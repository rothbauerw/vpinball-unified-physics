// license:GPLv3+
//
// The jet bumper (contactbumper.h) as a moving part of the contact system (contactmovingpart.h): one host with three surfaces (skirt, ring, centre body); the engine's per-ball
// shear states of the host's contact (SurfaceContact::st, stExtra) carry one state per surface.
//
// World frame in metres, z up; the part is built at the foot of its axis (x, y, playfield height).

#pragma once

#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>

#include "physics/contactmovingpart.h"
#include "physics/contactbumper.h"

namespace JetBumper
{
   class EnginePart : public ContactMovingPart
   {
   public:
      Bumper model;
      std::string spikeNote;   // set when a ball's contact force on the bumper exceeds kSpikeN: the host logs it and raises a notification (rule 8)
      int spikes = 0;
      V3 origin;               // the foot of the bumper's axis on the playfield, world (m)
      RubberContact::BallSI ball;

      void Build(const Params &prm, const V3 &foot)
      {
         model.Configure(prm);
         origin = foot;
         m_reach = 0.5 * prm.skirtOdMm * MM + 2.0e-3;
      }

      // ---- ContactMovingPart ----
      void SetBall(const RubberContact::BallSI &b) override { ball = b; }
      void SetGravity(const V3 &g) override { model.SetGravity(g); }
      V3 Origin() const override { return origin; }
      double Reach() const override { return m_reach; }
      // the ring: its speed now plus what the coil at full current in the closed position and the return spring can add in a cycle; the skirt: its speed (the rim moves faster than the centre)
      // plus what its springs can add. The bumper can start at any time (Solenoid), so even asleep the ring's acceleration counts.
      double SurfaceSpeedBound(double cycleS) const override
      {
         const Params &p = model.p;
         const double ring = std::abs(model.vr) + RingAccelBound() * cycleS;
         const double rim = 0.5 * p.skirtOdMm * MM;
         const double skirt = std::abs(model.vz) + std::sqrt(model.wx * model.wx + model.wy * model.wy) * rim + SkirtAccelBound() * cycleS;
         return std::max(ring, skirt);
      }

      // the nearest surface as the engine's probe: dist = the ball centre's distance from the surface, rho 0 (the surface is a sheet, its thickness is in the profile)
      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         RubberContact::TubeProbe pr;
         const V3 rel = c - origin;
         if (std::sqrt(rel.x * rel.x + rel.y * rel.y) > m_reach + ball.R + 0.05 || rel.z > 0.06 + ball.R)
         {
            pr.dist = INFINITY;
            pr.delta = -INFINITY;
            return pr;
         }
         const Nearest nr = model.NearestSurface(rel, ball.R);
         pr.seg = nr.surf;
         pr.normal = nr.n;
         pr.point = nr.point + origin;
         pr.dist = nr.g + ball.R;
         pr.rho = 0.0;
         pr.delta = -nr.g;
         pr.overlap = nr.g < 0.0;
         pr.reff = nr.reff;
         return pr;
      }

      RubberContact::ContactResult ContactItem(int /*ballIndex*/, SurfaceContact &it, const V3 &c, const V3 &v, const V3 &w, double h) override
      {
         RubberContact::ContactResult r;
         BallContactState st;
         st.skirt = it.st;
         st.ring = it.stExtra[0];
         st.body = it.stExtra[1];
         BallLoad L;
         model.BallForces(ball, c - origin, v, w, h, st, L);
         it.st = st.skirt;
         it.stExtra[0] = st.ring;
         it.stExtra[1] = st.body;
         const double fn = std::max(L.fnSkirt, std::max(L.fnRing, L.fnBody));
         if (fn <= 0.0)
            return r;
         if (fn > kSpikeBaseN + kSpikePerMps * (v.Length() + SurfaceSpeedBound(0.0)) && spikeNote.empty() && spikes < 20)   // relative closing speed: the ball's plus the moving surface's
         {
            spikes++;
            const V3 rel = c - origin;
            char buf[400];
            std::snprintf(buf, sizeof buf, "contact force %.0f N (skirt %.0f ring %.0f body %.0f), overlap %.0f um, ball at %.1f mm from the axis, %.1f mm high, speed %.2f m/s; ", fn, L.fnSkirt, L.fnRing, L.fnBody, L.delta * 1e6,
               std::sqrt(rel.x * rel.x + rel.y * rel.y) * 1e3, rel.z * 1e3, v.Length());
            spikeNote = std::string(buf) + StateNote();
         }
         r.inContact = true;
         r.force = L.force;
         r.torqueOnBall = L.torque;
         r.fn = fn;
         r.normal = L.normal;
         r.delta = L.delta;
         return r;
      }

      std::string StateNote() const override
      {
         char buf[300];
         std::snprintf(buf, sizeof buf, "bumper: ring %.2f mm down at %.2f m/s, spoon %.2f mm, skirt %.2f mm high tilt %.2f/%.2f deg, switch %d, coil %.2f A, asleep %d", model.dr * 1e3, model.vr, model.ds * 1e3, model.zc * 1e3,
            model.RotXDeg(), model.RotYDeg(), (int)model.switchClosed, model.iCoil, (int)model.asleep);
         return buf;
      }

      void SaveState(State &out) const override
      {
         out.used = 0;
         Put(out, steppedThisCycle);
         Put(out, steppedS);
         Bumper::Dyn d;
         std::memset(&d, 0, sizeof d);   // (the padding too: saved bytes must be the same for the same state)
         model.SaveDyn(d);
         Put(out, d);
      }
      void RestoreState(const State &in) override
      {
         StateReader r(in);
         r.Get(steppedThisCycle);
         r.Get(steppedS);
         Bumper::Dyn d;
         r.Get(d);
         model.LoadDyn(d);
         events.clear();
         spikeNote.clear();
      }

      void SurfaceVelocity(const V3 &, V3 &v, V3 &w) const override
      {
         v = V3();
         w = V3();
      }

      // The step a ball's contact needs (design 7): the reduced mass of the ball and the surface's effective mass (the ring assembly 105 g, the skirt 7 g, the body fixed), the Hertz
      // stiffness at the deepest compression the closing speed can reach, ending at the arrival; the ring may be moving to meet the ball (a bound on its speed and acceleration).
      double BallStep(const V3 &c, const V3 &v, double hMax) const override
      {
         const V3 rel = c - origin;
         const double rho = std::sqrt(rel.x * rel.x + rel.y * rel.y);
         if (rho > m_reach + ball.R + std::max(v.Length(), 0.1) * hMax + 0.01 || rel.z > 0.06 + ball.R)
            return hMax;
         const Nearest nr = model.NearestSurface(rel, ball.R);
         const double mEffSurf = nr.surf == 1 ? model.p.ringMass : nr.surf == 0 ? model.p.skirtMass : 1e3;
         const double mu = 1.0 / (1.0 / ball.m + 1.0 / mEffSurf);
         const double eStar = nr.surf == 1 ? kSteelEstar : model.p.plasticEMpa * 1e6 / (1.0 - model.p.plasticNu * model.p.plasticNu);
         const double k = (4.0 / 3.0) * eStar * std::sqrt(nr.reff > 0.0 ? nr.reff : ball.R);
         // closing speed: the ball's, plus the surface's own bound (the ring may already be on its way down at up to ~2 m/s: the coil stroke gives 1 - 1.75 m/s)
         const double vc = std::max(0.0, -v.Dot(nr.n)) + (nr.surf == 2 ? 0.0 : SurfaceSpeedBound(hMax));
         const double ac = 9.81 + (nr.surf == 1 ? RingAccelBound() : nr.surf == 0 ? SkirtAccelBound() : 0.0);
         const double delta = -nr.g;
         auto hContact = [&](double vImpact) {
            const double vi = std::max(vImpact, 0.0);
            const double dMax = std::max(delta, 0.0) + std::pow(1.25 * mu * vi * vi / k, 0.4);
            const double kt = 1.5 * k * std::sqrt(std::max(dMax, 1e-12));
            return 0.1 / std::sqrt(kt / mu);
         };
         double h = hMax;
         if (delta > 0.0)
            h = std::min(h, hContact(vc));
         else
         {
            const double tArr = (-vc + std::sqrt(vc * vc + 4.0 * ac * -delta)) / (2.0 * ac);
            if (tArr < h)
               h = std::max(tArr, hContact(vc + ac * tArr));
         }
         return std::max(h, kMinStep);
      }

      // one island substep h: the balls' reactions of this substep act during all of it, the bumper's bodies advance in their own steps within it
      void StepInIsland(double h) override
      {
         model.StepHeld(h, model.OwnStep());
         Collect();
         steppedThisCycle = true;
         steppedS += h;
      }
      // no ball is near: one physics cycle
      int AdvanceAlone(double dt) override
      {
         const int n = model.StepHeld(dt, model.OwnStep());
         if (n > 0)
            Collect();
         return n;
      }

   private:
      double RingAccelBound() const
      {
         const Params &p = model.p;
         return (model.coilD.strength + p.returnPreloadN + p.returnKNPerM * p.ringStrokeMm * MM) / p.ringMass;
      }
      double SkirtAccelBound() const
      {
         const Params &p = model.p;
         return (p.skirtPreloadN + p.skirtKNPerM * p.skirtClearanceMm * MM) / p.skirtMass;
      }
      double m_reach = 0.0;
      static constexpr double kMinStep = 1e-7;
      // a ball hitting the steel ring peaks at about 1.45 kN per m/s of closing speed (offline T6: 1.4 kN at 1 m/s, 4.4 kN at 3 m/s); far more than that for the ball's speed is a deep first contact
      static constexpr double kSpikeBaseN = 500.0, kSpikePerMps = 2500.0;

      // the model's events as the part's (Hit, Switch closed / opened); the pulse and stop events are for the log only
      void Collect()
      {
         for (const JetBumper::Event &e : model.events)
         {
            if (e.kind == EvHit)
               events.push_back({ ContactMovingPart::Event::Hit, -1, e.value });
            else if (e.kind == EvSwitchClosed)
               events.push_back({ ContactMovingPart::Event::SwitchClosed, -1, e.value });
            else if (e.kind == EvSwitchOpened)
               events.push_back({ ContactMovingPart::Event::SwitchOpened, -1, e.value });
         }
         model.events.clear();
      }
   };
}
