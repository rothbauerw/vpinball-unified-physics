// license:GPLv3+
//
// A rollover button of the contact system (to-do 11.2, runs/2026-11-10-rollover/design.md section 2S): a cap on a centre peg in a round hole in the playfield, held up by the leaf
// blade of its switch. Sink only (user 2026-11-11, from a slow-motion video: no tilt, no sideways move): one body on one coordinate (bodypart.h), a slide along the cap's axis,
// square to the playfield: s = the cap's sink, 0 at rest (the push nut on the bushing's underside), rest_mm flush with the wood (the top's underside on the bushing's top).
//   spring   the leaf blade under the peg, linear in sink: spring_n up at rest (preload), spring_full_n at full sink; the cap's weight along the axis is its own term
//   stops    up at s 0 and down at s rest_mm: BodyStop, rigid (no unknown while the body is on one), the impact law with a constant restitution stop_e, the rest rule (bodypart.h)
//   surface  a flat round top, a vertical side cap_h_mm down from the top, the corner rounded with corner_mm: JetBumper::Revolve with 3 vertices, carried along the slide (one column, down).
//            A Star rollover (runs/2026-11-11-star/design.md) is the same model with another top (Top::Star): the dome through the stock star mesh's ridge tops, its tips' top edge
//            rounded with corner_mm, its gentle kinks at the probe's edge minimum (StarOutline, StarRounding); default mass kStarMassG
//   switch   closed while s >= close_mm, open below it (a plain leaf: no differential travel). Hit on close, Unhit on open
//   material the bumper skirt's plastic set (contactbumper.h Params), ESTIMATED, for the ball's contact and the stops' Hertz stiffness (the rest rule's squash)
// No friction, no tilt (design 2S). The table side is in trigger.cpp. Pure SI, world frame in metres, z up (the playfield's normal).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "physics/bodypart.h"
#include "physics/bumpergroup.h"
#include "physics/contactmovingpart.h"

namespace ButtonCap
{
   using GroupStep::Hit;
   using GroupStep::V3;

   struct Settings
   {
      double restMm = 3.175;       // rest_mm: the cap's top above the playfield at rest, and its sink to flush (Williams A-7385 / A-7903: 1/8 in, user 2026-11-10)
      double closeMm = 2.38125;    // close_mm: the switch is closed while the cap has sunk at least this far, open below it (3/32 in, user 2026-11-11; PartData S9)
      double springN = 0.2;        // spring_n: the leaf's force up on the cap at rest, N (its preload; the cap's weight not included; user 2026-11-11 after watching the drive: 0.34 too strong)
      double springFullN = 0.25;   // spring_full_n: the leaf's force up at full sink (flush), N; linear between (user 2026-11-11)
      double cornerMm = 0.79375;   // corner_mm: the rounding of the edge between the top and the side (1/32 in, user 2026-11-10); at most the cap's radius and rest_mm (larger: a notice, and the largest that fits)
      double massG = 1.0;          // mass_g: the cap with its peg, g (user 2026-11-10); a Star's default is kStarMassG (the owner sets it before reading the ini)
      double stopE = 0.2;          // stop_e: restitution of both stops (push nut up, flush down), the same at every speed (ESTIMATED from the user's video: no visible bounce; design 2S)
      double capHMm = 3.97;        // cap_h_mm: the cap's side from its top down (5/32 in, ESTIMATED, design 1H): the touchable shape; reaches into the hole at rest (must exceed rest_mm)
   };

   // The cap's top: the Button's flat round top, or the Star's dome (runs/2026-11-11-star/design.md 1b). Nothing else of the model differs
   enum class Top { Button, Star };
   constexpr double kStarMassG = 2.5;   // mass_g's default for a Star, g (PartData S10: a star rollover 2-3 g, user estimate 2026-11-10)

   // The Star's hit surface (design 1b; APPROXIMATION 1c: the smooth dome, not the spokes): a body of revolution through the ridge tops of the stock star mesh (triggerStarMesh.h: 10 spokes
   // from a centre plateau to the tips), (radius, height) in mesh units: the radius x R in plan; the height x rest / kStarTop, so the top is rest_mm above the wood at rest whatever R
   // (x 0.806 of the mesh at R 13.5 mm). The last point is the bottom of the tips' wall, below the wood (in the hole)
   constexpr double kStarProfile[6][2] = { { 0.0, 0.292 }, { 0.19, 0.292 }, { 0.47, 0.256 }, { 0.79, 0.177 }, { 0.97, 0.105 }, { 0.97, -0.037 } };
   constexpr double kStarTop = 0.292;   // the mesh's top (its plateau), mesh units
   constexpr int kStarTip = 4;          // the tips' top edge (rounded corner_mm); corners 1-3 are the dome's kinks (rounded at the probe's edge minimum)

   // An outline in (radius from the axis, height), m, starting on the axis, with each inner corner rounded by its own radius, as a Revolve sweep. Every straight part of the outline is a
   // segment at a constant sweep radius parallel to it, every corner the circle of its radius tangent to both sides (the Revolve's normal is then the outline's: no taper on a
   // surface segment). Where the radius changes between two corners, a joint along the side's normal joins the two circles tangent to that side at the same point: each circle of the
   // joint lies inside the larger one, so it is never the nearest part for a point outside the cap. The start takes the radius of corner 1; the end is rounded as the Button's side
   // (its circle's bottom on the last point) with the radius of the corner before it
   struct RoundedOutline
   {
      struct P2 { double r, z; };
      std::vector<P2> p;             // the points
      std::vector<P2> d, n;          // each side's direction and its outward normal (left of the direction: up on a side going out)
      std::vector<double> len;       // each side's length
      std::vector<double> tanHalf;   // each inner corner's tan(turn / 2): its rounding radius times it is the length the rounding takes off each side (0 at the ends)

      explicit RoundedOutline(std::vector<P2> pts) : p(std::move(pts))
      {
         for (size_t i = 0; i + 1 < p.size(); i++)
         {
            const double dr = p[i + 1].r - p[i].r, dz = p[i + 1].z - p[i].z, l = std::sqrt(dr * dr + dz * dz);
            d.push_back({ dr / l, dz / l });
            n.push_back({ -dz / l, dr / l });
            len.push_back(l);
         }
         tanHalf.assign(p.size(), 0.0);
         for (size_t i = 1; i + 1 < p.size(); i++)
         {
            const double c = d[i - 1].r * d[i].r + d[i - 1].z * d[i].z, s = std::fabs(d[i - 1].r * d[i].z - d[i - 1].z * d[i].r);
            tanHalf[i] = s / (1.0 + c);
         }
      }
      // rad[i]: corner i's rounding (1 .. size - 2; the ends' entries are not used)
      std::vector<JetBumper::Revolve::Vtx> Sweep(const std::vector<double> &rad) const
      {
         const size_t m = p.size() - 1;
         std::vector<JetBumper::Revolve::Vtx> v;
         double a = rad[1];
         v.push_back({ p[0].r - a * n[0].r, p[0].z - a * n[0].z, a });
         for (size_t i = 1; i < m; i++)
         {
            const double b = rad[i], k = 1.0 + n[i - 1].r * n[i].r + n[i - 1].z * n[i].z;
            const P2 c { p[i].r - b * (n[i - 1].r + n[i].r) / k, p[i].z - b * (n[i - 1].z + n[i].z) / k };   // the circle tangent to both sides
            if (b != a)
               v.push_back({ c.r + (b - a) * n[i - 1].r, c.z + (b - a) * n[i - 1].z, a });   // the joint: the side's own radius over the corner's tangent point
            v.push_back({ c.r, c.z, b });
            a = b;
         }
         v.push_back({ p[m].r - a * (n[m - 1].r + d[m - 1].r), p[m].z - a * (n[m - 1].z + d[m - 1].z), a });
         return v;
      }
   };
   // The Star's outline at rest for radius R and rest height rest (m)
   inline RoundedOutline StarOutline(const double R, const double rest)
   {
      std::vector<RoundedOutline::P2> pts;
      for (const auto &q : kStarProfile)
         pts.push_back({ q[0] * R, q[1] / kStarTop * rest });
      return RoundedOutline(std::move(pts));
   }
   // The Star's roundings that fit its sides (design 3, a tiny star: the largest that fits, as the Button's corner check): the three kinks share kink (the probe's edge minimum) unless
   // their sides are too short; then the tips' edge takes corner (corner_mm) unless the side before it or the wall is too short. kinkFit, cornerFit: the largest that fit, m (the owner's notices)
   struct StarRound
   {
      std::vector<double> rad;
      double kinkFit = 0.0, cornerFit = 0.0;
   };
   inline StarRound StarRounding(const RoundedOutline &o, const double kink, const double corner)
   {
      const std::vector<double> &t = o.tanHalf, &L = o.len;
      StarRound s;
      s.kinkFit = std::min({ L[0] / t[1], L[1] / (t[1] + t[2]), L[2] / (t[2] + t[3]), L[3] / t[3] });
      const double k = std::min(kink, s.kinkFit);
      s.cornerFit = std::min((L[3] - k * t[3]) / t[kStarTip], L[kStarTip] / (t[kStarTip] + 1.0));   // the wall: the rounding's tangent length and the end circle's radius
      s.rad = { 0.0, k, k, k, std::min(corner, s.cornerFit), 0.0 };
      return s;
   }

   // its own state beside the body's (PartRecords)
   struct CapState
   {
      V3 g = V3(0.0, 0.0, -9.81);   // the gravity the part feels (the owner's, between cycles; kept so a replay runs on the same)
      bool switchClosed = false;
   };

   class Part : public GroupStep::BodyPart, public CapState
   {
   public:
      Part()
      {
         KeepState<CapState>(this);
         KeepQueue(&swLog);
      }
      static constexpr double kStopReffM = 0.4e-3;   // the stops' Hertz radius, m: ESTIMATE (the gate's and wire switch's stop edge, stop_edge_mm 0.4); sets only the rest rule's squash
      static constexpr double kEdgeMinM = 0.2e-3;    // the smallest rounding a probe counts (the bumper's edgeMinMm)
      RubberContact::RubberParams rp;   // ball contact law
      JetBumper::Revolve prof;          // the cap at rest, in (radius from the axis, height above the playfield) round its axis
      V3 centre;                        // the axis on the playfield (the trigger's surface), m
      double R = 0.0, rc = 0.0, travel = 0.0, closeAt = 0.0, f0 = 0.0, k = 0.0;   // radius, corner, full sink, switch point (m); preload (N), blade rate (N/m)
      int stopUp = -1, stopDown = -1;
      std::vector<char> swLog;          // what happened, in order, until the owner reads it: 1 switch closed, 0 opened

      // axis: the cap's axis on the playfield, m; radius: the cap's radius (the trigger's Radius), m; s: corner_mm and cap_h_mm already checked (the owner's notices); top: its shape
      void Build(const Settings &s, const V3 &axis, const double radius, const RubberContact::BallSI &b, const Top top = Top::Button)
      {
         const double mm = 1e-3;
         centre = axis;
         R = radius;
         rc = s.cornerMm * mm;
         travel = s.restMm * mm;
         closeAt = s.closeMm * mm;
         mass = s.massG * 1e-3;
         reach = 1.0;        // a slide: the rate is the surface's speed
         maxMove = travel;
         f0 = s.springN;
         k = (s.springFullN - s.springN) / travel;
         // (0, top - rc), (R - rc, top - rc), (R - rc, top - cap_h + rc), each swept by rc: the flat top, the torus corner, the side down cap_h from the top (below the wood: the hole)
         if (top == Top::Star)   // the dome through the star's ridge tops, its wall down into the hole (cap_h_mm not used): the tips' edge rounded rc, the kinks at the edge minimum
         {
            const RoundedOutline o = StarOutline(R, travel);
            prof.v = o.Sweep(StarRounding(o, kEdgeMinM, rc).rad);
         }
         else
            prof.v = { { 0.0, travel - rc, rc }, { R - rc, travel - rc, rc }, { R - rc, travel - s.capHMm * mm + rc, rc } };
         const JetBumper::Params pp;   // the skirt's plastic set: E, nu, restitution at 0.5 m/s, friction
         rp = RubberContact::RubberParams();
         rp.eMpa = pp.plasticEMpa;
         rp.nu = pp.plasticNu;
         rp.thicknessMm = 100.0;
         rp.betaLayer = 0.0;
         rp.mu = pp.plasticMu;
         rp.twist = true;
         rp.normalLoss = false;
         rp.gammaKK = RubberContact::CalibrateKK(rp, b.R, b.m, 0.5, pp.plasticE);
         // the stops: rigid; Hertz k = 4/3 E* sqrt(reff) of plastic on plastic (E* = E / (2 (1 - nu^2))) for the rest rule's squash; restitution stop_e at every speed (every entry of the
         // law's table the same: its interpolation gives stop_e)
         const double estar = pp.plasticEMpa * 1e6 / (2.0 * (1.0 - pp.plasticNu * pp.plasticNu));
         q = rate = 0.0;
         const auto stop = [&](const double at, const int side) {
            GroupStep::BodyStopPlace st;
            st.at = at;
            st.side = side;
            st.lever = 1.0;
            st.k = (4.0 / 3.0) * estar * std::sqrt(kStopReffM);
            for (double &e : st.law.eTab)
               e = s.stopE;
            return st;
         };
         stopUp = PlaceStop(stop(0.0, -1));
         stopDown = PlaceStop(stop(travel, 1));
         switchClosed = false;
         Settle();
      }
      // At rest on the up stop from the start when the leaf holds it there (its preload above the weight): no unknown, asleep (design 2S rest)
      void Settle()
      {
         q = rate = 0.0;
         pinStop = StoredForce(0.0) < 0.0 ? stopUp : -1;
         stops[stopUp].touching = pinStop == stopUp;
         CycleBegin();
      }
      double Sink() const { return q; }
      double Weight() const { return -mass * g.z; }   // the cap's weight along the axis (down, + s), N

      // ---- the body's terms: V(s) = (f0 + k s / 2) s - W s ----
      double StoredEnergy(const double at) const override { return (f0 + 0.5 * k * at) * at - Weight() * at; }
      double StoredChord(const double at, const double dq, double, double &slope, double &lossF, double &lossSlope) const override
      {
         slope = 0.5 * k;
         lossF = lossSlope = 0.0;
         return f0 + k * (at + 0.5 * dq) - Weight();
      }
      double StoredForce(const double at) const override { return Weight() - f0 - k * at; }
      double FrictionLimit(double, double, const V3 &) const override { return 0.0; }   // no friction (design 2S)
      void AfterMove() override
      {
         const bool c = q >= closeAt;
         if (c != switchClosed)
         {
            switchClosed = c;
            swLog.push_back(c ? 1 : 0);
         }
      }
      double MaxStoredForce() const { return f0 + k * travel + std::fabs(Weight()); }

      // The cap's surface for a ball of radius Rb centred at c, the cap moved down by dq more (the trial): the Revolve's closest point; the column: the slide, down.
      // A ball rolling over the rounded corner turns round the corner's centre line (the distance to it, as a tube); on the flat top and the side it does not turn
      void ProbeSurface(const V3 &c, const double Rb, const double dq, Hit &h) const
      {
         const double x = c.x - centre.x, y = c.y - centre.y, rho = std::sqrt(x * x + y * y);
         const double cphi = rho > 1e-9 ? x / rho : 1.0, sphi = rho > 1e-9 ? y / rho : 0.0;
         const JetBumper::Revolve::Near nr = prof.Closest(rho, c.z - centre.z + q + dq, kEdgeMinM);
         h.delta = Rb - nr.g;
         h.n = V3(nr.nr * cphi, nr.nr * sphi, nr.nz);
         h.reff = BumperGroup::KappaEff(Rb, nr);
         h.curvDist = nr.edge ? nr.g + nr.rr : 0.0;
         h.grad = 1.0;
         const V3 down(0.0, 0.0, -1.0);
         h.nc = 1;
         h.c[0].Set(0, h.n.Dot(down), down, V3());
      }
      std::string HistoryNote() const override
      {
         char line[200];
         snprintf(line, sizeof line, "sink %.6f mm rate %.5g m/s KE %.5g stored %.9g uJ pinStop %d switch %d", q * 1e3, rate, KineticJ() * 1e6, StoredEnergy(q) * 1e6, PinStop(), (int)switchClosed);
         return line;
      }
   };

   // the engine side (PartRecords: part of its state)
   struct EngineState
   {
      bool idle = false;
      GroupStep::RestWatch rest;   // the rest rule's watch, alone and over its last cycle in a group (RestPose)
   };

   // The cap as a moving part of the contact system: the probe for engagement and the overlap replay, the speed bound, its surface, its own steps when no ball is near, its switch events.
   // A disabled trigger's cap still moves and is still touched (a disconnected switch is still a bump, design 4e): the table side drops its events
   class Cap : public ContactMovingPart, public EngineState
   {
   public:
      Cap() { KeepState<EngineState>(this); }
      Part part;
      RubberContact::BallSI ball;
      long failedSteps = 0, aloneSteps = 0, energySteps = 0;   // the part alone: a step not solved (rule 8), steps taken (cost), of them taken on energy alone (rule 8)
      int bookMs = 0;

      void Build(const Settings &s, const V3 &axis, const double radius, const RubberContact::BallSI &b, const Top top = Top::Button)
      {
         ball = b;
         part.Build(s, axis, radius, b, top);
         idle = true;   // held on its up stop: nothing to step until a ball comes near or the gravity changes
      }
      // The switch's changes as the part's events. The adapter calls it before the host takes the events
      void Collect()
      {
         for (const char e : part.swLog)
            events.push_back({ e == 1 ? Event::SwitchClosed : Event::SwitchOpened, -1, 0.0 });
         part.swLog.clear();
         part.stopEvents.clear();
      }

      // ---- the moving-part interface (contactmovingpart.h) ----
      void SetBall(const RubberContact::BallSI &b) override { ball = b; }
      void SetGravity(const V3 &g) override
      {
         if ((g - part.g).Length() > 1e-6)
            idle = false;   // (a nudge wakes it)
         part.g = g;
      }
      V3 Origin() const override { return part.centre; }
      double Reach() const override { return part.R + 2.0 * part.travel + 2.0 * part.rc; }   // every point of the cap at any sink
      double MotionRange() const override { return part.travel; }
      double SurfaceSpeedBound(const double cycleS) const override { return std::fabs(part.rate) + part.MaxStoredForce() / part.mass * cycleS; }
      RubberContact::TubeProbe Probe(const V3 &c) const override
      {
         Hit h;
         part.ProbeSurface(c, ball.R, 0.0, h);
         RubberContact::TubeProbe pr;
         pr.seg = 0;
         pr.normal = h.n;
         pr.dist = ball.R - h.delta;
         pr.point = c - h.n * pr.dist;
         pr.rho = 0.0;
         pr.delta = h.delta;
         pr.overlap = h.delta > 0.0;
         pr.reff = h.reff;
         return pr;
      }
      void SurfaceVelocity(const V3 &, V3 &v, V3 &w) const override
      {
         v = V3(0.0, 0.0, -part.rate);
         w = V3();
      }

      // ---- the group's calls ----
      GroupStep::Part *GroupPart() override { return &part; }
      int GroupSurfaces() const override { return 1; }
      void GroupProbe(int, const V3 &c, const double R, const double *const trial, Hit &h) const override { part.ProbeSurface(c, R, trial != nullptr ? trial[0] : 0.0, h); }
      const RubberContact::RubberParams *GroupLaw(int, double &layerM) const override
      {
         layerM = 0.0;
         return &part.rp;
      }
      void CycleBegin() override
      {
         idle = false;
         part.CycleBegin();
      }
      void TouchBegan(int, int, double) override {}   // (the button's Hit is its switch closing, not a touch)
      bool Sleeping() const override { return idle; }
      bool RestNumbers(double &moveM, double &speedMS) const override
      {
         moveM = rest.CycleMoved(RestPose());
         speedMS = std::fabs(part.rate);
         return true;
      }
      bool AtRest() const override { return idle || rest.Rests(RestPose(), std::fabs(part.rate)); }   // the one rest rule: its move over its last cycle in a group, and its speed now
      GroupStep::RestWatch *Rest() override { return &rest; }
      GroupStep::RestWatch::Pose RestPose() const override { return GroupStep::RestWatch::Pose().Add(part.q, 1.0); }   // its one rest pose: the sink
      // No ball near: the same group with no balls (the cap rises after a ball leaves), asleep while nothing moves, its state kept as it is
      int AdvanceAlone(const double dt) override
      {
         if (dt <= 0.0 || idle)
            return 0;
         const bool was = idle;
         part.ClearLoads();
         part.CycleBegin();
         const int steps = GroupStep::AloneCycle(part, dt, bookMs, failedSteps, aloneSteps, energySteps);
         idle = rest.Step(RestPose(), std::fabs(part.rate), false, restTouched, dt);
         GroupStep::RestWatch::Book(1, was, idle);
         return steps;
      }

      std::string StateNote() const override
      {
         char buf[200];
         std::snprintf(buf, sizeof buf, "button: sink %.4f mm at %.4f m/s, on stop %d, switch %d, asleep %d", part.q * 1e3, part.rate, part.PinStop(), (int)part.switchClosed, (int)idle);
         return buf;
      }
   };
}
