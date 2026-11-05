// license:GPLv3+
//
// What every part on one horizontal axis shares as a part of the group step (a gate, a spinner; docs/GateSpinnerDesign.txt sections 3 and 6): its touchable shape (wire capsules and sheet triangles in body coordinates), its stored energy
// (gravity on the centre of mass, a switch blade), the pivot's friction. Its body is a BodyPart (bodypart.h); the moving is the group's. Pure SI, no VPX types.
//
// Frame: the axis through A along E (horizontal). N: the way a ball passes from the front; Z = E x N, up. Body coordinates: x along E, y down at theta = 0, z front at theta = 0; a positive theta swings the lower part toward N.
// Ball contact: an ordinary Item, the steel-on-steel Hertz + viscoelastic law, plus the bend: out of its plane the part bends like a cantilever from the axis, far softer than Hertz, a massless spring in series (Hit::kb), taken when a touch starts and held while it lasts.
// Pivot: Coulomb friction of the pin wire in its hole: |tau| <= mu r_pin |N|, N = the load on the pins (weight, balls, the centre of mass's centripetal force).

#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "physics/bodypart.h"
#include "physics/contacthingebody.h"
#include "physics/rubbercontact.h"

namespace HingeGroup
{
   using GroupStep::Hit;
   using GroupStep::V3;
   using RubberContact::BallSI;
   using RubberContact::RubberParams;
   using RubberContact::TubeProbe;

   class Part : public GroupStep::BodyPart
   {
   public:
      RubberParams rp;            // ball contact law
      V3 A, E, N, Z;              // axis point and frame
      V3 g = V3(0.0, 0.0, -9.81); // the gravity the part feels: gravity minus the cabinet acceleration
      HingeBody::MassProps mp;
      double rPin = 0.6e-3, muPivot = 0.3, cViscous = 0.0;
      bool pieceSum = false;      // [Contact] hinge_piece_sum (SetPieceSum): the contact as its union's terms, each stepped by its own exact difference (ProbeTerms)
      std::vector<HingeBody::Cam> cams;
      double kBlade = 0.0, aBlade = 0.0;   // a one-sided torsion spring: a blade pressed once theta is past aBlade
      std::vector<char> swLog;    // what happened, in order, until the owner reads it: 1 switch closed, 0 opened, 2 a spin
      bool switchClosed = false;
      double prevWrapped = 0.0;   // (a spinner with no switch: where it was against its switch angle)

      // ---- what the derived part supplies ----
      virtual double BendKb(double y) const = 0;           // the bend stiffness out of the part's plane at the distance y from the axis
      // the bend out of the part's plane, in series with a touch at pb (body coordinates) along the normal nb: the part's plane holds the axis (out of it: z, N)
      virtual double BendAlongNormal(const V3 &nb, const V3 &pb) const
      {
         const double oop = std::fabs(nb.z);
         return oop > 1e-6 ? BendKb(std::max(std::fabs(pb.y), 1e-4)) / (oop * oop) : 0.0;
      }
      virtual bool NearStart(double th) const = 0;         // near its start angle: a Hit fires only then
      virtual int LimitKind(int) const { return 0; }       // a touch of this stop as an event: +1 EOS, -1 BOS, 0 none
      virtual void UpdateSwitches() = 0;                   // after a step

      int Surfaces() const { return surfaces; }
      // Where gravity alone hangs it: on a sloped table not theta 0 (6 deg on a 6 deg table)
      double HangAngle() const { return std::atan2(GravB(), GravA()); }
      double PendulumArm() const { const V3 cm = mp.Cm(); return std::sqrt(cm.y * cm.y + cm.z * cm.z); }
      // The most the stored energy can push it with (the speed bound). The blade's largest torque is at the farthest angle the body can reach, k (theta_max - a): the
      // highest stop on the + side (a driven one at the farther of its place and target); a preloaded blade (a < 0, a wire switch) is pressed over its whole travel. No + stop: half a radian past a.
      double MaxStoredTorque() const
      {
         double thMax = aBlade + 0.5;
         bool bounded = false;
         for (const GroupStep::BodyStop &s : stops)
            if (s.side > 0)
            {
               const double stopAt = s.rate > 0.0 ? std::max(s.at, s.target) : s.at;
               thMax = bounded ? std::max(thMax, stopAt) : stopAt;
               bounded = true;
            }
         double tau = mp.m * g.Length() * PendulumArm() + kBlade * std::max(thMax - aBlade, 0.0);
         for (const HingeBody::Cam &cm : cams)
            tau += (cm.F0 + cm.k * cm.r) * cm.r;
         return tau;
      }
      // One contact's probe: only that side's wires and faces, the part at its angle now + dq
      void ProbeSurface(const int surface, const V3 &c, const double R, const double dq, Hit &h) const
      {
         const double th = q + dq;
         const V3 dn = Down(th), fr = Front(th), rel = c - A, cb(rel.Dot(E), rel.Dot(dn), rel.Dot(fr));
         BallSI b;
         b.R = R;
         if (pieceSum)
         {
            ProbeTerms(surface, cb, b, dn, fr, h);
            return;
         }
         TubeProbe best;
         best.dist = INFINITY;
         best.delta = -INFINITY;
         // nearest first, exactly: a wire's distance from the ball centre is at least the distance to its middle less half its length. The wire nearest by that bound is
         // probed first; then, in their own order, only those whose bound is within its distance (and rounding) can be the nearest: the same best as probing them all
         const int nw = (int)wires.size();
         segBound.resize((size_t)nw);
         int j = -1;
         for (int i = 0; i < nw; i++)
         {
            const Seg &w = wires[i];
            segBound[i] = w.surface == surface ? (cb - w.m).Length() - w.half : INFINITY;
            if (w.surface == surface && (j < 0 || segBound[i] < segBound[j]))
               j = i;
         }
         TubeProbe prJ;
         double cut = INFINITY;
         if (j >= 0)
         {
            const V3 pts[2] = { wires[j].a, wires[j].b };
            GroupStep::envStats.capsules++;
            prJ = RubberContact::ProbeTube(b, pts, 2, rho, cb);
            if (prJ.seg >= 0)
               cut = prJ.dist * (1.0 + 1e-9) + 1e-12;
         }
         int bestI = -1;
         for (int i = 0; i < nw; i++)
         {
            const Seg &w = wires[i];
            if (w.surface != surface || segBound[i] > cut)
               continue;
            const V3 pts[2] = { w.a, w.b };
            if (i != j)
               GroupStep::envStats.capsules++;
            TubeProbe pr = i == j ? prJ : RubberContact::ProbeTube(b, pts, 2, rho, cb);
            if (pr.seg < 0)
               continue;
            pr.delta = R + rho - pr.dist;
            if (pr.delta > best.delta)
            {
               best = pr;
               bestI = i;
            }
         }
         // diagnostics (VPX_GROUP_ROUND_TRACE): the nearest piece and the second deepest (read only: no state, no counts)
         const bool noteOn = ProbeNoteOn();
         char note2[200] = "";
         if (noteOn)
         {
            const auto sAlong = [&](const int i) { const V3 t = wires[i].b - wires[i].a; return (cb - wires[i].a).Dot(t) / t.Dot(t); };
            int i2 = -1;
            TubeProbe p2;
            for (int i = 0; i < nw; i++)
            {
               if (wires[i].surface != surface || i == bestI)
                  continue;
               const V3 pts[2] = { wires[i].a, wires[i].b };
               TubeProbe pr = RubberContact::ProbeTube(b, pts, 2, rho, cb);
               if (pr.seg < 0)
                  continue;
               pr.delta = R + rho - pr.dist;
               if (i2 < 0 || pr.delta > p2.delta)
               {
                  i2 = i;
                  p2 = pr;
               }
            }
            if (i2 < 0)
               snprintf(note2, sizeof note2, " | 2nd none");
            else
            {
               const double s2 = sAlong(i2);
               snprintf(note2, sizeof note2, " | 2nd piece %d s %.4f depth %.4g um", i2, s2, p2.delta * 1e6);
            }
            if (bestI >= 0)
               snprintf(probeNote, sizeof probeNote, "piece %d s %.4f depth %.4g um", bestI, sAlong(bestI), best.delta * 1e6);
            else
               snprintf(probeNote, sizeof probeNote, "piece none");
         }
         if (!tris.empty() && surface == 0)
         {
            const TubeProbe pr = RubberContact::ProbeTriangles(b, tris, cand, halfT, cb, nullptr, 0, halfT);
            if (pr.seg >= 0 && pr.delta > best.delta)
            {
               best = pr;
            }
         }
         const V3 nb = best.normal;
         h.delta = best.delta;
         h.n = E * nb.x + dn * nb.y + fr * nb.z;
         h.reff = best.reff > 0.0 ? best.reff : R;
         h.curvDist = 0.0;
         h.grad = 1.0;
         // the contact point (body coordinates), its speed per unit rate of theta (E x Down = Front, E x Front = -Down): the depth's slope with theta is the normal's part of it
         const V3 pb = cb - nb * (R - 0.5 * std::max(best.delta, 0.0));
         h.nc = 1;
         h.c[0].Set(0, nb.z * pb.y - nb.y * pb.z, fr * pb.y - dn * pb.z, E);
         h.kb = BendAlongNormal(nb, pb);
         if (noteOn)   // (what the probe returns: the deepest piece or a face, its body normal and slope)
         {
            const size_t n0 = strlen(probeNote);
            snprintf(probeNote + n0, sizeof probeNote - n0, " -> depth %.4g um n (%.3f %.3f %.3f) slope %.4g faces %d%s", h.delta * 1e6, nb.x, nb.y, nb.z, h.c[0].slope * h.grad,
                     (int)tris.size(), note2);
         }
      }

      // [Contact] hinge_piece_sum: the contact as its union's terms, each a smooth depth of its own: every touchable piece of this side (+1), every joint shared by
      // m of them (-(m - 1)), the faces as one term (+1). U = sum s_t U(max(d_t, 0)) at every place of the ball: a piece off the ball has d_i <= 0 and its joints
      // d_v <= d_i, so the fixed list is the union wherever the ball is. The group steps each term by its own exact difference (Hit::t); the Hit's own
      // depth, normal and slope are the union's (FoldSum; the deepest piece's while none presses), for the step choice and arrival.
      // Faces: one term per flat quad (AddPlate), minus a capsule along each line two touchable elements share (a fold edge, a wire along a plate edge;
      // FinishShape). Open: a wire end lying on a face has no overlap term (FinishShape's overlapWarn).
      const std::string &OverlapWarn() const { return overlapWarn; }   // (FinishShape: a wire end on a face, empty if none)
      bool SetPieceSum(const bool on)
      {
         BuildJoints();
         bool fits = true;
         for (int sf = 0; sf < surfaces; sf++)
         {
            int n = sf == 0 ? (int)(quads.size() + lines.size()) : 0;
            for (const Seg &w : wires)
               n += w.surface == sf ? 1 : 0;
            for (const Joint &j : joints)
               n += j.surface == sf ? 1 : 0;
            fits = fits && n <= GroupStep::kTermMax;
         }
         pieceSum = on && fits;
         return fits;
      }
      void ProbeTerms(const int surface, const V3 &cb, const BallSI &b, const V3 &dn, const V3 &fr, Hit &h) const
      {
         BuildJoints();
         const double R = b.R, S = R + rho;
         h.nt = 0;
         RubberContact::FoldSum sum;
         V3 bestN(0.0, 1.0, 0.0);
         double bestD = -INFINITY, bestReff = 0.0;
         const auto add = [&](const double d, const V3 &nb, const double s, const double reff) {
            const V3 pb = cb - nb * (R - 0.5 * std::max(d, 0.0));
            GroupStep::HitTerm &t = h.t[h.nt++];
            t.d = d;
            t.n = E * nb.x + dn * nb.y + fr * nb.z;
            t.slope = nb.z * pb.y - nb.y * pb.z;
            t.s = s;
            if (d > 0.0)
            {
               const double d15 = d * std::sqrt(d);
               sum.s25 += s * d * d15;
               sum.g = sum.g + nb * (s * d15);
            }
            if (s > 0.0 && d > bestD)
            {
               bestD = d;
               bestN = nb;
               bestReff = reff;
            }
         };
         for (const Seg &w : wires)
            if (w.surface == surface)
            {
               const V3 pts[2] = { w.a, w.b };
               GroupStep::envStats.capsules++;
               const TubeProbe pr = RubberContact::ProbeTube(b, pts, 2, rho, cb);
               if (pr.seg >= 0)
                  add(S - pr.dist, pr.normal, 1.0, pr.reff);
               else
                  add(-S, bestN, 1.0, 0.0);   // (a piece of no length: never touched, kept so the list is the same everywhere)
            }
         for (const Joint &j : joints)
            if (j.surface == surface)
            {
               const V3 cv = cb - j.v;
               const double dist = cv.Length();
               add(S - dist, dist > 1e-12 ? cv / dist : bestN, -(double)j.weight, 0.0);
            }
         if (surface == 0)
         {
            for (const Quad &qd : quads)   // each flat quad by its own two triangles (convex and flat: its nearest point is exact)
            {
               const TubeProbe pr = RubberContact::ProbeTriangles(b, tris, qd.cand, halfT, cb, nullptr, 0, halfT);
               add(pr.seg >= 0 ? pr.delta : -S, pr.seg >= 0 ? pr.normal : bestN, 1.0, pr.seg >= 0 ? pr.reff : 0.0);
            }
            for (const Line &l : lines)   // a shared line, taken back once: where both nearest points are on it, the smaller element's depth is its own and cancels
            {
               const V3 pts[2] = { l.a, l.b };
               GroupStep::envStats.capsules++;
               const TubeProbe pr = RubberContact::ProbeTube(b, pts, 2, l.r, cb);
               add(pr.seg >= 0 ? R + l.r - pr.dist : -S, pr.seg >= 0 ? pr.normal : bestN, -1.0, 0.0);
            }
         }
         TubeProbe u;
         u.delta = bestD;
         u.normal = bestN;
         u.grad = 1.0;
         const bool summed = sum.s25 > 0.0;
         if (summed)
            sum.Finish(u, 0.0);
         const V3 nb = u.normal;
         h.delta = u.delta;
         h.n = E * nb.x + dn * nb.y + fr * nb.z;
         h.reff = bestReff > 0.0 ? bestReff : R;
         h.curvDist = 0.0;
         h.grad = summed ? u.grad : 1.0;
         const V3 pb = cb - nb * (R - 0.5 * std::max(u.delta, 0.0));
         h.nc = 1;
         h.c[0].Set(0, nb.z * pb.y - nb.y * pb.z, fr * pb.y - dn * pb.z, E);
         h.kb = BendAlongNormal(nb, pb);
         GroupStep::envStats.terms += h.nt;
         if (ProbeNoteOn())
         {
            int on = 0;
            for (int t = 0; t < h.nt; t++)
               on += h.t[t].d > 0.0 ? 1 : 0;
            snprintf(probeNote, sizeof probeNote, "terms %d pressing %d deepest piece %.4g um -> depth %.4g um n (%.3f %.3f %.3f) slope %.4g summed %d", h.nt, on, bestD * 1e6,
                     h.delta * 1e6, nb.x, nb.y, nb.z, h.c[0].slope * h.grad, summed ? 1 : 0);
         }
      }

      // ---- the body's terms ----
      double Viscous() const override { return cViscous; }
      double StoredEnergy(const double th) const override
      {
         double e = GravV(th);
         const double x = th - aBlade;
         if (kBlade > 0.0 && x > 0.0)
            e += 0.5 * kBlade * x * x;
         for (const HingeBody::Cam &cm : cams)
            e += cm.V(cm.Delta(th));
         return e;
      }
      double StoredForce(const double th) const override
      {
         double tau = -GravVp(th);
         const double x = th - aBlade;
         if (kBlade > 0.0 && x > 0.0)
            tau -= kBlade * x;
         for (const HingeBody::Cam &cm : cams)
            tau += cm.F(cm.Delta(th)) * cm.r * std::sin(th - cm.phase);
         return tau;
      }
      double StoredChord(const double th, const double dq, const double h, double &slope, double &lossF, double &lossSlope) const override
      {
         // gravity: closed form, no cancellation (sin a - sin b = 2 cos((a + b) / 2) sin((a - b) / 2))
         const double half = 0.5 * dq, thm = th + half, sinc = std::fabs(half) > 1e-4 ? std::sin(half) / half : 1.0 - half * half / 6.0;
         double chord = GravVp(thm) * sinc;
         slope = 0.5 * GravVpp(thm);
         lossF = lossSlope = 0.0;
         if (kBlade > 0.0)
         {
            // the linear spring's chord in closed form, k (x0 + dq / 2), while both ends are past the free angle: (V1 - V0) / dq with x1 = x0 + dq loses ulp(x0) / dq to rounding,
            // which froze the solve of a ball resting on a wire switch (dq ~ 1e-12 rad). Across the free angle one end is 0, so the difference form has no cancellation.
            const double x0 = th - aBlade, x1 = x0 + dq;
            if (x0 > 0.0 && x1 > 0.0)
               chord += kBlade * (x0 + 0.5 * dq);
            else if (x0 > 0.0 || x1 > 0.0)
            {
               const double a0 = std::max(x0, 0.0), a1 = std::max(x1, 0.0);
               chord += 0.5 * kBlade * (a1 * a1 - a0 * a0) / dq;
            }
            if (x0 > 0.0 || x1 > 0.0)
               slope += 0.5 * kBlade;
         }
         for (const HingeBody::Cam &cm : cams)
         {
            const double d0 = cm.Delta(th), d1 = cm.Delta(th + dq);
            if (d0 <= 0.0 && d1 <= 0.0)
               continue;
            const double sn = std::sin(thm - cm.phase), ddt = -cm.r * sn;   // d delta / d theta at mid-step
            const double dd = ddt * sinc;                                    // (d1 - d0) / dq in closed form (cos a - cos b = -2 sin((a + b) / 2) sin((a - b) / 2))
            const double Fm = cm.F(std::max(0.5 * (d0 + d1), 0.0));
            chord += cm.MeanF(d0, d1) * dd;
            slope += 0.5 * (cm.k + (cm.F0 > 0.0 ? cm.F0 / cm.rampM : 0.0)) * ddt * ddt;
            if (cm.c > 0.0)
            {
               // the damper on the part of the move inside the zone: continuous in the move, and its work has one sign; never pulling.
               // Inside the zone at both ends the move in it is dd dq (closed form); entering or leaving, one end is 0 (no cancellation)
               const bool in = d0 > 0.0 && d1 > 0.0;
               const double p0 = std::max(d0, 0.0), p1 = std::max(d1, 0.0), dp = in ? dd * dq : p1 - p0;
               double Fd = cm.c * dp / h;
               if (Fd < -Fm)
                  Fd = -Fm;
               if (std::fabs(dq) > 1e-12)
               {
                  lossF += Fd * (in ? dd : dp / dq);
                  lossSlope += cm.c * ddt * ddt / h;
               }
            }
         }
         return chord;
      }
      double FrictionLimit(const double th, const double om, const V3 &ballForce) const override { return muPivot * rPin * PinLoad(th, om, ballForce); }
      // the pin on its hole wall: Mindlin's sideways spring 8 G* a of the Hertz patch, times r^2 (steel on steel, E* 115 GPa, G* 23.8 GPa)
      double PinStiffness(const double th, const double om, const V3 &ballForce) const override
      {
         const double N = PinLoad(th, om, ballForce), a = std::cbrt(3.0 * N * rPin / (4.0 * 115e9));
         return 8.0 * 23.8e9 * a * rPin * rPin;
      }
      double PinLoad(const double th, const double om, const V3 &ballForce) const   // the hole's force on the pin, across the axis, N
      {
         const V3 cm = mp.Cm(), rcm = Down(th) * cm.y + Front(th) * cm.z, load = g * mp.m + ballForce;
         V3 pin = rcm * (-mp.m * om * om) - load;
         pin = pin - E * pin.Dot(E);
         return pin.Length();
      }
      void AfterMove() override { UpdateSwitches(); }
      void SaveState() override
      {
         BodyPart::SaveState();
         svSwitch = switchClosed;
         svWrapped = prevWrapped;
         svSwLog = swLog.size();
      }
      void RestoreState() override
      {
         BodyPart::RestoreState();
         switchClosed = svSwitch;
         prevWrapped = svWrapped;
         swLog.resize(svSwLog);
      }
      std::string StartNote() const override
      {
         char s[200];
         snprintf(s, sizeof s, " switch %d wrapped %.17g g %.17g %.17g %.17g", (int)svSwitch, svWrapped, g.x, g.y, g.z);
         return BodyPart::StartNote() + s;
      }
      // diagnostics (VPX_GROUP_HISTORY_TRACE): the body's note and the blade's share of the stored energy (uJ)
      std::string HistoryNote() const override
      {
         const double x = q - aBlade;
         char s[120];
         snprintf(s, sizeof s, " blade %.9g uJ switch %d", kBlade > 0.0 && x > 0.0 ? 0.5 * kBlade * x * x * 1e6 : 0.0, (int)switchClosed);
         return BodyPart::HistoryNote() + s;
      }
      // diagnostics (VPX_GROUP_CONTACT_TRACE): the body's note, the blade's share of the stored torque at mid-move
      std::string ContactNote(const double *x, const int off, const double h, const double *Res) const override
      {
         const double qm = q + 0.5 * (NumDof() == 1 ? x[off] : 0.0);
         char s[200];
         snprintf(s, sizeof s, " | blade %.5g (k %.4g a %.4f deg) switch %d", kBlade > 0.0 && qm > aBlade ? -kBlade * (qm - aBlade) : 0.0, kBlade, aBlade * 57.29577951308232, (int)switchClosed);
         return BodyPart::ContactNote(x, off, h, Res) + s;
      }

   protected:
      struct Seg { V3 a, b; int surface; V3 m; double half; };   // m, half: its middle and half length (the probe's nearest-first bound)
      mutable std::vector<double> segBound;   // (ProbeSurface's per-wire bound, kept to save the allocation)
      mutable char probeNote[480] = "";       // diagnostics (VPX_GROUP_ROUND_TRACE): ProbeSurface's last probe, its pieces (LastProbeNote)
      static bool ProbeNoteOn() { static const bool on = []() { const char *const e = getenv("VPX_GROUP_ROUND_TRACE"); return e != nullptr && atoi(e) != 0; }(); return on; }

   public:
      std::string LastProbeNote() const override { return probeNote; }

   protected:
      void Begin(const V3 &axisPoint, const V3 &front, const V3 &gravity)
      {
         g = gravity;
         mp = HingeBody::MassProps();
         stops.clear();
         cams.clear();
         wires.clear();
         joints.clear();
         jointsBuilt = false;
         tris.clear();
         cand.clear();
         quads.clear();
         lines.clear();
         overlapWarn.clear();
         q = rate = 0.0;
         A = axisPoint;
         N = front / front.Length();
         Z = V3(0.0, 0.0, 1.0);
         E = N.Cross(Z);
         E = E / E.Length();
         N = Z.Cross(E);
         surfaces = 1;
         surfaceNow = 0;
         kBlade = aBlade = 0.0;
         switchClosed = false;
      }
      void AddWire(const V3 &p, const V3 &u, const double lambda, const bool touchable)
      {
         mp.AddWire(p, u, lambda);
         if (touchable)
            wires.push_back({ p, u, surfaceNow, (p + u) * 0.5, 0.5 * (u - p).Length() });
      }
      void AddPlate(const V3 &a, const V3 &b, const V3 &c, const V3 &d, const double sigma)
      {
         mp.AddQuad(a, b, c, d, sigma);
         tris.push_back(RubberContact::MakeSurfaceTri(a, b, c));
         tris.push_back(RubberContact::MakeSurfaceTri(a, c, d));
         Quad qd;
         qd.p[0] = a; qd.p[1] = b; qd.p[2] = c; qd.p[3] = d;
         qd.cand = { (int)tris.size() - 2, (int)tris.size() - 1 };
         quads.push_back(qd);
      }
      // the shape is complete: the faces' list, the body's inertia and reach
      void FinishShape(const double reachM)
      {
         cand.resize(tris.size());
         for (size_t i = 0; i < tris.size(); i++)
            cand[i] = (int)i;
         // hinge_piece_sum's shared lines: two quads with two corners in common (a fold edge, radius halfT); a touchable wire along a quad's edge (the stretch
         // they share, the smaller radius). A wire end lying on a face is not a line: its overlap would count twice, so it is named (overlapWarn).
         lines.clear();
         overlapWarn.clear();
         for (size_t i = 0; i < quads.size(); i++)
         {
            const Quad &qa = quads[i];
            for (size_t j = i + 1; j < quads.size(); j++)
            {
               V3 sh[2];
               int ns = 0;
               for (const V3 &p : qa.p)
                  for (const V3 &o : quads[j].p)
                     if ((p - o).Length() < 1e-9 && ns < 2)
                        sh[ns++] = p;
               if (ns == 2)
                  lines.push_back({ sh[0], sh[1], halfT });
            }
            for (const Seg &w : wires)
            {
               const V3 u = w.b - w.a;
               const double L2 = u.Dot(u);
               if (w.surface != 0 || L2 <= 0.0)
                  continue;
               for (int e = 0; e < 4; e++)
               {
                  const V3 &p = qa.p[e], &o = qa.p[(e + 1) % 4];
                  const double tp = (p - w.a).Dot(u) / L2, to = (o - w.a).Dot(u) / L2;
                  if ((p - w.a - u * tp).Length() > 1e-9 || (o - w.a - u * to).Length() > 1e-9)
                     continue;
                  const double t0 = std::max(0.0, std::min(tp, to)), t1 = std::min(1.0, std::max(tp, to));
                  if ((t1 - t0) * std::sqrt(L2) > 1e-9)
                     lines.push_back({ w.a + u * t0, w.a + u * t1, std::min(rho, halfT) });
               }
            }
         }
         for (const Seg &w : wires)
            for (const V3 &v : { w.a, w.b })
               for (const Quad &qd : quads)
               {
                  if (w.surface != 0)
                     continue;
                  BallSI pt;
                  pt.R = 0.0;
                  const TubeProbe pr = RubberContact::ProbeTriangles(pt, tris, qd.cand, halfT, v, nullptr, 0, halfT);
                  if (pr.seg < 0 || pr.dist >= rho + halfT)
                     continue;
                  bool onLine = false;
                  for (const Line &l : lines)
                  {
                     const V3 lu = l.b - l.a;
                     const double s = std::min(std::max((v - l.a).Dot(lu) / std::max(lu.Dot(lu), 1e-300), 0.0), 1.0);
                     onLine = onLine || (v - l.a - lu * s).Length() < 1e-9;
                  }
                  if (!onLine)
                  {
                     char s[200];
                     snprintf(s, sizeof s, "a wire end at (%.2f %.2f %.2f) mm is %.3f mm from a plate face: their overlap counts twice", v.x * 1e3, v.y * 1e3, v.z * 1e3, pr.dist * 1e3);
                     overlapWarn = s;
                  }
               }
         mass = mp.I;
         reach = reachM;
      }
      // the ball contact law: steel on steel (E* of the pair in the law's single-body form), restitution of the part at 1 m/s on a fixed surface
      void SetContactLaw(const double elasticity, const double friction, const BallSI &b)
      {
         rp = RubberParams();
         rp.nu = 0.3;
         rp.eMpa = GroupStep::kSteelEstar * (1.0 - 0.3 * 0.3) * 1e-6;
         rp.betaLayer = 0.0;
         rp.mu = friction;
         rp.twist = true;
         rp.tanDelta = 0.0;
         rp.normalLoss = false;
         rp.gammaKK = RubberContact::CalibrateKK(rp, b.R, b.m, 1.0, std::min(0.99, std::max(0.05, elasticity)));
      }
      std::vector<Seg> wires;
      struct Joint { V3 v; int surface, weight; };   // a point shared by weight + 1 touchable pieces of one side (hinge_piece_sum's negative terms)
      mutable std::vector<Joint> joints;
      mutable bool jointsBuilt = false;
      void BuildJoints() const
      {
         if (jointsBuilt)
            return;
         jointsBuilt = true;
         joints.clear();
         for (const Seg &w : wires)
            for (const V3 &v : { w.a, w.b })
            {
               int m = 0;
               for (const Seg &o : wires)
                  if (o.surface == w.surface && ((o.a - v).Length() < 1e-9 || (o.b - v).Length() < 1e-9))
                     m++;
               bool known = false;
               for (const Joint &j : joints)
                  known = known || (j.surface == w.surface && (j.v - v).Length() < 1e-9);
               if (m >= 2 && !known)
                  joints.push_back({ v, w.surface, m - 1 });
            }
      }
      std::vector<RubberContact::SurfaceTri> tris;
      std::vector<int> cand;
      struct Quad { V3 p[4]; std::vector<int> cand; };   // an AddPlate face: its corners and its two triangles (hinge_piece_sum's face terms)
      struct Line { V3 a, b; double r; };               // a line two touchable elements share (hinge_piece_sum's overlap terms)
      std::vector<Quad> quads;
      std::vector<Line> lines;
      std::string overlapWarn;                          // FinishShape: a wire end on a face (no overlap term; the owner warns with hinge_piece_sum on)
      double rho = 0.0, halfT = 0.0;
      int surfaces = 1, surfaceNow = 0;   // a wire shape's separate contacts: split where a ball can touch two wires at once (the V of a wire W): nearest-only has a force step there, no root

   private:
      V3 Down(const double th) const { return Z * -std::cos(th) + N * std::sin(th); }
      V3 Front(const double th) const { return N * std::cos(th) + Z * std::sin(th); }
      // gravity on the centre of mass, V = -m g . r(theta) = -m (A cos theta + B sin theta) + const
      double GravA() const { const V3 cm = mp.Cm(); return -cm.y * g.Dot(Z) + cm.z * g.Dot(N); }
      double GravB() const { const V3 cm = mp.Cm(); return cm.y * g.Dot(N) + cm.z * g.Dot(Z); }
      double GravV(const double th) const { return -mp.m * (GravA() * std::cos(th) + GravB() * std::sin(th)); }
      double GravVp(const double th) const { return -mp.m * (-GravA() * std::sin(th) + GravB() * std::cos(th)); }
      double GravVpp(const double th) const { return mp.m * (GravA() * std::cos(th) + GravB() * std::sin(th)); }
      bool svSwitch = false;
      double svWrapped = 0.0;
      size_t svSwLog = 0;
   };

   // One side of a part, a surface it carries (the offline bench's; the engine has its own, on the same probe)
   class Shape : public GroupStep::Shape
   {
   public:
      const Part *owner = nullptr;
      int surface = 0;
      void Probe(const V3 &c, const double R, Hit &h) const override { owner->ProbeSurface(surface, c, R, trial != nullptr ? trial[0] : 0.0, h); }
   };
}
