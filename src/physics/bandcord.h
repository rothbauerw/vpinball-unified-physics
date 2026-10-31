// license:GPLv3+
//
// A band (a rubber stretched round posts) and a sling (a band with a kicker plate behind it) as a massless stretched cord, stepped so that energy is kept by
// construction (docs/BandCordDesign.txt; ProjectAssessment-2026-10-27.txt item 0). Replaces the node chain of contactbandpart.h once proven.
//
// Laws:
//  - the cord: the shortest path round its posts, the balls pressing on it, the stop posts behind it and the plate's face. Stored energy EA / (2 L0) (L - L0)^2 of
//    its length L, one tension round the loop. The force on a ball (the torque on the plate) is the tension pulling round it: -V'(L) dL/dq.
//  - backstops: Hertz + Hunt-Crossley of the rubber where a ball squashes the cord against an own post, a stop post or the plate's face.
//  - friction, twist and rolling resistance at a ball's contact with the cord: Coulomb with stick on the mid-step slip; the cord's surface is still along its length.
//  - the plate: one body (plate and plunger through the rigid link), its return spring and rest stop as stored energy, the coil's torque as a drive, the coil stop
//    inelastic (the step is shortened to end on it).
// Step: the balls and the plate move together; every stored energy's force is its discrete gradient over the step (Gonzalez: exact energy differences over the move) and
// every loss acts on the mid-step rates, so once the step's equations are solved, energy end = start + coil work - losses at any step size. The step size is set by
// accuracy only (ChooseStep).
// Pure SI maths, no VPX types and no allocation. Its own solve (Step) is planar (the band's plane; spins are 3D): the offline test runs/2026-10-29-band-cord.
// In a group (groupstep.h; GroupStepDesign.txt 10) the group's 3D solve moves the band with the balls: a ball's height above the cord's centre line is in the cord's
// energy (round a ball the path keeps the clearance sqrt((R + r)^2 - dz^2); the squash is measured along the 3D line), so the band lifts a ball whose centre is
// above the cord and presses one below it into the playfield; the grip turns into the 3D contact plane. Offline test runs/2026-10-29-group-band.

#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <cstdio>
#include <cstring>

#include "physics/hitsling.h"
#include "physics/groupstep.h"

namespace BandCord
{
   constexpr double kPi = 3.14159265358979323846;
   constexpr int kMaxBalls = 8, kMaxPosts = 8, kMaxStops = 4, kMaxSwitches = 4;
   constexpr int kMaxPath = kMaxPosts + kMaxStops + 2 + 2 * kMaxBalls;
   constexpr int kMaxDof = 2 * kMaxBalls + 1;

   struct V2
   {
      double x = 0.0, y = 0.0;
      V2() = default;
      V2(double a, double b) : x(a), y(b) { }
      V2 operator+(const V2 &o) const { return V2(x + o.x, y + o.y); }
      V2 operator-(const V2 &o) const { return V2(x - o.x, y - o.y); }
      V2 operator*(double s) const { return V2(x * s, y * s); }
      double Dot(const V2 &o) const { return x * o.x + y * o.y; }
      double Cross(const V2 &o) const { return x * o.y - y * o.x; }
      double Length() const { return std::sqrt(x * x + y * y); }
   };
   inline V2 RightOf(const V2 &t) { return V2(t.y, -t.x); }   // the outside of a counter-clockwise loop

   struct V3
   {
      double x = 0.0, y = 0.0, z = 0.0;
      V3() = default;
      V3(double a, double b, double c) : x(a), y(b), z(c) { }
      V3 operator+(const V3 &o) const { return V3(x + o.x, y + o.y, z + o.z); }
      V3 operator-(const V3 &o) const { return V3(x - o.x, y - o.y, z - o.z); }
      V3 operator*(double s) const { return V3(x * s, y * s, z * s); }
      double Dot(const V3 &o) const { return x * o.x + y * o.y + z * o.z; }
      V3 Cross(const V3 &o) const { return V3(y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x); }
      double Length() const { return std::sqrt(x * x + y * y + z * z); }
   };

   struct Settings
   {
      double cordRadius = 2.69875e-3;           // m
      double tensionN = 20.0;                   // the cord's tension at rest, N
      double prestretch = 0.07;                 // its strain at rest
      double tanDelta = 0.15;                   // loss tangent of the stretch
      double refHz = 20.0;                      // the frequency tanDelta is given at
      double eMpa = 2.0;                        // rubber modulus of the backstops (the cord squashed between a ball and a post, a stop or the face)
      double alphaSPerM = 0.18;                 // Hunt-Crossley damping of the backstops, s/m (the post-rubber value, user 2026-10-29)
      double betaLayer = 0.33;                  // thin-layer coefficient of the squash (the shared rubber law's beta_layer; unchecked)
      double layerM = 0.0;                      // rubber thickness of the squash, m; 0 = the cord's diameter (as layer_mm 0)
      double mu = 0.8;                          // ball on cord, and cord on the surfaces behind it (black / white rubber; user, 2026-10-30)
      double playfieldWorldZ = 0.0;             // the playfield's height in the group's world, under the cord
      RubberContact::RubberParams ringRp;       // a rubber on a post (GroupStepDesign 13): the ring's own values when ringRpSet, else made from the values above
      bool ringRpSet = false;
      double cordTwistScale = 1.0;              // on the cord's twist stiffness G J (1 / L1 + 1 / L2), G = E / 3, J = pi r^4 / 2, L1 and L2 along the cord to its next post or ball each side (posts ASSUMED to hold the twist)
      double postMu = 0.8;                      // ASSUMED, the ball's: the cord on an own post (capstan friction, GroupStepDesign 11.24)
      double twistRadiusM = 1.0e-3;             // ASSUMED: the contact's radius for the twist torque (limit mu N a)
      double rollingResistance = 0.005;         // ASSUMED: rolling resistance coefficient on the cord (limit rr R N)
      double ballMassKg = 0.080;
      double ballRadiusM = 25.0 * 0.53975e-3;
      // the step (ChooseStep): accuracy rules, first guesses to be set from the offline sweep
      double omegaH = 0.6;                      // h <= omegaH / omega of the stiffest body (0.6 from the 2026-10-29 sweep: sling hit 26 us per contact ms, exit within 1.7 %)
      double turnMaxRad = 2.0 * kPi / 180.0;    // the cord's bend at a ball, and the plate, turn at most this per step
      double hMax = 1e-3;
      double hMin = 1e-9;                       // a step halved below this fails
      // the solver
      double tolV = 1e-12;                      // solved when a further round changes no ball's speed (or its friction's share of it) by more than this, m/s
      double tolW = 1e-9;                       // the plate's, rad/s
      int maxRounds = 60;
      SlingSettings sling;                      // the plate, its coil and its driver (the plate's link is taken as rigid: sling slop 0)
      double switchGapM = 1.5875e-3;            // a blade closes at this displacement of the cord toward it
   };

   struct Circle { V2 c; double r = 0.0; };    // centre, centreline clearance (an own post: post radius + cord radius)

   // The plate's face (as contactbandpart.h PlateGeom): at plate angle phi the face at band height is at originN + a / cos(t) + h tan(t) along n from the span's
   // rest centreline AB, t = phi + rest, over f0..f1 of AB
   struct PlateGeom
   {
      V2 A, B;              // the span's rest centreline
      V2 n;                 // unit normal from the plate toward the band and the ball
      double originN = 0.0, a = 0.0, h = 0.0, rest = 0.0;
      double f0 = 0.25, f1 = 0.75;
   };
   struct SwitchSpec { V2 P, b; bool coil = false; };   // a blade: rest point on the centreline, unit vector toward the blade; coil = wired to the driver
   struct WallSpec { V2 A, B; double halfThick = 0.0; };   // a wall or rail behind the band (inside the loop): its centre line and half its thickness
   constexpr int kMaxWalls = 4;

   struct Geometry
   {
      int nPosts = 0;
      Circle posts[kMaxPosts];   // counter-clockwise
      int nStops = 0;
      Circle stops[kMaxStops];   // posts behind the band (inside the loop)
      bool hasPlate = false;
      PlateGeom plate;
      int nSwitches = 0;
      SwitchSpec sw[kMaxSwitches];
      int nWalls = 0;            // walls and rails behind the band (group step only; GroupStepDesign 11.2)
      WallSpec walls[kMaxWalls];
      bool postTwin[kMaxPosts] = {};   // a post the cord passes twice (a post inside the loop touching its front and rear runs: listed once per run): the second listing has no squash of its own
   };

   // The energy books and counters (losses positive)
   struct Books
   {
      double postSlipJ = 0.0;   // (in frictionJ) the cord slipping round its posts
      double floorMaxM = 0.0;   // the deepest squash of the cord on the playfield. Any at all: the rubber is placed wrong (user, GroupStepDesign 11.24): shown, not tuned
      bool cordOnPlayfield = false;
      double coilJ = 0.0, cordLossJ = 0.0, backstopLossJ = 0.0, frictionJ = 0.0, twistJ = 0.0, rollJ = 0.0, plateLossJ = 0.0, stopJ = 0.0;
      long steps = 0, rounds = 0, halvings = 0;
      int maxRounds = 0, switchCloses = 0;
      double maxStretch = 0.0, maxOverlapM = 0.0, hMinUsed = 1.0;
      bool backTouched = false, through = false, failed = false;
      double Losses() const { return cordLossJ + backstopLossJ + frictionJ + twistJ + rollJ + plateLossJ + stopJ; }
   };

   enum Kind : int { kPost = 0, kStop, kFace, kBall, kBallBack, kWallEnd };   // kWallEnd id = 2 x wall + end (group step only)
   struct Elem { int kind = 0, id = 0; double wrapRef = 0.0; };
   struct Path { int n = 0; Elem e[kMaxPath]; };
   // z: a ball centre's height above the cord's centre line; held / mh: the cord's arc on that ball held in this plane (group step: GroupStepDesign 11.3)
   struct Config { V2 c[kMaxBalls]; double z[kMaxBalls] = {}; double phi = 0.0; bool held[kMaxBalls] = {}; V3 mh[kMaxBalls]; };
   struct PathGeom
   {
      V2 cen[kMaxPath];
      double rho[kMaxPath];   // signed clearance: + kept inside the loop (left of travel), - outside
      V2 t[kMaxPath];         // the run from pulley k to k + 1
      double l[kMaxPath];
      double wrap[kMaxPath];
      double L = 0.0;
   };

   class Part : public GroupStep::Part
   {
   public:
      Settings set;
      Geometry geo;
      Books books;
      int nb = 0;
      V2 c[kMaxBalls], v[kMaxBalls];
      double z[kMaxBalls] = {};                   // a ball centre's height above the cord's centre line (0 in the planar solve)
      V3 w[kMaxBalls];
      V2 g;                                       // gravity in the band's plane, m/s^2
      double phi = 0.0, omega = 0.0, iHigh = 0.0; // the plate: angle from rest toward the coil stop, speed, coil current as a fraction of V/R
      bool pinned = false;                        // the plate on its coil stop
      SlingDriver driver;
      bool switchClosed[kMaxSwitches] = {};
      Path path;                                  // the cord's path now
      // Solver diagnostics for the offline trace (off by default; when on, the last Solve call's rounds are kept and Step counts its tries)
      struct SolveDiag
      {
         static constexpr int kRounds = 64;
         struct Round
         {
            double moveV, moveW;        // the round's Newton move: balls m/s, plate rad/s
            double slip, twist, roll;   // the round's change of the friction / twist / rolling impulses, m/s
            int stick0, twistHeld0, rollHeld0;   // ball 0: friction sticking (1) or sliding (0), twist / rolling held within their limit (1) or at it (0); -1 = no contact
            double scale;               // how much of the Newton move was taken (1 = all)
         };
         bool on = false;
         int tries = 0;                 // Solve calls in the last Step (1 = taken first time; more = halved or the coil stop's search)
         int nRound = 0;                // rounds of the last Solve call
         int stepRule = 0;              // the rule that set the last ChooseStep (always kept): 0 hMax, 1 ball turn, 2 ball stiffness, 3 backstop arrival, 4 plate stiffness, 5 plate turn, 6 floor
         double kCord0 = 0.0, kSquash0 = 0.0;   // ball 0's stiffness in the last ChooseStep, N/m: the cord (stretch + bend) and the backstop squash parts
         Round r[kRounds];
      };
      mutable SolveDiag diag;
      PathGeom G;                                 // its geometry

      // After set and geo (and before balls): derived constants and the cord at rest. The rest path's spans are then readable with SpanFrame.
      void Configure()
      {
         m = set.ballMassKg;
         R = set.ballRadiusM;
         r = set.cordRadius;
         Ib = 0.4 * m * R * R;
         EA = set.tensionN / set.prestretch;
         const double eStar = set.eMpa * 1e6 / 0.75;
         const double reff = R * r / (R + r);   // the ball on the cord's round cross-section
         kH = 4.0 / 3.0 * eStar * std::sqrt(reff);
         hL = set.layerM > 0.0 ? set.layerM : 2.0 * r;
         cL = kH * set.betaLayer * reff * std::sqrt(reff);
         nb = 0;
         phi = omega = iHigh = 0.0;
         pinned = false;
         path.n = geo.nPosts;
         for (int k = 0; k < geo.nPosts; k++)
            path.e[k] = Elem { kPost, k, 0.0 };
         Config q;
         Build(path, q, G);
         for (int k = 0; k < path.n; k++)   // a counter-clockwise loop turns left at every post: the turn in (0, 2 pi)
         {
            double wr = G.wrap[k];
            while (wr <= 0.0)
               wr += 2.0 * kPi;
            path.e[k].wrapRef = wr;
         }
         Build(path, q, G);
         Lr = G.L / (1.0 + set.prestretch);
         cv = EA * set.tanDelta / (2.0 * kPi * std::max(set.refHz, 1e-3)) / Lr;
         if (geo.hasPlate)
         {
            FlipperModel::Derive(set.sling.mech, der);
            driver = SlingDriver();
            driver.SetPulse(set.sling.pulseMs * 1e-3);
            driver.SetDelay(set.sling.coilDelayMs * 1e-3);
            phiRest = RestPhi();
         }
         for (int s = 0; s < kMaxSwitches; s++)
            switchClosed[s] = false;
         books = Books();
         Lcur = G.L;
      }
      // span k (from own post k to post k + 1) at rest: its centreline ends and its outward unit normal
      void SpanFrame(int k, V2 &A, V2 &B, V2 &nOut) const
      {
         const int k1 = (k + 1) % path.n;
         nOut = RightOf(G.t[k]);
         A = G.cen[k] + nOut * G.rho[k];
         B = G.cen[k1] + nOut * G.rho[k1];
      }
      int AddBall(const V2 &pos, const V2 &vel)
      {
         if (nb >= kMaxBalls)
            return -1;
         c[nb] = pos;
         v[nb] = vel;
         w[nb] = V3();
         Pf[nb] = Mt[nb] = 0.0;
         Pz[nb] = 0.0;
         Mr[nb] = V3();
         nb++;
         Config q = Now();
         Settle(path, q, G);   // a ball placed inside the band (a contact found late) is wrapped at once
         Commit(path);
         return nb - 1;
      }

      double Energy() const
      {
         double e = CordV(Lcur) + PlateV(phi) + (geo.hasPlate ? 0.5 * Inertia(phi) * omega * omega : 0.0);
         const Config q = Now();
         for (int i = 0; i < nb; i++)
         {
            e += 0.5 * m * v[i].Dot(v[i]) + 0.5 * Ib * w[i].Dot(w[i]) - m * g.Dot(c[i]);
            for (int j = 0; j < NumBackstops(); j++)
               e += SquashV(Overlap(i, j, q, Dofs(), nullptr));
         }
         return e;
      }
      double CordEnergy() const { return groupMode ? SpanV(spanL) : CordV(Lcur); }
      double MaxOverlapNow() const
      {
         double mx = 0.0;
         const Config q = Now();
         for (int i = 0; i < nb; i++)
            for (int j = 0; j < NumBackstops(); j++)
               mx = std::max(mx, Overlap(i, j, q, Dofs(), nullptr));
         return mx;
      }
      double Length() const { return Lcur; }
      double Stretch() const
      {
         if (!groupMode)
            return Lcur / Lr - 1.0;
         double st = -1.0;
         for (int j = 0; j < geo.nPosts; j++)
            st = std::max(st, spanL[j] / spanS[j] - 1.0);
         return st;
      }
      double Tension() const   // (group mode: the tightest span's)
      {
         if (!groupMode)
            return Lcur > Lr ? EA * (Lcur - Lr) / Lr : 0.0;
         double tm = 0.0;
         for (int j = 0; j < geo.nPosts; j++)
            tm = std::max(tm, SpanT(j, spanL[j]));
         return tm;
      }
      double PlateDeg() const { return phi * 180.0 / kPi; }
      double CurrentA() const { return geo.hasPlate ? iHigh * set.sling.mech.vSupply / set.sling.mech.rP : 0.0; }
      bool InPath(int kind, int id) const { return groupMode ? Find3(path3, kind, id) >= 0 : Find(path, kind, id) >= 0; }
      bool Touching(int i) const
      {
         if (InPath(kBall, i))
            return true;
         const Config q = Now();
         for (int j = 0; j < NumBackstops(); j++)
            if (Overlap(i, j, q, Dofs(), nullptr) > 0.0)
               return true;
         if (groupMode)
            for (int j = 0; j < geo.nPosts + geo.nStops; j++)
               if (OverlapRaw(i, j, q, Dofs(), nullptr) - postState[i].s[j].zero > 0.0)
                  return true;
         return false;
      }

      // Once per 1 ms cycle: the switch blades, and the coil driver from its switches
      void CycleBegin()
      {
         bool coil = false;
         for (int s = 0; s < geo.nSwitches; s++)
         {
            const bool closed = (groupMode ? SwitchDisp3(s) : SwitchDisp(s)) >= set.switchGapM;
            if (closed && !switchClosed[s])
               books.switchCloses++;
            switchClosed[s] = closed;
            coil = coil || (closed && geo.sw[s].coil);
         }
         if (geo.hasPlate)
            driver.Update(coil);
      }
      // the cord's displacement toward blade s from its rest point: where the path first crosses the blade's line coming from outside
      double SwitchDisp(int s) const
      {
         const SwitchSpec &sw = geo.sw[s];
         double best = 1e30;
         const int n = path.n;
         for (int k = 0; k < n; k++)
         {
            const int k1 = (k + 1) % n;
            // the run k, then the arc on pulley k1 (in chords of up to 10 deg)
            V2 pts[40];
            int np = 0;
            pts[np++] = G.cen[k] + RightOf(G.t[k]) * G.rho[k];
            pts[np++] = G.cen[k1] + RightOf(G.t[k]) * G.rho[k1];
            const int na = std::min(36, std::max(1, (int)std::ceil(std::fabs(G.wrap[k1]) / (10.0 * kPi / 180.0))));
            const double a0 = std::atan2(G.t[k].y, G.t[k].x);
            const double sg = G.rho[k1] >= 0.0 ? 1.0 : -1.0;
            for (int a = 1; a <= na && np < 40; a++)
            {
               const double th = a0 + sg * G.wrap[k1] * a / na;
               pts[np++] = G.cen[k1] + RightOf(V2(std::cos(th), std::sin(th))) * G.rho[k1];
            }
            for (int p = 0; p + 1 < np; p++)
            {
               // X = sw.P + s b crossing the chord pts[p] -> pts[p + 1]
               const V2 e = pts[p + 1] - pts[p];
               const double den = sw.b.Cross(e);
               if (std::fabs(den) < 1e-15)
                  continue;
               const V2 d = pts[p] - sw.P;
               const double sPar = d.Cross(e) / den;
               const double u = d.Cross(sw.b) / den;
               if (u >= 0.0 && u <= 1.0)
                  best = std::min(best, sPar);
            }
         }
         return best < 1e29 ? best : 0.0;
      }

      // The engine's step choice from the state now (accuracy only)
      double ChooseStep() const
      {
         double h = set.hMax;
         int rule = 0;
         const auto lim = [&](const double c, const int r) { if (c < h) { h = c; rule = r; } };
         const Config q = Now();
         const double T = Tension(), kAx = EA / Lr;
         const double lMin = 1e-3;
         diag.kCord0 = diag.kSquash0 = 0.0;
         for (int i = 0; i < nb; i++)
         {
            double k = 0.0;
            const int e = Find(path, kBall, i);
            if (e >= 0)
            {
               const int n = path.n;
               const double lin = std::max(G.l[(e + n - 1) % n], lMin), lout = std::max(G.l[e], lMin);
               const V2 gk = G.t[(e + n - 1) % n] - G.t[e];
               k += kAx * gk.Dot(gk) + T * (1.0 / lin + 1.0 / lout);
               if (i == 0)
                  diag.kCord0 = k;
               const double sp = v[i].Length();
               if (sp > 1e-9)
                  lim(set.turnMaxRad * std::min(lin, lout) / sp, 1);
            }
            for (int j = 0; j < NumBackstops(); j++)
            {
               double gr[kMaxDof];
               const Dofs d = MakeDofs();
               const double del = Overlap(i, j, q, d, gr);
               const double dDot = gr[2 * i] * v[i].x + gr[2 * i + 1] * v[i].y + (d.phi >= 0 ? gr[d.phi] * omega : 0.0);   // closing speed, m/s
               // the depth this approach reaches (all the ball's kinetic energy into the squash): the stiffness while it lasts
               const double vc = std::max(std::fabs(dDot), 0.05);
               if (del > 0.0 || -del < vc * set.hMax)
               {
                  double lo = 0.0, hi = hL;
                  for (int it = 0; it < 40; it++)
                  {
                     const double mid = 0.5 * (lo + hi);
                     (SquashV(mid) < 0.5 * m * vc * vc ? lo : hi) = mid;
                  }
                  const double kExp = SquashK(std::max(del, lo));
                  k += kExp;
                  if (i == 0)
                     diag.kSquash0 += kExp;
                  if (del < 0.0 && dDot > 0.0)   // reached in about two steps, not overshot
                     lim(std::max(0.5 * -del / dDot, set.omegaH / std::sqrt(kExp / m)), 3);
               }
            }
            if (k > 0.0)
               lim(set.omegaH / std::sqrt(k / m), 2);
         }
         if (geo.hasPlate && !pinned)
         {
            const double I = Inertia(phi);
            double k = std::max(PlateVpp(phi), 0.0);
            for (int e = 0; e < path.n; e++)
               if (path.e[e].kind == kFace)
               {
                  const int n = path.n;
                  const double dd = FaceDD(phi);
                  k += T * dd * dd * (1.0 / std::max(G.l[(e + n - 1) % n], lMin) + 1.0 / std::max(G.l[e], lMin));
               }
            double dLdphi = 0.0;
            for (int e = 0; e < path.n; e++)
               if (path.e[e].kind == kFace)
                  dLdphi += (G.t[(e + path.n - 1) % path.n] - G.t[e]).Dot(geo.plate.n) * FaceDD(phi);
            k += kAx * dLdphi * dLdphi;
            if (k > 0.0)
               lim(set.omegaH / std::sqrt(k / I), 4);
            if (std::fabs(omega) > 1e-9)
               lim(set.turnMaxRad / std::fabs(omega), 5);
         }
         diag.stepRule = h < 1e-7 ? 6 : rule;
         return std::max(h, 1e-7);
      }

      // One step of at most hTry: shorter if the solver needs it (halving) or the plate reaches its coil stop. Returns the time taken, 0 = failed.
      double Step(double hTry)
      {
         double h = hTry;
         Trial T;
         diag.tries = 1;
         while (!Solve(h, T))
         {
            diag.tries++;
            h *= 0.5;
            books.halvings++;
            if (h < set.hMin)
            {
               books.failed = true;
               return 0.0;
            }
         }
         if (geo.hasPlate && !pinned && T.q1.phi > PhiStop())
         {
            // the step that ends on the coil stop (Illinois false position on the step length)
            double lo = 0.0, flo = phi - PhiStop(), hi = h, fhi = T.q1.phi - PhiStop();
            Trial Ts, Tbelow;
            double hBelow = 0.0;
            int side = 0;
            for (int it = 0; it < 80; it++)
            {
               double hm = (lo * fhi - hi * flo) / (fhi - flo);
               if (!(hm > lo && hm < hi))
                  hm = 0.5 * (lo + hi);
               diag.tries++;
               if (!Solve(hm, Ts))
               {
                  books.failed = true;
                  return 0.0;
               }
               const double f = Ts.q1.phi - PhiStop();
               if (f <= 0.0)
               {
                  Tbelow = Ts;
                  hBelow = hm;
                  if (f > -1e-10)
                     break;
                  lo = hm;
                  flo = f;
                  if (side == -1)
                     fhi *= 0.5;
                  side = -1;
               }
               else
               {
                  hi = hm;
                  fhi = f;
                  if (side == 1)
                     flo *= 0.5;
                  side = 1;
               }
            }
            if (hBelow <= 0.0)
            {
               books.failed = true;
               return 0.0;
            }
            T = Tbelow;
            h = hBelow;
            // the plunger stops on the steel coil stop: inelastic, the plate and plunger lose their speed
            T.stopJ = 0.5 * Inertia(T.q1.phi) * T.omega1 * T.omega1;
            T.omega1 = 0.0;
            T.pinNow = true;
         }
         Accept(T, h);
         return h;
      }

      // The group's balls changed (the engine's islands, GroupStepDesign 8): what this part keeps per ball follows them. oldOfNew[i]: the ball's index before, or -1 for a new
      // ball. A ball that has gone leaves the path.
      void RemapBalls(const int *oldOfNew, int nNew)
      {
         const auto Move = [&](auto *arr, auto zero) {
            std::decay_t<decltype(arr[0])> t[kMaxBalls];
            for (int i = 0; i < nNew; i++)
               t[i] = oldOfNew[i] >= 0 ? arr[oldOfNew[i]] : zero;
            for (int i = 0; i < kMaxBalls; i++)
               arr[i] = i < nNew ? t[i] : zero;
         };
         Move(c, V2()); Move(v, V2()); Move(z, 0.0); Move(w, V3());
         Move(gP, V3()); Move(gTq, V3()); Move(gPrate, V3()); Move(gTqRate, V3()); Move(gFn, V3()); Move(gFnTaken, V3());
         Move(holdOn, false); Move(mHold, V3(0.0, 0.0, 1.0)); Move(psi, 0.0); Move(kTw, 0.0); Move(cTw, 0.0); Move(gSurfRate, 0.0);
         Move(gRootOk, false);
         Move(bsZero, BsZero()); Move(bsNew, true); Move(postState, PostStates());
         int newOfOld[kMaxBalls];
         for (int i = 0; i < kMaxBalls; i++)
            newOfOld[i] = -1;
         for (int i = 0; i < nNew; i++)
            if (oldOfNew[i] >= 0)
               newOfOld[oldOfNew[i]] = i;
         int n = 0;
         for (int k = 0; k < path3.n; k++)
         {
            if (path3.e[k].kind == kBall || path3.e[k].kind == kBallBack)
            {
               if (newOfOld[path3.e[k].id] < 0)
                  continue;
               path3.e[k].id = newOfOld[path3.e[k].id];
            }
            path3.e[n] = path3.e[k];
            path3.ct[n] = path3.ct[k];
            n++;
         }
         path3.n = n;
         nb = nNew;
      }
      // an energy step (GroupStep::Group::EnergyStep)
      void SaveState() override
      {
         if (!snap)
            snap = std::make_shared<Part>(*this);
         else
            *snap = *this;
         snap->snap.reset();
      }
      void RestoreState() override
      {
         const std::shared_ptr<Part> keep = snap;
         *this = *keep;
         snap = keep;
      }
      void HoldArcs() override   // the held arcs stay where they are, as one held still by a surface (their turn is a solve of its own, and it may have no root either), and no pulley joins in this step (which of them are on it can go round in a ring near a pressed pair: the next step's solve sees to that)
      {
         for (int i = 0; i < nb; i++)
            gStill[i] = holdOn[i];
         gSamePath = true;
      }
      void ExplicitMove(const double *Res, const int off, const double h, double *x) const override
      {
         if (gPl)
            x[off] -= Res[off] * h / (2.0 * gI0);
      }
      void BallsMoved(const GroupStep::Ball *b_, const int n_) override   // the last resort of an energy step: the path follows the balls with the elements it has
      {
         gSamePath = true;
         Sync(b_, n_);
         gSamePath = false;
      }
      double KineticJ() const override { return geo.hasPlate ? 0.5 * Inertia(phi) * omega * omega : 0.0; }
      void ScaleSpeeds(const double s) override { omega *= s; }
      void AfterEnergyStep() override   // the push and the grip of a step with no solve are nothing to go on from
      {
         for (int i = 0; i < kMaxBalls; i++)
         {
            gFnTaken[i] = snap->gFnTaken[i];
            gPrate[i] = snap->gPrate[i];
            gTqRate[i] = snap->gTqRate[i];
            gSurfRate[i] = snap->gSurfRate[i];
         }
         gAlphaLast = snap->gAlphaLast;
         books = snap->books;
         gdiag.settleFails = snap->gdiag.settleFails;
      }
      double SwitchPush(int s) const { return SwitchDisp3(s); }
      double TwistOf(int i) const { return psi[i]; }
      bool SwitchIsClosed(int s) const { return switchClosed[s]; }
      bool PlateQuiet() const { return !geo.hasPlate || (std::fabs(phi - phiRest) < 1e-5 && std::fabs(omega) < 1e-3 && !driver.Coil() && driver.pending <= 0.0); }
      // where the plate comes to rest: its return spring's preload balanced by the rest stop (phi below 0, PlateV's slope 0)
      double RestPhi() const
      {
         if (!geo.hasPlate)
            return 0.0;
         double lo = -0.2, hi = 0.0;
         for (int it = 0; it < 40; it++)
         {
            const double mid = 0.5 * (lo + hi);
            (PlateVp(mid) > 0.0 ? hi : lo) = mid;
         }
         return 0.5 * (lo + hi);
      }
      // how far ball i is clear of the cord (the group's release rule): 0 while it is on the cord, less than 0 while it squashes the cord on a post, a stop or the plate's
      // face, whose stored energy is the part's: a ball let go there would lose it. A backstop it was found inside and has not pressed since holds nothing and does not count.
      // how fast ball i (velocity vw, world) closes on the cord's push on it, from the last step: the speed of a hit
      double ClosingOnCord(const int i, const GroupStep::V3 &vw) const
      {
         const V3 f = gdiag.cordF[i], u = ToBand(vw);
         const double l = f.Length();
         return l > 0.0 ? std::max(-u.Dot(f) / l, 0.0) : u.Length();
      }
      double CordGap(int i) const
      {
         if (Find3(path3, kBall, i) >= 0)
            return 0.0;
         const V3 C(c[i].x, c[i].y, z[i]);
         double gap = 1e30;
         const Config q = Now();
         for (int j = 0; j < NumBackstops(); j++)
         {
            const double del = Overlap(i, j, q, Dofs(), nullptr);
            if (del > 0.0 || bsZero[i].d[j] <= 0.0)
               gap = std::min(gap, -del);
         }
         if (groupMode)   // its contact on a post, a group contact: as a backstop above
            for (int j = 0; j < geo.nPosts + geo.nStops; j++)
            {
               const double zero = postState[i].s[j].zero, del = OverlapRaw(i, j, q, Dofs(), nullptr) - zero;
               if (del > 0.0 || zero <= 0.0)
                  gap = std::min(gap, -del);
            }
         for (int k = 0; k < path3.n; k++)
         {
            const V3 dc = C - path3.ct[k].X;
            const double u = std::min(std::max(dc.Dot(G3.t[k]), 0.0), G3.l[k]);
            gap = std::min(gap, (dc - G3.t[k] * u).Length() - (R + r));
         }
         return gap;
      }

      // ---- A ball pressing the cord on an own post or a stop post is a rubber on a post (GroupStepDesign 13, ModelStandards 2d): an ordinary group contact with the ring's law,
      // shear spring, twist, rolling and step rule, not this part's. The depth is the ring's: R + r less the distance to the cord's centre line at the post. The cord's own force and
      // grip go on acting: its grip's limit is mu x the cord's push, the contact's is mu x its own, so one Coulomb limit on the whole push. No blend: the force moves over as the squash grows.
      class PostShape : public GroupStep::Shape
      {
      public:
         const Part *part = nullptr;
         int j = 0;
         void Probe(const GroupStep::V3 &cw, double Rb, GroupStep::Hit &h) const override { part->ProbePost(j, cw, Rb, h); }
      };
      void ProbePost(const int j, const GroupStep::V3 &cw, const double Rb, GroupStep::Hit &h) const
      {
         const bool own = j < geo.nPosts;
         if (own && geo.postTwin[j])
         {
            h.delta = -1.0;
            h.n = ToWorld(V3(0.0, 0.0, 1.0));
            h.reff = Rb;
            h.curvDist = 0.0;
            return;
         }
         const V3 p = ToBand(cw - GroupStep::V3(frame.O.x, frame.O.y, frame.O.z));
         const Circle &cs = own ? geo.posts[j] : geo.stops[j - geo.nPosts];
         const V2 dv = V2(p.x, p.y) - cs.c;
         const double D = dv.Length(), aa = D - cs.r, zz = p.z, dist = std::sqrt(aa * aa + zz * zz);   // the groove at the cord's height, a stop's too: the rubber pinned on it does not slide (user, 2026-10-30)
         const V2 eh = D > 0.0 ? dv * (1.0 / D) : V2(1.0, 0.0);
         const double na = dist > 0.0 ? aa / dist : 1.0, nz = dist > 0.0 ? zz / dist : 0.0;
         h.delta = Rb + r - dist;
         h.n = ToWorld(V3(eh.x * na, eh.y * na, nz));
         const double R1 = 1.0 / (1.0 / Rb + 1.0 / r), kLong = RubberContact::TubeAlongCurvature(1.0 / cs.r, r, 1.0);   // the effective radius as ProbeTube
         h.reff = std::sqrt(R1 / std::max(1.0 / Rb + kLong, 0.05 / Rb));
         h.curvDist = D;
      }
      // Into the group, before its Configure
      void AddPostShapes(GroupStep::Group &g)
      {
         RubberContact::RubberParams rp = set.ringRp;
         if (!set.ringRpSet)
         {
            rp.eMpa = set.eMpa;
            rp.alphaSPerM = set.alphaSPerM;
            rp.betaLayer = set.betaLayer;
            rp.mu = set.mu;
            rp.tanDelta = set.tanDelta;
         }
         postFirst = g.ns;
         nPostShapes = 0;
         for (int j = 0; j < geo.nPosts + geo.nStops && g.ns < GroupStep::kMaxShapes; j++)
         {
            PostShape &s = postShape[j];
            s.part = this;
            s.j = j;
            s.rp = rp;
            s.layerM = hL;
            g.shape[g.ns++] = &s;
            nPostShapes++;
         }
      }
      // After its Configure, with the part attached: what each ball kept on each post
      void RestorePostItems(GroupStep::Group &g) const
      {
         for (int k = 0; k < g.ni; k++)
         {
            const GroupStep::Item &it = g.item[k];
            if (it.b >= 0 || it.s < postFirst || it.s >= postFirst + nPostShapes)
               continue;
            const PostState &st = postState[it.a].s[it.s - postFirst];
            g.SetShapeItem(k, st.touching, st.zero, st.xi);
         }
      }
      void AfterStep(const GroupStep::Group &g) override
      {
         for (int k = 0; k < g.ni; k++)
         {
            const GroupStep::Item &it = g.item[k];
            if (it.b >= 0 || it.s < postFirst || it.s >= postFirst + nPostShapes)
               continue;
            PostState &st = postState[it.a].s[it.s - postFirst];
            st.touching = it.touching;
            st.xi = it.xi;
            st.zero = it.dZero;
         }
      }

      // ---- the band in a group's solve (groupstep.h; GroupStepDesign.txt 10, 11) ----
      // The group moves the band with its balls: the balls are the group's (every ball of the group: the band decides which touch it), read in the band's frame
      // (its plane at the cord's height, ez up); the plate's angle is the band's one unknown. The cord is free in height (GroupStepDesign 11): its path is the
      // shortest path in 3D round its own posts (they hold its height), the balls (spheres) and the round objects behind it (it slides up and down on them),
      // solved to rounding at every evaluation. The balls' mass and radius are the band's settings (ballMassKg, ballRadiusM). The band's planar gravity g is not
      // used (the group applies gravity).
      using GV3 = GroupStep::V3;
      struct Frame
      {
         V3 O;                                                                        // on the cord's centre line, at its height
         V3 ex = V3(1.0, 0.0, 0.0), ey = V3(0.0, 1.0, 0.0), ez = V3(0.0, 0.0, 1.0);   // the band's plane, and its normal (up from the playfield)
      };
      Frame frame;
      // per ball, from the last Eval and UpdateGrip of a group step (the trace; band frame): the cord's and the squash's force, the grip's force (impulse / h), the
      // band's push the grip's limits use, the contact-point slip along / across the cord at mid-step without the grip's share, whether the friction held
      struct GroupDiag
      {
         V3 cordF[kMaxBalls], squashF[kMaxBalls], gripF[kMaxBalls];
         double push[kMaxBalls] = {}, wrap[kMaxBalls] = {}, clear[kMaxBalls] = {}, slip1[kMaxBalls] = {}, slip2[kMaxBalls] = {};
         int stick[kMaxBalls] = {};   // 1 held, 0 sliding, -1 no push
         double heldTq[kMaxBalls] = {};   // the cord's pull on a held ball's turn, N m (0 if not held)
         // the grip's parts from the last UpdateGrip (averaged over the step): friction along the cord (P1) and across it (P2), N; rolling resistance and twist,
         // N m; the friction's limit mu N, N; the surface's friction on a turning pressed arc, N m; and the contact frame (band frame: t1 along, t2 across, n push)
         double gripP1[kMaxBalls] = {}, gripP2[kMaxBalls] = {}, gripRoll[kMaxBalls] = {}, gripTwist[kMaxBalls] = {}, gripLim[kMaxBalls] = {}, gripSurf[kMaxBalls] = {};
         V3 gripT1[kMaxBalls], gripT2[kMaxBalls], gripN[kMaxBalls];
         int held[kMaxBalls] = {};         // 1 held, 0 not
         long settleFails = 0;        // Evals whose path would not settle (the step's end or mid-step), the last one's place and pulley count
         bool settleFailMid = false;
         int settleFailN = 0;
      };
      // the last actions of a path's settling (kept for the last failure): '-' a pulley left (val: its wrap, rad), '+' one joined (val: how far the run crossed it, m),
      // at its place in the path, the path's pulley count before the action
      struct SettleLog
      {
         static constexpr int kN = 24;
         struct Entry { char act = ' '; int kind = 0, id = 0, at = 0, pathN = 0; double val = 0.0; };
         Entry e[kN];
         int n = 0;   // actions logged (the last kN kept, oldest at n % kN)
      };
      SettleLog settleLog;
      // the 3D path solves so far (diagnostics), and by caller: settling the end and mid paths (with trial insertions), the length's derivative, the held arcs'
      // turn, commits; evaluations
      struct PassStats
      {
         long solves = 0, passes = 0, fails = 0;
         int maxPasses = 0;
         long bySettle = 0, byGrad = 0, byBeta = 0, byOther = 0, evals = 0;
         void Add(const PassStats &a, const PassStats &b)   // the counts between two readings
         {
            solves += a.solves - b.solves;
            passes += a.passes - b.passes;
            fails += a.fails - b.fails;
            bySettle += a.bySettle - b.bySettle;
            byGrad += a.byGrad - b.byGrad;
            byBeta += a.byBeta - b.byBeta;
            byOther += a.byOther - b.byOther;
            evals += a.evals - b.evals;
         }
      };
      mutable int solveCaller = 0;   // 0 other, 1 settle, 2 derivative, 3 held arcs
      mutable PassStats passStats;
      bool debugSolve = false;   // print the last passes of a path solve that does not settle (diagnostics)
      bool debugGrad = false;    // check the length's exact derivative against central differences in every evaluation: the largest differences (m/m, m/rad)
      mutable double gradErrQ = 0.0, gradErrPhi = 0.0, gradErrBeta = 0.0;
      GroupDiag gdiag;
      int groupRule = 0;   // what set the last group step choice: 0 hMax, 1 ball turn, 2 stiffness (cord + squash), 3 squash arrival, 4 plate stiffness, 5 plate turn, 7 a steady squash at omegaHRest

      // after Configure and the frame: the group's balls as they lie (the path wrapped round them, no grip carried)
      void AttachGroup(const GroupStep::Ball *b, int n)
      {
         nb = 0;
         for (int i = 0; i < kMaxBalls; i++)
         {
            gP[i] = gTq[i] = gPrate[i] = gTqRate[i] = gFnTaken[i] = V3();
            postState[i] = PostStates();
         }
         // the cord at rest round its own posts (as Configure's planar rest path), at its height
         groupMode = true;
         gFloorTop = (V3(0.0, 0.0, set.playfieldWorldZ) - frame.O).Dot(frame.ez) + r;   // where the cord's centre meets the playfield
         {
            const Config q0;
            Path pr;
            pr.n = geo.nPosts;
            for (int k = 0; k < geo.nPosts; k++)
               pr.e[k] = Elem { kPost, k, 0.0 };
            PathGeom gr;
            Build(pr, q0, gr);
            path3.n = geo.nPosts;
            for (int k = 0; k < geo.nPosts; k++)
            {
               const int km = (k + geo.nPosts - 1) % geo.nPosts;
               double wr = gr.wrap[k];
               while (wr <= 0.0)   // a counter-clockwise loop turns left at every post
                  wr += 2.0 * kPi;
               path3.e[k] = Elem { kPost, k, wr };
               path3.ct[k] = Contact3();
               const V2 E2 = gr.cen[k] + RightOf(gr.t[km]) * gr.rho[k], X2 = gr.cen[k] + RightOf(gr.t[k]) * gr.rho[k];
               path3.ct[k].E = V3(E2.x, E2.y, 0.0);
               path3.ct[k].X = V3(X2.x, X2.y, 0.0);
               path3.ct[k].m = V3(0.0, 0.0, 1.0);
               path3.ct[k].tOut = V3(gr.t[k].x, gr.t[k].y, 0.0);
            }
            // each span's unstretched length: the same stretch in all at rest (the ring's tension equal round the posts)
            Geom3 g0;
            if (Solve3(path3, q0, g0))
            {
               SpanMap sm;
               MakeSpanMap(path3, g0, sm);
               SpanLengths(path3, g0, sm, spanL);
               for (int j = 0; j < geo.nPosts; j++)
                  spanS[j] = spanL[j] / (1.0 + set.prestretch);
            }
         }
         Sync(b, n, true);
         for (int i = 0; i < nb; i++)
         {
            const int e = Find3(path3, kBall, i);
            holdOn[i] = e >= 0;
            mHold[i] = e >= 0 ? path3.ct[e].m : V3(0.0, 0.0, 1.0);
         }
         Commit3(path3);
         for (int i = 0; i < nb; i++)
            TwistJoin(i);
      }

      // the cord's energy with ball i's held arc turned by ang about the cord's direction there, from the state now (the force check; 0 if not held)
      double EnergyTurned(int i, double ang) const
      {
         const int e = holdOn[i] ? Find3(path3, kBall, i) : -1;
         if (e < 0)
            return SpanV(spanL);
         const V3 ch = path3.ct[e].X - path3.ct[e].E, t1 = ch.Length() > 1e-9 ? ch * (1.0 / ch.Length()) : G3.t[e];
         Config q = Now();
         q.mh[i] = Turn(t1, ang, q.mh[i]);
         q.mh[i] = q.mh[i] * (1.0 / q.mh[i].Length());
         Path3 pp = path3;
         Geom3 gg;
         if (!Solve3(pp, q, gg))
            return 0.0;
         SpanMap sm;
         double ls[kMaxPosts];
         MakeSpanMap(pp, gg, sm);
         SpanLengths(pp, gg, sm, ls);
         return SpanV(ls);
      }
      // where the cord sits on ball i now (group mode): the elevation of its arc's middle above the ball's centre (deg; + above), the arc's wrap (deg), the
      // arc plane's tilt from level (deg); false if the cord does not wrap the ball
      bool CordOnBall(int i, double &elevDeg, double &wrapDeg, double &tiltDeg) const
      {
         const int e = groupMode ? Find3(path3, kBall, i) : -1;
         if (e < 0)
            return false;
         const V3 C(c[i].x, c[i].y, z[i]), mid = (path3.ct[e].E + path3.ct[e].X) * 0.5 - C;
         elevDeg = std::atan2(mid.z, std::sqrt(mid.x * mid.x + mid.y * mid.y)) * 180.0 / kPi;
         wrapDeg = G3.wrap[e] * 180.0 / kPi;
         tiltDeg = std::acos(std::min(1.0, std::fabs(path3.ct[e].m.z))) * 180.0 / kPi;
         return true;
      }
      // the force check: the push of the last Eval taken as a step's, so that the grip's limits know it
      void TakePushForCheck()
      {
         for (int i = 0; i < nb; i++)
            gFnTaken[i] = gFn[i];
      }
      // the force check: ball i's held arc turned by ang (a held grip under load), the path re-solved there
      void TurnHoldForCheck(int i, double ang)
      {
         const int e = holdOn[i] ? Find3(path3, kBall, i) : -1;
         if (e < 0)
            return;
         const V3 ch = path3.ct[e].X - path3.ct[e].E, t1 = ch.Length() > 1e-9 ? ch * (1.0 / ch.Length()) : G3.t[e];
         mHold[i] = Turn(t1, ang, mHold[i]);
         mHold[i] = mHold[i] * (1.0 / mHold[i].Length());
         Commit3(path3);
      }
      // the 3D path now, one line per pulley (diagnostics)
      void DumpPath3(FILE *f) const
      {
         for (int k = 0; k < path3.n; k++)
            fprintf(f, "   %c%d wrap %.6f deg  E (%.4f %.4f %.4f) mm  X (%.4f %.4f %.4f) mm  m (%.4f %.4f %.4f)  run %.4f mm\n", "PSFBKW"[path3.e[k].kind], path3.e[k].id,
                    G3.wrap[k] * 180.0 / kPi, path3.ct[k].E.x * 1e3, path3.ct[k].E.y * 1e3, path3.ct[k].E.z * 1e3, path3.ct[k].X.x * 1e3, path3.ct[k].X.y * 1e3,
                    path3.ct[k].X.z * 1e3, path3.ct[k].m.x, path3.ct[k].m.y, path3.ct[k].m.z, G3.l[k] * 1e3);
         fprintf(f, "   L %.6f mm, ball 0 at (%.4f %.4f %.4f) mm\n", Lcur * 1e3, c[0].x * 1e3, c[0].y * 1e3, z[0] * 1e3);
      }
      int NumDof() const override { return geo.hasPlate && !pinned ? 1 : 0; }

      double Energy(const GroupStep::Ball *b, int n) const override
      {
         double e = SpanV(spanL) + gFloorV + TwistV() + PlateV(phi) + (geo.hasPlate ? 0.5 * Inertia(phi) * omega * omega : 0.0);   // the path is committed at the balls' positions
         Config q = Now();
         for (int i = 0; i < std::min(n, nb); i++)
         {
            V2 vv;
            Project(b[i], q.c[i], q.z[i], vv);
         }
         for (int i = 0; i < nb; i++)
            if (!bsNew[i])   // (a ball just joined has stored nothing: Overlap)
               for (int j = 0; j < NumBackstops(); j++)
                  e += SquashV(Overlap(i, j, q, Dofs(), nullptr));
         return e;
      }

      // the step choice (GroupStepDesign 10.3): the planar solve's rules on the 3D state, the steady rule on the squash
      double ChooseStep(const GroupStep::Ball *b, int n, const GroupStep::Settings &s, int &rule) override
      {
         Sync(b, n);
         double h = s.hMax;
         int rl = 0;
         const auto lim = [&](const double cc, const int r0) { if (cc < h) { h = cc; rl = r0; } };
         const Config q = Now();
         const Dofs d = MakeDofs();
         double sMinS = spanS[0];
         for (int j = 1; j < geo.nPosts; j++)
            sMinS = std::min(sMinS, spanS[j]);
         const double T = Tension(), kAx = EA / sMinS;
         const double lMin = 1e-3;
         SpanMap sMapNow;
         MakeSpanMap(path3, G3, sMapNow);
         double delNow[kMaxBalls];
         FloorState(path3, G3, delNow);
         for (int i = 0; i < nb; i++)
         {
            const V3 vb = ToBand(b[i].v), wb = ToBand(b[i].w);
            double kStrict = 0.0, kRest = 0.0;
            const int e = Find3(path3, kBall, i);
            if (e >= 0)
            {
               const int pn = path3.n;
               const double lin = std::max(G3.l[(e + pn - 1) % pn], lMin), lout = std::max(G3.l[e], lMin);
               const V3 gk = G3.t[(e + pn - 1) % pn] - G3.t[e];
               const int spanE = sMapNow.arcSpan[e];
               const double T = SpanT(spanE, spanL[spanE]), kAx = EA / spanS[spanE];   // its own span
               kStrict += kAx * gk.Dot(gk) + T * (1.0 / lin + 1.0 / lout);   // stretch and bend (in every direction: the cord is free in height)
               const double sp = vb.Length();
               if (sp > 1e-9)
                  lim(s.turnMaxRad * std::min(lin, lout) / sp, 1);
               if (holdOn[i])   // the held arc: the ball's turn against the cord's bend round it (GroupStepDesign 11.4)
               {
                  // (the cord's bend and its twist in series, GroupStepDesign 11.28)
                  const double kl = T * (R + r) * (R + r) * (1.0 / lin + 1.0 / lout), kw = kTw[i] / (r * r), kb = kl * kw * R * R / (kl + kw * (R + r) * (R + r));
                  lim(s.omegaH / std::sqrt(std::max(kb, 1e-30) / Ib), 2);
               }
            }
            for (int j = 0; j < NumBackstops(); j++)
            {
               double gr[kMaxDof], gz = 0.0;
               const double del = Overlap(i, j, q, d, gr, &gz);
               const V3 g3(gr[2 * i], gr[2 * i + 1], gz);   // the overlap's gradient in the ball's centre
               const double dDot = g3.Dot(vb) + (d.phi >= 0 ? gr[d.phi] * omega : 0.0);
               const double vc = std::max(std::fabs(dDot), s.vFloor);
               if (!(del > 0.0 || -del < vc * s.hMax))
                  continue;
               // the depth this approach reaches (all the ball's kinetic energy into the squash): the stiffness while it lasts
               double lo = 0.0, hi = hL;
               for (int it = 0; it < 40; it++)
               {
                  const double mid = 0.5 * (lo + hi);
                  (SquashV(mid) < 0.5 * m * vc * vc ? lo : hi) = mid;
               }
               const double kExp = SquashK(std::max(del, lo));
               bool steady = false;   // touching, with its compression rate and the slip of the ball's contact point within vSteady
               if (s.vSteady > 0.0 && del > 0.0 && std::fabs(dDot) <= s.vSteady)
               {
                  const double gl = g3.Length();
                  const V3 nn = gl > 0.0 ? g3 * (-1.0 / gl) : V3(0.0, 0.0, 1.0);
                  const V3 u = vb + wb.Cross(nn * -R);
                  steady = (u - nn * u.Dot(nn)).Length() <= s.vSteady;
               }
               (steady ? kRest : kStrict) += kExp;
               if (del < 0.0 && dDot > 0.0)   // reached in about two steps, not overshot
                  lim(std::max(0.5 * -del / dDot, s.omegaH / std::sqrt(kExp / m)), 3);
            }
            {
               // the cord under the ball, against the playfield: as a backstop above, at the rate its lowest point can fall
               const double del = delNow[i], vc = std::max(std::fabs(vb.z) + wb.Length() * (R + r), s.vFloor);
               if (e >= 0 && (del > 0.0 || -del < vc * s.hMax))
               {
                  double lo = 0.0, hi = hL;
                  for (int it = 0; it < 40; it++)
                  {
                     const double mid = 0.5 * (lo + hi);
                     (SquashV(mid) < 0.5 * m * vc * vc ? lo : hi) = mid;
                  }
                  const double kExp = SquashK(std::max(del, lo));
                  kStrict += kExp;
                  if (del < 0.0)
                     lim(std::max(0.5 * -del / vc, s.omegaH / std::sqrt(kExp / m)), 3);
               }
            }
            if (kStrict > 0.0)
               lim(s.omegaH / std::sqrt(kStrict / m), 2);
            if (kRest > 0.0)
               lim(s.omegaHRest / std::sqrt((kStrict + kRest) / m), 7);
         }
         if (geo.hasPlate && !pinned)
         {
            const double I = Inertia(phi);
            double k = std::max(PlateVpp(phi), 0.0), dLdphi = 0.0;
            for (int e = 0; e < path3.n; e++)
               if (path3.e[e].kind == kFace)
               {
                  const int pn = path3.n;
                  const double dd = FaceDD(phi);
                  const V3 gk = G3.t[(e + pn - 1) % pn] - G3.t[e];
                  k += T * dd * dd * (1.0 / std::max(G3.l[(e + pn - 1) % pn], lMin) + 1.0 / std::max(G3.l[e], lMin));
                  dLdphi += (gk.x * geo.plate.n.x + gk.y * geo.plate.n.y) * dd;
               }
            k += kAx * dLdphi * dLdphi;
            if (k > 0.0)
               lim(s.omegaH / std::sqrt(k / I), 4);
            if (std::fabs(omega) > 1e-9)
               lim(s.turnMaxRad / std::fabs(omega), 5);
         }
         groupRule = rule = rl;
         return std::max(h, 1e-7);
      }

      bool BeginStep(const GroupStep::Ball *b, int n, double h, double *x, int off) override
      {
         Sync(b, n);
         gb = b;
         gd = MakeDofs();
         gPl = gd.phi >= 0;
         gq0 = Now();
         gL0 = Lcur;
         gI0 = Inertia(phi);
         gp0 = gI0 * omega;
         const bool coilOn = geo.hasPlate && driver.Coil();
         gCoil = coilOn;
         gIHigh1 = geo.hasPlate ? CurrentAfter(h, coilOn) : 0.0;
         gI2 = 0.5 * (iHigh * iHigh + gIHigh1 * gIHigh1);
         for (int i = 0; i < nb; i++)
            for (int j = 0; j < NumBackstops(); j++)
               gDel0[i * kMaxBs + j] = Overlap(i, j, gq0, gd, nullptr);
         double posScale = 0.1;
         for (int k = 0; k < geo.nPosts; k++)
            posScale = std::max(posScale, std::fabs(geo.posts[k].c.x) + std::fabs(geo.posts[k].c.y));
         const double eps = 2.220446049250313e-16;
         gNoiseL = 16.0 * eps * (path.n + 2) * (posScale + gL0);
         gNoiseE = gNoiseL * (Tension() + 1.0);
         FloorState(path3, G3, gFloorDel0);
         gNoiseD = 16.0 * eps * posScale;
         gMemo.ok = gMemoSlide.ok = gGradB0Ok = false;
         gSamePath = false;
         gNdEv = 0;
         for (int i = 0; i < nb; i++)
            gSlideOn[i] = false;
         for (int i = 0; i < nb; i++)
            gRootOk[i] = false;
         for (int i = 0; i < nb; i++)   // grip warm start: the last step's over this one (the band's push too: the hold's limit)
         {
            gP[i] = gPrate[i] * h;
            gTq[i] = gTqRate[i] * h;
            gBb[i] = gBc[i] = gTqB[i] = 0.0;
            gSlid[i] = false;
            gN0[i] = gFnTaken[i].Length();
            gNs0[i] = 0.0;
            gMode[i] = 0;
            gTauB[i] = gTauS[i] = 0.0;
            gSurfImp[i] = gSurfRate[i] * h;   // the surface's friction on a turning arc: warm start
            for (int j = 0; j < NumBackstops(); j++)
               gNs0[i] += SquashF(gDel0[i * kMaxBs + j]);
            gNs0[i] += SquashF(gFloorDel0[i]);
            // the two grips in series on the pressed arc: the weaker one slides (limits at the step's start: the ball's mu N R, the surface's mu Ns (R + 2 r))
            gStill[i] = holdOn[i] && gNs0[i] > 0.0 && set.mu * gN0[i] * R < set.mu * gNs0[i] * (R + 2.0 * r);
            const int e = holdOn[i] ? Find3(path3, kBall, i) : -1;
            if (e >= 0)   // the cord's direction at the ball: the arc's chord, or its run out for a vanishing arc
            {
               const V3 ch = path3.ct[e].X - path3.ct[e].E;
               gT1[i] = ch.Length() > 1e-9 ? ch * (1.0 / ch.Length()) : G3.t[e];
            }
            else
               gq0.held[i] = false;
         }
         if (gPl)
            x[off] = (omega + 0.5 * gAlphaLast * h) * h;
         return true;
      }

      bool InRange(const double *x, int off) const override
      {
         double dq[kMaxGD];
         LocalMove(x, off, dq);
         const Config q1 = Moved(dq, 1.0), qm = Moved(dq, 0.5);
         const double lim = hL * (1.0 - 1e-9);
         for (int i = 0; i < nb; i++)
            for (int j = 0; j < NumBackstops(); j++)
               if (Overlap(i, j, q1, gd, nullptr) >= lim || Overlap(i, j, qm, gd, nullptr) >= lim)
                  return false;
         return true;
      }

      // the step's forces at x: the cord and the squash as exact differences over the move (Gonzalez), their losses on the mid-step rates, the plate as in the
      // planar solve; the grip's impulses (held from UpdateGrip) added on the balls. The held grip across the cord (GroupStepDesign 11.3): each held ball's
      // turn about the cord's direction over the step (bb, from its mid-step spin); the arc turns with it (bc = bb) while the cord's pull on that turn stays within
      // mu x push x R, else it slides to where the pull is at that limit. The cord's length change is split exactly: the arcs' turn at the start positions
      // (L0 -> Lb: a torque on each held ball), then the move at the turned arcs (Lb -> L1: the forces on the balls and the plate).
      bool Eval(const double *x, int off, double h, GV3 *F, GV3 *M, double *Res, double (*J)[GroupStep::kMaxUnk]) override
      {
         double dq[kMaxGD];
         const int nd = LocalMove(x, off, dq), iphi = 3 * nb;
         Config q1 = Moved(dq, 1.0), qm = Moved(dq, 0.5);
         int nh = 0, nTurn = 0;
         for (int i = 0; i < nb; i++)
         {
            gBb[i] = gBc[i] = gTqB[i] = 0.0;
            gSlid[i] = false;
            if (!gq0.held[i])
               continue;
            nh++;
            const V3 wm = ToBand(gb[i].w + GV3(x[6 * i + 3], x[6 * i + 4], x[6 * i + 5]) * 0.5);
            gBb[i] = wm.Dot(gT1[i]) * h;
            // Each held arc by its step's mode (BeginStep; GroupStepDesign 11.10). Pressed onto a surface whose grip is the stronger (gStill): the arc stays where it
            // is (bc = 0) and the ball rubs on it (the grip's friction against a still cord, UpdateGrip). Otherwise the arc turns with the ball while the cord's
            // pull stays within the ball's limit (the surface's weaker friction, if any, acting on that turn in UpdateGrip), else it slides on the ball the way the
            // cord pulls it until the pull is down to the limit (Coulomb).
            if (gStill[i])
               gMode[i] = 4;
            else
            {
               gBc[i] = gBb[i];
               gMode[i] = gNs0[i] > 0.0 ? 6 : 2;
               nTurn++;
            }
         }
         double Lb = gL0, gbt[kMaxBalls] = {}, Fel[kMaxBalls] = {}, gBe[kMaxBalls] = {}, LbS[kMaxPosts], LbD[kMaxBalls], Vb = SpanV(spanL) + gFloorV;
         for (int j = 0; j < geo.nPosts; j++)
            LbS[j] = spanL[j];
         for (int i = 0; i < nb; i++)
            LbD[i] = gFloorDel0[i];
         const auto hit = [&](const BetaMemo &M) {
            if (!M.ok)
               return false;
            for (int i = 0; i < nb; i++)
               if (M.bc[i] != gBc[i])
                  return false;
            return true;
         };
         const auto BetaPart = [&]() -> bool {   // the arcs' turn at the start positions: Lb, its discrete gradient over the turns, the elastic pull on each
            const BetaMemo *M = hit(gMemo) ? &gMemo : hit(gMemoSlide) ? &gMemoSlide : nullptr;
            if (!M)
            {
               struct Restore { int &c; ~Restore() { c = 0; } } rs { solveCaller };
               solveCaller = 3;
               Config qb = gq0;
               TurnHeld(qb, gBc);
               sPb = path3;
               if (!Solve3(sPb, qb, sGb))
                  return false;
               if (!gGradB0Ok)
               {
                  double w0[kMaxPosts], wf0[kMaxBalls];
                  for (int j = 0; j < geo.nPosts; j++)
                     w0[j] = SpanT(j, spanL[j]);
                  for (int i = 0; i < nb; i++)
                     wf0[i] = SquashF(gFloorDel0[i]);
                  if (!GradBeta(path3, G3, gq0, gGradB0, w0, wf0))
                     return false;
                  gGradB0Ok = true;
               }
               BetaMemo &N = gMemo;
               N.ok = false;
               MakeSpanMap(sPb, sGb, sMapB);
               SpanLengths(sPb, sGb, sMapB, N.LbS);
               double wb[kMaxPosts], wfb[kMaxBalls];
               FloorState(sPb, sGb, N.LbD);
               for (int i = 0; i < nb; i++)
               {
                  if (N.LbD[i] >= hL * (1.0 - 1e-9))
                     return false;
                  wfb[i] = SquashF(N.LbD[i]);
               }
               for (int j = 0; j < geo.nPosts; j++)
                  wb[j] = SpanT(j, N.LbS[j]);
               if (!GradBeta(sPb, sGb, qb, N.ge, wb, wfb))   // (the energy's derivative over the turn: each span's tension on its share)
                  return false;
               if (debugGrad)
               {
                  double gn[kMaxBalls], gx[kMaxBalls];
                  if (GradBetaNum(path3, qb, gBc, gn) && GradBeta(sPb, sGb, qb, gx))
                     for (int i = 0; i < nb; i++)
                        gradErrBeta = std::max(gradErrBeta, std::fabs(gx[i] - gn[i]));
               }
               N.Lb = sGb.L;
               double gm[kMaxBalls], Wb[kMaxBalls];
               for (int i = 0; i < nb; i++)
               {
                  N.bc[i] = gBc[i];
                  gm[i] = 0.5 * (gGradB0[i] + N.ge[i]);   // the arcs' derivative over the turn: its two ends (second order, as the mid-turn's, with no path solve of its own)
                  Wb[i] = 1.0 / Ib;
               }
               N.Vb = SpanV(N.LbS) + FloorVOf(N.LbD);
               Gonzalez(gm, N.Vb - SpanV(spanL) - gFloorV, gNoiseE, gBc, Wb, nb, N.gbt);
               for (int i = 0; i < nb; i++)
                  N.Fel[i] = -N.gbt[i];
               N.ok = true;
               M = &N;
            }
            Lb = M->Lb;
            for (int j = 0; j < geo.nPosts; j++)
               LbS[j] = M->LbS[j];
            Vb = M->Vb;
            for (int i = 0; i < nb; i++)
               LbD[i] = M->LbD[i];
            for (int i = 0; i < nb; i++)
            {
               gbt[i] = M->gbt[i];
               Fel[i] = M->Fel[i];
               gBe[i] = M->ge[i];
            }
            return true;
         };
         if (nh > 0)
         {
            // A round cord on a rolling ball (GroupStepDesign 11.28): at the contact the ball's surface moves R bb; the cord's centre moves S bc round the ball (S = R + r, bc the arc's turn), and
            // the cord rolls about its own axis by dpsi, against its twist: without slip S bc - r dpsi = R bb. The energy's stationary point gives one force at the contact: the pull on the arc
            // Fel = S f, the twist's torque M = r f (M the mid value of k psi), and the ball's torque R f. Past mu N the contact slides: f at the limit against the slip, which turns bc and psi
            // back until both match it: the loss is f x the slip, never negative. All the arcs on one cord solve together, as they pull on each other through the runs between them: a
            // quasi-Newton solve over the turns (matrix started from each arc's own slope, corrected by the moves made, Broyden); then which arcs slide is checked.
            if (nTurn > 0)
            {
               const double S = R + r, Sr = S / r, Rr = R / r;
               const auto Fric = [&](int i) { return set.mu * gN0[i]; };
               const auto Dpsi = [&](int i) { return (S * gBc[i] - R * gBb[i]) / r; };
               const auto Mid = [&](int i) { return kTw[i] * (psi[i] + 0.5 * Dpsi(i)) + cTw[i] / h * Dpsi(i); };   // the twist's torque with the rubber's own loss
               const auto Res = [&](int i) { return gMode[i] == 3 ? Fel[i] - gSlideF[i] : Fel[i] - Sr * Mid(i); };
               // the first guess: warm from this step's last answer, else from the pull at no turn and the cord's bend
               bool anyCold = false;
               for (int i = 0; i < nb; i++)
                  if (gq0.held[i] && !gStill[i])
                  {
                     gBc[i] = 0.0;
                     anyCold = anyCold || !gRootOk[i];
                  }
               if (anyCold && !BetaPart())
               {
                  gdiag.settleFails++;
                  return false;
               }
               for (int i = 0; i < nb; i++)
               {
                  if (!gq0.held[i] || gStill[i])
                     continue;
                  const int e = Find3(path3, kBall, i);
                  const int pn = path3.n;
                  const double lin = e >= 0 ? std::max(G3.l[(e + pn - 1) % pn], 1e-3) : 1.0, lout = e >= 0 ? std::max(G3.l[e], 1e-3) : 1.0;
                  const double Kb = 0.5 * Tension() * S * S * (1.0 / lin + 1.0 / lout), kq = (0.5 * kTw[i] + cTw[i] / h) / (r * r), den = Kb + kq * S * S;
                  gRootA[i] = den > 0.0 ? kq * S * R / den : R / S;
                  if (gRootOk[i])
                     gBc2[i] = gRoot[i] + gRootA[i] * (gBb[i] - gRootBb[i]);
                  else
                     gBc2[i] = den > 0.0 ? (Fel[i] - Sr * kTw[i] * psi[i] + kq * S * R * gBb[i]) / den : R / S * gBb[i];
               }
               for (int i = 0; i < nb; i++)
                  if (gq0.held[i] && !gStill[i])
                     gBc[i] = gBc2[i];
               if (!BetaPart())
               {
                  gdiag.settleFails++;
                  return false;
               }
               bool lawMet = false;
               for (int setPass = 0; setPass < kArcSetPasses && !lawMet; setPass++)
               {
                  int act[kMaxBalls], na = 0;
                  for (int i = 0; i < nb; i++)
                     if (gq0.held[i] && !gStill[i])
                        act[na++] = i;
                  const auto Met = [&]() {   // each arc's energy error in a step: far below the books' limit over a whole run
                     for (int a = 0; a < na; a++)
                        if (std::fabs(Res(act[a]) * gBc[act[a]]) > kArcTolJ)
                           return false;
                     return true;
                  };
                  double B[kMaxBalls][kMaxBalls] = {}, res[kMaxBalls], dx[kMaxBalls];
                  for (int a = 0; a < na; a++)
                  {
                     const int i = act[a];
                     const double slope = std::fabs(gBc[i]) > 1e-9 ? -(gBe[i] + Fel[i]) / gBc[i] : 0.0;
                     B[a][a] = (slope < 0.0 ? slope : -Tension() * S) - (gMode[i] == 3 ? 0.0 : Sr * Sr * (0.5 * kTw[i] + cTw[i] / h));
                     res[a] = Res(i);
                  }
                  bool innerOk = Met();
                  double brP = 0.0, brN = 0.0, resPrev = INFINITY;
                  bool hasP = false, hasN = false;
                  for (int it = 0; it < kArcIters && !innerOk; it++)
                  {
                     // B dx = -res
                     double A[kMaxBalls][kMaxBalls + 1];
                     for (int a = 0; a < na; a++)
                     {
                        for (int b = 0; b < na; b++)
                           A[a][b] = B[a][b];
                        A[a][na] = -res[a];
                     }
                     bool singular = false;
                     for (int c = 0; c < na && !singular; c++)
                     {
                        int pv = c;
                        for (int a = c + 1; a < na; a++)
                           if (std::fabs(A[a][c]) > std::fabs(A[pv][c]))
                              pv = a;
                        if (std::fabs(A[pv][c]) < 1e-30)
                        {
                           singular = true;
                           break;
                        }
                        if (pv != c)
                           for (int b = c; b <= na; b++)
                              std::swap(A[pv][b], A[c][b]);
                        for (int a = c + 1; a < na; a++)
                        {
                           const double f = A[a][c] / A[c][c];
                           for (int b = c; b <= na; b++)
                              A[a][b] -= f * A[c][b];
                        }
                     }
                     if (singular)
                        break;
                     double big = 0.0;
                     for (int a = na - 1; a >= 0; a--)
                     {
                        double v = A[a][na];
                        for (int b = a + 1; b < na; b++)
                           v -= A[a][b] * dx[b];
                        dx[a] = v / A[a][a];
                        big = std::max(big, std::fabs(dx[a]));
                     }
                     if (big > kArcMaxMoveRad)   // a far guess stays in the arc's neighbourhood
                        for (int a = 0; a < na; a++)
                           dx[a] *= kArcMaxMoveRad / big;
                     // The pull on an arc rises steeply as its turn crosses the point where the runs are straight (a ball close to a post: a step in the pull over a few um of turn), and a secant can circle such
                     // a step for ever. The residual is monotone there. One arc: the root stays inside the bracket its signs give. Several: a move that does not lower the residuals is shortened.
                     if (na == 1)
                     {
                        const double ra = res[0], bc0 = gBc[act[0]];
                        (ra > 0.0 ? brP : brN) = bc0;
                        (ra > 0.0 ? hasP : hasN) = true;
                        if (hasP && hasN)
                        {
                           const double t = bc0 + dx[0], lo = std::min(brP, brN), hi = std::max(brP, brN);
                           if (!(t > lo && t < hi) || std::fabs(ra) > 0.5 * resPrev)
                              dx[0] = 0.5 * (lo + hi) - bc0;
                        }
                        resPrev = std::fabs(ra);
                     }
                     double n0 = 0.0;
                     for (int a = 0; a < na; a++)
                        n0 += res[a] * res[a];
                     for (int bt = 0;; bt++)
                     {
                        for (int a = 0; a < na; a++)
                           gBc[act[a]] += dx[a];
                        if (!BetaPart())
                        {
                           gdiag.settleFails++;
                           return false;
                        }
                        if (na == 1 || bt >= 6)
                           break;
                        double nn = 0.0;
                        for (int a = 0; a < na; a++)
                        {
                           const double rn = Res(act[a]);
                           nn += rn * rn;
                        }
                        if (nn < n0)
                           break;
                        for (int a = 0; a < na; a++)
                        {
                           gBc[act[a]] -= dx[a];
                           dx[a] *= 0.5;
                        }
                     }
                     double dd = 0.0, Bdx[kMaxBalls], dres[kMaxBalls];
                     for (int a = 0; a < na; a++)
                     {
                        const double rn = Res(act[a]);
                        dres[a] = rn - res[a];
                        res[a] = rn;
                        dd += dx[a] * dx[a];
                        Bdx[a] = 0.0;
                        for (int b = 0; b < na; b++)
                           Bdx[a] += B[a][b] * dx[b];
                     }
                     if (dd > 0.0)
                        for (int a = 0; a < na; a++)
                           for (int b = 0; b < na; b++)
                              B[a][b] += (dres[a] - Bdx[a]) * dx[b] / dd;
                     innerOk = Met();
                  }
                  if (!innerOk)
                  {
                     gdiag.settleFails++;
                     return false;
                  }
                  // which arcs slide: one past its limit starts to; one whose slip has turned against its friction stops
                  lawMet = true;
                  for (int a = 0; a < na; a++)
                  {
                     const int i = act[a];
                     if (gMode[i] == 3)
                     {
                        const double dp = (gSlideM[i] - kTw[i] * psi[i]) / std::max(0.5 * kTw[i] + cTw[i] / h, 1e-300), slip = S * gBc[i] - r * dp - R * gBb[i];
                        if (slip * gArcDir[i] < 0.0)
                        {
                           gMode[i] = gNs0[i] > 0.0 ? 6 : 2;
                           lawMet = false;
                        }
                     }
                     else if (std::fabs(Mid(i)) > r * Fric(i))
                     {
                        gMode[i] = 3;
                        gArcDir[i] = Mid(i) >= 0.0 ? 1.0 : -1.0;
                        gSlideF[i] = gArcDir[i] * Fric(i) * S;
                        gSlideM[i] = gArcDir[i] * Fric(i) * r;
                        gSlideTq[i] = gArcDir[i] * Fric(i) * R;
                        lawMet = false;
                     }
                  }
               }
               if (!lawMet)
               {
                  gdiag.settleFails++;
                  return false;
               }
               for (int i = 0; i < nb; i++)
                  if (gq0.held[i] && !gStill[i])
                  {
                     gRoot[i] = gBc[i];
                     gRootBb[i] = gBb[i];
                     gRootOk[i] = true;
                     if (gMode[i] == 3)
                     {
                        gTwM[i] = gSlideM[i];
                        gPsi1[i] = psi[i] + (gSlideM[i] - kTw[i] * psi[i]) / std::max(0.5 * kTw[i] + cTw[i] / h, 1e-300);
                     }
                     else
                     {
                        gTwM[i] = Mid(i);
                        gPsi1[i] = psi[i] + Dpsi(i);
                     }
                  }
               if (hit(gMemo))
                  gMemoSlide = gMemo;
            }
            // every arc's turn is now set: the final arcs' part (the gradient and the pulls at the chosen turns; already in hand unless an arc slid)
            if (nTurn > 0 && !BetaPart())
            {
               gdiag.settleFails++;
               return false;
            }
            for (int i = 0; i < nb; i++)
               gSlid[i] = gq0.held[i] && gMode[i] != 2 && gMode[i] != 6;   // the spin Jacobian's held term: only an arc turning with the ball
            TurnHeld(q1, gBc);
            TurnHeld(qm, gBc);
         }
         gPath1 = path3;
         Path3 &pm = sPm;
         pm = path3;
         Geom3 &Gm = sGm;
         passStats.evals++;
         solveCaller = 1;
         const bool okEnd = gSamePath ? LeaveOnly(gPath1, q1, gG1) : Settle3(gPath1, q1, gG1), okMid = okEnd && (gSamePath ? LeaveOnly(pm, qm, Gm) : Settle3(pm, qm, Gm));
         solveCaller = 0;
         if (!okMid)
         {
            gdiag.settleFails++;
            gdiag.settleFailMid = okEnd;
            gdiag.settleFailN = okEnd ? pm.n : gPath1.n;
            Path3 pl = path3;   // the failing settle again, logged
            Geom3 gl;
            Settle3(pl, okEnd ? qm : q1, gl, &settleLog);
            return false;
         }
         const double L1 = gG1.L, dL = L1 - gL0;   // the whole change (the cord's stretch loss); the move's part is L1 - Lb
         double W[kMaxGD], Fl[kMaxGD], gLm[kMaxGD], gL[kMaxGD], dqc[kMaxGD];
         double H[kMaxGD][kMaxGD];
         for (int a = 0; a < nd; a++)
         {
            W[a] = a < iphi ? 1.0 / m : 1.0 / gI0;
            for (int b2 = 0; b2 < nd; b2++)
               H[a][b2] = 0.0;
         }
         // the cord: its gradient at mid-step, its discrete gradient over the balls it touches during the move (a ball off the cord takes none of the
         // correction) and the plate
         solveCaller = 2;
         MakeSpanMap(pm, Gm, sMapM);
         SpanLengths(pm, Gm, sMapM, sLsM);
         // the playfield under the cord, at mid-step and at the end: its pull, with the squash law's damping on the move
         double delM[kMaxBalls], del1[kMaxBalls], wFloor[kMaxBalls], gFl[kMaxGD], dLossF = 0.0, kFloor[kMaxBalls];
         FloorState(pm, Gm, delM);
         FloorState(gPath1, gG1, del1);
         for (int i = 0; i < nb; i++)
         {
            if (delM[i] >= hL * (1.0 - 1e-9) || del1[i] >= hL * (1.0 - 1e-9))
            {
               gdiag.settleFails++;
               return false;
            }
            const double fm = SquashF(delM[i]), cD = delM[i] > 0.0 ? 1.5 * set.alphaSPerM * fm : 0.0, dd = std::max(del1[i], 0.0) - std::max(LbD[i], 0.0), cf = cD * dd / h;
            wFloor[i] = fm + cf;
            dLossF += cf * dd;
            kFloor[i] = 0.5 * h * SquashK(delM[i]) + cD;
         }
         gFloorLoss = dLossF;
         const bool okGrad = GradL3(pm, Gm, qm, gPl, gLm, sGs, &sMapM, wFloor, gFl);
         if (debugGrad && okGrad)
         {
            double gn[kMaxGD];
            if (GradL3Num(pm, qm, gPl, gn))
               for (int a = 0; a < 3 * nb + 1; a++)
               {
                  const double er = std::fabs(gLm[a] - gn[a]);
                  if (a == 3 * nb)
                     gradErrPhi = std::max(gradErrPhi, er);
                  else
                     gradErrQ = std::max(gradErrQ, er);
               }
         }
         solveCaller = 0;
         if (!okGrad)
         {
            gdiag.settleFails++;
            return false;
         }
         for (int i = 0; i < nb; i++)
         {
            const bool on = Find3(path3, kBall, i) >= 0 || Find3(gPath1, kBall, i) >= 0 || Find3(pm, kBall, i) >= 0 || Find3(path3, kBallBack, i) >= 0 ||
                            Find3(gPath1, kBallBack, i) >= 0 || Find3(pm, kBallBack, i) >= 0;
            for (int a = 0; a < 3; a++)
               dqc[3 * i + a] = on ? dq[3 * i + a] : 0.0;
         }
         if (gPl)   // the cord's energy does not depend on a plate it is not on: as a ball off the cord, it takes none of the correction
         {
            bool onFace = false;
            for (int k = 0; k < path3.n && !onFace; k++)
               onFace = path3.e[k].kind == kFace;
            for (int k = 0; k < gPath1.n && !onFace; k++)
               onFace = gPath1.e[k].kind == kFace;
            for (int k = 0; k < pm.n && !onFace; k++)
               onFace = pm.e[k].kind == kFace;
            dqc[iphi] = onFace ? dq[iphi] : 0.0;
         }
         // one scalar over the spans, so one discrete gradient closes the energy exactly: the cord's stored energy change, plus each span's stretch loss written as its
         // frozen rate times its length (its change is the loss itself)
         MakeSpanMap(gPath1, gG1, sMap1);
         SpanLengths(gPath1, gG1, sMap1, sLs1);
         double gPhi[kMaxGD], dLossS = 0.0;
         for (int a = 0; a < nd; a++)
            gPhi[a] = 0.0;
         for (int j = 0; j < geo.nPosts; j++)
         {
            const double dl = sLs1[j] - LbS[j], cj = SpanCv(j) * dl / h, wj = SpanT(j, sLsM[j]) + cj;
            dLossS += cj * dl;
            for (int a = 0; a < nd; a++)
               gPhi[a] += wj * sGs[j][a];
         }
         for (int a = 0; a < nd; a++)
            gPhi[a] += gFl[a];
         Gonzalez(gPhi, SpanV(sLs1) + FloorVOf(del1) - Vb + dLossS + dLossF, gNoiseE, dqc, W, nd, gL);   // the stretch loss acts on the move's share of the length change (the arcs' turn takes none)
         for (int a = 0; a < nd; a++)
            Fl[a] = -gL[a];
         gCordLoss = dLossS;
         // The correction's own slope in the Newton matrix. Its force is -(e / qq) W dq, e = df - g . dq, qq = dq^T W dq (e changes only to third order in the move). Left out, a light body (the plate) overshoots
         // every round when e is large: the row's true slope is 2 I / h + h e / (I qq), the matrix had the first term only.
         {
            const double df = SpanV(sLs1) + FloorVOf(del1) - Vb + dLossS + dLossF;
            double gdq = 0.0, qq = 0.0, Wd[kMaxGD];
            for (int a = 0; a < nd; a++)
            {
               gdq += gPhi[a] * dqc[a];
               Wd[a] = W[a] * dqc[a];
               qq += dqc[a] * Wd[a];
            }
            if (qq > 1e-300 && std::fabs(df - gdq) > gNoiseE)
            {
               const double cf = (df - gdq) / qq;
               for (int a = 0; a < nd; a++)
                  if (dqc[a] != 0.0)
                  {
                     H[a][a] += h * cf * W[a];
                     for (int b2 = 0; b2 < nd; b2++)
                        H[a][b2] -= 2.0 * h * cf / qq * Wd[a] * Wd[b2];
                  }
            }
         }
         // the held balls: the cord's elastic pull on the arc's turn, a torque about the cord's direction; a slid arc's torque is the friction itself, at its limit
         // against the slip (GroupStepDesign 11.19), so its loss is never negative by construction: the search stops where the pull meets that limit, but at a kink
         // in the length (a run squashed to nothing on a post) the pull can jump across it
         gHoldLoss = 0.0;
         for (int i = 0; i < nb; i++)
            if (gq0.held[i])
            {
               // the cord's pull on the arc goes to the ball while the arc turns with it or slides on it (the friction is then that pull, against the slip); an
               // arc held still by the surface puts it on the surface, which does no work
               if (gStill[i])
               {
                  gTqB[i] = 0.0;
                  gPsi1[i] = psi[i];
                  continue;
               }
               gTqB[i] = gMode[i] == 3 ? gSlideTq[i] : gTwM[i] * R / r;
               gHoldLoss += -gTqB[i] * gBb[i] + Fel[i] * gBc[i] - 0.5 * kTw[i] * (gPsi1[i] * gPsi1[i] - psi[i] * psi[i]);   // what the ball gave, less what the cord took: f x the slip and the twist's own loss
            }
         {
            for (int j = 0; j < geo.nPosts; j++)
            {
               const double sj = 0.5 * h * (sLsM[j] > spanS[j] ? EA / spanS[j] : 0.0) + SpanCv(j);
               for (int a = 0; a < nd; a++)
                  for (int b2 = 0; b2 < nd; b2++)
                     H[a][b2] += sj * sGs[j][a] * sGs[j][b2];
            }
            for (int i = 0; i < nb; i++)
               H[3 * i + 2][3 * i + 2] += kFloor[i];
            double sRun[kMaxPath];
            for (int k = 0; k < pm.n; k++)
            {
               const int sp = sMapM.partSpan[k][0];
               sRun[k] = 0.5 * h * SpanT(sp, sLsM[sp]);
            }
            AddHessL3(pm, Gm, qm, gPl, sRun, H);
         }
         V3 Fn3[kMaxBalls];
         for (int i = 0; i < nb; i++)
         {
            Fn3[i] = V3(Fl[3 * i], Fl[3 * i + 1], Fl[3 * i + 2]);
            gdiag.cordF[i] = Fn3[i];
            gdiag.wrap[i] = 0.0;
            gdiag.clear[i] = 0.0;
            const int e = Find3(pm, kBall, i);
            if (e >= 0)
               gdiag.wrap[i] = Gm.wrap[e];
         }
         // the squash: ball i against backstop j, its discrete gradient over ball i (and the plate)
         gBsLoss = gFloorLoss;
         gMaxOv = 0.0;
         for (int i = 0; i < nb; i++)
            for (int j = 0; j < NumBackstops(); j++)
            {
               double g2[kMaxDof], gz = 0.0;
               const double d0 = gDel0[i * kMaxBs + j], d1 = Overlap(i, j, q1, gd, nullptr), dm = Overlap(i, j, qm, gd, g2, &gz);
               gMaxOv = std::max(gMaxOv, d1);
               if (d0 <= 0.0 && d1 <= 0.0 && dm <= 0.0)
                  continue;
               int idx[4] = { 3 * i, 3 * i + 1, 3 * i + 2, iphi };
               const int na = gPl ? 4 : 3;
               double gm4[4] = { g2[2 * i], g2[2 * i + 1], gz, gPl ? g2[gd.phi] : 0.0 }, dq4[4], W4[4], gd4[4];
               for (int a = 0; a < na; a++)
               {
                  dq4[a] = dq[idx[a]];
                  W4[a] = W[idx[a]];
               }
               Gonzalez(gm4, d1 - d0, gNoiseD, dq4, W4, na, gd4);
               const double cD = dm > 0.0 ? 1.5 * set.alphaSPerM * SquashF(dm) : 0.0;
               const double fv = SquashChord(d0, d1, dm) + cD * (d1 - d0) / h;
               for (int a = 0; a < na; a++)
                  Fl[idx[a]] -= fv * gd4[a];
               Fn3[i] = Fn3[i] - V3(gd4[0], gd4[1], gd4[2]) * fv;
               gBsLoss += cD * (d1 - d0) * (d1 - d0) / h;
               const double kd = SquashK(dm);
               for (int a = 0; a < na; a++)
                  for (int b2 = 0; b2 < na; b2++)
                     H[idx[a]][idx[b2]] += (0.5 * h * kd + cD) * gm4[a] * gm4[b2];
            }
         // the plate (as the planar solve)
         gPlateLoss = gCoilW = 0.0;
         if (gPl)
         {
            const double dphi = dq[iphi], phm = qm.phi, I1 = Inertia(q1.phi);
            const double Jb = 0.5 * (1.0 / gI0 + 1.0 / I1);
            const double pEnd = 2.0 * dphi / (h * Jb) - gp0;
            double dTdphi;
            if (std::fabs(dphi) > 1e-7)
               dTdphi = (1.0 / I1 - 1.0 / gI0) / dphi * 0.25 * (gp0 * gp0 + pEnd * pEnd);
            else
            {
               const double e = 1e-6;
               dTdphi = (1.0 / Inertia(phm + e) - 1.0 / Inertia(phm - e)) / (2.0 * e) * 0.25 * (gp0 * gp0 + pEnd * pEnd);
            }
            const double Vp = std::fabs(dphi) > 1e-7 ? (PlateV(q1.phi) - PlateV(gq0.phi)) / dphi : PlateVp(phm);
            double tauS, cR;
            RestStop(phm, tauS, cR);
            const double om = dphi / h;
            const double tauD = tauS > 0.0 ? std::max(-cR * om, -tauS) : 0.0;
            const double tauC = gCoil || iHigh > 0.0 || gIHigh1 > 0.0 ? CoilTorquePerI2(phm) * gI2 : 0.0;
            const double bd = set.sling.mech.bDamping;
            Fl[iphi] += -Vp + tauD - bd * om + tauC - dTdphi;
            gPlateLoss = -(tauD - bd * om) * dphi;
            gCoilW = tauC * dphi;
            H[iphi][iphi] += 2.0 / (h * Jb) + 0.5 * h * std::max(PlateVpp(phm), 0.0) + bd + cR;
            Res[off] = pEnd - gp0 - h * Fl[iphi];
         }
         // onto the balls in the group's frame, with the grip's impulses over the step
         for (int i = 0; i < nb; i++)
         {
            F[i] = F[i] + ToWorld(V3(Fl[3 * i], Fl[3 * i + 1], Fl[3 * i + 2])) + ToWorld(gP[i]) * (1.0 / h);
            M[i] = M[i] + ToWorld(gTq[i]) * (1.0 / h) + ToWorld(gT1[i] * gTqB[i]);
            gFn[i] = Fn3[i];
            gdiag.squashF[i] = Fn3[i] - gdiag.cordF[i];
            gdiag.heldTq[i] = gTqB[i];
            gdiag.held[i] = gq0.held[i] ? 1 : 0;
            gdiag.gripF[i] = gP[i] * (1.0 / h);
            gdiag.push[i] = Fn3[i].Length();
         }
         if (J)
         {
            // the band's rows and columns turned into the group's: a world move dc_a changes band coordinate p by E_p . e_a
            const V3 E[3] = { frame.ex, frame.ey, frame.ez };
            const auto Ec = [&](int p, int a) { return a == 0 ? E[p].x : a == 1 ? E[p].y : E[p].z; };
            for (int i = 0; i < nb; i++)
               for (int j = 0; j < nb; j++)
                  for (int a = 0; a < 3; a++)
                     for (int b2 = 0; b2 < 3; b2++)
                     {
                        double s2 = 0.0;
                        for (int p = 0; p < 3; p++)
                           for (int qq = 0; qq < 3; qq++)
                              s2 += Ec(p, a) * H[3 * i + p][3 * j + qq] * Ec(qq, b2);
                        J[6 * i + a][6 * j + b2] += s2;
                     }
            if (gPl)
            {
               for (int j = 0; j < nb; j++)
                  for (int b2 = 0; b2 < 3; b2++)
                  {
                     double sr = 0.0, sc = 0.0;
                     for (int p = 0; p < 3; p++)
                     {
                        sr += H[iphi][3 * j + p] * Ec(p, b2);
                        sc += Ec(p, b2) * H[3 * j + p][iphi];
                     }
                     J[off][6 * j + b2] += sr;
                     J[6 * j + b2][off] += sc;
                  }
               J[off][off] += H[iphi][iphi];
            }
            // a friction at its limit: mu x the push, so it changes with the move as the push does (-n . H / h, h from the impulse's h): its direction held
            for (int i = 0; i < nb; i++)
            {
               if (!gSlideOn[i])
                  continue;
               const double Nn = gFn[i].Length();
               if (Nn < 1e-9)
                  continue;
               const V3 nh = gFn[i] * (1.0 / Nn);
               double hv[kMaxGD];
               for (int b2 = 0; b2 < nd; b2++)
                  hv[b2] = nh.x * H[3 * i][b2] + nh.y * H[3 * i + 1][b2] + nh.z * H[3 * i + 2][b2];
               const V3 dirF = gSlideDir[i], dirT = (gN3[i] * -R).Cross(gSlideDir[i]);
               const GV3 fw = ToWorld(dirF), tw = ToWorld(dirT);
               const double fv[6] = { fw.x, fw.y, fw.z, tw.x, tw.y, tw.z };
               for (int a = 0; a < 6; a++)
               {
                  const double c0 = set.mu * fv[a];
                  for (int j = 0; j < nb; j++)
                     for (int b2 = 0; b2 < 3; b2++)
                     {
                        double sr = 0.0;
                        for (int pp = 0; pp < 3; pp++)
                           sr += hv[3 * j + pp] * Ec(pp, b2);
                        J[6 * i + a][6 * j + b2] += c0 * sr;
                     }
                  if (gPl)
                     J[6 * i + a][off] += c0 * hv[iphi];
               }
            }
            gNdEv = nd;
            for (int a = 0; a < nd; a++)
            {
               gDqEv[a] = dq[a];
               for (int b2 = 0; b2 < nd; b2++)
                  gHev[a][b2] = H[a][b2];
            }
            // a held ball's spin: the pull's change with the arc's turn (stretch and the cord's bend round the ball, about T S^2 (1/l_in + 1/l_out)), the turn
            // moving with half the spin change (bb = (w0 + dw / 2) . t1 h)
            const double S = R + r;   // (the held arcs carry no stretch loss: the cv term is left out)
            MakeSpanMap(path3, G3, sMapB);
            for (int i = 0; i < nb; i++)
            {
               const int e = gq0.held[i] && !gSlid[i] ? Find3(path3, kBall, i) : -1;
               if (e < 0)
                  continue;
               const int spanE = sMapB.arcSpan[e];
               const double Tm = SpanT(spanE, sLs1[spanE]), kAx = sLs1[spanE] > spanS[spanE] ? EA / spanS[spanE] : 0.0, gbL = Tm > 0.0 ? gbt[i] / Tm : 0.0;
               const int pn = path3.n;
               const double lin = std::max(G3.l[(e + pn - 1) % pn], 1e-3), lout = std::max(G3.l[e], 1e-3);
               const double kl = 0.5 * kAx * gbL * gbL + Tm * S * S * (1.0 / lin + 1.0 / lout), kw = (0.5 * kTw[i] + cTw[i] / h) / (r * r);
               const double kb = kl * kw * R * R / std::max(kl + kw * S * S, 1e-30) * h * 0.5 * h;   // the bend and the twist in series
               const GV3 tw = ToWorld(gT1[i]);
               const double tv[3] = { tw.x, tw.y, tw.z };
               for (int a = 0; a < 3; a++)
                  for (int b2 = 0; b2 < 3; b2++)
                     J[6 * i + 3 + a][6 * i + 3 + b2] += kb * tv[a] * tv[b2];
            }
         }
         return true;
      }

      // the grip (GroupStepDesign 10.2): the planar solve's stick / slide friction, twist and rolling resistance turned into the 3D contact plane of the band's
      // push on the ball (n from the push, t1 along the cord, t2 = n x t1 across it); rates at mid-step less the grip's own share, the free ball's compliance (the
      // fixed point needs no more: at it the impulses meet their limits on the true rates). Sticking means no slip at the step's END (GroupStepDesign 11.12): the
      // impulse is chosen on the end rate (mid-step rate less half the start rate, the slip taken as changing evenly over the step), which spends the start slip;
      // where that impulse would do positive work (other forces turning the slip over within the step) the mid-step rule is kept (no slip at mid-step: no work).
      // A grip that held only the mid-step slip at zero let the slip reverse every step at full size (MEASURED, notes-group-band-chatter-2026-10-29.md).
      double UpdateGrip(const double *x, double h) override
      {
         double ch = 0.0;
         const double k1 = 0.5 / m + 0.5 * R * R / Ib;
         // the band's push at this move: the last Eval's, and its change over the Newton move (the push the grip uses is then as new as the rates it uses)
         double dqNow[kMaxGD];
         LocalMove(x, 6 * nb, dqNow);
         for (int i = 0; i < nb; i++)
         {
            V3 FnNow = gFn[i];
            for (int a = 0; a < 3; a++)
            {
               double df = 0.0;
               for (int b2 = 0; b2 < gNdEv; b2++)
                  df -= gHev[3 * i + a][b2] * (dqNow[b2] - gDqEv[b2]) / h;
               (a == 0 ? FnNow.x : a == 1 ? FnNow.y : FnNow.z) += df;
            }
            const double N = FnNow.Length();
            V3 Pn, Tq;
            double Mtw = 0.0, Mr = 0.0;
            V3 n3, t2;
            gSlideOn[i] = false;
            if (N >= 1e-9)
            {
               n3 = FnNow * (1.0 / N);
               V3 t1(-n3.y, n3.x, 0.0);
               const double tl = t1.Length();
               if (tl > 1e-9)
                  t1 = t1 * (1.0 / tl);
               else   // pushed straight along the band's normal: any direction in the contact plane
               {
                  t1 = V3(1.0, 0.0, 0.0) - n3 * n3.x;
                  t1 = t1 * (1.0 / t1.Length());
               }
               t2 = n3.Cross(t1);
               const V3 rc = n3 * -R;
               const V3 vm = ToBand(GV3(x[6 * i], x[6 * i + 1], x[6 * i + 2])) * (1.0 / h) - gP[i] * (0.5 / m);
               const V3 wm = ToBand(gb[i].w + GV3(x[6 * i + 3], x[6 * i + 4], x[6 * i + 5]) * 0.5) - gTq[i] * (0.5 / Ib);
               const V3 sl = vm + wm.Cross(rc);
               // the slip and turn rates at the step's start (the ball as it began the step)
               const V3 ws = ToBand(gb[i].w), sls = ToBand(gb[i].v) + ws.Cross(rc);
               const bool across = !(holdOn[i] && !gStill[i]);   // a held ball: across the cord its held arc grips (11.3); an arc held still: the cord is still
               Grip g;
               g.s1 = sl.Dot(t1);
               g.s2 = across ? sl.Dot(t2) : 0.0;
               g.b = wm.Dot(t2);
               g.kPP = g.kPz = k1;
               g.kM = 0.5 / Ib;
               g.c = 0.5 * R / Ib;
               g.LP = set.mu * N * h;
               g.LM = set.rollingResistance * R * N * h;
               {
                  Grip ge = g;   // on the end rates
                  ge.s1 = g.s1 - 0.5 * sls.Dot(t1);
                  ge.s2 = across ? g.s2 - 0.5 * sls.Dot(t2) : 0.0;
                  ge.b = g.b - 0.5 * ws.Dot(t2);
                  SolveGrip(ge);
                  const double m1 = g.s1 + g.kPP * ge.P1 - g.c * ge.M, m2 = g.s2 + g.kPz * ge.P2, mb = g.b - g.c * ge.P1 + g.kM * ge.M;   // its mid-step rates
                  if (ge.P1 * m1 + ge.P2 * m2 + ge.M * mb <= 0.0)
                  {
                     g.P1 = ge.P1;
                     g.P2 = ge.P2;
                     g.M = ge.M;
                  }
                  else
                     SolveGrip(g);
               }
               gdiag.slip1[i] = g.s1;
               gdiag.slip2[i] = g.s2;
               gdiag.stick[i] = g.LP <= 0.0 || g.P1 * g.P1 + g.P2 * g.P2 < g.LP * g.LP * (1.0 - 1e-9) ? 1 : 0;
               if (holdOn[i])   // the held arc: 2 turns with the ball, 3 slides on it, 4 held still by the surface (the ball rubs), 6 turns with the ball, pressed
                  gdiag.stick[i] = gMode[i];
               const double tlim = set.mu * N * set.twistRadiusM * h;
               {
                  const double c0 = wm.Dot(n3), cE = c0 - 0.5 * ws.Dot(n3), mE = std::min(std::max(-2.0 * Ib * cE, -tlim), tlim);
                  Mtw = mE * (c0 + mE / (2.0 * Ib)) <= 0.0 ? mE : std::min(std::max(-2.0 * Ib * c0, -tlim), tlim);
               }
               Mr = g.M;
               Pn = t1 * g.P1 + t2 * g.P2;
               {
                  const double pl = Pn.Length();
                  gSlideOn[i] = g.LP > 0.0 && pl >= g.LP * (1.0 - 1e-9);
                  if (gSlideOn[i])
                     gSlideDir[i] = Pn * (1.0 / pl);
               }
               Tq = rc.Cross(Pn) + n3 * Mtw + t2 * Mr;
               // an arc turning with the ball while pressed onto a surface: the surface's friction on that turn (about the cord's direction, limit mu Ns (R + 2 r)
               // over the step), stick / slide on the ball's mid-step turn rate without this impulse (free ball: compliance 1 / (2 I))
               double si = 0.0;
               if (holdOn[i] && !gStill[i] && gNs0[i] > 0.0)
               {
                  const double ls = set.mu * gNs0[i] * (R + 2.0 * r) * h, u0 = wm.Dot(gT1[i]), uE = u0 - 0.5 * ws.Dot(gT1[i]);
                  const double sE = std::min(std::max(-2.0 * Ib * uE, -ls), ls);
                  si = sE * (u0 + sE / (2.0 * Ib)) <= 0.0 ? sE : std::min(std::max(-2.0 * Ib * u0, -ls), ls);
                  Tq = Tq + gT1[i] * si;
               }
               ch = std::max(ch, std::fabs(si - gSurfImp[i]) / (Ib / R));
               gSurfImp[i] = si;
               gdiag.gripP1[i] = g.P1 / h;
               gdiag.gripP2[i] = g.P2 / h;
               gdiag.gripRoll[i] = Mr / h;
               gdiag.gripTwist[i] = Mtw / h;
               gdiag.gripLim[i] = set.mu * N;
               gdiag.gripSurf[i] = si / h;
               gdiag.gripT1[i] = t1;
               gdiag.gripT2[i] = t2;
               gdiag.gripN[i] = n3;
            }
            if (N < 1e-9)
            {
               gdiag.stick[i] = -1;
               gdiag.gripP1[i] = gdiag.gripP2[i] = gdiag.gripRoll[i] = gdiag.gripTwist[i] = gdiag.gripLim[i] = gdiag.gripSurf[i] = 0.0;
               ch = std::max(ch, std::fabs(gSurfImp[i]) / (Ib / R));
               gSurfImp[i] = 0.0;
            }
            ch = std::max(ch, std::max((Pn - gP[i]).Length() / m, (Tq - gTq[i]).Length() / (Ib / R)));
            gP[i] = Pn;
            gTq[i] = Tq;
            gMtw[i] = Mtw;
            gMr[i] = Mr;
            gN3[i] = n3;
            gT2[i] = t2;
         }
         return ch;
      }

      double MoveScale(const double *stp, int off) const override { return gPl && std::fabs(stp[off]) > 0.05 ? 0.05 / std::fabs(stp[off]) : 1.0; }
      bool Settled(const double *stp, int off, double scale, double h) const override { return !gPl || std::fabs(stp[off]) * scale / h <= set.tolW; }
      // the plate's row (an angular momentum) in the line search's speed units: over its inertia, times an arm of 0.05 m (ASSUMED; it only weighs the rows)
      double ResNorm2(const double *Res, int off) const override
      {
         if (!gPl)
            return 0.0;
         const double s = Res[off] / gI0 * 0.05;
         return s * s;
      }
      double StopExcess(const double *x, int off) const override { return gPl ? phi + x[off] - PhiStop() : -INFINITY; }

      void Accept(const double *x, int off, double h, bool onStop, double &driveJ, double &lossJ, bool &through) override
      {
         double dq[kMaxGD];
         LocalMove(x, off, dq);
         // the grip's work over the step from the end rates (its impulses are the last UpdateGrip's)
         double fr = 0.0, tw = 0.0, ro = 0.0;
         for (int i = 0; i < nb; i++)
         {
            const V3 vm = ToBand(GV3(x[6 * i], x[6 * i + 1], x[6 * i + 2])) * (1.0 / h);
            const V3 wm = ToBand(gb[i].w + GV3(x[6 * i + 3], x[6 * i + 4], x[6 * i + 5]) * 0.5);
            const V3 rc = gN3[i] * -R;
            fr += -gP[i].Dot(vm + wm.Cross(rc));
            tw += -gMtw[i] * wm.Dot(gN3[i]);
            ro += -gMr[i] * wm.Dot(gT2[i]);
            fr += -gSurfImp[i] * wm.Dot(gT1[i]);   // the surface's friction on a turning pressed arc
            gSurfRate[i] = gSurfImp[i] * (1.0 / h);
            gPrate[i] = gP[i] * (1.0 / h);
            gTqRate[i] = gTq[i] * (1.0 / h);
         }
         double stopJ = 0.0;
         if (gPl)
         {
            const double phi1 = phi + dq[3 * nb], I1 = Inertia(phi1), Jb = 0.5 * (1.0 / gI0 + 1.0 / I1);
            double omega1 = (2.0 * dq[3 * nb] / (h * Jb) - gp0) / I1;
            gAlphaLast = (omega1 - omega) / h;
            if (onStop)   // the plunger stops on the steel coil stop: inelastic
            {
               stopJ = 0.5 * I1 * omega1 * omega1;
               omega1 = 0.0;
            }
            phi = phi1;
            omega = omega1;
         }
         iHigh = gIHigh1;
         for (int i = 0; i < nb; i++)
            gFnTaken[i] = gFn[i];
         // the balls where the group puts them (c + dc, v = 2 dc / h - v)
         for (int i = 0; i < nb; i++)
         {
            const GV3 dc(x[6 * i], x[6 * i + 1], x[6 * i + 2]);
            GroupStep::Ball bn = gb[i];
            bn.c = gb[i].c + dc;
            bn.v = dc * (2.0 / h) - gb[i].v;
            Project(bn, c[i], z[i], v[i]);
         }
         for (int i = 0; i < nb; i++)
            if (gq0.held[i])
            {
               mHold[i] = Turn(gT1[i], gBc[i], gq0.mh[i]);
               mHold[i] = mHold[i] * (1.0 / mHold[i].Length());
            }
         path3 = gPath1;
         for (int i = 0; i < nb; i++)
         {
            const int e = Find3(path3, kBall, i);
            if (e >= 0 && !holdOn[i])
               mHold[i] = path3.ct[e].m;   // joined: held at the plane it joined in
            gTwJoined[i] = e >= 0 && !holdOn[i];
            psi[i] = e >= 0 && holdOn[i] ? gPsi1[i] : psi[i];
            if (holdOn[i] && e < 0)   // left the cord: any twist still in it is let go (the friction limit has brought it to almost nothing as the push fell)
            {
               gTwDropJ += 0.5 * kTw[i] * gPsi1[i] * gPsi1[i];
               psi[i] = 0.0;
            }
            holdOn[i] = e >= 0;
         }
         Commit3(path3);
         for (int i = 0; i < nb; i++)
            if (gTwJoined[i])
               TwistJoin(i);
         gPostSlipJ = gTwDropJ;
         gTwDropJ = 0.0;
         PostSlip();
         if (onStop)
            pinned = true;
         else if (pinned && TauNow() < 0.0)
            pinned = false;
         if (geo.hasPlate)
            driver.Tick(h);
         books.steps++;
         books.hMinUsed = std::min(books.hMinUsed, h);
         books.coilJ += gCoilW;
         books.cordLossJ += gCordLoss;
         books.backstopLossJ += gBsLoss;
         books.frictionJ += fr + gHoldLoss + gPostSlipJ;
         books.postSlipJ += gPostSlipJ;
         books.twistJ += tw;
         books.rollJ += ro;
         books.plateLossJ += gPlateLoss;
         books.stopJ += stopJ;
         books.maxOverlapM = std::max(books.maxOverlapM, gMaxOv);
         driveJ = gCoilW;
         lossJ = gCordLoss + gBsLoss + fr + gHoldLoss + gPostSlipJ + tw + ro + gPlateLoss + stopJ;
         through = books.through;
      }

      double phiRest = 0.0;   // RestPhi of the settings (Configure)
   private:
      // derived
      double m = 0.08, R = 0.0135, r = 0.0027, Ib = 0.0, EA = 0.0, Lr = 1.0, cv = 0.0, kH = 0.0, hL = 0.0054, cL = 0.0, Lcur = 0.0;
      FlipperDerived der;
      // friction impulses of the last step per unit time (warm start): along the in-plane tangent, along z, twist torque, rolling torque
      double Pf[kMaxBalls] = {}, Pz[kMaxBalls] = {}, Mt[kMaxBalls] = {};
      V3 Mr[kMaxBalls];

      struct Dofs { int n = 0, phi = -1; };
      Dofs MakeDofs() const
      {
         Dofs d;
         d.n = 2 * nb;
         if (geo.hasPlate && !pinned)
            d.phi = d.n++;
         return d;
      }

      // ---- group-step scratch and helpers (the band in a group's solve) ----
      static constexpr int kMaxBs = kMaxPosts + kMaxStops + 1 + kMaxWalls;   // backstops
      struct BsZero { double d[kMaxBs] = {}; };
      BsZero bsZero[kMaxBalls];               // per ball and backstop: the depth its squash counts from (0 unless it was found inside)
      struct PostState { bool touching = false; GroupStep::V3 xi; double zero = 0.0; };
      struct PostStates { PostState s[kMaxPosts + kMaxStops]; };
      PostStates postState[kMaxBalls];         // per ball, what its group contact on each post keeps between cycles
      PostShape postShape[kMaxPosts + kMaxStops];
      int postFirst = 0, nPostShapes = 0;      // where the shapes are in the group's list
      bool bsNew[kMaxBalls] = { true, true, true, true };   // the ball has just joined: the next commit sets its zero
      static constexpr int kMaxGD = 3 * kMaxBalls + 1;                       // the band's coordinates in a group step: per ball x, y (its plane) and z (height), the plate
      const GroupStep::Ball *gb = nullptr;   // the group's balls at the step's start
      Dofs gd;                               // the planar helpers' layout (2 per ball, the plate): the backstops' gradients
      bool gPl = false, gCoil = false;
      Config gq0;
      double gL0 = 0.0, gI0 = 1.0, gp0 = 0.0, gIHigh1 = 0.0, gI2 = 0.0, gNoiseL = 0.0, gNoiseD = 0.0;
      double gDel0[kMaxBalls * kMaxBs] = {};
      // the grip: impulses over the step in the band's frame (friction, its torque with twist and rolling resistance), per unit time for the warm start, and the
      // last UpdateGrip's frame and parts
      std::shared_ptr<Part> snap;   // an energy step's starting state
      bool gSamePath = false;       // this attempt: the path keeps its elements (HoldArcs)
      V3 gP[kMaxBalls], gTq[kMaxBalls], gPrate[kMaxBalls], gTqRate[kMaxBalls], gFn[kMaxBalls], gFnTaken[kMaxBalls], gN3[kMaxBalls], gT2[kMaxBalls];   // gFnTaken: the band's push at the end of the last step taken (a failed attempt's push is no measure of anything)
      double gMtw[kMaxBalls] = {}, gMr[kMaxBalls] = {};
      double gCordLoss = 0.0, gBsLoss = 0.0, gMaxOv = 0.0, gPlateLoss = 0.0, gCoilW = 0.0;
      // the held grip across the cord (GroupStepDesign 11.3): per ball on the cord, the plane its arc lies in, fixed to the ball (state, band frame); in a step
      // the axis it turns about (the cord's direction at the ball), the ball's turn about it (bb), the arc's (bc: equal while held, less when it slides), the
      // cord's pull on that turn (a torque, N m), whether it slid, its loss
      bool holdOn[kMaxBalls] = {};
      V3 mHold[kMaxBalls];
      V3 gT1[kMaxBalls];
      double spanS[kMaxPosts] = {}, spanL[kMaxPosts] = {};   // each span's unstretched length (changed only by a slip round a post) and its length now
      double psi[kMaxBalls] = {}, kTw[kMaxBalls] = {}, cTw[kMaxBalls] = {}, gPsi1[kMaxBalls] = {}, gTwM[kMaxBalls] = {}, gSlideF[kMaxBalls] = {}, gSlideM[kMaxBalls] = {}, gRootBb[kMaxBalls] = {}, gRootA[kMaxBalls] = {}, gBc2[kMaxBalls] = {}, gTwDropJ = 0.0;
      bool gTwJoined[kMaxBalls] = {};
      double gPostSlipJ = 0.0, gNoiseE = 0.0, gFloorTop = -1e30, gFloorV = 0.0, gFloorLoss = 0.0, gFloorDel0[kMaxBalls] = {};
      mutable double linArc[kMaxPath] = {};
      static constexpr int kArcSetPasses = 8, kArcIters = 30;   // the several arcs' solve: how many times which arcs slide may change, and rounds in each
      static constexpr double kArcTolJ = 2e-13;       // a sliding arc's energy error in a step, (pull - limit) x turn: over 250 000 steps still far below the books' 1 uJ
      static constexpr double kArcMaxMoveRad = 0.2;   // the largest move of an arc's turn in one round
      double gArcDir[kMaxBalls] = {};
      double gBb[kMaxBalls] = {}, gBc[kMaxBalls] = {}, gTqB[kMaxBalls] = {}, gSlideTq[kMaxBalls] = {}, gHoldLoss = 0.0, gN0[kMaxBalls] = {};
      // the surface's grip on a held arc (GroupStepDesign 11.10): the squash force at the step's start (its limit), the two frictions on the arc (N m), the mode
      // (2 the arc turns with the ball, 3 slides on the ball (no surface), 4 held by the surface, 5 slides on both)
      double gNs0[kMaxBalls] = {}, gTauB[kMaxBalls] = {}, gTauS[kMaxBalls] = {}, gSurfImp[kMaxBalls] = {}, gSurfRate[kMaxBalls] = {};
      int gMode[kMaxBalls] = {};
      bool gStill[kMaxBalls] = {};   // this step: the arc held still by the surface it is pressed onto (its grip the stronger), the ball rubbing on it
      bool gSlid[kMaxBalls] = {};
      // the arcs' part of a step (Eval): a pure function of the arcs' turns (start positions, fixed limits), so kept: the last one worked out and the last slide's; the
      // arcs' derivative at the step's start; each arc's slide point. Reset by BeginStep.
      struct BetaMemo { bool ok = false; double bc[kMaxBalls], Lb, Vb, LbS[kMaxPosts], LbD[kMaxBalls], gbt[kMaxBalls], Fel[kMaxBalls], ge[kMaxBalls]; };
      BetaMemo gMemo, gMemoSlide;
      double gGradB0[kMaxBalls] = {}, gRoot[kMaxBalls] = {};
      bool gGradB0Ok = false, gRootOk[kMaxBalls] = {};
      // the last Eval's move and the band's stiffness there (band coordinates: H / h = -dF / dq): the push at a new move; and a friction at its limit, along which
      double gDqEv[kMaxGD] = {}, gHev[kMaxGD][kMaxGD] = {};
      int gNdEv = 0;
      double gAlphaLast = 0.0;   // the plate's angular acceleration over the last step (the next step's first guess)
      V3 gSlideDir[kMaxBalls];
      bool gSlideOn[kMaxBalls] = {};
      static V3 Turn(const V3 &ax, double ang, const V3 &v)   // v turned by ang about the unit axis ax (Rodrigues)
      {
         const double cs = std::cos(ang), sn = std::sin(ang);
         return v * cs + ax.Cross(v) * sn + ax * (ax.Dot(v) * (1.0 - cs));
      }
      // the configuration with every held ball's plane turned by its arc's turn bc (a step's end and middle positions keep the same arc turn)
      void TurnHeld(Config &q, const double *bc) const
      {
         for (int i = 0; i < nb; i++)
            if (q.held[i])
            {
               q.mh[i] = Turn(gT1[i], bc[i], gq0.mh[i]);
               q.mh[i] = q.mh[i] * (1.0 / q.mh[i].Length());
            }
      }

      // ---- the cord's path in 3D (GroupStepDesign 11.1 - 11.2) ----
      // A pulley is a circle in a plane: an own post's groove (horizontal, at the cord's height: kFixed: it holds the height, the cord leaving it at a slope bends
      // over its lip), a round object behind the band the cord slides up and down on (a stop post, the plate face's ends, a wall's ends: horizontal, at the height
      // where its runs' slopes balance: kSlide), or a ball (a great circle of the sphere of radius R + r, R + 3 r seen by the back span, in the plane through its
      // centre and its neighbours' touch points: kSphere). (A post taken as a sphere lets the cord tilt round it and slip past a ball: not a groove.) Signed radius
      // as the planar path (+ kept inside the loop, - outside). The path keeps, per pulley, the points where the cord meets (E) and leaves (X) it and the arc's
      // plane normal m. The pull on the balls and the plate is the exact derivative of this length (GradL3).
      enum PType : int { kFixed = 0, kSlide, kSphere };
      // tIn, tOut: the tangents at E and X (tOut also stands for a run squashed to no length). The inputs they were worked out from are kept with them (cOwn: this
      // pulley's centre; cIn, cOut: its neighbours'; rrIn, rrOut: their circles' radii seen in its plane; wr: the wrap they give): the same inputs give the same
      // tangents, so a pass or a solve only works out what has moved.
      struct Contact3
      {
         V3 E, X, m = V3(0.0, 0.0, 1.0), tOut, tIn;
         double zc = 0.0;
         V3 cOwn, cIn, cOut;
         double rrIn = 0.0, rrOut = 0.0, wr = 0.0;
         bool inOk = false, outOk = false, wrOk = false;
      };
      struct Path3
      {
         int n = 0;
         Elem e[kMaxPath];
         Contact3 ct[kMaxPath];
         Path3() = default;
         Path3(const Path3 &o) { *this = o; }
         Path3 &operator=(const Path3 &o)   // only the pulleys in use
         {
            n = o.n;
            for (int k = 0; k < n; k++)
            {
               e[k] = o.e[k];
               ct[k] = o.ct[k];
            }
            return *this;
         }
      };
      struct Geom3 { V3 t[kMaxPath]; double l[kMaxPath] = {}, wrap[kMaxPath] = {}, rho[kMaxPath] = {}; V3 cen[kMaxPath]; double L = 0.0; int passes = 0; bool gone[kMaxPath] = {}; };
      Path3 path3, gPath1;   // the cord's path now (group mode); the last Eval's end path
      Geom3 G3, gG1;
      // scratch of a step's evaluation (one group at a time, as the group's Jacobian): held here so that they are not built afresh for every path solve
      mutable Path3 sPb, sPm, sTr;
      mutable Geom3 sGb, sGm, sGt;
      mutable Contact3 sPrev[kMaxPath], sPrev2[kMaxPath], sLast[kMaxPath];
      bool groupMode = false;

      V2 WallEnd(int id) const { const WallSpec &wl = geo.walls[id / 2]; return id % 2 == 0 ? wl.A : wl.B; }
      void Pulley3(const Elem &e, const Config &q, int &type, V3 &cen, double &rho) const
      {
         switch (e.kind)
         {
         case kPost: type = kFixed; cen = V3(geo.posts[e.id].c.x, geo.posts[e.id].c.y, 0.0); rho = geo.posts[e.id].r; break;
         case kStop: type = kFixed; cen = V3(geo.stops[e.id].c.x, geo.stops[e.id].c.y, 0.0); rho = geo.stops[e.id].r; break;   // the cord pinned on it does not slide up or down it (user, 2026-10-30)
         case kFace: { const V2 f = FacePoint(e.id, q.phi); type = kSlide; cen = V3(f.x, f.y, 0.0); rho = r; break; }
         case kWallEnd: { const V2 f = WallEnd(e.id); type = kSlide; cen = V3(f.x, f.y, 0.0); rho = geo.walls[e.id / 2].halfThick + r; break; }
         case kBall: type = kSphere; cen = V3(q.c[e.id].x, q.c[e.id].y, q.z[e.id]); rho = -(R + r); break;
         default: type = kSphere; cen = V3(q.c[e.id].x, q.c[e.id].y, q.z[e.id]); rho = R + 3.0 * r; break;   // kBallBack
         }
      }
      // One tangent of a pulley in plane m: d the neighbour's centre in the plane, seen from where the run leaves, dr the radii's difference
      static void Tan3(const V3 &m, const V3 &d, double dr, V3 &t)
      {
         const double D2 = d.Dot(d);
         const V3 mxd = m.Cross(d);
         if (D2 > dr * dr)
         {
            const double l = std::sqrt(D2 - dr * dr), iD = 1.0 / D2;
            t = V3((d.x * l + mxd.x * dr) * iD, (d.y * l + mxd.y * dr) * iD, (d.z * l + mxd.z * dr) * iD);
         }
         else if (D2 > 0.0)   // the cord squashed between the two (a backstop carries it from here)
            t = mxd * ((dr >= 0.0 ? 1.0 : -1.0) / std::sqrt(D2));
         else
         {
            V3 pv = V3(1.0, 0.0, 0.0) - m * m.x;
            if (pv.Length() < 0.5)
               pv = V3(0.0, 1.0, 0.0) - m * m.y;
            t = m.Cross(pv * (1.0 / pv.Length())) * (dr >= 0.0 ? 1.0 : -1.0);
         }
      }
      // A pulley's touch points for the runs to its neighbours (the previous one ending at A, of centre cA and signed radius rA; the next starting at B, cB, rB),
      // from the plane in ct.m. A ball (kSphere): a great circle of its sphere, in the plane through A, B and the centre (its side kept by ct.m), or fixed (mFix: a
      // held arc); anything else a circle in the horizontal plane, a sliding one (kSlide) at the height where the straight line of the cord, unrolled, crosses it (its
      // runs' slopes balance). Each touch point is the planar common tangent with the neighbour's circle seen in this plane (its centre projected; its radius times
      // the cosine between the two planes, mA and mB their normals: a run touching both lies in both, so that is how far the other circle reaches sideways of the
      // run, its sign the side seen from this plane); defined when the two circles overlap (a ball squashing the cord onto a post: the run between them has no
      // length; a tangent from a point inside the other circle would not exist). The result: the largest squared move of E or X.
      double Touch3(int type, V3 &cen, double rho, const V3 &A, const V3 &B, const V3 &cA, double rA, const V3 &mA, const V3 &cB, double rB, const V3 &mB, Contact3 &ct,
                    const V3 *mFix = nullptr) const
      {
         const bool flat = type != kSphere;
         V3 mm = ct.m;
         if (!flat)
         {
            if (mFix)
               mm = *mFix;
            else
            {
               const V3 da = A - cen, db = B - cen, nrm = da.Cross(db);
               const double nl2 = nrm.Dot(nrm);
               if (nl2 > 1e-18 * da.Dot(da) * db.Dot(db))
               {
                  const double il = 1.0 / std::sqrt(nl2);
                  mm = nrm * (mm.Dot(nrm) < 0.0 ? -il : il);
               }
            }
         }
         const auto same = [flat](const V3 &a, const V3 &b) { return a.x == b.x && a.y == b.y && (flat || a.z == b.z); };
         const bool ownSame = same(cen, ct.cOwn) && mm.x == ct.m.x && mm.y == ct.m.y && mm.z == ct.m.z;
         const double rrI = rA * mm.Dot(mA), rrO = rB * mm.Dot(mB);
         bool moved = false;
         if (!(ct.inOk && ownSame && rrI == ct.rrIn && same(cA, ct.cIn)))
         {
            const V3 u = cA - cen;
            Tan3(mm, mm * mm.Dot(u) - u, rrI - rho, ct.tIn);
            ct.cIn = cA;
            ct.rrIn = rrI;
            ct.inOk = true;
            moved = true;
         }
         if (!(ct.outOk && ownSame && rrO == ct.rrOut && same(cB, ct.cOut)))
         {
            const V3 u = cB - cen;
            Tan3(mm, u - mm * mm.Dot(u), rho - rrO, ct.tOut);
            ct.cOut = cB;
            ct.rrOut = rrO;
            ct.outOk = true;
            moved = true;
         }
         double mv2 = 0.0;
         const V3 E0 = ct.E, X0 = ct.X;
         if (moved)
         {
            ct.m = mm;
            ct.cOwn = cen;
            ct.wrOk = false;
            const V3 tmI = ct.tIn.Cross(mm), tmO = ct.tOut.Cross(mm);
            ct.E = V3(cen.x + rho * tmI.x, cen.y + rho * tmI.y, cen.z + rho * tmI.z);
            ct.X = V3(cen.x + rho * tmO.x, cen.y + rho * tmO.y, cen.z + rho * tmO.z);
         }
         if (type == kSlide)   // the runs' plan lengths from this circle's touch points to the neighbours'
         {
            const double ax = ct.E.x - A.x, ay = ct.E.y - A.y, bx = B.x - ct.X.x, by = B.y - ct.X.y;
            const double la = std::sqrt(ax * ax + ay * ay), lb = std::sqrt(bx * bx + by * by);
            cen.z = la + lb > 0.0 ? A.z + (B.z - A.z) * la / (la + lb) : 0.5 * (A.z + B.z);
            ct.E.z = ct.X.z = cen.z;
            moved = true;
         }
         if (type != kSphere)
            ct.zc = cen.z;
         if (moved)
         {
            const V3 e0 = ct.E - E0, x0 = ct.X - X0;
            mv2 = std::max(e0.Dot(e0), x0.Dot(x0));
         }
         return mv2;
      }
      // the arc's wrap from its tangents, followed from wrapRef
      static double WrapOf(Contact3 &ct, double rho, double wrapRef)
      {
         if (!ct.wrOk)
         {
            const double wr = std::atan2(ct.tIn.Cross(ct.tOut).Dot(ct.m), ct.tIn.Dot(ct.tOut));
            ct.wr = rho < 0.0 ? -wr : wr;
            ct.wrOk = true;
         }
         const double e = wrapRef - ct.wr;
         return std::fabs(e) < kPi ? ct.wr : ct.wr + 2.0 * kPi * std::round(e / (2.0 * kPi));
      }
      // The path's touch points at q, from its present ones: block passes over the pulleys (each one's touch points exact for its neighbours') until no point
      // moves more than 2e-14 m (the length's error is second order in that), then the runs, wraps and the length. False: not settled in 400 passes.
      bool Solve3(Path3 &p, const Config &q, Geom3 &g) const
      {
         const int n = p.n;
         int type[kMaxPath];
         for (int k = 0; k < n; k++)
         {
            Pulley3(p.e[k], q, type[k], g.cen[k], g.rho[k]);
            if (type[k] == kSlide)
               g.cen[k].z = p.ct[k].zc;
         }
         bool conv = false;
         int it = 0;
         double mvPrev = 0.0;
         for (; it < 400 && !conv; it++)
         {
            const bool creepPass = it >= 8 && it % 3 == 2;
            if (creepPass)
               for (int k = 0; k < n; k++)
               {
                  sPrev2[k] = sLast[k];
                  sPrev[k] = p.ct[k];
               }
            else if (it >= 7 && it % 3 == 1)
               for (int k = 0; k < n; k++)
                  sLast[k] = p.ct[k];
            double mv2 = 0.0;
            for (int k = 0; k < n; k++)
            {
               const int km = k == 0 ? n - 1 : k - 1, kp = k + 1 == n ? 0 : k + 1;
               const V3 *mFix = p.e[k].kind == kBall && q.held[p.e[k].id] ? &q.mh[p.e[k].id] : nullptr;   // a held ball: its plane fixed
               mv2 = std::max(mv2, Touch3(type[k], g.cen[k], g.rho[k], p.ct[km].X, p.ct[kp].E, g.cen[km], g.rho[km], p.ct[km].m, g.cen[kp], g.rho[kp], p.ct[kp].m, p.ct[k], mFix));
               if (debugSolve && it >= 394)
               {
                  const Contact3 &c2 = p.ct[k];
                  fprintf(stderr, "pass %d %c%d mv %.3e E (%.9f %.9f %.9f) X (%.9f %.9f %.9f) mm\n", it, "PSFBKW"[p.e[k].kind], p.e[k].id, std::sqrt(mv2),
                          c2.E.x * 1e3, c2.E.y * 1e3, c2.E.z * 1e3, c2.X.x * 1e3, c2.X.y * 1e3, c2.X.z * 1e3);
               }
            }
            const double mv = std::sqrt(mv2);
            conv = mv <= 2e-14;
            // A pair of pulleys pressed together (a ball squashing the cord onto a post: the run between them has almost no length) makes the passes creep: each
            // moves the points by nearly the same step. Once the steps shrink by a steady ratio, the creep is taken to its limit (Aitken), then the passes go on.
            // The same for a pair that swaps between two places (the steps turn back on themselves: a gain near -1, with the same formula half way between them).
            if (!conv && creepPass && mvPrev > 0.0)
            {
               double turn = 0.0;
               for (int k = 0; k < n; k++)
                  turn += (p.ct[k].E - sPrev[k].E).Dot(sPrev[k].E - sPrev2[k].E) + (p.ct[k].X - sPrev[k].X).Dot(sPrev[k].X - sPrev2[k].X) + (p.ct[k].zc - sPrev[k].zc) * (sPrev[k].zc - sPrev2[k].zc);
               const double ratio = (turn < 0.0 ? -1.0 : 1.0) * mv / mvPrev;
               if ((ratio > 0.5 && ratio < 0.99999) || ratio < -0.5)
               {
                  const double f = ratio / (1.0 - ratio);
                  for (int k = 0; k < n; k++)
                  {
                     p.ct[k].E = p.ct[k].E + (p.ct[k].E - sPrev[k].E) * f;
                     p.ct[k].X = p.ct[k].X + (p.ct[k].X - sPrev[k].X) * f;
                     p.ct[k].zc += (p.ct[k].zc - sPrev[k].zc) * f;
                     p.ct[k].inOk = p.ct[k].outOk = false;
                  }
               }
            }
            mvPrev = mv;
         }
         g.passes = it;
         passStats.solves++;
         (solveCaller == 1 ? passStats.bySettle : solveCaller == 2 ? passStats.byGrad : solveCaller == 3 ? passStats.byBeta : passStats.byOther)++;
         passStats.passes += it;
         passStats.maxPasses = std::max(passStats.maxPasses, it);
         if (!conv)
         {
            passStats.fails++;
            if (!gSamePath)   // (an energy step takes the touch points the passes have reached: the length is measured from them, and the real energy closes the step)
               return false;
         }
         g.L = 0.0;
         for (int k = 0; k < n; k++)
         {
            const int kp1 = (k + 1) % n;
            V3 d = p.ct[kp1].E - p.ct[k].X;
            // two pulleys pressed together (a ball squashing the cord onto a post or the face): the touch points have passed each other along the line of the centres, and the run between them has no
            // length there: what is left of the cord's length is the arcs (the squash carries the force). A height difference between them stays.
            g.gone[k] = false;
            {
               const V3 cd = g.cen[kp1] - g.cen[k];
               const double a = d.Dot(cd), cc = cd.Dot(cd);
               if (a < 0.0 && cc > 0.0)
               {
                  d = d - cd * (a / cc);
                  g.gone[k] = true;
               }
            }
            const double l = d.Length();
            g.l[k] = l;
            g.t[k] = l > 1e-12 ? d * (1.0 / l) : p.ct[k].tOut;
            g.wrap[k] = WrapOf(p.ct[k], g.rho[k], p.e[k].wrapRef);
            g.L += l + std::fabs(g.rho[k]) * g.wrap[k];
         }
         return true;
      }
      // ---- the exact derivative of the path's length (GroupStepDesign 11.21) ----
      // The length, as Solve3 defines it: L = sum |E(k+1) - X(k)| + |rho(k)| wrap(k). On pulley k (plane normal m, centre c), a neighbour's centre projected on the
      // plane and its circle seen in the plane give each side's tangent tau (d the centre's projection, seen from where the run leaves; dr the radii's difference):
      //    tau = (d l + (m x d) dr) / |d|^2,  l = sqrt(|d|^2 - dr^2),  turning at  dtheta = dd . ((m x d) - d dr / l) / |d|^2 + ddr / l,
      //    and E (or X) = c + rho tau x m moves by  dc + rho dtheta tau + rho tau x dm;  the wrap's change is sign(rho) (dtheta_out - dtheta_in).
      // Added up over the path the large dtheta terms (1 / l as a run closes) cancel down to their remainders, rho (1 - t . tau) = rho |t - tau|^2 / 2, t the real run:
      //    dL = sum (t_in - t_out) . dc + rho [dtheta_out (1 - t_out . tau_out) - dtheta_in (1 - t_in . tau_in) + t_in . (tau_in x dm) - t_out . (tau_out x dm)]
      // (t_in - t_out is the classic pull; the rest is the run leaving a groove at a slope bending over its lip, zero for a level run.) The only feedback between the
      // pulleys is through a sliding pulley's height (where its runs' slopes balance) and a free ball's plane (through its centre and its neighbours' touch points):
      // those few tangents are the unknowns of one small linear solve, the linearised fixed point of Solve3's passes. No path solve is needed.
      static constexpr int kMaxLinY = 64;
      struct Lin3
      {
         struct Side { V3 tau, g; double c, s, w; };   // the tangent; dtheta = dd . g + ddr c; m . (neighbour - c); rho (1 - t . tau)
         struct Un { int k, at; bool slide, frozen; double f, sum, la, lb, dz, nl, sg; V3 ua, ub, da, db; };   // an unknown: a sliding pulley's height, a free ball's plane
         Side in[kMaxPath], out[kMaxPath];
         Un un[kMaxPath];
         int nun, ny;
         double A[kMaxLinY][kMaxLinY];   // I - G, factored
         int piv[kMaxLinY];
      };
      mutable Lin3 lin;
      mutable V3 linDc[kMaxPath], linDm[kMaxPath], linDE[kMaxPath], linDX[kMaxPath];
      // one side of pulley k, seen from its neighbour nbk
      void LinSide(const Path3 &p, const Geom3 &g, int k, int nbk, bool isIn, const V3 &tRun, Lin3::Side &S) const
      {
         const V3 &mk = p.ct[k].m;
         const V3 u = g.cen[nbk] - g.cen[k];
         S.s = mk.Dot(u);
         const V3 pu = u - mk * S.s, d = isIn ? pu * -1.0 : pu;
         const double rr = g.rho[nbk] * mk.Dot(p.ct[nbk].m), dr = isIn ? rr - g.rho[k] : g.rho[k] - rr, D2 = d.Dot(d);
         const V3 mxd = mk.Cross(d);
         if (D2 <= 1e-30)
         {
            S.tau = tRun;
            S.g = V3();
            S.c = 0.0;
         }
         else if (D2 > dr * dr)
         {
            const double l = std::sqrt(D2 - dr * dr);
            S.tau = (d * l + mxd * dr) * (1.0 / D2);
            S.g = (mxd - d * (dr / l)) * (1.0 / D2);
            S.c = 1.0 / l;
         }
         else   // squashed together: the run has no length, the tangent square to the centres' line
         {
            S.tau = mxd * ((dr >= 0.0 ? 1.0 : -1.0) / std::sqrt(D2));
            S.g = mxd * (1.0 / D2);
            S.c = 0.0;
         }
         const V3 e = tRun - S.tau;
         S.w = g.rho[k] * 0.5 * e.Dot(e);
      }
      // the length's change, and the touch points' (dE, dX), for centres moving dc and planes turning dm
      double LinSweep(const Path3 &p, const Geom3 &g) const
      {
         double dL = 0.0;
         const int n = p.n;
         for (int k = 0; k < n; k++)
         {
            const int km = k == 0 ? n - 1 : k - 1, kp = k + 1 == n ? 0 : k + 1;
            const Lin3::Side &I = lin.in[k], &O = lin.out[k];
            const V3 &mk = p.ct[k].m, &dm = linDm[k], &dc = linDc[k];
            const double rho = g.rho[k];
            const V3 wI = linDc[km] - dc - dm * I.s, wO = linDc[kp] - dc - dm * O.s;
            const double thI = -wI.Dot(I.g) + g.rho[km] * (dm.Dot(p.ct[km].m) + mk.Dot(linDm[km])) * I.c;
            const double thO = wO.Dot(O.g) - g.rho[kp] * (dm.Dot(p.ct[kp].m) + mk.Dot(linDm[kp])) * O.c;
            const V3 tmI = I.tau.Cross(dm) * rho, tmO = O.tau.Cross(dm) * rho;
            linDE[k] = dc + I.tau * (rho * thI) + tmI;
            linDX[k] = dc + O.tau * (rho * thO) + tmO;
            const V3 &tin = g.t[km], &tout = g.t[k];
            dL += (tin - tout).Dot(dc) + thO * O.w - thI * I.w + tin.Dot(tmI) - tout.Dot(tmO);
            linArc[k] = rho * (thO - thI);
         }
         for (int k = 0; k < n; k++)   // (a run with no length changes nothing)
            if (g.gone[k] && g.l[k] <= 1e-12)
               dL -= g.t[k].Dot(linDE[k + 1 == n ? 0 : k + 1] - linDX[k]);
         return dL;
      }
      // what the unknowns come out as from the touch points' changes (the passes' map, linearised)
      void LinOut(const Path3 &p, double *yo) const
      {
         const int n = p.n;
         for (int u = 0; u < lin.nun; u++)
         {
            const Lin3::Un &U = lin.un[u];
            const int km = U.k == 0 ? n - 1 : U.k - 1, kp = U.k + 1 == n ? 0 : U.k + 1;
            const V3 &dA = linDX[km], &dB = linDE[kp];
            if (U.slide)
            {
               double dz = dA.z * (1.0 - U.f) + dB.z * U.f;
               if (U.sum > 0.0)
               {
                  const V3 ea = linDE[U.k] - dA, eb = dB - linDX[U.k];
                  dz += U.dz * (U.ua.Dot(ea) * U.lb - U.la * U.ub.Dot(eb)) / (U.sum * U.sum);
               }
               yo[U.at] = dz;
            }
            else
            {
               V3 dw;
               if (!U.frozen)
               {
                  const V3 &mk = p.ct[U.k].m, &dc = linDc[U.k];
                  dw = (dA - dc).Cross(U.db) + U.da.Cross(dB - dc);
                  dw = (dw - mk * mk.Dot(dw)) * (U.sg / U.nl);
               }
               yo[U.at] = dw.x;
               yo[U.at + 1] = dw.y;
               yo[U.at + 2] = dw.z;
            }
         }
      }
      void LinClear(int n) const
      {
         for (int k = 0; k < n; k++)
            linDc[k] = linDm[k] = V3();
      }
      void LinAddY(const double *y) const
      {
         for (int u = 0; u < lin.nun; u++)
         {
            const Lin3::Un &U = lin.un[u];
            if (U.slide)
               linDc[U.k].z += y[U.at];
            else
               linDm[U.k] = linDm[U.k] + V3(y[U.at], y[U.at + 1], y[U.at + 2]);
         }
      }
      // the path p, solved at q, linearised: every side, and the unknowns' small solve factored
      bool Linearize3(const Path3 &p, const Geom3 &g, const Config &q) const
      {
         const int n = p.n;
         lin.nun = lin.ny = 0;
         for (int k = 0; k < n; k++)
         {
            const int km = k == 0 ? n - 1 : k - 1, kp = k + 1 == n ? 0 : k + 1;
            LinSide(p, g, k, km, true, g.t[km], lin.in[k]);
            LinSide(p, g, k, kp, false, g.t[k], lin.out[k]);
            const Elem &e = p.e[k];
            const bool slide = e.kind == kStop || e.kind == kFace || e.kind == kWallEnd;
            const bool freeBall = e.kind == kBallBack || (e.kind == kBall && !q.held[e.id]);
            if (!slide && !freeBall)
               continue;
            Lin3::Un &U = lin.un[lin.nun++];
            U.k = k;
            U.at = lin.ny;
            U.slide = slide;
            lin.ny += slide ? 1 : 3;
            const V3 &A = p.ct[km].X, &B = p.ct[kp].E;
            if (slide)
            {
               const V3 &E = p.ct[k].E, &X = p.ct[k].X;
               const V3 ea(E.x - A.x, E.y - A.y, 0.0), eb(B.x - X.x, B.y - X.y, 0.0);
               U.la = ea.Length();
               U.lb = eb.Length();
               U.sum = U.la + U.lb;
               U.f = U.sum > 0.0 ? U.la / U.sum : 0.5;
               U.dz = B.z - A.z;
               U.ua = U.la > 1e-15 ? ea * (1.0 / U.la) : V3();
               U.ub = U.lb > 1e-15 ? eb * (1.0 / U.lb) : V3();
            }
            else
            {
               U.da = A - g.cen[k];
               U.db = B - g.cen[k];
               const V3 nrm = U.da.Cross(U.db);
               U.nl = nrm.Length();
               U.frozen = !(U.nl > 1e-9 * U.da.Length() * U.db.Length());   // along one line: the plane stays where it was (TouchSphere)
               U.sg = nrm.Dot(p.ct[k].m) >= 0.0 ? 1.0 : -1.0;
            }
         }
         const int ny = lin.ny;
         if (ny == 0)
            return true;
         if (ny > kMaxLinY)
            return false;
         double y[kMaxLinY], yo[kMaxLinY];
         for (int a = 0; a < ny; a++)
            y[a] = 0.0;
         for (int c = 0; c < ny; c++)
         {
            LinClear(n);
            y[c] = 1.0;
            LinAddY(y);
            y[c] = 0.0;
            LinSweep(p, g);
            LinOut(p, yo);
            for (int r = 0; r < ny; r++)
               lin.A[r][c] = (r == c ? 1.0 : 0.0) - yo[r];
         }
         for (int c = 0; c < ny; c++)   // LU with row pivots
         {
            int pv = c;
            for (int r = c + 1; r < ny; r++)
               if (std::fabs(lin.A[r][c]) > std::fabs(lin.A[pv][c]))
                  pv = r;
            if (std::fabs(lin.A[pv][c]) < 1e-300)
               return false;
            lin.piv[c] = pv;
            if (pv != c)
               for (int a = 0; a < ny; a++)
                  std::swap(lin.A[c][a], lin.A[pv][a]);
            for (int r = c + 1; r < ny; r++)
            {
               const double f = lin.A[r][c] / lin.A[c][c];
               lin.A[r][c] = f;
               for (int a = c + 1; a < ny; a++)
                  lin.A[r][a] -= f * lin.A[c][a];
            }
         }
         return true;
      }
      // the length's change for the seed now in linDc and linDm (the unknowns solved for it)
      double LinDeriv(const Path3 &p, const Geom3 &g) const
      {
         const int ny = lin.ny;
         if (ny == 0)
            return LinSweep(p, g);
         double y[kMaxLinY];
         LinSweep(p, g);
         LinOut(p, y);
         for (int c = 0; c < ny; c++)   // (the factor's rows carry all the swaps: the right-hand side takes them first)
            std::swap(y[c], y[lin.piv[c]]);
         for (int c = 0; c < ny; c++)
            for (int r = c + 1; r < ny; r++)
               y[r] -= lin.A[r][c] * y[c];
         for (int c = ny - 1; c >= 0; c--)
         {
            for (int a = c + 1; a < ny; a++)
               y[c] -= lin.A[c][a] * y[a];
            y[c] /= lin.A[c][c];
         }
         LinAddY(y);
         return LinSweep(p, g);
      }
      // dL / d(arc turn) of the held balls: the path p solved at q, each arc's plane turning about its axis
      bool GradBeta(const Path3 &p, const Geom3 &g, const Config &q, double *gb2, const double *wSpan = nullptr, const double *wF = nullptr) const
      {
         for (int i = 0; i < nb; i++)
            gb2[i] = 0.0;
         if (!Linearize3(p, g, q))
            return false;
         SpanMap sm;
         if (wSpan)
            MakeSpanMap(p, g, sm);
         for (int i = 0; i < nb; i++)
         {
            const int e = q.held[i] ? Find3(p, kBall, i) : -1;
            if (e < 0)
               continue;
            LinClear(p.n);
            linDm[e] = gT1[i].Cross(q.mh[i]);
            gb2[i] = LinDeriv(p, g);
            if (wSpan)
            {
               double one[kMaxPosts][kMaxGD];
               for (int j = 0; j < geo.nPosts; j++)
                  one[j][0] = 0.0;
               LinSpans(p, g, sm, one, 0);
               gb2[i] = 0.0;
               for (int j = 0; j < geo.nPosts; j++)
                  gb2[i] += wSpan[j] * one[j][0];
               if (wF)
                  gb2[i] += LinFloor(p, g, wF);
            }
         }
         return true;
      }
      // the same by central differences of the path solve, 1e-3 rad (debugGrad: the check)
      bool GradBetaNum(const Path3 &p, const Config &q, const double *bc, double *gb2) const
      {
         for (int i = 0; i < nb; i++)
         {
            gb2[i] = 0.0;
            if (!q.held[i])
               continue;
            const double d = 1e-3;
            double Lp = 0.0, Lm = 0.0;
            for (int sgn = -1; sgn <= 1; sgn += 2)
            {
               double b2[kMaxBalls];
               for (int k = 0; k < nb; k++)
                  b2[k] = bc[k];
               b2[i] += sgn * d;
               Config qq = q;
               TurnHeld(qq, b2);
               Path3 pp = p;
               Geom3 gg;
               if (!Solve3(pp, qq, gg))
                  return false;
               (sgn < 0 ? Lm : Lp) = gg.L;
            }
            gb2[i] = (Lp - Lm) / (2.0 * d);
         }
         return true;
      }
      static int Find3(const Path3 &p, int kind, int id)
      {
         for (int k = 0; k < p.n; k++)
            if (p.e[k].kind == kind && p.e[k].id == id)
               return k;
         return -1;
      }
      // The path at q, from p: a pulley leaves where its wrap goes below zero; one joins only if, put in where a run cuts into it from its own side, its wrap comes
      // out positive (the same test both ways: no pulley can join and leave in turn). A ball has no fixed side: it may meet a run from any direction, and its arc's
      // plane starts with its centre on the bend's inside (GroupStepDesign 11.17). The back span: a run beyond the ball's own arc, seen from the ball. lg: the last
      // actions, for the trace.
      bool Settle3(Path3 &p, const Config &q, Geom3 &g, SettleLog *lg = nullptr) const
      {
         if (lg)
            lg->n = 0;
         const auto Log = [&](char act, const Elem &e, int at, double val, int pn) {
            if (!lg)
               return;
            SettleLog::Entry &en = lg->e[lg->n % SettleLog::kN];
            en.act = act;
            en.kind = e.kind;
            en.id = e.id;
            en.at = at;
            en.val = val;
            en.pathN = pn;
            lg->n++;
         };
         bool rejected[8 * 16] = {};
         for (int it = 0; it < 4 * kMaxPath; it++)
         {
            if (!Solve3(p, q, g))
               return false;
            int worst = -1;
            double wmin = -1e-12;
            for (int k = 0; k < p.n; k++)
               if (g.wrap[k] < wmin)
               {
                  wmin = g.wrap[k];
                  worst = k;
               }
            if (worst >= 0 && p.n > 2)
            {
               Log('-', p.e[worst], worst, wmin, p.n);
               for (int k = worst; k + 1 < p.n; k++)
               {
                  p.e[k] = p.e[k + 1];
                  p.ct[k] = p.ct[k + 1];
               }
               p.n--;
               continue;
            }
            // the deepest cut of a run into a pulley not in the path, from its own side, inside the run (not at its ends)
            int bestSeg = -1;
            Elem bestE;
            double bestDepth = 1e-12;
            const auto Try = [&](int kind, int id) {
               if (Find3(p, kind, id) >= 0 || rejected[kind * 16 + id])
                  return;
               Elem e { kind, id, 0.0 };
               int type;
               V3 cen;
               double rho;
               Pulley3(e, q, type, cen, rho);
               const int own = kind == kBallBack ? Find3(p, kBall, id) : -1;
               for (int k = 0; k < p.n; k++)
               {
                  if (own >= 0 && (k == own || (k + 1) % p.n == own))   // the runs to its own arc
                     continue;
                  const V3 P0 = p.ct[k].X;
                  V3 t = g.t[k];
                  double l = g.l[k];
                  V3 dc = cen - P0;
                  if (type != kSphere)   // a vertical circle: the run's plan view
                  {
                     const double th = std::sqrt(t.x * t.x + t.y * t.y);
                     if (th < 1e-12)
                        continue;
                     l *= th;
                     t = V3(t.x / th, t.y / th, 0.0);
                     dc.z = 0.0;
                  }
                  if (type != kSphere && (t.x * dc.y - t.y * dc.x) * rho <= 0.0)   // an upright surface behind the band: not on its side of this run
                     continue;
                  const double u = dc.Dot(t);
                  if (u <= 0.0 || u >= l)
                     continue;
                  if (own >= 0 && (t * u - dc).Dot((p.ct[own].E - cen) + (p.ct[own].X - cen)) <= 0.0)   // not beyond its own arc
                     continue;
                  const double depth = std::fabs(rho) - (dc - t * u).Length();
                  if (depth > bestDepth)
                  {
                     bestDepth = depth;
                     bestSeg = k;
                     bestE = e;
                  }
               }
            };
            for (int j = 0; j < geo.nPosts; j++)
               Try(kPost, j);
            for (int j = 0; j < geo.nStops; j++)
               Try(kStop, j);
            if (geo.hasPlate)
            {
               Try(kFace, 0);
               Try(kFace, 1);
            }
            for (int j = 0; j < 2 * geo.nWalls; j++)
               Try(kWallEnd, j);
            for (int i = 0; i < nb; i++)
            {
               Try(kBall, i);
               if (Find3(p, kBall, i) >= 0)
                  Try(kBallBack, i);
            }
            if (bestSeg < 0)
               return true;
            if (p.n >= kMaxPath)
               return false;
            Path3 &tr = sTr;
            tr = p;
            const int at = bestSeg + 1;
            for (int k = tr.n; k > at; k--)
            {
               tr.e[k] = tr.e[k - 1];
               tr.ct[k] = tr.ct[k - 1];
            }
            tr.e[at] = bestE;
            tr.n++;
            {
               int type;
               V3 cen;
               double rho;
               Pulley3(bestE, q, type, cen, rho);
               tr.ct[at] = Contact3();
               tr.ct[at].zc = 0.5 * (tr.ct[at - 1].X.z + tr.ct[(at + 1) % tr.n].E.z);
               int ta, tb;
               V3 ca, cb;
               double ra, rb;
               Pulley3(tr.e[at - 1], q, ta, ca, ra);
               Pulley3(tr.e[(at + 1) % tr.n], q, tb, cb, rb);
               if (ta == kSlide)
                  ca.z = tr.ct[at - 1].zc;
               if (tb == kSlide)
                  cb.z = tr.ct[(at + 1) % tr.n].zc;
               const V3 *mFix = bestE.kind == kBall && q.held[bestE.id] ? &q.mh[bestE.id] : nullptr;
               V3 m0(0.0, 0.0, 1.0);
               if (type == kSphere)   // the plane through its neighbours' points and its centre, the centre on the bend's inside (right of travel for rho < 0)
               {
                  const V3 nrm = (tr.ct[at - 1].X - cen).Cross(tr.ct[(at + 1) % tr.n].E - cen);
                  const double nl = nrm.Length();
                  if (nl > 0.0)
                     m0 = nrm * ((rho < 0.0 ? -1.0 : 1.0) / nl);
               }
               tr.ct[at].m = m0;
               Touch3(type, cen, rho, tr.ct[at - 1].X, tr.ct[(at + 1) % tr.n].E, ca, ra, tr.ct[at - 1].m, cb, rb, tr.ct[(at + 1) % tr.n].m, tr.ct[at], mFix);
            }
            if (Solve3(tr, q, sGt) && sGt.wrap[at] > 0.0)
            {
               Log('+', bestE, at, bestDepth, p.n);
               p = tr;
               continue;
            }
            Log('x', bestE, at, bestDepth, p.n);   // put in, it would not wrap: left out (its squash, if any, carries it)
            rejected[bestE.kind * 16 + bestE.id] = true;
         }
         return false;
      }
      // The path at q in an energy step: a pulley the cord has come off leaves, nothing joins (which of them are on it can go round in a ring: that waits for the next step's solve)
      bool LeaveOnly(Path3 &p, const Config &q, Geom3 &g) const
      {
         for (;;)
         {
            if (!Solve3(p, q, g))
               return false;
            int worst = -1;
            double wmin = -1e-12;
            for (int k = 0; k < p.n; k++)
               if (g.wrap[k] < wmin)
               {
                  wmin = g.wrap[k];
                  worst = k;
               }
            if (worst < 0 || p.n <= 2)
               return true;
            for (int k = worst; k + 1 < p.n; k++)
            {
               p.e[k] = p.e[k + 1];
               p.ct[k] = p.ct[k + 1];
            }
            p.n--;
         }
      }
      // Each span between two own posts has its own unstretched length and so its own tension (GroupStepDesign 11.24). A span is the runs and arcs from post j to post j + 1, and
      // half of each of the two posts' arcs. A post the cord has lifted off (no wrap, so no grip) still marks where its two spans meet: the nearest point of the run passing it.
      struct SpanMap
      {
         int arcSpan[kMaxPath];                 // an element that is not a post: its span (a post: -1)
         int nPart[kMaxPath];                   // run k (from element k to k + 1): its parts by span, in order, as fractions of its length
         int partSpan[kMaxPath][kMaxPosts + 1];
         double partFrac[kMaxPath][kMaxPosts + 1];
      };
      mutable SpanMap sMapB, sMapM, sMap1;
      mutable double sLsM[kMaxPosts] = {}, sLs1[kMaxPosts] = {}, sGs[kMaxPosts][kMaxGD] = {};
      void MakeSpanMap(const Path3 &p, const Geom3 &g, SpanMap &M) const
      {
         const int n = p.n, nS = geo.nPosts;
         int k0 = -1;
         for (int k = 0; k < n && k0 < 0; k++)
            if (p.e[k].kind == kPost)
               k0 = k;
         if (k0 < 0)   // no post on the path: all one span
         {
            for (int k = 0; k < n; k++)
            {
               M.arcSpan[k] = 0;
               M.nPart[k] = 1;
               M.partSpan[k][0] = 0;
               M.partFrac[k][0] = 1.0;
            }
            return;
         }
         int at = 0;
         while (at < n)
         {
            // the stretch from this post to the next post on the path
            const int ka = (k0 + at) % n, a = p.e[ka].id;
            int len = 1;
            while (len < n && p.e[(ka + len) % n].kind != kPost)
               len++;
            const int b = len < n ? p.e[(ka + len) % n].id : a;
            int missing = ((b - a - 1) % nS + nS) % nS;
            if (len >= n && nS > 1)
               missing = nS - 1;   // one post on the path: all the others are off it
            int mRun[kMaxPosts];
            double mFrac[kMaxPosts];
            for (int m = 0; m < missing; m++)
            {
               const Circle &c0 = geo.posts[(a + 1 + m) % nS];
               double best = 1e30;
               int bestRun = ka;
               double bestF = 0.0;
               for (int i = 0; i < len; i++)
               {
                  const int k = (ka + i) % n;
                  const double l = g.l[k];
                  const V3 d = V3(c0.c.x, c0.c.y, 0.0) - p.ct[k].X;
                  const double u = l > 0.0 ? std::min(std::max(d.x * g.t[k].x + d.y * g.t[k].y, 0.0), l) : 0.0;
                  const double dx = d.x - g.t[k].x * u, dy = d.y - g.t[k].y * u, dd = dx * dx + dy * dy;
                  if (dd < best)
                  {
                     best = dd;
                     bestRun = i;
                     bestF = l > 0.0 ? u / l : 0.0;
                  }
               }
               mRun[m] = bestRun;
               mFrac[m] = bestF;
            }
            int cur = a, mNext = 0;
            for (int i = 0; i < len; i++)
            {
               const int k = (ka + i) % n;
               M.arcSpan[k] = p.e[k].kind == kPost ? -1 : cur;
               int np = 0;
               double from = 0.0;
               while (mNext < missing && mRun[mNext] == i)
               {
                  const double f = std::max(mFrac[mNext], from);
                  M.partSpan[k][np] = cur;
                  M.partFrac[k][np] = f - from;
                  np++;
                  from = f;
                  cur = (a + 1 + mNext) % nS;
                  mNext++;
               }
               M.partSpan[k][np] = cur;
               M.partFrac[k][np] = 1.0 - from;
               M.nPart[k] = np + 1;
            }
            at += len;
         }
      }
      void SpanLengths(const Path3 &p, const Geom3 &g, const SpanMap &M, double *Ls) const
      {
         const int nS = geo.nPosts;
         for (int j = 0; j < nS; j++)
            Ls[j] = 0.0;
         for (int k = 0; k < p.n; k++)
         {
            const double arc = std::fabs(g.rho[k]) * g.wrap[k];
            if (p.e[k].kind == kPost)
            {
               Ls[(p.e[k].id + nS - 1) % nS] += 0.5 * arc;
               Ls[p.e[k].id] += 0.5 * arc;
            }
            else
               Ls[M.arcSpan[k]] += arc;
            for (int q2 = 0; q2 < M.nPart[k]; q2++)
               Ls[M.partSpan[k][q2]] += M.partFrac[k][q2] * g.l[k];
         }
      }
      // each span's share of the length's change for the seed just solved (LinDeriv): into gs[span][slot]
      void LinSpans(const Path3 &p, const Geom3 &g, const SpanMap &M, double (*gs)[kMaxGD], int slot) const
      {
         const int nS = geo.nPosts, n = p.n;
         for (int k = 0; k < n; k++)
         {
            const int kp = k + 1 == n ? 0 : k + 1;
            if (p.e[k].kind == kPost)
            {
               gs[(p.e[k].id + nS - 1) % nS][slot] += 0.5 * linArc[k];
               gs[p.e[k].id][slot] += 0.5 * linArc[k];
            }
            else
               gs[M.arcSpan[k]][slot] += linArc[k];
            const double dl = g.gone[k] && g.l[k] <= 1e-12 ? 0.0 : g.t[k].Dot(linDE[kp] - linDX[k]);
            for (int q2 = 0; q2 < M.nPart[k]; q2++)
               gs[M.partSpan[k][q2]][slot] += M.partFrac[k][q2] * dl;
         }
      }
      // The cord under the ball against the playfield (GroupStepDesign 11.24): the squash is the cord's radius less the height of the lowest point of its arc on the ball, the
      // bottom of the arc's great circle if the arc reaches it, else the lower of its ends. With the ball straight above that point this is a wall's law, R + 2 r less the centre's height.
      struct ArcLowInfo { double z = 0.0, len = 0.0; int how = 0; };   // how: 0 at E, 1 at X, 2 the circle's bottom
      ArcLowInfo ArcLow(const Path3 &p, const Geom3 &g, int e) const
      {
         ArcLowInfo a;
         const Contact3 &ct = p.ct[e];
         a.how = ct.E.z <= ct.X.z ? 0 : 1;
         a.z = std::min(ct.E.z, ct.X.z);
         const V3 &m3 = ct.m;
         const double len2 = 1.0 - m3.z * m3.z;
         if (len2 > 1e-12)
         {
            const double len = std::sqrt(len2), rr = std::fabs(g.rho[e]);
            const V3 dn = (V3(0.0, 0.0, 1.0) - m3 * m3.z) * (-1.0 / len);
            const V3 mid = (ct.E - g.cen[e]) + (ct.X - g.cen[e]);
            const double ml = mid.Length();
            if (ml > 0.0 && dn.Dot(mid) / ml >= std::cos(0.5 * g.wrap[e]))
            {
               a.how = 2;
               a.z = g.cen[e].z - rr * len;
               a.len = len;
            }
         }
         return a;
      }
      void FloorState(const Path3 &p, const Geom3 &g, double *del) const
      {
         for (int i = 0; i < nb; i++)
         {
            const int e = Find3(p, kBall, i);
            del[i] = e >= 0 ? gFloorTop - ArcLow(p, g, e).z : -1.0;
         }
      }
      double FloorVOf(const double *del) const
      {
         double e = 0.0;
         for (int i = 0; i < nb; i++)
            e += SquashV(del[i]);
         return e;
      }
      // the playfield's share of the energy's change for the seed just solved (LinDeriv): each ball's weight (its squash force, with damping) on its arc's lowest point
      double LinFloor(const Path3 &p, const Geom3 &g, const double *wF) const
      {
         if (!wF)
            return 0.0;
         double d = 0.0;
         for (int i = 0; i < nb; i++)
         {
            if (wF[i] <= 0.0)
               continue;
            const int e = Find3(p, kBall, i);
            if (e < 0)
               continue;
            const ArcLowInfo a = ArcLow(p, g, e);
            const double dz = a.how == 0 ? linDE[e].z : a.how == 1 ? linDX[e].z : linDc[e].z + std::fabs(g.rho[e]) * p.ct[e].m.z * linDm[e].z / a.len;
            d -= wF[i] * dz;
         }
         return d;
      }
      // the cord's twist (GroupStepDesign 11.28): each ball's own, from when it joins the cord; its stiffness fixed for the contact so that its energy is exact
      double TwistV() const
      {
         double e = 0.0;
         for (int i = 0; i < nb; i++)
            e += 0.5 * kTw[i] * psi[i] * psi[i];
         return e;
      }
      void TwistJoin(int i)
      {
         const int e = Find3(path3, kBall, i), pn = path3.n;
         psi[i] = 0.0;
         if (e < 0)
            return;
         // to the next ball, or part way into the next post (user, 2026-10-30: a twist travels partway into a post, not past it; half its wrap ASSUMED)
         const int em = (e + pn - 1) % pn, ep = (e + 1) % pn;
         const double a1 = path3.e[em].kind == kPost ? 0.5 * std::fabs(G3.rho[em]) * std::max(G3.wrap[em], 0.0) : 0.0, a2 = path3.e[ep].kind == kPost ? 0.5 * std::fabs(G3.rho[ep]) * std::max(G3.wrap[ep], 0.0) : 0.0;
         const double l1 = std::max(G3.l[em] + a1, 1e-3), l2 = std::max(G3.l[e] + a2, 1e-3), G = set.eMpa * 1e6 / 3.0, J = 0.5 * kPi * r * r * r * r;
         kTw[i] = set.cordTwistScale * G * J * (1.0 / l1 + 1.0 / l2);
         cTw[i] = kTw[i] * set.tanDelta / (2.0 * kPi * std::max(set.refHz, 1e-3));
      }
      double SpanT(int j, double L) const { return L > spanS[j] ? EA * (L - spanS[j]) / spanS[j] : 0.0; }
      double SpanV(const double *Ls) const
      {
         double e = 0.0;
         for (int j = 0; j < geo.nPosts; j++)
         {
            const double st = Ls[j] - spanS[j];
            if (st > 0.0)
               e += 0.5 * EA / spanS[j] * st * st;
         }
         return e;
      }
      double SpanCv(int j) const { return cv * Lr / spanS[j]; }
      // Capstan friction at the posts: past T_high / T_low = e^(mu theta) the cord slips round the post until the ratio is at that limit: unstretched cord moves from the slack span
      // to the tight one, the path unchanged, and the energy it gives up is the friction's (dV/d(moved) = (T_lo + T_lo^2 / 2 EA) - (T_hi + T_hi^2 / 2 EA) < 0: never a gain). Posts in passes,
      // since a slip changes the next post's tensions. No path solve, so no cost to the step.
      void PostSlip()
      {
         const int nS = geo.nPosts;
         if (nS < 2)
            return;
         double wr[kMaxPosts] = {};
         for (int k = 0; k < path3.n; k++)
            if (path3.e[k].kind == kPost)
               wr[path3.e[k].id] = std::max(G3.wrap[k], 0.0);
         const double V0 = SpanV(spanL);
         for (int pass = 0; pass < 40; pass++)
         {
            bool any = false;
            for (int j = 0; j < nS; j++)
            {
               const int a = (j + nS - 1) % nS, ta = SpanT(a, spanL[a]) >= SpanT(j, spanL[j]) ? a : j, lo = ta == a ? j : a;
               const double K = std::exp(std::min(set.postMu * wr[j], 40.0));
               if (SpanT(ta, spanL[ta]) <= K * SpanT(lo, spanL[lo]) * (1.0 + 1e-12) + 1e-9)
                  continue;
               const double s0h = spanS[ta], s0l = spanS[lo];
               double dA = 0.0, dB = 0.999 * s0l;
               for (int it = 0; it < 60; it++)
               {
                  const double d = 0.5 * (dA + dB);
                  spanS[ta] = s0h + d;
                  spanS[lo] = s0l - d;
                  (SpanT(ta, spanL[ta]) > K * SpanT(lo, spanL[lo]) ? dA : dB) = d;
               }
               spanS[ta] = s0h + dB;
               spanS[lo] = s0l - dB;
               any = true;
            }
            if (!any)
               break;
         }
         gPostSlipJ += V0 - SpanV(spanL);
      }
      // the path now becomes the step's start (group mode): its geometry at the balls' positions now, the wraps it is followed from, the books that read it
      void Commit3(Path3 &p)
      {
         const Config q = Now();
         Geom3 g;
         if (Solve3(p, q, g))
            G3 = g;
         for (int k = 0; k < p.n; k++)
         {
            p.e[k].wrapRef = G3.wrap[k];
            if (p.e[k].kind == kBallBack)
               books.backTouched = true;
         }
         Lcur = G3.L;
         for (int i = 0; i < nb; i++)   // (only ever down while the squash is 0: nothing stored changes)
         {
            for (int j = 0; j < NumBackstops(); j++)
            {
               const double raw = std::max(OverlapRaw(i, j, q, Dofs(), nullptr), 0.0);
               bsZero[i].d[j] = bsNew[i] || raw - bsZero[i].d[j] >= hL * (1.0 - 1e-9) ? raw : std::min(bsZero[i].d[j], raw);
            }
            bsNew[i] = false;
         }
         {
            SpanMap sm;
            MakeSpanMap(p, G3, sm);
            SpanLengths(p, G3, sm, spanL);
            double del[kMaxBalls];
            FloorState(p, G3, del);
            gFloorV = FloorVOf(del);
            for (int i = 0; i < nb; i++)
               if (del[i] > 0.0)
               {
                  books.cordOnPlayfield = true;
                  books.floorMaxM = std::max(books.floorMaxM, del[i]);
               }
         }
         books.maxStretch = std::max(books.maxStretch, Stretch());
         // through: a ball not on the cord that cuts into one of its runs deeper than its squash on a backstop there accounts for
         for (int i = 0; i < nb; i++)
         {
            if (Find3(p, kBall, i) >= 0)
               continue;
            double ov = 0.0;
            for (int j = 0; j < NumBackstops(); j++)
               ov = std::max(ov, OverlapRaw(i, j, q, Dofs(), nullptr));
            const V3 C(c[i].x, c[i].y, z[i]);
            for (int k = 0; k < p.n; k++)
            {
               const V3 dc = C - p.ct[k].X;
               const double u = std::min(std::max(dc.Dot(G3.t[k]), 0.0), G3.l[k]);
               if (R + r - (dc - G3.t[k] * u).Length() > ov + 1e-6)
                  books.through = true;
            }
         }
      }
      // dL/dq at the path p, solved at q (band coordinates: per ball x, y, z, then the plate when pl): the exact derivative of its length (above). Not only
      // t_in - t_out: a run leaving a groove at a slope bends over its lip, so the touch points are not free shortest-path points and the derivative has the lip's
      // share too. False: the unknowns' solve is singular.
      bool GradL3(const Path3 &p, const Geom3 &g, const Config &q, bool pl, double *gl, double (*gs)[kMaxGD] = nullptr, const SpanMap *M = nullptr, const double *wF = nullptr,
                  double *gF = nullptr) const
      {
         for (int a = 0; a < 3 * nb + 1; a++)
         {
            gl[a] = 0.0;
            if (gF)
               gF[a] = 0.0;
            if (gs)
               for (int j = 0; j < geo.nPosts; j++)
                  gs[j][a] = 0.0;
         }
         if (!Linearize3(p, g, q))
            return false;
         for (int i = 0; i < nb; i++)
         {
            const int e0 = Find3(p, kBall, i), e1 = Find3(p, kBallBack, i);
            if (e0 < 0 && e1 < 0)
               continue;
            for (int a = 0; a < 3; a++)
            {
               const V3 un(a == 0 ? 1.0 : 0.0, a == 1 ? 1.0 : 0.0, a == 2 ? 1.0 : 0.0);
               LinClear(p.n);
               if (e0 >= 0)
                  linDc[e0] = un;
               if (e1 >= 0)
                  linDc[e1] = un;
               gl[3 * i + a] = LinDeriv(p, g);
               if (gs)
                  LinSpans(p, g, *M, gs, 3 * i + a);
               if (gF)
                  gF[3 * i + a] = LinFloor(p, g, wF);
            }
         }
         if (pl)
         {
            bool face = false;
            LinClear(p.n);
            const V3 fv = V3(geo.plate.n.x, geo.plate.n.y, 0.0) * (geo.hasPlate ? FaceDD(q.phi) : 0.0);
            for (int k = 0; k < p.n; k++)
               if (p.e[k].kind == kFace)
               {
                  face = true;
                  linDc[k] = fv;
               }
            if (face)
            {
               gl[3 * nb] = LinDeriv(p, g);
               if (gs)
                  LinSpans(p, g, *M, gs, 3 * nb);
               if (gF)
                  gF[3 * nb] = LinFloor(p, g, wF);
            }
         }
         return true;
      }
      // the same by central differences of the path solve, 10 um, 1e-4 rad (debugGrad: the check)
      bool GradL3Num(const Path3 &p, const Config &q, bool pl, double *gl) const
      {
         for (int a = 0; a < 3 * nb + 1; a++)
            gl[a] = 0.0;
         const auto L = [&](const Config &qq, double &out) {
            Path3 pp = p;
            Geom3 gg;
            if (!Solve3(pp, qq, gg))
               return false;
            out = gg.L;
            return true;
         };
         for (int i = 0; i < nb; i++)
         {
            if (Find3(p, kBall, i) < 0 && Find3(p, kBallBack, i) < 0)
               continue;
            for (int a = 0; a < 3; a++)
            {
               const double d = 1e-5;
               Config qp = q, qm = q;
               if (a == 0)
               {
                  qp.c[i].x += d;
                  qm.c[i].x -= d;
               }
               else if (a == 1)
               {
                  qp.c[i].y += d;
                  qm.c[i].y -= d;
               }
               else
               {
                  qp.z[i] += d;
                  qm.z[i] -= d;
               }
               double Lp, Lm;
               if (!L(qp, Lp) || !L(qm, Lm))
                  return false;
               gl[3 * i + a] = (Lp - Lm) / (2.0 * d);
            }
         }
         bool face = false;
         for (int k = 0; k < p.n; k++)
            face = face || p.e[k].kind == kFace;
         if (pl && face)
         {
            const double d = 1e-4;
            Config qp = q, qm = q;
            qp.phi += d;
            qm.phi -= d;
            double Lp, Lm;
            if (!L(qp, Lp) || !L(qm, Lm))
               return false;
            gl[3 * nb] = (Lp - Lm) / (2.0 * d);
         }
         return true;
      }
      // s x (an approximate d2L/dq2: each run's bending stiffness (I - t t^T) / l between its two pulleys) added to H, band coordinates
      int ElemMap3(const Elem &e, double ph, bool pl, int *dof, V3 *col) const
      {
         if (e.kind == kBall || e.kind == kBallBack)
         {
            for (int a = 0; a < 3; a++)
               dof[a] = 3 * e.id + a;
            col[0] = V3(1.0, 0.0, 0.0);
            col[1] = V3(0.0, 1.0, 0.0);
            col[2] = V3(0.0, 0.0, 1.0);
            return 3;
         }
         if (e.kind == kFace && pl)
         {
            dof[0] = 3 * nb;
            col[0] = V3(geo.plate.n.x, geo.plate.n.y, 0.0) * FaceDD(ph);
            return 1;
         }
         return 0;
      }
      void AddHessL3(const Path3 &p, const Geom3 &g, const Config &q, bool pl, const double *sRun, double H[kMaxGD][kMaxGD]) const
      {
         const int n = p.n;
         const double lMin = 1e-3;
         for (int k = 0; k < n; k++)
         {
            const int k1 = (k + 1) % n;
            int da[3], db[3];
            V3 ca[3], cb[3];
            const int na = ElemMap3(p.e[k], q.phi, pl, da, ca), nbm = ElemMap3(p.e[k1], q.phi, pl, db, cb);
            if (na == 0 && nbm == 0)
               continue;
            const V3 t = g.t[k];
            const double il = sRun[k] / std::max(g.l[k], lMin);
            const auto Pm = [&](const V3 &x, const V3 &y) { return il * (x.Dot(y) - x.Dot(t) * y.Dot(t)); };
            for (int a = 0; a < na; a++)
               for (int b = 0; b < na; b++)
                  H[da[a]][da[b]] += Pm(ca[a], ca[b]);
            for (int a = 0; a < nbm; a++)
               for (int b = 0; b < nbm; b++)
                  H[db[a]][db[b]] += Pm(cb[a], cb[b]);
            for (int a = 0; a < na; a++)
               for (int b = 0; b < nbm; b++)
               {
                  const double x = Pm(ca[a], cb[b]);
                  H[da[a]][db[b]] -= x;
                  H[db[b]][da[a]] -= x;
               }
         }
      }
      // the cord's displacement toward blade s (group mode): the 3D path seen from above (each arc in chords of up to 10 deg)
      double SwitchDisp3(int s) const
      {
         const SwitchSpec &sw = geo.sw[s];
         double best = 1e30;
         const int n = path3.n;
         for (int k = 0; k < n; k++)
         {
            const int k1 = (k + 1) % n;
            V2 pts[40];
            int np = 0;
            pts[np++] = V2(path3.ct[k].X.x, path3.ct[k].X.y);
            pts[np++] = V2(path3.ct[k1].E.x, path3.ct[k1].E.y);
            const int na = std::min(36, std::max(1, (int)std::ceil(std::fabs(G3.wrap[k1]) / (10.0 * kPi / 180.0))));
            const V3 cen = G3.cen[k1], m3 = path3.ct[k1].m, e0 = path3.ct[k1].E - cen;
            const double sg = G3.rho[k1] >= 0.0 ? 1.0 : -1.0;
            for (int a = 1; a <= na && np < 40; a++)
            {
               const double th = sg * G3.wrap[k1] * a / na, cs = std::cos(th), sn = std::sin(th);
               const V3 pr = e0 * cs + m3.Cross(e0) * sn + m3 * (m3.Dot(e0) * (1.0 - cs));   // e0 turned th about m
               pts[np++] = V2(cen.x + pr.x, cen.y + pr.y);
            }
            for (int pp = 0; pp + 1 < np; pp++)
            {
               const V2 e = pts[pp + 1] - pts[pp];
               const double den = sw.b.Cross(e);
               if (std::fabs(den) < 1e-15)
                  continue;
               const V2 d = pts[pp] - sw.P;
               const double sPar = d.Cross(e) / den;
               const double u = d.Cross(sw.b) / den;
               if (u >= 0.0 && u <= 1.0)
                  best = std::min(best, sPar);
            }
         }
         return best < 1e29 ? best : 0.0;
      }

      V3 ToBand(const GV3 &u) const { const V3 a(u.x, u.y, u.z); return V3(a.Dot(frame.ex), a.Dot(frame.ey), a.Dot(frame.ez)); }
      GV3 ToWorld(const V3 &u) const { const V3 a = frame.ex * u.x + frame.ey * u.y + frame.ez * u.z; return GV3(a.x, a.y, a.z); }
      void Project(const GroupStep::Ball &b, V2 &cc, double &zz, V2 &vv) const
      {
         const V3 p = ToBand(b.c - GV3(frame.O.x, frame.O.y, frame.O.z)), u = ToBand(b.v);
         cc = V2(p.x, p.y);
         zz = p.z;
         vv = V2(u.x, u.y);
      }
      // the balls as the group has them; the path re-settled if any moved since the band last saw them
      void Sync(const GroupStep::Ball *b, int n, bool force = false)
      {
         n = std::min(n, kMaxBalls);
         force = force || n != nb;
         bool moved = force;
         nb = n;
         for (int i = 0; i < nb; i++)
         {
            V2 cc, vv;
            double zz;
            Project(b[i], cc, zz, vv);
            moved = moved || cc.x != c[i].x || cc.y != c[i].y || zz != z[i];
            c[i] = cc;
            z[i] = zz;
            v[i] = vv;
         }
         if (moved)
         {
            const Config q = Now();
            Geom3 g;
            const bool keep = gSamePath;
            gSamePath = !force;   // (the passes' touch points as they are reached: no element comes or goes)
            if (force ? !Settle3(path3, q, g) : !Solve3(path3, q, g))
               gdiag.settleFails++;
            Commit3(path3);
            gSamePath = keep;
         }
      }
      // the group's move x in the band's coordinates; the number of them
      int LocalMove(const double *x, int off, double *dq) const
      {
         for (int i = 0; i < nb; i++)
         {
            const V3 u = ToBand(GV3(x[6 * i], x[6 * i + 1], x[6 * i + 2]));
            dq[3 * i] = u.x;
            dq[3 * i + 1] = u.y;
            dq[3 * i + 2] = u.z;
         }
         int nd = 3 * nb;
         if (gPl)
            dq[nd++] = x[off];
         return nd;
      }
      Config Moved(const double *dq, double s) const
      {
         Config q = gq0;
         for (int i = 0; i < nb; i++)
         {
            q.c[i] = gq0.c[i] + V2(dq[3 * i], dq[3 * i + 1]) * s;
            q.z[i] = gq0.z[i] + dq[3 * i + 2] * s;
         }
         if (gPl)
            q.phi = gq0.phi + dq[3 * nb] * s;
         return q;
      }
      // the net torque on the plate now (whether it leaves the coil stop; group mode)
      double TauNow() const
      {
         double dLdphi = 0.0;
         const int n = path3.n;
         SpanMap smT;
         MakeSpanMap(path3, G3, smT);
         for (int e = 0; e < n; e++)
            if (path3.e[e].kind == kFace)
            {
               const V3 gk = G3.t[(e + n - 1) % n] - G3.t[e];
               const int sp = smT.arcSpan[e];
               dLdphi += SpanT(sp, spanL[sp]) * (gk.x * geo.plate.n.x + gk.y * geo.plate.n.y) * FaceDD(phi);
            }
         double tau = -dLdphi - PlateVp(phi) + CoilTorquePerI2(phi) * iHigh * iHigh;
         Dofs dp;
         dp.n = 2 * nb + 1;
         dp.phi = 2 * nb;
         const Config q = Now();
         for (int i = 0; i < nb; i++)
         {
            double gr[kMaxDof];
            const double del = Overlap(i, FaceBs(), q, dp, gr);
            tau -= SquashF(del) * gr[dp.phi];
         }
         return tau;
      }
      Config Now() const
      {
         Config q;
         for (int i = 0; i < nb; i++)
         {
            q.c[i] = c[i];
            q.z[i] = z[i];
            q.held[i] = holdOn[i];
            q.mh[i] = mHold[i];
         }
         q.phi = phi;
         return q;
      }

      // ---- the plate (one body: plate + plunger through the rigid link) ----
      double PhiStop() const { return set.sling.mech.stopDeg * FlipperModel::DEG; }
      double Inertia(double ph) const
      {
         if (!geo.hasPlate)
            return 1.0;
         const FlipperModel::Kin k = FlipperModel::Kinematics(set.sling.mech, der, der.theta0 + ph);
         return der.iFlipper + set.sling.mech.plungerKg * (k.lever * 1e-3) * (k.lever * 1e-3);
      }
      // return spring (on the plunger) + rest stop, J
      double PlateV(double ph) const
      {
         if (!geo.hasPlate)
            return 0.0;
         const FlipperSettings &s = set.sling.mech;
         const double u = der.xRest - FlipperModel::SliderX(s, der, der.theta0 + ph);   // mm the plunger has moved in from rest
         double e = 1e-3 * (s.preloadN * u + 0.5 * der.kSpringNPerMm * u * u);
         const double pen = -ph;
         if (pen > 0.0)
         {
            const double penRef = 1.4 * kPi / 180.0;
            e += s.restKNmPerRad * penRef * penRef / (s.restN + 1.0) * std::pow(pen / penRef, s.restN + 1.0);
         }
         return e;
      }
      double PlateVp(double ph) const { const double e = 1e-7; return (PlateV(ph + e) - PlateV(ph - e)) / (2.0 * e); }
      double PlateVpp(double ph) const { const double e = 1e-5; return (PlateV(ph + e) - 2.0 * PlateV(ph) + PlateV(ph - e)) / (e * e); }
      // rest stop's spring torque (out of the stop, + phi) and damping coefficient
      void RestStop(double ph, double &tauS, double &cD) const
      {
         tauS = cD = 0.0;
         const FlipperSettings &s = set.sling.mech;
         const double pen = -ph;
         if (pen <= 0.0)
            return;
         const double penRef = 1.4 * kPi / 180.0, rr = pen / penRef;
         tauS = s.restKNmPerRad * penRef * std::pow(rr, s.restN);
         cD = 2.0 * s.restZeta * std::sqrt(s.restKNmPerRad * Inertia(ph)) * std::pow(rr, 0.5 * (s.restN - 1.0));
      }
      // coil: force factor at the plunger place and the lever (torque = factor x i^2 x lever)
      double CoilTorquePerI2(double ph) const
      {
         const FlipperSettings &s = set.sling.mech;
         const FlipperModel::Kin k = FlipperModel::Kinematics(s, der, der.theta0 + ph);
         const double frac = (k.x - der.xMin) / s.bobbinMm;
         const double gHigh = FlipperModel::Sq(s.x0P / (frac + s.x0P));
         return der.strengthP * gHigh * k.lever * 1e-3;
      }
      // the coil current after h (fraction of V/R): the circuit's exponential with its resistance (and back-EMF) and inductance at the step's start
      double CurrentAfter(double h, bool coilOn) const
      {
         const FlipperSettings &s = set.sling.mech;
         const FlipperModel::Kin k = FlipperModel::Kinematics(s, der, der.theta0 + phi);
         const double frac = (k.x - der.xMin) / s.bobbinMm;
         const double dLdx = der.bP / (s.bobbinMm * FlipperModel::Sq(frac + s.x0P));   // mH/mm
         const double lP = der.aP + der.bP / (frac + s.x0P);                           // mH
         const double vPos = std::max(k.lever * omega * 1e-3, 0.0);                   // plunger closing speed, m/s (mm/ms)
         const double rEff = s.rP + dLdx * vPos;                                       // ohm
         const double issP = s.vSupply / s.rP;
         const double iss = (coilOn ? s.vSupply : 0.0) / rEff;
         const double i0 = iHigh * issP;
         const double i1 = iss + (i0 - iss) * std::exp(-h * 1e3 * rEff / lP);
         return std::max(i1, 0.0) / issP;
      }
      // the face's place along n and its rate with the plate angle
      mutable double faceAt = INFINITY, faceD = 0.0, faceDD = 0.0;   // the last angle asked for, its place and rate
      void FaceAt(double ph) const
      {
         if (ph == faceAt)
            return;
         const PlateGeom &p = geo.plate;
         const double t = ph + p.rest, ct = std::cos(t), st = std::sin(t);
         faceAt = ph;
         faceD = p.originN + p.a / ct + p.h * st / ct;
         faceDD = (p.a * st + p.h) / (ct * ct);
      }
      double FaceD(double ph) const { FaceAt(ph); return faceD; }
      double FaceDD(double ph) const { FaceAt(ph); return faceDD; }
      V2 FacePoint(int id, double ph) const
      {
         const PlateGeom &p = geo.plate;
         return p.A + (p.B - p.A) * (id == 0 ? p.f0 : p.f1) + p.n * FaceD(ph);
      }

      // ---- the cord ----
      double CordV(double L) const { const double s = L - Lr; return s > 0.0 ? 0.5 * EA / Lr * s * s : 0.0; }
      double CordVp(double L) const { const double s = L - Lr; return s > 0.0 ? EA / Lr * s : 0.0; }
      double CordChord(double L0, double L1) const
      {
         const double s0 = L0 - Lr, s1 = L1 - Lr;
         if (s0 >= 0.0 && s1 >= 0.0)
            return 0.5 * EA / Lr * (s0 + s1);
         if (s0 <= 0.0 && s1 <= 0.0)
            return 0.0;
         return (CordV(L1) - CordV(L0)) / (L1 - L0);
      }

      // round a ball whose centre is dz above the cord's centre line, the cord's centre line keeps this in-plane clearance from the ball's centre: the 3D distance
      // stays R + r (GroupStepDesign 10.1); 0 once the ball is clear above or below the cord
      double BallClear(double dz) const
      {
         const double S = R + r;
         return std::sqrt(std::max(S * S - dz * dz, 0.0));
      }
      void Pulley(const Elem &e, const Config &q, V2 &cen, double &rho) const
      {
         switch (e.kind)
         {
         case kPost: cen = geo.posts[e.id].c; rho = geo.posts[e.id].r; break;
         case kStop: cen = geo.stops[e.id].c; rho = geo.stops[e.id].r; break;
         case kFace: cen = FacePoint(e.id, q.phi); rho = r; break;
         case kBall: cen = q.c[e.id]; rho = -BallClear(q.z[e.id]); break;
         default: cen = q.c[e.id]; rho = BallClear(q.z[e.id]) + 2.0 * r; break;   // the ball with the dented front span round it, seen by the back span
         }
      }
      // the run leaving circle a (signed radius ra) and arriving at circle b: its direction and length (the two touch points are centre + rho x RightOf(t))
      static void Tangent(const V2 &a, double ra, const V2 &b, double rb, V2 &t, double &l)
      {
         const V2 d = b - a;
         const double D2 = d.Dot(d), dr = ra - rb;
         double beta;
         if (D2 > dr * dr)
         {
            l = std::sqrt(D2 - dr * dr);
            beta = std::atan2(dr, l);
         }
         else   // the cord squashed between the two (a backstop carries it from here)
         {
            l = 0.0;
            beta = dr >= 0.0 ? 0.5 * kPi : -0.5 * kPi;
         }
         const double th = std::atan2(d.y, d.x) + beta;
         t = V2(std::cos(th), std::sin(th));
      }
      void Build(const Path &p, const Config &q, PathGeom &g2) const
      {
         const int n = p.n;
         for (int k = 0; k < n; k++)
            Pulley(p.e[k], q, g2.cen[k], g2.rho[k]);
         for (int k = 0; k < n; k++)
         {
            const int k1 = (k + 1) % n;
            Tangent(g2.cen[k], g2.rho[k], g2.cen[k1], g2.rho[k1], g2.t[k], g2.l[k]);
         }
         g2.L = 0.0;
         for (int k = 0; k < n; k++)
         {
            const V2 &ti = g2.t[(k + n - 1) % n], &to = g2.t[k];
            double wr = std::atan2(ti.Cross(to), ti.Dot(to));   // the turn, left positive
            if (g2.rho[k] < 0.0)
               wr = -wr;
            wr += 2.0 * kPi * std::round((p.e[k].wrapRef - wr) / (2.0 * kPi));   // followed from the step's start
            g2.wrap[k] = wr;
            g2.L += g2.l[k] + std::fabs(g2.rho[k]) * wr;
         }
      }
      static int Find(const Path &p, int kind, int id)
      {
         for (int k = 0; k < p.n; k++)
            if (p.e[k].kind == kind && p.e[k].id == id)
               return k;
         return -1;
      }
      // the path at q, from p: a pulley leaves where its wrap goes below zero, joins where a run crosses it from its own side
      bool Settle(Path &p, const Config &q, PathGeom &g2, SettleLog *lg = nullptr) const
      {
         if (lg)
            lg->n = 0;
         const auto Log = [&](char act, const Elem &e, int at, double val) {
            if (!lg)
               return;
            SettleLog::Entry &en = lg->e[lg->n % SettleLog::kN];
            en.act = act;
            en.kind = e.kind;
            en.id = e.id;
            en.at = at;
            en.val = val;
            en.pathN = p.n;
            lg->n++;
         };
         for (int it = 0; it < 4 * kMaxPath; it++)
         {
            Build(p, q, g2);
            int worst = -1;
            double wmin = -1e-12;
            for (int k = 0; k < p.n; k++)
               if (g2.wrap[k] < wmin)
               {
                  wmin = g2.wrap[k];
                  worst = k;
               }
            if (worst >= 0 && p.n > 2)
            {
               Log('-', p.e[worst], worst, wmin);
               for (int k = worst; k + 1 < p.n; k++)
                  p.e[k] = p.e[k + 1];
               p.n--;
               continue;
            }
            int bestSeg = -1;
            Elem bestE;
            double bestDepth = 1e-12;
            const auto Try = [&](int kind, int id) {
               if (Find(p, kind, id) >= 0)
                  return;
               Elem e { kind, id, 0.0 };
               V2 cen;
               double rho;
               Pulley(e, q, cen, rho);
               for (int k = 0; k < p.n; k++)
               {
                  const V2 P0 = g2.cen[k] + RightOf(g2.t[k]) * g2.rho[k];
                  const V2 dc = cen - P0;
                  if (g2.t[k].Cross(dc) * rho <= 0.0)   // not on its side of this run
                     continue;
                  const double u = std::min(std::max(dc.Dot(g2.t[k]), 0.0), g2.l[k]);
                  const double depth = std::fabs(rho) - (dc - g2.t[k] * u).Length();
                  if (depth > bestDepth)
                  {
                     bestDepth = depth;
                     bestSeg = k;
                     bestE = e;
                  }
               }
            };
            for (int j = 0; j < geo.nPosts; j++)
               Try(kPost, j);
            for (int j = 0; j < geo.nStops; j++)
               Try(kStop, j);
            if (geo.hasPlate)
            {
               Try(kFace, 0);
               Try(kFace, 1);
            }
            for (int i = 0; i < nb; i++)
            {
               Try(kBall, i);
               if (Find(p, kBall, i) >= 0)
                  Try(kBallBack, i);
            }
            if (bestSeg < 0)
               return true;
            Log('+', bestE, bestSeg + 1, bestDepth);
            if (p.n >= kMaxPath)
               return false;
            for (int k = p.n; k > bestSeg + 1; k--)
               p.e[k] = p.e[k - 1];
            p.e[bestSeg + 1] = bestE;
            p.n++;
         }
         return false;
      }

      // the map from a pulley's centre to the unknowns: up to two (dof, column) pairs
      int ElemMap(const Elem &e, double ph, const Dofs &d, int *dof, V2 *col) const
      {
         if (e.kind == kBall || e.kind == kBallBack)
         {
            dof[0] = 2 * e.id;
            col[0] = V2(1.0, 0.0);
            dof[1] = 2 * e.id + 1;
            col[1] = V2(0.0, 1.0);
            return 2;
         }
         if (e.kind == kFace && d.phi >= 0)
         {
            dof[0] = d.phi;
            col[0] = geo.plate.n * FaceDD(ph);
            return 1;
         }
         return 0;
      }
      // dL/dq: a pulley's centre moved along d changes L by d . (t_in - t_out)
      void GradL(const Path &p, const PathGeom &g2, const Config &q, const Dofs &d, double *gr) const
      {
         for (int a = 0; a < d.n; a++)
            gr[a] = 0.0;
         const int n = p.n;
         for (int k = 0; k < n; k++)
         {
            int dof[2];
            V2 col[2];
            const int nm = ElemMap(p.e[k], q.phi, d, dof, col);
            const V2 gk = g2.t[(k + n - 1) % n] - g2.t[k];
            for (int a = 0; a < nm; a++)
               gr[dof[a]] += col[a].Dot(gk);
         }
      }
      // s x (an approximate d2L/dq2: each run's bending stiffness (I - t t^T) / l between its two pulleys) added to H
      void AddHessL(const Path &p, const PathGeom &g2, const Config &q, const Dofs &d, double s, double H[kMaxDof][kMaxDof]) const
      {
         const int n = p.n;
         const double lMin = 1e-3;
         for (int k = 0; k < n; k++)
         {
            const int k1 = (k + 1) % n;
            int da[2], db[2];
            V2 ca[2], cb[2];
            const int na = ElemMap(p.e[k], q.phi, d, da, ca), nbm = ElemMap(p.e[k1], q.phi, d, db, cb);
            if (na == 0 && nbm == 0)
               continue;
            const V2 t = g2.t[k];
            const double il = s / std::max(g2.l[k], lMin);
            const auto Pm = [&](const V2 &x, const V2 &y) { return il * (x.Dot(y) - x.Dot(t) * y.Dot(t)); };   // x^T (I - t t^T) y / l
            for (int a = 0; a < na; a++)
               for (int b = 0; b < na; b++)
                  H[da[a]][da[b]] += Pm(ca[a], ca[b]);
            for (int a = 0; a < nbm; a++)
               for (int b = 0; b < nbm; b++)
                  H[db[a]][db[b]] += Pm(cb[a], cb[b]);
            for (int a = 0; a < na; a++)
               for (int b = 0; b < nbm; b++)
               {
                  const double x = Pm(ca[a], cb[b]);
                  H[da[a]][db[b]] -= x;
                  H[db[b]][da[a]] -= x;
               }
         }
      }

      // ---- backstops: ball i against own post j, stop post, or the face (j counts posts, then stops, then the face) ----
      int NumBackstops() const { return geo.nPosts + geo.nStops + (geo.hasPlate ? 1 : 0) + geo.nWalls; }
      int FaceBs() const { return geo.nPosts + geo.nStops; }   // the plate face's backstop (when there is a plate); the walls follow it
      // The squash (the rubber's shared law, rubbercontact.h ContactCore, with the thin layer taken at the thickness left; docs/BandCordDesign.txt section 2):
      //    F = kH [delta^1.5 + beta reff^1.5 delta^3 / (h - delta)^3],  a = sqrt(reff delta), the layer term beta (a / (h - delta))^3 of Hertz
      // Rubber keeps its volume, so a layer pinched between steel only flows sideways, harder the thinner it is: the stiffness rises without limit as delta nears h.
      // Stored energy V = kH [0.4 delta^2.5 + beta reff^1.5 W(delta)], W = integral of x^3 / (h - x)^3 = 1.5 h + 3 h ln(h / u) + h^3 / (2 u^2) - 3 h^2 / u + u, u = h - delta.
      double LayerW(double del) const
      {
         const double u = hL - del;
         return 1.5 * hL + 3.0 * hL * std::log(hL / u) + hL * hL * hL / (2.0 * u * u) - 3.0 * hL * hL / u + u;
      }
      double SquashV(double del) const { return del > 0.0 ? 0.4 * kH * del * del * std::sqrt(del) + cL * LayerW(del) : 0.0; }
      double SquashF(double del) const
      {
         if (del <= 0.0)
            return 0.0;
         const double x = del / (hL - del);
         return kH * del * std::sqrt(del) + cL * x * x * x;
      }
      double SquashK(double del) const   // dF / d delta
      {
         if (del <= 0.0)
            return 0.0;
         const double u = hL - del, x = del / u;
         return 1.5 * kH * std::sqrt(del) + cL * 3.0 * x * x * (1.0 / u + del / (u * u));
      }
      // the move dq from q0 keeps every squash short of the rubber's thickness at the end and at mid-step
      bool InRange(const Config &q0, const double *dqv, const Dofs &d) const
      {
         Config q1 = q0, qm = q0;
         for (int i = 0; i < nb; i++)
         {
            q1.c[i] = q0.c[i] + V2(dqv[2 * i], dqv[2 * i + 1]);
            qm.c[i] = q0.c[i] + V2(dqv[2 * i], dqv[2 * i + 1]) * 0.5;
         }
         if (d.phi >= 0)
         {
            q1.phi = q0.phi + dqv[d.phi];
            qm.phi = q0.phi + 0.5 * dqv[d.phi];
         }
         const double lim = hL * (1.0 - 1e-9);
         for (int i = 0; i < nb; i++)
            for (int j = 0; j < NumBackstops(); j++)
               if (Overlap(i, j, q1, d, nullptr) >= lim || Overlap(i, j, qm, d, nullptr) >= lim)
                  return false;
         return true;
      }
      double SquashChord(double d0, double d1, double dm) const
      {
         if (d0 <= 0.0 && d1 <= 0.0)
            return 0.0;
         if (std::fabs(d1 - d0) < 1e-9)   // the midpoint force: the difference would be mostly rounding (energy error ~F'' dd^3)
            return SquashF(dm);
         double hz;
         if (d0 > 0.0 && d1 > 0.0)   // Hertz part: (s1^5 - s0^5) / (s1^2 - s0^2) with s = sqrt(delta), without the cancellation
         {
            const double s0 = std::sqrt(d0), s1 = std::sqrt(d1);
            const double s0q = s0 * s0, s1q = s1 * s1;
            hz = 0.4 * kH * (s1q * s1q + s1q * s1 * s0 + s1q * s0q + s1 * s0q * s0 + s0q * s0q) / (s1 + s0);
         }
         else
         {
            const double a = std::max(d0, 0.0), b = std::max(d1, 0.0);
            hz = 0.4 * kH * (b * b * std::sqrt(b) - a * a * std::sqrt(a)) / (d1 - d0);
         }
         const double lw = cL * ((d1 > 0.0 ? LayerW(d1) : 0.0) - (d0 > 0.0 ? LayerW(d0) : 0.0)) / (d1 - d0);
         return hz + lw;
      }
      // the overlap of ball i with backstop j at q (> 0 = squashing the cord), and its gradient in the unknowns (gr may be null; grz: its change with the ball's
      // height). The squash is R + r less the distance from the ball's centre to where the cord's centre line can lie against the backstop: on an own post its
      // groove at the cord's height (the 3D distance); on every other surface (a stop post, the plate's face, a wall) anywhere up or down it, as the cord slides
      // (the in-plane distance less r off the surface). One law for every surface (GroupStepDesign 11.2).
      // A squash grows from zero: a ball found already inside a backstop (over the top of a stop post, a contact found late) has stored nothing, so that depth is its zero
      // until it has backed out of it (bsZero, Commit3). No energy can appear when a ball joins (GroupStepDesign 12).
      double Overlap(int i, int j, const Config &q, const Dofs &d, double *gr, double *grz = nullptr) const
      {
         if (groupMode && j < geo.nPosts + geo.nStops)   // a group contact (AddPostShapes): the plate and the walls stay here
         {
            if (gr)
               for (int a = 0; a < d.n; a++)
                  gr[a] = 0.0;
            if (grz)
               *grz = 0.0;
            return -1.0;
         }
         return OverlapRaw(i, j, q, d, gr, grz) - bsZero[i].d[j];
      }
      double OverlapRaw(int i, int j, const Config &q, const Dofs &d, double *gr, double *grz = nullptr) const
      {
         if (gr)
            for (int a = 0; a < d.n; a++)
               gr[a] = 0.0;
         if (j < geo.nPosts && geo.postTwin[j])
            return -1.0;
         const V2 cb = q.c[i];
         double a;      // the ball centre's in-plane distance from where the cord's centre line lies against it
         V2 ga;         // its gradient in the ball's centre
         double gaPhi = 0.0;
         if (j < geo.nPosts || (groupMode && j < geo.nPosts + geo.nStops))   // its groove at the cord's height: the 3D distance to the groove's centre line (ProbePost)
         {
            const Circle &cs = j < geo.nPosts ? geo.posts[j] : geo.stops[j - geo.nPosts];
            const V2 dv = cb - cs.c;
            const double D = dv.Length(), aa = D - cs.r, d3 = std::sqrt(aa * aa + q.z[i] * q.z[i]);
            const double fa = d3 > 0.0 ? aa / d3 : 1.0, fz = d3 > 0.0 ? q.z[i] / d3 : 0.0;
            if (gr && D > 0.0)
            {
               gr[2 * i] = -fa * dv.x / D;
               gr[2 * i + 1] = -fa * dv.y / D;
            }
            if (grz)
               *grz = -fz;
            return R + r - d3;
         }
         if (j < geo.nPosts + geo.nStops)
         {
            const Circle &cs = geo.stops[j - geo.nPosts];
            const V2 dv = cb - cs.c;
            const double D = dv.Length();
            a = D - cs.r;
            ga = D > 0.0 ? dv * (1.0 / D) : V2();
         }
         else if (geo.hasPlate && j == FaceBs())
         {
            const PlateGeom &p = geo.plate;
            const V2 F0 = FacePoint(0, q.phi), F1 = FacePoint(1, q.phi), e = F1 - F0;
            const double u = (cb - F0).Dot(e) / e.Dot(e);
            const double dd = FaceDD(q.phi);
            double dist;
            if (u > 0.0 && u < 1.0)
            {
               const double sd = (cb - F0).Dot(p.n);   // the distance to the plate, from either side: a ball that passes its end from behind meets no step in its depth (a plane's side has one: it is deep inside at once)
               dist = std::fabs(sd);
               ga = sd >= 0.0 ? p.n : p.n * -1.0;
            }
            else
            {
               const V2 dv = cb - (u <= 0.0 ? F0 : F1);
               dist = dv.Length();
               ga = dist > 0.0 ? dv * (1.0 / dist) : p.n;
            }
            a = dist - r;
            gaPhi = -ga.Dot(p.n) * dd;
         }
         else
         {
            const WallSpec &wl = geo.walls[j - FaceBs() - (geo.hasPlate ? 1 : 0)];
            const V2 e = wl.B - wl.A;
            const double u = std::min(std::max((cb - wl.A).Dot(e) / e.Dot(e), 0.0), 1.0);
            const V2 dv = cb - (wl.A + e * u);
            const double dist = dv.Length();
            a = dist - wl.halfThick - r;
            ga = dist > 0.0 ? dv * (1.0 / dist) : V2();
         }
         if (gr)
         {
            gr[2 * i] = -ga.x;
            gr[2 * i + 1] = -ga.y;
            if (d.phi >= 0)
               gr[d.phi] = -gaPhi;
         }
         if (grz)
            *grz = 0.0;
         return R + r - a;
      }

      // A ball's grip on the cord over one step: friction impulse P (P1 along the cord, P2 across the band's plane) on the disc |P| <= LP and the rolling
      // resistance impulse M (about z: rolling along the cord; spin about the cord's direction is slip across the plane, which the friction takes) on |M| <= LM.
      // Mid-step rates u = u0 + K x, x = (P1, P2, M):
      //    s1 = s1_0 + kPP P1 - c M,   s2 = s2_0 + kPz P2,   b = b_0 - c P1 + kM M
      // Maximum dissipation: x minimises 1/2 x^T K x + u0 . x over the two limits, so each impulse at its limit opposes its rate and one within its limit holds its
      // rate at zero: the work x . u is never positive. KKT: x = -(K + D)^-1 u0, D = diag(bP, bP, bM), bP and bM >= 0 found by the limits.
      struct Grip
      {
         double s1 = 0.0, s2 = 0.0, b = 0.0;                     // rates with no impulse
         double kPP = 0.0, kPz = 0.0, kM = 0.0, c = 0.0;         // compliance
         double LP = 0.0, LM = 0.0;                              // limits (impulse over the step)
         double P1 = 0.0, P2 = 0.0, M = 0.0;                     // result
      };
      // x for given multipliers (infinite = that impulse is zero: a zero limit)
      static void GripAt(Grip &g, double bP, double bM, bool pZero, bool mZero)
      {
         if (pZero)
         {
            g.P1 = g.P2 = 0.0;
            g.M = mZero ? 0.0 : -g.b / (g.kM + bM);
            return;
         }
         g.P2 = -g.s2 / (g.kPz + bP);
         if (mZero)
         {
            g.P1 = -g.s1 / (g.kPP + bP);
            g.M = 0.0;
            return;
         }
         // [kPP + bP, -c; -c, kM + bM] [P1; M] = -[s1; b]
         const double a11 = g.kPP + bP, a22 = g.kM + bM, det = a11 * a22 - g.c * g.c;
         g.P1 = (-g.s1 * a22 - g.c * g.b) / det;
         g.M = (-g.b * a11 - g.c * g.s1) / det;
      }
      // the root of a decreasing f on [0, inf) where f(0) > 0 (bracket, then Illinois false position)
      template <class Fn> static double DecreasingRoot(const Fn &f)
      {
         double lo = 0.0, flo = f(0.0), hi = 1.0, fhi = f(hi);
         while (fhi > 0.0 && hi < 1e300)
         {
            lo = hi;
            flo = fhi;
            hi *= 8.0;
            fhi = f(hi);
         }
         int side = 0;
         for (int it = 0; it < 200; it++)
         {
            const double x = (lo * fhi - hi * flo) / (fhi - flo);
            const double fx = f(x);
            if (fx == 0.0 || (hi - lo) <= 1e-15 * hi)
               return x;
            if (fx > 0.0)
            {
               lo = x;
               flo = fx;
               if (side == 1)
                  fhi *= 0.5;
               side = 1;
            }
            else
            {
               hi = x;
               fhi = fx;
               if (side == -1)
                  flo *= 0.5;
               side = -1;
            }
            if (std::fabs(fx) < 1e-14 * (1.0 + std::fabs(flo)))
               return x;
         }
         return 0.5 * (lo + hi);
      }
      static void SolveGrip(Grip &g)
      {
         const bool pZero = g.LP <= 0.0, mZero = g.LM <= 0.0;
         // the rolling multiplier for a given friction multiplier: 0 if M is within its limit, else the one that puts it on the limit
         const auto BM = [&](double bP) {
            if (mZero)
               return 0.0;
            GripAt(g, bP, 0.0, pZero, false);
            if (std::fabs(g.M) <= g.LM)
               return 0.0;
            return DecreasingRoot([&](double bM) { GripAt(g, bP, bM, pZero, false); return std::fabs(g.M) - g.LM; });
         };
         double bP = 0.0;
         if (!pZero)
         {
            const auto PMag = [&](double b) { const double bm = BM(b); GripAt(g, b, bm, false, mZero); return std::sqrt(g.P1 * g.P1 + g.P2 * g.P2) - g.LP; };
            if (PMag(0.0) > 0.0)
               bP = DecreasingRoot(PMag);
         }
         const double bM = BM(bP);
         GripAt(g, bP, bM, pZero, mZero);
      }

      // the discrete gradient (Gonzalez): gm + (df - gm . dq) / (dq^T W dq) W dq, so that the result . dq = df exactly. A mismatch within dfNoise (the rounding of
      // f itself) carries no information and is left out: dividing it by a tiny move would only amplify the rounding (the energy error left is below f's rounding).
      static void Gonzalez(const double *gm, double df, double dfNoise, const double *dq, const double *W, int n, double *out)
      {
         double gdq = 0.0, qq = 0.0;
         for (int a = 0; a < n; a++)
         {
            gdq += gm[a] * dq[a];
            qq += dq[a] * W[a] * dq[a];
         }
         const double cf = qq > 1e-300 && std::fabs(df - gdq) > dfNoise ? (df - gdq) / qq : 0.0;
         for (int a = 0; a < n; a++)
            out[a] = gm[a] + cf * W[a] * dq[a];
      }
      static bool SolveLinear(double A[kMaxDof][kMaxDof], double *b, int n)
      {
         for (int col = 0; col < n; col++)
         {
            int piv = col;
            for (int rw = col + 1; rw < n; rw++)
               if (std::fabs(A[rw][col]) > std::fabs(A[piv][col]))
                  piv = rw;
            if (std::fabs(A[piv][col]) < 1e-300)
               return false;
            if (piv != col)
            {
               for (int k = 0; k < n; k++)
                  std::swap(A[col][k], A[piv][k]);
               std::swap(b[col], b[piv]);
            }
            for (int rw = col + 1; rw < n; rw++)
            {
               const double f = A[rw][col] / A[col][col];
               if (f == 0.0)
                  continue;
               for (int k = col; k < n; k++)
                  A[rw][k] -= f * A[col][k];
               b[rw] -= f * b[col];
            }
         }
         for (int rw = n - 1; rw >= 0; rw--)
         {
            double s = b[rw];
            for (int k = rw + 1; k < n; k++)
               s -= A[rw][k] * b[k];
            b[rw] = s / A[rw][rw];
         }
         return true;
      }

      // ---- the step ----
      struct Trial
      {
         Config q1;
         V2 v1[kMaxBalls];
         V3 w1[kMaxBalls];
         double omega1 = 0.0, iHigh1 = 0.0, L1 = 0.0;
         Path path1;
         PathGeom G1;
         double Pf[kMaxBalls], Pz[kMaxBalls], Mt[kMaxBalls];
         V3 Mr[kMaxBalls];
         double coilJ = 0.0, cordLossJ = 0.0, backstopLossJ = 0.0, frictionJ = 0.0, twistJ = 0.0, rollJ = 0.0, plateLossJ = 0.0, stopJ = 0.0, maxOverlap = 0.0;
         double tauEnd = 0.0;   // net torque on the plate at the end (its coil stop's release)
         int rounds = 0;
         bool pinNow = false;
      };

      // Solve one step of h from the state now; the end state into T (not taken)
      bool Solve(double h, Trial &T) const
      {
         const Dofs d = MakeDofs();
         const int n = d.n;
         const int nPairs = nb * NumBackstops();
         const Config q0 = Now();
         const bool pl = d.phi >= 0;
         const double I0 = Inertia(phi), p0 = I0 * omega;
         const bool coilOn = geo.hasPlate && driver.Coil();
         T.iHigh1 = geo.hasPlate ? CurrentAfter(h, coilOn) : 0.0;
         const double i2 = 0.5 * (iHigh * iHigh + T.iHigh1 * T.iHigh1);
         double W[kMaxDof];
         for (int a = 0; a < 2 * nb; a++)
            W[a] = 1.0 / m;
         if (pl)
            W[d.phi] = 1.0 / I0;
         const double L0 = Lcur;
         // rounding of the path's length and of an overlap: positions up to posScale from the origin, a few ulps per run
         double posScale = 0.1;
         for (int k = 0; k < geo.nPosts; k++)
            posScale = std::max(posScale, std::fabs(geo.posts[k].c.x) + std::fabs(geo.posts[k].c.y));
         const double eps = 2.220446049250313e-16;
         const double noiseL = 16.0 * eps * (path.n + 2) * (posScale + L0), noiseD = 16.0 * eps * posScale;
         double del0[kMaxBalls * (kMaxPosts + kMaxStops + 1)];
         for (int i = 0; i < nb; i++)
            for (int j = 0; j < NumBackstops(); j++)
               del0[i * NumBackstops() + j] = Overlap(i, j, q0, d, nullptr);

         // the unknowns are the moves over the step (not the end positions: a position 0.1 m from the origin holds only ~1e-17 m, too coarse for a 1 us step's speed)
         double dq[kMaxDof];
         for (int i = 0; i < nb; i++)
         {
            dq[2 * i] = v[i].x * h;
            dq[2 * i + 1] = v[i].y * h;
         }
         if (pl)
            dq[d.phi] = omega * h;
         for (int it = 0; !InRange(q0, dq, d); it++)   // a first guess past the squash's range is shortened
         {
            if (it == 60)
               return false;
            for (int a = 0; a < n; a++)
               dq[a] *= 0.5;
         }
         Config q1 = q0;
         for (int i = 0; i < nb; i++)   // friction warm start: last step's forces over this step
         {
            T.Pf[i] = Pf[i] * h;
            T.Pz[i] = Pz[i] * h;
            T.Mt[i] = Mt[i] * h;
            T.Mr[i] = Mr[i] * h;
         }

         double gLm[kMaxDof], gL[kMaxDof], gdm[kMaxDof], gd[kMaxDof], F[kMaxDof], Res[kMaxDof];
         double H[kMaxDof][kMaxDof];
         V2 Fn[kMaxBalls];
         Path pm;
         PathGeom Gm;
         bool converged = false;
         int round = 0;
         double cordLoss = 0.0, bsLoss = 0.0, plateLoss = 0.0, coilW = 0.0, maxOv = 0.0;
         for (round = 1; round <= set.maxRounds; round++)
         {
            Config qm = q0;
            for (int i = 0; i < nb; i++)
            {
               q1.c[i] = q0.c[i] + V2(dq[2 * i], dq[2 * i + 1]);
               qm.c[i] = q0.c[i] + V2(dq[2 * i], dq[2 * i + 1]) * 0.5;
            }
            if (pl)
            {
               q1.phi = q0.phi + dq[d.phi];
               qm.phi = q0.phi + 0.5 * dq[d.phi];
            }
            T.path1 = path;
            pm = path;
            if (!Settle(T.path1, q1, T.G1) || !Settle(pm, qm, Gm))
               return false;
            const double L1 = T.G1.L, dL = L1 - L0;
            GradL(pm, Gm, qm, d, gLm);
            Gonzalez(gLm, dL, noiseL, dq, W, n, gL);
            const double Vc = CordChord(L0, L1);
            for (int a = 0; a < n; a++)
            {
               F[a] = -(Vc + cv * dL / h) * gL[a];
               for (int b = 0; b < n; b++)
                  H[a][b] = 0.0;
            }
            cordLoss = cv * dL * dL / h;
            // stiffness for the Newton rounds: stretch, bend, viscous
            const double Tm = CordVp(Gm.L), kAx = Gm.L > Lr ? EA / Lr : 0.0;
            for (int a = 0; a < n; a++)
               for (int b = 0; b < n; b++)
                  H[a][b] += (0.5 * h * kAx + cv) * gLm[a] * gLm[b];
            AddHessL(pm, Gm, qm, d, 0.5 * h * Tm, H);
            // backstops
            bsLoss = 0.0;
            maxOv = 0.0;
            for (int i = 0; i < nb; i++)
               Fn[i] = V2(F[2 * i], F[2 * i + 1]);
            for (int pr = 0; pr < nPairs; pr++)
            {
               const int i = pr / NumBackstops(), j = pr % NumBackstops();
               const double d0 = del0[pr], d1 = Overlap(i, j, q1, d, nullptr), dm = Overlap(i, j, qm, d, gdm);
               maxOv = std::max(maxOv, d1);
               if (d0 <= 0.0 && d1 <= 0.0 && dm <= 0.0)
                  continue;
               Gonzalez(gdm, d1 - d0, noiseD, dq, W, n, gd);
               const double cD = dm > 0.0 ? 1.5 * set.alphaSPerM * SquashF(dm) : 0.0;
               const double fv = SquashChord(d0, d1, dm) + cD * (d1 - d0) / h;
               for (int a = 0; a < n; a++)
                  F[a] -= fv * gd[a];
               Fn[i] = Fn[i] - V2(gd[2 * i], gd[2 * i + 1]) * fv;
               bsLoss += cD * (d1 - d0) * (d1 - d0) / h;
               const double kd = SquashK(dm);
               for (int a = 0; a < n; a++)
                  for (int b = 0; b < n; b++)
                     H[a][b] += (0.5 * h * kd + cD) * gdm[a] * gdm[b];
            }
            for (int i = 0; i < nb; i++)
            {
               F[2 * i] += m * g.x;
               F[2 * i + 1] += m * g.y;
            }
            // the plate
            double pEnd = 0.0, I1 = I0;
            plateLoss = 0.0;
            coilW = 0.0;
            if (pl)
            {
               const double dphi = dq[d.phi], phm = qm.phi;
               I1 = Inertia(q1.phi);
               const double Jb = 0.5 * (1.0 / I0 + 1.0 / I1);
               pEnd = 2.0 * dphi / (h * Jb) - p0;
               double dTdphi;
               // exact differences over the turn; below 1e-7 rad the midpoint slope (the difference would be mostly rounding; the energy error left is ~V''' dphi^3)
               if (std::fabs(dphi) > 1e-7)
                  dTdphi = (1.0 / I1 - 1.0 / I0) / dphi * 0.25 * (p0 * p0 + pEnd * pEnd);
               else
               {
                  const double e = 1e-6;
                  dTdphi = (1.0 / Inertia(phm + e) - 1.0 / Inertia(phm - e)) / (2.0 * e) * 0.25 * (p0 * p0 + pEnd * pEnd);
               }
               const double Vp = std::fabs(dphi) > 1e-7 ? (PlateV(q1.phi) - PlateV(q0.phi)) / dphi : PlateVp(phm);
               double tauS, cR;
               RestStop(phm, tauS, cR);
               const double om = dphi / h;
               const double tauD = tauS > 0.0 ? std::max(-cR * om, -tauS) : 0.0;
               const double tauC = coilOn || iHigh > 0.0 || T.iHigh1 > 0.0 ? CoilTorquePerI2(phm) * i2 : 0.0;
               const double b = set.sling.mech.bDamping;
               F[d.phi] += -Vp + tauD - b * om + tauC - dTdphi;
               plateLoss = -(tauD - b * om) * dphi;
               coilW = tauC * dphi;
               H[d.phi][d.phi] += 2.0 / (h * Jb) + 0.5 * h * std::max(PlateVpp(phm), 0.0) + b + cR;
            }
            // residual and the Newton round
            for (int i = 0; i < nb; i++)
            {
               const V2 nrm = Fn[i].Length() > 1e-12 ? Fn[i] * (1.0 / Fn[i].Length()) : V2();
               const V2 t1(-nrm.y, nrm.x);
               Res[2 * i] = m * (2.0 * dq[2 * i] / h - 2.0 * v[i].x) - h * F[2 * i] - T.Pf[i] * t1.x;
               Res[2 * i + 1] = m * (2.0 * dq[2 * i + 1] / h - 2.0 * v[i].y) - h * F[2 * i + 1] - T.Pf[i] * t1.y;
            }
            if (pl)
               Res[d.phi] = pEnd - p0 - h * F[d.phi];
            for (int a = 0; a < 2 * nb; a++)
               H[a][a] += 2.0 * m / h;
            double step[kMaxDof];
            for (int a = 0; a < n; a++)
            {
               if (!std::isfinite(Res[a]))
                  return false;
               step[a] = -Res[a];
            }
            if (n > 0 && !SolveLinear(H, step, n))
               return false;
            // limit a round's move (2 mm, 0.05 rad)
            double scale = 1.0;
            for (int i = 0; i < nb; i++)
            {
               const double s = std::sqrt(step[2 * i] * step[2 * i] + step[2 * i + 1] * step[2 * i + 1]);
               if (s > 2e-3)
                  scale = std::min(scale, 2e-3 / s);
            }
            if (pl && std::fabs(step[d.phi]) > 0.05)
               scale = std::min(scale, 0.05 / std::fabs(step[d.phi]));
            // and keep it inside the squash's range (no rubber squashed to nothing)
            {
               double trial[kMaxDof];
               int it = 0;
               for (; it < 60; it++, scale *= 0.5)
               {
                  for (int a = 0; a < n; a++)
                     trial[a] = dq[a] + step[a] * scale;
                  if (InRange(q0, trial, d))
                     break;
               }
               if (it == 60)
                  return false;
            }
            double moveV = 0.0, moveW = 0.0;
            for (int i = 0; i < nb; i++)
            {
               dq[2 * i] += step[2 * i] * scale;
               dq[2 * i + 1] += step[2 * i + 1] * scale;
               moveV = std::max(moveV, std::sqrt(step[2 * i] * step[2 * i] + step[2 * i + 1] * step[2 * i + 1]) * scale / h);
            }
            if (pl)
            {
               dq[d.phi] += step[d.phi] * scale;
               moveW = std::fabs(step[d.phi]) * scale / h;
            }
            // the ball's grip on the new move: friction, rolling resistance and twist solved together (they act on the same spin; solved one after the other they
            // undo 5/7 of each other every round: runs/2026-10-29-band-cord rounds trace)
            double fricChange = 0.0, chSlip = 0.0, chTwist = 0.0, chRoll = 0.0;
            int stick0 = -1, twistHeld0 = -1, rollHeld0 = -1;
            for (int i = 0; i < nb; i++)
            {
               const double N = Fn[i].Length();
               if (N < 1e-9)
               {
                  chSlip = std::max(chSlip, (std::fabs(T.Pf[i]) + std::fabs(T.Pz[i])) / m);
                  chRoll = std::max(chRoll, T.Mr[i].Length() / (Ib / R));
                  chTwist = std::max(chTwist, std::fabs(T.Mt[i]) / (Ib / R));
                  T.Pf[i] = T.Pz[i] = T.Mt[i] = 0.0;
                  T.Mr[i] = V3();
                  continue;
               }
               const V3 n3(Fn[i].x / N, Fn[i].y / N, 0.0), t1(-n3.y, n3.x, 0.0), ez(0.0, 0.0, 1.0);
               // the mid-step rates with no grip impulse (the Newton move already carries this round's in-plane friction: take its share back out)
               const double V0 = dq[2 * i] / h * t1.x + dq[2 * i + 1] / h * t1.y - 0.5 * T.Pf[i] / m;
               const double A0 = w[i].Dot(t1), B0 = w[i].z, C0 = w[i].Dot(n3);
               Grip g;
               g.s1 = V0 - R * B0;   // slip along the cord (contact point at -R n: slip = v + w x rc)
               g.s2 = R * A0;        // slip across the band's plane
               g.b = B0;             // rolling spin along the cord (about z)
               g.kPP = 0.5 / m + 0.5 * R * R / Ib;
               g.kPz = 0.5 * R * R / Ib;
               g.kM = 0.5 / Ib;
               g.c = 0.5 * R / Ib;
               g.LP = set.mu * N * h;
               g.LM = set.rollingResistance * R * N * h;
               SolveGrip(g);
               const double M2 = g.M;
               // twist about the normal: on its own (no other impulse turns the ball about n), stick up to mu N a
               const double tl = set.mu * N * set.twistRadiusM * h;
               const double Mtw = std::min(std::max(-2.0 * Ib * C0, -tl), tl);
               if (i == 0)
               {
                  stick0 = g.P1 * g.P1 + g.P2 * g.P2 < g.LP * g.LP * (1.0 - 1e-9) || g.LP == 0.0 ? 1 : 0;
                  rollHeld0 = std::fabs(M2) < g.LM * (1.0 - 1e-9) || g.LM == 0.0 ? 1 : 0;
                  twistHeld0 = std::fabs(2.0 * Ib * C0) <= tl ? 1 : 0;
               }
               const double dP = std::sqrt((g.P1 - T.Pf[i]) * (g.P1 - T.Pf[i]) + (g.P2 - T.Pz[i]) * (g.P2 - T.Pz[i])) / m;
               const double dM = std::fabs(M2 - T.Mr[i].z) / (Ib / R), dT = std::fabs(Mtw - T.Mt[i]) / (Ib / R);
               chSlip = std::max(chSlip, dP);
               chRoll = std::max(chRoll, dM);
               chTwist = std::max(chTwist, dT);
               T.Pf[i] = g.P1;
               T.Pz[i] = g.P2;
               T.Mr[i] = ez * M2;
               T.Mt[i] = Mtw;
            }
            fricChange = std::max(chSlip, std::max(chRoll, chTwist));
            if (diag.on && round <= SolveDiag::kRounds)
            {
               SolveDiag::Round &rd = diag.r[round - 1];
               rd.moveV = moveV;
               rd.moveW = moveW;
               rd.slip = chSlip;
               rd.twist = chTwist;
               rd.roll = chRoll;
               rd.stick0 = stick0;
               rd.twistHeld0 = twistHeld0;
               rd.rollHeld0 = rollHeld0;
               rd.scale = scale;
               diag.nRound = round;
            }
            if (moveV <= set.tolV && moveW <= set.tolW && fricChange <= set.tolV)
            {
               converged = true;
               break;
            }
         }
         T.rounds = std::min(round, set.maxRounds);
         if (!converged)
            return false;

         // the end state
         for (int i = 0; i < nb; i++)
            q1.c[i] = q0.c[i] + V2(dq[2 * i], dq[2 * i + 1]);
         if (pl)
            q1.phi = q0.phi + dq[d.phi];
         Settle(T.path1, q1, T.G1);
         T.q1 = q1;
         T.L1 = T.G1.L;
         T.frictionJ = T.twistJ = T.rollJ = 0.0;
         for (int i = 0; i < nb; i++)
         {
            T.v1[i] = V2(dq[2 * i], dq[2 * i + 1]) * (2.0 / h) - v[i];
            const double N = Fn[i].Length();
            if (N < 1e-9)
            {
               T.w1[i] = w[i];
               continue;
            }
            const V3 n3(Fn[i].x / N, Fn[i].y / N, 0.0), t1(-n3.y, n3.x, 0.0), ez(0.0, 0.0, 1.0);
            const V3 rc = n3 * -R;
            const V3 P = t1 * T.Pf[i] + ez * T.Pz[i];
            const V3 dw = (rc.Cross(P) + n3 * T.Mt[i] + T.Mr[i]) * (1.0 / Ib);
            T.w1[i] = w[i] + dw;
            const V3 wm = w[i] + dw * 0.5;
            const V3 vm(dq[2 * i] / h, dq[2 * i + 1] / h, 0.0);
            T.frictionJ += -P.Dot(vm + wm.Cross(rc));
            T.twistJ += -T.Mt[i] * wm.Dot(n3);
            T.rollJ += -T.Mr[i].Dot(wm - n3 * wm.Dot(n3));
         }
         if (pl)
         {
            const double I1 = Inertia(q1.phi), Jb = 0.5 * (1.0 / I0 + 1.0 / I1);
            T.omega1 = (2.0 * dq[d.phi] / (h * Jb) - p0) / I1;
         }
         else
         {
            T.q1.phi = phi;
            T.omega1 = 0.0;
         }
         T.cordLossJ = cordLoss;
         T.backstopLossJ = bsLoss;
         T.plateLossJ = plateLoss;
         T.coilJ = coilW;
         T.maxOverlap = maxOv;
         // the net torque on the plate at the end (whether it leaves the coil stop)
         T.tauEnd = 0.0;
         if (geo.hasPlate)
         {
            double dLdphi = 0.0;
            for (int e = 0; e < T.path1.n; e++)
               if (T.path1.e[e].kind == kFace)
                  dLdphi += (T.G1.t[(e + T.path1.n - 1) % T.path1.n] - T.G1.t[e]).Dot(geo.plate.n) * FaceDD(T.q1.phi);
            double tau = -CordVp(T.L1) * dLdphi - PlateVp(T.q1.phi) + CoilTorquePerI2(T.q1.phi) * T.iHigh1 * T.iHigh1;
            Dofs dp;
            dp.n = 2 * nb + 1;
            dp.phi = 2 * nb;
            double gr[kMaxDof];
            for (int i = 0; i < nb; i++)
            {
               const double del = Overlap(i, FaceBs(), T.q1, dp, gr);
               tau -= SquashF(del) * gr[dp.phi];
            }
            T.tauEnd = tau;
         }
         return true;
      }

      void Accept(const Trial &T, double h)
      {
         for (int i = 0; i < nb; i++)
         {
            c[i] = T.q1.c[i];
            v[i] = T.v1[i];
            w[i] = T.w1[i];
            Pf[i] = T.Pf[i] / h;
            Pz[i] = T.Pz[i] / h;
            Mt[i] = T.Mt[i] / h;
            Mr[i] = T.Mr[i] * (1.0 / h);
         }
         phi = T.q1.phi;
         omega = T.omega1;
         iHigh = T.iHigh1;
         if (T.pinNow)
            pinned = true;
         else if (pinned && T.tauEnd < 0.0)
            pinned = false;
         path = T.path1;
         Commit(path);
         if (geo.hasPlate)
            driver.Tick(h);
         books.steps++;
         books.rounds += T.rounds;
         books.maxRounds = std::max(books.maxRounds, T.rounds);
         books.hMinUsed = std::min(books.hMinUsed, h);
         books.coilJ += T.coilJ;
         books.cordLossJ += T.cordLossJ;
         books.backstopLossJ += T.backstopLossJ;
         books.frictionJ += T.frictionJ;
         books.twistJ += T.twistJ;
         books.rollJ += T.rollJ;
         books.plateLossJ += T.plateLossJ;
         books.stopJ += T.stopJ;
         books.maxOverlapM = std::max(books.maxOverlapM, T.maxOverlap);
      }
      // the path now becomes the step's start: its geometry, the wraps it is followed from, the books that read it
      void Commit(Path &p)
      {
         const Config q = Now();
         Build(p, q, G);
         for (int k = 0; k < p.n; k++)
         {
            p.e[k].wrapRef = G.wrap[k];
            if (p.e[k].kind == kBallBack)
               books.backTouched = true;
         }
         Lcur = G.L;
         books.maxStretch = std::max(books.maxStretch, Stretch());
         // through: a ball not on the cord overlapping one of its runs (in 3D: within its clearance at its height)
         for (int i = 0; i < nb; i++)
         {
            const double cl = BallClear(z[i]);
            if (Find(p, kBall, i) >= 0 || cl <= 1e-6)
               continue;
            for (int k = 0; k < p.n; k++)
            {
               const V2 P0 = G.cen[k] + RightOf(G.t[k]) * G.rho[k];
               const V2 dc = c[i] - P0;
               const double u = std::min(std::max(dc.Dot(G.t[k]), 0.0), G.l[k]);
               if ((dc - G.t[k] * u).Length() < cl - 1e-6)
                  books.through = true;
            }
         }
      }
   };
}
