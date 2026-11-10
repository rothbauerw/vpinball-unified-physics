// license:GPLv3+

#include "core/stdafx.h"
#include "physics/playfieldcut.h"
#include "physics/hitrubber.h"
#include "physics/physicsini.h"
#include "physics/contactsystem.h"
#include "parts/pintable.h"
#include "core/player.h"
#include "ui/live/LiveUI.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

namespace
{
   using RubberContact::V3;
   using RubberContact::SurfaceTri;
   constexpr double U = ContactSystem::kMetresPerUnit;
   constexpr double kRampM = 7.5e-3;          // the drilled ramp: the outline narrows to a point this far from the end of the straight sides (user 2026-11-09; fixed, not a setting)
   constexpr double kEdgeRhoM = 0.25e-3;      // the rim's rounding radius (user 2026-11-09: ASSUMED, not measured)
   constexpr double kRampDepthM = 1.75e-3;    // the ramp floor's depth on the axis at the round end's edge (u = a), rising linearly to 0 at the tip (design section 2 revised: drill cross-section)
   constexpr double kThickM = 12.7e-3;        // the playfield's thickness: the through-hole's walls (1/2 in, user 2026-11-09)
   constexpr double kBandMarginM = 1.0e-3;    // the band reaches this far past the outline (flat there)
   constexpr double kTolM = 0.5e-6;           // the faces' largest distance from the shape on the top, rims and outline curves (design section 2: under the 0.8 um squash of a resting ball)
   constexpr int kFilletSeg = 16;             // facets on a rim's rounding (at most 90 deg: 5.6 deg each, 0.3 um inside the 0.25 mm circle)
   constexpr int kFloorSeg = 12;              // facets on each half of the ramp floor's arc (a ball reaches the floor only at the tip)
   constexpr double kMaxStepU = 0.2e-3;       // largest step between cross-sections along the ramp (the floor depth and rim radius change along it)
   constexpr double kWidthMaxMm = 5.03;       // width_mm default and largest: 0.198 in, the slot under a wire switch (PartData 7)
   constexpr double kWidthMinMm = 3.175;      // width_mm smallest: 1/8 in

   // The slot in its own frame (playfieldcut.h): a = half-width, hs = half the straight part, the round end to U1, the ramp from U1 (half-width wJ) to a point at L
   struct SlotShape
   {
      double a, hs, L, U1, wJ;
      SlotShape(const double a_, const double hs_) : a(a_), hs(hs_), L(kRampM)
      {
         U1 = L - std::sqrt(L * L - a * a);
         wJ = std::sqrt(a * a - U1 * U1);
      }
      double W(const double u) const   // the outline's half-width (where the top ends)
      {
         if (u <= 0.0)
            return a;
         if (u <= U1)
            return std::sqrt(std::max(0.0, a * a - u * u));
         if (u < L)
            return std::sqrt(std::max(0.0, wJ * wJ - 2.0 * U1 * (u - U1)));
         return 0.0;
      }
      double C(const double u) const { return u <= 0.0 ? a : u < a ? std::sqrt(a * a - u * u) : 0.0; }   // the through-hole's half-width
      bool Floor(const double u) const { return u > U1 && u < L; }
      double D(const double u) const { return kRampDepthM * (L - u) / (L - a); }   // the floor's depth on the axis
      double UofW(const double w) const { return U1 + (wJ * wJ - w * w) / (2.0 * U1); }   // on the ramp
   };

   struct PT
   {
      double t, z;
   };

   // Half a cross-section (t >= 0) at u, from the band's edge inward: the band edge, the rim's top tangent point, nf points round the rim (the last where it meets the wall or
   // floor), nfl points down the floor (to the hole's edge or the axis), the wall's foot. Always 3 + nf + nfl points (collapsed where a part is missing), so sections join face to face.
   void Profile(const SlotShape &sh, const double u, const double T, const int nf, const int nfl, std::vector<PT> &p)
   {
      p.clear();
      p.push_back({ T, 0.0 });
      const double w = sh.W(u);
      if (w <= 1e-12)   // the tip and past it: flat
      {
         p.insert(p.end(), (size_t)(2 + nf + nfl), PT { 0.0, 0.0 });
         return;
      }
      const double c = sh.C(u);
      const bool fl = sh.Floor(u);
      double r = kEdgeRhoM, D = 0.0, rho = 0.0, zc = 0.0;
      if (fl)
      {
         D = sh.D(u);
         r = std::min(r, D);   // the rounding is no deeper than the cut (near the tip)
         rho = (w * w + D * D) / (2.0 * D);   // the floor arc through (+-w, 0) and (0, -D): centre (0, zc)
         zc = rho - D;
      }
      // the rounding: a circle of radius r tangent to the top (centre (tf, -r)) and to the cut (the hole's wall or the floor's circle)
      const auto dA = [&](const double tf) { return c > 0.0 ? tf - c : 1e30; };
      const auto dB = [&](const double tf) { return fl ? std::hypot(tf, -r - zc) - rho : 1e30; };
      double lo = 0.0, hi = w + 2.0 * r + std::abs(zc) + 1e-3;
      for (int k = 0; k < 100; k++)
      {
         const double m = 0.5 * (lo + hi);
         (std::min(dA(m), dB(m)) < r ? lo : hi) = m;
      }
      const double tf = 0.5 * (lo + hi);
      const PT P { tf, -r };
      PT Q;
      bool onFloor = false;
      if (dA(tf) <= dB(tf))
         Q = { c, -r };
      else
      {
         const double l = std::hypot(tf, -r - zc);
         Q = { rho * tf / l, zc + rho * (-r - zc) / l };
         onFloor = true;
      }
      const double th1 = std::atan2(Q.z - P.z, Q.t - P.t);
      p.push_back({ tf, 0.0 });
      for (int i = 1; i <= nf; i++)
      {
         const double th = 0.5 * M_PI + (th1 - 0.5 * M_PI) * i / nf;
         p.push_back(i == nf ? Q : PT { P.t + r * std::cos(th), P.z + r * std::sin(th) });
      }
      if (onFloor && Q.t > c)
      {
         const double psQ = std::atan2(Q.z - zc, Q.t);
         const double psE = c > 0.0 ? std::atan2(-std::sqrt(std::max(0.0, rho * rho - c * c)), c) : -0.5 * M_PI;
         for (int i = 1; i <= nfl; i++)
         {
            const double ps = psQ + (psE - psQ) * i / nfl;
            p.push_back(i == nfl && c <= 0.0 ? PT { 0.0, -D } : PT { rho * std::cos(ps), zc + rho * std::sin(ps) });
         }
      }
      else
         p.insert(p.end(), (size_t)nfl, Q);
      const PT e = p.back();
      p.push_back(c > 0.0 ? PT { e.t, -kThickM } : e);
   }

   // The cross-sections past one end of the straight part (u from 0 to L, then the band's end): the outline's chord within kTolM of it
   std::vector<double> EndStations(const SlotShape &sh)
   {
      std::vector<double> u { 0.0 };
      const double th1 = std::asin(sh.U1 / sh.a);
      const int nr = std::max(1, (int)std::ceil(th1 * sh.a / std::sqrt(8.0 * kTolM * sh.a)));
      for (int k = 1; k <= nr; k++)
         u.push_back(k == nr ? sh.U1 : sh.a * std::sin(th1 * k / nr));
      for (double w = sh.wJ; w > 0.0;)   // the ramp, stepped in half-width: its curvature rises to 1 / U1 at the tip
      {
         const double q = w / sh.U1, kap = (1.0 / sh.U1) / std::pow(1.0 + q * q, 1.5);
         const double dw = std::min(std::sqrt(8.0 * kTolM / kap) / std::sqrt(1.0 + q * q), kMaxStepU * sh.U1 / w);
         w = std::max(0.0, w - dw);
         u.push_back(w > 0.0 ? sh.UofW(w) : sh.L);
      }
      u.push_back(sh.a);                                            // the hole closes
      u.push_back(sh.L - kEdgeRhoM * (sh.L - sh.a) / kRampDepthM);  // the rounding starts to follow the depth
      std::sort(u.begin(), u.end());
      u.erase(std::unique(u.begin(), u.end(), [](const double x, const double y) { return std::abs(x - y) < 1e-9; }), u.end());
      u.push_back(sh.L + kBandMarginM);
      return u;
   }

   void Notice(const std::string &msg)   // rule 8: on screen and in the log, once per table start (Setup runs once per start)
   {
      PLOGW << msg;
      if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(msg, 8000);
   }

   struct CutKeys
   {
      std::string name, shape;
      double x = 0.0, y = 0.0, angle = 0.0, lengthMm = -1.0, widthMm = kWidthMaxMm;
      bool xSet = false, ySet = false, lengthSet = false;
      std::vector<std::string> unknown;
   };

   std::vector<CutKeys> Parse(const std::filesystem::path &file)
   {
      std::vector<CutKeys> out;
      std::ifstream in(file);
      std::string line;
      CutKeys *cur = nullptr;
      while (std::getline(in, line))
      {
         const size_t c = line.find_first_of("#;");
         if (c != std::string::npos)
            line.resize(c);
         line = PhysicsIni::Trim(line);
         if (line.empty())
            continue;
         if (line.front() == '[' && line.back() == ']')
         {
            cur = nullptr;
            const std::string head = PhysicsIni::Trim(line.substr(1, line.size() - 2));
            if (PhysicsIni::Lower(head).rfind("playfieldcut:", 0) == 0)
            {
               out.emplace_back();
               cur = &out.back();
               cur->name = PhysicsIni::Trim(head.substr(13));
            }
            continue;
         }
         const size_t eq = line.find('=');
         if (cur == nullptr || eq == std::string::npos)
            continue;
         const std::string key = PhysicsIni::Lower(PhysicsIni::Trim(line.substr(0, eq)));
         const std::string sval = PhysicsIni::Lower(PhysicsIni::Trim(line.substr(eq + 1)));
         const double v = atof(sval.c_str());
         if (key == "shape") cur->shape = sval;
         else if (key == "x") { cur->x = v; cur->xSet = true; }
         else if (key == "y") { cur->y = v; cur->ySet = true; }
         else if (key == "angle") cur->angle = v;
         else if (key == "length_mm") { cur->lengthMm = v; cur->lengthSet = true; }
         else if (key == "width_mm") cur->widthMm = v;
         else cur->unknown.push_back(line);
      }
      return out;
   }

   // ---- The exact rim (physics, design section 2 third version) ----
   // The rim is a tube of radius kEdgeRhoM round a closed centre line at depth kEdgeRhoM under the top: the sharp outline offset outward by kEdgeRhoM everywhere (design section 2, REVISED
   // 11-10 after build v4). Straight part: t = a + rho. Round end: the concentric circle of radius a + rho, to the angle where the outline's round end meets the ramp (U1). Ramp: the outline
   // point plus rho along its outward normal; the parameter tau = sqrt(L - u) (outline u = L - tau^2, w = sqrt(2 U1) tau) runs smoothly through the tip, where the line reaches u = L + rho.
   // The outline is tangent-continuous and convex, so the line has no corner (its curvature steps at the straight part's end and at U1).
   constexpr int kRampSamples = 200;    // intervals of each ramp's centre line for the minima search, uniform in tau (about 20 um apart at the tip, 70 um near U1)
   constexpr int kRoundSamples = 8;     // intervals of each round end's arc (about 10 deg at 5.03 mm; a circle: g = (c - A).T changes sign at most once on it)
   constexpr double kRootTolM = 1e-12;  // a minimum is refined until it is known to this along the centre line (brief: 1e-12 m)

   struct Rim
   {
      const ContactSurface *host = nullptr;
      std::string name;
      double cx = 0, cy = 0, z0 = 0, lx = 0, ly = 1, px = -1, py = 0;   // the slot's frame in the world (m): centre, along (lx, ly), across (px, py), top z0
      double a = 0, hs = 0, L = kRampM, U1 = 0, wJ = 0;                  // half-width, half the straight part, the ramp's length, the round end's end (u), the half-width there
      double c = 0, tauJ = 0;                                            // outline w = c tau, tau at U1
      double rho = kEdgeRhoM, r = 0, th1 = 0, uRound = 0, uTip = 0;      // rim radius; centre line: the round end's radius a + rho, its end angle, its end u, the tip's u (L + rho)
      double lo[3] = {}, hi[3] = {};                                     // the centre line's world box
      struct Seg   // one interval of the centre line: piece, its parameter range, the sample at p0 (slot frame) and the unit tangent there, the arc length from the loop's start
      {
         int piece;
         double p0, p1, s, t, Ts, Tt, arc;
      };
      std::vector<Seg> seg;
      double total = 0.0;
      std::vector<double> rampU;   // the + end's ramp centre line at tau_j = tauJ j / kRampSamples: its u (decreasing in j), for the inside test
      double kapTip = 0.0, kapRound = 0.0;   // logged: the curvature (1/m, toward the hole) at the tip and on the round end

      // the ramp's centre line at tau (u past the end of the straight part, t across) and its derivative in tau: the outline point plus rho along the outline's outward normal
      void Ramp(const double tau, double &u, double &t, double &du, double &dt) const
      {
         const double m = std::sqrt(c * c + 4.0 * tau * tau), dm = 4.0 * tau / m;
         const double nu = c / m, nt = 2.0 * tau / m, dnu = -c * dm / (m * m), dnt = 2.0 / m - 2.0 * tau * dm / (m * m);   // the outline's outward normal
         u = L - tau * tau + rho * nu;
         t = c * tau + rho * nt;
         du = -2.0 * tau + rho * dnu;
         dt = c + rho * dnt;
      }
      // the closed centre line, piece by piece in loop order (the hole on the right): 0 side +, 1 + end round (side +), 2 + end ramp, 3 + end round (side -), 4 side -, 5 - end round (side -),
      // 6 - end ramp, 7 - end round (side +). Point (s, t) and its derivative in p.
      double PieceEnd(const int piece) const { return piece == 0 || piece == 4 ? 2.0 * hs : piece == 2 || piece == 6 ? 2.0 * tauJ : th1; }
      void Eval(const int piece, const double p, double &s, double &t, double &ds, double &dt) const
      {
         double u, w, du, dw;
         switch (piece)
         {
         case 0: s = -hs + p; t = r; ds = 1.0; dt = 0.0; break;
         case 1: s = hs + r * std::sin(p); t = r * std::cos(p); ds = r * std::cos(p); dt = -r * std::sin(p); break;
         case 2: Ramp(tauJ - p, u, w, du, dw); s = hs + u; t = w; ds = -du; dt = -dw; break;
         case 3: { const double th = th1 - p; s = hs + r * std::sin(th); t = -r * std::cos(th); ds = -r * std::cos(th); dt = -r * std::sin(th); } break;
         case 4: s = hs - p; t = -r; ds = -1.0; dt = 0.0; break;
         case 5: s = -hs - r * std::sin(p); t = -r * std::cos(p); ds = -r * std::cos(p); dt = r * std::sin(p); break;
         case 6: Ramp(-tauJ + p, u, w, du, dw); s = -(hs + u); t = w; ds = -du; dt = dw; break;
         default: { const double th = th1 - p; s = -hs - r * std::sin(th); t = r * std::cos(th); ds = r * std::cos(th); dt = r * std::sin(th); } break;
         }
      }
      // inside the centre line (slot frame): the region the flat top leaves out
      bool Inside(const double s, const double t) const
      {
         const double u = std::abs(s) - hs, T = std::abs(t);
         if (u <= 0.0)
            return T < r;
         if (u <= uRound)
            return T * T + u * u < r * r;
         if (u >= uTip)
            return false;
         // the ramp: the tau where the centre line's u is this u (decreasing in tau), by Newton inside the table's bracket
         const int K = kRampSamples;
         int j0 = 0, j1 = K;   // rampU[j0] >= u >= rampU[j1]
         while (j1 - j0 > 1)
         {
            const int jm = (j0 + j1) / 2;
            (rampU[jm] >= u ? j0 : j1) = jm;
         }
         double ta = tauJ * j0 / K, tb = tauJ * j1 / K, tau = ta + (tb - ta) * (rampU[j0] - u) / std::max(rampU[j0] - rampU[j1], 1e-300);
         double cu = 0.0, ct = 0.0, du, dt;
         for (int it = 0; it < 60; it++)
         {
            Ramp(tau, cu, ct, du, dt);
            const double f = cu - u;
            (f > 0.0 ? ta : tb) = tau;
            double next = du < 0.0 ? tau - f / du : 0.5 * (ta + tb);
            if (!(next > ta && next < tb))
               next = 0.5 * (ta + tb);
            if (std::abs(next - tau) * std::sqrt(du * du + dt * dt) < kRootTolM)
            {
               tau = next;
               break;
            }
            tau = next;
         }
         Ramp(tau, cu, ct, du, dt);
         return T < ct;
      }
      // a point of the playfield's faces that is not surface: inside the centre line, at the top or below it within the playfield's thickness (the through-hole; a thin playfield mesh's underside)
      bool Holds(const V3 &p) const
      {
         if (p.z > z0 + 1e-6 || p.z < z0 - kThickM - 1e-6)
            return false;
         const double dx = p.x - cx, dy = p.y - cy;
         return Inside(dx * lx + dy * ly, dx * px + dy * py);
      }
      double BoxDist(const V3 &q) const
      {
         const double v[3] = { q.x, q.y, q.z };
         double d2 = 0.0;
         for (int i = 0; i < 3; i++)
         {
            const double o = std::max(std::max(lo[i] - v[i], v[i] - hi[i]), 0.0);
            d2 += o * o;
         }
         return std::sqrt(d2);
      }
      // g(p) = (c - A(p)).T(p): > 0 while the distance falls along the loop, a minimum where it turns <= 0
      double G(const int piece, const double p, const double cs, const double ct) const
      {
         double s, t, ds, dt;
         Eval(piece, p, s, t, ds, dt);
         return ((cs - s) * ds + (ct - t) * dt) / std::sqrt(ds * ds + dt * dt);
      }
      // the signed curvature toward the hole at p (1/m): the unit tangent's turn over the parameter's neighbourhood within [pa, pb] (the piece's own side of a curvature step)
      double Curvature(const int piece, const double p, const double pa, const double pb) const
      {
         const double h = 1e-4 * (pb - pa), q0 = std::max(pa, p - h), q1 = std::min(pb, p + h);
         double s, t, ds0, dt0, ds1, dt1, ds, dt;
         Eval(piece, q0, s, t, ds0, dt0);
         Eval(piece, q1, s, t, ds1, dt1);
         Eval(piece, p, s, t, ds, dt);
         const double l0 = std::hypot(ds0, dt0), l1 = std::hypot(ds1, dt1), l = std::hypot(ds, dt);
         const double dTs = ds1 / l1 - ds0 / l0, dTt = dt1 / l1 - dt0 / l0;
         return (dTs * (dt / l) - dTt * (ds / l)) / (l * (q1 - q0));   // . n_in = (Tt, -Ts)
      }
   };

   // A rim contact: one true local minimum of the distance from the ball's centre to a cut's centre line, the rim's surface facing the ball (its upper inner part)
   struct RimHit
   {
      int cut;
      double delta, dist, reff, arc;   // depth R + rho - dist, the tube law's effective radius, the place on the loop (arc length)
      V3 point, normal;
   };
   // A local maximum of the same distance (the hump between two minima along the line): the place on the loop and the distance
   struct RimHump
   {
      int cut;
      double arc, dist;
   };

   std::vector<Rim> g_rims;   // this table's cuts, rebuilt at every table load (Setup)
   bool g_warnedMany = false; // three or more rim contacts pressing on one cut: logged once per table start

   // The root of g in [lo, hi] (g(lo) = glo and g(hi) = ghi of opposite signs, g > 0 counted with lo's side when glo > 0), by false position (Illinois) to kRootTolM along the line
   double Refine(const Rim &m, const int piece, double lo, double hi, double glo, double ghi, const double cs, const double ct)
   {
      if (ghi == 0.0)
         return hi;
      double sx, tx, dsx, dtx;
      m.Eval(piece, 0.5 * (lo + hi), sx, tx, dsx, dtx);
      const double speed = std::hypot(dsx, dtx) * 2.0;   // m per unit of p (an upper bound for the bracket's length)
      const bool loPos = glo > 0.0;
      double p = hi;
      int side = 0;
      for (int it = 0; it < 100; it++)
      {
         p = hi - ghi * (hi - lo) / (ghi - glo);
         if (!(p > lo && p < hi))
            p = 0.5 * (lo + hi);
         const double gp = m.G(piece, p, cs, ct);
         if ((gp > 0.0) == loPos)
         {
            lo = p;
            glo = gp;
            if (side == 1)
               ghi *= 0.5;
            side = 1;
         }
         else
         {
            hi = p;
            ghi = gp;
            if (side == -1)
               glo *= 0.5;
            side = -1;
         }
         if (gp == 0.0 || (hi - lo) * speed < kRootTolM)
            break;
      }
      return p;
   }

   // Every rim contact (true local minimum of the distance, facing the ball) of cut ci for a ball (radius R) at c, and every local maximum of the distance (humps): a sign change of
   // g = (c - A).T between an interval's two samples, refined inside it (the line is smooth, so g is continuous along it).
   void RimExtrema(const int ci, const double R, const V3 &c, std::vector<RimHit> &out, std::vector<RimHump> *humps)
   {
      const Rim &m = g_rims[ci];
      const double dx = c.x - m.cx, dy = c.y - m.cy, cs = dx * m.lx + dy * m.ly, ct = dx * m.px + dy * m.py, dz = c.z - (m.z0 - m.rho);
      const size_t n = m.seg.size();
      const auto arcOf = [&](const size_t i, const double p) {
         const Rim::Seg &q = m.seg[i];
         return q.arc + ((i + 1 < n ? m.seg[i + 1].arc : m.total) - q.arc) * (p - q.p0) / (q.p1 - q.p0);
      };
      const auto addHump = [&](const double arc, const double s, const double t) {
         if (humps != nullptr)
            humps->push_back({ ci, arc, std::sqrt((cs - s) * (cs - s) + (ct - t) * (ct - t) + dz * dz) });
      };
      const auto addMin = [&](const size_t i, const double p, const double s, const double t, const double Ts, const double Tt, const double kin) {
         if ((cs - s) * Tt - (ct - t) * Ts <= 0.0)   // the horizontal direction to the ball does not point into the hole: the rim's lower or outer side, not surface
            return;
         const V3 A(m.cx + m.lx * s + m.px * t, m.cy + m.ly * s + m.py * t, m.z0 - m.rho);
         const V3 dv = c - A;
         const double dist = dv.Length();
         if (dist <= 1e-12)
            return;
         RimHit h;
         h.cut = ci;
         h.point = A;
         h.dist = dist;
         h.normal = dv / dist;
         h.delta = R + m.rho - dist;
         h.arc = arcOf(i, p);
         // the tube law (rubbercontact.h ProbeTube): across the rim 1/R1 = 1/R + 1/rho; along it the centre line's curvature (TubeAlongCurvature), concave bends held within 20 ball radii
         double kLong = 0.0;
         const double nh = std::hypot(h.normal.x, h.normal.y);
         if (std::abs(kin) > 1e-9 && nh > 1e-9)
         {
            const double sg = kin > 0.0 ? 1.0 : -1.0, nis = Tt * sg, nit = -Ts * sg;   // toward the centre of curvature (slot frame)
            const double dnx = m.lx * nis + m.px * nit, dny = m.ly * nis + m.py * nit;
            kLong = RubberContact::TubeAlongCurvature(std::abs(kin), m.rho, -(h.normal.x * dnx + h.normal.y * dny) / nh);
         }
         const double R1 = 1.0 / (1.0 / R + 1.0 / m.rho), R2 = 1.0 / std::max(1.0 / R + kLong, 0.05 / R);
         h.reff = std::sqrt(R1 * R2);
         out.push_back(h);
      };
      for (size_t i = 0; i < n; i++)
      {
         const Rim::Seg &q = m.seg[i], &nx = m.seg[i + 1 < n ? i + 1 : 0];
         const double gA = (cs - q.s) * q.Ts + (ct - q.t) * q.Tt, gB = (cs - nx.s) * nx.Ts + (ct - nx.t) * nx.Tt;
         if ((gA > 0.0) != (gB > 0.0))   // a minimum (g from > 0 to <= 0) or a hump inside the interval
         {
            const double p = Refine(m, q.piece, q.p0, q.p1, gA, gB, cs, ct);
            double s, t, ds, dt;
            m.Eval(q.piece, p, s, t, ds, dt);
            if (gA > 0.0)
            {
               const double l = std::hypot(ds, dt);
               addMin(i, p, s, t, ds / l, dt / l, m.Curvature(q.piece, p, q.p0, q.p1));
            }
            else
               addHump(arcOf(i, p), s, t);
         }
      }
   }

   // the cuts of this host within reach of a ball at c
   void NearCuts(const ContactSurface *host, const V3 &c, const double reachM, std::vector<int> &out)
   {
      out.clear();
      for (size_t i = 0; i < g_rims.size(); i++)
         if (g_rims[i].host == host && g_rims[i].BoxDist(c) <= reachM)
            out.push_back((int)i);
   }

   // the faces whose nearest point to c is surface: one inside a near cut's centre line (at or under the top) is not (as CarvedHole::Holds for a turntable)
   void FacesLeft(const std::vector<SurfaceTri> &tris, const std::vector<int> &faces, const std::vector<int> &cuts, const V3 &c, std::vector<int> &out)
   {
      out.clear();
      for (const int i : faces)
      {
         int kind;
         const V3 p = RubberContact::ClosestOnTriangle(tris[i], c, kind);
         if (std::none_of(cuts.begin(), cuts.end(), [&](const int ci) { return g_rims[ci].Holds(p); }))
            out.push_back(i);
      }
   }

   // build cut m's centre line samples and log numbers (slot frame set, a, hs, rho)
   void BuildRim(Rim &m)
   {
      m.U1 = m.L - std::sqrt(m.L * m.L - m.a * m.a);
      m.wJ = std::sqrt(m.a * m.a - m.U1 * m.U1);
      m.c = std::sqrt(2.0 * m.U1);
      m.tauJ = std::sqrt(m.L - m.U1);
      m.r = m.a + m.rho;
      m.th1 = std::asin(m.U1 / m.a);
      m.uRound = m.r * std::sin(m.th1);
      m.uTip = m.L + m.rho;
      const int K = kRampSamples;
      m.rampU.assign(K + 1, 0.0);
      for (int j = 0; j <= K; j++)
      {
         double u, t, du, dt;
         m.Ramp(m.tauJ * j / K, u, t, du, dt);
         m.rampU[j] = u;
      }
      m.seg.clear();
      for (int piece = 0; piece < 8; piece++)
      {
         const double P = m.PieceEnd(piece);
         if (P <= 0.0)
            continue;
         std::vector<double> ps;
         const int nI = piece == 0 || piece == 4 ? 1 : piece == 2 || piece == 6 ? 2 * K : kRoundSamples;
         for (int i = 0; i <= nI; i++)
            ps.push_back(P * i / nI);
         for (size_t i = 0; i + 1 < ps.size(); i++)
         {
            Rim::Seg q;
            q.piece = piece;
            q.p0 = ps[i];
            q.p1 = ps[i + 1];
            double ds, dt;
            m.Eval(piece, q.p0, q.s, q.t, ds, dt);
            const double l = std::hypot(ds, dt);
            q.Ts = ds / l;
            q.Tt = dt / l;
            q.arc = m.total;
            double s0 = q.s, t0 = q.t;
            for (int j = 1; j <= 16; j++)   // arc length by chords
            {
               double s1, t1;
               m.Eval(piece, q.p0 + (q.p1 - q.p0) * j / 16, s1, t1, ds, dt);
               m.total += std::hypot(s1 - s0, t1 - t0);
               s0 = s1;
               t0 = t1;
            }
            m.seg.push_back(q);
         }
      }
      for (int i = 0; i < 3; i++)
      {
         m.lo[i] = 1e30;
         m.hi[i] = -1e30;
      }
      for (const Rim::Seg &q : m.seg)
      {
         const double v[3] = { m.cx + m.lx * q.s + m.px * q.t, m.cy + m.ly * q.s + m.py * q.t, m.z0 - m.rho };
         for (int i = 0; i < 3; i++)
         {
            m.lo[i] = std::min(m.lo[i], v[i]);
            m.hi[i] = std::max(m.hi[i], v[i]);
         }
      }
      // logged: the curvature at the tip and on the round end
      m.kapTip = m.Curvature(2, m.tauJ, m.tauJ - 1e-3 * m.tauJ, m.tauJ + 1e-3 * m.tauJ);
      m.kapRound = 1.0 / m.r;
   }
}

namespace PlayfieldCutMesh
{
   void BuildSlot(const Slot &s, std::vector<SurfaceTri> &out, Stats &st)
   {
      out.clear();
      st = Stats();
      const SlotShape sh(s.halfWidth, s.halfStraight);
      const std::vector<double> us = EndStations(sh);
      std::vector<double> ss, uu;   // the cross-sections along the slot
      for (size_t k = us.size(); k-- > 0;)
      {
         ss.push_back(-(sh.hs + us[k]));
         uu.push_back(us[k]);
      }
      for (const double u : us)
      {
         if (sh.hs + u - ss.back() < 1e-12)   // no straight part: the two ends share the middle section
            continue;
         ss.push_back(sh.hs + u);
         uu.push_back(u);
      }
      const double T = sh.a + kBandMarginM;
      st.sBand = sh.hs + sh.L + kBandMarginM;
      st.tBand = T;
      st.stations = (int)ss.size();
      const double px = -s.ly, py = s.lx;   // across
      const auto world = [&](const double sv, const double t, const double z) { return V3(s.cx + s.lx * sv + px * t, s.cy + s.ly * sv + py * t, s.z + z); };
      const auto strip = [&](const double s0, const std::vector<PT> &p0, const double s1, const std::vector<PT> &p1, const int side, std::vector<SurfaceTri> &o) {
         const auto add = [&](const V3 &a, const V3 &b, const V3 &c) {
            if ((b - a).Cross(c - a).Length() > 1e-18)   // a part of the section that is missing here gives zero-area faces: left out
               o.push_back(RubberContact::MakeSurfaceTri(a, b, c));
         };
         for (size_t j = 0; j + 1 < p0.size(); j++)
         {
            const V3 A = world(s0, side * p0[j].t, p0[j].z), B = world(s0, side * p0[j + 1].t, p0[j + 1].z);
            const V3 C = world(s1, side * p1[j + 1].t, p1[j + 1].z), D = world(s1, side * p1[j].t, p1[j].z);
            add(A, B, C);
            add(A, C, D);
         }
      };
      std::vector<std::vector<PT>> prof(ss.size());
      for (size_t i = 0; i < ss.size(); i++)
         Profile(sh, uu[i], T, kFilletSeg, kFloorSeg, prof[i]);
      for (size_t i = 0; i + 1 < ss.size(); i++)
         for (int side = 1; side >= -1; side -= 2)
            strip(ss[i], prof[i], ss[i + 1], prof[i + 1], side, out);
      st.faces = (int)out.size();
      // check: the true shape midway between sections, against the faces of that strip (one side; the other is its mirror)
      std::vector<PT> dense;
      std::vector<SurfaceTri> one;
      constexpr int nfD = 128, nflD = 64;
      for (size_t i = 0; i + 1 < ss.size(); i++)
      {
         const double sm = 0.5 * (ss[i] + ss[i + 1]);
         Profile(sh, std::abs(sm) - sh.hs, T, nfD, nflD, dense);
         one.clear();
         strip(ss[i], prof[i], ss[i + 1], prof[i + 1], 1, one);
         for (size_t q = 0; q + 1 < dense.size(); q++)   // (not the wall's foot)
         {
            const V3 P = world(sm, dense[q].t, dense[q].z);
            double best = 1e30;
            for (const SurfaceTri &t : one)
            {
               int kind;
               best = std::min(best, (P - RubberContact::ClosestOnTriangle(t, P, kind)).Length());
            }
            double &dev = q <= (size_t)(1 + nfD) ? st.devRim : st.devFloor;
            dev = std::max(dev, best);
         }
      }
   }

}

namespace PlayfieldCuts
{
   void Setup(PhysicsEngine *, PinTable *const table)
   {
      g_rims.clear();
      g_warnedMany = false;
      const std::filesystem::path file = PhysicsIni::FindUnified(table->m_filename);
      if (file.empty())
         return;
      const PhysicsIni::Features feat = PhysicsIni::Active(file);
      std::vector<CutKeys> cuts = Parse(file);
      if (cuts.empty())
         return;
      if (!feat.playfieldCuts)
      {
         PLOGI << "Playfield cuts: " << cuts.size() << " [PlayfieldCut:<name>] section(s) ignored: [Features] playfield_cuts = 0";
         return;
      }
      std::vector<HitRubber *> pf;   // the playfield's hosts (the plane's, or a playfield_mesh primitive's): a host with faces has a height under the cut
      for (ContactHost *const h : ContactHost::AllHosts())
         if (h->IsPlayfieldKind())
            if (HitRubber *const r = dynamic_cast<HitRubber *>(h); r != nullptr)
               pf.push_back(r);
      if (pf.empty())
      {
         Notice("Playfield cuts: " + std::to_string(cuts.size()) + " [PlayfieldCut:<name>] section(s) ignored: cuts need the playfield on the force model ([Features] playfield_contact = 1)");
         return;
      }
      const double R = 25.0 * U;   // a standard ball, for the logged paper numbers
      for (const CutKeys &k : cuts)
      {
         const std::string where = "Playfield cut " + k.name + ": ";
         for (const std::string &line : k.unknown)
            Notice(where + "unknown key, ignored line '" + line + "' (physics ini)");
         if (k.shape != "slot")
         {
            Notice(where + (k.shape.empty() ? std::string("no shape") : "shape '" + k.shape + "' is not known") + ": the cut is ignored. The only shape is slot (physics ini)");
            continue;
         }
         if (!k.xSet || !k.ySet || !k.lengthSet || k.lengthMm < 0.0)
         {
            Notice(where + "needs x, y and length_mm (0 or more): the cut is ignored (physics ini)");
            continue;
         }
         double widthMm = k.widthMm;
         if (widthMm > kWidthMaxMm || widthMm < kWidthMinMm)
         {
            const double lim = widthMm > kWidthMaxMm ? kWidthMaxMm : kWidthMinMm;
            std::ostringstream o;
            o << where << "width_mm " << widthMm << " is outside " << kWidthMinMm << " .. " << kWidthMaxMm << ": using " << lim << " (physics ini)";
            Notice(o.str());
            widthMm = lim;
         }
         HitRubber *host = nullptr;
         double zU = 0.0;
         for (HitRubber *const h : pf)
         {
            double z;
            if (h->SurfaceHeightAt(k.x, k.y, z) && (host == nullptr || z > zU))
            {
               host = h;
               zU = z;
            }
         }
         if (host == nullptr)
         {
            Notice(where + "no playfield surface at its centre: the cut is ignored (physics ini x, y)");
            continue;
         }
         Rim m;
         m.host = host;
         m.name = k.name;
         m.cx = k.x * U;
         m.cy = k.y * U;
         m.z0 = zU * U;
         const double rot = k.angle * (M_PI / 180.0);
         m.lx = -std::sin(rot);
         m.ly = std::cos(rot);
         m.px = -m.ly;
         m.py = m.lx;
         m.hs = 0.5 * k.lengthMm * 1e-3;
         m.a = 0.5 * widthMm * 1e-3;
         // the playfield must be flat over the cut and 1 mm round it (the rim's top tangent line is the top's edge): its top height on a 0.5 mm grid
         bool flat = true;
         const double sB = m.hs + m.L + 1e-3, tB = m.a + m.rho + 1e-3;
         for (double s = -sB; flat && s <= sB + 1e-9; s += 0.5e-3)
            for (double t = -tB; flat && t <= tB + 1e-9; t += 0.5e-3)
            {
               double z;
               flat = host->SurfaceHeightAt((m.cx + m.lx * s + m.px * t) / U, (m.cy + m.ly * s + m.py * t) / U, z) && std::abs(z - zU) * U <= 1e-6;
            }
         if (!flat)
         {
            Notice(where + "the playfield is not flat over the cut (a raised or curved playfield, or its edge): the cut is ignored (physics ini x, y)");
            continue;
         }
         BuildRim(m);
         const double ar = m.a + m.rho;
         const double hc = std::sqrt((R + m.rho) * (R + m.rho) - ar * ar) - m.rho;   // the ball's centre over the middle of the straight part, resting on both rims
         std::ostringstream o;
         o << "Playfield cut " << k.name << ": slot in the playfield " << host->Name() << " (exact rim), centre " << k.x << ", " << k.y << " (playfield z " << zU << "), angle " << k.angle
           << " deg (along " << m.lx << ", " << m.ly << "), straight part " << k.lengthMm << " mm, width " << widthMm << " mm; round end to u " << m.U1 * 1e3 << " mm, ramp to " << m.L * 1e3
           << " mm, rim radius " << m.rho * 1e3 << " mm (assumed). Rim centre line (outline + rho): " << m.seg.size() << " intervals, " << m.total * 1e3
           << " mm round; curvature toward the hole " << m.kapRound << " /m on the round end, " << m.kapTip << " /m at the tip. Paper, 27 mm ball: sinks " << (R - hc) * 1e3
           << " mm, rim normals " << std::asin(ar / (R + m.rho)) * 180.0 / M_PI << " deg from vertical";
         PLOGI << o.str();
         g_rims.push_back(std::move(m));
      }
   }

   bool Any(const ContactSurface *const host)
   {
      return std::any_of(g_rims.begin(), g_rims.end(), [&](const Rim &m) { return m.host == host; });
   }

   RubberContact::TubeProbe Probe(const ContactSurface *const host, const RubberContact::BallSI &bs, const std::vector<SurfaceTri> &tris, const std::vector<int> &faces, const double edgeRho,
                                  const V3 &c, const double reachM, const RubberContact::CarvedHole *const holes, const int nHoles, const bool patch, const double zero)
   {
      static std::vector<int> cuts, left;   // (the physics thread's; reused)
      static std::vector<RimHit> hits;
      static std::vector<RimHump> humps;
      NearCuts(host, c, reachM, cuts);
      if (cuts.empty())   // today's probe
      {
         RubberContact::TubeProbe pr = RubberContact::ProbeTriangles(bs, tris, faces, edgeRho, c, holes, nHoles);
         if (patch)
            RubberContact::PatchTriangles(bs, tris, faces, c, zero, pr, holes, nHoles);
         return pr;
      }
      FacesLeft(tris, faces, cuts, c, left);
      RubberContact::TubeProbe pr = RubberContact::ProbeTriangles(bs, tris, left, edgeRho, c, holes, nHoles);
      hits.clear();
      humps.clear();
      for (const int ci : cuts)
         RimExtrema(ci, bs.R, c, hits, patch ? &humps : nullptr);
      int best = -1;   // the deepest: a face (-1) or a rim contact
      double bestDelta = pr.seg >= 0 ? pr.delta : -1e30;
      for (size_t j = 0; j < hits.size(); j++)
         if (hits[j].delta > bestDelta)
         {
            bestDelta = hits[j].delta;
            best = (int)j;
         }
      if (best < 0)   // the flat (or nothing): the faces' fold sum as today; a rim facing a ball whose centre is outside the centre line is out of its reach
      {
         if (patch)
            RubberContact::PatchTriangles(bs, tris, left, c, zero, pr, holes, nHoles);
         return pr;
      }
      const RimHit &b = hits[best];
      pr = RubberContact::TubeProbe();
      pr.point = b.point;
      pr.normal = b.normal;
      pr.dist = b.dist;
      pr.delta = b.delta;
      pr.overlap = b.delta > 0.0;
      pr.reff = b.reff;
      pr.rho = g_rims[b.cut].rho;
      pr.seg = (int)tris.size() + b.cut;
      if (!patch || !pr.overlap || pr.delta <= zero)
         return pr;
      // The patch sum (design section 2, REVISED 11-10): every rim contact pressing the ball (d = delta - zero > 0) adds its Hertz energy, s25 = sum w d^2.5, gradient sum w d^1.5 n
      // (rubbercontact.h FoldSum), with a weight from the hump between it and each other contact of the same cut. For two contacts 1, 2 of one cut:
      //   ds = R + rho - (the hump's distance) - zero   how deep the hump (saddle) between them presses: the smaller of the largest distances along the two arcs of the line between them
      //   s  = clamp(ds / min(d1, d2), 0, 1)           0: free between (two contacts), 1: the shallower one is on the hump (one contact)
      //   p_i = d_i - ds,  w_i = (1 - s) + s p_i / (p1 + p2)
      // Symmetric in 1 and 2 (no contact fixed at weight 1): a ball on the axis gets no sideways force from which one is deeper. A minimum born with its hump (p = 0) comes in at weight 0, and
      // two minima merging leave one at full weight, so contacts come and go with no jump. Contacts on different cuts: weight 1 each. Three or more on one cut (not expected on the slot: one
      // minimum per side; logged once as a WARN with the ball's position): w_i is the product of its pair weights with each other one. The general rule for several contacts on one object
      // is to-do 1.33.
      static std::vector<double> wts;
      RubberContact::FoldSum sum;
      int nPress = 0;
      for (const RimHit &h : hits)
         nPress += h.delta - zero > 0.0;
      wts.assign(hits.size(), 1.0);
      for (const int ci : cuts)
      {
         int nCut = 0;
         for (const RimHit &h : hits)
            nCut += h.cut == ci && h.delta - zero > 0.0;
         if (nCut < 2)
            continue;
         if (nCut >= 3 && !g_warnedMany)
         {
            g_warnedMany = true;
            PLOGW << "Playfield cut " << g_rims[ci].name << ": " << nCut << " rim contacts on one cut at ball x " << c.x / U << " y " << c.y / U << " z " << c.z / U
                  << " (VPX units; once per table start): weighted pair by pair";
         }
         const double tot = g_rims[ci].total, top = bs.R + g_rims[ci].rho - zero;
         for (size_t i = 0; i < hits.size(); i++)
            for (size_t j = i + 1; j < hits.size(); j++)
            {
               const RimHit &h1 = hits[i], &h2 = hits[j];
               if (h1.cut != ci || h2.cut != ci || h1.delta - zero <= 0.0 || h2.delta - zero <= 0.0)
                  continue;
               // the largest distance on each arc between them (forward from h1 to h2, and the rest); the hump is the lower of the two
               const double fwd = std::fmod(h2.arc - h1.arc + 2.0 * tot, tot);
               double m1 = -1.0, m2 = -1.0;
               for (const RimHump &q : humps)
               {
                  if (q.cut != ci)
                     continue;
                  const double x = std::fmod(q.arc - h1.arc + 2.0 * tot, tot);
                  double &mm = x < fwd ? m1 : m2;
                  mm = std::max(mm, q.dist);
               }
               const double d1 = h1.delta - zero, d2 = h2.delta - zero;
               const double hump = m1 < 0.0 ? m2 : m2 < 0.0 ? m1 : std::min(m1, m2);
               const double ds = hump < 0.0 ? -1.0 : top - hump;   // (no hump found: free between)
               const double s = std::min(std::max(ds / std::min(d1, d2), 0.0), 1.0), p1 = d1 - ds, p2 = d2 - ds;
               const double f1 = p1 + p2 > 0.0 ? p1 / (p1 + p2) : 0.5;
               wts[i] *= (1.0 - s) + s * f1;
               wts[j] *= (1.0 - s) + s * (1.0 - f1);
            }
      }
      if (nPress < 2)
         return pr;
      for (size_t j = 0; j < hits.size(); j++)
      {
         const double d = hits[j].delta - zero;
         if (d <= 0.0 || wts[j] <= 0.0)
            continue;
         const double d15 = d * std::sqrt(d);
         sum.s25 += wts[j] * d * d15;
         sum.g = sum.g + hits[j].normal * (wts[j] * d15);
      }
      if (sum.s25 > 0.0)
         sum.Finish(pr, zero);
      return pr;
   }

   double ArrivalS(const ContactSurface *const host, const RubberContact::BallSI &bs, const std::vector<SurfaceTri> &tris, const std::vector<int> &faces, const int nearestSeg, const V3 &c,
                   const V3 &v, const double a, const double hMax, const double touchM, const double arriveM, const double reachM)
   {
      static std::vector<int> cuts, left;
      static std::vector<RimHit> hits;
      NearCuts(host, c, reachM, cuts);
      if (cuts.empty())   // today's bound
         return RubberContact::ArrivalTriangles(bs, tris, faces, nearestSeg, c, v, a, hMax, touchM, arriveM);
      const int nT = (int)tris.size();
      FacesLeft(tris, faces, cuts, c, left);
      double h = RubberContact::ArrivalTriangles(bs, tris, left, nearestSeg < nT ? nearestSeg : -1, c, v, a, hMax, touchM, arriveM);
      for (const int ci : cuts)
      {
         hits.clear();
         RimExtrema(ci, bs.R, c, hits, nullptr);
         int skip = -1;   // the probe's nearest on this cut: its deepest (the probe takes the first of equals)
         if (nearestSeg == nT + ci)
            for (size_t j = 0; j < hits.size(); j++)
               if (skip < 0 || hits[j].delta > hits[skip].delta)
                  skip = (int)j;
         for (size_t j = 0; j < hits.size(); j++)
         {
            if ((int)j == skip)
               continue;
            const double gap = hits[j].dist - (bs.R + g_rims[ci].rho);
            if (gap <= touchM)   // touching: the step's own law has it
               continue;
            h = RubberContact::ArrivalStep(gap - arriveM, std::max(-v.Dot(hits[j].normal), 0.0), a, h);
         }
      }
      return h;
   }
}
