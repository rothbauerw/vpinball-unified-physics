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
   constexpr double kWidthMaxMm = 5.03;       // width_mm default and largest: 0.198 in, the slot under a wire switch (PartData 7)
   constexpr double kWidthMinMm = 3.175;      // width_mm smallest: 1/8 in

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
   constexpr int kPieces = 8;           // the centre line's pieces (Rim::Eval); a place on the loop is piece + its fraction of the piece, 0 .. kPieces
   constexpr int kRampSamples = 200;    // the inside test's table of the ramp's centre line, uniform in tau (Rim::rampU)
   constexpr double kRootTolM = 1e-12;  // the inside test's ramp point is found to this along the centre line (brief: 1e-12 m)
   constexpr double kReachPadM = 1e-9;  // added to every reach: a point is passed over only when it is surely beyond it (rounding is ~1e-17 m)
   constexpr double kNoLimitM = 1e30;   // a reach that passes over nothing

   struct Rim
   {
      const ContactSurface *host = nullptr;
      std::string name;
      double cx = 0, cy = 0, z0 = 0, lx = 0, ly = 1, px = -1, py = 0;   // the slot's frame in the world (m): centre, along (lx, ly), across (px, py), top z0
      double a = 0, hs = 0, L = kRampM, U1 = 0, wJ = 0;                  // half-width, half the straight part, the ramp's length, the round end's end (u), the half-width there
      double c = 0, tauJ = 0;                                            // outline w = c tau, tau at U1
      double rho = kEdgeRhoM, r = 0, th1 = 0, uRound = 0, uTip = 0;      // rim radius; centre line: the round end's radius a + rho, its end angle, its end u, the tip's u (L + rho)
      double lo[3] = {}, hi[3] = {};                                     // the centre line's world box
      double js[kPieces] = {}, jt[kPieces] = {}, jTs[kPieces] = {}, jTt[kPieces] = {};   // the junctions: the start of each piece (slot frame) and the unit tangent there
      double total = 0.0;          // logged: the centre line's length
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
      // the ramp's centre line's curvature toward the hole at tau (1/m): the outline's radius of curvature (c^2 + 4 tau^2)^1.5 / (2 c) plus rho (the offset shares its centres of curvature)
      double RampCurvature(const double tau) const
      {
         const double m = std::sqrt(c * c + 4.0 * tau * tau);
         return 1.0 / (m * m * m / (2.0 * c) + rho);
      }
      // the closed centre line, piece by piece in loop order (the hole on the right): 0 side +, 1 + end round (side +), 2 + end ramp, 3 + end round (side -), 4 side -, 5 - end round (side -),
      // 6 - end ramp, 7 - end round (side +); p runs from 0 to 2 hs on a side, 2 tauJ on a ramp, th1 on a round end. Point (s, t) and its derivative in p.
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
   };

   // A rim contact: one true local minimum of the distance from the ball's centre to a cut's centre line, the rim's surface facing the ball (its upper inner part)
   struct RimHit
   {
      int cut;
      double delta, dist, reff, loop;   // depth R + rho - dist, the tube law's effective radius, the place on the loop (piece + fraction, 0 .. kPieces)
      V3 point, normal;
   };
   // A local maximum of the same distance (the hump between two minima along the line): the place on the loop and the distance
   struct RimHump
   {
      int cut;
      double loop, dist;
   };

   std::vector<Rim> g_rims;   // this table's cuts, rebuilt at every table load (Setup)
   const PinTable *g_rimsTable = nullptr;   // the table they were read for (the drawn playfield of another table is not cut)
   bool g_warnedMany = false; // three or more rim contacts pressing on one cut: logged once per table start

   // the ball's centre in cut m's frame: along, across, height over the centre line
   void SlotFrame(const Rim &m, const V3 &c, double &cs, double &ct, double &dz)
   {
      const double dx = c.x - m.cx, dy = c.y - m.cy;
      cs = dx * m.lx + dy * m.ly;
      ct = dx * m.px + dy * m.py;
      dz = c.z - (m.z0 - m.rho);
   }

   // The real roots of x^3 + P x + Q = 0, ascending, and how many: three (D < 0) by Viete's trigonometric form, one by Cardano's formula on its branch with no cancellation
   int CubicRoots(const double P, const double Q, double x[3])
   {
      const double D = 0.25 * Q * Q + P * P * P / 27.0;
      if (D < 0.0)   // (so P < 0)
      {
         const double m = 2.0 * std::sqrt(-P / 3.0), th = std::acos(std::clamp(1.5 * Q / P * std::sqrt(-3.0 / P), -1.0, 1.0)) / 3.0;
         x[0] = m * std::cos(th - 4.0 * M_PI / 3.0);
         x[1] = m * std::cos(th - 2.0 * M_PI / 3.0);
         x[2] = m * std::cos(th);
         return 3;
      }
      const double t1 = -std::copysign(std::cbrt(0.5 * std::abs(Q) + std::sqrt(D)), Q);
      x[0] = t1 != 0.0 ? t1 - P / (3.0 * t1) : 0.0;
      return 1;
   }

   // Every rim contact (true local minimum of the distance, facing the ball, within reachMin of the centre) of cut ci for a ball (radius R) at c, and every local maximum of the distance
   // (humps, when humps is given), in loop order. The distance's extrema along the line are its normal foot points: where the line's normal passes through the ball's centre (the line is
   // level, so the 3D distance's extrema are the plane's). g = (c - A).T > 0 while the distance falls along the loop: a minimum where g turns from > 0 to <= 0, a hump from <= 0 to > 0.
   // The line is tangent-continuous, so g is continuous along it: its sign at the junctions, shared by the two pieces that meet there, says whether a piece holds a root (minima and humps
   // alternate round the loop), and each root is placed in closed form:
   //   side        the foot of the normal (s = the ball's s); a line has no maximum
   //   round end   a circle (at most one root on an arc shorter than pi): the nearest point O + r e (a minimum) or the farthest O - r e (a hump), e the unit vector from its centre O
   //               to the ball's; the ball's centre on O (every point equally far): the arc's start
   //   ramp        the rim line is the outline offset along its normals, so they share their normal lines: the outline's foot points. In the end's frame (u along, w across; outline
   //               u = L - tau^2, w = c tau) with the ball's centre at (u0, w0), d/dtau of half the squared distance is F(tau) = 2 tau^3 + q tau - c w0, q = c^2 - 2 (L - u0), times a
   //               positive factor (g has the sign of F on the + end, of -F on the - end). F' = 6 tau^2 + q splits the ramp at tau = +-sqrt(-q / 6) (q < 0) into monotone pieces,
   //               each holding one root where the signs at its ends differ: the cubic's root in it.
   void RimExtrema(const int ci, const double R, const V3 &c, std::vector<RimHit> &out, std::vector<RimHump> *humps, const double reachMin)
   {
      const Rim &m = g_rims[ci];
      double cs, ct, dz;
      SlotFrame(m, c, cs, ct, dz);
      // one extremum: the point (s, t), the unit tangent along the loop, the curvature toward the hole (1/m)
      const auto add = [&](const bool isMin, const double loop, const double s, const double t, const double Ts, const double Tt, const double kin) {
         if (!isMin)
         {
            if (humps != nullptr)
               humps->push_back({ ci, loop, std::sqrt((cs - s) * (cs - s) + (ct - t) * (ct - t) + dz * dz) });
            return;
         }
         if ((cs - s) * Tt - (ct - t) * Ts <= 0.0)   // the horizontal direction to the ball does not point into the hole: the rim's lower or outer side, not surface
            return;
         const V3 A(m.cx + m.lx * s + m.px * t, m.cy + m.ly * s + m.py * t, m.z0 - m.rho);
         const V3 dv = c - A;
         const double dist = dv.Length();
         if (dist <= 1e-12 || dist > reachMin)   // (beyond the caller's reach: it cannot change the caller's answer)
            return;
         RimHit h;
         h.cut = ci;
         h.point = A;
         h.dist = dist;
         h.normal = dv / dist;
         h.delta = R + m.rho - dist;
         h.loop = loop;
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
      bool pos[kPieces];   // g > 0 at each junction
      for (int k = 0; k < kPieces; k++)
         pos[k] = (cs - m.js[k]) * m.jTs[k] + (ct - m.jt[k]) * m.jTt[k] > 0.0;
      for (int k = 0; k < kPieces; k++)
      {
         const int k1 = (k + 1) % kPieces;
         const bool pa = pos[k], pb = pos[k1];   // a root in the piece: a minimum if g > 0 at its start
         if (k == 0 || k == 4)   // a side, t = +-r along +-s
         {
            if (pa == pb)
               continue;
            const double sg = k == 0 ? 1.0 : -1.0, len = 2.0 * m.hs, p = std::clamp(m.hs + sg * cs, 0.0, len);
            add(pa, k + (len > 0.0 ? p / len : 0.0), sg * (p - m.hs), sg * m.r, sg, 0.0, 0.0);
         }
         else if (k & 1)   // a round end: the circle of radius r about O = (+-hs, 0), clockwise from junction k to junction k1
         {
            if (pa == pb)
               continue;
            const double O = k < 4 ? m.hs : -m.hs, as = m.js[k] - O, at = m.jt[k], bs = m.js[k1] - O, bt = m.jt[k1];   // the arc's ends seen from O
            double es = cs - O, et = ct;
            const double n = std::hypot(es, et);
            if (n > 0.0)
            {
               es *= (pa ? m.r : -m.r) / n;   // r e: the nearest point (a minimum) or the farthest (a hump), from O
               et *= (pa ? m.r : -m.r) / n;
            }
            if (n == 0.0 || as * et - at * es > 0.0)   // the ball's centre on O, or the point before the arc's start (rounding at the junction): the start
               add(pa, k, m.js[k], m.jt[k], m.jTs[k], m.jTt[k], 1.0 / m.r);
            else if (es * bt - et * bs > 0.0)          // past its end
               add(pa, k + 1, m.js[k1], m.jt[k1], m.jTs[k1], m.jTt[k1], 1.0 / m.r);
            else   // (the place on the loop: sin of the angle from the start over sin th1, which orders the same way as the angle)
               add(pa, k + (as * et - at * es) / (as * bt - at * bs), O + es, et, et / m.r, -es / m.r, 1.0 / m.r);
         }
         else   // a ramp: k 2 the + end (tau from tauJ down to -tauJ along the loop), k 6 the - end (tau from -tauJ up to tauJ)
         {
            const double sg = k == 2 ? 1.0 : -1.0, q = m.c * m.c - 2.0 * (m.L - (sg * cs - m.hs)), r0 = m.c * ct;
            const auto F = [&](const double x) { return 2.0 * x * x * x + q * x - r0; };
            double tau[4], x[3];   // the monotone pieces' ends in loop order; the cubic's roots
            bool gp[4];            // g > 0 there
            int nb = 0, nr = 0;
            tau[nb] = sg * m.tauJ;
            gp[nb++] = pa;
            const double tc = q < 0.0 ? std::sqrt(-q / 6.0) : m.tauJ;   // F's turning points
            if (tc < m.tauJ)
               for (const double xc : { sg * tc, -sg * tc })
               {
                  tau[nb] = xc;
                  gp[nb++] = sg * F(xc) > 0.0;
               }
            tau[nb] = -sg * m.tauJ;
            gp[nb++] = pb;
            for (int i = 0; i + 1 < nb; i++)
            {
               if (gp[i] == gp[i + 1])
                  continue;
               if (nr == 0)
                  nr = CubicRoots(0.5 * q, -0.5 * r0, x);
               // the root in this piece (of rounding's several there, the one nearest F = 0; with none there, the nearest, at the piece's end)
               const double lo = std::min(tau[i], tau[i + 1]), hi = std::max(tau[i], tau[i + 1]);
               int j = -1;
               for (int jj = 0; jj < nr; jj++)
                  if (x[jj] >= lo && x[jj] <= hi && (j < 0 || std::abs(F(x[jj])) < std::abs(F(x[j]))))
                     j = jj;
               double xr;
               if (j >= 0)
                  xr = x[j];
               else
               {
                  j = 0;
                  for (int jj = 1; jj < nr; jj++)
                     if (std::min(std::abs(x[jj] - lo), std::abs(x[jj] - hi)) < std::min(std::abs(x[j] - lo), std::abs(x[j] - hi)))
                        j = jj;
                  xr = std::clamp(x[j], lo, hi);
               }
               double u, w, du, dw;
               m.Ramp(xr, u, w, du, dw);
               const double l = std::hypot(du, dw);
               add(gp[i], k + (m.tauJ - sg * xr) / (2.0 * m.tauJ), sg * (m.hs + u), w, -du / l, -sg * dw / l, m.RampCurvature(xr));
            }
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

   // build cut m's centre line: its numbers, junctions, world box and the inside test's table (slot frame set, a, hs, rho)
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
      for (int k = 0; k < kPieces; k++)
      {
         double ds, dt;
         m.Eval(k, 0.0, m.js[k], m.jt[k], ds, dt);
         const double l = std::hypot(ds, dt);
         m.jTs[k] = ds / l;
         m.jTt[k] = dt / l;
      }
      // the world box: the line lies within |s| <= hs + L + rho (the tip), |t| <= r
      for (int i = 0; i < 3; i++)
      {
         m.lo[i] = 1e30;
         m.hi[i] = -1e30;
      }
      for (int corner = 0; corner < 4; corner++)
      {
         const double s = (corner & 1 ? 1.0 : -1.0) * (m.hs + m.uTip), t = (corner & 2 ? 1.0 : -1.0) * m.r;
         const double v[3] = { m.cx + m.lx * s + m.px * t, m.cy + m.ly * s + m.py * t, m.z0 - m.rho };
         for (int i = 0; i < 3; i++)
         {
            m.lo[i] = std::min(m.lo[i], v[i]);
            m.hi[i] = std::max(m.hi[i], v[i]);
         }
      }
      // logged: the length (the outline's, the ramps' by the parabola's arc length, plus rho times the loop's whole turn 2 pi) and the curvature at the tip and on the round end
      const double mJ = std::sqrt(m.c * m.c + 4.0 * m.tauJ * m.tauJ);
      m.total = 4.0 * m.hs + 4.0 * m.a * m.th1 + 2.0 * (m.tauJ * mJ + 0.5 * m.c * m.c * std::asinh(2.0 * m.tauJ / m.c)) + 2.0 * M_PI * m.rho;
      m.kapTip = m.RampCurvature(0.0);
      m.kapRound = 1.0 / m.r;
   }
}

namespace PlayfieldCuts
{
   void Setup(PhysicsEngine *, PinTable *const table)
   {
      g_rims.clear();
      g_rimsTable = table;
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
           << " mm, rim radius " << m.rho * 1e3 << " mm (assumed). Rim centre line (outline + rho): " << m.total * 1e3
           << " mm round, extrema in closed form; curvature toward the hole " << m.kapRound << " /m on the round end, " << m.kapTip << " /m at the tip. Paper, 27 mm ball: sinks " << (R - hc) * 1e3
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
      // The reach (beyond S = R + rho): a rim minimum matters here only if it is deeper than the faces' nearest (it becomes the probe) or presses (d = delta - zero > 0 within S - zero,
      // the patch sum), so with a face left the rim is searched to S + max(the faces' gap, -zero (patch), 0). With no face the deepest minimum is the probe wherever it is: no limit.
      // Every extremum is in closed form, so one search costs the same whatever its reach. Every hump is kept at its distance (one beyond the top gives s = 0).
      const bool face = pr.seg >= 0;
      const double extra = face ? std::max(std::max(-pr.delta, patch ? -zero : 0.0), 0.0) : kNoLimitM;
      hits.clear();
      humps.clear();
      for (const int ci : cuts)
         RimExtrema(ci, bs.R, c, hits, patch ? &humps : nullptr, bs.R + g_rims[ci].rho + extra + kReachPadM);
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
         const double tot = kPieces, top = bs.R + g_rims[ci].rho - zero;
         for (size_t i = 0; i < hits.size(); i++)
            for (size_t j = i + 1; j < hits.size(); j++)
            {
               const RimHit &h1 = hits[i], &h2 = hits[j];
               if (h1.cut != ci || h2.cut != ci || h1.delta - zero <= 0.0 || h2.delta - zero <= 0.0)
                  continue;
               // the largest distance on each arc between them (forward from h1 to h2, and the rest); the hump is the lower of the two
               const double fwd = std::fmod(h2.loop - h1.loop + 2.0 * tot, tot);
               double m1 = -1.0, m2 = -1.0;
               for (const RimHump &q : humps)
               {
                  if (q.cut != ci)
                     continue;
                  double &mm = std::fmod(q.loop - h1.loop + 2.0 * tot, tot) < fwd ? m1 : m2;
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
      // the reach: a minimum farther than this leaves h as it is (ArrivalStep: d > vc h + a h^2, with vc <= |v| and h only falling from here)
      const double h0 = h, beyondM = v.Length() * h0 + std::max(a, 0.0) * h0 * h0 + arriveM;
      for (const int ci : cuts)
      {
         hits.clear();
         RimExtrema(ci, bs.R, c, hits, nullptr, bs.R + g_rims[ci].rho + beyondM + kReachPadM);
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

   std::vector<DrawnSlot> Drawn(const PinTable *const table)
   {
      std::vector<DrawnSlot> out;
      if (table != g_rimsTable)
         return out;
      for (const Rim &m : g_rims)
      {
         DrawnSlot d;
         d.name = m.name;
         d.unit = U;
         d.cx = m.cx;
         d.cy = m.cy;
         d.z0 = m.z0;
         d.lx = m.lx;
         d.ly = m.ly;
         d.px = m.px;
         d.py = m.py;
         d.a = m.a;
         d.hs = m.hs;
         d.L = m.L;
         d.U1 = m.U1;
         d.wJ = m.wJ;
         d.c = m.c;
         d.tauJ = m.tauJ;
         d.th1 = m.th1;
         d.thick = kThickM;
         out.push_back(d);
      }
      return out;
   }

   double DrawnSlot::W(const double u) const
   {
      if (u <= 0.0)
         return a;
      if (u <= U1)
         return std::sqrt(std::max(0.0, a * a - u * u));
      if (u < L)
         return std::sqrt(std::max(0.0, wJ * wJ - 2.0 * U1 * (u - U1)));
      return 0.0;
   }

   double DrawnSlot::C(const double u) const { return u <= 0.0 ? a : u < a ? std::sqrt(a * a - u * u) : 0.0; }

   void DrawnSlot::FloorArc(const double u, double &rho, double &zc) const
   {
      const double D = kRampDepthM * (L - u) / (L - a), w = W(u);
      rho = (w * w + D * D) / (2.0 * D);
      zc = rho - D;
   }
}
