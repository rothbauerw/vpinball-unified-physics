// license:GPLv3+

#include "core/stdafx.h"
#include "physics/playfieldcut.h"
#include "physics/contactsurfacehost.h"
#include "physics/physicsini.h"
#include "physics/contactsystem.h"
#include "physics/contactprof.h"
#include <intrin.h>
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
   constexpr double kGrooveRadiusM = 0.625 * 25.4e-3;   // the plunger lane groove's floor radius (5/8 in, the cutter's; fixed, not a setting: paper runs/2026-11-10-plunger-groove/paper.md section 1)
   constexpr double kGrooveDepthM = 0.188 * 25.4e-3;    // its depth at the plunger end (0.188 in), falling linearly to 0 at the far tip (fixed, not a setting: paper section 1)
   constexpr double kGrooveLenMinMm = 12.0 * 25.4;      // length_mm smallest and largest for a groove: the 12 - 20 in the paper checked (section 3: brute force, cases)
   constexpr double kGrooveLenMaxMm = 20.0 * 25.4;

   void Notice(const std::string &msg)   // rule 8: on screen and in the log, once per table start (Setup runs once per start)
   {
      PLOGW << msg;
      if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
         g_pplayer->m_liveUI->PushNotification(msg, 8000);
   }

   struct CutKeys
   {
      std::string name, shape;
      double x = 0.0, y = 0.0, angle = 0.0, lengthMm = -1.0, widthMm = kWidthMaxMm, cornerMm = 0.0;
      bool xSet = false, ySet = false, lengthSet = false, widthSet = false, cornerSet = false;
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
         else if (key == "width_mm") { cur->widthMm = v; cur->widthSet = true; }
         else if (key == "corner_mm") { cur->cornerMm = v; cur->cornerSet = true; }
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

   // ---- The plunger lane groove (paper runs/2026-11-10-plunger-groove/paper.md sections 1 - 2) ----
   // Frame: X along the groove from the plunger end (X = 0) to the far tip, Y across, Z up from the playfield's top. Three exact pieces, joined tangent:
   //   floor  a cylinder of radius Rg about a straight axis through (0, 0, Pz), Pz = Rg / cos a - D0, direction (cos a, 0, sin a), tan a = D0 / Lg: the vertical depth is D0 (1 - X / Lg)
   //   rim    the rolling-ball fillet between the top and the floor: exactly a tube of radius rho about its spine, the plane Z = -rho cut with the cylinder of radius B = Rg + rho,
   //          a level ellipse X = (B cos phi - u0) / sin a, Y = B sin phi, u0 = (Pz + rho) cos a; its exposed arc runs from +Z (tangent to the top) to -r (r the unit direction from
   //          the axis to the spine point; tangent to the floor)
   //   top    the playfield's own faces outside the spine's outline (FacesLeft, as for the slot)
   // The plunger end is CLOSED (lead 11-10, design runs/2026-11-10-cut-junction section 1: the wood is the slab minus the cut, so the removed volume ends in the vertical plane X = 0): an
   // end face (that plane inside the floor's section), its rim with the top (a straight tube), its edge with the floor (a tube about an ellipse in the plane X = -rho) and the two corner
   // spheres where the three spines meet at (-rho, +-yh, -rho) (GrooveExtrema's end). The groove rim's arc runs on to X = -rho to meet them. A cut covering the end removes them (Gone).
   struct GrooveShape
   {
      double Lg = 0, a = 0, s = 0, c = 1, Rg = kGrooveRadiusM, D0 = kGrooveDepthM, rho = kEdgeRhoM;   // length, tilt (sin, cos), floor radius, depth at X = 0, rim radius
      double Pz = 0, B = 0, K = 0, u0 = 0, phiE = 0, tE = 0, xv = 0, yh = 0;   // axis height at X = 0, spine radius, K = B cos^2 a, u0, the arc's end (X = 0) and tan(phiE / 2), the vertex X, half-width at X = 0
      double floorR1 = 0;   // the floor's Hertz radius across (concave): 1 / (1 / R - 1 / Rg) for the ball it was built for (logged)

      // the height Z (groove frame) where the vertical line through (X, Y) meets the cylinder of radius B about the floor's axis, from below (|Y| < B): (c wz - X s)^2 = B^2 - Y^2,
      // wz = Z - Pz, the lower root. Where a through cut's wall meets the floor, the edge's spine (rho from both) lies on it over the wall's centre line (Join kind 1). dZ/dX, dZ/dY.
      double ZonB(const double X, const double Y, double &zx, double &zy) const
      {
         const double q = std::sqrt(std::max(B * B - Y * Y, 1e-30));
         zx = s / c;
         zy = Y / (q * c);
         return Pz + (X * s - q) / c;
      }
      // the outward radial unit vector from the floor's axis through p (groove frame), and p's distance from the axis
      double Radial(const double X, const double Y, const double Z, double &rx, double &ry, double &rz) const
      {
         const double wz = Z - Pz, al = X * c + wz * s;
         rx = X - al * c;
         ry = Y;
         rz = wz - al * s;
         const double l = std::sqrt(rx * rx + ry * ry + rz * rz);
         if (l > 0.0)
         {
            rx /= l;
            ry /= l;
            rz /= l;
         }
         return l;
      }

      void Build(const double length)
      {
         Lg = length;
         a = std::atan(D0 / Lg);
         s = std::sin(a);
         c = std::cos(a);
         Pz = Rg / c - D0;
         B = Rg + rho;
         K = B * c * c;
         u0 = (Pz + rho) * c;
         phiE = std::acos((u0 - rho * s) / B);   // the spine at X = -rho (the closed end's corner)
         tE = std::tan(0.5 * phiE);
         xv = (B - u0) / s;
         yh = B * std::sin(phiE);
      }
      // the outline the top leaves out: inside the spine's ellipse, on the groove's side of the open plunger end
      bool InsideOutline(const double X, const double Y) const
      {
         if (X < -rho)   // (the closed end: the top ends at the end rim's tangent line)
            return false;
         const double u = u0 + X * s;
         return std::abs(u) < B && Y * Y < B * B - u * u;
      }
      // h(t), t = tan(phi / 2): the sign of d/dphi of the squared distance from the ball's centre to the spine (paper section 1); P = u0 + X s, Q = Y s^2
      double H(const double t, const double P, const double Q) const { return ((Q * t + 2.0 * (K + P)) * t * t - 2.0 * (K - P)) * t - Q; }
      double Hp(const double t, const double P, const double Q) const { return (4.0 * Q * t + 6.0 * (K + P)) * t * t - 2.0 * (K - P); }
   };

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
      bool groove = false;   // a plunger lane groove (gr; the frame above with (cx, cy) its plunger end on the axis, along toward the tip), not a slot
      GrooveShape gr;
      bool rect = false;     // a rectangle (RectExtrema): hs, a its half length and half width, k its corner radius; the frame above with (cx, cy) its centre, along its length
      double k = 0.0;
      // cuts that meet (design runs/2026-11-10-cut-junction, found once at load, Meets): meets, the other cuts on this host whose boxes come within rho of this one's (each of this cut's
      // pieces counts only where it bounds wood, Gone); joins, the rims and corners where this cut's surfaces meet another's that this cut's probe carries (g_joins). Both empty: alone,
      // every piece as before
      std::vector<int> meets, joins;

      // p strictly within rho of this cut's removed volume: a through cut's centre line's inside (the outline offset by rho), through the playfield's thickness as Holds; a groove's
      // cylinder of radius B = Rg + rho about the floor's axis, on the groove's side of the open end. Another cut's piece there bounds no wood (design section 1)
      bool Near(const V3 &p) const
      {
         if (p.z > z0 + rho || p.z < z0 - kThickM - rho)
            return false;
         const double dx = p.x - cx, dy = p.y - cy, s = dx * lx + dy * ly, t = dx * px + dy * py;
         if (!groove)
            return rect ? InsideRect(s, t) : Inside(s, t);
         const double wz = p.z - z0 - gr.Pz, al = s * gr.c + wz * gr.s;
         return s > -gr.rho && s * s + t * t + wz * wz - al * al < gr.B * gr.B;   // (within rho of the closed end's plane X = 0)
      }

      // inside a rectangle's rim centre line (the outline offset outward by rho: half sides hs + rho, a + rho, corner radius k + rho about the corner centres (+-(hs - k), +-(a - k)))
      bool InsideRect(const double s, const double t) const
      {
         const double S = std::abs(s), T = std::abs(t), ds = std::max(S - (hs - k), 0.0), dt = std::max(T - (a - k), 0.0);
         return S < hs + rho && T < a + rho && ds * ds + dt * dt < (k + rho) * (k + rho);
      }

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
      // a point of the playfield's faces that is not surface: inside the centre line (a groove: its spine), at the top or below it within the playfield's thickness (the through-hole or the
      // groove; a thin playfield mesh's underside)
      bool Holds(const V3 &p) const
      {
         if (p.z > z0 + 1e-6 || p.z < z0 - kThickM - 1e-6)
            return false;
         const double dx = p.x - cx, dy = p.y - cy;
         const double s = dx * lx + dy * ly, t = dx * px + dy * py;
         return groove ? gr.InsideOutline(s, t) : rect ? InsideRect(s, t) : Inside(s, t);
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

   // A cut's contact: one true local minimum of the distance from the ball's centre to a rim's centre line, the rim's surface facing the ball (its upper inner part), or to a groove's floor
   struct RimHit
   {
      int cut;
      double delta, dist, reff, loop;   // depth R + rho - dist, the law's effective radius, the place on the loop (slot: piece + fraction, 0 .. kPieces; groove rim: the spine's phi)
      double rho;                       // the piece's tube radius: the rim's rho; 0 for a groove's floor (dist is then to the surface)
      bool floor;                       // a groove's floor (no hump lies between it and another contact)
      V3 point, normal;
      int edge = -1;                    // a floor-wall edge's (Join kind 1): the pair's key (groove << 16 | through cut), loop the through cut's piece + fraction; else -1
   };
   // A local maximum of the same distance (the hump between two minima along the line): the place on the loop and the distance
   struct RimHump
   {
      int cut;
      double loop, dist;
      int edge = -1;   // as RimHit::edge
   };

   std::vector<Rim> g_rims;   // this table's cuts, rebuilt at every table load (Setup)
   uint64_t g_rimsVer = 0;    // shared_finder P4b: changes at every Setup (a ball's cut query made before it is stale)
   std::vector<const ContactSurface *> g_anyHosts;   // P4b (key on): the hosts holding cuts, found once per Setup (Any)
   uint64_t g_anyVer = 0;
   PlayfieldCuts::CutStats g_cutStats;   // P4b counters (diagnostics only)
   bool CutCheckOn()   // P4b check mode (VPX_FINDER_CHECK=1): every NearCuts answer against today's scan (diagnostics only)
   {
      static const bool on = [] { const char *e = std::getenv("VPX_FINDER_CHECK"); return e != nullptr && e[0] == '1'; }();
      return on;
   }
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

   // Cuts that meet (design runs/2026-11-10-cut-junction sections 1 - 2): the wood is the slab minus the union of the cuts, and the ball meets only its boundary. A piece of cut m (a rim's
   // spine point, or a surface's foot point moved rho into the wood) bounds wood only if that point is not within rho of another cut's removed volume (Rim::Near; x1, x2: cuts left out,
   // the two a junction lies on). Every piece's minimum is tested so: a minimum inside another cut is not a contact, and the piece's distance then rises from its clipped end, where
   // the junction's rim or corner sphere (Join) meets it tangent and takes over, so the switch is continuous. A cut that meets none: always true (every piece as before).
   bool Gone(const Rim &m, const V3 &p, const int x1 = -1, const int x2 = -1)
   {
      for (const int j : m.meets)
         if (j != x1 && j != x2 && g_rims[j].Near(p))
            return true;
      return false;
   }

   // A rim where two cuts' surfaces meet, or a corner where three do (design sections 1 - 2), found once at load (Meets), carried by its owner's probe (JoinExtrema); each is the 0.25 mm
   // rounding of a convex edge of the wood, one Hertz law:
   //   kind 0  where two through cuts' outlines cross: a vertical tube of radius rho about the line through V = the crossing of their centre lines (rho from both walls), from the top
   //           rim's spine height to the underside rim's; n[0], n[1] the two walls' normals into their holes (horizontal). The wood there is outside both: a wedge under 180 deg, convex.
   //   kind 1  where a groove's floor meets a through cut's wall (b): a tube of radius rho about the spine rho from both, the wall's centre line P(p) (piece, p0 .. p1 of cut b, where
   //           it is inside the groove's spine outline) lifted onto the floor's cylinder of radius B (GrooveShape::ZonB): on a straight wall an ellipse, on a curved wall (a round end,
   //           a corner) a space curve. A slot inside the groove: its whole loop.
   //   kind 2  a sphere of radius rho about S where three surfaces meet (the ends of kind 0 and kind 1 rims, where they meet the top or underside rims): n[3] the surfaces' normals out of
   //           the wood (the rounded corner's region is their cone)
   struct Join
   {
      int kind = 0, owner = 0, a = 0, b = 0, piece = 0;   // kind 0: the two through cuts; kind 1: a the groove, b the through cut
      double p0 = 0.0, p1 = 0.0;                           // kind 1: the part of the piece
      double s2 = 0.0, a2 = 0.0;                           // kind 1: bounds on |E'|^2 and |E''| over p0 .. p1 (EdgeBounds), for the root search's exclusion test
      double rad = 0.0;                                    // kind 1: every spine point lies within rad of S (groove frame; EdgeBounds), the probe's first reach test
      V3 S, n[3];
   };
   std::vector<Join> g_joins;   // rebuilt at every table load (Meets)

   // a through cut's rim centre line (the outline offset outward by rho) in the world, by piece in loop order, the hole on the right (the slot's Rim::Eval; a rectangle's: side j at
   // 2 j, p from 0 to 2 (its straight half length), then corner j, p the angle from side j's normal to side j + 1's): the point (x, y) and its derivative in p
   double PieceLen(const Rim &m, const int k)
   {
      if (m.rect)
         return k & 1 ? 0.5 * M_PI : 2.0 * ((k / 2) & 1 ? m.a - m.k : m.hs - m.k);
      return k == 0 || k == 4 ? 2.0 * m.hs : k == 2 || k == 6 ? 2.0 * m.tauJ : m.th1;
   }
   void LineAt(const Rim &m, const int k, const double p, double &x, double &y, double &dx, double &dy)
   {
      double s, t, ds, dt;
      if (m.rect)
      {
         static const double Ns[4] = { 0, 1, 0, -1 }, Nt[4] = { 1, 0, -1, 0 }, Ts[4] = { 1, 0, -1, 0 }, Tt[4] = { 0, -1, 0, 1 }, Os[4] = { 1, 1, -1, -1 }, Ot[4] = { 1, -1, -1, 1 };
         const int j = k / 2, j1 = (j + 1) % 4;
         const double off = j & 1 ? m.hs : m.a, ext = j & 1 ? m.a - m.k : m.hs - m.k;
         if (!(k & 1))
         {
            const double xx = p - ext;
            s = xx * Ts[j] + (off + m.rho) * Ns[j];
            t = xx * Tt[j] + (off + m.rho) * Nt[j];
            ds = Ts[j];
            dt = Tt[j];
         }
         else
         {
            const double r = m.k + m.rho, cp = std::cos(p), sp = std::sin(p);
            s = Os[j] * (m.hs - m.k) + r * (Ns[j] * cp + Ns[j1] * sp);
            t = Ot[j] * (m.a - m.k) + r * (Nt[j] * cp + Nt[j1] * sp);
            ds = r * (Ns[j1] * cp - Ns[j] * sp);
            dt = r * (Nt[j1] * cp - Nt[j] * sp);
         }
      }
      else
         m.Eval(k, p, s, t, ds, dt);
      x = m.cx + m.lx * s + m.px * t;
      y = m.cy + m.ly * s + m.py * t;
      dx = m.lx * ds + m.px * dt;
      dy = m.ly * ds + m.py * dt;
   }

   // the same centre line's second derivative in p (world), closed form: sides 0; round ends and corners (p the angle) -(P - O), O the arc's centre; a slot's ramp from Rim::Ramp's
   // u'' = -2 - 4 rho c (m^2 - 12 tau^2) / m^5, t'' = -24 rho c^2 tau / m^5 (m = sqrt(c^2 + 4 tau^2)), with the piece's signs (Eval: tau = tauJ - p on 2, -tauJ + p on 6)
   void LineAt2(const Rim &m, const int k, const double p, double &ddx, double &ddy)
   {
      double dds = 0.0, ddt = 0.0;
      if (m.rect)
      {
         if (k & 1)
         {
            static const double Ns[4] = { 0, 1, 0, -1 }, Nt[4] = { 1, 0, -1, 0 };
            const int j = k / 2, j1 = (j + 1) % 4;
            const double r = m.k + m.rho, cp = std::cos(p), sp = std::sin(p);
            dds = -r * (Ns[j] * cp + Ns[j1] * sp);
            ddt = -r * (Nt[j] * cp + Nt[j1] * sp);
         }
      }
      else if (k == 2 || k == 6)
      {
         const double tau = k == 2 ? m.tauJ - p : -m.tauJ + p, m2 = m.c * m.c + 4.0 * tau * tau, m5 = m2 * m2 * std::sqrt(m2);
         const double uu = -2.0 - 4.0 * m.rho * m.c * (m2 - 12.0 * tau * tau) / m5, tt = -24.0 * m.rho * m.c * m.c * tau / m5;
         dds = k == 2 ? uu : -uu;
         ddt = tt;
      }
      else if (k != 0 && k != 4)
      {
         double s, t, ds, dt;
         m.Eval(k, p, s, t, ds, dt);
         dds = -(s - (k < 4 ? m.hs : -m.hs));
         ddt = -t;
      }
      ddx = m.lx * dds + m.px * ddt;
      ddy = m.ly * dds + m.py * ddt;
   }

   // A floor-wall edge's spine (Join kind 1) in the groove g's frame: the through cut w's centre line, piece k at p, lifted onto the floor's cylinder of radius B (GrooveShape::ZonB): the
   // point E, its derivative E' in p, and a, the wall's horizontal normal into the through cut; dd given: E'' too, closed form (Z = Pz + (X s - q) / c, q = sqrt(B^2 - Y^2):
   // Z'' = (s X'' + (Y'^2 + Y Y'') / q + Y^2 Y'^2 / q^3) / c)
   void EdgeSpine(const Rim &g, const Rim &w, const int k, const double p, double e[3], double d[3], double an[2], double *dd = nullptr)
   {
      double x, y, dx, dy, zx, zy;
      LineAt(w, k, p, x, y, dx, dy);
      x -= g.cx;
      y -= g.cy;
      e[0] = x * g.lx + y * g.ly;
      e[1] = x * g.px + y * g.py;
      d[0] = dx * g.lx + dy * g.ly;
      d[1] = dx * g.px + dy * g.py;
      e[2] = g.gr.ZonB(e[0], e[1], zx, zy);
      d[2] = zx * d[0] + zy * d[1];
      const double l = std::hypot(d[0], d[1]);
      an[0] = d[1] / l;
      an[1] = -d[0] / l;
      if (dd != nullptr)
      {
         double ddx, ddy;
         LineAt2(w, k, p, ddx, ddy);
         dd[0] = ddx * g.lx + ddy * g.ly;
         dd[1] = ddx * g.px + ddy * g.py;
         const GrooveShape &G = g.gr;
         const double q2 = std::max(G.B * G.B - e[1] * e[1], 1e-30), q = std::sqrt(q2);
         dd[2] = (G.s * dd[0] + (d[1] * d[1] + e[1] * dd[1]) / q + e[1] * e[1] * d[1] * d[1] / (q2 * q)) / G.c;
      }
   }
   // At load: s2 >= |E'|^2 and a2 >= |E''| over the edge's part p0 .. p1, for the root search's exclusion test (JoinExtrema). A stated bound: the largest of 2048 even samples (E'' by
   // the central difference of the closed-form E' over 1e-7 m), times 1.5 (the pieces are smooth: between two samples, 1/2048 of the part, they change by far less than half).
   void EdgeBounds(const Rim &g, const Rim &w, Join &J)
   {
      constexpr int kS = 2048;
      double s2 = 0.0, a2 = 0.0, ec[3], rad = 0.0;
      {
         double d[3], an[2];
         EdgeSpine(g, w, J.piece, 0.5 * (J.p0 + J.p1), ec, d, an);
      }
      for (int i = 0; i <= kS; i++)
      {
         const double p = J.p0 + (J.p1 - J.p0) * i / kS;
         double e[3], d[3], an[2], e1[3], d1[3], e2[3], d2[3];
         EdgeSpine(g, w, J.piece, p, e, d, an);
         rad = std::max(rad, std::sqrt((e[0] - ec[0]) * (e[0] - ec[0]) + (e[1] - ec[1]) * (e[1] - ec[1]) + (e[2] - ec[2]) * (e[2] - ec[2])));
         const double sp2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2], h = 1e-7 / std::sqrt(std::max(sp2, 1e-300));
         EdgeSpine(g, w, J.piece, p + h, e1, d1, an);
         EdgeSpine(g, w, J.piece, p - h, e2, d2, an);
         const double ax = (d1[0] - d2[0]) / (2.0 * h), ay = (d1[1] - d2[1]) / (2.0 * h), az = (d1[2] - d2[2]) / (2.0 * h);
         s2 = std::max(s2, sp2);
         a2 = std::max(a2, std::sqrt(ax * ax + ay * ay + az * az));
      }
      J.s2 = 1.5 * s2;
      J.a2 = 1.5 * a2;
      // (between two samples a point is within |E'| (p1 - p0) / kS / 2 of one)
      J.S = V3(ec[0], ec[1], ec[2]);
      J.rad = rad + std::sqrt(J.s2) * (J.p1 - J.p0) / kS * 0.5 + 1e-9;
   }

   // n in the cone of m[0..2] (n = sum of nonnegative multiples of them; Cramer's rule): the corner sphere's own region
   bool InCone(const V3 &n, const V3 m[3])
   {
      const auto det = [](const V3 &a, const V3 &b, const V3 &c) { return a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) + a.z * (b.x * c.y - b.y * c.x); };
      const double D = det(m[0], m[1], m[2]);
      if (std::abs(D) < 1e-12)
         return false;
      const double s = D > 0.0 ? 1.0 : -1.0;
      return s * det(n, m[1], m[2]) >= 0.0 && s * det(m[0], n, m[2]) >= 0.0 && s * det(m[0], m[1], n) >= 0.0;
   }

   // The minima of the distance from the ball's centre to cut ci's joins (Join), as RimHits of cut ci, weight 1 in the patch sum (floor flag: each meets its neighbours tangent, so it is
   // never a second minimum beside one of them at the same place). The law: a tube's across 1/R1 = 1/R + 1/rho, along it the spine's curvature by Meusnier (TubeAlongCurvature), concave
   // bends held within 20 ball radii; a sphere's 1/Reff = 1/R + 1/rho.
   //   kind 0  the vertical line's foot (the ball's height between the spines): a contact when the horizontal direction to the ball lies between the two walls' normals
   //   kind 1  the spine's minima: f(p) = (E - C).E' changes sign from < 0 to >= 0, and its humps (when humps is given) from >= 0 to < 0. A stated solve that finds every one more than
   //           1 um apart along the spine: a bracket lo .. hi is dropped when no point of it is within reach, or when f keeps one sign and cannot reach zero, |f(lo)| + |f(hi)| > L (hi - lo)
   //           with L >= |f'| = |E'^2 + (E - C).E''| on it (EdgeBounds: s2 + (the nearer end's distance + the bracket's length) a2); otherwise it is halved, down to 1 um of spine, and a
   //           sign change there is bisected to 1e-13 m. (Measured 11-10, runs/2026-11-10-cut-junction results_fix.txt: the earlier 17 even samples missed one of the two minima a slot's
   //           tip splits into, a minimum and a hump in one bracket, so the contact count jumped 1 -> 2.) A minimum at a piece's end is the piece's that ends there (f continuous across: the
   //           line is tangent-continuous). A contact when the direction spine -> centre lies between the wall's normal a (into the through cut) and the floor's b = -r (both normal to
   //           the spine). Several minima on one edge are weighted by the hump between them (Probe, as the slot's rim), so a minimum splitting in two is continuous. E'' for the curvature:
   //           central difference of E' (closed form) over 1e-7 m along (relative error ~1e-6, in the Hertz radius only).
   //   kind 2  the sphere's centre: a contact when the direction lies in the cone of the three normals
   void JoinExtrema(const int ci, const double R, const V3 &c, std::vector<RimHit> &out, std::vector<RimHump> *humps, const double reachGap)
   {
      const double rho = g_rims[ci].rho, R1 = 1.0 / (1.0 / R + 1.0 / rho);
      const auto hit = [&](const V3 &S, const V3 &dv, const double dist, const double reff) {
         RimHit h;
         h.cut = ci;
         h.dist = dist;
         h.delta = R + rho - dist;
         h.rho = rho;
         h.floor = true;
         h.loop = 0.0;
         h.normal = dv / dist;
         h.point = S;
         h.reff = reff;
         out.push_back(h);
      };
      for (const int ji : g_rims[ci].joins)
      {
         const Join &J = g_joins[ji];
         if (J.kind == 2)
         {
            const V3 dv = c - J.S;
            const double d = dv.Length();
            if (d > 1e-12 && d - (R + rho) <= reachGap && InCone(dv, J.n))
               hit(J.S, dv, d, R1);
            continue;
         }
         if (J.kind == 0)
         {
            const Rim &A = g_rims[J.a];
            if (c.z >= A.z0 - rho || c.z <= A.z0 - kThickM + rho)
               continue;
            const V3 S(J.S.x, J.S.y, c.z), dv = c - S;
            const double d = dv.Length(), D = J.n[0].x * J.n[1].y - J.n[0].y * J.n[1].x;
            if (d <= 1e-12 || d - (R + rho) > reachGap || std::abs(D) < 1e-12)
               continue;
            const double al = (dv.x * J.n[1].y - dv.y * J.n[1].x) / D, be = (J.n[0].x * dv.y - J.n[0].y * dv.x) / D;
            if (al >= 0.0 && be >= 0.0 && !Gone(g_rims[ci], S, J.a, J.b))
               hit(S, dv, d, std::sqrt(R1 * R));
            continue;
         }
         // kind 1, in the groove's frame
         const Rim &g = g_rims[J.a], &w = g_rims[J.b];
         const GrooveShape &G = g.gr;
         const double cx = c.x - g.cx, cy = c.y - g.cy, X = cx * g.lx + cy * g.ly, Y = cx * g.px + cy * g.py, Z = c.z - g.z0;
         // every point of the part within rad of S (EdgeBounds): out of reach as a whole
         if (std::sqrt((X - J.S.x) * (X - J.S.x) + (Y - J.S.y) * (Y - J.S.y) + (Z - J.S.z) * (Z - J.S.z)) - J.rad > R + rho + reachGap)
            continue;
         // the spine E(p) and E'(p), and a, the wall's normal into the through cut
         const auto spine = [&](const double p, double e[3], double d[3], double an[2]) { EdgeSpine(g, w, J.piece, p, e, d, an); };
         // f(p) and the distance |E - C|
         const auto ev = [&](const double p, double &fv, double &dv) {
            double e[3], d[3], an[2];
            spine(p, e, d, an);
            const double ex = e[0] - X, ey = e[1] - Y, ez = e[2] - Z;
            fv = ex * d[0] + ey * d[1] + ez * d[2];
            dv = std::sqrt(ex * ex + ey * ey + ez * ez);
         };
         const int key = J.a << 16 | J.b;
         const double reach = R + rho + reachGap, sq = std::sqrt(J.s2), len = PieceLen(w, J.piece);
         struct Br
         {
            double lo, hi, flo, fhi, dlo, dhi;
         };
         constexpr int kStart = 8;
         Br stack[kStart + 64];
         int ns = 0;
         {
            double f0, d0;
            ev(J.p1, f0, d0);
            for (int i = kStart - 1; i >= 0; i--)   // (pushed from the far end: popped in order along the line)
            {
               Br b;
               b.hi = J.p0 + (J.p1 - J.p0) * (i + 1) / kStart;
               b.fhi = f0;
               b.dhi = d0;
               b.lo = J.p0 + (J.p1 - J.p0) * i / kStart;
               ev(b.lo, f0, d0);
               b.flo = f0;
               b.dlo = d0;
               stack[ns++] = b;
            }
         }
         while (ns > 0)
         {
            Br b = stack[--ns];
            const double wd = b.hi - b.lo, al = wd * sq, dm = std::min(b.dlo, b.dhi);
            if (dm - al > reach)   // no point of it within reach
               continue;
            const bool nlo = b.flo < 0.0, nhi = b.fhi < 0.0;
            if (nlo == nhi && std::abs(b.flo) + std::abs(b.fhi) > (J.s2 + (dm + al) * J.a2) * wd)   // one sign, and f cannot reach zero in it
               continue;
            if (al > 1e-6 && ns + 2 <= (int)(sizeof stack / sizeof stack[0]))   // halve
            {
               Br l = b, r = b;
               l.hi = r.lo = 0.5 * (b.lo + b.hi);
               ev(l.hi, l.fhi, l.dhi);
               r.flo = l.fhi;
               r.dlo = l.dhi;
               stack[ns++] = r;
               stack[ns++] = l;
               continue;
            }
            if (nlo == nhi)   // (two roots within 1 um of spine: a minimum and a hump being born, not resolved)
               continue;
            // the sign change's root: Newton on f, f' = |E'|^2 + (E - C).E'' (E'' closed form), held in the bracket (a step out of it bisects), to 1e-13 m of spine
            double lo = b.lo, hi = b.hi, p = 0.5 * (lo + hi);
            double e[3], d[3], an[2], dd[3];
            for (int it = 0; it < 60; it++)
            {
               EdgeSpine(g, w, J.piece, p, e, d, an, dd);
               const double ex = e[0] - X, ey = e[1] - Y, ez = e[2] - Z, fm = ex * d[0] + ey * d[1] + ez * d[2];
               const double fp = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + ex * dd[0] + ey * dd[1] + ez * dd[2];
               ((fm < 0.0) == nlo ? lo : hi) = p;
               double next = fp != 0.0 ? p - fm / fp : 0.5 * (lo + hi);
               if (!(next > lo && next < hi))
                  next = 0.5 * (lo + hi);
               const bool done = std::abs(next - p) * sq < 1e-13 || (hi - lo) * sq < 1e-13;
               p = next;
               if (done)
                  break;
            }
            EdgeSpine(g, w, J.piece, p, e, d, an, dd);
            const double ex = X - e[0], ey = Y - e[1], ez = Z - e[2], dist = std::sqrt(ex * ex + ey * ey + ez * ez);
            if (!nlo)   // a hump (f from >= 0 to < 0)
            {
               if (humps != nullptr)
                  humps->push_back({ ci, J.piece + p / len, dist, key });
               continue;
            }
            if (dist <= 1e-12 || dist - (R + rho) > reachGap)
               continue;
            double rx, ry, rz;
            G.Radial(e[0], e[1], e[2], rx, ry, rz);
            const double nx = ex / dist, ny = ey / dist, nz = ez / dist, ab = -(an[0] * rx + an[1] * ry), na = nx * an[0] + ny * an[1], nb = -(nx * rx + ny * ry + nz * rz);
            if (na - ab * nb < 0.0 || nb - ab * na < 0.0)
               continue;
            const V3 S(g.cx + g.lx * e[0] + g.px * e[1], g.cy + g.ly * e[0] + g.py * e[1], g.z0 + e[2]);
            if (Gone(g, S, J.a, J.b))
               continue;
            // the curvature vector (E'' less its part along E') / |E'|^2, E'' closed form (EdgeSpine)
            const double sp2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            double kv[3], ddt = 0.0;
            for (int i3 = 0; i3 < 3; i3++)
               ddt += dd[i3] * d[i3];
            for (int i3 = 0; i3 < 3; i3++)
               kv[i3] = (dd[i3] - ddt / sp2 * d[i3]) / sp2;
            const double k = std::sqrt(kv[0] * kv[0] + kv[1] * kv[1] + kv[2] * kv[2]);
            double kLong = 0.0;
            if (k > 1e-9)
               kLong = RubberContact::TubeAlongCurvature(k, rho, -(nx * kv[0] + ny * kv[1] + nz * kv[2]) / k);
            const double R2 = 1.0 / std::max(1.0 / R + kLong, 0.05 / R);
            hit(S, c - S, dist, std::sqrt(R1 * R2));
            out.back().edge = key;
            out.back().loop = J.piece + p / len;
         }
      }
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

   // Every rim contact (true local minimum of the distance, facing the ball, its gap dist - (R + rho) within reachGap) of cut ci for a ball (radius R) at c, and every local maximum of the distance
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
   void GrooveExtrema(const int ci, const double R, const V3 &c, std::vector<RimHit> &out, std::vector<RimHump> *humps, const double reachGap);
   void RectExtrema(const int ci, const double R, const V3 &c, std::vector<RimHit> &out, std::vector<RimHump> *humps, const double reachGap);
   void RimExtrema(const int ci, const double R, const V3 &c, std::vector<RimHit> &out, std::vector<RimHump> *humps, const double reachGap)
   {
      const Rim &m = g_rims[ci];
      if (!m.joins.empty())
         JoinExtrema(ci, R, c, out, humps, reachGap);
      if (m.groove)
      {
         GrooveExtrema(ci, R, c, out, humps, reachGap);
         return;
      }
      if (m.rect)
      {
         RectExtrema(ci, R, c, out, humps, reachGap);
         return;
      }
      const double reachMin = R + m.rho + reachGap;
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
         if (!m.meets.empty() && Gone(m, A))   // inside another cut: no wood there (Gone)
            return;
         RimHit h;
         h.cut = ci;
         h.point = A;
         h.dist = dist;
         h.normal = dv / dist;
         h.delta = R + m.rho - dist;
         h.loop = loop;
         h.rho = m.rho;
         h.floor = false;
         // the tube law (rubbercontact.h ProbeTube): across the rim 1/R1 = 1/R + 1/rho; along it the centre line's curvature (TubeAlongCurvature, Meusnier: the 3D cosine between the normal
         // and the line's in-plane normal), concave bends held within 20 ball radii
         double kLong = 0.0;
         if (std::abs(kin) > 1e-9)
         {
            const double sg = kin > 0.0 ? 1.0 : -1.0, nis = Tt * sg, nit = -Ts * sg;   // toward the centre of curvature (slot frame)
            const double dnx = m.lx * nis + m.px * nit, dny = m.ly * nis + m.py * nit;
            kLong = RubberContact::TubeAlongCurvature(std::abs(kin), m.rho, -(h.normal.x * dnx + h.normal.y * dny));
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

   // A groove's contacts (paper section 1, "Extrema, per piece"): every true local minimum of the distance from the ball's centre to its floor and to its rim, each piece in closed form, and
   // (humps given) the rim's humps; gap = dist - (R + rho) within reachGap. The top is the playfield's faces (FacesLeft). Same rules as the slot: every true local minimum is a contact.
   //   floor  the foot on the ray from the axis through the ball's centre: dist = Rg - r (r the centre's distance from the axis); a contact when the foot lies on the floor, within
   //          +-theta_r of straight down in its section (cos theta_r = e / B, e = (the axis' height there + rho) / cos a) and on the groove's side of the plunger end
   //   rim    the spine is level, so its 3D extrema are the plane's: the real roots of h(t) on the arc |phi| <= phiE (GrooveShape::H, no squaring: every real root is a true extremum).
   //          h's critical points (exact for Q = 0, 3 Newton steps) split the arc into monotone brackets; a sign change holds one root, seeded by the depressed cubic with Q t^4 dropped
   //          (CubicRoots) and polished by 2 Newton steps held in its bracket (paper: 1.2e-15 mm from the exact roots). A minimum is a contact when the direction spine -> centre lies in the
   //          exposed arc (from +Z to -r, both in the spine's normal plane).
   void GrooveExtrema(const int ci, const double R, const V3 &c, std::vector<RimHit> &out, std::vector<RimHump> *humps, const double reachGap)
   {
      const Rim &m = g_rims[ci];
      const GrooveShape &g = m.gr;
      const double dx = c.x - m.cx, dy = c.y - m.cy;
      const double X = dx * m.lx + dy * m.ly, Y = dx * m.px + dy * m.py, Z = c.z - m.z0;   // the ball's centre in the groove's frame
      const auto world = [&](const double x, const double y, const double z) { return V3(m.cx + m.lx * x + m.px * y, m.cy + m.ly * x + m.py * y, m.z0 + z); };
      const auto worldDir = [&](const double x, const double y, const double z) { return V3(m.lx * x + m.px * y, m.ly * x + m.py * y, z); };
      // the floor: w = C - (axis point at X = 0), sax its part along the axis, v the rest (perpendicular to the axis), up = (-sin a, 0, cos a)
      {
         const double wz = Z - g.Pz, sax = X * g.c + wz * g.s;
         const double vx = X - sax * g.c, vy = Y, vz = wz - sax * g.s, rc = std::sqrt(vx * vx + vy * vy + vz * vz);
         const double down = g.s * vx - g.c * vz;   // -v.up: rc cos(theta), theta the foot's angle from straight down in its section
         const double e = (g.Pz + sax * g.s + g.rho) / g.c;
         // on the floor: past the plunger end (open: sax >= 0), before the tip (e < B), within the edge: cos theta >= cos theta_r = e / B (theta_r < 90 deg, e > 0); on wood: the foot
         // moved rho into the wood (the point B from the axis on its ray) not inside another cut (Gone; where it meets a through cut's wall, the floor ends on the edge tube's tangent line)
         if (rc > 1e-12 && sax * g.c + g.B * vx / rc >= -g.rho && e < g.B && down * g.B >= e * rc)   // (closed end: the foot moved rho into the wood up to the end edge's spine)
         {
            const double dist = g.Rg - rc;
            if (dist - R <= reachGap && (m.meets.empty() || !Gone(m, world(sax * g.c + g.B * vx / rc, g.B * vy / rc, g.Pz + sax * g.s + g.B * vz / rc))))
            {
               RimHit h;
               h.cut = ci;
               h.dist = dist;
               h.delta = R - dist;
               h.rho = 0.0;
               h.floor = true;
               h.loop = 0.0;
               h.normal = worldDir(-vx / rc, -vy / rc, -vz / rc);   // surface -> centre: toward the axis
               h.point = world(sax * g.c + g.Rg * vx / rc, g.Rg * vy / rc, g.Pz + sax * g.s + g.Rg * vz / rc);
               // the law for a concave surface (rubbercontact.h EffectiveRadius(R, -Rg)): across 1/R1 = 1/R - 1/Rg, held within 20 ball radii; along the straight generator R2 = R
               const double R1 = 1.0 / std::max(1.0 / R - 1.0 / g.Rg, 0.05 / R);
               h.reff = std::sqrt(R1 * R);
               out.push_back(h);
            }
         }
      }
      // the rim (where another cut meets it, it ends there: Gone)
      const double tEnd = g.tE;
      const double P = g.u0 + X * g.s, Q = Y * g.s * g.s, KP = g.K + P;
      if (KP <= 0.0)   // (the centre a metre behind the plunger end, past the ellipse's far end: no extremum on the arc)
         return;
      double br[4];   // the monotone brackets' ends, ascending in t
      int nb = 0;
      br[nb++] = -tEnd;
      if (g.K > P)
         for (const double sg : { -1.0, 1.0 })
         {
            double tc = sg * std::sqrt((g.K - P) / (3.0 * KP));   // h's critical points, exact for Q = 0; 3 Newton steps on h'
            for (int it = 0; it < 3; it++)
            {
               const double d2 = 12.0 * Q * tc * tc + 12.0 * KP * tc;
               if (d2 == 0.0)
                  break;
               tc -= g.Hp(tc, P, Q) / d2;
            }
            if (tc > -tEnd && tc < tEnd && tc > br[nb - 1])
               br[nb++] = tc;
         }
      br[nb++] = tEnd;
      double seed[3];
      int ns = 0;
      for (int i = 0; i + 1 < nb; i++)
      {
         const double lo = br[i], hi = br[i + 1], hl = g.H(lo, P, Q), hh = g.H(hi, P, Q);
         if ((hl > 0.0) == (hh > 0.0))   // (a root exactly on a bracket's end counts as h <= 0, once)
            continue;
         if (ns == 0)
            ns = CubicRoots(-(g.K - P) / KP, -Q / (2.0 * KP), seed);
         double t = 0.5 * (lo + hi);
         for (int j = 0; j < ns; j++)
            if (seed[j] >= lo && seed[j] <= hi)
            {
               t = seed[j];
               break;
            }
         for (int it = 0; it < 2; it++)
         {
            const double d = g.Hp(t, P, Q);
            if (d != 0.0)
               t = std::clamp(t - g.H(t, P, Q) / d, lo, hi);
         }
         const double t2 = t * t, cp = (1.0 - t2) / (1.0 + t2), sp = 2.0 * t / (1.0 + t2), phi = 2.0 * std::atan(t);
         const double Sx = (g.B * cp - g.u0) / g.s, Sy = g.B * sp, Sz = -g.rho;   // the spine point
         const double ex = X - Sx, ey = Y - Sy, ez = Z - Sz, d3 = std::sqrt(ex * ex + ey * ey + ez * ez);
         if (hh <= 0.0)   // h from > 0 to <= 0: the distance rises then falls, a hump
         {
            if (humps != nullptr)
               humps->push_back({ ci, phi, d3 });
            continue;
         }
         if (d3 <= 1e-12 || d3 - (R + g.rho) > reachGap)
            continue;
         if (!m.meets.empty() && Gone(m, world(Sx, Sy, Sz)))   // inside another cut: no wood there (Gone)
            continue;
         // the exposed arc: n between +Z and m = -r = (-s cos phi, -sin phi, c cos phi) in the spine's normal plane, nin = m's horizontal part made unit (toward the groove, the spine's
         // centre of curvature): n.nin >= 0 and n no farther from +Z than m
         const double mh = std::sqrt(g.s * g.s * cp * cp + sp * sp), ninx = -g.s * cp / mh, niny = -sp / mh;   // (mh: m's horizontal length, the sine of the arc's end angle)
         const double nx = ex / d3, ny = ey / d3, nz = ez / d3, nn = nx * ninx + ny * niny;
         if (nn < 0.0 || nn * (g.c * cp) - nz * mh > 0.0)
            continue;
         RimHit h;
         h.cut = ci;
         h.dist = d3;
         h.delta = R + g.rho - d3;
         h.rho = g.rho;
         h.floor = false;
         h.loop = phi;
         h.normal = worldDir(nx, ny, nz);
         h.point = world(Sx, Sy, Sz);
         // the tube law (as the slot's): across 1/R1 = 1/R + 1/rho; along it the spine's curvature k = (B^2 / sin a) / |S'|^3 (the ellipse's), by Meusnier with the 3D cosine between n and
         // the spine's outward in-plane normal (-nin): -nn (the inner side: concave), concave bends held within 20 ball radii
         const double sp2 = g.B * g.B * (sp * sp / (g.s * g.s) + cp * cp), k = g.B * g.B / g.s / (sp2 * std::sqrt(sp2));
         const double kLong = RubberContact::TubeAlongCurvature(k, g.rho, -nn);
         const double R1 = 1.0 / (1.0 / R + 1.0 / g.rho), R2 = 1.0 / std::max(1.0 / R + kLong, 0.05 / R);
         h.reff = std::sqrt(R1 * R2);
         out.push_back(h);
      }
      // The plunger end, CLOSED (GrooveShape comment): four exact pieces, each counted only where it bounds wood (Gone: a cut covering the end removes them).
      //   face       the plane X = 0: the foot (0, Y, Z) for a centre on the groove's side (X > 0), when the foot moved rho into the wood, (-rho, Y, Z), lies within the floor's spine
      //              cylinder (radius B) and below the top rim's spine (Z < -rho); flat: Reff = R
      //   top rim    the tube about the line X = -rho, Z = -rho, |Y| < yh: the foot on the line, facing when the direction lies between +Z (top) and +X (face)
      //   floor rim  the tube about the plane X = -rho cut with the cylinder of radius B: the ellipse E = E0 + U sin(phi) + V cos(phi), E0 = (-rho, 0, Pz - rho s / c), U = (0, B, 0),
      //              V = (0, 0, -B / c), |phi| <= phiE. d/dphi |E - C|^2 / 2 = A cos - Bb sin + k sin cos (W = E0 - C, A = W.U, Bb = W.V, k = |U|^2 - |V|^2, U.V = 0); times (1 + t^2)^2,
      //              t = tan(phi / 2), the quartic -A t^4 - 2 (Bb + k) t^3 + 2 (k - Bb) t + A. Its critical points (the cubic -4A t^3 - 6 (Bb + k) t^2 + 2 (k - Bb), CubicRoots, 2 Newton
      //              steps) split |t| <= tE into monotone brackets; a sign change holds one root, placed by Newton held in its bracket to rounding. Facing: between the face's normal +X
      //              and the floor's -r. Curvature: the ellipse's, closed form. Two minima weighted by their hump (edge key ci << 16 | 0xFFFF), as a floor-wall edge
      //   corners    spheres of radius rho at (-rho, +-yh, -rho), where the three spines meet: a contact when the direction lies in the cone of +Z, +X and -r there; Reff 1/R + 1/rho
      {
         const double rho = g.rho, R1 = 1.0 / (1.0 / R + 1.0 / rho);
         const auto push = [&](const double px, const double py, const double pz, const double dist, const double tube, const double reff, const int edge, const double loop) {
            RimHit h;
            h.cut = ci;
            h.dist = dist;
            h.delta = R + tube - dist;
            h.rho = tube;
            h.floor = true;
            h.loop = loop;
            h.normal = worldDir((X - px) / dist, (Y - py) / dist, (Z - pz) / dist);
            h.point = world(px, py, pz);
            h.reff = reff;
            h.edge = edge;
            out.push_back(h);
         };
         const bool cut = !m.meets.empty();
         // face
         if (X > 1e-12 && X - R <= reachGap && Z < -rho)
         {
            const double wz = Z - g.Pz, al = -rho * g.c + wz * g.s;
            if (rho * rho + Y * Y + wz * wz - al * al < g.B * g.B && (!cut || !Gone(m, world(-rho, Y, Z))))
               push(0.0, Y, Z, X, 0.0, R, -1, 0.0);
         }
         // top rim
         if (std::abs(Y) < g.yh && X + rho >= 0.0 && Z + rho >= 0.0)
         {
            const double dist = std::hypot(X + rho, Z + rho);
            if (dist > 1e-12 && dist - (R + rho) <= reachGap && (!cut || !Gone(m, world(-rho, Y, -rho))))
               push(-rho, Y, -rho, dist, rho, std::sqrt(R1 * R), -1, 0.0);
         }
         // corners
         for (const double sg : { -1.0, 1.0 })
         {
            const double sy = sg * g.yh, dx = X + rho, dy = Y - sy, dz = Z + rho, dist = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (dist <= 1e-12 || dist - (R + rho) > reachGap)
               continue;
            double rx, ry, rz;
            g.Radial(-rho, sy, -rho, rx, ry, rz);
            const V3 n3[3] = { V3(0.0, 0.0, 1.0), V3(1.0, 0.0, 0.0), V3(-rx, -ry, -rz) };
            if (InCone(V3(dx, dy, dz), n3) && (!cut || !Gone(m, world(-rho, sy, -rho))))
               push(-rho, sy, -rho, dist, rho, R1, -1, 0.0);
         }
         // floor rim
         const double E0z = g.Pz - rho * g.s / g.c, Vz = -g.B / g.c, Wy = -Y, Wz = E0z - Z;
         const double A = Wy * g.B, Bb = Wz * Vz, k = g.B * g.B - Vz * Vz;
         const double c4 = -A, c3 = -2.0 * (Bb + k), c1 = 2.0 * (k - Bb), c0 = A;
         const auto Qf = [&](const double t) { return ((c4 * t + c3) * t * t + c1) * t + c0; };
         const auto Qd = [&](const double t) { return (4.0 * c4 * t + 3.0 * c3) * t * t + c1; };
         // a quick reach test: every point of the ellipse is within B / c of its centre (0 at phi = 0 lies at (-rho, 0, E0z))
         if (std::sqrt((X + rho) * (X + rho) + Y * Y + Wz * Wz) - g.B / g.c - (R + rho) > reachGap)
            return;
         double crit[3];
         int nc = 0;
         {
            const double a3 = 4.0 * c4, a2 = 3.0 * c3, a0 = c1;   // the cubic a3 t^3 + a2 t^2 + a0 (no linear term)
            if (std::abs(a3) * 1e3 > std::abs(a2))
            {
               const double sh = a2 / (3.0 * a3), P = -a2 * a2 / (3.0 * a3 * a3), Qc = 2.0 * sh * sh * sh + a0 / a3;   // t = x - sh: x^3 + P x + Qc
               double xr[3];
               const int n = CubicRoots(P, Qc, xr);
               for (int i = 0; i < n; i++)
                  crit[nc++] = xr[i] - sh;
            }
            else if (a2 != 0.0 && -a0 / a2 > 0.0)   // (the third root lies beyond |a2 / a3| > 1000: off the arc)
            {
               crit[nc++] = -std::sqrt(-a0 / a2);
               crit[nc++] = std::sqrt(-a0 / a2);
            }
            for (int i = 0; i < nc; i++)
               for (int it = 0; it < 2; it++)
               {
                  const double d2 = 3.0 * a3 * crit[i] * crit[i] + 2.0 * a2 * crit[i];
                  if (d2 == 0.0)
                     break;
                  crit[i] -= Qd(crit[i]) / d2;
               }
            std::sort(crit, crit + nc);
         }
         double br[5];
         int nb = 0;
         br[nb++] = -g.tE;
         for (int i = 0; i < nc; i++)
            if (crit[i] > br[nb - 1] && crit[i] < g.tE)
               br[nb++] = crit[i];
         br[nb++] = g.tE;
         const int key = ci << 16 | 0xFFFF;
         for (int i = 0; i + 1 < nb; i++)
         {
            const double lo = br[i], hi = br[i + 1], ql = Qf(lo), qh = Qf(hi);
            if ((ql < 0.0) == (qh < 0.0))
               continue;
            const bool isMin = ql < 0.0;   // (the distance falls, then rises along t, i.e. along Y)
            double a = lo, b = hi, t = 0.5 * (lo + hi);
            for (int it = 0; it < 60; it++)
            {
               const double q = Qf(t);
               if (q == 0.0)
                  break;
               ((q < 0.0) == isMin ? a : b) = t;
               const double dq = Qd(t);
               double next = dq != 0.0 ? t - q / dq : 0.5 * (a + b);
               if (!(next > a && next < b))
                  next = 0.5 * (a + b);
               const bool done = std::abs(next - t) <= 1e-16;
               t = next;
               if (done)
                  break;
            }
            const double phi = 2.0 * std::atan(t), sp = std::sin(phi), cp = std::cos(phi);
            const double ey = g.B * sp, ez = E0z + Vz * cp, ex = -rho;
            const double vx = X - ex, vy = Y - ey, vz = Z - ez, dist = std::sqrt(vx * vx + vy * vy + vz * vz);
            if (!isMin)
            {
               if (humps != nullptr)
                  humps->push_back({ ci, phi, dist, key });
               continue;
            }
            if (dist <= 1e-12 || dist - (R + rho) > reachGap)
               continue;
            double rx, ry, rz;
            g.Radial(ex, ey, ez, rx, ry, rz);
            const double nx = vx / dist, ny = vy / dist, nz = vz / dist, ab = -rx, na = nx, nb2 = -(nx * rx + ny * ry + nz * rz);
            if (na - ab * nb2 < 0.0 || nb2 - ab * na < 0.0)
               continue;
            if (cut && Gone(m, world(ex, ey, ez)))
               continue;
            // the ellipse's curvature vector: (E'' less its part along E') / |E'|^2, E' = (0, B cos, -Vz sin), E'' = (0, -B sin, -Vz cos)
            const double d1y = g.B * cp, d1z = -Vz * sp, d2y = -g.B * sp, d2z = -Vz * cp, sp2 = d1y * d1y + d1z * d1z, dd = (d2y * d1y + d2z * d1z) / sp2;
            const double ky = (d2y - dd * d1y) / sp2, kz = (d2z - dd * d1z) / sp2, kk = std::hypot(ky, kz);
            const double kLong = kk > 1e-9 ? RubberContact::TubeAlongCurvature(kk, rho, -(ny * ky + nz * kz) / kk) : 0.0;
            const double R2 = 1.0 / std::max(1.0 / R + kLong, 0.05 / R);
            push(ex, ey, ez, dist, rho, std::sqrt(R1 * R2), key, phi);
         }
      }
   }

   // A rectangle's contacts (playfieldcut.h): every true local minimum of the distance from the ball's centre to each exact piece, in closed form. Frame: s along the length, t across, half
   // sides hs, a, corner radius k, corner centres O = (+-(hs - k), +-(a - k)).
   //   rims   two tubes of radius rho, the top's (spine at z0 - rho) and the underside's (z0 - thick + rho), about the same level centre line: the outline offset outward by rho, 4 straight
   //          lines (side j at distance off_j + rho from the centre) and 4 quarter circles of radius k + rho about the O (k = 0: radius rho about the sharp corner, the Minkowski offset of
   //          the corner point, so the line stays tangent-continuous and the tube closes round the corner as a horn torus). The line is level, so the 3D extrema are the plane's: a side's
   //          normal foot (within the side), a quarter circle's nearest point O + r e (a minimum) or farthest O - r e (a hump), e the unit vector from O to the ball's centre, when it lies in
   //          the quarter. A minimum is a contact on its exposed quarter of the tube: facing the hole (the ball's centre inside the line, horizontally) and above the top rim's spine (below
   //          the underside rim's). Loop place: piece (side j = 2 j, corner j = 2 j + 1) + fraction, the underside rim's 8 more.
   //   walls  vertical: the outline's 4 planes and (k > 0) 4 concave quarter cylinders of radius k, between the two rims' spine heights (the rims join them tangent there): the foot on the
   //          plane, dist = off_j - (the centre's distance across), or the cylinder's foot on the ray from O, dist = k - |e|, with the centre inside the outline.
   void RectExtrema(const int ci, const double R, const V3 &c, std::vector<RimHit> &out, std::vector<RimHump> *humps, const double reachGap)
   {
      const Rim &m = g_rims[ci];
      const double dx = c.x - m.cx, dy = c.y - m.cy;
      const double cs = dx * m.lx + dy * m.ly, ct = dx * m.px + dy * m.py;
      const double As = m.hs - m.k, At = m.a - m.k, r = m.k + m.rho;
      const auto world = [&](const double s, const double t, const double z) { return V3(m.cx + m.lx * s + m.px * t, m.cy + m.ly * s + m.py * t, z); };
      const auto worldDir = [&](const double s, const double t) { return V3(m.lx * s + m.px * t, m.ly * s + m.py * t, 0.0); };
      // side j: outward normal N, tangent along the loop T, the outline's distance from the centre, the half extent of its straight part; corner j (after side j): signs of O
      static const double Ns[4] = { 0, 1, 0, -1 }, Nt[4] = { 1, 0, -1, 0 }, Ts[4] = { 1, 0, -1, 0 }, Tt[4] = { 0, -1, 0, 1 }, Os[4] = { 1, 1, -1, -1 }, Ot[4] = { 1, -1, -1, 1 };
      const double off[4] = { m.a, m.hs, m.a, m.hs }, ext[4] = { As, At, As, At };
      const double zTop = m.z0 - m.rho, zUnder = m.z0 - kThickM + m.rho;
      const auto rimHit = [&](const int level, const double loop, const double s, const double t, const double zl, const double kin, const double ds, const double dt) {
         // a minimum on the rim at level (0 top, 1 underside); kin the line's curvature toward the hole there, (ds, dt) the unit direction toward its centre of curvature
         const double dz = c.z - zl;
         if (level == 0 ? dz <= 0.0 : dz >= 0.0)   // below the top rim's spine (above the underside's): the wall's, not this tube's
            return;
         const V3 A = world(s, t, zl), dv = c - A;
         const double dist = dv.Length();
         if (dist <= 1e-12 || dist - (R + m.rho) > reachGap)
            return;
         if (!m.meets.empty() && Gone(m, A))   // inside another cut: no wood there (Gone)
            return;
         RimHit h;
         h.cut = ci;
         h.point = A;
         h.dist = dist;
         h.normal = dv / dist;
         h.delta = R + m.rho - dist;
         h.loop = loop + 8.0 * level;
         h.rho = m.rho;
         h.floor = false;
         // the tube law (as the slot's): across 1/R1 = 1/R + 1/rho; along it the line's curvature by Meusnier (TubeAlongCurvature), concave bends held within 20 ball radii
         double kLong = 0.0;
         if (kin > 0.0)
         {
            const V3 dn = worldDir(ds, dt);
            kLong = RubberContact::TubeAlongCurvature(kin, m.rho, -(h.normal.x * dn.x + h.normal.y * dn.y));
         }
         const double R1 = 1.0 / (1.0 / R + 1.0 / m.rho), R2 = 1.0 / std::max(1.0 / R + kLong, 0.05 / R);
         h.reff = std::sqrt(R1 * R2);
         out.push_back(h);
      };
      const auto wallHit = [&](const double dist, const double s, const double t, const double ns, const double nt, const double R1) {
         // a wall contact: dist from the centre to the wall, the foot (s, t) at the centre's height, the unit normal toward the centre, the law's radius across (along: vertical, R)
         if (dist <= 1e-12 || dist - R > reachGap)
            return;
         if (!m.meets.empty() && Gone(m, world(s - m.rho * ns, t - m.rho * nt, c.z)))   // the foot moved rho into the wood inside another cut: no wood there (Gone)
            return;
         RimHit h;
         h.cut = ci;
         h.point = world(s, t, c.z);
         h.dist = dist;
         h.normal = worldDir(ns, nt);
         h.delta = R - dist;
         h.loop = 0.0;
         h.rho = 0.0;
         h.floor = true;   // (weight 1 in the patch sum: a wall and a rim meet tangent, never both minima at one place)
         h.reff = std::sqrt(R1 * R);
         out.push_back(h);
      };
      const bool wall = c.z >= zUnder && c.z <= zTop;
      for (int j = 0; j < 4; j++)
      {
         // side j: along x, across y
         const double x = cs * Ts[j] + ct * Tt[j], y = cs * Ns[j] + ct * Nt[j];
         if (std::abs(x) <= ext[j])
         {
            const double fr = ext[j] > 0.0 ? 0.5 * (x + ext[j]) / ext[j] : 0.5, ss = x * Ts[j] + (off[j] + m.rho) * Ns[j], st = x * Tt[j] + (off[j] + m.rho) * Nt[j];
            if (y < off[j] + m.rho)   // facing the hole
               for (int level = 0; level < 2; level++)
                  rimHit(level, 2 * j + fr, ss, st, level == 0 ? zTop : zUnder, 0.0, 0.0, 0.0);
            if (wall && y < off[j])
               wallHit(off[j] - y, x * Ts[j] + off[j] * Ns[j], x * Tt[j] + off[j] * Nt[j], -Ns[j], -Nt[j], R);
         }
         // corner j: e from O
         const double os = Os[j] * As, ot = Ot[j] * At, es = cs - os, et = ct - ot, n = std::hypot(es, et);
         if (n <= 0.0)
            continue;
         const double us = es / n, ut = et / n;
         const auto frac = [&](const double vs, const double vt) { return std::atan2(std::abs(Ns[j] * vt - Nt[j] * vs), Ns[j] * vs + Nt[j] * vt) / (0.5 * M_PI); };   // from side j's normal
         if (es * Os[j] > 0.0 && et * Ot[j] > 0.0)   // the nearest point in the quarter: a minimum
         {
            if (n < r)
               for (int level = 0; level < 2; level++)
                  rimHit(level, 2 * j + 1 + frac(us, ut), os + r * us, ot + r * ut, level == 0 ? zTop : zUnder, 1.0 / r, -us, -ut);
            if (wall && m.k > 0.0 && n < m.k)
               wallHit(m.k - n, os + m.k * us, ot + m.k * ut, -us, -ut, 1.0 / std::max(1.0 / R - 1.0 / m.k, 0.05 / R));
         }
         else if (humps != nullptr && -es * Os[j] > 0.0 && -et * Ot[j] > 0.0)   // the farthest point in the quarter: a hump, on each rim
         {
            const double hs2 = (n + r) * (n + r);
            for (int level = 0; level < 2; level++)
            {
               const double dz = c.z - (level == 0 ? zTop : zUnder);
               humps->push_back({ ci, 2 * j + 1 + frac(-us, -ut) + 8.0 * level, std::sqrt(hs2 + dz * dz) });
            }
         }
      }
   }

   // the cuts of this host within reach of a ball at c, in index order. cq null (shared_finder = 0): today's scan of every cut.
   // shared_finder P4b (cq set): the ball's query holds every cut whose box is within cq->r of cq->c0; a cut within reachM of c is within |c - c0| + reachM of c0 (triangle inequality),
   // so the query's list is complete when that is inside its volume. Otherwise the cuts alone are queried again around c0 with reachM + 2 |c - c0| (design.md section 4's inside check:
   // the factor 2 bounds repeats), never scanned. No query yet or a stale one: queried around c within reachM.
   void NearCuts(const ContactSurface *host, const V3 &c, const double reachM, std::vector<int> &out, PlayfieldCuts::CutQuery *const cq)
   {
      const bool prof = ContactProf::CheckOn() || CutCheckOn();
      const uint64_t t0 = prof ? __rdtsc() : 0;
      out.clear();
      size_t tested = 0;
      if (cq == nullptr)
      {
         for (size_t i = 0; i < g_rims.size(); i++)
            if (g_rims[i].host == host && g_rims[i].BoxDist(c) <= reachM)
               out.push_back((int)i);
         tested = g_rims.size();
      }
      else
      {
         const bool valid = cq->r >= 0.0 && cq->ver == g_rimsVer;
         const double stray = valid ? (c - cq->c0).Length() : 0.0;
         if (!valid || stray + reachM > cq->r)
         {
            if (valid)
               ContactSurfaceHost::FinderCuts(*cq, cq->c0, reachM + 2.0 * stray);
            else
               ContactSurfaceHost::FinderCuts(*cq, c, reachM);
            if (prof)
               g_cutStats.requery++;
         }
         for (const int i : cq->cuts)
            if (g_rims[i].host == host && g_rims[i].BoxDist(c) <= reachM)
               out.push_back(i);
         tested = cq->cuts.size();
      }
      if (prof)
      {
         g_cutStats.calls++;
         g_cutStats.tested += tested;
         g_cutStats.ticks += __rdtsc() - t0;
      }
      if (cq != nullptr && CutCheckOn())   // check mode: today's scan, the same order
      {
         static std::vector<int> full;
         full.clear();
         for (size_t i = 0; i < g_rims.size(); i++)
            if (g_rims[i].host == host && g_rims[i].BoxDist(c) <= reachM)
               full.push_back((int)i);
         g_cutStats.checks++;
         if (full != out)
         {
            g_cutStats.mismatch++;
            for (const int i : full)
               if (std::find(out.begin(), out.end(), i) == out.end())
                  g_cutStats.miss++;
            static int logged = 0;
            if (logged++ < 20)
               PLOGW << "Finder cut check: NearCuts differs from the scan at c " << c.x << " " << c.y << " " << c.z << ", reach " << reachM << ": " << out.size() << " against " << full.size()
                     << " cuts; query c0 " << cq->c0.x << " " << cq->c0.y << " " << cq->c0.z << " r " << cq->r;
         }
      }
   }

   // the faces whose nearest point to c is surface: one inside a near cut's centre line (at or under the top) is not (as CarvedHole::Holds for a turntable)
   // shared_finder P3a: fr reads the faces' measurements from the ball's table (faces is then the query list), and outPos gets each kept face's position in it
   void FacesLeft(const std::vector<SurfaceTri> &tris, const std::vector<int> &faces, const std::vector<int> &cuts, const V3 &c, std::vector<int> &out, const RubberContact::FaceRead *fr = nullptr,
                  std::vector<int> *outPos = nullptr)
   {
      out.clear();
      if (outPos != nullptr)
         outPos->clear();
      for (size_t j = 0; j < faces.size(); j++)
      {
         const int i = faces[j];
         int kind;
         const V3 p = fr != nullptr ? fr->At(j, tris[i], c).p : RubberContact::ClosestOnTriangle(tris[i], c, kind);
         if (std::none_of(cuts.begin(), cuts.end(), [&](const int ci) { return g_rims[ci].Holds(p); }))
         {
            out.push_back(i);
            if (outPos != nullptr)
               outPos->push_back(fr != nullptr && fr->pos != nullptr ? fr->pos[j] : (int)j);
         }
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

   // A [PlayfieldCut:<name>] shape = groove: its keys (playfieldcut.h), its playfield host and the groove's numbers; R a standard ball for the logged paper numbers
   void SetupGroove(const CutKeys &k, const std::vector<ContactSurfaceHost *> &pf, const double R)
   {
      const std::string where = "Playfield cut " + k.name + ": ";
      if (k.widthSet)
         Notice(where + "width_mm is not a groove key (its radius is fixed): ignored (physics ini)");
      double lengthMm = k.lengthMm;
      if (lengthMm > kGrooveLenMaxMm || lengthMm < kGrooveLenMinMm)
      {
         const double lim = lengthMm > kGrooveLenMaxMm ? kGrooveLenMaxMm : kGrooveLenMinMm;
         std::ostringstream o;
         o << where << "groove length_mm " << lengthMm << " is outside " << kGrooveLenMinMm << " .. " << kGrooveLenMaxMm << ": using " << lim << " (physics ini)";
         Notice(o.str());
         lengthMm = lim;
      }
      ContactSurfaceHost *host = nullptr;   // the highest playfield surface at the plunger end
      double zU = 0.0;
      for (ContactSurfaceHost *const h : pf)
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
         Notice(where + "no playfield surface at its plunger end: the groove is ignored (physics ini x, y)");
         return;
      }
      Rim m;
      m.groove = true;
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
      m.rho = kEdgeRhoM;
      GrooveShape &g = m.gr;
      g.Build(lengthMm * 1e-3);
      // the playfield must be flat over the groove and 1 mm round it (the rim's top tangent line is the top's edge): its top height on a 1 mm grid
      bool flat = true;
      for (double s = -1e-3; flat && s <= g.xv + 1e-3 + 1e-9; s += 1e-3)
         for (double t = -(g.yh + 1e-3); flat && t <= g.yh + 1e-3 + 1e-9; t += 1e-3)
         {
            double z;
            flat = host->SurfaceHeightAt((m.cx + m.lx * s + m.px * t) / U, (m.cy + m.ly * s + m.py * t) / U, z) && std::abs(z - zU) * U <= 1e-6;
         }
      if (!flat)
      {
         Notice(where + "the playfield is not flat over the groove (a raised or curved playfield, or its edge): the groove is ignored (physics ini x, y, angle, length_mm)");
         return;
      }
      // the world box: X from the plunger end to the spine's vertex, |Y| <= the half-width at X = 0, Z from the floor's lowest point (X = 0) to the top
      for (int i = 0; i < 3; i++)
      {
         m.lo[i] = 1e30;
         m.hi[i] = -1e30;
      }
      for (int corner = 0; corner < 8; corner++)
      {
         const double s = corner & 1 ? g.xv : -2.0 * g.rho, t = (corner & 2 ? 1.0 : -1.0) * g.yh, z = corner & 4 ? 0.0 : -g.D0;   // (the closed end's pieces reach X = -rho)
         const double v[3] = { m.cx + m.lx * s + m.px * t, m.cy + m.ly * s + m.py * t, m.z0 + z };
         for (int i = 0; i < 3; i++)
         {
            m.lo[i] = std::min(m.lo[i], v[i]);
            m.hi[i] = std::max(m.hi[i], v[i]);
         }
      }
      // logged paper numbers (paper sections 1, 2, 3a): a ball resting on the floor's bottom 20 mm from the plunger end has its centre Rg - R from the axis, straight below it
      g.floorR1 = 1.0 / (1.0 / R - 1.0 / g.Rg);
      const double sax = (20e-3 - (g.Rg - R) * g.s) / g.c, zc = g.Pz + sax * g.s - (g.Rg - R) * g.c;
      std::ostringstream o;
      o << "Playfield cut " << k.name << ": groove in the playfield " << host->Name() << " (exact floor and rim), plunger end " << k.x << ", " << k.y << " (playfield z " << zU << "), angle " << k.angle
        << " deg (toward the tip " << m.lx << ", " << m.ly << "), length " << lengthMm << " mm, depth " << g.D0 * 1e3 << " mm at the plunger end (open) to 0 at the tip, radius " << g.Rg * 1e3
        << " mm, rim radius " << m.rho * 1e3 << " mm (assumed); tilt " << g.a * 180.0 / M_PI << " deg, opening " << 2.0 * g.yh * 1e3 << " mm at the plunger end, spine vertex at "
        << g.xv * 1e3 << " mm, its radius " << g.B * g.s * 1e6 << " um. Paper, 27 mm ball: floor Hertz radius " << std::sqrt(g.floorR1 * R) * 1e3 << " mm; resting 20 mm from the plunger end, centre "
        << zc * 1e3 << " mm above the top";
      PLOGI << o.str();
      g_rims.push_back(std::move(m));
   }

   // A [PlayfieldCut:<name>] shape = rectangle: its keys (playfieldcut.h), its playfield host and its numbers; R a standard ball for the logged paper numbers
   void SetupRect(const CutKeys &k, const std::vector<ContactSurfaceHost *> &pf, const double R)
   {
      const std::string where = "Playfield cut " + k.name + ": ";
      if (!k.widthSet || k.widthMm <= 0.0 || k.lengthMm <= 0.0)
      {
         Notice(where + "a rectangle needs length_mm and width_mm above 0: the cut is ignored (physics ini)");
         return;
      }
      double cornerMm = std::max(k.cornerMm, 0.0);
      const double cMax = 0.5 * std::min(k.lengthMm, k.widthMm);
      if (cornerMm > cMax || k.cornerMm < 0.0)
      {
         std::ostringstream o;
         o << where << "corner_mm " << k.cornerMm << " is outside 0 .. " << cMax << " (half the shorter side): using " << std::clamp(k.cornerMm, 0.0, cMax) << " (physics ini)";
         Notice(o.str());
         cornerMm = std::clamp(k.cornerMm, 0.0, cMax);
      }
      ContactSurfaceHost *host = nullptr;   // the highest playfield surface at the centre
      double zU = 0.0;
      for (ContactSurfaceHost *const h : pf)
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
         return;
      }
      Rim m;
      m.rect = true;
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
      m.rho = kEdgeRhoM;
      m.hs = 0.5 * k.lengthMm * 1e-3;
      m.a = 0.5 * k.widthMm * 1e-3;
      m.k = cornerMm * 1e-3;
      // the playfield must be flat over the cut and 1 mm round it (the rim's top tangent line is the top's edge): its top height on a 0.5 mm grid
      bool flat = true;
      const double sB = m.hs + m.rho + 1e-3, tB = m.a + m.rho + 1e-3;
      for (double s = -sB; flat && s <= sB + 1e-9; s += 0.5e-3)
         for (double t = -tB; flat && t <= tB + 1e-9; t += 0.5e-3)
         {
            double z;
            flat = host->SurfaceHeightAt((m.cx + m.lx * s + m.px * t) / U, (m.cy + m.ly * s + m.py * t) / U, z) && std::abs(z - zU) * U <= 1e-6;
         }
      if (!flat)
      {
         Notice(where + "the playfield is not flat over the cut (a raised or curved playfield, or its edge): the cut is ignored (physics ini x, y, angle, length_mm, width_mm)");
         return;
      }
      // the world box: the rim centre line's corners (|s| <= hs + rho, |t| <= a + rho), from the underside to the top
      for (int i = 0; i < 3; i++)
      {
         m.lo[i] = 1e30;
         m.hi[i] = -1e30;
      }
      for (int corner = 0; corner < 8; corner++)
      {
         const double s = (corner & 1 ? 1.0 : -1.0) * (m.hs + m.rho), t = (corner & 2 ? 1.0 : -1.0) * (m.a + m.rho), z = corner & 4 ? 0.0 : -kThickM;
         const double v[3] = { m.cx + m.lx * s + m.px * t, m.cy + m.ly * s + m.py * t, m.z0 + z };
         for (int i = 0; i < 3; i++)
         {
            m.lo[i] = std::min(m.lo[i], v[i]);
            m.hi[i] = std::max(m.hi[i], v[i]);
         }
      }
      // logged paper numbers: a ball narrower than the hole's width (2 (a + rho) < 2 (R + rho)) rests on the two long rims, its centre hc over the top (as the slot)
      const double ar = m.a + m.rho;
      std::ostringstream o;
      o << "Playfield cut " << k.name << ": rectangle in the playfield " << host->Name() << " (exact rims, walls and underside rim), centre " << k.x << ", " << k.y << " (playfield z " << zU
        << "), angle " << k.angle << " deg (along " << m.lx << ", " << m.ly << "), " << k.lengthMm << " x " << k.widthMm << " mm, corner " << cornerMm << " mm, through " << kThickM * 1e3
        << " mm, rim radius " << m.rho * 1e3 << " mm (assumed). Paper, 27 mm ball: ";
      if (ar < R + m.rho)
         o << "rests on the long rims, centre " << (std::sqrt((R + m.rho) * (R + m.rho) - ar * ar) - m.rho) * 1e3 << " mm above the top, rim normals " << std::asin(ar / (R + m.rho)) * 180.0 / M_PI
           << " deg from vertical";
      else
         o << "falls through (the hole is " << 2.0 * (m.a - R) * 1e3 << " mm wider than the ball)";
      PLOGI << o.str();
      g_rims.push_back(std::move(m));
   }

   // the places (piece, p) along through cut m's centre line where in(x, y) changes, in loop order; returns in's value at the line's start. The line sampled every 0.05 mm along (at
   // least 16 samples per piece), each change bisected to 1e-13 m along the line (a stated solve: two changes closer than 0.05 mm along the line, or a tangent touch, are not seen)
   struct Cross
   {
      int piece;
      double p;
   };
   template <class F> bool Walk(const Rim &m, const F &in, std::vector<Cross> &out)
   {
      out.clear();
      double x, y, dx, dy;
      LineAt(m, 0, 0.0, x, y, dx, dy);
      const bool first = in(x, y);
      bool prev = first;
      for (int k = 0; k < kPieces; k++)
      {
         const double len = PieceLen(m, k);
         if (len <= 0.0)
            continue;
         double sp = 0.0;   // the largest speed |dP/dp| of three on the piece
         for (const double fr : { 0.0, 0.5, 1.0 })
         {
            LineAt(m, k, fr * len, x, y, dx, dy);
            sp = std::max(sp, std::hypot(dx, dy));
         }
         const int n = std::max(16, (int)std::ceil(len * sp / 0.05e-3));
         double pp = 0.0;
         for (int i = 1; i <= n; i++)
         {
            const double p = len * i / n;
            LineAt(m, k, p, x, y, dx, dy);
            const bool cur = in(x, y);
            if (cur != prev)
            {
               double lo = pp, hi = p;
               for (int it = 0; it < 200 && (hi - lo) * sp > 1e-13; it++)
               {
                  const double mid = 0.5 * (lo + hi);
                  LineAt(m, k, mid, x, y, dx, dy);
                  (in(x, y) == prev ? lo : hi) = mid;
               }
               out.push_back({ k, 0.5 * (lo + hi) });
            }
            prev = cur;
            pp = p;
         }
      }
      return first;
   }

   // Cuts that meet (design runs/2026-11-10-cut-junction sections 1 - 2), found once at load, one rule for every pair (no per-pair code): every cut's meets (Rim::meets: the cuts on its
   // host whose boxes come within rho of its own; its pieces then count only on wood, Gone) and the rims and corners where the cuts' surfaces meet (Join, carried by the owner's probe):
   //   two through cuts (slot, rectangle: vertical walls, so their overlap is the union of their outlines): where their centre lines cross, V (rho from both walls), a vertical edge
   //            (kind 0) and its two corner spheres (kind 2: top and the two walls at the top rim's spine height; underside and the two walls at the underside rim's)
   //   a groove and a through cut: the edge where the floor meets the wall (kind 1) along the part of the wall's centre line inside the groove's spine outline (there the floor lies under
   //            the top's rim height); where the centre line crosses the spine outline, a corner sphere (top, wall, floor)
   //   two grooves (a floor meeting a floor): geometry the rule does not describe (a notice; each solved alone there)
   // A junction inside a third cut is left out. A cut that meets no other has no meets and no joins: every piece as before.
   void Meets()
   {
      g_joins.clear();
      const int n = (int)g_rims.size();
      const auto inThird = [&](const V3 &p, const int i, const int j) {
         for (int k = 0; k < n; k++)
            if (k != i && k != j && g_rims[k].host == g_rims[i].host && g_rims[k].Near(p))
               return true;
         return false;
      };
      const auto grow = [](Rim &m, const V3 &p) {
         const double v[3] = { p.x, p.y, p.z };
         for (int i = 0; i < 3; i++)
         {
            m.lo[i] = std::min(m.lo[i], v[i]);
            m.hi[i] = std::max(m.hi[i], v[i]);
         }
      };
      const auto add = [&](Join J) {
         J.owner = J.kind == 1 ? J.a : g_rims[J.b].rect && !g_rims[J.a].rect && !g_rims[J.a].groove ? J.b : J.a;
         g_rims[J.owner].joins.push_back((int)g_joins.size());
         g_joins.push_back(J);
      };
      std::vector<Cross> ca, cb;
      for (int i = 0; i < n; i++)
         for (int j = i + 1; j < n; j++)
         {
            Rim &p = g_rims[i], &q = g_rims[j];
            const double pad = p.rho + q.rho;
            if (p.host != q.host || p.lo[0] > q.hi[0] + pad || q.lo[0] > p.hi[0] + pad || p.lo[1] > q.hi[1] + pad || q.lo[1] > p.hi[1] + pad)
               continue;
            p.meets.push_back(j);
            q.meets.push_back(i);
            if (p.groove && q.groove)
            {
               Notice("Playfield cuts " + p.name + " and " + q.name + ": two grooves meeting (a floor meeting a floor) is not described by the cuts-meet rule: each is solved alone there (physics ini)");
               continue;
            }
            int nEdge = 0, nCorner = 0;
            if (!p.groove && !q.groove)   // two through cuts: their centre lines' crossings
            {
               const bool a0 = Walk(p, [&](const double x, const double y) { return q.Near(V3(x, y, q.z0 - q.rho)); }, ca);
               Walk(q, [&](const double x, const double y) { return p.Near(V3(x, y, p.z0 - p.rho)); }, cb);
               (void)a0;
               for (const Cross &x : ca)
               {
                  double vx, vy, dx, dy, bx, by, bdx, bdy, best = 1e30;
                  LineAt(p, x.piece, x.p, vx, vy, dx, dy);
                  V3 nb;
                  for (const Cross &y : cb)
                  {
                     LineAt(q, y.piece, y.p, bx, by, bdx, bdy);
                     const double d = std::hypot(bx - vx, by - vy), l = std::hypot(bdx, bdy);
                     if (d < best)
                     {
                        best = d;
                        nb = V3(bdy / l, -bdx / l, 0.0);
                     }
                  }
                  const double la = std::hypot(dx, dy);
                  const V3 na(dy / la, -dx / la, 0.0), top(vx, vy, p.z0 - p.rho), under(vx, vy, p.z0 - kThickM + p.rho);
                  if (best > 1e-7 || inThird(top, i, j))   // (no matching crossing on the other line within 0.1 um: not a crossing of the two)
                     continue;
                  Join J;
                  J.a = i;
                  J.b = j;
                  J.kind = 0;
                  J.S = top;
                  J.n[0] = na;
                  J.n[1] = nb;
                  add(J);
                  J.kind = 2;
                  J.n[0] = V3(0.0, 0.0, 1.0);
                  J.n[1] = na;
                  J.n[2] = nb;
                  add(J);
                  J.S = under;
                  J.n[0] = V3(0.0, 0.0, -1.0);
                  add(J);
                  grow(g_rims[g_joins.back().owner], V3(vx, vy, p.z0));
                  grow(g_rims[g_joins.back().owner], V3(vx, vy, p.z0 - kThickM));
                  nEdge++;
                  nCorner += 2;
               }
            }
            else   // a groove and a through cut: the wall's centre line inside the groove's spine outline
            {
               const int gi = p.groove ? i : j, wi = p.groove ? j : i;
               Rim &g = g_rims[gi];
               const Rim &w = g_rims[wi];
               const auto inG = [&](const double x, const double y) {
                  const double dx = x - g.cx, dy = y - g.cy;
                  return g.gr.InsideOutline(dx * g.lx + dy * g.ly, dx * g.px + dy * g.py);
               };
               bool in = Walk(w, inG, ca);
               size_t ic = 0;
               for (int k = 0; k < kPieces; k++)
               {
                  const double len = PieceLen(w, k);
                  if (len <= 0.0)
                     continue;
                  double start = 0.0;
                  for (;; ic++)
                  {
                     const bool more = ic < ca.size() && ca[ic].piece == k;
                     const double end = more ? ca[ic].p : len;
                     if (in && end > start)
                     {
                        Join J;
                        J.kind = 1;
                        J.a = gi;
                        J.b = wi;
                        J.piece = k;
                        J.p0 = start;
                        J.p1 = end;
                        EdgeBounds(g, w, J);
                        add(J);
                        nEdge++;
                     }
                     if (!more)
                        break;
                     // the crossing: a corner sphere at (P, -rho), normals +Z (top), a (the wall's, into the through cut), b = -r (the floor's, toward the axis)
                     double x, y, dx, dy, rx, ry, rz;
                     LineAt(w, k, ca[ic].p, x, y, dx, dy);
                     const double l = std::hypot(dx, dy), X = (x - g.cx) * g.lx + (y - g.cy) * g.ly, Y = (x - g.cx) * g.px + (y - g.cy) * g.py;
                     g.gr.Radial(X, Y, -g.rho, rx, ry, rz);
                     const V3 S(x, y, g.z0 - g.rho);
                     if (!inThird(S, gi, wi))
                     {
                        Join J;
                        J.kind = 2;
                        J.a = gi;
                        J.b = wi;
                        J.S = S;
                        J.n[0] = V3(0.0, 0.0, 1.0);
                        J.n[1] = V3(dy / l, -dx / l, 0.0);
                        J.n[2] = V3(-(g.lx * rx + g.px * ry), -(g.ly * rx + g.py * ry), -rz);
                        J.owner = gi;
                        g.joins.push_back((int)g_joins.size());
                        g_joins.push_back(J);
                        nCorner++;
                     }
                     in = !in;
                     start = ca[ic].p;
                  }
               }
               if (nEdge > 0)
                  g.lo[2] -= g.rho;   // the world box reaches the edge tube's bottom, rho under the floor's
            }
            std::ostringstream o;
            o << "Playfield cuts " << p.name << " and " << q.name << " meet (cuts-meet rule): " << nEdge << (p.groove || q.groove ? " floor-wall edge piece(s), " : " vertical edge(s), ") << nCorner
              << " corner sphere(s); each cut's rims, walls and floor end where the other has removed the wood";
            PLOGI << o.str();
         }
   }

}

namespace PlayfieldCuts
{
   void Setup(PhysicsEngine *, PinTable *const table)
   {
      g_rims.clear();
      g_rimsVer++;   // shared_finder P4b: every ball's cut query is stale
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
      std::vector<ContactSurfaceHost *> pf;   // the playfield's hosts (the plane's, or a playfield_mesh primitive's): a host with faces has a height under the cut
      for (ContactHost *const h : ContactHost::AllHosts())
         if (h->IsPlayfieldKind())
            if (ContactSurfaceHost *const r = dynamic_cast<ContactSurfaceHost *>(h); r != nullptr)
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
         if (k.shape != "slot" && k.shape != "groove" && k.shape != "rectangle")
         {
            Notice(where + (k.shape.empty() ? std::string("no shape") : "shape '" + k.shape + "' is not known") + ": the cut is ignored. The shapes are slot, groove and rectangle (physics ini)");
            continue;
         }
         if (!k.xSet || !k.ySet || !k.lengthSet || k.lengthMm < 0.0)
         {
            Notice(where + "needs x, y and length_mm (0 or more): the cut is ignored (physics ini)");
            continue;
         }
         if (k.shape == "groove")
         {
            SetupGroove(k, pf, R);
            continue;
         }
         if (k.shape == "rectangle")
         {
            SetupRect(k, pf, R);
            continue;
         }
         if (k.cornerSet)
            Notice(where + "corner_mm is a rectangle key: ignored (physics ini)");
         double widthMm = k.widthMm;
         if (widthMm > kWidthMaxMm || widthMm < kWidthMinMm)
         {
            const double lim = widthMm > kWidthMaxMm ? kWidthMaxMm : kWidthMinMm;
            std::ostringstream o;
            o << where << "width_mm " << widthMm << " is outside " << kWidthMinMm << " .. " << kWidthMaxMm << ": using " << lim << " (physics ini)";
            Notice(o.str());
            widthMm = lim;
         }
         ContactSurfaceHost *host = nullptr;
         double zU = 0.0;
         for (ContactSurfaceHost *const h : pf)
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
      Meets();
   }

   bool Any(const ContactSurface *const host)
   {
      if (!ContactSurfaceHost::FinderDrives())   // shared_finder = 0: today's scan
         return std::any_of(g_rims.begin(), g_rims.end(), [&](const Rim &m) { return m.host == host; });
      if (g_anyVer != g_rimsVer)   // P4b: the hosts holding cuts, once per Setup
      {
         g_anyHosts.clear();
         for (const Rim &m : g_rims)
            if (std::find(g_anyHosts.begin(), g_anyHosts.end(), m.host) == g_anyHosts.end())
               g_anyHosts.push_back(m.host);
         g_anyVer = g_rimsVer;
      }
      return std::find(g_anyHosts.begin(), g_anyHosts.end(), host) != g_anyHosts.end();
   }

   uint64_t Version() { return g_rimsVer; }
   int Count() { return (int)g_rims.size(); }
   void Box(const int i, V3 &lo, V3 &hi)
   {
      lo = V3(g_rims[i].lo[0], g_rims[i].lo[1], g_rims[i].lo[2]);
      hi = V3(g_rims[i].hi[0], g_rims[i].hi[1], g_rims[i].hi[2]);
   }
   CutStats &Stats() { return g_cutStats; }

   RubberContact::TubeProbe Probe(const ContactSurface *const host, const RubberContact::BallSI &bs, const std::vector<SurfaceTri> &tris, const std::vector<int> &faces, const double edgeRho,
                                  const V3 &c, const double reachM, const RubberContact::CarvedHole *const holes, const int nHoles, const bool patch, const double zero, RubberContact::FaceTable *const ft, CutQuery *const cq)
   {
      static std::vector<int> cuts, left, leftPos;   // (the physics thread's; reused)
      const RubberContact::FaceRead frF { ft, nullptr };   // shared_finder P3a: the ball's table (null: measured as before)
      const RubberContact::FaceRead *const pF = ft != nullptr ? &frF : nullptr;
      static std::vector<RimHit> hits;
      static std::vector<RimHump> humps;
      NearCuts(host, c, reachM, cuts, cq);
      // shared_finder P3b: the near pass's faces within R (ContactSurfaceHost::FinderTable), in load order with their query-list positions; the probe with none within R reads the whole list with the
      // box bound (its nearest wherever it is, as today). A turntable's hole (nHoles) moves a face's nearest point, so there the whole list as P3a.
      const bool nearOn = ft != nullptr && ft->useNear && nHoles == 0;
      const std::vector<int> &in = nearOn ? ft->inI : faces;
      const RubberContact::FaceRead frI { ft, nearOn ? ft->inJ.data() : nullptr };
      const RubberContact::FaceRead *const pI = ft != nullptr ? &frI : nullptr;
      if (cuts.empty())   // today's probe
      {
         RubberContact::TubeProbe pr;
         if (nearOn && in.empty())
         {
            const RubberContact::FaceRead frP { ft, nullptr, true };
            pr = RubberContact::ProbeTriangles(bs, tris, faces, edgeRho, c, holes, nHoles, 0.0, nullptr, &frP);
         }
         else
            pr = RubberContact::ProbeTriangles(bs, tris, in, edgeRho, c, holes, nHoles, 0.0, nullptr, pI);
         if (patch)
            RubberContact::PatchTriangles(bs, tris, in, c, zero, pr, holes, nHoles, pI);   // (the patch reads only faces within R)
         return pr;
      }
      FacesLeft(tris, in, cuts, c, left, pI, ft != nullptr ? &leftPos : nullptr);   // P3b: FacesLeft filters the near pass's set
      const bool wide = nearOn && left.empty();   // P3b: no face left within R: the probe's nearest of every face left, with the box bound
      if (wide)
         FacesLeft(tris, faces, cuts, c, left, pF, &leftPos);
      const RubberContact::FaceRead frL { ft, leftPos.data(), wide };   // the faces left, read at their query-list positions
      const RubberContact::FaceRead *const pL = ft != nullptr ? &frL : nullptr;
      RubberContact::TubeProbe pr = RubberContact::ProbeTriangles(bs, tris, left, edgeRho, c, holes, nHoles, 0.0, nullptr, pL);
      // The reach (beyond S = R + rho): a rim minimum matters here only if it is deeper than the faces' nearest (it becomes the probe) or presses (d = delta - zero > 0 within S - zero,
      // the patch sum), so with a face left the rim is searched to S + max(the faces' gap, -zero (patch), 0). With no face the deepest minimum is the probe wherever it is: no limit.
      // Every extremum is in closed form, so one search costs the same whatever its reach. Every hump is kept at its distance (one beyond the top gives s = 0).
      const bool face = pr.seg >= 0;
      const double extra = face ? std::max(std::max(-pr.delta, patch ? -zero : 0.0), 0.0) : kNoLimitM;
      hits.clear();
      humps.clear();
      for (const int ci : cuts)
         RimExtrema(ci, bs.R, c, hits, patch ? &humps : nullptr, extra + kReachPadM);
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
            RubberContact::PatchTriangles(bs, tris, left, c, zero, pr, holes, nHoles, pL);
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
      pr.rho = b.rho;
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
               double hump = -1.0;
               if (h1.edge >= 0 || h2.edge >= 0)   // a floor-wall edge (Join kind 1): two minima on one edge as the slot's (the through cut's loop, its humps); with any other piece weight 1
               {                                   // (they meet tangent)
                  if (h1.edge != h2.edge)
                     continue;
                  const double fwd = std::fmod(h2.loop - h1.loop + 2.0 * tot, tot);
                  double m1 = -1.0, m2 = -1.0;
                  for (const RimHump &q : humps)
                  {
                     if (q.cut != ci || q.edge != h1.edge)
                        continue;
                     double &mm = std::fmod(q.loop - h1.loop + 2.0 * tot, tot) < fwd ? m1 : m2;
                     mm = std::max(mm, q.dist);
                  }
                  hump = m1 < 0.0 ? m2 : m2 < 0.0 ? m1 : std::min(m1, m2);
               }
               else if (g_rims[ci].rect)   // a rectangle: a wall contact keeps weight 1 (a wall and a rim meet tangent), as do contacts on its two different rims; two on one rim: as the slot's,
               {                      // on that rim's loop (8 pieces each, the underside's 8 on)
                  if (h1.floor || h2.floor || (h1.loop >= 8.0) != (h2.loop >= 8.0))
                     continue;
                  const double base = h1.loop >= 8.0 ? 8.0 : 0.0, l1 = h1.loop - base, fwd = std::fmod(h2.loop - h1.loop + 2.0 * tot, tot);
                  double m1 = -1.0, m2 = -1.0;
                  for (const RimHump &q : humps)
                  {
                     if (q.cut != ci || q.edge >= 0 || (q.loop >= 8.0) != (base > 0.0))
                        continue;
                     double &mm = std::fmod(q.loop - base - l1 + 2.0 * tot, tot) < fwd ? m1 : m2;
                     mm = std::max(mm, q.dist);
                  }
                  hump = m1 < 0.0 ? m2 : m2 < 0.0 ? m1 : std::min(m1, m2);
               }
               else if (g_rims[ci].groove)   // a groove: a floor contact keeps weight 1 (its pieces meet tangent: no hump between, and never two pressing at once, paper section 3e); two rim
               {                             // contacts: the largest distance on the open spine arc between them
                  if (h1.floor || h2.floor)
                     continue;
                  for (const RimHump &q : humps)
                     if (q.cut == ci && q.edge < 0 && q.loop > std::min(h1.loop, h2.loop) && q.loop < std::max(h1.loop, h2.loop))
                        hump = std::max(hump, q.dist);
               }
               else   // a slot: the largest distance on each arc between them (forward from h1 to h2, and the rest); the hump is the lower of the two (a join's contact: weight 1)
               {
                  if (h1.floor || h2.floor)
                     continue;
                  const double fwd = std::fmod(h2.loop - h1.loop + 2.0 * tot, tot);
                  double m1 = -1.0, m2 = -1.0;
                  for (const RimHump &q : humps)
                  {
                     if (q.cut != ci || q.edge >= 0)
                        continue;
                     double &mm = std::fmod(q.loop - h1.loop + 2.0 * tot, tot) < fwd ? m1 : m2;
                     mm = std::max(mm, q.dist);
                  }
                  hump = m1 < 0.0 ? m2 : m2 < 0.0 ? m1 : std::min(m1, m2);
               }
               const double d1 = h1.delta - zero, d2 = h2.delta - zero;
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
                   const V3 &v, const double a, const double hMax, const double touchM, const double arriveM, const double reachM, RubberContact::FaceTable *const ft, CutQuery *const cq)
   {
      static std::vector<int> cuts, left, leftPos;
      static std::vector<RimHit> hits;
      const RubberContact::FaceRead frF { ft, nullptr };   // shared_finder P3a: the ball's table (null: measured as before)
      const RubberContact::FaceRead *const pF = ft != nullptr ? &frF : nullptr;
      NearCuts(host, c, reachM, cuts, cq);
      if (cuts.empty())   // today's bound
         return RubberContact::ArrivalTriangles(bs, tris, faces, nearestSeg, c, v, a, hMax, touchM, arriveM, pF);
      const int nT = (int)tris.size();
      FacesLeft(tris, faces, cuts, c, left, pF, ft != nullptr ? &leftPos : nullptr);
      const RubberContact::FaceRead frL { ft, leftPos.data() };
      double h = RubberContact::ArrivalTriangles(bs, tris, left, nearestSeg < nT ? nearestSeg : -1, c, v, a, hMax, touchM, arriveM, ft != nullptr ? &frL : nullptr);
      // the reach: a minimum farther than this leaves h as it is (ArrivalStep: d > vc h + a h^2, with vc <= |v| and h only falling from here)
      const double h0 = h, beyondM = v.Length() * h0 + std::max(a, 0.0) * h0 * h0 + arriveM;
      for (const int ci : cuts)
      {
         hits.clear();
         RimExtrema(ci, bs.R, c, hits, nullptr, beyondM + kReachPadM);
         int skip = -1;   // the probe's nearest on this cut: its deepest (the probe takes the first of equals)
         if (nearestSeg == nT + ci)
            for (size_t j = 0; j < hits.size(); j++)
               if (skip < 0 || hits[j].delta > hits[skip].delta)
                  skip = (int)j;
         for (size_t j = 0; j < hits.size(); j++)
         {
            if ((int)j == skip)
               continue;
            const double gap = hits[j].dist - (bs.R + hits[j].rho);
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
         if (m.groove || m.rect)   // (the groove and the rectangle are not drawn yet: a later piece)
            continue;
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
