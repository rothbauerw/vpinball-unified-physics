#pragma once
// Playfield cuts: holes cut in the playfield (design runs/2026-11-09-slot-hole/design.md, section 2 third version). Each cut is a [PlayfieldCut:<name>] section of the table's physics ini;
// [Features] playfield_cuts = 0 ignores every one. Two shapes: the plunger lane groove (below) and the slot a wire switch pokes through:
//   [PlayfieldCut:<name>]
//   shape = slot        (or groove, below); any other value: the cut is ignored (on-screen notice and WARN)
//   x, y                centre of the straight part, VPX units
//   angle               direction of the slot, degrees, the trigger's rotation convention: along (-sin angle, cos angle) (wireshape.h Place)
//   length_mm           the straight part, between the round ends' centres
//   width_mm            default and largest 5.03 (0.198 in, PartData 7), smallest 3.175 (1/8 in); outside: the nearest limit, with a notice and a WARN
// Playfield only (no z key): the cut is at the playfield's height at its centre, and the playfield must be flat over the cut's band (else the cut is ignored, with a notice).
//
// The plunger lane groove (paper runs/2026-11-10-plunger-groove/paper.md; physics only, not drawn yet):
//   [PlayfieldCut:<name>]
//   shape = groove      a round-bottomed groove routed along the plunger lane: radius 0.625 in, depth 0.188 in at the plunger end falling linearly to 0 at the far tip (both fixed)
//   x, y                the plunger end, on the groove's axis, VPX units. The end is open: the groove runs on out of the cut's area there, the ball is held by the plunger tip
//   angle               direction from the plunger end toward the tip, degrees, the slot's convention: along (-sin angle, cos angle)
//   length_mm           plunger end to tip; typical and limit 304.8 .. 508 (12 .. 20 in, the range the paper checked); outside: the nearest limit, with a notice and a WARN
// Its contacts are exact: the floor (a tilted cylinder), the rim (a 0.25 mm tube about a level ellipse) and the playfield's own faces outside the rim (playfieldcut.cpp GrooveShape).
//
// A generic rectangle hole through the playfield (brief runs/2026-11-10-rect-cut; physics only, not drawn yet):
//   [PlayfieldCut:<name>]
//   shape = rectangle
//   x, y                its centre, VPX units
//   angle               direction of its length, degrees, the slot's convention: along (-sin angle, cos angle)
//   length_mm, width_mm the side lengths (above 0)
//   corner_mm           the corners' radius, 0 (sharp) .. half the shorter side (default 0); outside: the nearest limit, with a notice and a WARN
// Through the playfield's thickness (12.7 mm). Its contacts are exact, one Hertz law (playfieldcut.cpp RectExtrema): the top rim and the underside rim (0.25 mm tubes about the outline offset
// outward by their radius: 4 lines and 4 quarter circles of radius corner + 0.25 mm, so a sharp corner's tube closes round the corner point), and the 4 vertical walls with (corner > 0)
// 4 concave quarter cylinders between them. A hole narrower than the ball holds it on the rims; a wider one lets it fall through, the walls and underside rim in contact on the way.
//
// Cuts that meet (design runs/2026-11-10-cut-junction sections 1 - 2; found once at load, playfieldcut.cpp Meets): the wood is the slab minus the union of the cuts, and the ball meets
// only its boundary. One rule for every pair: each cut's rims, walls and floor count only where they bound wood (not within 0.25 mm of another cut's removed volume, Gone); where two
// through cuts' outlines cross, a rounded vertical edge; where a groove's floor meets any through cut's wall (straight or curved, or a slot inside the groove), a 0.25 mm tube about
// the wall's centre line lifted onto the floor's offset cylinder; a 0.25 mm sphere wherever three surfaces meet (Join). Two grooves meeting is not described (a notice). A cut that
// meets no other is unchanged.
//
// Shape, in the slot's frame (s along the slot, t across, z up from the playfield top; a = width / 2, u = |s| - length / 2 past each end of the straight part):
//   outline (where the top ends)  straight sides t = +-a; the drilled round end w^2 = a^2 - u^2 up to U1 = L - sqrt(L^2 - a^2); then the drill tip's ramp w^2 = a^2 - U1^2 - 2 U1 (u - U1),
//                                 tangent at U1, to a point at u = L (7.5 mm).
//   through-hole                  the straight part and the round ends: open, vertical walls through the playfield.
//   ramp floor                    drawing only (DrawnSlot::FloorArc): depth on the axis D = 1.75 mm (L - u) / (L - a); in each cross-section u = const a circular arc through the two
//                                 outline points and the floor centre (0, -D). Not a contact: the ball meets only the rim.
//   rim                           rounded, radius rho = 0.25 mm.
// The contact is exact (no faces): the rim is a tube of radius rho round a smooth closed centre line at depth rho, the outline offset outward by rho everywhere (straight part, round ends, ramp,
// tip). The line has no corner: the outline is smooth (the ramp joins the round end tangent at U1) and convex, so its outward offset is too. The playfield's own faces are left out inside the
// centre line (PlayfieldCuts::Probe); the rim gives one contact per true local minimum of the distance from the ball's centre to the centre line.

#include <memory>
#include <string>
#include <vector>
#include "physics/rubbercontact.h"

class PhysicsEngine;
class PinTable;
class ContactSurface;

namespace PlayfieldCuts
{
   // Reads the cuts of the table's physics ini and gives each to the playfield host it lies on. Called when the physics engine is built, after the playfield host and the turntables.
   void Setup(PhysicsEngine *physics, PinTable *table);
   // Cuts are held for this playfield host (false: the host's probe, patch and step bound are today's, unchanged)
   bool Any(const ContactSurface *host);
   // shared_finder P4b ([Features] shared_finder = 1; runs/2026-11-11-finder/design.md section 3): a ball's cut pieces found by its per-cycle query (ContactSurfaceHost::FinderRun): every cut whose
   // box (Rim::BoxDist's) lies within r of c0, in index order. r < 0: none yet. ver: the cuts' Version it was made for (a later Setup makes it stale). NearCuts reads it in place of a scan.
   struct CutQuery
   {
      RubberContact::V3 c0;
      double r = -1.0;
      uint64_t ver = 0;
      std::vector<int> cuts;
   };
   uint64_t Version();   // changes at every Setup (the cuts read anew)
   int Count();          // the cuts held
   void Box(int i, RubberContact::V3 &lo, RubberContact::V3 &hi);   // cut i's world box, m (the box NearCuts tests)
   // P4b counters (VPX_CHECK_PROFILE=1 or VPX_FINDER_CHECK=1; diagnostics only): NearCuts calls, cuts tested, ticks (TSC), re-queries; cut queries and their ticks; check-mode calls, mismatches, misses
   struct CutStats
   {
      uint64_t calls = 0, tested = 0, ticks = 0, requery = 0, queries = 0, queryTicks = 0, checks = 0, mismatch = 0, miss = 0;
   };
   CutStats &Stats();
   // The playfield's faces (tris, candidates faces) with the cuts of this host within reachM of the ball: the faces' nearest point (ProbeTriangles, the turntable holes as there) where it is not
   // inside a cut, or the deepest rim contact. patch: the group step's patch sum (the faces' fold sum, or the rim contacts' symmetric sum, playfieldcut.cpp Probe; zero = the touch's depth
   // zero). A rim contact's seg is tris.size() + the cut's index.
   RubberContact::TubeProbe Probe(const ContactSurface *host, const RubberContact::BallSI &bs, const std::vector<RubberContact::SurfaceTri> &tris, const std::vector<int> &faces, double edgeRho,
                                  const RubberContact::V3 &c, double reachM, const RubberContact::CarvedHole *holes, int nHoles, bool patch, double zero, RubberContact::FaceTable *ft = nullptr, CutQuery *cq = nullptr);   // ft: shared_finder P3a, the ball's measurement table; cq: P4b, the ball's cut pieces (null: key off, the scan)
   // The step bound (ContactSurface::ArrivalS) of the same surface: the faces not inside a cut (ArrivalTriangles), and every rim contact but the probe's nearest (ArrivalStep each)
   double ArrivalS(const ContactSurface *host, const RubberContact::BallSI &bs, const std::vector<RubberContact::SurfaceTri> &tris, const std::vector<int> &faces, int nearestSeg,
                   const RubberContact::V3 &c, const RubberContact::V3 &v, double a, double hMax, double touchM, double arriveM, double reachM, RubberContact::FaceTable *ft = nullptr, CutQuery *cq = nullptr);
   // A cut as the drawing needs it (parts/playfieldcutdraw.cpp, design section 8): the same numbers as its rim. Frame: centre (cx, cy) at the playfield's top z0, along (lx, ly), across
   // (px, py), metres (unit: metres per VPX unit). Shape (the slot's frame, comment above): a, hs, L, U1, wJ; the ramp's outline u = L - tau^2, w = c tau (tau from 0 at the tip to tauJ at
   // U1); th1 the round end's angle at U1; thick the through-hole's depth (the physics' through-hole, Rim::Holds).
   struct DrawnSlot
   {
      std::string name;
      double unit = 1.0, cx = 0.0, cy = 0.0, z0 = 0.0, lx = 0.0, ly = 1.0, px = -1.0, py = 0.0;
      double a = 0.0, hs = 0.0, L = 0.0, U1 = 0.0, wJ = 0.0, c = 0.0, tauJ = 0.0, th1 = 0.0, thick = 0.0;
      double W(double u) const;   // the outline's half-width at u past the straight part's end (0 past the tip)
      double C(double u) const;   // the through-hole's half-width (0 past the round end, u >= a)
      // the ramp floor's cross-section at u (U1 < u < L): the circular arc t^2 + (z - zc)^2 = rho^2 through (+-W, 0) and (0, -depth), z up from the top
      void FloorArc(double u, double &rho, double &zc) const;
   };
   // The cuts read for this table (none: no cut active, or the cuts are another table's)
   std::vector<DrawnSlot> Drawn(const PinTable *table);
}
