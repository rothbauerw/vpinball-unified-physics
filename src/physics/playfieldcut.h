#pragma once
// Playfield cuts: holes cut in the playfield (design runs/2026-11-09-slot-hole/design.md, section 2 third version). Each cut is a [PlayfieldCut:<name>] section of the table's physics ini;
// [Features] playfield_cuts = 0 ignores every one. The only shape so far is the slot a wire switch pokes through:
//   [PlayfieldCut:<name>]
//   shape = slot        the only shape; any other value: the cut is ignored (on-screen notice and WARN)
//   x, y                centre of the straight part, VPX units
//   angle               direction of the slot, degrees, the trigger's rotation convention: along (-sin angle, cos angle) (wireshape.h Place)
//   length_mm           the straight part, between the round ends' centres
//   width_mm            default and largest 5.03 (0.198 in, PartData 7), smallest 3.175 (1/8 in); outside: the nearest limit, with a notice and a WARN
// Playfield only (no z key): the cut is at the playfield's height at its centre, and the playfield must be flat over the cut's band (else the cut is ignored, with a notice).
//
// Shape, in the slot's frame (s along the slot, t across, z up from the playfield top; a = width / 2, u = |s| - length / 2 past each end of the straight part):
//   outline (where the top ends)  straight sides t = +-a; the drilled round end w^2 = a^2 - u^2 up to U1 = L - sqrt(L^2 - a^2); then the drill tip's ramp w^2 = a^2 - U1^2 - 2 U1 (u - U1),
//                                 tangent at U1, to a point at u = L (7.5 mm).
//   through-hole                  the straight part and the round ends: open, vertical walls through the playfield.
//   ramp floor                    drawing only (PlayfieldCutMesh::BuildSlot): depth on the axis D = 1.75 mm (L - u) / (L - a); in each cross-section u = const a circular arc through the two
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

// The mesh builder, kept apart from the physics host for the drawing (design section 8); the physics does not use it.
namespace PlayfieldCutMesh
{
   struct Slot   // world frame, metres
   {
      double cx = 0.0, cy = 0.0, z = 0.0;   // centre of the straight part, at the playfield's top
      double lx = 0.0, ly = 1.0;            // along the slot (unit)
      double halfStraight = 0.0;            // half the straight part
      double halfWidth = 0.0;               // a
   };
   struct Stats
   {
      int faces = 0;                        // faces of the band
      int stations = 0;                     // cross-sections along the slot
      double sBand = 0.0, tBand = 0.0;      // the band: |s| <= sBand, |t| <= tBand (m)
      double devRim = 0.0, devFloor = 0.0;  // largest distance of the true shape from the faces (m): top, rims and outline / ramp floor (midway between cross-sections)
   };
   // The faces of the band around the slot (world frame, metres; a, b, c per face)
   void BuildSlot(const Slot &s, std::vector<RubberContact::SurfaceTri> &out, Stats &st);
}

namespace PlayfieldCuts
{
   // Reads the cuts of the table's physics ini and gives each to the playfield host it lies on. Called when the physics engine is built, after the playfield host and the turntables.
   void Setup(PhysicsEngine *physics, PinTable *table);
   // Cuts are held for this playfield host (false: the host's probe, patch and step bound are today's, unchanged)
   bool Any(const ContactSurface *host);
   // The playfield's faces (tris, candidates faces) with the cuts of this host within reachM of the ball: the faces' nearest point (ProbeTriangles, the turntable holes as there) where it is not
   // inside a cut, or the deepest rim contact. patch: the group step's patch sum (the faces' fold sum, or the rim contacts' symmetric sum, playfieldcut.cpp Probe; zero = the touch's depth
   // zero). A rim contact's seg is tris.size() + the cut's index.
   RubberContact::TubeProbe Probe(const ContactSurface *host, const RubberContact::BallSI &bs, const std::vector<RubberContact::SurfaceTri> &tris, const std::vector<int> &faces, double edgeRho,
                                  const RubberContact::V3 &c, double reachM, const RubberContact::CarvedHole *holes, int nHoles, bool patch, double zero);
   // The step bound (ContactSurface::ArrivalS) of the same surface: the faces not inside a cut (ArrivalTriangles), and every rim contact but the probe's nearest (ArrivalStep each)
   double ArrivalS(const ContactSurface *host, const RubberContact::BallSI &bs, const std::vector<RubberContact::SurfaceTri> &tris, const std::vector<int> &faces, int nearestSeg,
                   const RubberContact::V3 &c, const RubberContact::V3 &v, double a, double hMax, double touchM, double arriveM, double reachM);
}
