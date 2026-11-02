// license:GPLv3+
//
// What a part on one horizontal axis is made of (a gate, a spinner; docs/GateSpinnerDesign.txt section 3): its mass properties from its shape, its switch blade. The moving is the group's (hingegroup.h, bodypart.h). Pure SI maths with no VPX types.
//
// Frame: the axis passes through A along E (unit, horizontal). N is the direction a ball moves when it passes from the front, Z = E x N is up, E = N x Z.
// theta = 0: the body as built (a gate hangs straight down); a positive theta swings the lower part toward N (a rotation about +E). A point of the body is given in
// body coordinates (e along E, d "down" at theta = 0, f "front" at theta = 0) and sits at A + E e + Down(theta) d + Front(theta) f,
// Down(theta) = -Z cos theta + N sin theta, Front(theta) = N cos theta + Z sin theta.

#pragma once

#include <cmath>
#include <algorithm>
#include <vector>

#include "physics/rubbercontact.h"

namespace HingeBody
{
   using RubberContact::V3;

   // Mass properties of the body about its axis, built from its shape (body coordinates, metres: x = e along the axis, y = d down, z = f front)
   struct MassProps
   {
      double m = 0.0;     // kg
      V3 mc;              // sum of m * position (first moment), kg m
      double I = 0.0;     // about the axis, kg m^2
      V3 Cm() const { return m > 0.0 ? mc / m : V3(); }

      // A straight wire from p to q, linear density lambda (kg/m). Exact for a thin rod: the squared distance to the axis is quadratic along it.
      void AddWire(const V3 &p, const V3 &q, double lambda)
      {
         const double L = (q - p).Length();
         const double dm = lambda * L;
         m += dm;
         mc = mc + (p + q) * (0.5 * dm);
         // r_perp^2 = d^2 + f^2: integral over the segment of (a + t (b - a))^2 = L (a^2 + a b + b^2) / 3 per component
         auto sq = [](double a, double b) { return (a * a + a * b + b * b) / 3.0; };
         I += dm * (sq(p.y, q.y) + sq(p.z, q.z));
      }
      // A flat triangle (a, b, c), areal density sigma (kg/m^2). Exact for a thin plate (second moment of a triangle).
      void AddPlate(const V3 &a, const V3 &b, const V3 &c, double sigma)
      {
         const double A = 0.5 * (b - a).Cross(c - a).Length();
         const double dm = sigma * A;
         m += dm;
         mc = mc + (a + b + c) * (dm / 3.0);
         auto sq = [](double x, double y, double z) { return (x * x + y * y + z * z + x * y + y * z + z * x) / 6.0; };
         I += dm * (sq(a.y, b.y, c.y) + sq(a.z, b.z, c.z));   // integral of y^2 over a triangle = A/6 (sum of squares + cross terms of the vertices)
      }
      void AddQuad(const V3 &a, const V3 &b, const V3 &c, const V3 &d, double sigma)
      {
         AddPlate(a, b, c, sigma);
         AddPlate(a, c, d, sigma);
      }
   };

   // A crank on the shaft pressing a switch blade (a spinner's leaf or micro switch): the crank of radius r points along the angle phase; the blade sits gap from
   // the shaft, so it is pushed in by delta = r cos(theta - phase) - gap while that is positive. The blade pushes back with F = F0 + k delta + a damper c on the part of the move inside its zone,
   // never pulling. Stored energy: the blade's spring takes energy while it is pressed in and gives it back as it is released: the speed-up / slow-down of a spinner over its switch.
   // A preload F0 comes on over the blade's first rampM of travel (the lever's own give): a force that jumps with the position has no root at a slow step (the bat's damper fault), and a real lever has give.
   struct Cam
   {
      double phase = 0.0;   // rad
      double r = 3e-3;      // crank radius, m
      double gap = 0.5e-3;  // blade distance from the shaft, m
      double k = 0.0;       // N/m
      double F0 = 0.0;      // N (preload while pressed)
      double rampM = 1e-4;  // m (see above); UNVERIFIED, only with a preload
      double c = 0.0;       // N s/m
      double Delta(const double th) const { return r * std::cos(th - phase) - gap; }
      double V(const double d) const
      {
         if (d <= 0.0)
            return 0.0;
         const double vk = 0.5 * k * d * d;
         if (F0 <= 0.0)
            return vk;
         return vk + (d < rampM ? 0.5 * F0 * d * d / rampM : F0 * (d - 0.5 * rampM));
      }
      double F(const double d) const
      {
         if (d <= 0.0)
            return 0.0;
         return k * d + (F0 <= 0.0 ? 0.0 : d < rampM ? F0 * d / rampM : F0);
      }
   };
}
