// license:GPLv3+
// Wire switch shapes A - E as the real parts (PartData S1, S7; docs/unified-physics/wire-switches/wire_A..E_installed.md; GroupStepDesign 23.4).
// Plain math only (no engine types) so offline benches and the engine can both include it. Built once at load; no runtime cost.
// Each shape is a fixed centreline in mm in the switch frame: x along the wire's plane from the lever's pivot away from the switch, z the height above the playfield
// top at rest with the switch at its agreed mount (pivot 19.3 mm below the playfield top). D's arch is a 14.7 mm arc; every other bend a 1.3 mm centreline-radius
// arc tangent to both straight runs (runs/2026-11-09-wire-sizes/radii/radii.md); arcs sampled every 5 deg. The lever end, the peak and the knee are named points
// (the wire files' corners), not searched for. No mesh tracing: the drawn tube is built from the same centreline (BuildTube), so drawn = touched.
#pragma once

#include <vector>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <algorithm>

namespace WireShape
{

struct V3
{
   double x = 0., y = 0., z = 0.;
   V3() = default;
   V3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) { }
   V3 operator+(const V3& o) const { return { x + o.x, y + o.y, z + o.z }; }
   V3 operator-(const V3& o) const { return { x - o.x, y - o.y, z - o.z }; }
   V3 operator*(double s) const { return { x * s, y * s, z * s }; }
   double Dot(const V3& o) const { return x * o.x + y * o.y + z * o.z; }
   double Len() const { return std::sqrt(Dot(*this)); }
};

struct Capsule
{
   V3 a, b;      // centreline end points
   double r = 0.; // tube radius
};

enum class Shape { WireA, WireB, WireC, WireD, WireE };

constexpr double kMmPerVpu = 25.4 * 1.0625 / 50.; // same as VPUTOMM
constexpr double kPivotZMm = -19.3;                // the lever's pivot below the playfield top at the agreed mount (S7: bracket flush under the playfield)
constexpr double kBendRMm = 1.3;                   // centreline radius of every bend but D's arch (measured 1.2 - 1.4 mm on every bend of 50 deg or more)
constexpr double kArchRMm = 14.7;                  // D's arch radius (user 2026-11-09)
constexpr double kArcStepDeg = 5.;                 // arcs sampled at most this many degrees per capsule
constexpr double kMaxWireDMm = 2. * kBendRMm;      // a wire this thick or thicker folds on the inside of a 1.3 mm bend
constexpr double kPi = 3.14159265358979323846;

inline double DistPointSeg(const V3& p, const V3& a, const V3& b)
{
   const V3 ab = b - a;
   const double l2 = ab.Dot(ab);
   double t = l2 > 0. ? (p - a).Dot(ab) / l2 : 0.;
   t = std::clamp(t, 0., 1.);
   return (p - (a + ab * t)).Len();
}

// Distance from p to the chain's centreline (minus r gives the distance to the tube surface).
inline double DistToChain(const V3& p, const std::vector<Capsule>& caps)
{
   double d = 1e30;
   for (const Capsule& c : caps)
      d = std::min(d, DistPointSeg(p, c.a, c.b));
   return d;
}

// A point in the switch frame, mm: x from the pivot away from the switch, z above the playfield top (heightoffset 0)
struct P2
{
   double x = 0., z = 0.;
   P2 operator+(const P2& o) const { return { x + o.x, z + o.z }; }
   P2 operator-(const P2& o) const { return { x - o.x, z - o.z }; }
   P2 operator*(double s) const { return { x * s, z * s }; }
   double Len() const { return std::sqrt(x * x + z * z); }
};

struct Centreline
{
   std::vector<P2> pts;          // sampled, from the wire's start on the lever to its free end
   P2 pivot { 0., kPivotZMm };
   P2 leverEnd, peak, knee;      // named points (the wire file's corners); knee only where the file names one (A, B)
   bool hasKnee = false;
   double spanMid = 0.;          // x of the middle of the part above the playfield: the trigger's x, y
};

// Arc about c of radius r from angle a0 to a1 (rad), at most kArcStepDeg per piece; the first point is skipped when out already ends there
inline void AppendArc(std::vector<P2>& out, const P2& c, double r, double a0, double a1)
{
   const int n = std::max(1, (int)std::ceil(std::fabs(a1 - a0) / (kArcStepDeg * kPi / 180.) - 1e-9));
   for (int k = 0; k <= n; k++)
   {
      const double a = a0 + (a1 - a0) * k / n;
      const P2 p { c.x + r * std::cos(a), c.z + r * std::sin(a) };
      if (k == 0 && !out.empty() && (out.back() - p).Len() < 1e-6)
         continue;
      out.push_back(p);
   }
}

// Straight runs through c[0] .. c.back(), each inner corner rounded by an arc of centreline radius r tangent to both runs
inline void AppendFilleted(std::vector<P2>& out, const std::vector<P2>& c, double r)
{
   if (out.empty() || (out.back() - c[0]).Len() > 1e-6)
      out.push_back(c[0]);
   for (size_t i = 1; i + 1 < c.size(); i++)
   {
      const P2 in = (c[i] - c[i - 1]) * (1. / (c[i] - c[i - 1]).Len()), ou = (c[i + 1] - c[i]) * (1. / (c[i + 1] - c[i]).Len());
      const double cr = in.x * ou.z - in.z * ou.x, turn = std::atan2(std::fabs(cr), in.x * ou.x + in.z * ou.z);
      if (turn < 1e-9)
      {
         out.push_back(c[i]);
         continue;
      }
      const double side = cr > 0. ? 1. : -1.;   // a left (counter-clockwise) turn has its centre on the left of the run in
      const P2 a = c[i] - in * (r * std::tan(0.5 * turn));
      const P2 ctr = a + P2 { -in.z * side, in.x * side } * r;
      const double a0 = std::atan2(a.z - ctr.z, a.x - ctr.x);
      AppendArc(out, ctr, r, a0, a0 + side * turn);
   }
   out.push_back(c.back());
}

// p turned up (counter-clockwise, the wire rising) by deg about the pivot
inline P2 RaiseLever(const P2& p, double deg)
{
   const double a = deg * (kPi / 180.), c = std::cos(a), sn = std::sin(a);
   const P2 r = p - P2 { 0., kPivotZMm };
   return { r.x * c - r.z * sn, r.x * sn + r.z * c + kPivotZMm };
}

// The agreed rest pose of each shape (the wire files' points; heights from the playfield top at heightoffset 0). A, B and C: the lever raised 2026-11-09 (user) so the
// rounded centreline's top is 6.0 mm (A, B) and C's free tip 8.0 mm; the files' points as photographed are turned up by that angle about the pivot
inline Centreline MakeCentreline(Shape s)
{
   Centreline cl;
   const auto up = [](std::vector<P2> v, double deg) { for (P2& p : v) p = RaiseLever(p, deg); return v; };
   switch (s)
   {
   case Shape::WireA:   // Cherry DB5: foot on the lever, up leg 22.5, span 22.8, down leg 12.3 mm
   {
      constexpr double raiseDeg = 1.90;   // rounded top 6.0 mm
      const std::vector<P2> c = up({ { 7.4, -17.6 }, { 22.3, -14.3 }, { 33.4, 5.3 }, { 53.8, -4.8 }, { 51.0, -16.7 } }, raiseDeg);
      cl.leverEnd = c[1]; cl.peak = c[2]; cl.knee = c[3]; cl.hasKnee = true;
      AppendFilleted(cl.pts, c, kBendRMm);
   }
      break;
   case Shape::WireB:   // V wire on a basic switch: lever run 19.3, rise 27.4, fall 30.8, tail 12.0 mm
   {
      constexpr double raiseDeg = 0.97;   // rounded top 6.0 mm
      const std::vector<P2> c = up({ { 1.6, -18.0 }, { 20.5, -14.3 }, { 39.0, 5.9 }, { 61.3, -15.4 }, { 72.9, -18.6 } }, raiseDeg);
      cl.leverEnd = c[1]; cl.peak = c[2]; cl.knee = c[3]; cl.hasKnee = true;
      AppendFilleted(cl.pts, c, kBendRMm);
   }
      break;
   case Shape::WireC:   // one-way: lever run, a leaning rise, a diagonal up to the tip (the peak)
   {
      constexpr double raiseDeg = 2.06;   // free tip 8.0 mm
      const std::vector<P2> c = up({ { 2.7, -18.0 }, { 27.2, -14.9 }, { 32.0, -1.5 }, { 56.2, 6.0 } }, raiseDeg);
      cl.leverEnd = c[1]; cl.peak = c[3];
      AppendFilleted(cl.pts, c, kBendRMm);
   }
      break;
   case Shape::WireD:   // the arch: file x runs from the arch top toward the pivot (41.2), here x = 41.2 - file x
   {
      const double xp = 41.2;
      cl.leverEnd = { xp - 19.1, -14.3 }; cl.peak = { xp, 5.8 };
      const P2 ctr { xp, -8.9 }, side { xp + 14.7, -8.9 }, corner { xp + 14.7, -14.8 }, tip { xp + 22.5, -14.8 };
      // the wire runs 15.9 mm along the lever from its end toward the pivot
      const P2 toPivot = (cl.pivot - cl.leverEnd) * (1. / (cl.pivot - cl.leverEnd).Len());
      const P2 start = cl.leverEnd + toPivot * 15.9;
      // the up leg leaves the arch on the tangent through the lever end
      const P2 dl = cl.leverEnd - ctr;
      const double aT = std::atan2(dl.z, dl.x) - std::acos(kArchRMm / dl.Len());
      const P2 tan { ctr.x + kArchRMm * std::cos(aT), ctr.z + kArchRMm * std::sin(aT) };
      AppendFilleted(cl.pts, { start, cl.leverEnd, tan }, kBendRMm);
      AppendArc(cl.pts, ctr, kArchRMm, aT + (aT < 0. ? 2. * kPi : 0.), 0.);
      AppendFilleted(cl.pts, { side, corner, tip }, kBendRMm);
      break;
   }
   case Shape::WireE:   // Cherry DB: along the lever 17.4, up leg 15.3, diagonal 17.9, down leg 30.8 mm
      cl.leverEnd = { 25.8, -14.3 }; cl.peak = { 45.7, 9.4 };
      AppendFilleted(cl.pts, { { 8.5, -15.7 }, cl.leverEnd, { 30.2, 0.4 }, cl.peak, { 45.5, -21.5 } }, kBendRMm);
      break;
   }
   // the part above the playfield: from the first to the last crossing of z = 0
   double x0 = 1e30, x1 = -1e30;
   for (size_t i = 0; i < cl.pts.size(); i++)
   {
      const P2 &p = cl.pts[i];
      if (p.z >= 0.) { x0 = std::min(x0, p.x); x1 = std::max(x1, p.x); }
      if (i > 0 && (cl.pts[i - 1].z < 0.) != (p.z < 0.))
      {
         const P2 &q = cl.pts[i - 1];
         const double xc = q.x + (p.x - q.x) * (0. - q.z) / (p.z - q.z);
         x0 = std::min(x0, xc); x1 = std::max(x1, xc);
      }
   }
   cl.spanMid = x1 >= x0 ? 0.5 * (x0 + x1) : 0.;
   return cl;
}

// The table: the trigger's x, y (VPU) at the middle of the part above the playfield; rotation 0 puts the wire's plane along y with the switch toward the top of the
// table (-y) and the wire toward the player (+y), turned as VPX turns a mesh by the trigger's rotation; the playfield top = the trigger's Surface (VPU); heightOffsetMm
// raises (+) or lowers the whole switch from the agreed mount
struct Placement
{
   double x = 0., y = 0., rotDeg = 0., surfaceZ = 0., heightOffsetMm = 0.;
};

struct WireBody
{
   std::vector<Capsule> caps;   // world, VPU
   std::vector<V3> nodes;       // the sampled centreline, world, VPU
   V3 pivot, lane;              // the lever's pivot; lane: unit, horizontal, from the switch toward the wire
   V3 leverEnd, peak, knee;     // the named points, world, VPU
   double r = 0.;               // tube radius, VPU
};

inline V3 PlaceMm(const Centreline& cl, const Placement& pl, const P2& p)
{
   const double rot = pl.rotDeg * (kPi / 180.);
   const V3 lane(-std::sin(rot), std::cos(rot), 0.);   // +y turned as Trigger::GenerateMesh turns a mesh: (x c - y s, x s + y c)
   const double s = (p.x - cl.spanMid) / kMmPerVpu;
   return { pl.x + lane.x * s, pl.y + lane.y * s, pl.surfaceZ + (p.z + pl.heightOffsetMm) / kMmPerVpu };
}

inline WireBody Place(const Centreline& cl, const Placement& pl, double wireDMm)
{
   WireBody wb;
   const double rot = pl.rotDeg * (kPi / 180.);
   wb.lane = V3(-std::sin(rot), std::cos(rot), 0.);
   wb.r = 0.5 * wireDMm / kMmPerVpu;
   for (const P2& p : cl.pts)
      wb.nodes.push_back(PlaceMm(cl, pl, p));
   for (size_t i = 0; i + 1 < wb.nodes.size(); i++)
   {
      Capsule c; c.a = wb.nodes[i]; c.b = wb.nodes[i + 1]; c.r = wb.r;
      wb.caps.push_back(c);
   }
   wb.pivot = PlaceMm(cl, pl, cl.pivot);
   wb.leverEnd = PlaceMm(cl, pl, cl.leverEnd);
   wb.peak = PlaceMm(cl, pl, cl.peak);
   wb.knee = PlaceMm(cl, pl, cl.knee);
   return wb;
}

// The drawn wire: a tube of radius r around the nodes (a planar chain whose plane holds the unit horizontal lane), `sides` vertices per ring, a half sphere at each
// end, so every vertex lies on the touched capsules (within r (1 - cos 2.5 deg) on the inside of a 5 deg sample). Vertex layout x, y, z, nx, ny, nz, tu, tv.
inline void BuildTube(const std::vector<V3>& nodes, const V3& lane, double r, int sides, std::vector<float>& verts, std::vector<uint16_t>& idx)
{
   verts.clear();
   idx.clear();
   const size_t n = nodes.size();
   if (n < 2)
      return;
   const V3 e(lane.y, -lane.x, 0.);   // across the wire's plane
   const auto unit = [](const V3& v) { const double l = v.Len(); return l > 0. ? v * (1. / l) : v; };
   const auto cross = [](const V3& a, const V3& b) { return V3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); };
   double total = 0.;
   for (size_t i = 0; i + 1 < n; i++)
      total += (nodes[i + 1] - nodes[i]).Len();
   const auto ring = [&](const V3& c, const V3& t, double rr, double off, double tu) {
      const V3 w = unit(cross(t, e));
      for (int k = 0; k < sides; k++)
      {
         const double ph = 2. * kPi * k / sides;
         const V3 radial = e * std::cos(ph) + w * std::sin(ph);
         const V3 p = c + radial * rr + t * off;
         const V3 nn = unit(p - c);
         const float v[8] = { (float)p.x, (float)p.y, (float)p.z, (float)nn.x, (float)nn.y, (float)nn.z, (float)tu, (float)(ph / (2. * kPi)) };
         verts.insert(verts.end(), v, v + 8);
      }
   };
   const auto join = [&](int r0, int r1) {
      for (int k = 0; k < sides; k++)
      {
         const int k1 = (k + 1) % sides;
         const uint16_t a = (uint16_t)(r0 * sides + k), b = (uint16_t)(r0 * sides + k1), c = (uint16_t)(r1 * sides + k), d = (uint16_t)(r1 * sides + k1);
         idx.insert(idx.end(), { a, c, b, b, c, d });
      }
   };
   int rings = 0;
   // start cap: rings at 60 and 30 deg from the end's axis, then the body rings
   const V3 t0 = unit(nodes[1] - nodes[0]), t1 = unit(nodes[n - 1] - nodes[n - 2]);
   for (const double a : { 60., 30. })
      ring(nodes[0], t0, r * std::sin(a * kPi / 180.), -r * std::cos(a * kPi / 180.), 0.), rings++;
   double s = 0.;
   for (size_t i = 0; i < n; i++)
   {
      if (i > 0)
         s += (nodes[i] - nodes[i - 1]).Len();
      const V3 t = i == 0 ? t0 : i + 1 == n ? t1 : unit(unit(nodes[i] - nodes[i - 1]) + unit(nodes[i + 1] - nodes[i]));
      ring(nodes[i], t, r, 0., total > 0. ? s / total : 0.);
      rings++;
   }
   for (const double a : { 30., 60. })
      ring(nodes[n - 1], t1, r * std::sin(a * kPi / 180.), r * std::cos(a * kPi / 180.), 1.), rings++;
   for (int i = 0; i + 1 < rings; i++)
      join(i, i + 1);
   // the poles
   const auto pole = [&](const V3& p, const V3& nn, int rg, bool first) {
      const uint16_t pi = (uint16_t)(verts.size() / 8);
      const float v[8] = { (float)p.x, (float)p.y, (float)p.z, (float)nn.x, (float)nn.y, (float)nn.z, first ? 0.f : 1.f, 0.f };
      verts.insert(verts.end(), v, v + 8);
      for (int k = 0; k < sides; k++)
      {
         const uint16_t a = (uint16_t)(rg * sides + k), b = (uint16_t)(rg * sides + (k + 1) % sides);
         if (first) idx.insert(idx.end(), { pi, a, b });
         else idx.insert(idx.end(), { pi, b, a });
      }
   };
   pole(nodes[0] - t0 * r, t0 * -1., 0, true);
   pole(nodes[n - 1] + t1 * r, t1, rings - 1, false);
}

} // namespace WireShape
