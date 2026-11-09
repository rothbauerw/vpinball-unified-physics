// license:GPLv3+
// Wire switch collision shape: a chain of capsules (centreline + radius) traced from a wire trigger's drawn mesh.
// Plain math only (no engine types) so offline benches and the engine can both include it. Built once at load; no runtime cost.
// Law (design-wire.md topic 2): the wire is a round tube; its centreline is traced from the mesh, joined into capsules, and
// near-collinear capsules are merged (tolerance mergeTolMm). Scale X sizes the whole wire evenly (1 = 1 in span),
// placed at the author's x, y and rotation with the wire top peakMm above the surface; the author's z is ignored.
#pragma once

#include <vector>
#include <cmath>
#include <cstddef>
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

enum class Shape { WireA, WireB, WireC, WireD };

constexpr double kMmPerVpu = 25.4 * 1.0625 / 50.; // same as VPUTOMM
constexpr double kStdSpanMm = 25.4;               // span at Scale X = 1 (the pictured wire: 1 in)

// Mesh-space tilt of each shape, as Trigger::GenerateMesh poses it (rotation about X before the Z rotation).
inline double TiltDeg(Shape s) { return s == Shape::WireB ? -23. : s == Shape::WireC ? 140. : 0.; }

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

// Centreline of a planar tube mesh, in mesh space (VPU). Positions and outward normals are read with a stride in floats
// (Vertex3D_NoTex2: x,y,z,nx,ny,nz,tu,tv -> stride 8). Method: per vertex, centre = p - R n (n averaged over duplicate
// positions, so seams and flat facets give the radial direction); R from the tube's half width across its plane
// (R = |dx| / |nx|, median); centres projected to the wire plane, clustered into nodes along the tube, ordered end to
// end, then merged where the chain stays within mergeTol of the dropped nodes (Douglas-Peucker).
struct Centreline
{
   std::vector<V3> nodes; // ordered, mesh space
   double r = 0.;         // tube radius, mesh units
};

inline void SimplifyDP(const std::vector<V3>& in, size_t i0, size_t i1, double tol, std::vector<char>& keep)
{
   if (i1 <= i0 + 1)
      return;
   double worst = -1.;
   size_t iw = i0;
   for (size_t i = i0 + 1; i < i1; i++)
   {
      const double d = DistPointSeg(in[i], in[i0], in[i1]);
      if (d > worst) { worst = d; iw = i; }
   }
   if (worst > tol)
   {
      keep[iw] = 1;
      SimplifyDP(in, i0, iw, tol, keep);
      SimplifyDP(in, iw, i1, tol, keep);
   }
}

inline Centreline TraceCentreline(const float* verts, size_t nVerts, size_t stride, double mergeTol)
{
   Centreline out;
   if (nVerts == 0)
      return out;
   // 1. unique positions, normals averaged over duplicates
   std::vector<V3> pos, nrm;
   const double eps2 = 1e-6;
   for (size_t i = 0; i < nVerts; i++)
   {
      const float* v = verts + i * stride;
      const V3 p(v[0], v[1], v[2]), n(v[3], v[4], v[5]);
      size_t j = 0;
      for (; j < pos.size(); j++)
         if ((pos[j] - p).Dot(pos[j] - p) < eps2)
            break;
      if (j == pos.size()) { pos.push_back(p); nrm.push_back(n); }
      else nrm[j] = nrm[j] + n;
   }
   for (V3& n : nrm) { const double l = n.Len(); if (l > 0.) n = n * (1. / l); }
   // 2. wire plane x = x0 (the meshes are flat tubes across x), tube radius from the across-plane extent
   double x0 = 0.;
   for (const V3& p : pos) x0 += p.x;
   x0 /= (double)pos.size();
   std::vector<double> rs;
   for (size_t i = 0; i < pos.size(); i++)
      if (std::fabs(nrm[i].x) > 0.3)
         rs.push_back(std::fabs(pos[i].x - x0) / std::fabs(nrm[i].x));
   if (rs.empty())
      return out;
   std::nth_element(rs.begin(), rs.begin() + rs.size() / 2, rs.end());
   const double R = rs[rs.size() / 2];
   out.r = R;
   // 3. per-vertex ring centre, projected to the wire plane
   std::vector<V3> c(pos.size());
   for (size_t i = 0; i < pos.size(); i++) { c[i] = pos[i] - nrm[i] * R; c[i].x = x0; }
   // 4. cluster centres into nodes (one per ring, rings are at least ~R apart along the tube)
   std::vector<V3> nodes; std::vector<int> cnt;
   const double clus = 0.5 * R;
   for (const V3& p : c)
   {
      size_t j = 0;
      for (; j < nodes.size(); j++)
         if ((nodes[j] * (1. / cnt[j]) - p).Len() < clus)
            break;
      if (j == nodes.size()) { nodes.push_back(p); cnt.push_back(1); }
      else { nodes[j] = nodes[j] + p; cnt[j]++; }
   }
   for (size_t j = 0; j < nodes.size(); j++) nodes[j] = nodes[j] * (1. / cnt[j]);
   // 5. order end to end: start at the node farthest from node 0 (an end), then nearest unvisited each step
   size_t start = 0; double farthest = -1.;
   for (size_t j = 0; j < nodes.size(); j++) { const double d = (nodes[j] - nodes[0]).Len(); if (d > farthest) { farthest = d; start = j; } }
   std::vector<char> used(nodes.size(), 0);
   std::vector<V3> chain; chain.reserve(nodes.size());
   size_t cur = start; used[cur] = 1; chain.push_back(nodes[cur]);
   for (size_t k = 1; k < nodes.size(); k++)
   {
      double best = 1e30; size_t nb = cur;
      for (size_t j = 0; j < nodes.size(); j++)
         if (!used[j]) { const double d = (nodes[j] - nodes[cur]).Len(); if (d < best) { best = d; nb = j; } }
      cur = nb; used[cur] = 1; chain.push_back(nodes[cur]);
   }
   // 6. merge near-collinear capsules
   std::vector<char> keep(chain.size(), 0);
   keep.front() = keep.back() = 1;
   SimplifyDP(chain, 0, chain.size() - 1, mergeTol, keep);
   for (size_t j = 0; j < chain.size(); j++)
      if (keep[j]) out.nodes.push_back(chain[j]);
   return out;
}

// Pose: author's x, y, rotation (deg) and Scale X; surface height at (x, y); all in VPU. wireDmm > 0 overrides the radius.
struct Pose
{
   double x = 0., y = 0., rotDeg = 0., scaleX = 1., surfaceZ = 0.;
   double peakMm = 6.;
   double wireDmm = 0.;
   double thicknessVpu = 0.; // the trigger's WireThickness (added to the radius, as stock adds it to the drawn tube)
};

// Mesh space -> world, uniform scale k and z shift dz (shared by the drawn mesh and the capsules).
struct Transform
{
   double tiltRad = 0., rotRad = 0., k = 1., x = 0., y = 0., dz = 0.;
   V3 Apply(const V3& p) const
   {
      const double ct = std::cos(tiltRad), st = std::sin(tiltRad);
      const V3 q(p.x, ct * p.y - st * p.z, st * p.y + ct * p.z);
      const double cr = std::cos(rotRad), sr = std::sin(rotRad);
      return { (cr * q.x - sr * q.y) * k + x, (sr * q.x + cr * q.y) * k + y, q.z * k + dz };
   }
   V3 ApplyDir(const V3& n) const
   {
      const double ct = std::cos(tiltRad), st = std::sin(tiltRad);
      const V3 q(n.x, ct * n.y - st * n.z, st * n.y + ct * n.z);
      const double cr = std::cos(rotRad), sr = std::sin(rotRad);
      return { cr * q.x - sr * q.y, sr * q.x + cr * q.y, q.z };
   }
};

struct WireBody
{
   std::vector<Capsule> caps; // world, VPU
   Transform xf;              // also poses the drawn mesh's vertices
   double spanMeshVpu = 0.;   // horizontal end-to-end length of the centreline before scaling
};

inline WireBody Build(Shape shape, const float* verts, size_t nVerts, size_t stride, const Pose& pose, double mergeTolMm = 0.05)
{
   WireBody wb;
   Transform& xf = wb.xf;
   xf.tiltRad = TiltDeg(shape) * (3.14159265358979323846 / 180.);
   xf.rotRad = pose.rotDeg * (3.14159265358979323846 / 180.);
   // span: horizontal extent of the tilted centreline (a first trace with a fine tolerance in mesh units)
   const Centreline pre = TraceCentreline(verts, nVerts, stride, 1e-3);
   if (pre.nodes.size() < 2)
      return wb;
   double ymin = 1e30, ymax = -1e30;
   const double ct = std::cos(xf.tiltRad), st = std::sin(xf.tiltRad);
   for (const V3& n : pre.nodes) { const double y = ct * n.y - st * n.z; ymin = std::min(ymin, y); ymax = std::max(ymax, y); }
   wb.spanMeshVpu = ymax - ymin;
   xf.k = pose.scaleX * (kStdSpanMm / kMmPerVpu) / wb.spanMeshVpu;
   // final trace with the merge tolerance mapped back to mesh units
   const Centreline cl = TraceCentreline(verts, nVerts, stride, mergeTolMm / (kMmPerVpu * xf.k));
   xf.x = pose.x; xf.y = pose.y; xf.dz = 0.;
   double topZ = -1e30;
   std::vector<V3> w;
   for (const V3& n : cl.nodes) { w.push_back(xf.Apply(n)); topZ = std::max(topZ, w.back().z); }
   // tube radius: midrange of the drawn vertices' distances to the chain (minimises the worst surface error;
   // hand-modelled tubes are not exactly round)
   for (size_t i = 0; i + 1 < w.size(); i++) { Capsule c; c.a = w[i]; c.b = w[i + 1]; wb.caps.push_back(c); }
   double dmin = 1e30, dmax = 0.;
   for (size_t i = 0; i < nVerts; i++)
   {
      const float* v = verts + i * stride;
      const double d = DistToChain(xf.Apply(V3(v[0], v[1], v[2])), wb.caps);
      dmin = std::min(dmin, d); dmax = std::max(dmax, d);
   }
   wb.caps.clear();
   const double r = pose.wireDmm > 0. ? 0.5 * pose.wireDmm / kMmPerVpu : 0.5 * (dmin + dmax) + pose.thicknessVpu;
   xf.dz = pose.surfaceZ + pose.peakMm / kMmPerVpu - (topZ + r);
   for (size_t i = 0; i + 1 < w.size(); i++)
   {
      Capsule c; c.a = w[i]; c.b = w[i + 1]; c.a.z += xf.dz; c.b.z += xf.dz; c.r = r;
      wb.caps.push_back(c);
   }
   return wb;
}

// The drawn wire (requirements.md section 4: the drawn switch is the collidable switch): each mesh vertex posed by the body's transform, then put on the
// tube around the nearest capsule (radius r, the wire_d_mm the physics uses), its normal the radial direction. Positions and normals with a stride in floats
// (Vertex3D_NoTex2: 8) in and out. Returns the worst |distance to the chain - r| of the result, VPU (round-off).
inline double PoseDrawn(const WireBody& wb, const float* verts, size_t nVerts, size_t stride, float* out, size_t outStride)
{
   double worst = 0.;
   for (size_t i = 0; i < nVerts; i++)
   {
      const float* v = verts + i * stride;
      const V3 p = wb.xf.Apply(V3(v[0], v[1], v[2]));
      const V3 n0 = wb.xf.ApplyDir(V3(v[3], v[4], v[5]));
      double best = 1e30; V3 c;
      for (const Capsule& cp : wb.caps)
      {
         const V3 ab = cp.b - cp.a;
         const double l2 = ab.Dot(ab);
         const double t = l2 > 0. ? std::clamp((p - cp.a).Dot(ab) / l2, 0., 1.) : 0.;
         const V3 q = cp.a + ab * t;
         const double d = (p - q).Len();
         if (d < best) { best = d; c = q; }
      }
      const double r = wb.caps.empty() ? 0. : wb.caps[0].r;
      const V3 n = best > 1e-9 ? (p - c) * (1. / best) : n0;   // (a vertex on the centreline: an end cap's centre keeps its own normal)
      const V3 w = c + n * r;
      float* o = out + i * outStride;
      o[0] = (float)w.x; o[1] = (float)w.y; o[2] = (float)w.z;
      o[3] = (float)n.x; o[4] = (float)n.y; o[5] = (float)n.z;
      worst = std::max(worst, std::fabs(DistToChain(V3(o[0], o[1], o[2]), wb.caps) - r));
   }
   return worst;
}

} // namespace WireShape
