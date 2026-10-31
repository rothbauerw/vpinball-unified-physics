// license:GPLv3+

#include "core/stdafx.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>

#include "imgui/imgui.h"
#include "parts/ball.h"
#include "physics/hitrubber.h"
#include "physics/cordband.h"
#include "physics/turntable.h"
#include "physics/contactsystem.h"
#include "physics/hitball.h"
#include "physics/magnet.h"
#include "physics/hitflipper.h"
#include "physics/physicsini.h"
#include "ui/live/LiveUI.h"
#include "physics/cabinet/NudgeHandler.h"
#include "physics/PhysicsEngine.h"
#include "parts/rubber.h"
#include "parts/primitive.h"

bool ContactTraceDetail()   // VPX_CONTACT_TRACE=2: also the per-cycle lines (ContactTrace, ContactPair): huge logs, the log file rotates; 1 = engage / release / part hit lines only
{
   static const bool on = []() { const char *const e = ContactSystem::EnvVar("VPX_CONTACT_TRACE", "VPX_RING_TRACE"); return e != nullptr && atoi(e) >= 2; }();
   return on;
}

bool ContactTraceOn()   // read once: getenv in a per-cycle path costs microseconds
{
   static const bool on = ContactSystem::EnvVar("VPX_CONTACT_TRACE", "VPX_RING_TRACE") != nullptr;
   return on;
}

namespace
{
   constexpr double U = ContactSystem::kMetresPerUnit;
   constexpr double VS = ContactSystem::kVelocityScale;
   constexpr double kMassUnitKg = ContactSystem::kMassUnitKg;
   using PhysicsIni::Trim;

   // Earliest time t >= 0 at which a point at p0 moving with velocity v is at distance r from the point q (-1 = none)
   double SphereTime(const RubberContact::V3 &p0, const RubberContact::V3 &v, const RubberContact::V3 &q, double r)
   {
      const RubberContact::V3 w = p0 - q;
      const double a = v.Dot(v);
      if (a < 1e-18)
         return -1.0;
      const double b = w.Dot(v);
      const double c = w.Dot(w) - r * r;
      const double disc = b * b - a * c;
      if (disc < 0.0)
         return -1.0;
      const double t = (-b - std::sqrt(disc)) / a;
      return t >= 0.0 ? t : -1.0;
   }
}

bool HitRubber::SurfaceHeightAt(const double xU, const double yU, double &zU) const
{
   if (m_surface == nullptr)
      return false;
   const double x = xU * U, y = yU * U;
   bool found = false;
   double best = -1e30;
   for (const RubberContact::SurfaceTri &t : m_surface->tris)
   {
      if (std::abs(t.n.z) < 0.5 || x < t.lo.x || x > t.hi.x || y < t.lo.y || y > t.hi.y)
         continue;
      // point in triangle (xy), then the height of its plane there
      const double d = (t.b.y - t.c.y) * (t.a.x - t.c.x) + (t.c.x - t.b.x) * (t.a.y - t.c.y);
      if (std::abs(d) < 1e-18)
         continue;
      const double l1 = ((t.b.y - t.c.y) * (x - t.c.x) + (t.c.x - t.b.x) * (y - t.c.y)) / d;
      const double l2 = ((t.c.y - t.a.y) * (x - t.c.x) + (t.a.x - t.c.x) * (y - t.c.y)) / d;
      const double l3 = 1.0 - l1 - l2;
      if (l1 < -1e-9 || l2 < -1e-9 || l3 < -1e-9)
         continue;
      const double z = l1 * t.a.z + l2 * t.b.z + l3 * t.c.z;
      if (z > best)
      {
         best = z;
         found = true;
      }
   }
   zU = best / U;
   return found;
}

HitRubber::HitRubber(IEditable *editable, const std::vector<Vertex3Ds> &pts, const float radius, const RubberObjectSettings &settings, const std::filesystem::path &iniPath, const std::string &name)
   : ContactHost(editable, settings, iniPath, name, ""), m_pts(pts), m_radius(radius)
{
   if (!m_pts.empty())
      m_pts.push_back(m_pts.front()); // close the loop
   for (const Vertex3Ds &p : m_pts)
      m_ptsSI.emplace_back((double)p.x * U, (double)p.y * U, (double)p.z * U);
   m_radiusSI = (double)m_radius * U;
   RubberContact::BuildTubePre(m_tubePre, m_ptsSI.data(), (int)m_ptsSI.size());
   m_cumSI.assign(m_ptsSI.size(), 0.0);
   for (size_t i = 1; i < m_ptsSI.size(); i++)
      m_cumSI[i] = m_cumSI[i - 1] + (m_ptsSI[i] - m_ptsSI[i - 1]).Length();
   if (m_set.band)
   {
      DetectBand();
      // A band's law needs a straight span. A rubber without one wraps its post all the way round: it is a ring, and a ball on it takes the ring's law (ModelStandards 2d), never stock.
      if (m_bandSpans.empty())
      {
         const std::string msg = "Rubber " + m_name + ": the physics ini makes it a band or sling, but it has no straight span. Running it as a ring: check the [Sling:] or [Rubber:] section";
         PLOGW << msg;
         if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
            g_pplayer->m_liveUI->PushNotification(msg, 5000);
         m_set.band = false;
         m_set.slingPlate.clear();
      }
   }
}

void RubberSurface::Add(const Vertex3Ds &a, const Vertex3Ds &b, const Vertex3Ds &c)
{
   using RubberContact::V3;
   const Vertex3Ds *p[3] = { &a, &b, &c };
   if (tris.empty())
      for (int k = 0; k < 3; k++)
         lo[k] = hi[k] = (&a.x)[k];
   for (const Vertex3Ds *v : p)
      for (int k = 0; k < 3; k++)
      {
         lo[k] = std::min(lo[k], (&v->x)[k]);
         hi[k] = std::max(hi[k], (&v->x)[k]);
      }
   tris.push_back(RubberContact::MakeSurfaceTri(V3((double)a.x * U, (double)a.y * U, (double)a.z * U), V3((double)b.x * U, (double)b.y * U, (double)b.z * U),
                                                V3((double)c.x * U, (double)c.y * U, (double)c.z * U)));
}

HitRubber::HitRubber(IEditable *editable, std::shared_ptr<RubberSurface> surface, const RubberObjectSettings &settings, const std::filesystem::path &iniPath, const std::string &name,
                     const char *kind)
   : ContactHost(editable, settings, iniPath, name, kind), m_surface(std::move(surface)), m_radius(0.f)
{
}

HitRubber::~HitRubber()
{
   // A ball can outlive its physics engine (a script still references it when the table closes): it must not keep a pointer to a band that is gone (the contact system lets
   // go of the host's balls in ~ContactHost)
   for (const BandBall &e : m_engaged)
      if (e.ball != nullptr && e.ball->m_bandOwner == this)
         e.ball->m_bandOwner = nullptr;
}

static std::vector<RubberContact::SurfaceHole> g_playfieldHoles;   // metres
bool PlayfieldHolesExist() { return !g_playfieldHoles.empty(); }

void HitRubber::AddPlayfieldHole(const float x, const float y, const float r)
{
   g_playfieldHoles.push_back({ (double)x * U, (double)y * U, (double)r * U });
}

void HitRubber::ClearPlayfieldHoles()
{
   g_playfieldHoles.clear();
}

// A ball whose centre is inside a non-legacy kicker's circle belongs to the kicker: VPX implements the kicker's funnel, capture and release together with the stock floor
// contact (the kicker's own impulses against the floor's embed correction), so the playfield model leaves it to them.
bool InPlayfieldHole(const Vertex3Ds &pos)
{
   for (const RubberContact::SurfaceHole &h : g_playfieldHoles)
   {
      const double dx = (double)pos.x * U - h.x, dy = (double)pos.y * U - h.y;
      if (dx * dx + dy * dy < h.r * h.r)
         return true;
   }
   return false;
}

// Band detection (step 1 of the band model: find the posts and the free spans; the contact itself is still the rigid tube).
// A VPX rubber stores only a closed centreline, so the posts are recovered from it: a vertex whose local bend radius is below
// band_curve_radius_mm is on a post wrap, a run of such vertices is one post (a circle is fitted to it), and the straight runs between
// consecutive posts are the spans. band_posts = "x,y,r | ..." replaces the fit with given circles (vertices within 15% of r of a circle are on it).
void HitRubber::DetectBand()
{
   m_bandPosts.clear();
   m_bandSpans.clear();
   const int n = (int)m_pts.size() - 1;   // the loop closed by repeating the first vertex
   if (n < 4)
   {
      PLOGW << "Band " << m_name << ": only " << n << " centreline points, no band";
      return;
   }
   const double mmPerU = U * 1000.0;
   auto at = [&](int i) -> const Vertex3Ds & { return m_pts[((i % n) + n) % n]; };

   // posts given by the user, if any
   struct Circle { double x, y, r; };
   std::vector<Circle> given;
   {
      const std::string rest = m_set.bandPosts;
      size_t pos = 0;
      while (pos < rest.size())
      {
         size_t bar = rest.find('|', pos);
         if (bar == std::string::npos)
            bar = rest.size();
         const std::string item = rest.substr(pos, bar - pos);
         double a = 0, b = 0, c = 0;
         if (sscanf(item.c_str(), " %lf , %lf , %lf", &a, &b, &c) == 3 && c > 0.0)
            given.push_back({ a, b, c });
         else if (!Trim(item).empty())
            PLOGW << "Band " << m_name << ": ignored band_posts entry '" << Trim(item) << "' (want x,y,r)";
         pos = bar + 1;
      }
   }

   // which vertices are on a post wrap; owner = index into `given` (override mode), or -1
   std::vector<char> onPost(n, 0);
   std::vector<int> owner(n, -1);
   const double rCurve = m_set.bandCurveRadiusMm / mmPerU;   // VPX units
   for (int i = 0; i < n; i++)
   {
      if (!given.empty())
      {
         for (size_t k = 0; k < given.size(); k++)
         {
            const double d = std::hypot((double)at(i).x - given[k].x, (double)at(i).y - given[k].y);
            if (std::fabs(d - given[k].r) <= std::max(0.15 * given[k].r, 1.0))
            {
               onPost[i] = 1;
               owner[i] = (int)k;
               break;
            }
         }
      }
      else
      {
         const double ax = at(i).x - at(i - 1).x, ay = at(i).y - at(i - 1).y;
         const double bx = at(i + 1).x - at(i).x, by = at(i + 1).y - at(i).y;
         const double la = std::hypot(ax, ay), lb = std::hypot(bx, by);
         if (la < 1e-6 || lb < 1e-6)
            continue;
         const double turn = std::atan2(ax * by - ay * bx, ax * bx + ay * by);   // signed turning angle at the vertex
         const double k = std::fabs(turn) / (0.5 * (la + lb));                    // curvature, 1/U
         onPost[i] = (k > 1.0 / rCurve) ? 1 : 0;
      }
   }

   // runs of post vertices around the loop, starting just after a free vertex
   int start = -1;
   for (int i = 0; i < n; i++)
      if (!onPost[i])
      {
         start = i;
         break;
      }
   if (start < 0)
   {
      PLOGW << "Band " << m_name << ": every centreline vertex bends tighter than " << m_set.bandCurveRadiusMm << " mm, so there are no straight spans (a ring, not a band)";
      return;
   }
   for (int k = 1; k < n;)
   {
      const int i = start + k;
      if (!onPost[i % n])
      {
         k++;
         continue;
      }
      BandPost bp;
      bp.first = i;
      int j = i;
      while (j - i < n && onPost[(j + 1) % n] && (given.empty() || owner[(j + 1) % n] == owner[i % n]))
         j++;
      bp.last = j;
      const int cnt = j - i + 1;
      if (!given.empty())
      {
         bp.x = given[owner[i % n]].x;
         bp.y = given[owner[i % n]].y;
         bp.r = given[owner[i % n]].r;
      }
      else
      {
         // algebraic circle fit (Kasa) to the run, with the vertex on each side (the tangent points) included
         std::vector<std::pair<double, double>> q;
         for (int m = i - 1; m <= j + 1; m++)
            q.emplace_back(at(m).x, at(m).y);
         double mx = 0, my = 0;
         for (const auto &p : q) { mx += p.first; my += p.second; }
         mx /= (double)q.size();
         my /= (double)q.size();
         double suu = 0, suv = 0, svv = 0, suuu = 0, svvv = 0, suvv = 0, svuu = 0;
         for (const auto &p : q)
         {
            const double u = p.first - mx, v = p.second - my;
            suu += u * u; suv += u * v; svv += v * v; suuu += u * u * u; svvv += v * v * v; suvv += u * v * v; svuu += v * u * u;
         }
         const double det = suu * svv - suv * suv;
         bp.fitted = cnt >= 3 && std::fabs(det) > 1e-9 * (suu + svv) * (suu + svv);
         if (bp.fitted)
         {
            const double b1 = 0.5 * (suuu + suvv), b2 = 0.5 * (svvv + svuu);
            const double uc = (b1 * svv - b2 * suv) / det, vc = (b2 * suu - b1 * suv) / det;
            bp.x = mx + uc;
            bp.y = my + vc;
            bp.r = std::sqrt(uc * uc + vc * vc + (suu + svv) / (double)q.size());
         }
         else
         {
            bp.x = mx;
            bp.y = my;
            double rr = 0;
            for (const auto &p : q)
               rr += std::hypot(p.first - mx, p.second - my);
            bp.r = rr / (double)q.size();
         }
      }
      m_bandPosts.push_back(bp);
      k = j - start + 1;
   }
   if (m_bandPosts.size() < 2)
   {
      PLOGW << "Band " << m_name << ": found " << m_bandPosts.size() << " post(s), a band needs at least 2";
      return;
   }

   // spans: the vertices between one post run and the next
   const size_t np = m_bandPosts.size();
   for (size_t a = 0; a < np; a++)
   {
      const BandPost &pa = m_bandPosts[a];
      const BandPost &pb = m_bandPosts[(a + 1) % np];
      BandSpan sp;
      sp.postA = (int)a;
      sp.postB = (int)((a + 1) % np);
      sp.first = pa.last;                                         // the last post vertex is where the straight run leaves the wrap
      sp.last = pb.first + (pb.first <= pa.last ? n : 0);         // and the first post vertex of the next post is where it arrives
      double len = 0.0;
      for (int m = sp.first; m < sp.last; m++)
         len += std::hypot((double)at(m + 1).x - at(m).x, (double)at(m + 1).y - at(m).y);
      sp.lengthMm = len * mmPerU;
      sp.free = sp.lengthMm >= m_set.bandMinSpanMm;
      m_bandSpans.push_back(sp);
   }

   PLOGI << "Band " << m_name << ": " << (given.empty() ? "fitted" : "given") << ' ' << np << " posts, " << m_bandSpans.size() << " spans, tube radius "
         << m_radius * mmPerU << " mm, " << n << " centreline points";
   for (size_t a = 0; a < np; a++)
      PLOGI << "  post " << a << ": centre (" << m_bandPosts[a].x << ", " << m_bandPosts[a].y << ") VPX units, centreline radius " << m_bandPosts[a].r << " U = "
            << m_bandPosts[a].r * mmPerU << " mm, post radius ~" << (m_bandPosts[a].r - m_radius) * mmPerU << " mm, vertices " << m_bandPosts[a].first << ".." << m_bandPosts[a].last
            << (m_bandPosts[a].fitted ? "" : " (too few vertices to fit a circle: centre is the mean)");
   for (size_t a = 0; a < m_bandSpans.size(); a++)
      PLOGI << "  span " << a << ": post " << m_bandSpans[a].postA << " -> post " << m_bandSpans[a].postB << ", vertices " << m_bandSpans[a].first << ".." << m_bandSpans[a].last
            << ", length " << m_bandSpans[a].lengthMm << " mm, " << (m_bandSpans[a].free ? "free" : "short (stays rigid)");
}

// ---------------------------------------------------------------------------
// Sling: the plate and the switch blades are the table's own primitives. Nothing is written to the table file; the plate's pose is changed in memory.
void HitRubber::SlingAttach(Primitive *const plate, const std::vector<Primitive *> &switches)
{
   if (m_bandSpans.empty())
      return;   // a ring (see the constructor): nothing to put a plate behind
   m_slingPlate = plate;
   m_slingSwitches = switches;
   if (plate == nullptr)
      return;
   m_slingBase = plate->m_d.m_vPosition;
   m_slingShownRotX = plate->m_d.m_aRotAndTra[0];
   // the ball must pass through these: only the band touches it
   plate->put_Collidable(VARIANT_FALSE);
   m_switches.clear();
   for (Primitive *const sw : switches)
   {
      sw->put_Collidable(VARIANT_FALSE);
      BandSwitch ss;
      ss.prim = sw;
      ss.baseTransX = sw->m_d.m_aRotAndTra[3];   // the table's value is the blade's resting place; the engine slides it with the band
      ss.shownTransX = ss.baseTransX;
      ss.coil = true;
      ss.legacySign = true;
      m_switches.push_back(ss);
   }
   m_plate.Configure(m_set.sling, true);
   m_plateConfigured = true;
   m_driver.SetPulse(m_set.sling.pulseMs * 1e-3);
   m_driver.SetDelay(m_set.sling.coilDelayMs * 1e-3);
   PLOGI << "Sling " << m_name << ": plate " << plate->GetName() << " (pivot " << m_slingBase.x << ", " << m_slingBase.y << ", " << m_slingBase.z << "), " << switches.size()
         << " switch blades, collision off on all of them";
   SlingPlace();
}

void HitRubber::SlingPlace() const
{
   Primitive *const pl = m_slingPlate;
   if (pl == nullptr || m_bandSpans.empty())
      return;
   const int n = (int)m_pts.size() - 1;
   if (n < 4)
      return;
   // measure from the table's own pose, with the plate at rest (the engine may have it swung)
   const float shownRotX = pl->m_d.m_aRotAndTra[0];
   pl->m_d.m_aRotAndTra[0] = 0.f;
   pl->m_d.m_vPosition = m_slingBase;
   const Matrix3D &M = pl->RecalculateMatrices();
   std::vector<Vertex3Ds> w(pl->m_mesh.NumVertices());
   double cx = 0.0, cy = 0.0;
   for (size_t i = 0; i < w.size(); i++)
   {
      w[i] = M * pl->m_mesh.m_vertices[i];
      cx += w[i].x;
      cy += w[i].y;
   }
   if (w.empty())
      return;
   cx /= (double)w.size();
   cy /= (double)w.size();

   // the span the face looks at: for every free span, the face is the plate's surface at the band's height that lies along the span; its distance from the back of the tube is the gap. The span with the smallest gap wins.
   double bestGap = 1e30, nfx = 0.0, nfy = 0.0;
   int bestSpan = -1;
   std::vector<std::pair<double, double>> facePts, bestPts;   // (place along the chord, distance from the rest centreline) of the plate's surface at band height
   for (size_t k = 0; k < m_bandSpans.size(); k++)
   {
      const BandSpan &sp = m_bandSpans[k];
      if (!sp.free)
         continue;
      const Vertex3Ds &A = m_pts[((sp.first % n) + n) % n];
      const Vertex3Ds &B = m_pts[((sp.last % n) + n) % n];
      double tx = (double)B.x - A.x, ty = (double)B.y - A.y;
      const double len = std::sqrt(tx * tx + ty * ty);
      if (len < 1e-6)
         continue;
      tx /= len;
      ty /= len;
      double fx = ty, fy = -tx;                                   // a normal of the span
      if (((cx - A.x) * fx + (cy - A.y) * fy) > 0.0)              // make it point from the plate toward the band
      {
         fx = -fx;
         fy = -fy;
      }
      const double bandZ = 0.5 * ((double)A.z + B.z);
      double face = -1e30;
      // the plate's surface at the band's height: where the mesh edges cross the plane z = bandZ (the mesh may have vertices only at the ends of the face,
      // and the hood, above the band, never reaches this plane)
      facePts.clear();
      auto consider = [&](const Vertex3Ds &v) {
         const double u = ((double)v.x - A.x) * tx + ((double)v.y - A.y) * ty;
         if (u < 0.0 || u > len)
            return;
         const double dn = ((double)v.x - A.x) * fx + ((double)v.y - A.y) * fy;
         facePts.emplace_back(u, dn);
         face = std::max(face, dn);
      };
      const size_t ni = pl->m_mesh.NumIndices();
      for (size_t i = 0; i + 2 < ni; i += 3)
         for (int e = 0; e < 3; e++)
         {
            const Vertex3Ds &p0 = w[pl->m_mesh.m_indices[i + e]];
            const Vertex3Ds &p1 = w[pl->m_mesh.m_indices[i + (e + 1) % 3]];
            const double d0 = (double)p0.z - bandZ, d1 = (double)p1.z - bandZ;
            if (d0 * d1 > 0.0)
               continue;
            if (d0 == d1)
            {
               consider(p0);
               consider(p1);
               continue;
            }
            const double t = d0 / (d0 - d1);
            consider(Vertex3Ds(p0.x + (float)(t * (p1.x - p0.x)), p0.y + (float)(t * (p1.y - p0.y)), (float)bandZ));
         }
      if (face < -1e29 || face > 0.0)
         continue;
      const double gap = -face - m_radius;
      if (gap < bestGap)
      {
         bestGap = gap;
         bestSpan = (int)k;
         nfx = fx;
         nfy = fy;
         bestPts = facePts;
      }
   }
   if (bestSpan < 0)
   {
      PLOGW << "Sling " << m_name << ": plate " << pl->GetName() << " has no face behind a free span of the band at the band's height; the plate is left where the table has it";
      pl->m_d.m_aRotAndTra[0] = shownRotX;
      pl->RecalculateMatrices();
      return;
   }
   constexpr double mmPerU = U * 1000.0;
   m_slingGapTableMm = bestGap * mmPerU;
   const double shiftU = bestGap;   // the plate goes where its face touches the band at RotX = 0; the gap at rest is made by turning the plate back (below)
   pl->m_d.m_vPosition.x = m_slingBase.x + (float)(nfx * shiftU);
   pl->m_d.m_vPosition.y = m_slingBase.y + (float)(nfy * shiftU);
   pl->m_d.m_aRotAndTra[0] = 0.f;
   pl->RecalculateMatrices();
   // the switches sense the band where their blades are: the place along the span's chord
   {
      const BandSpan &sp = m_bandSpans[bestSpan];
      const Vertex3Ds &A = m_pts[((sp.first % n) + n) % n];
      const Vertex3Ds &B = m_pts[((sp.last % n) + n) % n];
      const double abx = (double)B.x - A.x, aby = (double)B.y - A.y, l2 = abx * abx + aby * aby;
      m_slingSpan = bestSpan;
      m_slingNx = nfx;
      m_slingNy = nfy;
      // the plate's geometry for the contact: the face (its surface within half a unit of the front) and the pivot axis, which runs through the plate's position
      const double faceRest = -m_radius;   // VPX units, from the rest centreline: after the move the face touches the back of the band at RotX = 0
      double uLo = 1e30, uHi = -1e30;
      for (const auto &pt : bestPts)
         if (pt.second + shiftU >= faceRest - 0.6)
         {
            uLo = std::min(uLo, pt.first);
            uHi = std::max(uHi, pt.first);
         }
      if (uHi > uLo && std::sqrt(abx * abx + aby * aby) > 1e-6)
      {
         const double len = std::sqrt(abx * abx + aby * aby);
         const double originN = (((double)pl->m_d.m_vPosition.x - A.x) * nfx + ((double)pl->m_d.m_vPosition.y - A.y) * nfy);
         m_plateOriginN = originN * U;
         m_plateA = (faceRest - originN) * U;
         m_plateH = (0.5 * ((double)A.z + B.z) - (double)pl->m_d.m_vPosition.z) * U;
         m_plateF0 = uLo / len;
         m_plateF1 = uHi / len;
         m_plateGeomValid = m_plateA > 0.0 && m_plateH > 0.0;
         // the gap at rest: the plate stands turned back from the touching pose by the angle at which the face is slingGapMm from the band
         //   a (1 / cos t - 1) + h tan t = -gap   for t < 0   (monotonic in t as h > a)
         m_plateRest = 0.0;
         if (m_plateGeomValid && m_set.slingGapMm > 0.0)
         {
            const double gapM = m_set.slingGapMm * 1e-3;
            double lo = -1.2, hi = 0.0;
            for (int it = 0; it < 60; it++)
            {
               const double mid = 0.5 * (lo + hi);
               const double g = m_plateA * (1.0 / std::cos(mid) - 1.0) + m_plateH * std::tan(mid);
               if (g > -gapM)
                  hi = mid;
               else
                  lo = mid;
            }
            m_plateRest = 0.5 * (lo + hi);
         }
         const float rx = (float)((m_plate.PhiDeg() + m_plateRest * 180.0 / M_PI));
         pl->m_d.m_aRotAndTra[0] = rx;
         m_slingShownRotX = rx;
         pl->RecalculateMatrices();
      }
      for (BandSwitch &sw : m_switches)
      {
         if (!sw.coil)
            continue;
         sw.span = bestSpan;
         sw.bx = -nfx;   // behind the band is toward the plate
         sw.by = -nfy;
         // the blade's place in the table: its position moved by the table's TransX along its own X axis (turned by its RotZ)
         const double rz = (double)sw.prim->m_d.m_aRotAndTra[2] * M_PI / 180.0;
         const double px = (double)sw.prim->m_d.m_vPosition.x + std::cos(rz) * sw.baseTransX, py = (double)sw.prim->m_d.m_vPosition.y + std::sin(rz) * sw.baseTransX;
         sw.f = std::min(std::max(((px - A.x) * abx + (py - A.y) * aby) / std::max(l2, 1e-12), 0.0), 1.0);
      }
   }
   if (std::fabs(m_slingPlacedGapMm - m_set.slingGapMm) > 1e-9)
   {
      m_slingPlacedGapMm = m_set.slingGapMm;
      PLOGI << "Sling " << m_name << ": plate " << pl->GetName() << " faces span " << bestSpan << ", face-to-band gap in the table " << m_slingGapTableMm << " mm, set to " << m_set.slingGapMm
            << " mm: the plate is moved " << std::fabs(shiftU) * mmPerU << (shiftU >= 0.0 ? " mm toward" : " mm away from") << " the band (" << nfx << ", " << nfy << ") so that its face touches the band at RotX = 0, pivot now ("
            << pl->m_d.m_vPosition.x << ", " << pl->m_d.m_vPosition.y << "), and it rests turned back to RotX = " << m_plateRest * 180.0 / M_PI << " deg for the gap";
   }
}

// How far the band is pushed from the ball toward the blade at a switch, m (> 0 = toward the blade, 0 while the band is at rest).
double HitRubber::SwitchDisplacementM(const BandSwitch &sw) const
{
   using namespace RubberContact;
   if (!m_dyn.active || sw.span < 0)
      return 0.0;
   int ci = -1;
   for (size_t k = 0; k < m_dyn.span.size(); k++)
      if (m_dyn.span[k] == sw.span)
         ci = (int)k;
   if (ci < 0 || m_dyn.chains[ci].Count() == 0)
      return 0.0;
   const BandChain &ch = m_dyn.chains[ci];
   const double f = sw.f;
   int j = 0;   // the leg that holds the material fraction f
   while (j < ch.Count() && ch.PtF(j + 1) < f)
      j++;
   const double f0 = ch.PtF(j), f1 = ch.PtF(j + 1);
   const double u = f1 > f0 ? std::min(std::max((f - f0) / (f1 - f0), 0.0), 1.0) : 0.0;
   const V3 P = ch.PtPos(j) + (ch.PtPos(j + 1) - ch.PtPos(j)) * u;
   const V3 rest = ch.A + (ch.B - ch.A) * f;
   return (P.x - rest.x) * sw.bx + (P.y - rest.y) * sw.by;
}

// Switch blades behind the band that only the script reads: collision off, the span they sense and the side they are on found from where they stand.
void HitRubber::BandSwitchesAttach(const std::vector<Primitive *> &blades)
{
   if (m_bandSpans.empty())
      return;
   for (Primitive *const prim : blades)
   {
      prim->put_Collidable(VARIANT_FALSE);
      BandSwitch ss;
      ss.prim = prim;
      ss.baseTransX = prim->m_d.m_aRotAndTra[3];
      ss.shownTransX = ss.baseTransX;
      SwitchPlace(ss);
      m_switches.push_back(ss);
   }
   PLOGI << "Band " << m_name << ": " << blades.size() << " switch blades for the script, " << m_switches.size() << " in all, collision off on them";
}

void HitRubber::SwitchPlace(BandSwitch &sw) const
{
   const int n = (int)m_pts.size() - 1;
   if (n < 4 || m_bandSpans.empty())
      return;
   // the blade's place in the table: its position moved by the table's TransX along its own X axis (turned by its RotZ)
   const double rz = (double)sw.prim->m_d.m_aRotAndTra[2] * M_PI / 180.0;
   const double px = (double)sw.prim->m_d.m_vPosition.x + std::cos(rz) * sw.baseTransX, py = (double)sw.prim->m_d.m_vPosition.y + std::sin(rz) * sw.baseTransX;
   int best = -1;
   double bestD = 1e30, bestT = 0.0, bnx = 0.0, bny = 0.0, bestS = 0.0;
   for (size_t k = 0; k < m_bandSpans.size(); k++)
   {
      const BandSpan &sp = m_bandSpans[k];
      if (!sp.free)
         continue;
      const Vertex3Ds &A = m_pts[((sp.first % n) + n) % n];
      const Vertex3Ds &B = m_pts[((sp.last % n) + n) % n];
      const double abx = (double)B.x - A.x, aby = (double)B.y - A.y, l2 = abx * abx + aby * aby;
      if (l2 < 1e-9)
         continue;
      const double t = std::min(std::max(((px - A.x) * abx + (py - A.y) * aby) / l2, 0.0), 1.0);
      const double dx = px - (A.x + abx * t), dy = py - (A.y + aby * t);
      const double d = std::sqrt(dx * dx + dy * dy);
      if (d < bestD)
      {
         const double len = std::sqrt(l2);
         best = (int)k;
         bestD = d;
         bestT = t;
         bnx = aby / len;
         bny = -abx / len;
         bestS = (px - A.x) * bnx + (py - A.y) * bny;
      }
   }
   if (best < 0)
   {
      PLOGW << "Band " << m_name << ": switch blade " << sw.prim->GetName() << " has no free span to sense";
      return;
   }
   const double side = bestS >= 0.0 ? 1.0 : -1.0;
   if (std::fabs(bestS) < 1.0)
      PLOGW << "Band " << m_name << ": switch blade " << sw.prim->GetName() << " stands on the band's chord (" << bestS * U * 1e3 << " mm off it): which side is behind is a guess";
   sw.span = best;
   sw.f = bestT;
   sw.bx = bnx * side;
   sw.by = bny * side;
   sw.axisProj = std::cos(rz) * sw.bx + std::sin(rz) * sw.by;
   PLOGI << "Band " << m_name << ": switch blade " << sw.prim->GetName() << " senses span " << best << " at " << bestT << " of its chord, " << std::fabs(bestS) * U * 1e3
         << " mm behind it, back direction (" << sw.bx << ", " << sw.by << "), blade X axis along it " << sw.axisProj;
}

// Once per physics cycle: every blade reads the band; a switch closes when the band is pushed back by the switch gap, opens again when it no longer is, but not before
// its minimum closed time. The script is told of every change. Returns true when a sling (coil) switch is closed.
bool HitRubber::SwitchesUpdate(const double dt)
{
   std::vector<std::pair<int, bool>> changes;
   bool anyCoil = false;
   for (size_t i = 0; i < m_switches.size(); i++)
   {
      BandSwitch &sw = m_switches[i];
      const double gapM = (sw.coil ? m_set.sling.switchGapMm : m_set.switchGapMm) * 1e-3;
      const double minS = (sw.coil ? m_set.sling.minClosedMs : m_set.switchMinClosedMs) * 1e-3;
      const bool was = sw.closed;
      sw.dispM = SwitchDisplacementM(sw);
      const bool raw = sw.dispM >= gapM;
      if (raw && !sw.closed)
      {
         sw.closed = true;
         sw.closedForS = 0.0;
      }
      else if (sw.closed)
      {
         sw.closedForS += dt;
         if (!raw && sw.closedForS >= minS)
            sw.closed = false;
      }
      if (sw.closed != was)
         changes.emplace_back((int)i, sw.closed);
      if (sw.coil)
         anyCoil = anyCoil || sw.closed;
   }
   if (m_editable != nullptr)
      for (const auto &c : changes)   // the script handler may do anything: nothing of ours is in flight (the vector is not touched by it)
         static_cast<Rubber *>(m_editable)->FireVoidEventParm(c.second ? DISPID_RubberEvents_SwitchClosed : DISPID_RubberEvents_SwitchOpened, c.first);
   return anyCoil;
}

// The blades follow the band: pushed back toward a blade, the blade slides with it; when the band stretches toward the ball the blade stays at its resting place.
void HitRubber::SwitchesFinish()
{
   for (BandSwitch &sw : m_switches)
   {
      const double push = std::max(sw.dispM, 0.0) / U;   // VPX units
      const float tx = (!sw.legacySign && std::fabs(sw.axisProj) >= 0.2) ? sw.baseTransX + (float)(push / sw.axisProj) : sw.baseTransX - (float)push;
      if (std::fabs(tx - sw.shownTransX) > 1e-3f)
      {
         sw.shownTransX = tx;
         sw.prim->m_d.m_aRotAndTra[3] = tx;
         sw.prim->RecalculateMatrices();
      }
   }
}

// Once per physics cycle: the switches read the band, the driver turns a closing into a pulse, the plate body steps with the coil on while it pulses, and
// the table's primitives follow (the plate's RotX, the blades' TransX). The plate does not touch the band yet.
void HitRubber::SlingAdvance(const float dtime)
{
   if (m_slingPlate == nullptr || !m_plateConfigured)
      return;
   const double dt = (double)dtime * DEFAULT_STEPTIME_S;
   if (dt <= 0.0)
      return;
   const SlingSettings &ss = m_set.sling;
   m_driver.SetPulse(ss.pulseMs * 1e-3);
   m_driver.SetDelay(ss.coilDelayMs * 1e-3);

   // a plate that has left the model's range (not a number, or far past its coil stop) is put back at rest, and the band with it, instead of staying broken
   {
      const auto &st = m_plate.st;
      const double sum = st.phi + st.omega + st.xp + st.vp + st.iHigh;
      const double limit = 3.0 * ss.mech.stopDeg * FlipperModel::DEG + 0.5;
      if (!std::isfinite(sum) || std::fabs(st.phi) > limit)
      {
         PLOGW << "Sling " << m_name << ": the plate left the model's range (angle " << st.phi / FlipperModel::DEG << " deg), put back at rest";
         m_plate.Configure(ss, true);
         m_driver.Cancel();
         m_slingTriggerReq = false;
         BandReleaseAll("sling reset");
         m_dyn = BandDyn();
         BandPublish();
         m_slingIdle = false;
      }
   }

   // switches (the sling's coil switches drive the pulse; the script is told of every switch)
   const bool any = SwitchesUpdate(dt);
   m_slingAnyClosed = any;
   const int firedBefore = m_driver.fired;
   m_slingTriggered = false;
   if (m_slingEnabled)
   {
      if (m_slingTriggerReq)
         m_slingTriggered = m_driver.Trigger();
      m_driver.Update(any);
   }
   else
      m_driver.Cancel();
   m_slingTriggerReq = false;
   if (m_driver.fired != firedBefore)
   {
      m_slingIdle = false;
      m_slingInfo.fireTimeMs = (int)g_pplayer->m_time_msec;
      PLOGI << "Sling " << m_name << (m_slingTriggered ? ": script trigger" : ": switch closed") << " at t=" << g_pplayer->m_time_msec << " ms (band " << (m_switches.empty() ? 0.0 : m_switches[0].dispM * 1e3) << ", "
            << (m_switches.size() > 1 ? m_switches[1].dispM * 1e3 : 0.0) << " mm back), coil pulse " << ss.pulseMs << " ms";
      m_slingInfo.pulses = m_driver.fired;
      static_cast<Rubber *>(m_editable)->FireGroupEvent(DISPID_RubberEvents_SlingFire);   // the script handler may do anything: nothing of ours is in flight yet this cycle
   }

   // does the plate work on the band this cycle? Yes once the face is about to reach the cord (a cycle ahead, and half a millimetre of margin) or while the
   // band is deflected; then the plate steps together with the band (SlingBandStep inside BandAdvanceCore). Otherwise on its own.
   const double phiAhead = m_plate.st.phi + 1.2 * m_plate.st.omega * dt;
   m_slingInBand = m_plateGeomValid && (m_dyn.active || SlingFaceN(phiAhead) + m_radiusSI > -0.5e-3);
   if (m_slingInBand)
   {
      m_slingIdle = false;
      SlingEnsureNodes();
      return;
   }

   // the plate alone: a whole quiet step (coil off, nothing moved) is skipped from then on until something changes
   const bool coilNow = m_driver.Coil();
   const double hMax = 2e-5;
   if (!(m_slingIdle && !coilNow && m_plate.st == m_slingIdleState))
   {
      const SlingPlate::State s0 = m_plate.st;
      const int n = std::max(1, (int)std::ceil(dt / hMax));
      const double h = dt / n;
      for (int i = 0; i < n; i++)
         SlingStepPlate(h, 0.0);
      const SlingPlate::State &s1 = m_plate.st;
      m_slingIdle = dt >= 0.9 * PHYSICS_STEPTIME_S && !coilNow && m_driver.remaining <= 0.0 && std::fabs(s1.phi - s0.phi) < 1e-12 && std::fabs(s1.omega - s0.omega) < 1e-9
         && std::fabs(s1.xp - s0.xp) < 1e-12 && std::fabs(s1.vp - s0.vp) < 1e-9 && std::fabs(s1.iHigh - s0.iHigh) < 1e-12 && s1.rigid == s0.rigid;
      if (m_slingIdle)
         m_slingIdleState = s1;
   }
}

// one plate substep with the driver's coil state; logs the first arrival at the coil stop of a stroke
void HitRubber::SlingStepPlate(const double h, const double torqueNm)
{
   m_plate.Step(h, m_driver.Coil(), torqueNm);
   m_driver.Tick(h);
   if (m_plate.HitStop() && !m_slingStopLogged)
   {
      m_slingStopLogged = true;
      PLOGI << "Sling " << m_name << ": plate at the coil stop at t=" << g_pplayer->m_time_msec << " ms, angle " << m_plate.PhiDeg() << " deg, speed " << m_plate.st.omega << " rad/s";
      {
         // log only: the shape of the span the plate faces at that moment (place along the chord, displacement along the normal, mm, + toward the ball)
         const int ci = SlingChain();
         if (m_dyn.active && ci >= 0)
         {
            std::ostringstream o;
            for (const int id : m_dyn.chains[ci].order)
            {
               const RubberContact::BandChain::Node &nd = m_dyn.chains[ci].nodes[id];
               o << " (" << (double)((int)(nd.f * 100.0 + 0.5)) / 100.0 << ": " << (double)((int)(((nd.P.x - nd.P0.x) * m_slingNx + (nd.P.y - nd.P0.y) * m_slingNy) * 1e4 + 0.5 * 1.0)) / 10.0 << ")";
            }
            PLOGI << "Sling " << m_name << ": span shape at the stop, f along the chord: displacement mm:" << o.str();
         }
      }
   }
   if (m_slingStopLogged && m_driver.remaining <= 0.0 && m_plate.PhiDeg() < 1.0)
      m_slingStopLogged = false;
}

// the table's primitives follow the plate and the switches
void HitRubber::SlingFinish()
{
   m_slingInfo.enabled = m_slingEnabled;
   m_slingInfo.plateAngleDeg = (float)m_plate.PhiDeg();
   m_slingInfo.switchClosed = m_slingAnyClosed;
   m_slingInfo.coilOn = m_driver.Coil();
   m_slingInfo.coilCurrentA = (float)m_plate.CurrentA();
   m_slingInfo.pulses = m_driver.fired;
   {
      // log only: how far the span the plate faces is displaced, back (toward the plate) and forward (toward the ball), over one deflection of the band
      const int ci = SlingChain();
      if (m_dyn.active && ci >= 0)
         for (const int id : m_dyn.chains[ci].order)
         {
            const RubberContact::BandChain::Node &nd = m_dyn.chains[ci].nodes[id];
            const double xn = ((nd.P.x - nd.P0.x) * m_slingNx + (nd.P.y - nd.P0.y) * m_slingNy) * 1e3;
            m_slingPkFwdMm = std::max(m_slingPkFwdMm, xn);
            m_slingPkBackMm = std::max(m_slingPkBackMm, -xn);
         }
      if (m_slingWasActive && !m_dyn.active)
      {
         PLOGI << "Sling " << m_name << ": band back to rest at t=" << g_pplayer->m_time_msec << " ms; largest displacement of the span: " << m_slingPkBackMm << " mm toward the plate, "
               << m_slingPkFwdMm << " mm toward the ball (past its rest line)";
         m_slingPkBackMm = m_slingPkFwdMm = 0.0;
      }
      m_slingWasActive = m_dyn.active;
   }
   const float rotX = (float)(m_plate.PhiDeg() + m_plateRest * 180.0 / M_PI);   // RotX = 0 is the pose in which the face touches the band
   if (std::fabs(rotX - m_slingShownRotX) > 1e-4f)
   {
      m_slingShownRotX = rotX;
      m_slingPlate->m_d.m_aRotAndTra[0] = rotX;
      m_slingPlate->RecalculateMatrices();
   }
   SwitchesFinish();
}

double HitRubber::SlingFaceN(const double phi) const
{
   const double t = phi + m_plateRest;   // the angle from the touching pose (RotX)
   const double c = std::cos(t);
   return m_plateOriginN + m_plateA / c + m_plateH * std::tan(t);
}

double HitRubber::SlingFaceDN(const double phi) const
{
   const double t = phi + m_plateRest;
   const double c = std::cos(t);
   return (m_plateA * std::sin(t) + m_plateH) / (c * c);
}

int HitRubber::SlingChain() const
{
   for (size_t k = 0; k < m_dyn.span.size(); k++)
      if (m_dyn.span[k] == m_slingSpan)
         return (int)k;
   return -1;
}

bool HitRubber::SlingInFootprint(const int chain, const int node) const
{
   if (!m_slingInBand || chain != SlingChain() || chain < 0)
      return false;
   const RubberContact::BandChain &ch = m_dyn.chains[chain];
   const double tol = m_set.bandFollowMm * 1e-3 / std::max(ch.Ls, 1e-6);   // a node within a node spacing of the face carries its share (nodes closer than that are unstable)
   const double f = ch.nodes[node].f;
   return f >= m_plateF0 - tol && f <= m_plateF1 + tol;
}

// Nodes at the two ends of the face, so the cord between them can lie flat on it (a face narrower than 6 mm gets one node in the middle).
void HitRubber::SlingEnsureNodes()
{
   if (!BandEnsureDyn())
      return;
   const int ci = SlingChain();
   if (ci < 0)
      return;
   RubberContact::BandChain &ch = m_dyn.chains[ci];
   const double widthM = (m_plateF1 - m_plateF0) * ch.Ls;
   double fs[2] = { m_plateF0, m_plateF1 };
   int nf = 2;
   if (widthM < 6e-3)
   {
      fs[0] = 0.5 * (m_plateF0 + m_plateF1);
      nf = 1;
   }
   bool changed = false;
   for (int k = 0; k < nf; k++)
   {
      bool have = false;
      for (const int id : ch.order)
         have = have || std::fabs(ch.nodes[id].f - fs[k]) * ch.Ls < m_set.bandFollowMm * 1e-3;   // no node closer than the follow spacing: a ball's node there serves
      if (!have && ch.Count() < 12)
      {
         ch.Add(fs[k]);
         changed = true;
      }
   }
   if (changed)
      BandRebuild();
}

// One substep of the plate working on the band: the face presses on every node under it (a spring and damper from the rubber's own modulus, one tributary
// width each), the plate steps with the reaction torque, and the driver's time runs.
void HitRubber::SlingBandStep(std::vector<std::vector<RubberContact::V3>> &fNode, const double h)
{
   using namespace RubberContact;
   double torque = 0.0;
   const int ci = SlingChain();
   if (ci >= 0 && ci < (int)m_dyn.chains.size() && m_dyn.chains[ci].Count() > 0)
   {
      BandChain &ch = m_dyn.chains[ci];
      const double phi = m_plate.st.phi;
      const double face = SlingFaceN(phi), dFaceDPhi = SlingFaceDN(phi), vFace = dFaceDPhi * m_plate.st.omega;
      int count = 0;
      for (const int id : ch.order)
         count += SlingInFootprint(ci, id) ? 1 : 0;
      if (count > 0)
      {
         const double width = std::max((m_plateF1 - m_plateF0) * ch.Ls, 3e-3) / count;
         const double E = m_set.eMpa * 1e6;
         for (const int id : ch.order)
         {
            if (!SlingInFootprint(ci, id))
               continue;
            const BandChain::Node &nd = ch.nodes[id];
            const double xn = (nd.P.x - nd.P0.x) * m_slingNx + (nd.P.y - nd.P0.y) * m_slingNy;   // + toward the ball
            const double pen = face - xn + m_radiusSI;                                           // the face against the back of the cord
            if (pen <= 0.0)
               continue;
            const double vn = nd.vp.x * m_slingNx + nd.vp.y * m_slingNy;
            const double k = E * width, c = 2.0 * m_set.sling.plateContactZeta * std::sqrt(k * nd.m);
            const double F = k * pen + c * (vFace - vn);
            if (F <= 0.0)
               continue;
            fNode[ci][id] = fNode[ci][id] + V3(m_slingNx, m_slingNy, 0.0) * F;
            torque -= F * dFaceDPhi;
         }
      }
   }
   SlingStepPlate(h, torque);
}

void HitRubber::CalcHitBBox()
{
   if (m_surface != nullptr)
   {
      const float m = 1.0f;
      m_hitBBox = FRect3D(m_surface->lo[0] - m, m_surface->hi[0] + m, m_surface->lo[1] - m, m_surface->hi[1] + m, m_surface->lo[2] - m, m_surface->hi[2] + m);
      if (const char *nearArg = ContactSystem::EnvVar("VPX_CONTACT_NEAR", "VPX_RING_NEAR"); nearArg != nullptr && m_set.contact)   // "x,y,r" (VPX units): the faces within r, as a 4 U grid of cells (count, z range)
      {
         double cx = 0, cy = 0, cr = 0;
         if (sscanf(nearArg, "%lf,%lf,%lf", &cx, &cy, &cr) == 3)
         {
            std::map<std::pair<int, int>, std::array<double, 3>> cells;   // count, zmin, zmax
            constexpr double U = 0.53975e-3;
            for (const RubberContact::SurfaceTri &t : m_surface->tris)
            {
               const double x = (t.a.x + t.b.x + t.c.x) / (3.0 * U), y = (t.a.y + t.b.y + t.c.y) / (3.0 * U);
               if (std::hypot(x - cx, y - cy) > cr)
                  continue;
               auto &c = cells.try_emplace({ (int)std::floor(x / 4.0), (int)std::floor(y / 4.0) }, std::array<double, 3>{ 0.0, 1e9, -1e9 }).first->second;
               c[0] += 1.0;
               c[1] = std::min(c[1], std::min({ t.a.z, t.b.z, t.c.z }) / U);
               c[2] = std::max(c[2], std::max({ t.a.z, t.b.z, t.c.z }) / U);
            }
            for (const RubberContact::SurfaceTri &t : m_surface->tris)   // DIAG: the vertices (VPX units) of every face whose bounding box reaches the circle
               if (t.lo.x / U <= cx + cr && t.hi.x / U >= cx - cr && t.lo.y / U <= cy + cr && t.hi.y / U >= cy - cr)
                  PLOGI << "RubberNearTri " << m_name << " a " << t.a.x / U << ' ' << t.a.y / U << ' ' << t.a.z / U << " b " << t.b.x / U << ' ' << t.b.y / U << ' ' << t.b.z / U << " c " << t.c.x / U << ' ' << t.c.y / U << ' ' << t.c.z / U;
            for (const auto &kv : cells)
               PLOGI << "RubberNear " << m_name << " cell " << kv.first.first * 4 << ' ' << kv.first.second * 4 << " faces " << kv.second[0] << " z " << kv.second[1] << ' ' << kv.second[2];
         }
      }
      if (ContactTraceOn() && m_set.contact)
         PLOGI << "RubberGeom " << m_name << " " << m_surfaceKind << " x " << m_surface->lo[0] << " " << m_surface->hi[0] << " y " << m_surface->lo[1] << " " << m_surface->hi[1] << " z " << m_surface->lo[2] << " " << m_surface->hi[2];
      return;
   }
   if (m_pts.empty())
      return;
   float x0 = FLT_MAX, x1 = -FLT_MAX, y0 = FLT_MAX, y1 = -FLT_MAX, z0 = FLT_MAX, z1 = -FLT_MAX;
   for (const Vertex3Ds &p : m_pts)
   {
      x0 = std::min(x0, p.x); x1 = std::max(x1, p.x);
      y0 = std::min(y0, p.y); y1 = std::max(y1, p.y);
      z0 = std::min(z0, p.z); z1 = std::max(z1, p.z);
   }
   m_hitBBox = FRect3D(x0 - m_radius, x1 + m_radius, y0 - m_radius, y1 + m_radius, z0 - m_radius, z1 + m_radius);
   if (ContactTraceOn() && m_set.contact)
      PLOGI << "RubberGeom " << m_name << " ring x " << x0 << " " << x1 << " y " << y0 << " " << y1 << " z " << z0 << " " << z1 << " radius " << m_radius << " points " << m_pts.size();
}

float HitRubber::HitTest(const BallS &ball, const float dtime, CollisionEvent &coll) const
{
   if (!m_enabled || m_pts.size() < 2)
      return -1.0f;
   if (CordModel())
      return -1.0f;   // the cord never meets the stock collision: the contact system takes the ball before it touches (CordNear; GroupStepDesign 15)
   if (BandIsEngaged(ball) || ContactOwns(ball))
      return -1.0f;   // the band or contact model owns this ball until it is released

   using RubberContact::V3;
   // a sling's band is tested at the shape it has now (the plate may have pushed it far past its rest line, and a ball arriving then must meet it there);
   // every other rubber is tested at its rest shape
   const bool useDyn = m_slingPlate != nullptr && m_dyn.active && m_dyn.pts.size() >= 3 && !BandDiverged();
   const size_t np = useDyn ? m_dyn.pts.size() : m_pts.size();
   auto pt = [&](size_t k) -> V3 {
      if (useDyn)
         return V3(m_dyn.pts[k].x / U, m_dyn.pts[k].y / U, m_dyn.pts[k].z / U);
      return V3(m_pts[k].x, m_pts[k].y, m_pts[k].z);
   };
   const double rt = (double)ball.m_radius + (double)m_radius;      // ball centre to centreline at touch, VPX units
   const V3 c(ball.m_pos.x, ball.m_pos.y, ball.m_pos.z);
   const V3 v(ball.m_vel.x, ball.m_vel.y, ball.m_vel.z);

   double bestT = -1.0;
   bool bestContact = false;
   double bestBnv = 0.0, bestBnd = 0.0;

   for (size_t i = 0; i + 1 < np; i++)
   {
      const V3 A = pt(i), B = pt(i + 1);
      const V3 seg = B - A;
      const double L2 = seg.Dot(seg);
      if (L2 <= 1e-12)
         continue;
      const double L = std::sqrt(L2);
      const V3 u = seg / L;

      // distance now, to the nearest point of this segment
      const double s0 = std::min(std::max((c - A).Dot(seg) / L2, 0.0), 1.0);
      const V3 P0 = A + seg * s0;
      const double d0 = (c - P0).Length();
      if (d0 <= 1.0e-6)
         continue;
      const double bnv = (c - P0).Dot(v) / d0;                    // ball normal velocity
      if (bnv > (double)C_CONTACTVEL)
         continue;                                                // clearly receding
      const double bnd = d0 - rt;                                 // ball distance to the rubber surface

      double t;
      bool isContact = false;
      if (bnd < (double)PHYS_TOUCH)                               // already touching
      {
         if (std::fabs(bnv) <= (double)C_CONTACTVEL)
         {
            isContact = true;
            t = 0.0;
         }
         else
            t = std::max(0.0, -bnd / bnv);
      }
      else
      {
         // swept test of the ball centre against the capsule of radius rt around the segment: side of the cylinder, then both end caps
         t = -1.0;
         const V3 w = c - A;
         const V3 wp = w - u * w.Dot(u), vp = v - u * v.Dot(u);
         const double a = vp.Dot(vp);
         if (a > 1e-18)
         {
            const double b = wp.Dot(vp), cc = wp.Dot(wp) - rt * rt, disc = b * b - a * cc;
            if (disc >= 0.0)
            {
               const double tc = (-b - std::sqrt(disc)) / a;
               if (tc >= 0.0)
               {
                  const double s = ((w + v * tc).Dot(u)) / L;
                  if (s >= 0.0 && s <= 1.0)
                     t = tc;
               }
            }
         }
         for (const V3 *q : { &A, &B })
         {
            const double te = SphereTime(c, v, *q, rt);
            if (te >= 0.0 && (t < 0.0 || te < t))
               t = te;
         }
         if (t < 0.0)
            continue;
      }
      if (!std::isfinite(t) || t > (double)dtime)
         continue;
      if (bestT < 0.0 || t < bestT)
      {
         bestT = t;
         bestContact = isContact;
         bestBnv = bnv;
         bestBnd = bnd;
      }
   }
   if (bestT < 0.0)
      return -1.0f;

   // normal at the hit position: from the nearest centreline point
   const V3 hp = c + v * bestT;
   double bd = 1e30;
   V3 bp;
   for (size_t i = 0; i + 1 < np; i++)
   {
      const V3 A = pt(i), B = pt(i + 1);
      const V3 seg = B - A;
      const double L2 = seg.Dot(seg);
      if (L2 <= 1e-12)
         continue;
      const double s = std::min(std::max((hp - A).Dot(seg) / L2, 0.0), 1.0);
      const V3 P = A + seg * s;
      const double d = (hp - P).Length();
      if (d < bd)
      {
         bd = d;
         bp = P;
      }
   }
   if (bd <= 1e-9)
      return -1.0f;
   const V3 n = (hp - bp) / bd;
   coll.m_hitnormal = Vertex3Ds((float)n.x, (float)n.y, (float)n.z);
   coll.m_isContact = bestContact;
   if (bestContact)
      coll.m_hit_org_normalvelocity = (float)bestBnv;
   coll.m_hitdistance = (float)bestBnd;
   return (float)bestT;
}

void HitRubber::Reload() const
{
   RubberObjectSettings s;
   if (m_surface != nullptr && m_surfaceKind == "playfield")
   {
      if (!LoadPlayfieldSettings(s, m_ini, m_tableMu))
         s.contact = false;
      m_set = s;
      PLOGI << "Playfield settings reloaded from " << m_ini.string() << ": contact " << (s.contact ? 1 : 0) << ", e_mpa " << s.eMpa << ", alpha " << s.alphaSPerM << ", tan_delta " << s.tanDelta << ", mu " << s.mu;
      return;
   }
   if (m_surface != nullptr)
   {
      if (m_hardSurface)
         LoadHardSurfaceSettings(s, m_ini, m_surfaceKind.c_str(), m_name, m_partE, m_partMu);
      else
         LoadSurfaceSettings(s, m_ini, m_surfaceKind.c_str(), m_name);   // contact follows the object's  rubber  flag, so it can be switched off while running
      m_set = s;
      PLOGI << "Rubber surface settings for " << m_name << " reloaded from " << m_ini.string() << ": rubber " << (s.contact ? 1 : 0) << ", alpha " << s.alphaSPerM << ", mu " << s.mu;
      return;
   }
   LoadSettings(s, m_ini, m_name);
   if (m_bandSpans.empty())   // a rubber with no straight span stays a ring whatever a reloaded section says
   {
      s.band = false;
      s.slingPlate.clear();
   }
   m_set = s;   // contact on/off is decided when the table is loaded; the values apply immediately
   m_set.contact = true;
   SlingPlace();
   if (m_plateConfigured)
      const_cast<HitRubber *>(this)->m_plate.Configure(m_set.sling, false);
   PLOGI << "Rubber settings for " << m_name << " reloaded from " << m_ini.string() << ": alpha " << s.alphaSPerM << ", mu " << s.mu << ", tan delta " << s.tanDelta;
}

static constexpr float kFlipperHorizonVpt = 1.0f;   // trace only: VPX time units (1 VPT = 10 ms) a bat is looked ahead for in the flip_clear numbers

void HitRubber::Collide(const CollisionEvent &coll)
{
   HitBall *const pball = coll.m_ball;
   const Vertex3Ds &n = coll.m_hitnormal;
   const float dot = -(n.Dot(pball->m_d.m_vel));
   if (dot < C_LOWNORMVEL)
      return;   // receding or sliding: nothing to solve (a resting ball is handled by Contact)

   ReloadIfChanged();
   const bool bandTook = m_set.band && BandEngage(pball);   // checks for a flipper itself (band_flipper_margin_u)
   if (m_set.band && !bandTook)
      BandNote(pball);
   if (bandTook || (!m_set.band && ContactEngage(pball)))
      g_pplayer->m_physics->OnBallWallHit(*pball, n, dot);   // the force-model contact is stepped by the engine from here on (no impulse here)
   else
   {
      // the old node band that could not take the ball (BandNote says why): the stock impulse, until the band work closes it. A post or ring never gets here.
      if (ContactTraceOn())
         PLOGI << "StockImpulse " << m_name << " t=" << g_pplayer->m_time_msec << " held " << (pball->m_contact != nullptr ? 1 : 0) << " owns " << (ContactOwns(pball->m_d) ? 1 : 0) << " dot " << dot;
      pball->Collide3DWall(n, m_elasticity, m_elasticityFalloff, m_friction, m_scatter);
   }

   if (m_obj && m_fe && dot >= m_threshold)
   {
      m_obj->m_currentHitThreshold = dot;
      FireHitEvent(pball);
   }
}

HitObject *HitRubber::HitEventObject() const
{
   if (m_surface != nullptr)
      return m_eventCarrier;
   return m_set.band ? nullptr : const_cast<HitRubber *>(this);   // a rubber post / ring fires its own hit event (Collide above)
}

void HitRubber::SideGather(const BallS &ball, const RubberContact::V3 &c, const double ballR, const double reachM, RubberSideSet &set)
{
   static const bool noSide = ContactSystem::EnvVar("VPX_CONTACT_NOSIDE", "VPX_RING_NOSIDE") != nullptr;   // diagnostics: the flipper step without the rubbers around the ball (stock contact for them)
   if (noSide)
   {
      set.items.clear();
      set.depthM.clear();
      return;
   }
   std::vector<SurfaceContact> items;
   const float reach = ball.m_radius + 20.0f;   // VPX units, as for a engagement
   for (ContactHost *const o : Registry())
   {
      if (!o->m_enabled || !o->ContactSettings().contact)
         continue;
      const FRect3D &bb = o->m_hitBBox;
      if (ball.m_pos.x < bb.left - reach || ball.m_pos.x > bb.right + reach || ball.m_pos.y < bb.top - reach || ball.m_pos.y > bb.bottom + reach
          || ball.m_pos.z < bb.zlow - reach || ball.m_pos.z > bb.zhigh + reach)
         continue;
      o->ReloadIfChanged();
      SurfaceContact it;
      bool known = false;
      for (const SurfaceContact &old : set.items)
         if (old.r == o)
         {
            it = old;   // the shear state stays
            known = true;
            break;
         }
      it.r = o;
      it.rp = o->ContactSettings().ToParams();
      it.hLayer = o->LayerM();
      if (!known)
         it.st.Reset();
      o->Refresh(it, ballR, reachM, c);
      items.push_back(std::move(it));
   }
   set.items.swap(items);
   set.depthM.assign(set.items.size(), 0.0);
}

void HitRubber::SideStep(RubberSideSet &set, const RubberContact::BallSI &bs, const RubberContact::V3 &c, const RubberContact::V3 &v, const RubberContact::V3 &w, const double h,
                         RubberContact::V3 &F, RubberContact::V3 &Mb)
{
   using namespace RubberContact;
   for (size_t i = 0; i < set.items.size(); i++)
   {
      SurfaceContact &it = set.items[i];
      const TubeProbe pr = it.r->Probe(bs, it, c);
      const ContactResult r = EvaluateTube(it.rp, bs, pr, it.hLayer, v, w, h, it.st);
      if (r.inContact)
      {
         F = F + r.force;
         Mb = Mb + r.torqueOnBall;
         set.depthM[i] = std::max(set.depthM[i], r.delta);
      }
   }
}

double HitRubber::LayerM() const
{
   if (m_set.layerMm > 0.0 || m_surface != nullptr)
      return ContactHost::LayerM();
   return 2.0 * m_radiusSI;
}

void HitRubber::Refresh(SurfaceContact &it, const double ballR, const double reachM, const RubberContact::V3 &c) const
{
   if (m_surface == nullptr)
      return;
   // A list built with some margin stays good until the ball has moved farther than the margin allows: scanning every face of a playfield mesh for every ball on every
   // call was the main cost of a table with several balls on it.
   const double need = ballR + reachM;
   if (!(it.candR >= 0.0 && (c - it.candC).Length() + need <= it.candR))
   {
      constexpr double kMargin = 0.008;   // m
      RubberContact::SurfaceCandidates(m_surface->tris, c, need + kMargin, it.candV);
      it.candC = c;
      it.candR = need + kMargin;
   }
   // The faces the ball can touch in this cycle: those whose box is within `need` of the ball (a face farther than the ball's radius plus its reach cannot touch it, and the nearest-face
   // query only has to find faces that can). The probes of the cycle's substeps test this short list, not the whole Verlet list. A surface with carved holes keeps the whole list (the hole moves the nearest point).
   const RubberContact::CarvedHole *holes = nullptr;
   int nHoles = 0;
   if (m_surfaceKind == "playfield")
      Turntables::Holes(holes, nHoles);
   it.cand.clear();
   it.tightC = c;
   if (nHoles > 0)
   {
      it.cand = it.candV;
      return;
   }
   it.tightC = c;
   const double need2 = (need + SurfaceContact::kTightSkin) * (need + SurfaceContact::kTightSkin);
   for (const int i : it.candV)
   {
      const RubberContact::SurfaceTri &t = m_surface->tris[i];
      const double dx = std::max(std::max(t.lo.x - c.x, c.x - t.hi.x), 0.0), dy = std::max(std::max(t.lo.y - c.y, c.y - t.hi.y), 0.0), dz = std::max(std::max(t.lo.z - c.z, c.z - t.hi.z), 0.0);
      if (dx * dx + dy * dy + dz * dz <= need2)
         it.cand.push_back(i);
   }
}

RubberContact::TubeProbe HitRubber::Probe(const RubberContact::BallSI &bs, const SurfaceContact &it, const RubberContact::V3 &c) const
{
   if (m_surface != nullptr)
   {
      if (m_surfaceKind == "playfield")
      {
         const RubberContact::CarvedHole *holes = nullptr;
         int nHoles = 0;
         Turntables::Holes(holes, nHoles);   // the playfield has a round hole under every turntable disc
         return RubberContact::ProbeTriangles(bs, m_surface->tris, it.Faces(c), m_set.edgeRadiusMm * 1e-3, c, holes, nHoles);
      }
      return RubberContact::ProbeTriangles(bs, m_surface->tris, it.Faces(c), m_set.edgeRadiusMm * 1e-3, c);   // (no carved hole: a ball inside a non-legacy kicker's circle is handed to the stock contact instead, see InPlayfieldHole; a carved hole let balls fall through a plane playfield)
   }
   return RubberContact::ProbeTube(bs, m_ptsSI.data(), (int)m_ptsSI.size(), m_radiusSI, c, &m_tubePre);
}

double HitRubber::ArrivalS(const RubberContact::BallSI &bs, const SurfaceContact &it, const int nearestSeg, const RubberContact::V3 &c, const RubberContact::V3 &v, const double a, const double hMax, const double touchM, const double arriveM) const
{
   if (m_surface != nullptr)
      return RubberContact::ArrivalTriangles(bs, m_surface->tris, it.Faces(c), nearestSeg, c, v, a, hMax, touchM, arriveM);
   return RubberContact::ArrivalTube(bs, m_ptsSI.data(), (int)m_ptsSI.size(), m_radiusSI, m_tubePre, nearestSeg, c, v, a, hMax, touchM, arriveM);
}

// ---- Engine-stepped contact ----
//
// A ball hitting a static rubber (a post ring, a wall rubber) used to be solved to the end inside the one Collide() call: 900 - 1500 substeps of 5 us
// (400 - 600 us, half a physics cycle, for every hit; two balls in the same cycle would not fit). Now the ball is engaged: every physics cycle RingAdvance()
// advances the same force model by that cycle's substeps (about 80 us), and the ball is handed back to the engine when it is clear of every rubber and moving away.
// Only the rubber's own contact force and torque act in the substeps; gravity, the playfield and every other object are left to the engine (as for a band).
// Every force-model rubber within reach of the ball is part of the contact, each with its own parameters and shear state, so a ball between two rubbers feels both.

bool HitRubberContactHolds(const HitBall *ball);

std::vector<BallPairState> g_ballPairs;

// A ball the pair model can take part in: not a flipper's, not in a kicker, not held by a band; a contact-held ball or a free one.
bool PairEligible(const HitBall *x)
{
   // Every ball takes part in pair contact except one held by a kicker (a state of its own, until the kicker is redesigned). Unified: the bat and the band are stepped
   // with the island, so their balls are members like any other; the legacy engine steps them apart and keeps its balls out.
   return !x->m_d.m_lockedInKicker && (ContactSystem::Unified() || (x->m_rubberOwner == nullptr && x->m_bandOwner == nullptr));
}

bool HitRubber::ContactPairModel(const HitBall *a, const BallS &b)
{
   if (!ContactSettings::Get().ballContact || !PairEligible(a))
      return false;
   const HitBall *hb = nullptr;
   for (const auto *const ball : g_pplayer->m_vball)
      if (&ball->m_hitBall.m_d == &b)
         hb = &ball->m_hitBall;
   // Unified: the contact system owns every pair (two free balls too). The legacy engine only solves pairs with a held ball.
   return hb != nullptr && PairEligible(hb) && (ContactSystem::Unified() || HitRubberContactHolds(a) || HitRubberContactHolds(hb));
}

bool HitRubberContactPairModel(const HitBall *a, const BallS &b)
{
   return HitRubber::ContactPairModel(a, b);
}

// ---- Engine-stepped band contact ----
//
// A ball hitting a free span of a band cannot be solved in one Collide() call: the dwell is 40 - 150 ms (a ball rolling slowly up the slope can rest on
// the band), during which the rest of the table acts on the ball. So the ball is engaged: every physics cycle BandAdvance() integrates the force model for
// the ball and the band's node in small SI steps (like the flipper's rubber model) and moves the ball for the whole cycle; gravity, playfield and every
// other object are left to the engine. The band's state is kept in m_dyn, so after the ball leaves the band rings down on its own.

bool BandMover::AddToList() const
{
   return (m_owner->m_set.band && !m_owner->m_bandSpans.empty()) || m_owner->m_slingPlate != nullptr || m_owner->m_set.contact;
}

void BandMover::PreUpdateDisplacements(const float dtime)
{
   m_owner->BandAdvance(dtime);
}

bool HitRubberContactHolds(const HitBall *ball)
{
   return ContactSystem::Holds(ball);
}

void HitRubberContactRelease(HitBall *ball, const char *why)
{
   ContactSystem::Release(ball, why);
}

void HitRubberBandBallDeleted(HitBall *ball)
{
   if (ball->m_bandOwner != nullptr)
      ball->m_bandOwner->BandRelease(ball, "ball deleted");
   ContactSystem::Release(ball, "ball deleted");
}

bool HitRubber::BandIsEngaged(const BallS &ball) const
{
   for (const BandBall &e : m_engaged)
      if (&e.ball->m_d == &ball)
         return true;
   return false;
}

bool HitRubber::BandWants(const HitBall *ball) const
{
   if (!m_set.band || m_bandSpans.empty() || ball->m_bandOwner != nullptr || ball->m_rubberOwner != nullptr)   // a contact (the playfield's, above all) is handed over to the band
      return false;
   if (CordModel())
      return true;   // no cooling down and no giving way to a flipper: there is no stock collision to leave the ball to
   return g_pplayer->m_time_msec >= m_bandCooldownUntilMs && !g_pplayer->m_physics->FlipperNear(ball->m_d, (float)m_set.bandFlipperMarginU);
}

// The deformed centreline: the ring with the span of every chain that has nodes replaced by A - P_1 - ... - P_M - B. Rebuilt whenever a node is added.
void HitRubber::BandRebuild()
{
   const int n = (int)m_pts.size() - 1;
   BandDyn &d = m_dyn;
   std::vector<char> removed(n, 0);
   std::vector<int> insertAfter(n, -1);
   for (size_t c = 0; c < d.chains.size(); c++)
      if (d.chains[c].Count() > 0)
      {
         const BandSpan &sp = m_bandSpans[d.span[c]];
         for (int m = sp.first + 1; m < sp.last; m++)
            removed[m % n] = 1;
         insertAfter[sp.first % n] = (int)c;
      }
   d.pts.clear();
   d.base.assign(d.chains.size(), -1);
   for (int i = 0; i < n; i++)
      if (!removed[i])
      {
         d.pts.push_back(m_ptsSI[i]);
         if (insertAfter[i] >= 0)
         {
            const int c = insertAfter[i];
            d.base[c] = (int)d.pts.size() - 1;   // this point is the chain's A; its nodes follow, then B
            for (int k = 1; k <= d.chains[c].Count(); k++)
               d.pts.push_back(d.chains[c].PtPos(k));
         }
      }
   d.pts.push_back(d.pts[0]);   // closes the loop
}

// Copy the nodes' positions into the deformed centreline
void HitRubber::BandSyncPts()
{
   BandDyn &d = m_dyn;
   for (size_t c = 0; c < d.chains.size(); c++)
      for (int k = 1; k <= d.chains[c].Count(); k++)
         d.pts[d.base[c] + k] = d.chains[c].PtPos(k);
}

// Every free span of the band is a chain of nodes. The span under the ball gets a node at the ball's place along it (or the ball takes a node that is
// already within 8 mm of its place, still ringing from an earlier hit); the other free spans get a passive node at the point nearest to it, so the
// loaded span can push on them (a loop's back span). Several balls can be on the band at once, on different spans or on the same one (each loads its own
// node); a ball is refused when another ball already loads the node it would take, or the span has too many nodes. chain, node = what the ball loads.
// The deformable band: one chain per free span, none of them with nodes yet.
bool HitRubber::BandEnsureDyn()
{
   using namespace RubberContact;
   if (m_dyn.active)
      return true;
   const int n = (int)m_pts.size() - 1;
   BandDyn d;
   d.bp.eMpa = m_set.eMpa;               // the same rubber as the rings: one modulus and one loss tangent
   d.bp.prestretch = m_set.bandPrestretchPct * 0.01;
   d.bp.tensionN = m_set.bandTensionN;
   d.bp.tanDelta = m_set.tanDelta;
   d.bp.refHz = m_set.bandRefHz;
   d.bp.nodeMassKg = m_set.bandNodeMassG * 1e-3;
   d.bp.lockZ = m_set.bandLockZ;
   d.bp.zeta = m_set.bandNodeZeta;
   {
      const std::string &rest = m_set.bandStops;
      size_t pos = 0;
      while (pos < rest.size())
      {
         size_t bar = rest.find('|', pos);
         if (bar == std::string::npos)
            bar = rest.size();
         const std::string item = rest.substr(pos, bar - pos);
         double a = 0, b = 0, c = 0;
         if (sscanf(item.c_str(), " %lf , %lf , %lf", &a, &b, &c) == 3 && c > 0.0)
            d.bp.stops.push_back({ a * U, b * U, c * U + m_radiusSI });
         else if (!Trim(item).empty())
         {
            // a table object by name, optionally with a radius: "Post7" or "Post7, 9"
            std::string nm = Trim(item);
            bool hasR = false;
            double rIn = 0.0;
            const size_t comma = nm.find(',');
            if (comma != std::string::npos)
            {
               hasR = sscanf(nm.c_str() + comma + 1, " %lf", &rIn) == 1;
               nm = Trim(nm.substr(0, comma));
            }
            double x = 0, y = 0, r = 0;
            if (m_stopResolver && m_stopResolver(nm, hasR, rIn, x, y, r) && r > 0.0)
            {
               d.bp.stops.push_back({ x * U, y * U, r * U + m_radiusSI });
               if (m_stopWarned.insert("ok:" + nm + ":" + std::to_string((int)(x * 10.0)) + "," + std::to_string((int)(y * 10.0))).second)
                  PLOGI << "Band " << m_name << ": band_stops '" << nm << "' is at (" << x << ", " << y << ") VPX units, radius " << r << " U";
            }
            else if (m_stopWarned.insert(nm).second)
               PLOGW << "Band " << m_name << ": band_stops entry '" << Trim(item) << "' is neither x,y,r nor a table object that has a radius (rubber, bumper, kicker; any other part needs 'name, r')";
         }
         pos = bar + 1;
      }
   }
   for (size_t k = 0; k < m_bandSpans.size(); k++)
   {
      const BandSpan &sp = m_bandSpans[k];
      if (!sp.free || sp.last - sp.first < 1 || sp.last - sp.first >= n)
         continue;
      BandChain ch;
      ch.Init(m_ptsSI[sp.first % n], m_ptsSI[sp.last % n], m_radiusSI, d.bp);
      ch.PinStops();   // a post the band lies against becomes a fixed corner of the cord
      if (m_set.bandNodeMm > 0.0 && m_slingPlate == nullptr)
         ch.FillNodes(m_set.bandNodeMm * 1e-3);
      if (ch.stops.size() < d.bp.stops.size())
         PLOGI << "Band " << m_name << ": span " << k << " ignores " << d.bp.stops.size() - ch.stops.size() << " band_stops entr" << (d.bp.stops.size() - ch.stops.size() == 1 ? "y" : "ies")
               << " (it contains an anchor of the span)";
      d.chains.push_back(ch);
      d.span.push_back((int)k);
   }
   if (d.chains.empty())
      return false;
   d.active = true;
   m_dyn = d;
   BandRebuild();
   return true;
}

bool HitRubber::BandBuild(const RubberContact::V3 &c, const double ballR, int &chain, int &node, const double sameSpot)
{
   using namespace RubberContact;
   constexpr int kMaxNodes = 12;               // per span
   constexpr double kCompanionSpot = 0.010;    // m: a passive node on the other spans is not added if one is this close
   const int n = (int)m_pts.size() - 1;
   int best = -1;
   double bestD = 1e30, bestF = 0.5;
   for (size_t k = 0; k < m_bandSpans.size(); k++)
   {
      const BandSpan &sp = m_bandSpans[k];
      if (!sp.free || sp.last - sp.first < 1 || sp.last - sp.first >= n)
         continue;
      const V3 &a = m_ptsSI[sp.first % n], &b = m_ptsSI[sp.last % n];
      const V3 ab = b - a;
      const double l2 = ab.Dot(ab);
      if (l2 < 1e-12)
         continue;
      const double f = std::min(std::max((c - a).Dot(ab) / l2, 0.0), 1.0);
      const double d = (c - (a + ab * f)).Length();
      if (d < bestD)
      {
         bestD = d;
         best = (int)k;
         bestF = f;
      }
   }
   if (best < 0 || bestD > ballR + m_radiusSI + 0.03)
      return false;   // the ball is not near the middle of a free span (a post wrap or a short span): that part is rigid
   double f = std::min(std::max(bestF, 0.05), 0.95);

   if (!m_dyn.active && !BandEnsureDyn())
      return false;
   int ci = -1;
   for (size_t k = 0; k < m_dyn.span.size(); k++)
      if (m_dyn.span[k] == best)
         ci = (int)k;
   if (ci < 0)
      return false;
   BandChain &ch = m_dyn.chains[ci];
   const double Ls = ch.Ls;
   if (ch.Count() > 0)
   {
      // The ball's place on the cord as it is now: the material fraction of the point of the deformed cord nearest to the ball. The fraction along the rest line
      // is only the same while the cord is hardly deflected; under a heavy ball (a deep V) it is not, and a node given the rest-line fraction would get an
      // unstretched length on one side that is far too short (an unnaturally stiff leg that pulls the ball toward the nearer post).
      double bestDist = 1e30, fProj = f;
      for (int j = 0; j <= ch.Count(); j++)
      {
         const V3 p = ch.PtPos(j), q = ch.PtPos(j + 1);
         const double ex = q.x - p.x, ey = q.y - p.y, l2 = ex * ex + ey * ey;
         if (l2 < 1e-12)
            continue;
         const double u = std::min(std::max(((c.x - p.x) * ex + (c.y - p.y) * ey) / l2, 0.0), 1.0);
         const double dd = std::hypot(c.x - (p.x + ex * u), c.y - (p.y + ey * u));
         if (dd < bestDist)
         {
            bestDist = dd;
            fProj = ch.PtF(j) + u * (ch.PtF(j + 1) - ch.PtF(j));
         }
      }
      f = std::min(std::max(fProj, 0.05), 0.95);
   }

   if (BandFixedNodes() && ch.Count() > 0)
   {
      // permanent nodes: the ball is tied to the nearest one (only for its start-up bookkeeping; the contact loads the legs under it)
      int nn = ch.order[0];
      for (int id : ch.order)
         if (std::fabs(ch.nodes[id].f - f) < std::fabs(ch.nodes[nn].f - f))
            nn = id;
      chain = ci;
      node = nn;
      return true;
   }
   int nid = -1;
   for (int id : ch.order)
      if (std::fabs(ch.nodes[id].f - f) * Ls < sameSpot)   // a ball this close to a node loads that node
         nid = id;
   if (nid < 0 && m_slingPlate != nullptr && ci == SlingChain())
   {
      // under the face of a sling's plate the ball loads the nearest of the plate's own nodes: a node of its own between the two at the ends of the face is
      // about 6 mm from both, and that pair of stiff short legs under the pushing plate went unstable (the plate's angle ran away to infinity, the sling stopped)
      const double tol = m_set.bandFollowMm * 1e-3 / std::max(Ls, 1e-6);
      if (f >= m_plateF0 - tol && f <= m_plateF1 + tol)
      {
         double bd = 1e30;
         for (int id : ch.order)
            if (SlingInFootprint(ci, id) && std::fabs(ch.nodes[id].f - f) < bd)
            {
               bd = std::fabs(ch.nodes[id].f - f);
               nid = id;
            }
      }
   }
   bool changed = false;
   if (nid >= 0)
   {
      for (const BandBall &e : m_engaged)
         if (e.chain == ci && e.node == nid)
            return false;   // another ball already loads this node
   }
   else
   {
      if (ch.Count() >= kMaxNodes)
         return false;
      nid = ch.Add(f);
      changed = true;
   }
   const V3 P0 = ch.nodes[nid].P0;
   for (size_t cj = 0; cj < m_dyn.chains.size(); cj++)
   {
      if ((int)cj == ci)
         continue;
      BandChain &oc = m_dyn.chains[cj];
      const V3 ab = oc.B - oc.A;
      const double l2 = ab.Dot(ab);
      if (l2 < 1e-12 || oc.Count() >= kMaxNodes)
         continue;
      const double fo = std::min(std::max((P0 - oc.A).Dot(ab) / l2, 0.05), 0.95);
      bool close = false;
      for (int id : oc.order)
         close = close || std::fabs(oc.nodes[id].f - fo) * oc.Ls < kCompanionSpot;
      if (!close)
      {
         oc.Add(fo);
         changed = true;
      }
   }
   if (changed)
      BandRebuild();
   chain = ci;
   node = nid;
   return true;
}

bool HitRubber::BandEngage(HitBall *const pb)
{
   using namespace RubberContact;
   if (!BandWants(pb))
      return false;
   CordEnsure();
   if (CordOn())   // the group step takes it from here: the cord finds the ball itself
   {
      if (!ContactSystem::EngageBand(this, pb))
         return false;
      BandBall e;
      e.ball = pb;
      e.chain = -1;
      e.node = -1;
      e.vin = V3((double)pb->m_d.m_vel.x * VS, (double)pb->m_d.m_vel.y * VS, (double)pb->m_d.m_vel.z * VS);
      m_engaged.push_back(e);
      pb->m_bandOwner = this;
      return true;
   }
   BallSI bs;
   bs.R = (double)pb->m_d.m_radius * U;
   const Vertex3Ds p = pb->m_d.m_pos;
   V3 c((double)p.x * U, (double)p.y * U, (double)p.z * U);
   V3 v((double)pb->m_d.m_vel.x * VS, (double)pb->m_d.m_vel.y * VS, (double)pb->m_d.m_vel.z * VS);
   int chain = 0, node = 0;
   if (!BandBuild(c, bs.R, chain, node))
      return false;

   // The ball starts the engagement touching or a little inside the band (the hit was found within a step). Go back to the instant it first touched:
   // along its own velocity, by depth / closing speed, at most one physics step; if it is not closing, seat it at touching distance along the normal.
   // Nothing is stored as compression.
   constexpr double kMaxRewindS = PHYSICS_STEPTIME_S;
   double rewound = 0.0, delta0 = 0.0;
   const V3 c0 = c;
   for (int pass = 0; pass < 3; pass++)
   {
      BandSyncPts();
      const TubeProbe p0 = ProbeTube(bs, m_dyn.pts.data(), (int)m_dyn.pts.size(), m_radiusSI, c);
      if (!p0.overlap || p0.delta <= 1e-6)
         break;
      delta0 = std::max(delta0, p0.delta);
      const V3 vrel = v - m_dyn.chains[chain].nodes[node].vp;
      const double closing = -vrel.Dot(p0.normal);
      if (closing > 0.01)
      {
         const double t = std::min(p0.delta / closing, kMaxRewindS);
         c = c - v * t;
         rewound = std::max(rewound, t);
      }
      else
         c = c + p0.normal * p0.delta;
      ContactSystem::SeatKeepsOthers(pb, this, bs, c0, c);
   }
   pb->m_d.m_pos = Vertex3Ds((float)(c.x / U), (float)(c.y / U), (float)(c.z / U));

   if (ContactSystem::Unified())
   {
      if (!ContactSystem::EngageBand(this, pb))   // the ball's contact (the playfield, posts) stays and the band joins it
         return false;
   }
   else
      HitRubberContactRelease(pb, "band");   // the playfield's (or a post's) force-model contact makes way for the band
   if (m_engaged.empty())
      m_dyn.peakFn = m_dyn.peakBB = 0.0;
   BandBall e;
   e.ball = pb;
   e.chain = chain;
   e.node = node;
   e.st.Reset();
   e.vin = v;
   {
      const BandChain &sc = m_dyn.chains[chain];
      const V3 ab = sc.B - sc.A;
      const double s = (c - sc.A).Dot(V3(ab.y, -ab.x, 0.0));
      e.side = s >= 0.0 ? 1.0 : -1.0;
   }
   m_engaged.push_back(e);
   pb->m_bandOwner = this;
   const BandChain &ch = m_dyn.chains[chain];
   PLOGI << "Band " << m_name << ": engaged ball at t=" << g_pplayer->m_time_msec << " ms, span " << m_dyn.span[chain] << " at " << ch.nodes[node].f << " of its length, speed "
         << v.Length() << " m/s, start overlap " << delta0 * 1e3 << " mm, rewound " << rewound * 1e3 << " ms, " << ch.Count() << " nodes on the span of " << ch.nodes[node].m * 1e3
         << " g, T0 " << ch.EA * m_set.bandPrestretchPct * 0.01 << " N, balls on the band " << m_engaged.size();
   return true;
}

void HitRubber::BandRelease(HitBall *const ball, const char *why)
{
   size_t i = 0;
   while (i < m_engaged.size() && m_engaged[i].ball != ball)
      i++;
   if (i == m_engaged.size())
      return;
   const BandBall e = m_engaged[i];
   // The ball must leave in front of the rest line, on the side it came from. A fast ball on a band that a plate holds forward, or a trip of the model, can hand it back to the
   // engine already on the far side of the rest line, and then it simply passes through (it ended up inside the sling). Put it back in front of the
   // line, touching it, and give it a damped bounce off it.
   if (e.side != 0.0 && std::strcmp(why, "ball deleted") != 0 && std::strcmp(why, "band switched off") != 0
       && e.chain >= 0 && e.chain < (int)m_dyn.chains.size())
   {
      using namespace RubberContact;
      const BandChain &ch = m_dyn.chains[e.chain];
      const V3 ab = ch.B - ch.A;
      const double L2 = ab.Dot(ab);
      if (L2 > 1e-12)
      {
         const double L = std::sqrt(L2);
         const V3 n(ab.y / L * e.side, -ab.x / L * e.side, 0.0);   // from the band toward the side the ball came from
         V3 c((double)ball->m_d.m_pos.x * U, (double)ball->m_d.m_pos.y * U, (double)ball->m_d.m_pos.z * U);
         const double R = (double)ball->m_d.m_radius * U;
         const double along = (c - ch.A).Dot(ab) / L2;
         const double s = (c - ch.A).Dot(n);
         if (along > 0.0 && along < 1.0 && s < R + m_radiusSI - 0.5e-3)
         {
            c = c + n * (R + m_radiusSI - s);
            ball->m_d.m_pos = Vertex3Ds((float)(c.x / U), (float)(c.y / U), (float)(c.z / U));
            V3 v((double)ball->m_d.m_vel.x * VS, (double)ball->m_d.m_vel.y * VS, (double)ball->m_d.m_vel.z * VS);
            const double vn = v.Dot(n);
            if (vn < 0.0)
            {
               v = v - n * (vn * 1.3);   // the normal speed comes back as 0.3 of what it was
               ball->m_d.m_vel = Vertex3Ds((float)(v.x / VS), (float)(v.y / VS), (float)(v.z / VS));
            }
            PLOGW << "Band " << m_name << ": ball released (" << why << ") inside the band's rest position (its centre " << s * 1e3 << " mm in front of the rest line, less than the " << (R + m_radiusSI) * 1e3 << " mm it needs), put back touching it";
         }
      }
   }
   const Vertex3Ds &vel = ball->m_d.m_vel;
   const double speedOut = std::hypot(std::hypot(vel.x, vel.y), vel.z) * VS;
   const bool nodeOk = e.chain >= 0 && e.chain < (int)m_dyn.chains.size() && e.node >= 0 && e.node < (int)m_dyn.chains[e.chain].nodes.size();
   double peakSpan0 = 0.0, peakSpan1 = 0.0;   // the largest deflection reached on the first two spans, for the log
   for (size_t c = 0; c < m_dyn.chains.size() && c < 2; c++)
      for (const RubberContact::BandChain::Node &nd : m_dyn.chains[c].nodes)
         if (nd.alive)
            (c == 0 ? peakSpan0 : peakSpan1) = std::max(c == 0 ? peakSpan0 : peakSpan1, nd.peak);
   PLOGI << "Band " << m_name << ": released ball (" << why << ") at t=" << g_pplayer->m_time_msec << " ms after " << e.dwellS * 1e3 << " ms, peak deflection "
         << e.peak * 1e3 << " mm (spans " << peakSpan0 * 1e3 << ", " << peakSpan1 * 1e3 << "), peak force " << m_dyn.peakFn << " N, span overlap " << m_dyn.peakBB * 1e3
         << " mm at " << (nodeOk ? m_dyn.chains[e.chain].nodes[e.node].f : 0.0) << " of the span, speed in " << e.vin.Length() << " m/s, out " << speedOut << " m/s, vel out ("
         << vel.x * VS << ", " << vel.y * VS << ", " << vel.z * VS << ") m/s, balls left on the band " << m_engaged.size() - 1;
   if (ball->m_contact != nullptr && ball->m_contact->band == this)
      ball->m_contact->band = nullptr;   // (unified: the ball may stay in its other contacts)
   ball->m_bandOwner = nullptr;
   m_engaged.erase(m_engaged.begin() + (long)i);
   m_info.peakMm = (float)(e.peak * 1e3);
   if (nodeOk)
   {
      const double f = m_dyn.chains[e.chain].nodes[e.node].f;
      m_info.hitSpan = m_dyn.span[e.chain];
      m_info.hitPosSpan = (float)f;
      m_info.hitPosRubber = BandRubberPos(m_info.hitSpan, f);
   }
   m_info.hitTimeMs = (int)g_pplayer->m_time_msec;
   m_info.speedIn = (float)e.vin.Length();
   m_info.speedOut = (float)speedOut;
   BandPublish();
   if (std::strcmp(why, "ball deleted") != 0 && std::strcmp(why, "band switched off") != 0 && m_editable != nullptr)
      static_cast<Rubber *>(m_editable)->FireVoidEventParm(DISPID_RubberEvents_BandHit, m_info.peakMm);   // the script handler may do anything: nothing of ours is in flight
}

void HitRubber::BandReleaseAll(const char *why)
{
   while (!m_engaged.empty())
      BandRelease(m_engaged.back().ball, why);
}

// Refresh what scripts can read: the most deflected node's deflection, its place along its span and along the whole rubber, and its position.
void HitRubber::BandPublish()
{
   const BandDyn &d = m_dyn;
   int bc = -1, bn = -1;
   double best = -1.0;
   for (size_t c = 0; c < d.chains.size(); c++)
      for (const int k : d.chains[c].order)
         if (d.chains[c].nodes[k].Deflection() > best)
         {
            best = d.chains[c].nodes[k].Deflection();
            bc = (int)c;
            bn = k;
         }
   m_info.active = d.active;
   if (!d.active || bc < 0)
   {
      m_info.deflectionMm = 0.f;
      return;
   }
   const RubberContact::BandChain::Node &nd = d.chains[bc].nodes[bn];
   m_info.span = d.span[bc];
   m_info.deflectionMm = (float)(nd.Deflection() * 1e3);
   m_info.posSpan = (float)nd.f;
   m_info.nodeX = (float)(nd.P.x / U);
   m_info.nodeY = (float)(nd.P.y / U);
   m_info.nodeZ = (float)(nd.P.z / U);
   m_info.posRubber = BandRubberPos(d.span[bc], nd.f);
}

float HitRubber::BandRubberPos(const int span, const double f) const
{
   if (span < 0 || span >= (int)m_bandSpans.size())
      return 0.f;
   const BandSpan &sp = m_bandSpans[span];
   const int n = (int)m_pts.size() - 1;
   double len = 0.0;   // length of the span along the centreline
   for (int m = sp.first; m < sp.last; m++)
      len += (m_ptsSI[(m + 1) % n] - m_ptsSI[m % n]).Length();
   const double total = m_cumSI[n];
   return total > 0.0 ? (float)(std::fmod(m_cumSI[sp.first % n] + f * len, total) / total) : 0.f;
}

// True if a node of the deformable band has left the range the model is built for (not a number, or deflected more than 0.4 of its span). The collision test
// and the mesh must never see such a band: its legs run from a post to a point that may be metres away, and a ball near the post is thrown by them.
bool HitRubber::BandDiverged() const
{
   for (const RubberContact::BandChain &ch : m_dyn.chains)
      for (const int id : ch.order)
      {
         const RubberContact::BandChain::Node &nd = ch.nodes[id];
         if (!std::isfinite(nd.P.x + nd.P.y + nd.P.z + nd.vp.x + nd.vp.y + nd.vp.z) || !(nd.Deflection() <= 0.4 * ch.Ls))
            return true;
      }
   return false;
}

void HitRubber::BandBeforeCore(const float dtime)
{
   if (m_slingPlate != nullptr)
      SlingAdvance(dtime);
}

void HitRubber::BandAfterCore(const float dtime)
{
   if (m_slingPlate == nullptr)
   {
      if (m_dyn.active && BandDiverged())
      {
         PLOGW << "Band " << m_name << ": the band diverged, reset to rest";
         BandReleaseAll("diverged");
         m_dyn = BandDyn();
         BandPublish();
      }
      if (!m_switches.empty())
      {
         SwitchesUpdate((double)dtime * DEFAULT_STEPTIME_S);
         SwitchesFinish();
      }
      return;
   }
   if (m_dyn.active && BandDiverged())
   {
      PLOGW << "Sling " << m_name << ": the band diverged, reset to rest (the plate is reset too)";
      BandReleaseAll("diverged");
      m_dyn = BandDyn();
      BandPublish();
      m_plate.Configure(m_set.sling, true);
      m_driver.Cancel();
      m_slingIdle = false;
   }
   if (m_dyn.active)
      BandSyncPts();   // the collision test of the next cycle meets the band where it is now
   SlingFinish();
}

void HitRubber::BandAdvance(const float dtime)
{
   CordEnsure();
   if (CordOn())
   {
      if (!m_cordStepped)
         ContactSystem::StepCordAlone(this, dtime);
      m_cordStepped = false;
      CordFinish(dtime);
      return;
   }
   if (m_bandIslandStepped)   // unified engine: the contact system did the sling's driver and the band's step with its balls
   {
      m_bandIslandStepped = false;
      BandAfterCore(dtime);
      return;
   }
   BandBeforeCore(dtime);
   BandAdvanceCore(dtime);
   BandAfterCore(dtime);
}

// ---- The cord in the group step (docs/GroupStepDesign.txt 8, 11) ----
void HitRubber::CordEnsure()
{
   if (m_cord != nullptr || !m_set.band || m_set.bandModel == 0.0 || m_bandPosts.size() < 2)
      return;
   auto cb = std::make_shared<CordBand>();
   BandCord::Part &p = cb->part;
   BandCord::Settings &s = p.set;
   const double rc = m_radiusSI;
   s.cordRadius = rc;
   s.eMpa = m_set.eMpa;
   s.tanDelta = m_set.tanDelta;
   s.refHz = m_set.bandRefHz;
   s.alphaSPerM = m_set.alphaSPerM;
   s.mu = m_set.mu;
   s.ringRp = m_set.ToParams();   // a rubber on a post: the ring's own values (GroupStepDesign 13)
   s.ringRpSet = true;
   s.prestretch = m_set.bandPrestretchPct * 0.01;
   s.tensionN = s.eMpa * 1e6 * M_PI * rc * rc * s.prestretch;
   s.sling = m_set.sling;
   s.switchGapM = m_set.sling.switchGapMm * 1e-3;
   BandCord::Geometry &g = p.geo;
   const int np = std::min((int)m_bandPosts.size(), (int)BandCord::kMaxPosts);
   double area = 0.0;
   for (int i = 0; i < np; i++)
      area += m_bandPosts[i].x * m_bandPosts[(i + 1) % np].y - m_bandPosts[(i + 1) % np].x * m_bandPosts[i].y;
   g.nPosts = np;
   for (int i = 0; i < np; i++)
   {
      const BandPost &bp = m_bandPosts[area < 0.0 ? np - 1 - i : i];
      g.posts[i] = BandCord::Circle { BandCord::V2(bp.x * U, bp.y * U), bp.r * U };
   }
   // the posts behind the band (band_stops): "x,y,r" in VPX units, or a table object by name
   {
      const std::string &rest = m_set.bandStops;
      size_t pos = 0;
      while (pos < rest.size() && g.nStops < BandCord::kMaxStops)
      {
         size_t bar = rest.find('|', pos);
         if (bar == std::string::npos)
            bar = rest.size();
         const std::string item = rest.substr(pos, bar - pos);
         double a = 0, b = 0, c = 0;
         if (sscanf(item.c_str(), " %lf , %lf , %lf", &a, &b, &c) == 3 && c > 0.0)
            g.stops[g.nStops++] = BandCord::Circle { BandCord::V2(a * U, b * U), c * U + rc };
         else if (!Trim(item).empty() && m_stopResolver)
         {
            std::string nm = Trim(item);
            bool hasR = false;
            double rIn = 0.0, x = 0, y = 0, r = 0;
            const size_t comma = nm.find(',');
            if (comma != std::string::npos)
            {
               hasR = sscanf(nm.c_str() + comma + 1, " %lf", &rIn) == 1;
               nm = Trim(nm.substr(0, comma));
            }
            if (m_stopResolver(nm, hasR, rIn, x, y, r) && r > 0.0)
               g.stops[g.nStops++] = BandCord::Circle { BandCord::V2(x * U, y * U), r * U + rc };
         }
         pos = bar + 1;
      }
   }
   // the plate and its switches, as the sling's own setup measured them (SlingPlace)
   const int n = (int)m_pts.size() - 1;
   if (m_slingPlate != nullptr && m_plateGeomValid && m_slingSpan >= 0 && n >= 4)
   {
      const BandSpan &sp = m_bandSpans[m_slingSpan];
      const Vertex3Ds &A = m_pts[((sp.first % n) + n) % n], &B = m_pts[((sp.last % n) + n) % n];
      const BandCord::V2 a(A.x * U, A.y * U), bb(B.x * U, B.y * U);
      g.hasPlate = true;
      BandCord::PlateGeom &pg = g.plate;
      pg.A = a;
      pg.B = bb;
      pg.n = BandCord::V2(m_slingNx, m_slingNy);
      pg.originN = m_plateOriginN;
      pg.a = m_plateA;
      pg.h = m_plateH;
      pg.rest = m_plateRest;
      pg.f0 = m_plateF0;
      pg.f1 = m_plateF1;
   }
   g.nSwitches = 0;
   for (const BandSwitch &sw : m_switches)
      if (sw.span >= 0 && sw.span < (int)m_bandSpans.size() && g.nSwitches < BandCord::kMaxSwitches && n >= 4)
      {
         const BandSpan &sp = m_bandSpans[sw.span];
         const Vertex3Ds &A = m_pts[((sp.first % n) + n) % n], &B = m_pts[((sp.last % n) + n) % n];
         BandCord::SwitchSpec &ss = g.sw[g.nSwitches++];
         ss.P = BandCord::V2((A.x + (B.x - A.x) * sw.f) * U, (A.y + (B.y - A.y) * sw.f) * U);
         ss.b = BandCord::V2(sw.bx, sw.by);
         ss.coil = sw.coil;
      }
   p.Configure();
   if (m_set.bandRingIn > 0.0)   // the stretch at rest from the ring, the tension from the rubber
   {
      const double L = p.Length(), Lring = M_PI * (m_set.bandRingIn * 25.4e-3 + 2.0 * rc);
      s.prestretch = L / Lring - 1.0;
      if (s.prestretch <= 0.0)
      {
         PLOGW << "Band " << m_name << ": a " << m_set.bandRingIn << " in ring is longer than the loop on its posts (" << L * 1e3 << " mm): no stretch";
         s.prestretch = 0.001;
      }
      s.tensionN = s.eMpa * 1e6 * M_PI * rc * rc * s.prestretch;
      p.Configure();
   }
   double zc = 0.0;
   for (const Vertex3Ds &v : m_pts)
      zc += v.z;
   zc = m_pts.empty() ? 0.0 : zc / (double)m_pts.size();
   p.frame.O = BandCord::V3(0.0, 0.0, zc * U);
   s.playfieldWorldZ = 0.0;
   m_cord = cb;
   m_cordFired = 0;
   PLOGI << "Band " << m_name << ": cord model, " << np << " posts, " << g.nStops << " stops, " << (g.hasPlate ? "plate, " : "") << g.nSwitches << " switches, cord "
         << 2e3 * rc << " mm, E " << s.eMpa << " MPa, loop " << p.Length() * 1e3 << " mm, stretch " << s.prestretch * 100.0 << " %, tension " << s.tensionN << " N";
}

void HitRubber::CordBegin()
{
   BandCord::Part &p = m_cord->part;
   if (!m_slingEnabled)
      p.driver.Cancel();
   if (p.geo.hasPlate)
   {
      if (m_slingEnabled && m_slingTriggerReq)
         p.driver.Trigger();
      m_slingTriggerReq = false;
   }
   if (m_slingEnabled)
      p.CycleBegin();
   else
   {
      const bool keep = false;
      (void)keep;
   }
}

void HitRubber::CordFinish(const float /*dtime*/)
{
   BandCord::Part &p = m_cord->part;
   for (const auto &h : m_cord->hits)   // the script's Hit event, as Collide fires it for a stock hit
      if (m_obj && m_fe && h.second >= m_threshold)
      {
         m_obj->m_currentHitThreshold = h.second;
         FireHitEvent(h.first);
      }
   m_cord->hits.clear();
   std::vector<std::pair<int, bool>> changes;
   const int ns = std::min((int)m_switches.size(), p.geo.nSwitches);
   bool anyClosed = false;
   for (int i = 0; i < ns; i++)
   {
      BandSwitch &sw = m_switches[i];
      sw.dispM = m_cord->attached ? p.SwitchPush(i) : 0.0;
      const bool closed = p.SwitchIsClosed(i);
      if (closed != sw.closed)
         changes.emplace_back(i, closed);
      sw.closed = closed;
      anyClosed = anyClosed || (closed && sw.coil);
   }
   if (m_editable != nullptr)
      for (const auto &c : changes)
         static_cast<Rubber *>(m_editable)->FireVoidEventParm(c.second ? DISPID_RubberEvents_SwitchClosed : DISPID_RubberEvents_SwitchOpened, c.first);
   if (m_slingPlate != nullptr && p.geo.hasPlate)
   {
      m_plate.st.phi = p.phi;
      m_plate.st.omega = p.omega;
      m_slingAnyClosed = anyClosed;
      if (p.driver.fired != m_cordFired)
      {
         m_cordFired = p.driver.fired;
         m_slingInfo.fireTimeMs = (int)g_pplayer->m_time_msec;
         static_cast<Rubber *>(m_editable)->FireGroupEvent(DISPID_RubberEvents_SlingFire);
      }
      SlingFinish();
   }
   else if (!m_switches.empty())
      SwitchesFinish();
}

void HitRubber::CordSolveFailed(const char *why)
{
   PLOGW << "Band " << m_name << ": " << why;
   if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
      g_pplayer->m_liveUI->PushNotification("Band " + m_name + ": " + why, 5000);
   m_cord->attached = false;
   m_cord->slots.clear();
   BandReleaseAll("solve failed");
}

// ---- One cycle of the band, in phases ----
// BandAdvanceCore drives them for the legacy step (balls integrated here); the contact system drives the same phases for the balls it holds (unified engine), so that the
// band's nodes are stepped in the island's substeps next to the balls and their other contacts:
//   BandCycleBegin   once per cycle: the substep, the rubber, the buffers; a ball that has rolled away from its node takes a new one. False = nothing to step with balls.
//   per substep:     BandSubstepBegin, BandBallForce for each ball (the reaction goes into the nodes), BandSubstepEnd (nodes step), BandBallPeak for each ball
//   BandCycleEnd     once per cycle: idle nodes merged away, a band that left the model's range reset (the balls are handed back untouched)

void HitRubber::BandZeroForces()
{
   BandCycle &bc = m_bc;
   for (size_t c = 0; c < m_dyn.chains.size(); c++)
      bc.fNode[c].assign(m_dyn.chains[c].nodes.size(), RubberContact::V3());
}

// the spans pushing on each other (a ball presses the front span of a loop onto the back one)
void HitRubber::BandChainContacts()
{
   BandCycle &bc = m_bc;
   BandDyn &d = m_dyn;
   const size_t nc = d.chains.size();
   for (size_t i = 0; i < nc; i++)
      for (size_t j = i + 1; j < nc; j++)
         if (d.chains[i].Count() > 0 && d.chains[j].Count() > 0)
            BandChainContact(bc.rp, m_radiusSI, d.chains[i], d.chains[j], bc.fNode[i], bc.fNode[j], &d.peakBB);
}

// by node id: loaded by a ball (not to be merged away or taken over)
void HitRubber::BandMarkLoaded()
{
   BandCycle &bc = m_bc;
   BandDyn &d = m_dyn;
   const size_t nc = d.chains.size();
   bc.loaded.assign(nc, std::vector<char>());
   for (size_t c = 0; c < nc; c++)
      bc.loaded[c].assign(d.chains[c].nodes.size(), 0);
   for (const BandBall &e : m_engaged)
      if (e.chain >= 0 && e.chain < (int)nc && e.node >= 0 && e.node < (int)bc.loaded[e.chain].size())
         bc.loaded[e.chain][e.node] = 1;
   if (m_slingInBand)
      for (size_t c = 0; c < nc; c++)
         for (const int id : d.chains[c].order)
            if (SlingInFootprint((int)c, id))
               bc.loaded[c][id] = 1;
}

// Nothing on the band: it rings down (band_ringing = 0: critically damped, so it just returns; band_ring_zeta: the damping of the ringing).
// Done when every node is back at rest (band_rest_mm).
void HitRubber::BandRingDown()
{
   BandCycle &bc = m_bc;
   BandDyn &d = m_dyn;
   const size_t nc = d.chains.size();
   const double zRing = !m_set.bandRinging ? 1.0 : m_set.bandRingZeta;   // < 0 = the nodes' own damping ratio
   for (int i = 0; i < bc.steps; i++)
   {
      BandZeroForces();
      if (m_slingInBand)
         SlingBandStep(bc.fNode, bc.h);
      BandChainContacts();
      for (size_t c = 0; c < nc; c++)
         d.chains[c].Step(bc.fNode[c], bc.h, zRing);
      bool added = false;
      for (size_t c = 0; c < nc; c++)
         added = d.chains[c].ApplyStops() || added;
      if (added)
         BandRebuild();
   }
   bool rest = true;
   for (const RubberContact::BandChain &ch : d.chains)
      rest = rest && ch.AtRest(m_set.bandRestMm * 1e-3, m_set.bandRestMm * 1.0);   // m and m/s: 0.05 mm and 0.05 m/s (a node's own vibration is about 1000 rad/s)
   if (rest)
      d = BandDyn();
   BandPublish();
}

bool HitRubber::BandCycleBegin(const float dtime)
{
   using namespace RubberContact;
   BandCycle &bc = m_bc;
   bc.active = false;
   bc.nodeTripped = false;
   if (!m_dyn.active && m_engaged.empty())
      return false;
   if (!m_set.band)   // switched off in the ini while a ball is engaged
   {
      BandReleaseAll("band switched off");
      m_dyn = BandDyn();
      BandPublish();
      return false;
   }
   bc.dt = (double)dtime * DEFAULT_STEPTIME_S;
   if (bc.dt <= 0.0)
      return false;
   const double hMax = std::min(50.0, m_set.bandSubstepUs > 0.0 ? m_set.bandSubstepUs : m_set.substepUs) * 1e-6;
   bc.steps = std::max(1, (int)std::ceil(bc.dt / hMax));
   bc.h = bc.dt / bc.steps;
   BandDyn &d = m_dyn;
   const size_t nc = d.chains.size();
   bc.rp = m_set.ToParams();
   bc.rp.shearScale = 1.0;   // (the old node band stays as it was)
   if (m_set.bandShearZeta >= 0.0)
      bc.rp.shearTan = 2.0 * m_set.bandShearZeta;   // a ball resting on a band must not rock on the shear spring of its contact
   bc.fNode.resize(nc);
   if (m_engaged.empty())
   {
      BandRingDown();
      return false;
   }

   if (m_set.bandFollow)
   {
      // a ball that has rolled away from the node it loads takes a node under itself (placed on the current shape, so nothing is lost); the old one is idle
      // from now on, relaxes, and is merged away when it lies on the line again
      const double thr = m_set.bandFollowMm * 1e-3;
      for (size_t k = 0; k < m_engaged.size(); k++)
      {
         BandBall &e = m_engaged[k];
         if (e.chain < 0 || e.chain >= (int)d.chains.size() || e.node < 0 || e.node >= (int)d.chains[e.chain].nodes.size())
            continue;
         const BandChain &ch0 = d.chains[e.chain];
         const Vertex3Ds p = e.ball->m_d.m_pos;
         const V3 c((double)p.x * U, (double)p.y * U, (double)p.z * U);
         const V3 ab = ch0.B - ch0.A;
         const double fb = (c - ch0.A).Dot(ab) / std::max(ab.Dot(ab), 1e-12);
         const double fn = (ch0.nodes[e.node].P - ch0.A).Dot(ab) / std::max(ab.Dot(ab), 1e-12);   // where the node is now along the rest line (not its material fraction)
         if (!BandFixedNodes() && std::fabs(fb - fn) * ch0.Ls <= thr)
            continue;
         const int oldChain = e.chain, oldNode = e.node;
         e.chain = -1;   // not "another ball" to itself
         int chain = oldChain, node = oldNode;
         if (BandBuild(c, (double)e.ball->m_d.m_radius * U, chain, node, thr))
         {
            e.chain = chain;
            e.node = node;
         }
         else
         {
            e.chain = oldChain;
            e.node = oldNode;
         }
      }
   }
   bc.hLayer = m_set.layerMm > 0.0 ? m_set.layerMm * 1e-3 : 2.0 * m_radiusSI;
   BandMarkLoaded();
   bc.active = true;
   return true;
}

// One substep, first part: the shape as the nodes have it, the forces on the nodes cleared, the sling plate pressing on its span (and stepped)
void HitRubber::BandSubstepBegin()
{
   BandSyncPts();
   BandZeroForces();
   if (m_slingInBand)
      SlingBandStep(m_bc.fNode, m_bc.h);
}

// The ball's contact with the band at its state (c, v, w); the reaction goes into the nodes at the ends of the leg it touches. The caller applies force and torque to the ball.
void HitRubber::BandBallForce(BandBall &e, const RubberContact::BallSI &bs, const RubberContact::V3 &c, const RubberContact::V3 &v, const RubberContact::V3 &w, BandBallContact &out)
{
   using namespace RubberContact;
   BandCycle &bc = m_bc;
   BandDyn &d = m_dyn;
   const size_t nc = d.chains.size();
   const TubeProbe pr = ProbeTube(bs, d.pts.data(), (int)d.pts.size(), m_radiusSI, c);
   // the surface moves with the leg the ball touches: its velocity at the contact, from the velocities of the leg's two ends
   int hc = -1, hj = 0;
   double s = 0.0;
   if (pr.overlap)
      for (size_t k = 0; k < nc && hc < 0; k++)
      {
         const BandChain &ch = d.chains[k];
         if (ch.Count() > 0 && pr.seg >= d.base[k] && pr.seg <= d.base[k] + ch.Count())
         {
            hc = (int)k;
            hj = pr.seg - d.base[k];
            const V3 p = ch.PtPos(hj), q = ch.PtPos(hj + 1);
            const double len = (q - p).Length();
            s = len > 1e-9 ? std::min(1.0, std::max(0.0, (pr.point - p).Length() / len)) : 0.0;
         }
      }
   const V3 vs = hc >= 0 ? d.chains[hc].PtVel(hj) * (1.0 - s) + d.chains[hc].PtVel(hj + 1) * s : V3();
   const ContactResult r = EvaluateTube(bc.rp, bs, pr, bc.hLayer, v, w, bc.h, e.st, vs);
   out.inContact = r.inContact;
   out.vs = vs;
   out.normal = pr.normal;
   out.gap = pr.dist - (bs.R + m_radiusSI);
   if (r.inContact)
   {
      out.force = r.force;
      out.torque = r.torqueOnBall;
      out.delta = r.delta;
      d.peakFn = std::max(d.peakFn, r.fn);
      if (hc >= 0)
      {
         const V3 reaction = r.force * -1.0;   // on the band, shared between the two ends of the leg (an anchor's share is carried by the post)
         const int n0 = d.chains[hc].PtNode(hj), n1 = d.chains[hc].PtNode(hj + 1);
         if (n0 >= 0)
            bc.fNode[hc][n0] = bc.fNode[hc][n0] + reaction * (1.0 - s);
         if (n1 >= 0)
            bc.fNode[hc][n1] = bc.fNode[hc][n1] + reaction * s;
      }
   }
}

// One substep, last part: the spans push on each other, the nodes step, the stops. True if a node left the model's range.
bool HitRubber::BandSubstepEnd()
{
   using namespace RubberContact;
   BandCycle &bc = m_bc;
   BandDyn &d = m_dyn;
   const size_t nc = d.chains.size();
   BandChainContacts();
   for (size_t c = 0; c < nc; c++)
   {
      d.chains[c].Step(bc.fNode[c], bc.h);
      for (const int k : d.chains[c].order)
         if (!(d.chains[c].nodes[k].Deflection() <= 0.4 * d.chains[c].Ls))   // written so that a NaN trips too
            bc.nodeTripped = true;
   }
   bool added = false;
   for (size_t c = 0; c < nc; c++)
      added = d.chains[c].ApplyStops() || added;
   if (added)
      BandRebuild();
   return bc.nodeTripped;
}

// After the nodes stepped: the peak deflection of the ball's node; returns the sum of the node's coordinates (for the finite check of the ball's state)
double HitRubber::BandBallPeak(BandBall &e)
{
   const RubberContact::BandChain::Node &ln = m_dyn.chains[e.chain].nodes[e.node];
   e.peak = std::max(e.peak, ln.Deflection());
   return ln.P.x + ln.P.y + ln.P.z;
}

// Once per cycle after the substeps. False = a span deflected beyond the model's range: the band is reset and the balls released (untouched).
bool HitRubber::BandCycleEnd()
{
   BandCycle &bc = m_bc;
   BandDyn &d = m_dyn;
   const size_t nc = d.chains.size();
   for (BandBall &e : m_engaged)
      e.dwellS += bc.dt;
   if (m_set.bandFollow && !BandFixedNodes())
   {
      // idle nodes (their ball has left or moved on) that have relaxed onto the line between their neighbours are merged away
      bool merged = false;
      for (size_t c = 0; c < nc; c++)
      {
         std::vector<RubberContact::V3> guards;   // where balls on the other spans load the band: the nodes here are companions of those contacts
         for (const BandBall &e : m_engaged)
            if (e.chain >= 0 && e.chain < (int)nc && (size_t)e.chain != c && e.node >= 0 && e.node < (int)d.chains[e.chain].nodes.size())
               guards.push_back(d.chains[e.chain].nodes[e.node].P);
         merged = d.chains[c].MergeIdle(bc.loaded[c], std::max(m_set.bandRestMm, 0.3) * 1e-3, 3.0, guards, 0.03) || merged;
      }   // the speed limit is generous: a light idle node vibrates along the cord for a long time, with well under a mJ in it
      if (merged)
         BandRebuild();
   }
   if (bc.nodeTripped)
   {
      PLOGW << "Band " << m_name << ": a span deflected beyond the model's range, using the stock collision";
      BandReleaseAll("tripped");
      d = BandDyn();
      BandPublish();
      m_bandCooldownUntilMs = g_pplayer->m_time_msec + 20;
      return false;   // the balls are untouched; the engine moves them as usual
   }
   return true;
}

// The release test for a ball on the band: clear of the band by more than the gap and moving away from it, for 1.5 ms, with no contact in between
bool HitRubber::BandBallClear(BandBall &e, const RubberContact::BallSI &bs, const RubberContact::V3 &c, const bool contact, const double lastGap, const double relNormalVel)
{
   const double releaseGap = 0.2e-3;   // m, the same hysteresis as the flipper's rubber_release_gap_mm default
   bool clear = !contact && lastGap > releaseGap && relNormalVel > 0.0;
   if (clear)
   {
      // not while the ball still sits in the band's rest position: the engine's collision test knows only the rest shape, so it would take the ball straight
      // back (a heavy ball following a deeply deflected band is clear of the band but inside its rest position)
      const RubberContact::TubeProbe rest = RubberContact::ProbeTube(bs, m_ptsSI.data(), (int)m_ptsSI.size(), m_radiusSI, c);
      if (rest.dist - (bs.R + m_radiusSI) <= releaseGap)
         clear = false;
   }
   e.clearS = clear ? e.clearS + m_bc.dt : 0.0;
   return e.clearS >= 1.5e-3;
}

// The legacy driver: the balls on the band are integrated here, each against the shared deformed shape.
void HitRubber::BandAdvanceCore(const float dtime)
{
   using namespace RubberContact;
   if (!BandCycleBegin(dtime))
      return;
   BandCycle &bc = m_bc;
   const double dt = bc.dt, h = bc.h;
   struct Work
   {
      BandBall *e;
      BallSI bs;
      Vertex3Ds posStart;
      V3 c, v, w;
      double lScale = 0.0, vCap = 0.0;
      bool contact = false, tripped = false;
      double lastGap = 0.0, relNormalVel = 0.0, deepest = 0.0;
   };
   std::vector<Work> ws;
   for (BandBall &e : m_engaged)
   {
      Work wk;
      HitBall *const pb = e.ball;
      wk.e = &e;
      wk.bs.R = (double)pb->m_d.m_radius * U;
      wk.bs.m = kMassUnitKg * (double)pb->m_d.m_mass;
      const double inertiaU = (double)pb->Inertia();
      const double wScale = 1.0 / (inertiaU * DEFAULT_STEPTIME_S);
      wk.lScale = inertiaU * DEFAULT_STEPTIME_S;
      wk.posStart = pb->m_d.m_pos;
      wk.c = V3((double)wk.posStart.x * U, (double)wk.posStart.y * U, (double)wk.posStart.z * U);
      wk.v = V3((double)pb->m_d.m_vel.x * VS, (double)pb->m_d.m_vel.y * VS, (double)pb->m_d.m_vel.z * VS);
      wk.w = V3((double)pb->m_angularmomentum.x * wScale, (double)pb->m_angularmomentum.y * wScale, (double)pb->m_angularmomentum.z * wScale);
      wk.vCap = std::max(50.0, 4.0 * wk.v.Length() + 5.0);   // m/s
      ws.push_back(wk);
   }

   for (int i = 0; i < bc.steps && !bc.nodeTripped; i++)
   {
      BandSubstepBegin();
      for (Work &wk : ws)
      {
         if (wk.tripped)
            continue;
         BandBallContact bcn;
         BandBallForce(*wk.e, wk.bs, wk.c, wk.v, wk.w, bcn);
         if (bcn.inContact)
         {
            wk.contact = true;
            wk.deepest = std::max(wk.deepest, bcn.delta);
            wk.v = wk.v + bcn.force * (1.0 / wk.bs.m) * h;
            wk.w = wk.w + bcn.torque * (1.0 / wk.bs.Inertia()) * h;
         }
         wk.lastGap = bcn.gap;
         wk.relNormalVel = (wk.v - bcn.vs).Dot(bcn.normal);
      }
      BandSubstepEnd();
      for (Work &wk : ws)
      {
         if (wk.tripped)
            continue;
         const double nodeSum = BandBallPeak(*wk.e);
         wk.c = wk.c + wk.v * h;
         if (!std::isfinite(wk.c.x + wk.c.y + wk.c.z + wk.v.x + wk.v.y + wk.v.z + wk.w.x + wk.w.y + wk.w.z + nodeSum) || wk.v.Length() > wk.vCap || wk.deepest > 0.8 * wk.bs.R)
            wk.tripped = true;
      }
   }
   if (!BandCycleEnd())
      return;

   static const char *const logFile = getenv("VPX_BAND_LOG");   // diagnostics: one row per cycle and ball with the ball and every node (VPX units)
   static FILE *const f = (logFile && *logFile) ? fopen(logFile, "w") : nullptr;

   BandDyn &d = m_dyn;
   std::vector<std::pair<HitBall *, const char *>> releases;
   for (size_t b = 0; b < ws.size(); b++)
   {
      Work &wk = ws[b];
      HitBall *const pb = wk.e->ball;
      if (wk.tripped)
      {
         PLOGW << "Band " << m_name << ": ball step left the model's range (depth " << wk.deepest * 1e3 << " mm), using the stock collision for it";
         m_bandCooldownUntilMs = g_pplayer->m_time_msec + 20;
         releases.emplace_back(pb, "tripped");   // the ball is untouched; the engine moves it as usual
         continue;
      }
      pb->m_d.m_pos = Vertex3Ds((float)(wk.c.x / U), (float)(wk.c.y / U), (float)(wk.c.z / U));
      pb->m_d.m_vel = Vertex3Ds((float)(wk.v.x / VS), (float)(wk.v.y / VS), (float)(wk.v.z / VS));
      pb->m_angularmomentum = Vertex3Ds((float)(wk.w.x * wk.lScale), (float)(wk.w.y * wk.lScale), (float)(wk.w.z * wk.lScale));
      const Vertex3Ds ds = pb->m_d.m_pos - wk.posStart;
      pb->m_lastEventSqrDist += ds.LengthSquared();
#ifdef C_DYNAMIC
      pb->m_drsq = ds.LengthSquared();
#endif
      pb->m_rubberAdvanced = true;   // the ball has been moved for the whole cycle
      if (f)
      {
         fprintf(f, "%u,%u,%.3f,%.4f,%.4f,%.4f,%.4f", g_pplayer->m_time_msec, (unsigned)b, wk.e->dwellS * 1e3, wk.c.x / U, wk.c.y / U, wk.v.x, wk.v.y);
         for (const BandChain &ch : d.chains)
            for (int k = 1; k <= ch.Count(); k++)
               fprintf(f, ",%.4f,%.4f,%.4f", ch.PtPos(k).x / U, ch.PtPos(k).y / U, ch.nodes[ch.PtNode(k)].Deflection() * 1e3);
         fprintf(f, ",%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f", d.peakFn, d.peakBB * 1e3, wk.c.z / U, wk.v.z, wk.w.x, wk.w.y, wk.w.z);   // then: ball z (U), vz (m/s), spin (rad/s)
         fprintf(f, "%c", 10);
         fflush(f);
      }
      if (BandBallClear(*wk.e, wk.bs, wk.c, wk.contact, wk.lastGap, wk.relNormalVel))
         releases.emplace_back(pb, "clear");
   }
   (void)dt;
   for (const auto &r : releases)
      BandRelease(r.first, r.second);
   BandPublish();
}

// The mesh animation: where the rendered tube goes when the band is deflected. The ring of the mesh centred at c (table units, in the plane of the table
// what matters) is matched to the nearest point of the rest centreline; if that is on a span with nodes, the ring follows the cord: its centre moves from
// the rest chord to the cord's place at the same material fraction, and it turns by the angle between the rest chord and the leg there.
bool HitRubber::BandDeformAt(const Vertex3Ds &c, Vertex3Ds &shift, float rot[9]) const
{
   using namespace RubberContact;
   if (!m_dyn.active || m_pts.size() < 3)
      return false;
   const int n = (int)m_pts.size() - 1;
   int bi = -1;
   double bd = 1e30;
   for (int i = 0; i < n; i++)
   {
      const double ax = m_pts[i].x, ay = m_pts[i].y, bx = m_pts[i + 1].x - ax, by = m_pts[i + 1].y - ay;
      const double l2 = bx * bx + by * by;
      const double t = l2 > 1e-12 ? std::min(std::max(((c.x - ax) * bx + (c.y - ay) * by) / l2, 0.0), 1.0) : 0.0;
      const double dx = c.x - (ax + bx * t), dy = c.y - (ay + by * t);
      if (dx * dx + dy * dy < bd)
      {
         bd = dx * dx + dy * dy;
         bi = i;
      }
   }
   if (bi < 0)
      return false;
   int sk = -1;
   for (size_t k = 0; k < m_bandSpans.size() && sk < 0; k++)
   {
      const BandSpan &sp = m_bandSpans[k];
      if (sp.free && ((bi - sp.first % n) % n + n) % n < sp.last - sp.first)
         sk = (int)k;
   }
   if (sk < 0)
      return false;
   int ci = -1;
   for (size_t k = 0; k < m_dyn.span.size(); k++)
      if (m_dyn.span[k] == sk)
         ci = (int)k;
   if (ci < 0 || m_dyn.chains[ci].Count() == 0)
      return false;
   const BandChain &ch = m_dyn.chains[ci];
   const double px = (double)c.x * U, py = (double)c.y * U;
   const double abx = ch.B.x - ch.A.x, aby = ch.B.y - ch.A.y, l2 = abx * abx + aby * aby;
   if (l2 < 1e-12)
      return false;
   const double f = std::min(std::max(((px - ch.A.x) * abx + (py - ch.A.y) * aby) / l2, 0.0), 1.0);
   int j = 0;   // the leg that holds the material fraction f
   while (j < ch.Count() && ch.PtF(j + 1) < f)
      j++;
   const double f0 = ch.PtF(j), f1 = ch.PtF(j + 1);
   const double u = f1 > f0 ? std::min(std::max((f - f0) / (f1 - f0), 0.0), 1.0) : 0.0;
   const V3 P = ch.PtPos(j) + (ch.PtPos(j + 1) - ch.PtPos(j)) * u;
   const V3 rest = ch.A + (ch.B - ch.A) * f;
   shift = Vertex3Ds((float)((P.x - rest.x) / U), (float)((P.y - rest.y) / U), (float)((P.z - rest.z) / U));
   const V3 t0 = (ch.B - ch.A) / (ch.B - ch.A).Length();
   const V3 d1 = ch.PtPos(j + 1) - ch.PtPos(j);
   const V3 t1 = d1 / std::max(d1.Length(), 1e-12);
   const V3 ax = t0.Cross(t1);
   const double sn = ax.Length(), cs = t0.Dot(t1);
   double R[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
   if (sn > 1e-9)
   {
      const V3 k = ax / sn;   // Rodrigues: R = I + sin K + (1 - cos) K^2, K the cross-product matrix of the unit axis
      const double K[9] = { 0, -k.z, k.y, k.z, 0, -k.x, -k.y, k.x, 0 };
      double K2[9];
      for (int r = 0; r < 3; r++)
         for (int q = 0; q < 3; q++)
            K2[r * 3 + q] = K[r * 3] * K[q] + K[r * 3 + 1] * K[3 + q] + K[r * 3 + 2] * K[6 + q];
      for (int i = 0; i < 9; i++)
         R[i] += sn * K[i] + (1.0 - cs) * K2[i];
   }
   for (int i = 0; i < 9; i++)
      rot[i] = (float)R[i];
   return true;
}

// A ball touches this band but the band model does not take it: the stock (rigid) collision runs instead. Says why, at most twice a second.
void HitRubber::BandNote(const HitBall *pb) const
{
   const long long now = (long long)g_pplayer->m_time_msec;
   if (now - m_bandNoteMs < 500)
      return;
   m_bandNoteMs = now;
   const char *why = "the band model refused it (the node limit, or the node it would load is loaded by another ball)";
   if (m_bandSpans.empty())
      why = "this rubber has no free span";
   else if (g_pplayer->m_time_msec < m_bandCooldownUntilMs)
      why = "cooling down after a span was deflected beyond the model's range";
   else if (pb->m_bandOwner != nullptr || pb->m_contact != nullptr)
      why = "the ball is held by another band";
   else if (pb->m_rubberOwner != nullptr)
      why = "the ball is held by a force-model or flipper rubber contact";
   else if (g_pplayer->m_physics->FlipperNear(pb->m_d, (float)m_set.bandFlipperMarginU))
      why = "a flipper is within band_flipper_margin_u of the ball";
   else
   {
      using namespace RubberContact;
      const V3 c((double)pb->m_d.m_pos.x * U, (double)pb->m_d.m_pos.y * U, (double)pb->m_d.m_pos.z * U);
      const int n = (int)m_pts.size() - 1;
      double bestD = 1e30;
      for (const BandSpan &sp : m_bandSpans)
      {
         if (!sp.free || sp.last - sp.first < 1 || sp.last - sp.first >= n)
            continue;
         const V3 &a = m_ptsSI[sp.first % n], &b = m_ptsSI[sp.last % n];
         const V3 ab = b - a;
         const double l2 = ab.Dot(ab);
         if (l2 < 1e-12)
            continue;
         const double f = std::min(std::max((c - a).Dot(ab) / l2, 0.0), 1.0);
         bestD = std::min(bestD, (c - (a + ab * f)).Length());
      }
      if (bestD > (double)pb->m_d.m_radius * U + m_radiusSI + 0.03)
         why = "the ball is not at a free span (it is on a post wrap or a short span)";
   }
   PLOGI << "Band " << m_name << ": ball at (" << pb->m_d.m_pos.x << ", " << pb->m_d.m_pos.y << ") touches the band but is not carried by the band model, stock collision: " << why;
}

void HitRubber::Contact(CollisionEvent &coll, const float dtime)
{
   if (m_set.band)
   {
      if (BandEngage(coll.m_ball))
         return;   // a resting or slow touch: the band's force model carries the ball from here on
      BandNote(coll.m_ball);
   }
   else if (SlowContact(coll.m_ball))
      return;   // a post rubber: the force model carries the resting or slow touch
   HitObject::Contact(coll, dtime);
}

void HitRubber::DrawUI(std::function<Vertex2D(Vertex3Ds)> project, ImDrawList *drawList, bool fill) const
{
   if (!m_enabled)
      return;
   const ImU32 lCol = ImGui::GetColorU32(ImGuiCol_PlotLines);
   for (size_t i = 0; i + 1 < m_pts.size(); i++)
   {
      const Vertex2D p0 = project(m_pts[i]);
      const Vertex2D p1 = project(m_pts[i + 1]);
      if (p0.x != FLT_MAX && p1.x != FLT_MAX)
         drawList->AddLine(ImVec2(p0.x, p0.y), ImVec2(p1.x, p1.y), lCol);
   }
}
