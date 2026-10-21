// license:GPLv3+

#include "core/stdafx.h"
#include "physics/magnet.h"
#include "physics/physicsini.h"
#include "physics/physconst.h"
#include "physics/magnetcore.h"
#include "physics/magnetballs.h"
#include "parts/pintable.h"
#include "core/player.h"
#include "parts/ball.h"
#include "ui/live/LiveUI.h"
#include "renderer/Renderer.h"
#include "renderer/RenderDevice.h"
#include "renderer/MeshBuffer.h"
#include "renderer/VertexBuffer.h"
#include "renderer/IndexBuffer.h"
#include "renderer/Shader.h"
#include <cmath>
#include <map>
#include <atomic>
#include <chrono>
#include <mutex>
#include <thread>
#include <sstream>
#include <unordered_map>
#include <array>

namespace Magnets
{
   namespace
   {
      constexpr double MU0 = 4.0e-7 * 3.14159265358979323846;
      constexpr double PI = 3.14159265358979323846;
      constexpr int NPTS = 24;
      constexpr double REACH = 0.127;   // m: beyond 5 in from the pole face the pull is below rolling resistance (under 0.005 m/s2 for a 9 A coil), so it is skipped
                                        // (air-core and circuit magnets; core_model = cells: its tables reach as far as the field is above MagnetCore::B_GATE)

      struct Magnet
      {
         std::string name;
         // ini
         double x = 0, y = 0, z = 0;           // centre of the near pole face, table units (U)
         double ax = 0, ay = 0, az = 1;        // axis, pointing out of the near face towards the ball
         double turns = 2000;                  // N
         double wireMm = 0.45;                 // wire diameter including insulation, mm (0.405 bare = 26 AWG)
         double rOhm = 0;                      // 0 = computed from the winding
         bool aircore = false;                 // type = aircore: an air-core coil the ball rolls through (no iron core, no pole face); x, y, z are the coil centre
         double coilLenMm = 40;                // aircore: winding length
         double boreMm = 32;                   // aircore: inner diameter of the winding
         bool ballInteraction = true;          // the balls this magnet magnetises attract and repel each other (induced dipoles)
         bool tube = true;                     // aircore: a guide tube inside the bore confines the ball to the axis (the tube's outer diameter is the bore)
         double tubeWallMm = 1.5875;           // aircore: tube wall thickness (1/16 in); the ball runs in the bore less a wall on each side
         double tubeE = 0.5;                   // aircore: restitution of the ball on the tube and on the stop (steel on plastic)
         double tubeMu = 0.2;                  // aircore: friction coefficient of the ball on the tube and on the stop
         bool tubeVisible = false;             // aircore: draw the tube and the stop (translucent) in the playfield view
         std::shared_ptr<MeshBuffer> mesh;     // the drawn tube, built on first render
         bool stop = false;                    // aircore: a stop across the tube ends the ball's travel
         double stopPosMm = 1e9;               // aircore: stop face position along the axis from the coil centre (signed; default = the end of the bore, half the coil length)
         double tubeRU = 0, tubeHalfU = 0, stopU = 0;   // derived, table units: tube inner radius, half the tube length, stop face position
         double housingRU = 0, housingHalfU = 0;   // derived, table units: winding housing outer radius (winding + a tube wall) and half length (end flanges a tube wall
                                               // beyond the winding); solid with the tube, which runs through the flanges
         double solRadius[3] = { 0, 0, 0 };    // aircore: radii of the three winding layers (m), each with a third of the turns
         double solLen = 0;
         double windIn = 0, windOut = 0;       // aircore: inner and outer radius of the winding block (m)
         std::vector<float> tab;               // aircore: per node of a grid over (z along the axis >= 0, r from the axis): the gradient of the ball-energy orders Q_1..Q_4 per A^2, then Bz, Br per A
         mutable bool overlapWarned = false;   // aircore: "ball inside the winding" raised once per table start
         double coreDiaMm = 19;                // iron core diameter
         double coreLenMm = 38;                // core length = winding length
         double returnDistMm = 0;              // > 0: the iron yoke returns the flux to flange poles in the face plane; distance from the core axis to the flanges' inner edge
         double returnWidthMm = 12.7;          // radial width of a flange
         double returnLengthMm = 57;           // length of a flange along its edge
         double returnDepthMm = 0;             // how far below the pole face plane the flanges sit (a core that reaches the playfield surface: the playfield thickness)
         double returnCount = 2;               // number of flanges (1 = one-sided bracket, 2 = opposite sides, 4 = all round), spread evenly
         double returnAngleDeg = 0;            // direction of the first flange in the face plane (0 = +x for an axis pointing up, counter-clockwise)
         double returnFraction = 0;            // 0 = bare rod (the far pole is open to space); 1 = a closed iron return (no far-pole reluctance)
         double bSatT = 1.6;                   // core saturation flux density
         double eddySigma = 4.5e6;             // ball steel conductivity, S/m (chrome steel 52100: resistivity about 2.2e-7 ohm m)
         double eddyScale = 1.0;               // ball_eddy: scale of the eddy-current braking of a spinning ball (0 = off)
         double msAm = 1.35e6;                 // ball saturation magnetisation, A/m (chrome steel ~1.7 T)
         // iron core, core_model = cells (2026-10-18, magnetcore.h): the core's answer to the ball from load-time tables
         bool cells = true;                    // core_model: cells (default) or circuit (the old point-charge circuit; the only one that uses return_* yokes)
         double coreExposedMm = 6;             // core_exposed_mm: core length above the winding, up to the pole face
         std::vector<double> operatingVolts;   // operating_volts: the coil voltages the tables are built for (default: rated_volts)
         std::shared_ptr<MagnetCore::Model> core;   // built at load (cells)
         std::shared_ptr<MagnetBalls::MagnetDrive> drive;   // cells + ball_interaction: the Option A ball-ball drive tables (magnetballs.h), built at load
         mutable bool radiusWarned = false;    // cells: "a ball of another size" raised once per table start
         // derived
         double r = 0, l = 0, tau = 0, spread = 0, area = 0, radius = 0, length = 0;
         double nx = 0, ny = 0, nz = 1, ux = 1, uy = 0, uz = 0, vx = 0, vy = 1, vz = 0;
         double spreadY = 0;                   // spreading conductance of the flange poles together
         int nY = 0;
         double ypts[16][3];                   // flange pole point positions (m)
         double pts[NPTS][3];                  // point charge positions relative to the face centre (m)
         // state
         double ratedVolts = 48;               // rated_volts: the supply the script's MagnetDuty is a fraction of
         double volts = 0, amps = 0;
      };

      std::vector<Magnet> g_magnets;
      struct OptAWarnings { bool mixed = false, cost = false; } g_optAWarn;   // Option A ball-ball notices, once per table start
      MagnetBalls::Solver g_optASolver;
      // this physics step's Option A parts per ball (build step B7); cleared at the start of every UpdateBallInteraction
      struct BallParts { double massKg = 0; std::vector<BallBallPart> parts; };
      std::unordered_map<const HitBall *, BallParts> g_parts;
      std::unordered_map<const HitBall *, std::array<double, 3>> g_partStart;   // the balls' centres (m) when the step's solve ran   // Option A solver: warm state and the last result (gate (b)); its last result is cleared at table start
      std::map<std::string, size_t> g_byName;
      const void *g_owner = nullptr;
      bool g_loaded = false, g_anyCurrent = false;
      std::filesystem::path g_file;             // the ini the magnets were read from
      std::filesystem::file_time_type g_fileTime{};
      uint32_t g_nextCheckMs = 0;
      std::map<std::string, std::pair<bool, float>> g_pending;   // what the script set before the magnets were loaded (e.g. a call at the top of the script): (isDuty, volts or duty)

      constexpr int TN = 128;              // table nodes per direction
      constexpr double TH = 1.0e-3;        // table spacing, m (the table spans 127 mm, the cutoff of the magnet force)

      // Field of one circular current loop (radius a, 1 A) at axial distance z and radial distance rho, Smythe's elliptic integral form (T per A).
      inline void LoopField(double a, double z, double rho, double &bz, double &br)
      {
         const double z2 = z * z;
         if (rho < 1e-7)
         {
            bz = MU0 * a * a / (2 * std::pow(a * a + z2, 1.5));
            br = 0;
            return;
         }
         const double ap = (a + rho) * (a + rho) + z2, am = (a - rho) * (a - rho) + z2;
         const double m = std::min(4 * a * rho / ap, 1 - 1e-12);
         // complete elliptic integrals K(m), E(m) by the arithmetic-geometric mean
         double A = 1.0, G = std::sqrt(1 - m), sum = 0.5 * m, pw = 0.5;
         for (int i = 0; i < 20; ++i)
         {
            const double An = 0.5 * (A + G), c = 0.5 * (A - G);
            G = std::sqrt(A * G);
            A = An;
            pw *= 2;
            sum += pw * c * c;
            if (c < 1e-15)
               break;
         }
         const double K = PI / (2 * A), E = K * (1 - sum);
         const double pre = MU0 / (2 * PI) / std::sqrt(ap);
         bz = pre * (K + (a * a - rho * rho - z2) / am * E);
         br = pre * z / rho * (-K + (a * a + rho * rho + z2) / am * E);
      }

      // Air-core coil: the lone ball's energy in the coil field, order by order (the multipole response of a sphere, to order 4).
      // For a steel ball (mu >> 1: M = 3H in a uniform field) E = sum_l a^(2l+1) Q_l, Q_l = 2 pi mu0 (2l+1) <p_l^2>, where p_l is the degree-l part of the coil's
      // scalar potential about the ball centre (H = -grad phi) and <> the mean over the unit sphere. Q_l belongs to the field alone, so one table serves every
      // ball size. Q_1 = 2 pi B^2 / mu0 is the dipole term (the old law, F = V M grad|B|). Q_2..Q_4 come from the coil field on a sphere round the node,
      // projected on spherical harmonics: the radial field there is -sum_l l s^(l-1) Y_l, so Q_l = mu0 (2l+1) sum_m H_lm^2 / (2 l^2 s^(2l-2)) for a sphere of radius s.
      // Per node (z along the axis >= 0, r from the axis): grad Q_1..Q_4 (d/dz, d/dr) per A^2, then Bz, Br per A.
      constexpr int AC_L = 4;                   // multipole order of the ball's response
      constexpr int AC_NF = 2 * AC_L + 2;       // floats per table node
      constexpr int AC_RINGS = 8;               // rings of loops per winding layer
      constexpr double AC_SREF = 13.5e-3;       // m: projection sphere radius (a standard ball); smaller next to the winding
      constexpr double AC_NEAR = 12e-3;         // m: beyond AC_SREF + this from the winding, Q_l / Q_1 is smooth and read from a coarse lattice
      constexpr int AC_STRIDE = 4;              // coarse lattice spacing, nodes
      constexpr double AC_TOL = 1e-4;           // quadrature: aliasing estimate (s / nearest loop)^(2 n - 4) below this

      inline void CoilField(const Magnet &m, double z, double rho, double &bz, double &br)   // T per A, (z, rho) from the coil centre
      {
         bz = br = 0;
         const double w = m.turns / (3.0 * AC_RINGS);
         for (int k = 0; k < 3; ++k)
            for (int j = 0; j < AC_RINGS; ++j)
            {
               double lz, lr;
               LoopField(m.solRadius[k], z - (-0.5 * m.solLen + (j + 0.5) * m.solLen / AC_RINGS), rho, lz, lr);
               bz += w * lz;
               br += w * lr;
            }
      }

      // distance from (z, r) to the winding block (r windIn..windOut, |z| <= half the length) in the coil's meridian plane
      inline double WindingGap(const Magnet &m, double z, double r)
      {
         const double dz = std::max(0.0, std::fabs(z) - 0.5 * m.solLen), dr = std::max(0.0, std::max(m.windIn - r, r - m.windOut));
         return std::sqrt(dz * dz + dr * dr);
      }

      // Sphere quadrature: Gauss-Legendre in cos(theta) x midpoints in phi over [0, pi] (the field is even about the meridian plane), with the real
      // orthonormal spherical harmonics Y_lm (cos part, l = 1..AC_L, m = 0..l) at each point.
      struct SphereRule { int n = 0; std::vector<double> ct, st, cp, sp, w, Y; };
      constexpr int AC_NLM = (AC_L + 1) * (AC_L + 2) / 2 - 1;
      const std::vector<SphereRule> &SphereRules()
      {
         static std::vector<SphereRule> rules;
         static std::once_flag once;
         std::call_once(once, []
         {
            for (const int n : { 6, 8, 10, 12, 16, 20, 24, 32 })
            {
               SphereRule q;
               q.n = n;
               std::vector<double> x(n), wx(n);
               for (int i = 0; i < n; ++i)   // Gauss-Legendre nodes by Newton on P_n
               {
                  double t = std::cos(PI * (i + 0.75) / (n + 0.5)), dp = 1;
                  for (int it = 0; it < 100; ++it)
                  {
                     double p1 = 1, p2 = 0;
                     for (int j = 1; j <= n; ++j) { const double p3 = p2; p2 = p1; p1 = ((2 * j - 1) * t * p2 - (j - 1) * p3) / j; }
                     dp = n * (t * p1 - p2) / (t * t - 1);
                     const double d = p1 / dp;
                     t -= d;
                     if (std::fabs(d) < 1e-15) break;
                  }
                  x[i] = t; wx[i] = 2 / ((1 - t * t) * dp * dp);
               }
               for (int i = 0; i < n; ++i)
                  for (int j = 0; j < n; ++j)
                  {
                     const double ph = (j + 0.5) * PI / n, st = std::sqrt(1 - x[i] * x[i]);
                     q.ct.push_back(x[i]); q.st.push_back(st); q.cp.push_back(std::cos(ph)); q.sp.push_back(std::sin(ph));
                     q.w.push_back(wx[i] * 2 * PI / n);
                     for (int l = 1; l <= AC_L; ++l)
                        for (int mm = 0; mm <= l; ++mm)
                        {
                           double pmm = 1;   // associated Legendre P_l^m(x), no Condon-Shortley phase (it squares away)
                           for (int k = 1; k <= mm; ++k) pmm *= (2 * k - 1) * st;
                           double plm = pmm;
                           if (l > mm)
                           {
                              double pm1 = x[i] * (2 * mm + 1) * pmm;
                              plm = pm1;
                              for (int ll = mm + 2; ll <= l; ++ll) { plm = (x[i] * (2 * ll - 1) * pm1 - (ll + mm - 1) * pmm) / (ll - mm); pmm = pm1; pm1 = plm; }
                           }
                           double f = 1;
                           for (int k = l - mm + 1; k <= l + mm; ++k) f *= k;
                           q.Y.push_back(std::sqrt((2 * l + 1) / (4 * PI) / f) * (mm ? std::sqrt(2.0) : 1.0) * plm * std::cos(mm * ph));
                        }
                  }
               rules.push_back(std::move(q));
            }
         });
         return rules;
      }

      // Q_2..Q_AC_L at node (z, r) per A^2 (Q[0] = Q_1 is not touched)
      void ProjectOrders(const Magnet &m, double z, double r, double Q[AC_L])
      {
         for (int l = 2; l <= AC_L; ++l) Q[l - 1] = 0;
         const double s = std::min(AC_SREF, 0.9 * WindingGap(m, z, r));
         if (s < 0.5e-3)
            return;   // the node is in or on the winding: any ball here overlaps it (the force keeps the dipole term there)
         double rhoMin = 1e9;
         for (int k = 0; k < 3; ++k)
            for (int j = 0; j < AC_RINGS; ++j)
               rhoMin = std::min(rhoMin, std::hypot(r - m.solRadius[k], z - (-0.5 * m.solLen + (j + 0.5) * m.solLen / AC_RINGS)));
         const auto &rules = SphereRules();
         const double q = s / rhoMin;
         size_t pick = rules.size() - 1;
         for (size_t i = 0; i < rules.size(); ++i)
            if (std::pow(q, 2 * rules[i].n - 4) < AC_TOL) { pick = i; break; }
         const SphereRule &R = rules[pick];
         double f[AC_NLM] = {};
         for (size_t k = 0; k < R.w.size(); ++k)
         {
            const double nx = R.st[k] * R.cp[k], ny = R.st[k] * R.sp[k], nz = R.ct[k];
            const double px = r + s * nx, py = s * ny, pz = z + s * nz, rho = std::sqrt(px * px + py * py);
            double bz, br;
            CoilField(m, pz, rho, bz, br);
            const double bn = (rho > 1e-12 ? br * (px * nx + py * ny) / rho : 0.0) + bz * nz;   // B.n = mu0 H_r
            const double *Y = &R.Y[k * AC_NLM];
            for (int i = 0; i < AC_NLM; ++i) f[i] += R.w[k] * bn * Y[i];
         }
         int i = 0;
         for (int l = 1; l <= AC_L; ++l)
         {
            double sum = 0;
            for (int mm = 0; mm <= l; ++mm, ++i) sum += f[i] * f[i];
            if (l >= 2)
               Q[l - 1] = (2 * l + 1) * sum / (2 * MU0 * l * l * std::pow(s, 2 * l - 2));   // f = mu0 H_lm: mu0 (2l+1) sum H_lm^2 / (2 l^2 s^(2l-2))
         }
      }

      // The table, built once at load. Q_1 and the field at every node; Q_2..Q_4 by projection near the winding and on a coarse lattice elsewhere
      // (far nodes take Q_l / Q_1 bilinear from the lattice: the ratio is smooth there). Spread over threads.
      void BuildAircoreTable(Magnet &m)
      {
         const size_t NN = (size_t)TN * TN;
         std::vector<double> Q(NN * AC_L, 0.0), BZ(NN), BR(NN);
         std::vector<char> exact(NN, 0);
         std::vector<int> work;
         auto onLattice = [](int i) { return i % AC_STRIDE == 0 || i == TN - 1; };
         for (int iz = 0; iz < TN; ++iz)
            for (int ir = 0; ir < TN; ++ir)
               if (WindingGap(m, iz * TH, ir * TH) < AC_SREF + AC_NEAR || (onLattice(iz) && onLattice(ir)))
               {
                  work.push_back(iz * TN + ir);
                  exact[(size_t)iz * TN + ir] = 1;
               }
         std::atomic<int> next { 0 };
         const int total = (int)NN + (int)work.size();
         auto worker = [&]()
         {
            for (int w; (w = next.fetch_add(1)) < total;)
            {
               if (w < (int)NN)   // field and Q_1 at every node
               {
                  double bz, br;
                  CoilField(m, (w / TN) * TH, (w % TN) * TH, bz, br);
                  BZ[w] = bz; BR[w] = br;
                  Q[(size_t)w * AC_L] = 2 * PI * (bz * bz + br * br) / MU0;
               }
               else
               {
                  const int n = work[w - (int)NN];
                  double q[AC_L];
                  ProjectOrders(m, (n / TN) * TH, (n % TN) * TH, q);
                  for (int l = 1; l < AC_L; ++l) Q[(size_t)n * AC_L + l] = q[l];
               }
            }
         };
         const int nth = std::max(1, std::min(16, (int)std::thread::hardware_concurrency()));
         std::vector<std::thread> pool;
         for (int i = 1; i < nth; ++i)
            pool.emplace_back(worker);
         worker();
         for (auto &t : pool)
            t.join();
         // far nodes: Q_l / Q_1 bilinear from the lattice
         auto lo = [](int i) { return std::min(i / AC_STRIDE * AC_STRIDE, TN - 1); };
         auto hi = [](int i) { return std::min(i / AC_STRIDE * AC_STRIDE + AC_STRIDE, TN - 1); };
         auto ratio = [&](int iz, int ir, int l) { const size_t n = (size_t)iz * TN + ir; return Q[n * AC_L] > 0 ? Q[n * AC_L + l] / Q[n * AC_L] : 0.0; };
         for (int iz = 0; iz < TN; ++iz)
            for (int ir = 0; ir < TN; ++ir)
            {
               const size_t n = (size_t)iz * TN + ir;
               if (exact[n])
                  continue;
               const int z0 = lo(iz), z1 = hi(iz), r0 = lo(ir), r1 = hi(ir);
               const double wz = z1 > z0 ? double(iz - z0) / (z1 - z0) : 0.0, wr = r1 > r0 ? double(ir - r0) / (r1 - r0) : 0.0;
               for (int l = 1; l < AC_L; ++l)
                  Q[n * AC_L + l] = Q[n * AC_L] * ((1 - wz) * ((1 - wr) * ratio(z0, r0, l) + wr * ratio(z0, r1, l)) + wz * ((1 - wr) * ratio(z1, r0, l) + wr * ratio(z1, r1, l)));
            }
         m.tab.assign(NN * AC_NF, 0.f);
         for (int iz = 0; iz < TN; ++iz)
            for (int ir = 0; ir < TN; ++ir)
            {
               auto at = [&](int a, int b, int l) { return Q[((size_t)std::min(std::max(a, 0), TN - 1) * TN + std::min(std::max(b, 0), TN - 1)) * AC_L + l]; };
               float *t = &m.tab[((size_t)iz * TN + ir) * AC_NF];
               for (int l = 0; l < AC_L; ++l)
               {
                  // z mirrors about the centre (the coil is symmetric), r is zero on the axis
                  t[2 * l] = iz == 0 ? 0.f : (float)((at(iz + 1, ir, l) - at(iz - 1, ir, l)) / (2 * TH));
                  t[2 * l + 1] = ir == 0 ? 0.f : (float)((at(iz, ir + 1, l) - at(iz, ir - 1, l)) / (2 * TH));
               }
               t[2 * AC_L] = (float)BZ[(size_t)iz * TN + ir];
               t[2 * AC_L + 1] = (float)BR[(size_t)iz * TN + ir];
            }
      }

      void Derive(Magnet &m)
      {
         const double nl = std::sqrt(m.ax * m.ax + m.ay * m.ay + m.az * m.az);
         m.nx = m.ax / nl; m.ny = m.ay / nl; m.nz = m.az / nl;
         // u = t x n for a helper vector t not parallel to n, v = n x u
         const double tx = std::fabs(m.nx) < 0.9 ? 1 : 0, ty = tx ? 0 : 1;
         m.ux = ty * m.nz; m.uy = -tx * m.nz; m.uz = tx * m.ny - ty * m.nx;
         const double ul = std::sqrt(m.ux * m.ux + m.uy * m.uy + m.uz * m.uz);
         m.ux /= ul; m.uy /= ul; m.uz /= ul;
         m.vx = m.ny * m.uz - m.nz * m.uy; m.vy = m.nz * m.ux - m.nx * m.uz; m.vz = m.nx * m.uy - m.ny * m.ux;
         if (m.aircore)
         {
            // layers of wire on a bore: turns per layer = length / wire diameter; the winding is three layers of equal turns at the thirds of its build
            const double len = 1e-3 * m.coilLenMm, rIn = 0.5e-3 * m.boreMm, dw = 1e-3 * m.wireMm;
            const double perLayer = std::max(1.0, std::floor(len / dw));
            const double build = std::ceil(m.turns / perLayer) * dw;
            for (int k = 0; k < 3; ++k)
               m.solRadius[k] = rIn + build * (k + 0.5) / 3;
            const double am = rIn + 0.5 * build;
            m.solLen = len;
            m.windIn = rIn;
            m.windOut = rIn + build;
            m.r = m.rOhm > 0 ? m.rOhm : 1.72e-8 * m.turns * 2 * PI * am / (PI * 0.25 * dw * dw * 0.8);
            m.l = MU0 * PI * am * am * m.turns * m.turns / (len + 0.9 * am);   // air-core coil (Wheeler type estimate)
            m.tau = m.l / m.r;
            constexpr double UM = 0.53975e-3;   // m per table unit
            m.tubeRU = std::max(0.0, 0.5e-3 * (m.boreMm - 2 * m.tubeWallMm) / UM);
            m.tubeHalfU = 0.5 * len / UM;
            m.housingRU = (rIn + build + 1e-3 * m.tubeWallMm) / UM;
            m.housingHalfU = (0.5 * len + 1e-3 * m.tubeWallMm) / UM;
            m.stopU = (m.stopPosMm > 1e8 ? 0.5 * len : 1e-3 * m.stopPosMm) / UM;
            const auto t0 = std::chrono::steady_clock::now();
            BuildAircoreTable(m);
            PLOGI << "Magnet " << m.name << ": air-core field table (ball response to order " << AC_L << ") built in " << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() << " ms";
            return;
         }
         m.radius = 0.5e-3 * m.coreDiaMm;
         m.length = 1e-3 * m.coreLenMm;
         const double area = PI * m.radius * m.radius;
         // reluctance of a disc pole radiating into space: 1 / (8 mu0 a) (the "capacitance of a disc" analogy); the far pole has the same, less the iron return
         m.spread = 8 * MU0 * m.radius;
         m.area = area;
         if (m.returnDistMm > 0)
         {
            const int cnt = std::max(1, std::min(4, (int)std::lround(m.returnCount)));
            const double d = 1e-3 * m.returnDistMm, w = 1e-3 * m.returnWidthMm, len = 1e-3 * m.returnLengthMm;
            // e1 = v, e2 = n x v: +x and +y for an axis pointing up
            const double e2x = m.ny * m.vz - m.nz * m.vy, e2y = m.nz * m.vx - m.nx * m.vz, e2z = m.nx * m.vy - m.ny * m.vx;
            m.spreadY = 8 * MU0 * std::sqrt(cnt * len * w / PI);   // flanges of the total area radiate like a disc of that area
            m.nY = cnt * 4;
            for (int k = 0; k < cnt; ++k)
            {
               const double th = (m.returnAngleDeg + 360.0 * k / cnt) * PI / 180, c = std::cos(th), sn = std::sin(th);
               for (int j = 0; j < 4; ++j)
               {
                  const double along = (j - 1.5) / 4 * len, rad = d + 0.5 * w;
                  double *q = m.ypts[k * 4 + j];
                  q[0] = rad * (c * m.vx + sn * e2x) + along * (-sn * m.vx + c * e2x);
                  q[1] = rad * (c * m.vy + sn * e2y) + along * (-sn * m.vy + c * e2y);
                  q[2] = rad * (c * m.vz + sn * e2z) + along * (-sn * m.vz + c * e2z);
                  for (int a = 0; a < 3; ++a)
                     q[a] -= 1e-3 * m.returnDepthMm * (a == 0 ? m.nx : a == 1 ? m.ny : m.nz);
               }
            }
         }
         const double farR = m.returnDistMm > 0 ? 1 / m.spreadY : (1 - m.returnFraction) / m.spread;
         m.l = m.turns * m.turns / (farR + 1 / m.spread);   // unloaded
         if (m.rOhm > 0)
            m.r = m.rOhm;
         else
         {
            // layers of wire on the core: turns per layer = winding length / wire diameter; the mean turn sits half the winding build outside the core
            const double dw = 1e-3 * m.wireMm;
            const double perLayer = std::max(1.0, std::floor(m.length / dw));
            const double build = std::ceil(m.turns / perLayer) * dw;
            const double meanTurn = 2 * PI * (m.radius + 0.5 * build);
            m.r = 1.72e-8 * m.turns * meanTurn / (PI * 0.25 * dw * dw * 0.8);   // copper; 0.8 = copper share of the insulated wire's cross section
         }
         m.tau = m.l / m.r;
         m.core.reset();
         m.drive.reset();
         if (m.cells)
         {
            // core_model = cells (magnetcore.h): winding on a 1 mm bobbin round the core, its build from the wire and the turns (as the resistance above);
            // the core reaches core_exposed_mm above the winding to the pole face. Tables per operating level (operating_volts, default rated_volts).
            const double dw = 1e-3 * m.wireMm, perLayer = std::max(1.0, std::floor(m.length / dw)), build = std::ceil(m.turns / perLayer) * dw;
            MagnetCore::Geometry g;
            g.coreRadius = m.radius; g.windLen = m.length; g.exposed = 1e-3 * m.coreExposedMm;
            g.windIn = m.radius + 1e-3; g.windOut = g.windIn + build; g.turns = m.turns; g.bSat = m.bSatT;
            std::vector<double> amps;
            for (const double v : (m.operatingVolts.empty() ? std::vector<double> { m.ratedVolts } : m.operatingVolts))
               amps.push_back(v / m.r);
            auto model = std::make_shared<MagnetCore::Model>();
            if (model->Build(g, amps))
            {
               m.core = model;
               PLOGI << "Magnet " << m.name << ": core_model cells, winding " << 1e3 * g.windIn << ".." << 1e3 * g.windOut << " mm x " << 1e3 * g.windLen << " mm, core exposed " << m.coreExposedMm << " mm\n" << model->log;
            }
            else
            {
               const std::string msg = "Magnet " + m.name + ": core tables could not be built (" + model->log + "); this magnet uses core_model = circuit.";
               PLOGW << msg;
               if (g_pplayer && g_pplayer->m_liveUI)
                  g_pplayer->m_liveUI->PushNotification(msg, 8000);
            }
            if (m.returnDistMm > 0 || m.returnFraction != 0)
            {
               const std::string msg = "Magnet " + m.name + ": return_* (yoke) keys are used by core_model = circuit only; core_model = cells ignores them (a bare rod core).";
               PLOGW << msg;
               if (g_pplayer && g_pplayer->m_liveUI)
                  g_pplayer->m_liveUI->PushNotification(msg, 8000);
            }
            // Option A ball-ball model (magnetballs.h): drive tables built from the core tables; any axis (the drive is rotated from the magnet's frame)
            if (m.core && m.ballInteraction)
            {
               auto d = std::make_shared<MagnetBalls::MagnetDrive>();
               if (d->Build(*m.core))
               {
                  m.drive = d;
                  PLOGI << "Magnet " << m.name << ": ball-ball model Option A, " << d->log;
               }
            }
         }
         for (int i = 0; i < NPTS; ++i)   // uniform over the disc area (sunflower pattern)
         {
            const double rho = m.radius * std::sqrt((i + 0.5) / NPTS), ph = i * 2.399963229728653;
            const double a = rho * std::cos(ph), b = rho * std::sin(ph);
            m.pts[i][0] = a * m.ux + b * m.vx;
            m.pts[i][1] = a * m.uy + b * m.vy;
            m.pts[i][2] = a * m.uz + b * m.vz;
         }
      }

      void Load()
      {
         g_loaded = true;
         g_magnets.clear();
         g_optAWarn = OptAWarnings();
         g_optASolver.ClearLast();
         g_parts.clear();
         g_partStart.clear();
         g_byName.clear();
         g_anyCurrent = false;
         if (g_pplayer == nullptr || g_pplayer->m_ptable == nullptr)
            return;
         const std::filesystem::path file = PhysicsIni::FindUnified(g_pplayer->m_ptable->m_filename);
         g_file = file;
         std::error_code tec;
         if (!file.empty())
            g_fileTime = std::filesystem::last_write_time(file, tec);
         if (file.empty())
         {
            PLOGI << "Magnets: no vpx_physics.ini found next to " << g_pplayer->m_ptable->m_filename.string() << " (and VPX_PHYSICS_INI is not set), no magnets";
            return;
         }
         if (!PhysicsIni::Active(file).magnets)
         {
            PLOGI << "Magnets: [Features] magnets = 0 in " << file.string() << " or the stock engine runs (no unified = 1 in the physics ini), no magnets";
            return;
         }
         std::ifstream in(file);
         std::string line;
         Magnet *cur = nullptr;
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
               if (PhysicsIni::Lower(head).rfind("magnet:", 0) == 0)
               {
                  g_magnets.emplace_back();
                  cur = &g_magnets.back();
                  cur->name = PhysicsIni::Trim(head.substr(7));
               }
               continue;
            }
            const size_t eq = line.find('=');
            if (cur == nullptr || eq == std::string::npos)
               continue;
            const std::string key = PhysicsIni::Lower(PhysicsIni::Trim(line.substr(0, eq)));
            const double v = atof(PhysicsIni::Trim(line.substr(eq + 1)).c_str());
            if (key == "x") cur->x = v;
            else if (key == "y") cur->y = v;
            else if (key == "z") cur->z = v;
            else if (key == "axis_x") cur->ax = v;
            else if (key == "axis_y") cur->ay = v;
            else if (key == "axis_z") cur->az = v;
            else if (key == "type") { const std::string t = PhysicsIni::Lower(PhysicsIni::Trim(line.substr(eq + 1))); cur->aircore = t == "aircore" || t == "solenoid"; }   // "solenoid" is the old name
            else if (key == "ball_interaction") cur->ballInteraction = v != 0;
            else if (key == "tube") cur->tube = v != 0;
            else if (key == "tube_visible") cur->tubeVisible = v != 0;
            else if (key == "tube_wall_mm") cur->tubeWallMm = v;
            else if (key == "tube_restitution") cur->tubeE = v;
            else if (key == "tube_friction") cur->tubeMu = v;
            else if (key == "stop") cur->stop = v != 0;
            else if (key == "stop_position_mm") cur->stopPosMm = v;
            else if (key == "coil_length_mm") cur->coilLenMm = v;
            else if (key == "bore_diameter_mm") cur->boreMm = v;
            else if (key == "turns") cur->turns = v;
            else if (key == "wire_mm") cur->wireMm = v;
            else if (key == "resistance_ohm") cur->rOhm = v;
            else if (key == "core_diameter_mm") cur->coreDiaMm = v;
            else if (key == "core_length_mm") cur->coreLenMm = v;
            else if (key == "return_fraction") cur->returnFraction = v;
            else if (key == "return_pole_distance_mm") cur->returnDistMm = v;
            else if (key == "return_pole_depth_mm") cur->returnDepthMm = v;
            else if (key == "return_pole_width_mm") cur->returnWidthMm = v;
            else if (key == "return_pole_length_mm") cur->returnLengthMm = v;
            else if (key == "return_pole_count") cur->returnCount = v;
            else if (key == "return_pole_angle_deg") cur->returnAngleDeg = v;
            else if (key == "b_sat_t") cur->bSatT = v;
            else if (key == "ball_ms_a_per_m") cur->msAm = v;
            else if (key == "ball_conductivity_s_per_m") cur->eddySigma = v;
            else if (key == "rated_volts" || key == "supply_volts") cur->ratedVolts = v;
            else if (key == "ball_eddy") cur->eddyScale = v;
            else if (key == "core_model") { const std::string t = PhysicsIni::Lower(PhysicsIni::Trim(line.substr(eq + 1))); cur->cells = t != "circuit"; if (t != "circuit" && t != "cells") PLOGW << "Magnet " << cur->name << ": core_model '" << t << "' unknown, using cells"; }
            else if (key == "core_exposed_mm") cur->coreExposedMm = v;
            else if (key == "ball_model") { const std::string msg = "Magnet " + cur->name + ": key ball_model is no longer used (the multipole ball-ball model was removed 2026-10-19); core_model = cells magnets use the ball-ball model Option A (2026-10-20), the others the induced dipoles."; PLOGW << msg; if (g_pplayer && g_pplayer->m_liveUI) g_pplayer->m_liveUI->PushNotification(msg, 8000); }
            else if (key == "operating_volts")
            {
               cur->operatingVolts.clear();
               std::string list = line.substr(eq + 1);
               for (char &ch : list) if (ch == ',') ch = ' ';
               std::istringstream ls(list); double vv;
               while (ls >> vv) cur->operatingVolts.push_back(vv);
            }
            else PLOGW << "Magnet " << cur->name << ": unknown key, ignored line '" << line << "'";
         }
         if (g_magnets.empty())
            PLOGI << "Magnets: " << file.string() << " has no [Magnet:<name>] section";
         for (size_t i = 0; i < g_magnets.size(); ++i)
         {
            Magnet &m = g_magnets[i];
            Derive(m);
            g_byName[PhysicsIni::Lower(m.name)] = i;
            const auto pend = g_pending.find(PhysicsIni::Lower(m.name));
            if (pend != g_pending.end())
               m.volts = pend->second.first ? std::max(-1.0, std::min(1.0, (double)pend->second.second)) * m.ratedVolts : pend->second.second;
            PLOGI << "Magnet " << m.name << (m.aircore ? " (air-core, coil centre at " + std::to_string(m.x) + ", " + std::to_string(m.y) + ", " + std::to_string(m.z) + ")" : std::string()) << ": R " << m.r << " ohm, L " << m.l * 1e3 << " mH, tau " << m.tau * 1e3 << " ms";
            if (!m.aircore)
               PLOGI << "Magnet " << m.name << ": " << m.bSatT << " T saturation, unloaded flux " << m.turns / ((m.returnDistMm > 0 ? 1 / m.spreadY : (1 - m.returnFraction) / m.spread) + 1.0 / m.spread) * 1e6 << " uWb per A";
         }
      }

      // A new play session (also a second Play from the editor, which can reuse the Player's address): the game clock starts again at 0, so a clock that went backwards
      // means everything loaded for the last session is stale. Called first by every entry point.
      void CheckSession()
      {
         static uint32_t lastMs = 0;
         if (g_pplayer == nullptr)
            return;
         const uint32_t now = g_pplayer->m_time_msec;
         if (g_owner != g_pplayer || now < lastMs)
         {
            g_loaded = false;
            g_owner = g_pplayer;
            g_magnets.clear();
            g_byName.clear();
            g_pending.clear();
            g_anyCurrent = false;
         }
         lastMs = now;
      }

      // The magnets load from the physics step: a script call made while the script itself is being loaded finds no table yet.
      void Ensure()
      {
         CheckSession();
         if (g_pplayer == nullptr || g_pplayer->m_ptable == nullptr)
            return;
         if (!g_loaded)
         {
            g_owner = g_pplayer;
            Load();
            g_pending.clear();
            g_nextCheckMs = g_pplayer->m_time_msec + 1000;
         }
         else if (!g_file.empty() && g_pplayer->m_time_msec >= g_nextCheckMs)
         {
            // like the other settings: re-read the ini when it has been saved (checked once a second of game time); a coil keeps its voltage and current
            g_nextCheckMs = g_pplayer->m_time_msec + 1000;
            std::error_code ec;
            const auto t = std::filesystem::last_write_time(g_file, ec);
            if (!ec && t != g_fileTime)
            {
               std::map<std::string, std::pair<double, double>> keep;
               for (const Magnet &m : g_magnets)
                  keep[PhysicsIni::Lower(m.name)] = { m.volts, m.amps };
               PLOGI << "Magnets: " << g_file.string() << " changed, reloading";
               Load();
               for (Magnet &m : g_magnets)
               {
                  const auto it = keep.find(PhysicsIni::Lower(m.name));
                  if (it != keep.end())
                  {
                     m.volts = it->second.first;
                     m.amps = it->second.second;
                  }
               }
            }
         }
      }
   }

   namespace
   {
      // The coil field (T) of every energised magnet at p (table units) for a ball of radius radiusU: the same fields as the force, without the gradient.
      void ExternalField(const double p[3], double radiusU, double B[3])
      {
         B[0] = B[1] = B[2] = 0.0;
         for (const Magnet &m : g_magnets)
         {
            if (m.amps == 0)
               continue;
            const double px = VPUTOM(p[0] - m.x), py = VPUTOM(p[1] - m.y), pz = VPUTOM(p[2] - m.z);
            if (!m.core && px * px + py * py + pz * pz > REACH * REACH)
               continue;
            if (m.aircore)
            {
               const double zc = px * m.nx + py * m.ny + pz * m.nz;
               const double rx = px - zc * m.nx, ry = py - zc * m.ny, rz = pz - zc * m.nz, rr = std::sqrt(rx * rx + ry * ry + rz * rz);
               const double fz = std::fabs(zc) / TH, fr = rr / TH;
               if (fz >= TN - 1 || fr >= TN - 1)
                  continue;
               const int iz = (int)fz, ir = (int)fr;
               const double wz = fz - iz, wr = fr - ir;
               double v[2];
               for (int q = 0; q < 2; ++q)
               {
                  const float *t00 = &m.tab[((size_t)iz * TN + ir) * AC_NF + 2 * AC_L + q], *t01 = t00 + AC_NF, *t10 = t00 + (size_t)TN * AC_NF, *t11 = t10 + AC_NF;
                  v[q] = (1 - wz) * ((1 - wr) * *t00 + wr * *t01) + wz * ((1 - wr) * *t10 + wr * *t11);
               }
               const double bz = v[0] * m.amps, br = v[1] * m.amps * (zc < 0 ? -1.0 : 1.0), ur = rr > 1e-9 ? 1.0 / rr : 0.0;   // Br is odd about the coil centre
               B[0] += bz * m.nx + br * rx * ur; B[1] += bz * m.ny + br * ry * ur; B[2] += bz * m.nz + br * rz * ur;
               continue;
            }
            if (m.core)
            {
               // core_model = cells: the bare field (coil + coil-only core) from the table; behind the pole face plane, nothing
               const double zc = px * m.nx + py * m.ny + pz * m.nz;
               const double rx = px - zc * m.nx, ry = py - zc * m.ny, rz = pz - zc * m.nz, rr = std::sqrt(rx * rx + ry * ry + rz * rz);
               double br, bz;
               if (zc <= 0 || !m.core->Field(rr, zc, m.amps, br, bz))
                  continue;
               const double ur = rr > 1e-9 ? 1.0 / rr : 0.0;
               B[0] += bz * m.nx + br * rx * ur; B[1] += bz * m.ny + br * ry * ur; B[2] += bz * m.nz + br * rz * ur;
               continue;
            }
            const double gap = std::max(0.5e-3, std::sqrt(px * px + py * py + pz * pz) - VPUTOM(radiusU));
            const double farR = m.returnDistMm > 0 ? 1.0 / m.spreadY : (1 - m.returnFraction) / m.spread;
            const double reluctance = farR + 1.0 / (m.spread + MU0 * m.area / gap);
            const double bc = std::min(m.bSatT, m.turns * std::fabs(m.amps) / reluctance / m.area);
            const double q = (m.amps < 0 ? -1.0 : 1.0) * bc * PI * m.radius * m.radius / NPTS / (4 * PI);
            auto add = [&](double qq, double sx, double sy, double sz)
            {
               const double rx = px - sx, ry = py - sy, rz = pz - sz;
               const double r2 = rx * rx + ry * ry + rz * rz, ir3 = qq / (r2 * std::sqrt(r2));
               B[0] += rx * ir3; B[1] += ry * ir3; B[2] += rz * ir3;
            };
            for (const auto &pt : m.pts)
               add(q, pt[0], pt[1], pt[2]);
            if (m.returnDistMm > 0)
               for (int i = 0; i < m.nY; ++i)
                  add(-NPTS * q / m.nY, m.ypts[i][0], m.ypts[i][1], m.ypts[i][2]);
            else
               add(-NPTS * q, -m.length * m.nx, -m.length * m.ny, -m.length * m.nz);
         }
      }

      // ---- ball-ball forces (2026-10-20). A ball takes part (is gated) when (a) the energised magnets' field at its centre is at least MagnetCore::B_GATE (0.665 mT:
      // a touching pair in that field pulls ~8 mN, ~1 % of the ball's weight, REPORT.md "GATE"), or (b) its surface is within half a ball diameter of a gated ball and
      // the field at its centre from the magnets plus the balls of the last ball-ball solve (Option A only) is at least B_GATE; (b) repeats until no ball is added.
      // Every gated magnetic ball in a field is an induced dipole, and its field acts on the other balls: in line
      // with the field they pull together, side by side across it they push apart. The moments follow the total field (coils + the neighbours' dipoles), solved
      // by three sweeps; the force on each pair is the dipole-dipole force, equal and opposite. Stored per ball as this step's velocity change (m_magPair).
      using MagnetCore::B_GATE;
      void UpdateBallInteraction()
      {
         struct BallM
         {
            HitBall *hb;
            double p[3], r, vol, mass, Bext[3], m[3];
            bool gated = false;
         };
         static std::vector<BallM> bs;
         bs.clear();
         g_parts.clear();
         g_partStart.clear();
         for (const auto ball : g_pplayer->m_vball)
         {
            HitBall &hb = ball->m_hitBall;
            hb.m_magPair[0] = hb.m_magPair[1] = hb.m_magPair[2] = 0.f;
            if (!hb.m_magnetic || hb.m_d.m_lockedInKicker)
               continue;
            BallM b = {};
            b.hb = &hb;
            b.p[0] = hb.m_d.m_pos.x; b.p[1] = hb.m_d.m_pos.y; b.p[2] = hb.m_d.m_pos.z;
            b.r = VPUTOM(hb.m_d.m_radius);
            b.vol = 4.0 / 3.0 * PI * b.r * b.r * b.r;
            b.mass = 7800.0 * b.vol * std::max(0.001f, hb.m_d.m_mass);
            bs.push_back(b);
         }
         const size_t n = bs.size();
         bool pairOn = false;
         for (const Magnet &m : g_magnets)
            pairOn = pairOn || (m.amps != 0 && m.ballInteraction);
         if (n < 2 || !pairOn)
         {
            g_optASolver.ClearLast();
            return;
         }
         int ngated = 0;
         for (BallM &b : bs)   // gate (a)
         {
            ExternalField(b.p, b.hb->m_d.m_radius, b.Bext);
            b.gated = std::sqrt(b.Bext[0] * b.Bext[0] + b.Bext[1] * b.Bext[1] + b.Bext[2] * b.Bext[2]) >= B_GATE;
            ngated += b.gated ? 1 : 0;
         }
         for (bool added = ngated > 0; added;)   // gate (b): the balls' own fields carry the magnetisation on through a chain (one cycle late: the last solve)
         {
            added = false;
            for (BallM &b : bs)
            {
               if (b.gated)
                  continue;
               bool touchesGated = false;
               for (const BallM &o : bs)
                  if (o.gated)
                  {
                     const double dx = VPUTOM(b.p[0] - o.p[0]), dy = VPUTOM(b.p[1] - o.p[1]), dz = VPUTOM(b.p[2] - o.p[2]);
                     touchesGated = touchesGated || std::sqrt(dx * dx + dy * dy + dz * dz) - b.r - o.r <= b.r;
                  }
               MagnetBalls::V3 Bb;
               if (!touchesGated || !g_optASolver.FieldFromBalls({ VPUTOM(b.p[0]), VPUTOM(b.p[1]), VPUTOM(b.p[2]) }, (uint64_t)(uintptr_t)b.hb, Bb))
                  continue;
               const double tx = b.Bext[0] + Bb.x, ty = b.Bext[1] + Bb.y, tz = b.Bext[2] + Bb.z;
               if (std::sqrt(tx * tx + ty * ty + tz * tz) >= B_GATE)
               {
                  b.gated = true;
                  ++ngated;
                  added = true;
               }
            }
         }
         if (ngated < 2)
         {
            g_optASolver.ClearLast();
            return;
         }
         const double K = MS2TOVPUVPT2(1.0f) * PHYS_FACTOR;

         // Option A (magnetballs.h) when every energised ball-ball magnet has its drive tables; otherwise the induced dipoles below
         {
            bool allA = true, anyA = false;
            for (const Magnet &m : g_magnets)
               if (m.amps != 0 && m.ballInteraction)
               {
                  allA = allA && m.drive != nullptr;
                  anyA = anyA || m.drive != nullptr;
               }
            if (anyA && !allA && !g_optAWarn.mixed)
            {
               g_optAWarn.mixed = true;
               const std::string msg = "Magnets: some energised magnets use the ball-ball model Option A and some the induced dipoles (air-core or core_model = circuit); while they are on together all balls use the dipoles.";
               PLOGW << msg;
               if (g_pplayer && g_pplayer->m_liveUI)
                  g_pplayer->m_liveUI->PushNotification(msg, 8000);
            }
            if (allA && anyA)
            {
               MagnetBalls::Solver &solver = g_optASolver;
               std::vector<MagnetBalls::Source> src;
               std::vector<int> srcMag;   // source -> magnet index
               for (const Magnet &m : g_magnets)
                  if (m.amps != 0 && m.ballInteraction && m.drive)
                  {
                     srcMag.push_back((int)(&m - &g_magnets[0]));
                     MagnetBalls::Source s;
                     s.drive = m.drive.get();
                     s.amps = m.amps;
                     s.origin = { VPUTOM(m.x), VPUTOM(m.y), VPUTOM(m.z) };
                     s.axis = { m.nx, m.ny, m.nz };
                     src.push_back(s);
                  }
               std::vector<int> gidx;
               std::vector<MagnetBalls::V3> pos, dF;
               std::vector<uint64_t> ids;
               for (size_t i = 0; i < n; ++i)
                  if (bs[i].gated)
                  {
                     gidx.push_back((int)i);
                     pos.push_back({ VPUTOM(bs[i].p[0]), VPUTOM(bs[i].p[1]), VPUTOM(bs[i].p[2]) });
                     ids.push_back((uint64_t)(uintptr_t)bs[i].hb);
                  }
               MagnetBalls::Stats st;
               const auto t0 = std::chrono::steady_clock::now();
               std::vector<MagnetBalls::Part> parts;
               solver.Evaluate(src, pos, dF, &st, &ids, &parts);
               const double us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
               for (size_t k = 0; k < gidx.size(); ++k)
               {
                  BallM &b = bs[gidx[k]];
                  b.hb->m_magPair[0] = (float)(dF[k].x / b.mass * K);
                  b.hb->m_magPair[1] = (float)(dF[k].y / b.mass * K);
                  b.hb->m_magPair[2] = (float)(dF[k].z / b.mass * K);
                  g_partStart[b.hb] = { pos[k].x, pos[k].y, pos[k].z };
                  g_parts[b.hb].massKg = b.mass;
               }
               for (const MagnetBalls::Part &pt : parts)   // keep the parts with their reference: the pair's line, or the ball's direction round the magnet axis
               {
                  BallBallPart sp;
                  sp.a = bs[gidx[pt.a]].hb; sp.F[0] = pt.F.x; sp.F[1] = pt.F.y; sp.F[2] = pt.F.z;
                  sp.pa0[0] = pos[pt.a].x; sp.pa0[1] = pos[pt.a].y; sp.pa0[2] = pos[pt.a].z;
                  if (pt.b >= 0)
                  {
                     sp.b = bs[gidx[pt.b]].hb;
                     sp.pb0[0] = pos[pt.b].x; sp.pb0[1] = pos[pt.b].y; sp.pb0[2] = pos[pt.b].z;
                     const double d[3] = { pos[pt.a].x - pos[pt.b].x, pos[pt.a].y - pos[pt.b].y, pos[pt.a].z - pos[pt.b].z }, l = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                     for (int c = 0; c < 3; c++) sp.u0[c] = l > 0 ? d[c] / l : 0.0;
                  }
                  else
                  {
                     const Magnet &m = g_magnets[srcMag[pt.src]];
                     sp.n[0] = m.nx; sp.n[1] = m.ny; sp.n[2] = m.nz; sp.o[0] = VPUTOM(m.x); sp.o[1] = VPUTOM(m.y); sp.o[2] = VPUTOM(m.z);
                     double d[3] = { pos[pt.a].x - sp.o[0], pos[pt.a].y - sp.o[1], pos[pt.a].z - sp.o[2] }; const double dn = d[0] * m.nx + d[1] * m.ny + d[2] * m.nz;
                     for (int c = 0; c < 3; c++) d[c] -= dn * sp.n[c];
                     const double l = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                     if (l > 1e-6) for (int c = 0; c < 3; c++) sp.r0[c] = d[c] / l;
                     else sp.n[0] = sp.n[1] = sp.n[2] = 0;   // a ball on the axis: no direction round it, the part is held
                  }
                  g_parts[bs[gidx[pt.ball]].hb].parts.push_back(sp);
               }
               if (us > 500.0 && !g_optAWarn.cost)
               {
                  g_optAWarn.cost = true;
                  std::ostringstream os;
                  os << "Magnets: the ball-ball model took " << (int)us << " us in one physics step (" << st.balls << " balls, " << st.pairFix << " pairs within the pair-fix range); its budget is 500 us.";
                  PLOGW << os.str();
                  if (g_pplayer && g_pplayer->m_liveUI)
                     g_pplayer->m_liveUI->PushNotification(os.str(), 8000);
               }
               return;
            }
         }

         // induced dipoles between the gated balls
         g_optASolver.ClearLast();
         std::vector<int> gi;
         for (size_t i = 0; i < n; ++i)
            if (bs[i].gated)
               gi.push_back((int)i);
         const double Ms = g_magnets.front().msAm;
         for (const int i : gi)
         {
            BallM &b = bs[i];
            const double bm = std::sqrt(b.Bext[0] * b.Bext[0] + b.Bext[1] * b.Bext[1] + b.Bext[2] * b.Bext[2]);
            const double M = bm > 1e-12 ? std::min(3 * bm / MU0, Ms) * b.vol / bm : 0.0;
            for (int a = 0; a < 3; ++a)
               b.m[a] = M * b.Bext[a];
         }
         auto rel = [&](const BallM &a, const BallM &b, double d[3])   // from b to a, m
         {
            for (int k = 0; k < 3; ++k)
               d[k] = VPUTOM(a.p[k] - b.p[k]);
         };
         for (int sweep = 0; sweep < 3; ++sweep)
            for (const int i : gi)
            {
               double B[3] = { bs[i].Bext[0], bs[i].Bext[1], bs[i].Bext[2] };
               for (const int j : gi)
               {
                  if (j == i)
                     continue;
                  double d[3];
                  rel(bs[i], bs[j], d);
                  const double dist = std::max(std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]), bs[i].r + bs[j].r);
                  const double mr = (bs[j].m[0] * d[0] + bs[j].m[1] * d[1] + bs[j].m[2] * d[2]) / (dist * dist);   // m.r_hat / dist
                  const double c = MU0 / (4 * PI * dist * dist * dist);
                  for (int k = 0; k < 3; ++k)
                     B[k] += c * (3 * mr * d[k] - bs[j].m[k]);
               }
               const double bm = std::sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
               const double M = bm > 1e-12 ? std::min(3 * bm / MU0, Ms) * bs[i].vol / bm : 0.0;
               for (int k = 0; k < 3; ++k)
                  bs[i].m[k] = M * B[k];
            }
         for (size_t ii = 0; ii < gi.size(); ++ii)
            for (size_t jj = ii + 1; jj < gi.size(); ++jj)
            {
               const int i = gi[ii], j = gi[jj];
               double d[3];   // from i to j
               rel(bs[j], bs[i], d);
               const double dl = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
               if (dl < 1e-9)
                  continue;
               const double dist = std::max(dl, bs[i].r + bs[j].r);
               const double rh[3] = { d[0] / dl, d[1] / dl, d[2] / dl };
               const double* m1 = bs[i].m; const double* m2 = bs[j].m;
               const double m1r = m1[0] * rh[0] + m1[1] * rh[1] + m1[2] * rh[2], m2r = m2[0] * rh[0] + m2[1] * rh[1] + m2[2] * rh[2], m12 = m1[0] * m2[0] + m1[1] * m2[1] + m1[2] * m2[2];
               const double c = 3 * MU0 / (4 * PI * dist * dist * dist * dist);
               double Fj[3];   // force on j
               for (int k = 0; k < 3; ++k)
                  Fj[k] = c * (m1r * m2[k] + m2r * m1[k] + (m12 - 5 * m1r * m2r) * rh[k]);
               for (int k = 0; k < 3; ++k)
               {
                  bs[j].hb->m_magPair[k] += (float)(Fj[k] / bs[j].mass * K);
                  bs[i].hb->m_magPair[k] -= (float)(Fj[k] / bs[i].mass * K);
               }
            }
      }
   }

   void Step()
   {
      Ensure();
      if (g_magnets.empty())
         return;
      bool any = false;
      for (Magnet &m : g_magnets)
      {
         m.amps += (m.volts / m.r - m.amps) * (1.0 - std::exp(-PHYSICS_STEPTIME_S / m.tau));
         if (m.volts == 0 && std::fabs(m.amps) < 1e-6)
            m.amps = 0;
         any = any || m.amps != 0;
      }
      g_anyCurrent = any;
      UpdateBallInteraction();
   }

   // Draws the guide tubes (outer and inner wall), the winding housings and the stops as translucent surfaces, in table units, after the table's own dynamic parts.
   void RenderGuides(Renderer *renderer)
   {
      for (Magnet &m : g_magnets)
      {
         if (!m.aircore || !m.tubeVisible || m.tubeRU <= 0)
            continue;
         RenderDevice *rd = renderer->m_renderDevice;
         if (!m.mesh)
         {
            constexpr int SEG = 48;
            constexpr double UM = 0.53975e-3;
            const double rOut = 0.5e-3 * m.boreMm / UM, rIn = m.tubeRU, h = m.tubeHalfU;
            std::vector<Vertex3D_NoTex2> v;
            std::vector<unsigned int> idx;
            auto pt = [&](double rad, double ax, double ang, double nrad, double nax)
            {
               const double ca = std::cos(ang), sa = std::sin(ang);
               const double ex = ca * m.ux + sa * m.vx, ey = ca * m.uy + sa * m.vy, ez = ca * m.uz + sa * m.vz;
               Vertex3D_NoTex2 q = {};
               q.x = (float)(m.x + rad * ex + ax * m.nx); q.y = (float)(m.y + rad * ey + ax * m.ny); q.z = (float)(m.z + rad * ez + ax * m.nz);
               q.nx = (float)(nrad * ex + nax * m.nx); q.ny = (float)(nrad * ey + nax * m.ny); q.nz = (float)(nrad * ez + nax * m.nz);
               v.push_back(q);
            };
            auto cyl = [&](double rad, double nrad, double hl)   // a cylinder wall from -hl to +hl
            {
               const unsigned int base = (unsigned int)v.size();
               for (int k = 0; k <= SEG; ++k)
               {
                  pt(rad, -hl, 2 * PI * k / SEG, nrad, 0);
                  pt(rad, hl, 2 * PI * k / SEG, nrad, 0);
               }
               for (int k = 0; k < SEG; ++k)
               {
                  const unsigned int a = base + 2 * k;
                  idx.insert(idx.end(), { a, a + 1, a + 2, a + 1, a + 3, a + 2 });
               }
            };
            const double hh = m.tube ? m.housingHalfU : h;
            cyl(rIn, -1, hh);
            if (!m.tube)
               cyl(rOut, 1, h);
            else   // the winding housing: outer wall and the two end flanges from the tube's inside to the housing's outside
            {
               cyl(m.housingRU, 1, hh);
               for (const double sg : { -1.0, 1.0 })
               {
                  const unsigned int base = (unsigned int)v.size();
                  for (int k = 0; k <= SEG; ++k)
                  {
                     pt(rIn, sg * hh, 2 * PI * k / SEG, 0, sg);
                     pt(m.housingRU, sg * hh, 2 * PI * k / SEG, 0, sg);
                  }
                  for (int k = 0; k < SEG; ++k)
                  {
                     const unsigned int a = base + 2 * k;
                     idx.insert(idx.end(), { a, a + 1, a + 2, a + 1, a + 3, a + 2 });
                  }
               }
            }
            if (m.stop)   // the stop: a disc across the bore
            {
               const unsigned int base = (unsigned int)v.size();
               const double sg = m.stopU >= 0 ? 1.0 : -1.0;
               pt(0, m.stopU, 0, 0, sg);
               for (int k = 0; k <= SEG; ++k)
                  pt(rIn, m.stopU, 2 * PI * k / SEG, 0, sg);
               for (int k = 0; k < SEG; ++k)
                  idx.insert(idx.end(), { base, base + 1 + k, base + 2 + k });
            }
            auto vb = std::make_shared<VertexBuffer>(rd, (unsigned int)v.size(), (const float *)v.data());
            auto ib = std::make_shared<IndexBuffer>(rd, idx);
            m.mesh = std::make_shared<MeshBuffer>("Magnet guide"s, vb, ib, true);
         }
         rd->ResetRenderState();
         rd->SetRenderState(RenderState::ZENABLE, RenderState::RS_TRUE);
         rd->SetRenderState(RenderState::ZWRITEENABLE, RenderState::RS_FALSE);
         rd->SetRenderState(RenderState::CULLMODE, RenderState::CULL_NONE);
         rd->SetRenderState(RenderState::ALPHABLENDENABLE, RenderState::RS_TRUE);
         rd->m_basicShader->SetTechnique(ShaderTechnique::unshaded_without_texture);
         rd->m_basicShader->SetVector(ShaderUniform::staticColor_Alpha, 0.55f, 0.75f, 1.0f, 0.25f);
         rd->DrawMesh(rd->m_basicShader, true, Vertex3Ds((float)m.x, (float)m.y, (float)m.z), 0.f, m.mesh, RenderDevice::TRIANGLELIST, 0, m.mesh->m_ib->m_count);
         rd->m_basicShader->SetVector(ShaderUniform::staticColor_Alpha, 1.f, 1.f, 1.f, 1.f);
      }
   }

   void ReleaseRender()
   {
      for (Magnet &m : g_magnets)
         m.mesh.reset();   // while the render device is still alive: a mesh destroyed with the magnet list of the next session would touch a freed device
   }

   bool HasGuides()
   {
      for (const Magnet &m : g_magnets)
         if (m.aircore && (m.tube || m.stop) && m.tubeRU > 0)
            return true;
      return false;
   }

   // Impulse contact of the ball with the tube, the winding housing and the stop (table units, U and VPT; a sphere of inertia 2/5 m R^2). The tube and the housing
   // are one solid (a bobbin): in the coil's meridian plane the block r = tube inside .. housing outside, |z| <= the housing half length; the ball touches it where the distance from
   // its centre to the block is under its radius (inner wall, outer wall, end rings and their edges alike). Positions are projected out of the wall,
   // the normal velocity is reflected with the restitution (none below 0.05 m/s, so a ball pressed on the wall by the coil's pull does not chatter), and the contact
   // point's sliding velocity is removed up to the Coulomb limit, which also spins the ball up to rolling.
   bool ConstrainBall(float pos[3], float vel[3], float angmom[3], float radiusU, float mass, float inertia)
   {
      bool changed = false;
      for (const Magnet &m : g_magnets)
      {
         if (!m.aircore || !(m.tube || m.stop) || m.tubeRU <= 0)
            continue;
         const double R = radiusU;
         const double px = pos[0] - m.x, py = pos[1] - m.y, pz = pos[2] - m.z;
         const double zc = px * m.nx + py * m.ny + pz * m.nz;
         if (std::fabs(zc) > m.housingHalfU + R && std::fabs(zc) > std::fabs(m.stopU) + R)
            continue;
         const double rx = px - zc * m.nx, ry = py - zc * m.ny, rz = pz - zc * m.nz, rr = std::sqrt(rx * rx + ry * ry + rz * rz);
         if (rr > m.housingRU + R)
            continue;
         // contact with the unit normal (cx, cy, cz) pointing from the ball to the wall, penetration pen
         auto contact = [&](double cx, double cy, double cz, double pen)
         {
            pos[0] -= (float)(cx * pen); pos[1] -= (float)(cy * pen); pos[2] -= (float)(cz * pen);
            changed = true;
            const double vn = vel[0] * cx + vel[1] * cy + vel[2] * cz;
            if (vn <= 0)
               return;
            const double vMin = 0.05 / (0.53975e-3 / DEFAULT_STEPTIME_S);   // 0.05 m/s in U/VPT
            const double e = vn > vMin ? m.tubeE : 0.0;
            const double jn = (1 + e) * mass * vn;
            vel[0] -= (float)(jn / mass * cx); vel[1] -= (float)(jn / mass * cy); vel[2] -= (float)(jn / mass * cz);
            // sliding velocity of the contact point (at +R along the normal): v + w x (R n)
            const double wx = angmom[0] / inertia, wy = angmom[1] / inertia, wz = angmom[2] / inertia;
            double cvx = vel[0] + wy * R * cz - wz * R * cy, cvy = vel[1] + wz * R * cx - wx * R * cz, cvz = vel[2] + wx * R * cy - wy * R * cx;
            const double cn = cvx * cx + cvy * cy + cvz * cz;
            cvx -= cn * cx; cvy -= cn * cy; cvz -= cn * cz;
            const double vt = std::sqrt(cvx * cvx + cvy * cvy + cvz * cvz);
            if (vt < 1e-12)
               return;
            const double mt = mass / (1.0 + mass * R * R / inertia);   // effective mass of the contact point against a tangential impulse
            const double jt = std::min(mt * vt, m.tubeMu * jn);
            const double tx = -cvx / vt * jt, ty = -cvy / vt * jt, tz = -cvz / vt * jt;   // impulse on the ball
            vel[0] += (float)(tx / mass); vel[1] += (float)(ty / mass); vel[2] += (float)(tz / mass);
            angmom[0] += (float)(R * (cy * tz - cz * ty)); angmom[1] += (float)(R * (cz * tx - cx * tz)); angmom[2] += (float)(R * (cx * ty - cy * tx));
         };
         if (m.tube && rr > 1e-9)
         {
            const double ex = rx / rr, ey = ry / rr, ez = rz / rr;   // radial unit vector
            const double h = m.housingHalfU, a = m.tubeRU, b = m.housingRU;
            const double az = std::fabs(zc), sz = zc < 0 ? -1.0 : 1.0;
            if (az < h && rr > a && rr < b)
            {
               // centre inside the block (a fast ball, or one placed there): out through the nearest face
               const double dIn = rr - a, dOut = b - rr, dEnd = h - az;
               if (dIn <= dOut && dIn <= dEnd) contact(ex, ey, ez, dIn + R);
               else if (dOut <= dEnd) contact(-ex, -ey, -ez, dOut + R);
               else contact(-sz * m.nx, -sz * m.ny, -sz * m.nz, dEnd + R);
            }
            else
            {
               const double qz = std::min(az, h), qr = std::min(std::max(rr, a), b);   // nearest point of the block
               const double dz = qz - az, dr = qr - rr, d = std::sqrt(dz * dz + dr * dr);
               if (d < R && d > 1e-9)
               {
                  const double cz = sz * dz / d, cr = dr / d;
                  contact(cr * ex + cz * m.nx, cr * ey + cz * m.ny, cr * ez + cz * m.nz, R - d);
               }
            }
         }
         if (m.stop && rr < m.tubeRU)
         {
            const double sg = m.stopU >= 0 ? 1.0 : -1.0, d = sg * zc, st = std::fabs(m.stopU);
            if (d > st - R && d < st + R)
               contact(sg * m.nx, sg * m.ny, sg * m.nz, d - (st - R));
         }
      }
      return changed;
   }

   bool Active()
   {
      return g_anyCurrent;
   }

   bool HasBallBallParts(const HitBall *hb)
   {
      return g_parts.find(hb) != g_parts.end();
   }

   const std::vector<BallBallPart> *BallBallParts(const HitBall *hb, double &massKg)
   {
      const auto it = g_parts.find(hb);
      if (it == g_parts.end())
         return nullptr;
      massKg = it->second.massKg;
      return &it->second.parts;
   }

   bool FieldAt(const float pos[3], float radiusU, double B[3])
   {
      B[0] = B[1] = B[2] = 0.0;
      if (!g_anyCurrent)
         return false;
      const double p[3] = { pos[0], pos[1], pos[2] };
      ExternalField(p, radiusU, B);
      // the field INSIDE the steel ball, which is what its eddy currents see: B + 2/3 mu0 M, i.e. 3 B for a ball that is not saturated (a permeable sphere gathers the field)
      const double bm = std::sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
      if (bm > 1e-12)
      {
         const double f = 1.0 + 2.0 * std::min(1.0, MU0 * g_magnets.front().msAm / (3.0 * bm));
         B[0] *= f; B[1] *= f; B[2] *= f;
      }
      return true;
   }

   double EddyConductivity()
   {
      return g_magnets.empty() ? 0.0 : g_magnets.front().eddySigma * g_magnets.front().eddyScale;
   }

   void AddVelocityIncrement(const float pos[3], float radiusU, float mass, float dv[3])
   {
      if (!g_anyCurrent)
         return;
      for (const Magnet &m : g_magnets)
      {
         if (m.amps == 0)
            continue;
         const double px = VPUTOM(pos[0] - (float)m.x), py = VPUTOM(pos[1] - (float)m.y), pz = VPUTOM(pos[2] - (float)m.z);   // ball relative to the face centre, m
         if (!m.core && px * px + py * py + pz * pz > REACH * REACH)
            continue;
         if (m.aircore)
         {
            // Air-core coil: F = I^2 sum_l a^(2l+1) grad Q_l (the ball's energy to order 4, see BuildAircoreTable), bilinear in the cell; z mirrors about the centre.
            const double zc = px * m.nx + py * m.ny + pz * m.nz;
            const double rx = px - zc * m.nx, ry = py - zc * m.ny, rz = pz - zc * m.nz, rr = std::sqrt(rx * rx + ry * ry + rz * rz);
            const double fz = std::fabs(zc) / TH, fr = rr / TH;
            if (fz >= TN - 1 || fr >= TN - 1)
               continue;
            const int iz = (int)fz, ir = (int)fr;
            const double wz = fz - iz, wr = fr - ir;
            double v[AC_NF];
            for (int q = 0; q < AC_NF; ++q)
            {
               const float *t00 = &m.tab[((size_t)iz * TN + ir) * AC_NF + q], *t01 = t00 + AC_NF, *t10 = t00 + (size_t)TN * AC_NF, *t11 = t10 + AC_NF;
               v[q] = (1 - wz) * ((1 - wr) * *t00 + wr * *t01) + wz * ((1 - wr) * *t10 + wr * *t11);
            }
            const double a = VPUTOM(radiusU), i2 = m.amps * m.amps;
            // The series holds while the ball is clear of the winding. A ball inside the winding block (nothing solid stops it there) keeps the dipole term only.
            int orders = AC_L;
            if (WindingGap(m, zc, rr) < a)
            {
               orders = 1;
               if (!m.overlapWarned)
               {
                  m.overlapWarned = true;
                  const std::string msg = "Magnet " + m.name + ": a ball is inside the air-core winding (no tube or track keeps it out); its pull there uses the dipole term only. Give the coil a tube or the ball a track.";
                  PLOGW << msg;
                  if (g_pplayer && g_pplayer->m_liveUI)
                     g_pplayer->m_liveUI->PushNotification(msg, 8000);
               }
            }
            double gz = 0, gr = 0, a2l1 = a * a * a;
            for (int l = 0; l < orders; ++l, a2l1 *= a * a)
            {
               gz += a2l1 * v[2 * l];
               gr += a2l1 * v[2 * l + 1];
            }
            // ball saturation: the moment is capped at Ms (M = 3H below it); applied to the whole response, exact for the dipole term
            const double B = std::sqrt(v[2 * AC_L] * v[2 * AC_L] + v[2 * AC_L + 1] * v[2 * AC_L + 1]) * std::fabs(m.amps);
            const double sat = 3 * B > MU0 * m.msAm ? MU0 * m.msAm / (3 * B) : 1.0;
            const double s = i2 * sat / (7800.0 * 4.0 / 3.0 * PI * a * a * a * std::max(0.001f, mass)) * MS2TOVPUVPT2(1.0f) * PHYS_FACTOR;
            const double Fz = zc < 0 ? -gz : gz, ur = rr > 1e-9 ? 1.0 / rr : 0.0;
            dv[0] += (float)(s * (Fz * m.nx + gr * rx * ur));
            dv[1] += (float)(s * (Fz * m.ny + gr * ry * ur));
            dv[2] += (float)(s * (Fz * m.nz + gr * rz * ur));
            continue;
         }
         if (m.core)
         {
            // core_model = cells (magnetcore.h): the force table of the operating levels around the coil current: the order-4 ball in the bare field, the core's
            // answer (on-axis reaction carried by the tube equivalent height), the ball's saturation and the rim rule. Built for the standard ball (radius 25 U).
            const double zc = px * m.nx + py * m.ny + pz * m.nz;
            if (zc <= 0)
               continue;   // the ball centre behind the pole face plane: no table there
            const double rx = px - zc * m.nx, ry = py - zc * m.ny, rz = pz - zc * m.nz, rr = std::sqrt(rx * rx + ry * ry + rz * rz);
            const double a = VPUTOM(radiusU);
            if (std::fabs(a - m.core->g.ballR) > 0.02 * m.core->g.ballR && !m.radiusWarned)
            {
               m.radiusWarned = true;
               const std::string msg = "Magnet " + m.name + ": a ball of radius " + std::to_string(a * 1e3) + " mm; core_model = cells tables are built for the standard 13.49 mm ball, so its pull is approximate.";
               PLOGW << msg;
               if (g_pplayer && g_pplayer->m_liveUI)
                  g_pplayer->m_liveUI->PushNotification(msg, 8000);
            }
            double down, side;
            if (!m.core->Force(rr, zc - a, m.amps, down, side))
               continue;
            // down: toward the face (-n); side: toward the axis (-r); a = F / (7800 kg/m^3 * V * mass)
            const double s = 1.0 / (7800.0 * 4.0 / 3.0 * PI * a * a * a * std::max(0.001f, mass)) * MS2TOVPUVPT2(1.0f) * PHYS_FACTOR;
            const double ur = rr > 1e-9 ? side / rr : 0.0;
            dv[0] += (float)(-s * (down * m.nx + ur * rx));
            dv[1] += (float)(-s * (down * m.ny + ur * ry));
            dv[2] += (float)(-s * (down * m.nz + ur * rz));
            continue;
         }
         // magnetic circuit: far pole (open to space unless there is an iron return) in series with the near pole, whose conductance is the
         // spreading into space plus a parallel-plate term mu0 A / gap that takes over as the ball nears the face (the ball closes the circuit)
         const double gap = std::max(0.5e-3, std::sqrt(px * px + py * py + pz * pz) - VPUTOM(radiusU));
         const double farR = m.returnDistMm > 0 ? 1.0 / m.spreadY : (1 - m.returnFraction) / m.spread;
         const double reluctance = farR + 1.0 / (m.spread + MU0 * m.area / gap);
         const double bc = std::min(m.bSatT, m.turns * std::fabs(m.amps) / reluctance / m.area);
         const double q = (m.amps < 0 ? -1.0 : 1.0) * bc * PI * m.radius * m.radius / NPTS / (4 * PI);   // pole strength of one point / 4 pi
         double Bx = 0, By = 0, Bz = 0, G[6] = { 0, 0, 0, 0, 0, 0 };   // gradient tensor xx yy zz xy xz yz
         auto add = [&](double qq, double sx, double sy, double sz)
         {
            const double rx = px - sx, ry = py - sy, rz = pz - sz;
            const double r2 = rx * rx + ry * ry + rz * rz, ir3 = qq / (r2 * std::sqrt(r2)), ir5 = ir3 / r2;
            Bx += rx * ir3; By += ry * ir3; Bz += rz * ir3;
            G[0] += ir3 - 3 * rx * rx * ir5; G[1] += ir3 - 3 * ry * ry * ir5; G[2] += ir3 - 3 * rz * rz * ir5;
            G[3] -= 3 * rx * ry * ir5; G[4] -= 3 * rx * rz * ir5; G[5] -= 3 * ry * rz * ir5;
         };
         for (const auto &p : m.pts)
            add(q, p[0], p[1], p[2]);
         if (m.returnDistMm > 0)
            for (int i = 0; i < m.nY; ++i)
               add(-NPTS * q / m.nY, m.ypts[i][0], m.ypts[i][1], m.ypts[i][2]);   // the yoke's flange poles, in the face plane
         else
            add(-NPTS * q, -m.length * m.nx, -m.length * m.ny, -m.length * m.nz);   // the far face as one charge behind the near face
         const double b = std::sqrt(Bx * Bx + By * By + Bz * Bz);
         if (b < 1e-9)
            continue;
         // grad|B| = G B / |B|;  F = V M grad|B|, M = min(3 |B| / mu0, Ms);  a = F / (7800 kg/m^3 * V * mass)
         const double gx = (G[0] * Bx + G[3] * By + G[4] * Bz) / b, gy = (G[3] * Bx + G[1] * By + G[5] * Bz) / b, gz = (G[4] * Bx + G[5] * By + G[2] * Bz) / b;
         const double M = std::min(3 * b / MU0, m.msAm);
         const double a = M / (7800.0 * std::max(0.001f, mass));   // m/s^2 per (T/m)
         const double s = a * MS2TOVPUVPT2(1.0f) * PHYS_FACTOR;
         dv[0] += (float)(s * gx); dv[1] += (float)(s * gy); dv[2] += (float)(s * gz);
      }
   }

   bool SetVoltage(const char *name, float volts)
   {
      CheckSession();
      if (!g_loaded)
      {
         g_pending[PhysicsIni::Lower(name)] = { false, volts };   // applied when the physics step loads the magnets
         return true;
      }
      const auto it = g_byName.find(PhysicsIni::Lower(name));
      if (it == g_byName.end())
         return false;
      g_magnets[it->second].volts = volts;
      return true;
   }

   // PWM drive: the coil sees the mean voltage duty x rated volts (its L / R time constant, 5 - 30 ms, smooths a PWM of a kilohertz or more)
   bool SetDuty(const char *name, float duty)
   {
      CheckSession();
      if (!g_loaded)
      {
         g_pending[PhysicsIni::Lower(name)] = { true, duty };
         return true;
      }
      const auto it = g_byName.find(PhysicsIni::Lower(name));
      if (it == g_byName.end())
         return false;
      g_magnets[it->second].volts = std::max(-1.0, std::min(1.0, (double)duty)) * g_magnets[it->second].ratedVolts;
      return true;
   }

   bool GetCurrent(const char *name, float &amps)
   {
      amps = 0.f;
      CheckSession();
      if (!g_loaded)
         return false;
      const auto it = g_byName.find(PhysicsIni::Lower(name));
      if (it == g_byName.end())
         return false;
      amps = (float)g_magnets[it->second].amps;
      return true;
   }
}
