// license:GPLv3+
// Iron-core magnet: the core's answer to the ball, built at load (see magnetcore.h). Ported from the offline benches in runs/2026-10-14-mirror:
//   bench_bare3.cpp (cells, dense rings, axis law, smooth read-out), bench_force.cpp (order-4 ball on an inner sphere), bench_axial.cpp (on-axis cell + ball solve),
//   bench_geom.cpp (tube equivalent height with the visibility boundary term), face_neck.txt (ball saturation), rim_bound2.py (rim rule).

#include "physics/magnetcore.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <thread>
#include <sstream>

namespace MagnetCore
{
   namespace
   {
      constexpr double PI = 3.14159265358979323846, MU0 = 4e-7 * PI;
      constexpr double CHI0 = 500.0, MUB = 200.0;            // core steel initial susceptibility; ball relative permeability (linear)
      constexpr double GAP_MIN = 0.1e-3;                     // contact gap convention: the tables start at 0.1 mm (surface finish / paint / mylar), as every reference run
      constexpr double RIM_IN = 5.5e-3, RIM_OUT = 2.5e-3;    // rim rule trusted points: this far inside / outside the face edge (user, fixed for every core)
      // geometry + force grid: offset x 0.5 mm steps to 30 mm, 1 mm steps beyond; gap h log-spaced (49 steps 0.1 -> 30 mm), the same ratio continued. Both run out to the
      // furthest bare table, so a ball anywhere in the magnet's reach is a lookup (2026-10-20: the direct order-4 evaluation beyond 30 mm cost 8.4 us per call and stalled 6 balls)
      constexpr double X_KNEE = 30e-3, DX_NEAR = 0.5e-3, DX_FAR = 1e-3;
      constexpr int NX_NEAR = 60, NH_NEAR = 50;
      inline double HStepLog() { static const double v = std::log(30e-3 / GAP_MIN) / (NH_NEAR - 1); return v; }
      // (a 5 % tolerance for the slightly rising unsaturated rim was tried 2026-10-18 and removed by the user: NI 2000 at 8-10 mm reads 0.57-0.59 of the FEM, accepted)
      double Now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

      void EllipKE(double m, double &K, double &E)   // AGM, m = k^2
      {
         double a = 1, b = std::sqrt(1 - m), c2sum = 0.5 * m, p = 1;
         for (int i = 0; i < 30; i++)
         {
            const double an = 0.5 * (a + b), bn = std::sqrt(a * b), cn = 0.5 * (a - b);
            p *= 2; c2sum += 0.5 * p * cn * cn; a = an; b = bn;
            if (std::fabs(cn) < 1e-16 * a) break;
         }
         K = PI / (2 * a); E = K * (1 - c2sum);
      }
      inline void RingH(double Q, double a, double zr, double rho, double z, double &Hr, double &Hz)   // ring of magnetic charge Q
      {
         const double d = z - zr, al2 = (rho + a) * (rho + a) + d * d, be2 = (rho - a) * (rho - a) + d * d;
         double K, E; EllipKE(std::min(4 * a * rho / al2, 1 - 1e-15), K, E); const double al = std::sqrt(al2), c = Q / (4 * PI) * 2 / PI;
         Hz += c * d * E / (al * be2);
         if (rho > 1e-12) Hr += c / (2 * rho * al) * (K - (a * a - rho * rho + d * d) * E / be2);
      }
      inline void LoopH(double I, double a, double zr, double rho, double z, double &Hr, double &Hz)   // current loop, H = B / mu0
      {
         const double d = z - zr, al2 = (rho + a) * (rho + a) + d * d, be2 = (rho - a) * (rho - a) + d * d;
         double K, E; EllipKE(std::min(4 * a * rho / al2, 1 - 1e-15), K, E); const double al = std::sqrt(al2), c = I / (2 * PI);
         Hz += c / al * (K + (a * a - rho * rho - d * d) * E / be2);
         if (rho > 1e-12) Hr += c * d / (rho * al) * (-K + (a * a + rho * rho + d * d) * E / be2);
      }
      void LegGauss(int n, std::vector<double> &x, std::vector<double> &w)
      {
         x.resize(n); w.resize(n);
         for (int i = 0; i < n; i++)
         {
            double t = std::cos(PI * (i + 0.75) / (n + 0.5)), dp = 0;
            for (int it = 0; it < 100; it++)
            {
               double p0 = 1, p1 = t;
               for (int k = 2; k <= n; k++) { const double p2 = ((2 * k - 1) * t * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
               dp = n * (t * p1 - p0) / (t * t - 1); const double dt = p1 / dp; t -= dt; if (std::fabs(dt) < 1e-16) break;
            }
            x[n - 1 - i] = t; w[n - 1 - i] = 2 / ((1 - t * t) * dp * dp);
         }
      }
      bool LU(std::vector<double> &A, std::vector<int> &piv, int n)
      {
         for (int k = 0; k < n; k++)
         {
            int p = k; for (int i = k + 1; i < n; i++) if (std::fabs(A[i * n + k]) > std::fabs(A[p * n + k])) p = i;
            piv[k] = p; if (p != k) for (int j = 0; j < n; j++) std::swap(A[k * n + j], A[p * n + j]);
            const double d = A[k * n + k]; if (d == 0) return false;
            for (int i = k + 1; i < n; i++) { const double f = (A[i * n + k] /= d); if (f != 0) for (int j = k + 1; j < n; j++) A[i * n + j] -= f * A[k * n + j]; }
         }
         return true;
      }
      void LUSolve(const std::vector<double> &A, const std::vector<int> &piv, int n, std::vector<double> &b)
      {
         for (int k = 0; k < n; k++) std::swap(b[k], b[piv[k]]);
         for (int k = 0; k < n; k++) for (int i = k + 1; i < n; i++) b[i] -= A[i * n + k] * b[k];
         for (int k = n - 1; k >= 0; k--) { for (int j = k + 1; j < n; j++) b[k] -= A[k * n + j] * b[j]; b[k] /= A[k * n + k]; }
      }
      void Legendre(int L, double x, std::vector<double> &P, std::vector<double> &dP)
      {
         P.assign(L + 1, 0); dP.assign(L + 1, 0); P[0] = 1; if (L >= 1) { P[1] = x; dP[1] = 1; }
         for (int n = 1; n < L; n++) { P[n + 1] = ((2 * n + 1) * x * P[n] - n * P[n - 1]) / (n + 1); dP[n + 1] = dP[n - 1] + (2 * n + 1) * P[n]; }
      }

      // Gauss rules on [0, 1], per thread (the builders run on several threads)
      const std::pair<std::vector<double>, std::vector<double>> &Gauss01(int n)
      {
         thread_local std::map<int, std::pair<std::vector<double>, std::vector<double>>> cache;
         auto it = cache.find(n); if (it != cache.end()) return it->second;
         std::vector<double> x, w; LegGauss(n, x, w); for (int i = 0; i < n; i++) { x[i] = 0.5 * (x[i] + 1); w[i] *= 0.5; }
         return cache[n] = { x, w };
      }

      // ---- core cells (bench_bare3): elements of the per-cell magnetisation, dense rings near a field point, axis law in the centre band
      struct Elem { int kind; double r0, r1, z0, z1; };      // 0 axial face, 1 radial face, 2 volume, 3 volume with the axis law
      constexpr double KD = 4.0; constexpr int VCAP = 48;
      inline int NRule(const Elem &e, double rho, double z)
      {
         const double L = std::max(e.r1 - e.r0, e.z1 - e.z0);
         const double d = std::max(std::hypot(std::max(std::max(e.r0 - rho, 0.0), rho - e.r1), std::max(std::max(e.z0 - z, 0.0), z - e.z1)), 1e-7);
         const int nd = (int)std::ceil(KD * L / d); if (nd > 6) return std::min(nd, 600);
         const double dc = std::hypot(rho - 0.5 * (e.r0 + e.r1), z - 0.5 * (e.z0 + e.z1));
         return dc < 3 * L ? 6 : (dc < 8 * L ? 2 : 1);
      }
      void ElemH(const Elem &e, double s, double rho, double z, double &Hr, double &Hz)
      {
         int n = NRule(e, rho, z); if (e.kind >= 2) n = std::min(n, VCAP);
         const auto &g = Gauss01(n); const std::vector<double> &x = g.first, &w = g.second;
         if (e.kind == 0) for (int i = 0; i < n; i++) { const double r = e.r0 + x[i] * (e.r1 - e.r0); RingH(s * 2 * PI * r * (e.r1 - e.r0) * w[i], r, e.z0, rho, z, Hr, Hz); }
         else if (e.kind == 1) { if (e.r0 < 1e-12) return; for (int i = 0; i < n; i++) RingH(s * 2 * PI * e.r0 * (e.z1 - e.z0) * w[i], e.r0, e.z0 + x[i] * (e.z1 - e.z0), rho, z, Hr, Hz); }
         else
         {
            const double rc = 0.5 * (e.r0 + e.r1);
            for (int i = 0; i < n; i++) for (int k = 0; k < n; k++)
            {
               const double r = e.r0 + x[k] * (e.r1 - e.r0); double Q = s * 2 * PI * (e.r1 - e.r0) * w[k] * (e.z1 - e.z0) * w[i];
               if (e.kind == 3) Q *= 2 * r / rc;
               RingH(Q, r, e.z0 + x[i] * (e.z1 - e.z0), rho, z, Hr, Hz);
            }
         }
      }
      struct Block { double r0, r1, z0, z1, I; int nr, nz; };   // coil block per unit ampere-turn (I = share of 1 At)
      inline void BlockH(const Block &b, double NI, double rho, double z, double &Hr, double &Hz)
      {
         const double rc = 0.5 * (b.r0 + b.r1), zc = 0.5 * (b.z0 + b.z1), L = std::max(b.r1 - b.r0, b.z1 - b.z0), d = std::hypot(rho - rc, z - zc);
         const int mr = d < 3 * L ? b.nr : (d < 8 * L ? 2 : 1), mz = d < 3 * L ? b.nz : (d < 8 * L ? 2 : 1); const double I = NI * b.I * b.nr * b.nz / (mr * mz);
         for (int i = 0; i < mr; i++) for (int k = 0; k < mz; k++) LoopH(I, b.r0 + (i + 0.5) * (b.r1 - b.r0) / mr, b.z0 + (k + 0.5) * (b.z1 - b.z0) / mz, rho, z, Hr, Hz);
      }

      struct Cells
      {
         int nr = 4, nz = 30, N = 0, n2 = 0; double a = 0, zb = 0, Ms = 0, aa = 0, axs = 2;
         std::vector<double> re, ze, rc, zc, Nm, Hc1;          // Hc1: coil field at the cells per ampere-turn
         std::vector<Block> coil;
         void Setup(const Geometry &g)
         {
            a = g.coreRadius; zb = -(g.exposed + g.windLen); Ms = g.bSat / MU0; aa = Ms / CHI0;
            re.resize(nr + 1); ze.resize(nz + 1);
            for (int i = 0; i <= nr; i++) re[i] = a * std::sin(0.5 * PI * i / nr);
            for (int i = 0; i <= nz; i++) ze[i] = zb + (0 - zb) * 0.5 * (1 - std::cos(PI * i / nz));
            N = nr * nz; n2 = 2 * N; rc.resize(N); zc.resize(N);
            for (int i = 0; i < nr; i++) for (int k = 0; k < nz; k++) { rc[i * nz + k] = 0.5 * (re[i] + re[i + 1]); zc[i * nz + k] = 0.5 * (ze[k] + ze[k + 1]); }
            const int cr = 8, cz = 60, br = 4, bz = 10; const double c0 = -g.exposed - g.windLen, c1 = -g.exposed, I = 1.0 / (cr * cz);
            coil.clear();
            for (int i = 0; i < cr / br; i++) for (int k = 0; k < cz / bz; k++)
               coil.push_back({ g.windIn + i * br * (g.windOut - g.windIn) / cr, g.windIn + (i + 1) * br * (g.windOut - g.windIn) / cr, c0 + k * bz * (c1 - c0) / cz, c0 + (k + 1) * bz * (c1 - c0) / cz, I, br, bz });
            axs = re[1] / (0.5 * (re[0] + re[1]));
            Nm.assign((size_t)n2 * n2, 0.0); Hc1.assign(n2, 0.0);
            for (int j = 0; j < N; j++) for (int i = 0; i < N; i++)
            {
               double hr = 0, hz = 0; CellH(j, false, rc[i], zc[i], hr, hz); Nm[(size_t)i * n2 + j] = hr; Nm[(size_t)(N + i) * n2 + j] = hz;
               hr = hz = 0; CellH(j, true, rc[i], zc[i], hr, hz); Nm[(size_t)i * n2 + N + j] = hr; Nm[(size_t)(N + i) * n2 + N + j] = hz;
            }
            for (int i = 0; i < N; i++) { double hr = 0, hz = 0; for (auto &b : coil) BlockH(b, 1.0, rc[i], zc[i], hr, hz); Hc1[i] = hr; Hc1[N + i] = hz; }
         }
         void CellH(int j, bool zcomp, double rho, double z, double &hr, double &hz) const   // unit M_r (or M_z) of cell j
         {
            const int i0 = j / nz, k0 = j % nz; const double r0 = re[i0], r1 = re[i0 + 1], z0 = ze[k0], z1 = ze[k0 + 1];
            if (!zcomp) { const Elem out { 1, r1, r1, z0, z1 }, in { 1, r0, r0, z0, z1 }, vol { i0 == 0 ? 3 : 2, r0, r1, z0, z1 };
               ElemH(out, i0 == 0 ? axs : 1.0, rho, z, hr, hz); if (i0 > 0) ElemH(in, -1, rho, z, hr, hz); ElemH(vol, -1, rho, z, hr, hz); }
            else { const Elem tp { 0, r0, r1, z1, z1 }, bt { 0, r0, r1, z0, z0 }; ElemH(tp, 1, rho, z, hr, hz); ElemH(bt, -1, rho, z, hr, hz); }
         }
      };

      // Newton on the cells with a frozen Jacobian (re-factored every CHORD steps or after a halved step); one solver per thread
      struct Newton
      {
         const Cells *c = nullptr; int CHORD = 50, since = 1 << 30, nfact = 0; bool halved = false;
         std::vector<double> H, Rv, J, dM, Mt, Rn, Hn; std::vector<int> piv;
         void Init(const Cells &cells) { c = &cells; const int n2 = cells.n2; H.resize(n2); Rv.resize(n2); J.resize((size_t)n2 * n2); dM.resize(n2); Mt.resize(n2); Rn.resize(n2); Hn.resize(n2); piv.resize(n2); since = 1 << 30; }
         double Resid(const std::vector<double> &Hs, const std::vector<double> &Mx, std::vector<double> &R_, std::vector<double> &Hx)
         {
            const int n2 = c->n2, N = c->N; double q = 0;
            for (int i = 0; i < n2; i++) { double s = Hs[i]; const double *row = &c->Nm[(size_t)i * n2]; for (int j = 0; j < n2; j++) s += row[j] * Mx[j]; Hx[i] = s; }
            for (int i = 0; i < N; i++) { const double h = std::max(std::hypot(Hx[i], Hx[N + i]), 1e-300), f = c->Ms * h / (h + c->aa);
               R_[i] = Mx[i] - f * Hx[i] / h; R_[N + i] = Mx[N + i] - f * Hx[N + i] / h; q += R_[i] * R_[i] + R_[N + i] * R_[N + i]; }
            return std::sqrt(q);
         }
         int Solve(const std::vector<double> &Hs, std::vector<double> &M, int maxit)
         {
            const int n2 = c->n2, N = c->N; const double Ms = c->Ms, aa = c->aa;
            double nR = Resid(Hs, M, Rv, H); int it;
            for (it = 1; it <= maxit; it++)
            {
               if (since >= CHORD || halved)
               {
                  since = 0; nfact++;
                  for (int i = 0; i < N; i++)
                  {
                     const double h = std::max(std::hypot(H[i], H[N + i]), 1e-300), er = H[i] / h, ez = H[N + i] / h, f = Ms * h / (h + aa), s = f / h, fp = Ms * aa / ((h + aa) * (h + aa)), d = fp - s;
                     const double Drr = s + d * er * er, Drz = d * er * ez, Dzz = s + d * ez * ez;
                     const double *a_ = &c->Nm[(size_t)i * n2], *b_ = &c->Nm[(size_t)(N + i) * n2]; double *jr = &J[(size_t)i * n2], *jz = &J[(size_t)(N + i) * n2];
                     for (int j = 0; j < n2; j++) { jr[j] = -(Drr * a_[j] + Drz * b_[j]); jz[j] = -(Drz * a_[j] + Dzz * b_[j]); }
                     jr[i] += 1; jz[N + i] += 1;
                  }
                  LU(J, piv, n2);
               }
               since++; for (int i = 0; i < n2; i++) dM[i] = -Rv[i]; LUSolve(J, piv, n2, dM);
               double lam = 1, nRn;
               while (true) { for (int i = 0; i < n2; i++) Mt[i] = M[i] + lam * dM[i]; nRn = Resid(Hs, Mt, Rn, Hn); if (nRn < nR || lam < 1e-3) break; lam *= 0.5; }
               double mx = 0; for (int i = 0; i < n2; i++) mx = std::max(mx, std::fabs(lam * dM[i]));
               halved = lam < 1; M = Mt; Rv = Rn; H = Hn; nR = nRn; if (mx / Ms < 1e-9) break;
            }
            return std::min(it, maxit);
         }
      };

      // ---- smooth read-out of the coil-only state into the bare table (bench_bare3 'smooth')
      std::vector<double> AxisGrid(double lo, double hi, const std::vector<double> &fine, double hf, double hc, double gr)
      {
         std::set<double> p = { lo, hi };
         for (double f : fine) for (int sg = -1; sg <= 1; sg += 2) { double x = f, h = hf; while (x >= lo && x <= hi) { p.insert(x); x += sg * h; h = std::min(h * gr, x > 0.03 ? 4 * hc : hc); } }
         std::vector<double> g; for (double x : p) if (g.empty() || x - g.back() > 0.25 * hf) g.push_back(x); return g;
      }
      // reach < 0: first find this state's reach (largest distance from the face centre, rays 0..90 deg from the axis, where |B| >= B_GATE) and return it in
      // reach; the table then covers reach + one ball diameter in r and z (a ball's surface, and the drive projection around its centre, stay inside)
      void BuildBare(const Cells &C, const Geometry &geo, double NI, const std::vector<double> &M, BareTable &T, double &reach)
      {
         const int nr = C.nr, nz = C.nz, N = C.N; const double a = C.a, zb = C.zb;
         std::vector<double> rn = { 0.0 }, zn = { zb };
         for (int i = 0; i < nr; i++) rn.push_back(0.5 * (C.re[i] + C.re[i + 1])); rn.push_back(a);
         for (int k = 0; k < nz; k++) zn.push_back(0.5 * (C.ze[k] + C.ze[k + 1])); zn.push_back(0.0);
         const int NR = (int)rn.size(), NZ = (int)zn.size(); std::vector<double> nMr(NR * NZ), nMz(NR * NZ);
         for (int j = 0; j < NR; j++) for (int l = 0; l < NZ; l++) { const int i = std::min(std::max(j - 1, 0), nr - 1), k = std::min(std::max(l - 1, 0), nz - 1);
            nMr[j * NZ + l] = j == 0 ? 0.0 : M[i * nz + k]; nMz[j * NZ + l] = M[N + i * nz + k]; }
         auto Mat = [&](double r, double z, double &mr, double &mz, double &dmr_dr, double &dmz_dz)
         {
            int j = 0; while (j < NR - 2 && r > rn[j + 1]) j++; int l = 0; while (l < NZ - 2 && z > zn[l + 1]) l++;
            const double t = (r - rn[j]) / (rn[j + 1] - rn[j]), u = (z - zn[l]) / (zn[l + 1] - zn[l]);
            auto bl = [&](const std::vector<double> &v, double &val, double &dr, double &dz) {
               const double v00 = v[j * NZ + l], v10 = v[(j + 1) * NZ + l], v01 = v[j * NZ + l + 1], v11 = v[(j + 1) * NZ + l + 1];
               val = (1 - t) * (1 - u) * v00 + t * (1 - u) * v10 + (1 - t) * u * v01 + t * u * v11;
               dr = ((1 - u) * (v10 - v00) + u * (v11 - v01)) / (rn[j + 1] - rn[j]); dz = ((1 - t) * (v01 - v00) + t * (v11 - v10)) / (zn[l + 1] - zn[l]); };
            double d1, d2; bl(nMr, mr, dmr_dr, d1); bl(nMz, mz, d2, dmz_dz);
         };
         struct SE { int kind; double r0, r1, z0, z1; };
         std::vector<SE> se;
         for (int j = 0; j < NR - 1; j++) { se.push_back({ 0, rn[j], rn[j + 1], 0.0, 0.0 }); se.push_back({ 1, rn[j], rn[j + 1], zb, zb }); }
         for (int l = 0; l < NZ - 1; l++) se.push_back({ 2, a, a, zn[l], zn[l + 1] });
         for (int j = 0; j < NR - 1; j++) for (int l = 0; l < NZ - 1; l++) se.push_back({ 3, rn[j], rn[j + 1], zn[l], zn[l + 1] });
         auto seH = [&](const SE &e, double rho, double z, double &hr, double &hz)
         {
            const Elem g { e.kind == 3 ? 2 : (e.kind == 2 ? 1 : 0), e.r0, e.r1, e.z0, e.z1 };
            int n = NRule(g, rho, z); if (e.kind == 3) n = std::min(n, VCAP);
            const auto &G = Gauss01(n); const std::vector<double> &x = G.first, &w = G.second; double mr, mz, dr, dz;
            if (e.kind <= 1) for (int i = 0; i < n; i++) { const double r = e.r0 + x[i] * (e.r1 - e.r0); Mat(r, e.z0, mr, mz, dr, dz);
               RingH((e.kind == 0 ? mz : -mz) * 2 * PI * r * (e.r1 - e.r0) * w[i], r, e.z0, rho, z, hr, hz); }
            else if (e.kind == 2) for (int i = 0; i < n; i++) { const double zz = e.z0 + x[i] * (e.z1 - e.z0); Mat(a, zz, mr, mz, dr, dz);
               RingH(mr * 2 * PI * a * (e.z1 - e.z0) * w[i], a, zz, rho, z, hr, hz); }
            else for (int i = 0; i < n; i++) for (int k = 0; k < n; k++) { const double r = e.r0 + x[k] * (e.r1 - e.r0), zz = e.z0 + x[i] * (e.z1 - e.z0); Mat(r, zz, mr, mz, dr, dz);
               const double rho_v = -((r > 1e-12 ? mr / r : 0) + dr + dz); RingH(rho_v * 2 * PI * r * (e.r1 - e.r0) * w[k] * (e.z1 - e.z0) * w[i], r, zz, rho, z, hr, hz); }
         };
         auto Bmag = [&](double rho, double z) { double hr = 0, hz = 0; for (auto &b : C.coil) BlockH(b, NI, rho, z, hr, hz); for (auto &e : se) seH(e, rho, z, hr, hz); return MU0 * std::hypot(hr, hz); };
         if (reach < 0)
         {
            reach = 0;
            for (int k = 0; k <= 9; k++)   // march out in 1 mm steps; the reach is the last point above the gate, with 20 mm below it confirming (the field near the core is not monotone)
            {
               const double t = 0.5 * PI * k / 9, ur = std::sin(t), uz = std::max(std::cos(t), 1e-3);
               double last = 0; int below = 0;
               for (double s = 5e-3; s < 1.0 && below < 20; s += 1e-3) { if (Bmag(s * ur, std::max(s * uz, 0.05e-3)) >= B_GATE) { last = s; below = 0; } else below++; }
               reach = std::max(reach, last);
            }
         }
         const double ext = reach + 2 * geo.ballR;
         T.r = AxisGrid(0.0, ext, { 0.0, a }, 0.25e-3, 2e-3, 1.2);
         T.z = AxisGrid(0.05e-3, ext, { 0.05e-3 }, 0.1e-3, 2e-3, 1.2);
         const size_t np = T.r.size() * T.z.size();
         for (int c = 0; c < 2; c++) T.F[c].assign(np, 0.0);
         for (size_t i = 0; i < T.r.size(); i++) for (size_t k = 0; k < T.z.size(); k++)
         {
            double hr = 0, hz = 0; for (auto &b : C.coil) BlockH(b, NI, T.r[i], T.z[k], hr, hz);
            for (auto &e : se) seH(e, T.r[i], T.z[k], hr, hz);
            T.F[0][i * T.z.size() + k] = hr; T.F[1][i * T.z.size() + k] = hz;
         }
         T.Build();
      }

      // ---- order-4 linear ball on an inner sphere (bench_force, cheap rule (L+2) x (2L+4) at 0.5 R)
      struct BallRule
      {
         int L = 4, nth = 0, nph = 0; double R = 0, Rs = 0; std::vector<double> nx, ny, nz_, wt, P, cosm;
         void Init(double ballR)
         {
            R = ballR; Rs = 0.5 * ballR; nth = L + 2; nph = 2 * L + 4; std::vector<double> mu, w; LegGauss(nth, mu, w);
            P.assign((size_t)nth * (L + 1) * (L + 1), 0.0); cosm.resize((size_t)(L + 1) * nph); nx.clear(); ny.clear(); nz_.clear(); wt.clear();
            for (int j = 0; j < nph; j++) for (int m = 0; m <= L; m++) cosm[(size_t)m * nph + j] = std::cos(m * 2 * PI * j / nph);
            for (int i = 0; i < nth; i++)
            {
               const double x = mu[i], sx = std::sqrt(std::max(0.0, 1 - x * x)); double pmm = std::sqrt(1 / (4 * PI));
               for (int m = 0; m <= L; m++)
               {
                  if (m > 0) pmm *= std::sqrt((2.0 * m + 1) / (2.0 * m)) * sx;
                  double pm2 = 0, pm1 = pmm; auto put = [&](int n, double v) { if (n <= L) P[((size_t)i * (L + 1) + m) * (L + 1) + n] = v * (m > 0 ? std::sqrt(2.0) : 1.0); };
                  put(m, pmm);
                  if (m + 1 <= L) { double p = std::sqrt(2.0 * m + 3) * x * pmm; put(m + 1, p); pm2 = pmm; pm1 = p;
                     for (int n = m + 2; n <= L; n++) { const double aa = std::sqrt((4.0 * n * n - 1) / (1.0 * n * n - 1.0 * m * m)), b = std::sqrt(((n - 1.0) * (n - 1) - 1.0 * m * m) / (4.0 * (n - 1) * (n - 1) - 1));
                        const double p2 = aa * (x * pm1 - b * pm2); put(n, p2); pm2 = pm1; pm1 = p2; } }
               }
               for (int j = 0; j < nph; j++) { const double ph = 2 * PI * j / nph; nx.push_back(sx * std::cos(ph)); ny.push_back(sx * std::sin(ph)); nz_.push_back(x); wt.push_back(w[i] * 2 * PI / nph); }
            }
         }
         // force on the ball centred at offset xc, height zc above the face: Fx (+ = away from the axis), Fz (+ = away from the face)
         void Force(const BareTable &T, double xc, double zc, double &Fx, double &Fz) const
         {
            const int np = nth * nph; double g0[128], gx[128], gz[128];
            for (int p = 0; p < np; p++)
            {
               const double X = xc + Rs * nx[p], Y = Rs * ny[p], Z = zc + Rs * nz_[p], rho = std::hypot(X, Y); double v[2], vr[2], vz[2]; T.Eval(rho, Z, v, vr, vz);
               const double cx = rho > 1e-12 ? X / rho : 0, cy = rho > 1e-12 ? Y / rho : 0, inv = rho > 1e-12 ? 1 / rho : 0, rad = cx * nx[p] + cy * ny[p];
               g0[p] = v[0] * rad + v[1] * nz_[p];
               gx[p] = vr[0] * cx * rad + v[0] * (cy * cy * inv * nx[p] - cx * cy * inv * ny[p]) + vr[1] * cx * nz_[p];
               gz[p] = vz[0] * rad + vz[1] * nz_[p];
            }
            double sx = 0, sz = 0, f0[64], fx[64], fz[64];
            for (int i = 0; i < nth; i++) for (int m = 0; m <= L; m++) { const double *cm = &cosm[(size_t)m * nph]; const size_t o = (size_t)i * nph; double a0 = 0, ax = 0, az = 0;
               for (int j = 0; j < nph; j++) { a0 += g0[o + j] * cm[j]; ax += gx[o + j] * cm[j]; az += gz[o + j] * cm[j]; }
               const double w = wt[o]; f0[i * (L + 1) + m] = a0 * w; fx[i * (L + 1) + m] = ax * w; fz[i * (L + 1) + m] = az * w; }
            for (int m = 0; m <= L; m++) for (int n = std::max(m, 1); n <= L; n++)
            {
               double p0 = 0, px = 0, pz = 0;
               for (int i = 0; i < nth; i++) { const double y = P[((size_t)i * (L + 1) + m) * (L + 1) + n]; const int q = i * (L + 1) + m; p0 += f0[q] * y; px += fx[q] * y; pz += fz[q] * y; }
               const double sc = std::pow(R / Rs, n - 1); p0 *= sc; px *= sc; pz *= sc;
               const double wn = (2 * n + 1) * (n * (MUB - 1) / (n * MUB + n + 1)) * R * R * R / (1.0 * n * n);
               sx += wn * p0 * px; sz += wn * p0 * pz;
            }
            Fx = MU0 * sx; Fz = MU0 * sz;
         }
      };

      // ---- ball saturation (face_neck.txt): F_sat / F_lin of a ball on the face vs the linear gap field Bgap at the gap midpoint, per gap; 3 cores x 3 currents
      struct SatRow { double gap; std::vector<std::pair<double, double>> pts; };
      const std::vector<SatRow> &SatData()
      {
         static const std::vector<SatRow> d = {
            { 0.05e-3, { { 2.01, 0.890 }, { 2.06, 0.891 }, { 2.09, 0.884 }, { 2.39, 0.875 }, { 2.50, 0.868 }, { 2.54, 0.845 }, { 2.58, 0.842 }, { 2.71, 0.828 }, { 2.72, 0.826 } } },
            { 0.10e-3, { { 1.91, 0.913 }, { 1.97, 0.912 }, { 2.02, 0.903 }, { 2.35, 0.886 }, { 2.47, 0.879 }, { 2.52, 0.857 }, { 2.56, 0.854 }, { 2.69, 0.842 }, { 2.71, 0.841 } } },
            { 0.25e-3, { { 1.26, 0.971 }, { 1.39, 0.967 }, { 1.46, 0.961 }, { 2.09, 0.918 }, { 2.24, 0.906 }, { 2.30, 0.888 }, { 2.34, 0.884 }, { 2.50, 0.867 }, { 2.52, 0.865 } } },
            { 0.50e-3, { { 0.77, 0.989 }, { 0.86, 0.988 }, { 0.90, 0.986 }, { 1.66, 0.961 }, { 1.88, 0.945 }, { 1.92, 0.932 }, { 1.98, 0.926 }, { 2.18, 0.905 }, { 2.21, 0.902 } } },
            { 1.00e-3, { { 0.46, 0.995 }, { 0.52, 0.995 }, { 0.54, 0.994 }, { 1.08, 0.987 }, { 1.28, 0.977 }, { 1.30, 0.981 }, { 1.36, 0.975 }, { 1.64, 0.959 }, { 1.67, 0.957 } } } };
         return d;
      }
      double SatFactor(double gap, double Bgap)
      {
         const auto &d = SatData();
         auto row = [&](const SatRow &r) {   // least-squares line in Bgap, clamped to [0.5, 1]
            double n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0; for (auto &p : r.pts) { n++; sx += p.first; sy += p.second; sxx += p.first * p.first; sxy += p.first * p.second; }
            const double b = (n * sxy - sx * sy) / (n * sxx - sx * sx), a0 = (sy - b * sx) / n; return std::min(1.0, std::max(0.5, a0 + b * Bgap)); };
         if (gap >= 3e-3) return 1.0;
         if (gap <= d.front().gap) return row(d.front());
         for (size_t i = 0; i + 1 < d.size(); i++)
            if (gap <= d[i + 1].gap) { const double u = std::log(gap / d[i].gap) / std::log(d[i + 1].gap / d[i].gap); return (1 - u) * row(d[i]) + u * row(d[i + 1]); }
         const double u = std::log(gap / d.back().gap) / std::log(3e-3 / d.back().gap); return (1 - u) * row(d.back()) + u;   // 1 mm -> 3 mm: toward 1 (FEM >= 0.987 from 2 mm)
      }

      // ---- on-axis cell + ball solve at the 6 reference gaps (bench_axial, chord Newton, 2 steps per outer iteration, warm start far -> near)
      void OnAxis(const Cells &C, const Geometry &geo, double NI, const std::vector<double> &M0, Newton &nw, std::vector<double> &hR, std::vector<double> &Rlin, std::vector<double> &sat,
                  std::vector<std::vector<double>> &plugB, std::vector<double> &Mplug, std::ostringstream &lg)
      {
         const int L = 12, N = C.N, n2 = C.n2; const double R = geo.ballR;
         const double GAPS[] = { 10.1e-3, 3.1e-3, 1.1e-3, 0.6e-3, 0.3e-3, 0.1e-3 };
         std::vector<double> beta(L + 2, 0.0); for (int n = 1; n <= L; n++) beta[n] = n * (MUB - 1) / (n * MUB + n + 1);
         const int nq = std::max(4 * L, 80); std::vector<double> mu, wq; LegGauss(nq, mu, wq);
         std::vector<double> Hc(n2); for (int i = 0; i < n2; i++) Hc[i] = NI * C.Hc1[i];
         std::vector<double> Mw = M0, abw(L + 2, 0.0); bool have = false;
         hR.clear(); Rlin.clear(); sat.clear(); plugB.clear();
         for (double g : GAPS)
         {
            const double zcen = R + g, rho = R + 0.5 * g;
            std::vector<double> Hcoil(2 * nq), E((size_t)2 * nq * n2, 0.0);
            for (int q = 0; q < nq; q++)
            {
               const double st = std::sqrt(std::max(0.0, 1 - mu[q] * mu[q])), pr = rho * st, pz = zcen + rho * mu[q];
               double hr = 0, hz = 0; for (auto &b : C.coil) BlockH(b, NI, pr, pz, hr, hz); Hcoil[q] = hr; Hcoil[nq + q] = hz;
               for (int j = 0; j < N; j++)
               {
                  hr = hz = 0; C.CellH(j, false, pr, pz, hr, hz); E[(size_t)q * n2 + j] = hr; E[(size_t)(nq + q) * n2 + j] = hz;
                  hr = hz = 0; C.CellH(j, true, pr, pz, hr, hz); E[(size_t)q * n2 + N + j] = hr; E[(size_t)(nq + q) * n2 + N + j] = hz;
               }
            }
            std::vector<double> Bc((size_t)L * n2), P, dP;
            for (int i = 0; i < N; i++)
            {
               const double dz = C.zc[i] - zcen, r = std::hypot(C.rc[i], dz), ct = dz / r, st = C.rc[i] / r; Legendre(L, ct, P, dP);
               for (int n = 1; n <= L; n++) { const double f = std::pow(R / r, n + 2) / R, Hr_ = (n + 1) * f * P[n], Ht = f * st * dP[n];
                  Bc[(size_t)(n - 1) * n2 + i] = Hr_ * st + Ht * ct; Bc[(size_t)(n - 1) * n2 + N + i] = Hr_ * ct - Ht * st; }
            }
            std::vector<std::vector<double>> PL(nq); for (int q = 0; q < nq; q++) { std::vector<double> d_; Legendre(L + 1, mu[q], PL[q], d_); }
            auto coeffs = [&](const std::vector<double> &M, std::vector<double> &b)
            {
               b.assign(L + 2, 0.0);
               for (int q = 0; q < nq; q++)
               {
                  double hr = Hcoil[q], hz = Hcoil[nq + q]; const double *er = &E[(size_t)q * n2], *ez = &E[(size_t)(nq + q) * n2];
                  for (int j = 0; j < n2; j++) { hr += er[j] * M[j]; hz += ez[j] * M[j]; }
                  const double st = std::sqrt(std::max(0.0, 1 - mu[q] * mu[q])), Hrad = hr * st + hz * mu[q];
                  for (int n = 1; n <= L + 1; n++) b[n] += wq[q] * Hrad * PL[q][n];
               }
               for (int n = 1; n <= L + 1; n++) b[n] = -(2 * n + 1) / (2.0 * n * std::pow(rho, n - 1)) * b[n] * std::pow(R, n);
            };
            // field of the truncated expansions at direction ct on the sphere (H radial, H theta), and the Maxwell force (+ toward the core)
            auto fieldAt = [&](const std::vector<double> &b, const std::vector<double> &ab, double ct, double &Hr_, double &Ht)
            {
               const double st = std::sqrt(std::max(0.0, 1 - ct * ct)); std::vector<double> P2, dP2; Legendre(L + 1, ct, P2, dP2); Hr_ = Ht = 0;
               for (int n = 1; n <= L + 1; n++) { const double bn = b[n] / std::pow(R, n), rn1 = std::pow(rho, n - 1); Hr_ += -n * bn * rn1 * P2[n]; Ht += bn * rn1 * st * dP2[n]; }
               for (int n = 1; n <= L; n++) { const double f = ab[n] * std::pow(R / rho, n + 2) / R; Hr_ += (n + 1) * f * P2[n]; Ht += f * st * dP2[n]; }
            };
            auto force = [&](const std::vector<double> &b, const std::vector<double> &ab)
            {
               const int nf = 2 * L + 6; std::vector<double> x, w; LegGauss(nf, x, w); double F = 0;
               for (int q = 0; q < nf; q++) { const double ct = x[q], st = std::sqrt(std::max(0.0, 1 - ct * ct)); double Hr_, Ht; fieldAt(b, ab, ct, Hr_, Ht);
                  const double Bz = MU0 * (Hr_ * ct - Ht * st), Bn = MU0 * Hr_, B2 = MU0 * MU0 * (Hr_ * Hr_ + Ht * Ht);
                  F += (Bz * Bn / MU0 - B2 / (2 * MU0) * ct) * 2 * PI * rho * rho * w[q]; }
               return -F;
            };
            std::vector<double> b, ab(L + 2, 0.0), abn(L + 2, 0.0), M = M0, Hs(n2);
            coeffs(M0, b); for (int n = 1; n <= L; n++) ab[n] = -beta[n] * b[n]; const double F0 = force(b, ab);
            if (have) { M = Mw; ab = abw; }
            int outer;
            for (outer = 1; outer <= 300; outer++)
            {
               for (int i = 0; i < n2; i++) { double s = Hc[i]; for (int n = 1; n <= L; n++) s += ab[n] * Bc[(size_t)(n - 1) * n2 + i]; Hs[i] = s; }
               nw.Solve(Hs, M, 2); coeffs(M, b);
               double ch = 0, mx = 0; for (int n = 1; n <= L; n++) { abn[n] = -beta[n] * b[n]; ch = std::max(ch, std::fabs(abn[n] - ab[n])); mx = std::max(mx, std::fabs(abn[n])); }
               if (ch / std::max(mx, 1e-30) < 1e-7) { ab = abn; break; }
               for (int n = 1; n <= L; n++) ab[n] = 0.3 * ab[n] + 0.7 * abn[n];
            }
            for (int i = 0; i < n2; i++) { double s = Hc[i]; for (int n = 1; n <= L; n++) s += ab[n] * Bc[(size_t)(n - 1) * n2 + i]; Hs[i] = s; }
            nw.Solve(Hs, M, 50); coeffs(M, b); const double F2 = force(b, ab);
            double Hr_, Ht; fieldAt(b, ab, -1.0, Hr_, Ht);           // the gap midpoint (bottom of the sphere, on the axis): Bgap = mu0 |H|
            const double Bgap = MU0 * std::hypot(Hr_, Ht);
            Mw = M; abw = ab; have = true;
            hR.push_back(g); Rlin.push_back(F2 - F0); sat.push_back(SatFactor(g, Bgap));
            // drive of the plug for the ball-ball model: b[n] here is the Legendre-basis coefficient scaled by R^n; orthonormal Y_n0 = sqrt((2n+1)/4pi) P_n
            std::vector<double> bz(L + 2, 0.0); for (int n = 1; n <= L + 1; n++) bz[n] = b[n] / std::sqrt((2 * n + 1) / (4 * PI));
            plugB.push_back(bz); if (g == GAPS[5]) Mplug = M;
            lg << "    gap " << g * 1e3 << " mm: M0 " << F0 << " N, M2 " << F2 << " N, Bgap " << Bgap << " T, ball saturation " << sat.back() << (outer > 300 ? " (NOT CONVERGED)" : "") << "\n";
         }
         std::reverse(hR.begin(), hR.end()); std::reverse(Rlin.begin(), Rlin.end()); std::reverse(sat.begin(), sat.end()); std::reverse(plugB.begin(), plugB.end());
      }

      // log-log interpolation in the gap with a power-law tail beyond the last point (R falls as gap^-2.8 there); sat interpolated in log gap, 1 beyond
      double InterpR(const std::vector<double> &h, const std::vector<double> &R, double g)
      {
         g = std::max(g, h.front());
         if (g >= h.back()) { const double k = std::log(std::max(R[R.size() - 1], 1e-30) / std::max(R[R.size() - 2], 1e-30)) / std::log(h.back() / h[h.size() - 2]); return R.back() * std::pow(g / h.back(), k); }
         size_t i = 0; while (i + 2 < h.size() && g > h[i + 1]) i++;
         const double u = std::log(g / h[i]) / std::log(h[i + 1] / h[i]);
         if (R[i] > 0 && R[i + 1] > 0) return std::exp((1 - u) * std::log(R[i]) + u * std::log(R[i + 1]));
         return (1 - u) * R[i] + u * R[i + 1];
      }
      double InterpSat(const std::vector<double> &h, const std::vector<double> &s, double g)
      {
         g = std::max(g, h.front()); if (g >= h.back()) return s.back();
         size_t i = 0; while (i + 2 < h.size() && g > h[i + 1]) i++;
         const double u = std::log(g / h[i]) / std::log(h[i + 1] / h[i]); return (1 - u) * s[i] + u * s[i + 1];
      }
   }

   // ---- bare table: bicubic Hermite with spline slopes (bench_force Table)
   void BareTable::Build()
   {
      auto spline_d = [](const double *x, const double *y, int n, int stride, double *d, int dstride)
      {
         std::vector<double> u(n), y2(n); y2[0] = u[0] = 0;
         for (int i = 1; i < n - 1; i++) { const double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]), p = sig * y2[i - 1] + 2; y2[i] = (sig - 1) / p;
            const double dy = (y[(i + 1) * stride] - y[i * stride]) / (x[i + 1] - x[i]) - (y[i * stride] - y[(i - 1) * stride]) / (x[i] - x[i - 1]); u[i] = (6 * dy / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p; }
         y2[n - 1] = 0; for (int k = n - 2; k >= 0; k--) y2[k] = y2[k] * y2[k + 1] + u[k];
         for (int i = 0; i < n; i++) { const int a = i < n - 1 ? i : n - 2; const double h = x[a + 1] - x[a], s = (y[(a + 1) * stride] - y[a * stride]) / h;
            d[i * dstride] = i < n - 1 ? s - h * (2 * y2[a] + y2[a + 1]) / 6 : s + h * (y2[a] + 2 * y2[a + 1]) / 6; }
      };
      const int nr = (int)r.size(), nz = (int)z.size();
      for (int c = 0; c < 2; c++)
      {
         Fr[c].assign(nr * nz, 0); Fz[c].assign(nr * nz, 0); Frz[c].assign(nr * nz, 0);
         for (int i = 0; i < nr; i++) spline_d(z.data(), &F[c][i * nz], nz, 1, &Fz[c][i * nz], 1);
         for (int k = 0; k < nz; k++) spline_d(r.data(), &F[c][k], nr, nz, &Fr[c][k], nz);
         for (int i = 0; i < nr; i++) spline_d(z.data(), &Fr[c][i * nz], nz, 1, &Frz[c][i * nz], 1);
      }
   }
   void BareTable::Eval(double rho, double zz, double v[2], double vr[2], double vz[2]) const
   {
      const int nr = (int)r.size(), nz = (int)z.size();
      auto find = [](const std::vector<double> &g, double x) { int lo = 0, hi = (int)g.size() - 2; if (x <= g[0]) return 0; if (x >= g[hi + 1]) return hi;
         while (hi - lo > 1) { const int m = (lo + hi) / 2; if (g[m] <= x) lo = m; else hi = m; } return lo; };
      rho = std::min(std::max(rho, r[0]), r[nr - 1]); zz = std::min(std::max(zz, z[0]), z[nz - 1]);
      const int i = find(r, rho), k = find(z, zz); const double dr = r[i + 1] - r[i], dz = z[k + 1] - z[k], t = (rho - r[i]) / dr, u = (zz - z[k]) / dz;
      const double A[2] = { 2 * t * t * t - 3 * t * t + 1, -2 * t * t * t + 3 * t * t }, B[2] = { t * t * t - 2 * t * t + t, t * t * t - t * t };
      const double Ad[2] = { (6 * t * t - 6 * t) / dr, (-6 * t * t + 6 * t) / dr }, Bd[2] = { 3 * t * t - 4 * t + 1, 3 * t * t - 2 * t };
      const double Au[2] = { 2 * u * u * u - 3 * u * u + 1, -2 * u * u * u + 3 * u * u }, Bu[2] = { u * u * u - 2 * u * u + u, u * u * u - u * u };
      const double Aud[2] = { (6 * u * u - 6 * u) / dz, (-6 * u * u + 6 * u) / dz }, Bud[2] = { 3 * u * u - 4 * u + 1, 3 * u * u - 2 * u };
      for (int c = 0; c < 2; c++)
      {
         double s = 0, sr = 0, sz = 0;
         for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++)
         {
            const int idx = (i + a) * nz + k + b;
            const double f = F[c][idx], fr = Fr[c][idx] * dr, fz = Fz[c][idx] * dz, frz = Frz[c][idx] * dr * dz;
            s += f * A[a] * Au[b] + fr * B[a] * Au[b] + fz * A[a] * Bu[b] + frz * B[a] * Bu[b];
            sr += f * Ad[a] * Au[b] + fr * Bd[a] / dr * Au[b] + fz * Ad[a] * Bu[b] + frz * Bd[a] / dr * Bu[b];
            sz += f * A[a] * Aud[b] + fr * B[a] * Aud[b] + fz * A[a] * Bud[b] / dz + frz * B[a] * Bud[b] / dz;
         }
         v[c] = s; vr[c] = sr; vz[c] = sz;
      }
   }

   // The level's stored energy from its force table (Model::Energy). V = 0 on the far edges and beyond the level's bare table. Each grid edge asks V's difference along it to be the table's work along it
   // (exact for the table's interpolation); the residuals, in force units, are minimised in least squares by a band Cholesky over the grid (one row of the band per gap).
   static void BuildEnergy(const std::vector<double> &xg, const std::vector<double> &hg, const Geometry &g, Level &L)
   {
      const int NX = (int)xg.size(), NH = (int)hg.size(), n = NX * NH, W1 = NH + 1;
      std::vector<char> fixed(n);
      for (int i = 0; i < NX; i++) for (int k = 0; k < NH; k++) fixed[i * NH + k] = i == NX - 1 || k == NH - 1 || !L.bare.Inside(xg[i] + g.ballR, hg[k] + 2 * g.ballR);
      std::vector<double> A((size_t)n * W1, 0.0), b(n, 0.0);   // A[r * W1 + (r - c)]: the lower band
      auto edge = [&](int p, int q, double work, double w)      // V[q] - V[p] = work, weight w
      {
         if (!fixed[p]) { A[(size_t)p * W1] += w; b[p] -= w * work; }
         if (!fixed[q]) { A[(size_t)q * W1] += w; b[q] += w * work; }
         if (!fixed[p] && !fixed[q]) A[(size_t)std::max(p, q) * W1 + std::abs(p - q)] -= w;
      };
      for (int i = 0; i < NX; i++) for (int k = 0; k < NH; k++)
      {
         const int p = i * NH + k;
         if (i + 1 < NX) { const double dx = xg[i + 1] - xg[i]; edge(p, p + NH, 0.5 * (L.side[p] + L.side[p + NH]) * dx, 1.0 / (dx * dx)); }
         if (k + 1 < NH)
         {
            const double h0 = hg[k], h1 = hg[k + 1], dh = h1 - h0, la = std::log(h1 / h0);   // down is linear in ln h along the edge
            edge(p, p + 1, L.down[p] * dh + (L.down[p + 1] - L.down[p]) * (h1 - dh / la), 1.0 / (dh * dh));
         }
      }
      for (int p = 0; p < n; p++) if (fixed[p]) { A[(size_t)p * W1] = 1.0; b[p] = 0.0; }
      for (int r = 0; r < n; r++)
         for (int c = std::max(0, r - NH); c <= r; c++)
         {
            double s = A[(size_t)r * W1 + (r - c)];
            for (int m = std::max(0, r - NH); m < c; m++) if (c - m <= NH) s -= A[(size_t)r * W1 + (r - m)] * A[(size_t)c * W1 + (c - m)];
            A[(size_t)r * W1 + (r - c)] = r == c ? std::sqrt(s) : s / A[(size_t)c * W1];
         }
      for (int r = 0; r < n; r++) { double s = b[r]; for (int c = std::max(0, r - NH); c < r; c++) s -= A[(size_t)r * W1 + (r - c)] * b[c]; b[r] = s / A[(size_t)r * W1]; }
      for (int r = n - 1; r >= 0; r--) { double s = b[r]; for (int c = r + 1; c <= std::min(n - 1, r + NH); c++) s -= A[(size_t)c * W1 + (c - r)] * b[c]; b[r] = s / A[(size_t)r * W1]; }
      L.V = b; L.Vx.assign(n, 0.0); L.Vs.assign(n, 0.0); L.Vxs.assign(n, 0.0);
      for (int p = 0; p < n; p++) if (!fixed[p]) { L.Vx[p] = L.side[p]; L.Vs[p] = L.down[p] * hg[p % NH]; }
      for (int i = 1; i < NX; i++) for (int k = 0; k < NH; k++)   // (on the axis the energy is even in x: its slope in x is zero at every height, so the mixed slope is too)
      {
         const int p = i * NH + k;
         if (fixed[p]) continue;
         const int i0 = i - 1, i1 = std::min(NX - 1, i + 1), k0 = std::max(0, k - 1), k1 = std::min(NH - 1, k + 1);
         const double dVsDx = (L.Vs[i1 * NH + k] - L.Vs[i0 * NH + k]) / (xg[i1] - xg[i0]), dVxDs = (L.Vx[i * NH + k1] - L.Vx[i * NH + k0]) / std::log(hg[k1] / hg[k0]);
         L.Vxs[p] = 0.5 * (dVsDx + dVxDs);
      }
   }

   bool Model::Build(const Geometry &geom, const std::vector<double> &ampsIn)
   {
      g = geom; levels.clear(); xg.clear(); hg.clear();
      std::ostringstream lg; const double t0 = Now();
      std::vector<double> amps; for (double a : ampsIn) if (std::fabs(a) > 1e-9) amps.push_back(std::fabs(a));
      std::sort(amps.begin(), amps.end()); amps.erase(std::unique(amps.begin(), amps.end(), [](double a, double b) { return std::fabs(a - b) < 1e-6 * std::max(a, b); }), amps.end());
      if (amps.empty()) { log = "no operating level with a current"; return false; }
      Cells C; C.Setup(g); const double tC = Now();
      auto parallel = [](size_t n, const std::function<void(size_t)> &f) { std::vector<std::thread> pool; for (size_t i = 1; i < n; i++) pool.emplace_back(f, i); if (n > 0) f(0); for (auto &t : pool) t.join(); };
      // ---- per level, one thread each: coil-only Newton, bare table, on-axis reaction, plug state
      levels.resize(amps.size());
      std::vector<std::ostringstream> llog(amps.size()); std::vector<double> lt(amps.size()), lf(amps.size()), tbs(amps.size()), tas(amps.size()); std::vector<int> it0s(amps.size()), nfacts(amps.size());
      parallel(amps.size(), [&](size_t li)
      {
         const double tl = Now(); Level &L = levels[li]; L.amps = amps[li]; const double NI = g.turns * L.amps;
         Newton nw; nw.Init(C); std::vector<double> M0(C.n2, 0.0), Hc(C.n2); for (int i = 0; i < C.n2; i++) Hc[i] = NI * C.Hc1[i];
         nw.CHORD = 1; it0s[li] = nw.Solve(Hc, M0, 50); nw.CHORD = 50; nw.since = 1 << 30;   // coil-only from cold: full Newton (a frozen Jacobian from M = 0 stalled at low current)
         L.reach = -1; BuildBare(C, g, NI, M0, L.bare, L.reach); const double tb = Now();
         std::vector<double> Mplug; OnAxis(C, g, NI, M0, nw, L.hR, L.R, L.sat, L.plugB, Mplug, llog[li]); const double ta0 = Now();
         if (!Mplug.empty()) { double rp = L.reach; BuildBare(C, g, NI, Mplug, L.plug, rp); }   // the core with the plug on the face, for the ball-ball model (same extent)
         const double ta = Now(); llog[li] << "    plug-state table " << (ta - ta0) * 1e3 << " ms\n";
         tbs[li] = tb - tl; tas[li] = ta - tb; nfacts[li] = nw.nfact; lt[li] = ta - tl;
      });
      const double tL = Now();
      // ---- geometry table (once per model, over the extent of the furthest bare table): equivalent height of the linear tubes and its slopes (bench_geom, analytic + visibility boundary term)
      double xmax = X_KNEE, hmax = 30e-3;
      for (const Level &L : levels) if (!L.bare.r.empty()) { xmax = std::max(xmax, L.bare.r.back() - g.ballR); hmax = std::max(hmax, L.bare.z.back() - 2 * g.ballR); }
      int nClamped = 0;
      {
         const double a = g.coreRadius, top = g.exposed, R = g.ballR, DS = R / MUB; const int nr = 60, nt = 96, nz = 60, nth = nt / 2;
         std::vector<double> px, py, pz, nx, ny, nzv, w;
         for (int i = 0; i < nr; i++) for (int k = 0; k < nth; k++) { const double r = (i + 0.5) / nr * a, t = (k + 0.5) / nt * 2 * PI;
            px.push_back(r * std::cos(t)); py.push_back(r * std::sin(t)); pz.push_back(0); nx.push_back(0); ny.push_back(0); nzv.push_back(1); w.push_back(2 * r * (a / nr) * (2 * PI / nt)); }
         for (int i = 0; i < nz; i++) for (int k = 0; k < nth; k++) { const double u = (i + 0.5) / nz, z = -top * (1 - std::cos(0.5 * PI * u)), dz = top * 0.5 * PI * std::sin(0.5 * PI * u) / nz, t = (k + 0.5) / nt * 2 * PI;
            px.push_back(a * std::cos(t)); py.push_back(a * std::sin(t)); pz.push_back(z); nx.push_back(std::cos(t)); ny.push_back(std::sin(t)); nzv.push_back(0); w.push_back(2 * dz * a * (2 * PI / nt)); }
         const size_t NE = w.size();
         auto Pgrad = [&](double cx, double cz, double &P, double &Px, double &Pz)
         {
            P = Px = Pz = 0;
            for (size_t e = 0; e < NE; e++) { const double dx = cx - px[e], dy = -py[e], dz = cz - pz[e]; if (dx * nx[e] + dy * ny[e] + dz * nzv[e] <= 0) continue;
               const double d = std::sqrt(dx * dx + dy * dy + dz * dz), q = 1.0 / (d - R + DS), q2 = w[e] * q * q / d; P += w[e] * q; Px -= q2 * dx; Pz -= q2 * dz; }
            if (cx > a && top > 0) { const double t0_ = std::acos(a / cx), dt0 = a / (cx * cx * std::sin(t0_)), sx = a * std::cos(t0_), sy = a * std::sin(t0_); double acc = 0; const int nb = 200;
               for (int i = 0; i < nb; i++) { const double z = -top * (i + 0.5) / nb, d = std::sqrt((cx - sx) * (cx - sx) + sy * sy + (cz - z) * (cz - z)); acc += (top / nb) / (d - R + DS); }
               Px += 2 * a * dt0 * acc; }
         };
         const int NA = 400; std::vector<double> ha(NA), Pa(NA), dPa(NA);
         for (int i = 0; i < NA; i++) { ha[i] = 1e-6 * std::pow(0.3 / 1e-6, i / (NA - 1.0)); double P, Px, Pz; Pgrad(0, R + ha[i], P, Px, Pz); Pa[i] = P; dPa[i] = Pz; }
         auto invert = [&](double P, double &hs, double &dP, bool &clamped) { const double lp = std::log(P); int i = 0; while (i < NA - 2 && std::log(Pa[i + 1]) > lp) i++;
            double u = (lp - std::log(Pa[i])) / (std::log(Pa[i + 1]) - std::log(Pa[i])); clamped = u < 0 || u > 1; u = std::min(std::max(u, 0.0), 1.0);
            hs = std::exp(std::log(ha[i]) + u * (std::log(ha[i + 1]) - std::log(ha[i]))); dP = dPa[i] + u * (dPa[i + 1] - dPa[i]); };
         for (int i = 0; i <= NX_NEAR; i++) xg.push_back(DX_NEAR * i);
         while (xg.back() < xmax - 1e-9) xg.push_back(std::min(xmax, X_KNEE + DX_FAR * (double)(xg.size() - NX_NEAR)));
         for (int k = 0; ; k++) { const double h = GAP_MIN * std::exp(HStepLog() * k); if (h >= hmax - 1e-9) { hg.push_back(hmax); break; } hg.push_back(h); }
         const int NX = (int)xg.size(), NH = (int)hg.size();
         hs.assign(NX * NH, 0); dhx.assign(NX * NH, 0); dhz.assign(NX * NH, 0);
         const size_t nth_ = std::max<size_t>(1, std::min<size_t>(8, std::thread::hardware_concurrency())); std::vector<int> clamp(nth_, 0);
         parallel(nth_, [&](size_t t) { for (int i = (int)t; i < NX; i += (int)nth_) for (int k = 0; k < NH; k++) { double P, Px, Pz, h_, dP; bool cl; Pgrad(xg[i], R + hg[k], P, Px, Pz); invert(P, h_, dP, cl);
            clamp[t] += cl; hs[i * NH + k] = h_; dhx[i * NH + k] = Px / dP; dhz[i * NH + k] = Pz / dP; } });
         for (int c_ : clamp) nClamped += c_;
      }
      const double tG = Now();
      // ---- per level, one thread each: force table over the geometry grid, rim rule. Beyond a level's bare table: no force (the engine's reach)
      BallRule ball; ball.Init(g.ballR);
      std::vector<int> nfixes(amps.size());
      parallel(amps.size(), [&](size_t li)
      {
         const double tf = Now(); Level &L = levels[li];
         const int NX = (int)xg.size(), NH = (int)hg.size(); L.down.assign(NX * NH, 0); L.side.assign(NX * NH, 0);
         for (int i = 0; i < NX; i++) for (int k = 0; k < NH; k++)
         {
            if (!L.bare.Inside(xg[i] + g.ballR, hg[k] + 2 * g.ballR)) continue;
            double Fx, Fz; ball.Force(L.bare, xg[i], g.ballR + hg[k], Fx, Fz);
            const double h_ = hs[i * NH + k], Rr = InterpR(L.hR, L.R, h_), s = InterpSat(L.hR, L.sat, h_);
            L.down[i * NH + k] = s * (-Fz + Rr * dhz[i * NH + k]);
            L.side[i * NH + k] = s * (-Fx + Rr * dhx[i * NH + k]);
         }
         // rim rule on the down force, per height row
         const double xin = std::max(0.0, g.coreRadius - RIM_IN), xout = g.coreRadius + RIM_OUT;
         auto atx = [&](int k, double x) { int i = 0; while (i < NX - 2 && xg[i + 1] < x) i++; const double u = (x - xg[i]) / (xg[i + 1] - xg[i]); return (1 - u) * L.down[i * NH + k] + u * L.down[(i + 1) * NH + k]; };
         int nfix = 0;
         for (int k = 0; k < NH && g.rimRule; k++)
         {
            double xp = xin, fp = atx(k, xin); const double fo = atx(k, xout);
            for (int i = 0; i < NX; i++)
            {
               if (xg[i] <= xin || xg[i] >= xout) continue;
               double &f = L.down[i * NH + k];
               if (f < std::min(fp, fo) || f > std::max(fp, fo)) { f = fp + (fo - fp) * (xg[i] - xp) / (xout - xp); nfix++; }
               xp = xg[i]; fp = f;
            }
         }
         BuildEnergy(xg, hg, g, L);
         nfixes[li] = nfix; lf[li] = Now() - tf;
      });
      for (size_t li = 0; li < levels.size(); li++)
         llog[li] << "    coil-only Newton " << it0s[li] << " steps, field reach " << levels[li].reach * 1e3 << " mm (B_GATE " << B_GATE * 1e3 << " mT), bare table " << levels[li].bare.r.size() << " x " << levels[li].bare.z.size() << " (" << tbs[li] * 1e3 << " ms), on-axis "
                  << tas[li] * 1e3 << " ms, force table " << lf[li] * 1e3 << " ms, rim rule corrected " << nfixes[li] << " points, Jacobian factorisations " << nfacts[li] << "\n";
      lg << "core cells " << C.nr << " x " << C.nz << " (" << (tC - t0) * 1e3 << " ms), levels " << (tL - tC) * 1e3 << " ms, geometry table " << xg.size() << " x " << hg.size() << " to "
         << xg.back() * 1e3 << " x " << hg.back() * 1e3 << " mm (" << (tG - tL) * 1e3 << " ms" << (nClamped ? ", " + std::to_string(nClamped) + " points outside the on-axis range" : std::string()) << ")\n";
      for (size_t li = 0; li < levels.size(); li++)
         lg << "  level " << levels[li].amps << " A (" << g.turns * levels[li].amps << " At): " << (lt[li] + lf[li]) * 1e3 << " ms\n" << llog[li].str();
      lg << "total " << (Now() - t0) * 1e3 << " ms";
      log = lg.str();
      return true;
   }

   bool Model::Force(double x, double h, double I, double &down, double &side) const
   {
      down = side = 0;
      I = std::fabs(I); if (levels.empty() || I == 0) return false;
      h = std::max(h, GAP_MIN);
      // the two levels around I (0 A = no force; above the top level: linear extrapolation from the top two, or proportional with one level)
      size_t hi = 0; while (hi < levels.size() && levels[hi].amps < I) hi++;
      const Level *lo_ = nullptr, *hi_ = nullptr; double Ilo = 0, Ihi = 0;
      if (hi == 0) { hi_ = &levels[0]; Ihi = hi_->amps; }
      else if (hi == levels.size()) { hi_ = &levels.back(); Ihi = hi_->amps; if (levels.size() > 1) { lo_ = &levels[levels.size() - 2]; Ilo = lo_->amps; } }
      else { lo_ = &levels[hi - 1]; Ilo = lo_->amps; hi_ = &levels[hi]; Ihi = hi_->amps; }
      const int NX = (int)xg.size(), NH = (int)hg.size();
      if (NX < 2 || NH < 2 || x > xg.back() || h > hg.back())
         return true;   // beyond the tables (the furthest bare table): no force
      // the grid cell by index (x: 0.5 mm steps to 30 mm, 1 mm beyond; h: one log ratio)
      const int i = std::min(NX - 2, x <= X_KNEE ? (int)(x / DX_NEAR) : NX_NEAR + (int)((x - X_KNEE) / DX_FAR));
      const int k = std::min(NH - 2, std::max(0, (int)(std::log(h / GAP_MIN) / HStepLog())));
      const double u = (x - xg[i]) / (xg[i + 1] - xg[i]), v = std::log(h / hg[k]) / std::log(hg[k + 1] / hg[k]);
      auto one = [&](const Level &L, double &d, double &s)
      {
         auto bil = [&](const std::vector<double> &T) { return (1 - u) * (1 - v) * T[i * NH + k] + u * (1 - v) * T[(i + 1) * NH + k] + (1 - u) * v * T[i * NH + k + 1] + u * v * T[(i + 1) * NH + k + 1]; };
         d = bil(L.down); s = bil(L.side);
      };
      double dh, sh; one(*hi_, dh, sh);
      if (!lo_) { const double f = I / Ihi; down = f * dh; side = f * sh; return true; }
      double dl, sl; one(*lo_, dl, sl);
      const double w = (I - Ilo) / (Ihi - Ilo); down = (1 - w) * dl + w * dh; side = (1 - w) * sl + w * sh;
      return true;
   }

   void Model::Energy(double x, double h, double I, Pull &p) const
   {
      p = Pull();
      I = std::fabs(I);
      const int NX = (int)xg.size(), NH = (int)hg.size();
      if (levels.empty() || I == 0 || NX < 2 || NH < 2 || x > xg.back() || h > hg.back())
         return;
      // the two levels around I and their shares, as Force
      size_t hi = 0; while (hi < levels.size() && levels[hi].amps < I) hi++;
      const Level *lo_ = nullptr, *hi_ = nullptr; double wl = 0, wh = 1;
      if (hi == 0) { hi_ = &levels[0]; wh = I / hi_->amps; }
      else if (hi == levels.size()) { hi_ = &levels.back(); if (levels.size() > 1) { lo_ = &levels[levels.size() - 2]; wh = (I - lo_->amps) / (hi_->amps - lo_->amps); wl = 1 - wh; } else wh = I / hi_->amps; }
      else { lo_ = &levels[hi - 1]; hi_ = &levels[hi]; wh = (I - lo_->amps) / (hi_->amps - lo_->amps); wl = 1 - wh; }
      const double hg0 = std::max(h, GAP_MIN);
      const int i = std::min(NX - 2, x <= X_KNEE ? (int)(x / DX_NEAR) : NX_NEAR + (int)((x - X_KNEE) / DX_FAR));
      const int k = std::min(NH - 2, std::max(0, (int)(std::log(hg0 / GAP_MIN) / HStepLog())));
      const double dx = xg[i + 1] - xg[i], ds = std::log(hg[k + 1] / hg[k]), u = (x - xg[i]) / dx, v = std::log(hg0 / hg[k]) / ds;
      // the cubic Hermite shapes: A value at each end, C slope at each end, with their first and second derivatives
      auto shapes = [](double t, double A[2], double C[2], double dA[2], double dC[2], double ddA[2], double ddC[2])
      {
         const double t2 = t * t, t3 = t2 * t;
         A[0] = 2 * t3 - 3 * t2 + 1; A[1] = 3 * t2 - 2 * t3; C[0] = t3 - 2 * t2 + t; C[1] = t3 - t2;
         dA[0] = 6 * t2 - 6 * t; dA[1] = 6 * t - 6 * t2; dC[0] = 3 * t2 - 4 * t + 1; dC[1] = 3 * t2 - 2 * t;
         ddA[0] = 12 * t - 6; ddA[1] = 6 - 12 * t; ddC[0] = 6 * t - 4; ddC[1] = 6 * t - 2;
      };
      double Au[2], Cu[2], dAu[2], dCu[2], ddAu[2], ddCu[2], Av[2], Cv[2], dAv[2], dCv[2], ddAv[2], ddCv[2];
      shapes(u, Au, Cu, dAu, dCu, ddAu, ddCu);
      shapes(v, Av, Cv, dAv, dCv, ddAv, ddCv);
      double e = 0, eu = 0, ev = 0, euu = 0, euv = 0, evv = 0;
      auto one = [&](const Level &L, double w)
      {
         for (int c = 0; c < 2; c++) for (int d = 0; d < 2; d++)
         {
            const int q = (i + c) * NH + k + d;
            const double a = w * L.V[q], bx = w * L.Vx[q] * dx, bs = w * L.Vs[q] * ds, bm = w * L.Vxs[q] * dx * ds;
            e += Au[c] * Av[d] * a + Cu[c] * Av[d] * bx + Au[c] * Cv[d] * bs + Cu[c] * Cv[d] * bm;
            eu += dAu[c] * Av[d] * a + dCu[c] * Av[d] * bx + dAu[c] * Cv[d] * bs + dCu[c] * Cv[d] * bm;
            ev += Au[c] * dAv[d] * a + Cu[c] * dAv[d] * bx + Au[c] * dCv[d] * bs + Cu[c] * dCv[d] * bm;
            euu += ddAu[c] * Av[d] * a + ddCu[c] * Av[d] * bx + ddAu[c] * Cv[d] * bs + ddCu[c] * Cv[d] * bm;
            euv += dAu[c] * dAv[d] * a + dCu[c] * dAv[d] * bx + dAu[c] * dCv[d] * bs + dCu[c] * dCv[d] * bm;
            evv += Au[c] * ddAv[d] * a + Cu[c] * ddAv[d] * bx + Au[c] * ddCv[d] * bs + Cu[c] * ddCv[d] * bm;
         }
      };
      one(*hi_, wh);
      if (lo_) one(*lo_, wl);
      const double Vs = ev / ds, Vss = evv / (ds * ds);
      p.V = e; p.Vx = eu / dx; p.Vh = Vs / hg0; p.Vxx = euu / (dx * dx); p.Vxh = euv / (dx * ds) / hg0; p.Vhh = (Vss - Vs) / (hg0 * hg0);
      if (h < GAP_MIN)   // below the tables' first gap: that gap's down force goes on
      {
         const double dh = h - GAP_MIN;
         p.V += p.Vh * dh; p.Vx += p.Vxh * dh; p.Vhh = 0;
      }
   }

   bool Model::Field(double x, double z, double I, double &Br, double &Bz) const
   {
      Br = Bz = 0; const double sg = I < 0 ? -1.0 : 1.0; I = std::fabs(I);
      if (levels.empty() || I == 0 || !levels.back().bare.Inside(x, z)) return false;   // the top level reaches furthest
      size_t hi = 0; while (hi < levels.size() && levels[hi].amps < I) hi++;
      auto one = [&](const Level &L, double &r_, double &z_) { r_ = z_ = 0; if (!L.bare.Inside(x, std::max(z, L.bare.z.front()))) return;   // beyond a lower level's reach: under B_GATE there, taken as 0
         double v[2], vr[2], vz[2]; L.bare.Eval(x, std::max(z, L.bare.z.front()), v, vr, vz); r_ = MU0 * v[0]; z_ = MU0 * v[1]; };
      double rh, zh;
      if (hi == 0 || levels.size() == 1) { const Level &L = levels[std::min(hi, levels.size() - 1)]; one(L, rh, zh); const double f = I / L.amps; Br = sg * f * rh; Bz = sg * f * zh; return true; }
      if (hi == levels.size()) hi = levels.size() - 1;
      double rl, zl; one(levels[hi - 1], rl, zl); one(levels[hi], rh, zh);
      const double u = (I - levels[hi - 1].amps) / (levels[hi].amps - levels[hi - 1].amps);
      Br = sg * ((1 - u) * rl + u * rh); Bz = sg * ((1 - u) * zl + u * zh);
      return true;
   }
}
