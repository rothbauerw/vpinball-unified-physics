// license:GPLv3+
// Ball-ball forces through the balls' magnetisation (see magnetballs.h). Ported from the offline work in runs/2026-10-14-mirror (optA_recipe.py) and
// runs/2026-10-13-cpp-proto (lean.h pair frames: rotate -> coaxial translation -> rotate back; bench_recipe.cpp).
// Conventions (ballmag.h): potential psi, H = -grad psi; incident psi = sum b_nm r^n Y_nm (orthonormal, Condon-Shortley), scattered psi = sum a_nm r^-(n+1) Y_nm,
// a = -beta_n R^(2n+1) b; scaled b_bar = b R^n, a_bar = a / R^(n+1); force on ball i = -mu0 R Re[ sum a_i^H W dT a_l / R + db_bar^H W a_bar_i ] with W = 2n + 1.

#include "physics/magnetballs.h"
#include <algorithm>
#include <cmath>
#include <mutex>
#include <unordered_map>
#include <sstream>
#include <chrono>

namespace MagnetBalls
{
   namespace
   {
      constexpr double PI = 3.14159265358979323846, MU0 = 4e-7 * PI, MUB = 200.0;
      constexpr double BALL_R = 25 * 0.5397e-3;                     // the standard ball the tables are built for (m)
      constexpr int L4 = 4, LP = 8, L5 = 5, L9 = 9, LH = 24;
      constexpr int K4 = L4 * (L4 + 2), KP = LP * (LP + 2), K5 = L5 * (L5 + 2), K9 = L9 * (L9 + 2);
      constexpr double PAIR_FIX_MAX = 12e-3;                         // largest surface gap of the pair operator table (pair fix range, Settings::pairFixGap)
      constexpr double GMIN = 0.1e-3;                                // contact gap convention (magnetcore: tables start at 0.1 mm)
      constexpr double PLUG_AXIS = 1e-3;                             // the plug gets the on-axis order-8 drive within this offset from the axis
      inline int kidx(int n, int m) { return n * n - 1 + m + n; }
      inline int lmi(int l, int m) { return l * (l + 1) + m; }
      inline double beta(int n) { return n * (MUB - 1) / (n * MUB + n + 1); }
      inline double len(const V3 &a) { return std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }
      inline V3 sub(const V3 &a, const V3 &b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }

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

      // orthonormal Y_lm, l <= Lmax (ballmag.h Harmonics)
      void Ylm(int Lmax, double ct, double st, double phi, cd *Y)
      {
         double p[16][16];
         p[0][0] = std::sqrt(1.0 / (4.0 * PI));
         for (int m = 1; m <= Lmax; m++) p[m][m] = -std::sqrt((2.0 * m + 1.0) / (2.0 * m)) * st * p[m - 1][m - 1];
         for (int m = 0; m < Lmax; m++) p[m + 1][m] = std::sqrt(2.0 * m + 3.0) * ct * p[m][m];
         for (int m = 0; m <= Lmax; m++)
            for (int l = m + 2; l <= Lmax; l++)
            {
               const double a = std::sqrt((4.0 * l * l - 1.0) / ((double)l * l - (double)m * m)), b = std::sqrt((((double)l - 1) * (l - 1) - (double)m * m) / (4.0 * ((double)l - 1) * (l - 1) - 1.0));
               p[l][m] = a * (ct * p[l - 1][m] - b * p[l - 2][m]);
            }
         const cd e1(std::cos(phi), std::sin(phi)); cd em(1.0, 0.0);
         for (int m = 0; m <= Lmax; m++) { for (int l = m; l <= Lmax; l++) { const cd y = p[l][m] * em; Y[lmi(l, m)] = y; if (m > 0) Y[lmi(l, -m)] = ((m & 1) ? -1.0 : 1.0) * std::conj(y); } em *= e1; }
      }

      // pair frames (lean.h, rotation convention 1 as validated there): Wigner rotation, coaxial translation, force sandwich in the frame
      struct Lean
      {
         int L, K; double R; std::vector<int> n_, m_; std::vector<double> W_, f_, co_, cz_, cp_, cm_; std::vector<std::vector<int>> rowsByM_;
         // Wigner small-d terms per (n, m', m), precomputed (build step B6, as runs/2026-10-13-cpp-proto/lean_fast.h): coefficient (sign, square roots, factorials),
         // powers of cos(beta/2) and sin(beta/2); a frame is then multiply-adds only
         std::vector<double> tc_; std::vector<unsigned char> ta_, tb_; std::vector<int> ts_;
         Lean(int Lg, double Rb) : L(Lg), K(Lg * (Lg + 2)), R(Rb)
         {
            n_.resize(K); m_.resize(K); W_.resize(K);
            for (int n = 1; n <= L; n++) for (int m = -n; m <= n; m++) { const int k = kidx(n, m); n_[k] = n; m_[k] = m; W_[k] = 2.0 * n + 1; }
            f_.resize(4 * L + 8); f_[0] = 1; for (size_t i = 1; i < f_.size(); i++) f_[i] = f_[i - 1] * (double)i;
            co_.assign((size_t)K * K, 0); cz_.assign((size_t)K * K, 0); cp_.assign((size_t)K * K, 0); cm_.assign((size_t)K * K, 0);
            for (int r = 0; r < K; r++)
               for (int c = 0; c < K; c++)
               {
                  const int n1 = n_[r], m1 = m_[r], n = n_[c], m = m_[c], N = n + n1, M = m - m1;
                  if (std::abs(M) > 1) continue;   // only M = 0 (translation, axial force) and M = +-1 (transverse force) are used on the axis
                  const double kap = std::sqrt(f_[N + M] * f_[N - M] / (f_[n + m] * f_[n - m] * f_[n1 + m1] * f_[n1 - m1]));
                  const double cf = std::sqrt((2.0 * n + 1.0) / (2.0 * n1 + 1.0)) * kap * (((n1 + m1) & 1) ? -1.0 : 1.0) * std::sqrt(4.0 * PI / (2.0 * N + 1.0));
                  const double qn = (2.0 * N + 1.0) / (2.0 * N + 3.0), y1 = std::sqrt((2.0 * N + 3.0) / (4 * PI));
                  const size_t e = (size_t)r * K + c;
                  if (M == 0) { co_[e] = cf * std::sqrt((2.0 * N + 1.0) / (4 * PI)); cz_[e] = -cf * std::sqrt(qn * ((N + 1.0) * (N + 1.0))) * y1; }
                  if (M == -1) cp_[e] = cf * std::sqrt(qn * (N + M + 1.0) * (N + M + 2.0)) * y1;
                  if (M == 1) cm_[e] = -cf * std::sqrt(qn * (N - M + 1.0) * (N - M + 2.0)) * y1;
               }
            rowsByM_.resize(2 * L + 1); for (int k = 0; k < K; k++) rowsByM_[m_[k] + L].push_back(k);
            ts_.push_back(0);
            for (int n = 1; n <= L; n++) for (int mp = -n; mp <= n; mp++) for (int m = -n; m <= n; m++)
            {
               const int kmin = std::max(0, m - mp), kmax = std::min(n + m, n - mp); const double pre = std::sqrt(f_[n + mp] * f_[n - mp] * f_[n + m] * f_[n - m]);
               for (int k = kmin; k <= kmax; k++)
               {
                  tc_.push_back(pre * (((k - m + mp) & 1) ? -1.0 : 1.0) / (f_[n + m - k] * f_[k] * f_[mp - m + k] * f_[n - mp - k]));
                  ta_.push_back((unsigned char)(2 * n + m - mp - 2 * k)); tb_.push_back((unsigned char)(mp - m + 2 * k));
               }
               ts_.push_back((int)tc_.size());
            }
         }
         // the rotation of a frame: D_{m'm} = e^{-i m' alpha} d_{m'm}(beta), kept as the real small-d blocks per degree and the row phases
         struct Frame { std::vector<double> d; cd ph[64]; double up[64]; double Rm[9]; };
         int dOff(int n) const { return (n - 1) * (4 * (n - 1) * (n - 1) + 12 * (n - 1) + 11) / 3; }   // sum of (2k+1)^2 for k < n
         void frame(const V3 &Dv, Frame &F) const
         {
            const double rho = len(Dv), ct = Dv.z / rho, bt = std::acos(std::max(-1.0, std::min(1.0, ct))), al = std::atan2(Dv.y, Dv.x);
            const double c = std::cos(bt / 2), s = std::sin(bt / 2); double cpw[64], spw[64]; cpw[0] = spw[0] = 1; for (int k = 1; k < 2 * L + 2; k++) { cpw[k] = cpw[k - 1] * c; spw[k] = spw[k - 1] * s; }
            F.d.resize(dOff(L + 1));
            for (int m = -L; m <= L; m++) F.ph[m + L] = std::polar(1.0, -m * al);
            int e = 0;
            for (int n = 1; n <= L; n++)
            {
               const int w = 2 * n + 1; double *B = &F.d[dOff(n)];
               for (int mp = -n; mp <= n; mp++)
                  for (int m = -n; m <= n; m++, e++)
                  {
                     double sum = 0; for (int t = ts_[e]; t < ts_[e + 1]; t++) sum += tc_[t] * cpw[ta_[t]] * spw[tb_[t]];
                     B[(mp + n) * w + (m + n)] = sum;
                  }
            }
            const double u = R / rho; F.up[0] = 1; for (int k = 1; k < 2 * L + 4; k++) F.up[k] = F.up[k - 1] * u;
            const double sa = std::sin(al), ca = std::cos(al), sb = std::sin(bt), cb = std::cos(bt);
            const double M3[9] = { ca * cb, -sa, ca * sb, sa * cb, ca, sa * sb, -sb, 0, cb }; for (int i = 0; i < 9; i++) F.Rm[i] = M3[i];
         }
         void rot(const Frame &F, const cd *in, cd *out, bool back) const   // global -> frame (back = false) or frame -> global
         {
            for (int n = 1; n <= L; n++)
            {
               const int w = 2 * n + 1, o = n * n - 1; const double *B = &F.d[dOff(n)]; const cd *P = &F.ph[L - n];
               if (back)   // D in: out_i = ph_i sum_j d_ij in_j
                  for (int i = 0; i < w; i++) { cd acc = 0; const double *row = B + i * w; for (int j = 0; j < w; j++) acc += row[j] * in[o + j]; out[o + i] = P[i] * acc; }
               else        // D^H in: out_i = sum_j d_ji conj(ph_j) in_j
               {
                  cd x[64]; for (int j = 0; j < w; j++) x[j] = std::conj(P[j]) * in[o + j];
                  for (int i = 0; i < w; i++) { cd acc = 0; for (int j = 0; j < w; j++) acc += B[j * w + i] * x[j]; out[o + i] = acc; }
               }
            }
         }
         void coax(const Frame &F, const cd *in, cd *out, int dir) const
         {
            for (int m = -L; m <= L; m++)
            {
               const auto &rows = rowsByM_[m + L];
               for (int r : rows) { cd acc = 0; for (int c : rows) { const int N = n_[r] + n_[c]; const double sg = (dir < 0 && (N & 1)) ? -1.0 : 1.0; acc += sg * co_[(size_t)r * K + c] * F.up[N + 1] * in[c]; } out[r] += acc; }
            }
         }
         // acc += T(c_i - c_j) a_j; frame built for D = c_i - c_j (dir +1) or c_j - c_i (dir -1)
         void translate(const Frame &F, int dir, const cd *aj, cd *acc, cd *t1, cd *t2) const
         {
            rot(F, aj, t1, false); for (int k = 0; k < K; k++) t2[k] = 0; coax(F, t1, t2, dir); rot(F, t2, t1, true); for (int k = 0; k < K; k++) acc[k] += t1[k];
         }
         // force sandwich (ai, al in the frame; D = c_i - c_l along +z * dir): global-frame (x, y, z) of sum conj(ai) W dT al / R
         void sandwich(const Frame &F, int dir, const cd *ai, const cd *al, double out[3]) const
         {
            cd sx = 0, sy = 0, sz = 0;
            for (int r = 0; r < K; r++)
            {
               cd dz = 0, dp = 0, dm = 0; const int m1 = m_[r];
               for (int c : rowsByM_[m1 + L]) { const int N = n_[r] + n_[c]; const double sg = (dir < 0 && ((N + 1) & 1)) ? -1.0 : 1.0; dz += sg * cz_[(size_t)r * K + c] * F.up[N + 2] * al[c]; }
               if (m1 - 1 >= -L) for (int c : rowsByM_[m1 - 1 + L]) { const int N = n_[r] + n_[c]; const double sg = (dir < 0 && ((N + 1) & 1)) ? -1.0 : 1.0; dp += sg * cp_[(size_t)r * K + c] * F.up[N + 2] * al[c]; }
               if (m1 + 1 <= L) for (int c : rowsByM_[m1 + 1 + L]) { const int N = n_[r] + n_[c]; const double sg = (dir < 0 && ((N + 1) & 1)) ? -1.0 : 1.0; dm += sg * cm_[(size_t)r * K + c] * F.up[N + 2] * al[c]; }
               const cd cr = std::conj(ai[r]) * W_[r];
               sx += cr * 0.5 * (dp + dm); sy += cr * (dp - dm) * cd(0, -0.5); sz += cr * dz;
            }
            const double fx = sx.real() / R, fy = sy.real() / R, fz = sz.real() / R; const double *M = F.Rm;
            out[0] = M[0] * fx + M[1] * fy + M[2] * fz; out[1] = M[3] * fx + M[4] * fy + M[5] * fz; out[2] = M[6] * fx + M[7] * fy + M[8] * fz;
         }
      };

      // ---- pair operators (per gap): two balls on the frame axis, ball 1 at the origin, ball 2 at +z. Input x_m = drive of ball 1 then ball 2, n = max(1,|m|)..8.
      //   Op_m : (hi - lo) correction of both balls, n <= 8 (hi = order 24 inside the pair, lo = order-4 exchange, each ball's own response full);
      //   Qz_m, Qp_m (to x_{m-1}), Qm_m (to x_{m+1}): force on ball 2 of the hi pair, sum conj(a2) W dT a1 / R in the frame (m-coupled for the transverse part);
      //   kz   : axial field H_z at the contact midpoint (incident + both balls, hi), from x_0.
      struct PairOps
      {
         std::vector<double> gaps; std::vector<std::vector<double>> data;   // one flat block per gap, layout below
         int sP[2 * LP + 1], s[2 * LP + 1]; size_t oOp[2 * LP + 1], oQz[2 * LP + 1], oQp[2 * LP + 1], oQm[2 * LP + 1], oKz = 0, size = 0;
         void Layout()
         {
            size_t o = 0;
            for (int m = -LP; m <= LP; m++) { const int i = m + LP; sP[i] = LP - std::max(1, std::abs(m)) + 1; s[i] = 2 * sP[i]; }
            for (int m = -LP; m <= LP; m++)
            {
               const int i = m + LP; oOp[i] = o; o += (size_t)s[i] * s[i]; oQz[i] = o; o += (size_t)s[i] * s[i];
               oQp[i] = o; if (m - 1 >= -LP) o += (size_t)s[i] * s[i - 1];
               oQm[i] = o; if (m + 1 <= LP) o += (size_t)s[i] * s[i + 1];
            }
            oKz = o; o += s[LP]; size = o;
         }
         void Build(double R, const std::vector<double> &gapList = {})   // gapList: offline tests only (empty = the engine table)
         {
            Layout(); Lean lh(LH, R);
            gaps = gapList;
            if (gaps.empty())   // 51 gaps, log-spaced 0.1 .. 12 mm (x1.10): coarser tables (13 gaps, x1.45) erred by 10-40 mN at 2-8 mm, as much as the fix (B2 test)
               for (int k = 0; k <= 50; k++)
                  gaps.push_back(GMIN * std::pow(PAIR_FIX_MAX / GMIN, k / 50.0));
            data.assign(gaps.size(), std::vector<double>(size, 0.0));
            for (size_t gi = 0; gi < gaps.size(); gi++)
            {
               std::vector<double> &d = data[gi]; const double D = 2 * R + gaps[gi], u = R / D; double up[64]; up[0] = 1; for (int k = 1; k < 64; k++) up[k] = up[k - 1] * u;
               std::vector<std::vector<double>> A1(2 * LP + 1), A2(2 * LP + 1); std::vector<int> hH(2 * LP + 1);
               for (int m = -LP; m <= LP; m++)
               {
                  const int mi = m + LP, nlo = std::max(1, std::abs(m)), h = LH - nlo + 1, nP = sP[mi], S = s[mi]; hH[mi] = h;
                  const auto &rows = lh.rowsByM_[m + LH];
                  auto solve = [&](bool lo, std::vector<double> &X)   // X (2h x S) = M^-1 RHS
                  {
                     const int n2 = 2 * h; std::vector<double> M((size_t)n2 * n2, 0.0); X.assign((size_t)n2 * S, 0.0);
                     for (int r = 0; r < h; r++)
                     {
                        M[(size_t)r * n2 + r] = 1; M[(size_t)(h + r) * n2 + h + r] = 1; const double br = beta(lh.n_[rows[r]]);
                        for (int c = 0; c < h; c++)
                        {
                           const int N = lh.n_[rows[r]] + lh.n_[rows[c]]; if (lo && (lh.n_[rows[r]] > L4 || lh.n_[rows[c]] > L4)) continue;
                           const double t = lh.co_[(size_t)rows[r] * lh.K + rows[c]] * up[N + 1];
                           M[(size_t)r * n2 + h + c] += br * t * ((N & 1) ? -1.0 : 1.0);   // field at ball 1 from ball 2: D = c1 - c2 along -z (dir -1)
                           M[(size_t)(h + r) * n2 + c] += br * t;                          // field at ball 2 from ball 1: dir +1
                        }
                     }
                     for (int q = 0; q < nP; q++) { X[(size_t)q * S + q] = -beta(nlo + q); X[(size_t)(h + q) * S + nP + q] = -beta(nlo + q); }
                     std::vector<int> piv(n2);
                     for (int k = 0; k < n2; k++)
                     {
                        int p = k; for (int i = k + 1; i < n2; i++) if (std::fabs(M[(size_t)i * n2 + k]) > std::fabs(M[(size_t)p * n2 + k])) p = i;
                        if (p != k) { for (int j = 0; j < n2; j++) std::swap(M[(size_t)k * n2 + j], M[(size_t)p * n2 + j]); for (int j = 0; j < S; j++) std::swap(X[(size_t)k * S + j], X[(size_t)p * S + j]); }
                        const double dk = M[(size_t)k * n2 + k];
                        for (int i = k + 1; i < n2; i++) { const double f = M[(size_t)i * n2 + k] / dk; if (f == 0) continue; for (int j = k; j < n2; j++) M[(size_t)i * n2 + j] -= f * M[(size_t)k * n2 + j]; for (int j = 0; j < S; j++) X[(size_t)i * S + j] -= f * X[(size_t)k * S + j]; }
                     }
                     for (int k = n2 - 1; k >= 0; k--) for (int j = 0; j < S; j++) { double v = X[(size_t)k * S + j]; for (int i = k + 1; i < n2; i++) v -= M[(size_t)k * n2 + i] * X[(size_t)i * S + j]; X[(size_t)k * S + j] = v / M[(size_t)k * n2 + k]; }
                  };
                  std::vector<double> Xh, Xl; solve(false, Xh); solve(true, Xl);
                  double *Op = &d[oOp[mi]];
                  for (int b = 0; b < 2; b++) for (int q = 0; q < nP; q++) for (int j = 0; j < S; j++) Op[(size_t)(b * nP + q) * S + j] = Xh[(size_t)(b * h + q) * S + j] - Xl[(size_t)(b * h + q) * S + j];
                  A1[mi].assign(Xh.begin(), Xh.begin() + (size_t)h * S); A2[mi].assign(Xh.begin() + (size_t)h * S, Xh.end());
               }
               // force forms: Q = A2^T W D A1' (dir +1, u^(N+2))
               for (int m = -LP; m <= LP; m++)
               {
                  const int mi = m + LP, h = hH[mi], S = s[mi]; const auto &rows = lh.rowsByM_[m + LH];
                  auto form = [&](int m2, const std::vector<double> &cf, double *Q)
                  {
                     const int mj = m2 + LP, h2 = hH[mj], S2 = s[mj]; const auto &rows2 = lh.rowsByM_[m2 + LH];
                     std::vector<double> DA((size_t)h * S2, 0.0);   // W D A1(m2)
                     for (int r = 0; r < h; r++) for (int c = 0; c < h2; c++)
                     {
                        const size_t e = (size_t)rows[r] * lh.K + rows2[c]; const double v = cf[e]; if (v == 0) continue; const int N = lh.n_[rows[r]] + lh.n_[rows2[c]];
                        const double t = lh.W_[rows[r]] * v * up[N + 2]; for (int j = 0; j < S2; j++) DA[(size_t)r * S2 + j] += t * A1[mj][(size_t)c * S2 + j];
                     }
                     for (int i = 0; i < S; i++) for (int j = 0; j < S2; j++) { double acc = 0; for (int r = 0; r < h; r++) acc += A2[mi][(size_t)r * S + i] * DA[(size_t)r * S2 + j]; Q[(size_t)i * S2 + j] = acc; }
                  };
                  form(m, lh.cz_, &d[oQz[mi]]);
                  if (m - 1 >= -LP) form(m - 1, lh.cp_, &d[oQp[mi]]);
                  if (m + 1 <= LP) form(m + 1, lh.cm_, &d[oQm[mi]]);
               }
               // axial field at the midpoint (frame z = R + g/2 from ball 1; ball 2 below it), from x_0
               {
                  const int mi = LP, h = hH[mi], S = s[mi], nP = sP[mi]; const double rr = R + 0.5 * gaps[gi]; double *kz = &d[oKz];
                  for (int j = 0; j < S; j++) kz[j] = 0;
                  for (int q = 0; q < nP; q++) { const int n = q + 1; kz[q] += -n * std::pow(rr / R, n - 1) / R * std::sqrt((2 * n + 1) / (4 * PI)); }   // incident on ball 1
                  for (int r = 0; r < h; r++)
                  {
                     const int n = r + 1; const double y0 = std::sqrt((2 * n + 1) / (4 * PI)), f = (n + 1) * std::pow(R / rr, n + 2) / R * y0;
                     for (int j = 0; j < S; j++) kz[j] += f * A1[mi][(size_t)r * S + j] - f * ((n & 1) ? -1.0 : 1.0) * A2[mi][(size_t)r * S + j];
                  }
               }
            }
         }
      };
      const PairOps &Ops()
      {
         static PairOps ops; static std::once_flag once;
         std::call_once(once, [] { ops.Build(BALL_R); });
         return ops;
      }

      // ball-ball contact saturation F_sat / F_lin vs the linear gap field (T) per gap: FEM, uniform-field pairs (neck_drive.json 'U', runs/2026-10-14-mirror)
      struct BBRow { double gap; std::vector<std::pair<double, double>> pts; };
      const std::vector<BBRow> &BBData()
      {
         static const std::vector<BBRow> d = {
            { 0.05e-3, { { 1.454, 0.8646 }, { 2.908, 0.6488 }, { 4.362, 0.5052 }, { 7.270, 0.3488 }, { 14.539, 0.1951 }, { 29.079, 0.0968 }, { 43.618, 0.0581 } } },
            { 0.10e-3, { { 0.972, 0.9528 }, { 1.945, 0.8430 }, { 2.917, 0.7103 }, { 4.862, 0.5202 }, { 9.724, 0.3027 }, { 19.449, 0.1550 }, { 29.173, 0.0956 } } },
            { 0.25e-3, { { 0.498, 0.9873 }, { 0.995, 0.9684 }, { 1.493, 0.9377 }, { 2.488, 0.8196 }, { 4.975, 0.5404 }, { 9.951, 0.2908 }, { 14.926, 0.1811 } } },
            { 0.50e-3, { { 0.287, 0.9944 }, { 0.574, 0.9871 }, { 0.860, 0.9782 }, { 1.434, 0.9501 }, { 2.868, 0.7692 }, { 5.735, 0.4552 }, { 8.603, 0.2919 } } },
            { 1.00e-3, { { 0.165, 0.9972 }, { 0.330, 0.9938 }, { 0.495, 0.9900 }, { 0.825, 0.9806 }, { 1.650, 0.9334 }, { 3.301, 0.6743 }, { 4.951, 0.4579 } } } };
         return d;
      }
      double RowPhi(const BBRow &r, double B)   // log-B interpolation, linear extrapolation in log B at both ends (as sat_rules.py)
      {
         const auto &p = r.pts; B = std::max(B, 1e-6); const double lb = std::log(B);
         size_t i = 0; if (lb > std::log(p.front().first)) { i = 0; while (i + 2 < p.size() && lb > std::log(p[i + 1].first)) i++; }
         const double x0 = std::log(p[i].first), x1 = std::log(p[i + 1].first), t = (lb - x0) / (x1 - x0);
         return std::min(1.0, std::max(0.02, p[i].second + t * (p[i + 1].second - p[i].second)));
      }
      double PhiBB(double gap, double B)
      {
         const auto &d = BBData(); gap = std::max(gap, d.front().gap);
         if (gap >= 3e-3) return 1.0;
         for (size_t i = 0; i + 1 < d.size(); i++)
            if (gap <= d[i + 1].gap) { const double u = std::log(gap / d[i].gap) / std::log(d[i + 1].gap / d[i].gap); return (1 - u) * RowPhi(d[i], B) + u * RowPhi(d[i + 1], B); }
         const double u = std::log(gap / d.back().gap) / std::log(3e-3 / d.back().gap); return (1 - u) * RowPhi(d.back(), B) + u;
      }

      // log interpolation helpers over the on-axis gaps (as magnetcore InterpR / InterpSat)
      double InterpR(const std::vector<double> &h, const std::vector<double> &R, double g)
      {
         g = std::max(g, h.front());
         if (g >= h.back()) { const double k = std::log(std::max(R[R.size() - 1], 1e-30) / std::max(R[R.size() - 2], 1e-30)) / std::log(h.back() / h[h.size() - 2]); return R.back() * std::pow(g / h.back(), k); }
         size_t i = 0; while (i + 2 < h.size() && g > h[i + 1]) i++;
         const double u = std::log(g / h[i]) / std::log(h[i + 1] / h[i]);
         if (R[i] > 0 && R[i + 1] > 0) return std::exp((1 - u) * std::log(R[i]) + u * std::log(R[i + 1]));
         return (1 - u) * R[i] + u * R[i + 1];
      }
      double InterpLin(const std::vector<double> &h, const std::vector<double> &s, double g)
      {
         g = std::max(g, h.front()); if (g >= h.back()) return s.back();
         size_t i = 0; while (i + 2 < h.size() && g > h[i + 1]) i++;
         const double u = std::log(g / h[i]) / std::log(h[i + 1] / h[i]); return (1 - u) * s[i] + u * s[i + 1];
      }

      // the levels around current I and their weights (as MagnetCore::Model::Force: proportional below the first level, linear above the top two)
      int LevelWeights(const MagnetCore::Model &m, double I, int idx[2], double w[2])
      {
         I = std::fabs(I); const auto &L = m.levels; if (L.empty() || I == 0) return 0;
         size_t hi = 0; while (hi < L.size() && L[hi].amps < I) hi++;
         if (hi == 0) { idx[0] = 0; w[0] = I / L[0].amps; return 1; }
         if (hi == L.size()) { if (L.size() == 1) { idx[0] = 0; w[0] = I / L[0].amps; return 1; } hi = L.size() - 1; }
         const double u = (I - L[hi - 1].amps) / (L[hi].amps - L[hi - 1].amps); idx[0] = (int)hi - 1; w[0] = 1 - u; idx[1] = (int)hi; w[1] = u; return 2;
      }

      // gradient of a drive about the centre: d b_bar_(n-1) / dc from b_bar_n (regular ladder, ballmag.h gz/gp/gm), up to degree nmax
      void Ladder(const cd *b, int nmax, double R, cd *db /* 3 x KP */)
      {
         for (int k = 0; k < 3 * KP; k++) db[k] = 0;
         std::vector<cd> dp(KP, 0.0), dm(KP, 0.0);
         for (int n = 2; n <= nmax + 1; n++)
            for (int m = -n; m <= n; m++)
            {
               const cd v = b[kidx(n, m)] / R; if (v == cd(0, 0)) continue; const double qn = (2.0 * n + 1.0) / (2.0 * n - 1.0);
               if (std::abs(m) <= n - 1) db[2 * KP + kidx(n - 1, m)] += std::sqrt(qn * (double)(n * n - m * m)) * v;
               if (m + 1 <= n - 1) dp[kidx(n - 1, m + 1)] += std::sqrt(qn * (double)(n - m) * (n - m - 1)) * v;
               if (m - 1 >= -(n - 1)) dm[kidx(n - 1, m - 1)] += -std::sqrt(qn * (double)(n + m) * (n + m - 1)) * v;
            }
         for (int k = 0; k < KP; k++) { db[k] = 0.5 * (dp[k] + dm[k]); db[KP + k] = (dp[k] - dm[k]) * cd(0, -0.5); }
      }
   }

   V3 ReaimPair(const V3 &v, const V3 &u0, const V3 &u1)
   {
      const V3 k = { u0.y * u1.z - u0.z * u1.y, u0.z * u1.x - u0.x * u1.z, u0.x * u1.y - u0.y * u1.x };
      const double c = u0.x * u1.x + u0.y * u1.y + u0.z * u1.z, s2 = k.x * k.x + k.y * k.y + k.z * k.z;
      if (s2 <= 1e-24) return v;   // (an exact reversal of the line inside one step does not happen: the balls would pass through each other)
      const double kv = k.x * v.x + k.y * v.y + k.z * v.z, f = (1 - c) / s2;
      const V3 kx = { k.y * v.z - k.z * v.y, k.z * v.x - k.x * v.z, k.x * v.y - k.y * v.x };
      return { v.x * c + kx.x + k.x * kv * f, v.y * c + kx.y + k.y * kv * f, v.z * c + kx.z + k.z * kv * f };
   }
   V3 ReaimMagnet(const V3 &v, const V3 &n, double a) { return ReaimMagnetCS(v, n, std::cos(a), std::sin(a)); }
   V3 ReaimMagnetCS(const V3 &v, const V3 &n, double c, double s)
   {
      const double nv = n.x * v.x + n.y * v.y + n.z * v.z;
      const V3 nx = { n.y * v.z - n.z * v.y, n.z * v.x - n.x * v.z, n.x * v.y - n.y * v.x };
      return { v.x * c + nx.x * s + n.x * nv * (1 - c), v.y * c + nx.y * s + n.y * nv * (1 - c), v.z * c + nx.z * s + n.z * nv * (1 - c) };
   }

   // ---- drive tables -------------------------------------------------------------------------------------------------------------------------------------------
   void DriveTable::Build(const MagnetCore::BareTable &T, double R)
   {
      x.clear(); h.clear();
      // grid out to the table's edge (MagnetCore: the level's field reach + one ball diameter): the projection sphere and the ball stay inside
      const double xmax = std::max(30e-3, T.r.back() - R), hmax = std::max(10 * GMIN, T.z.back() - 2 * R);
      for (int i = 0; i <= 60; i++) x.push_back(0.5e-3 * i);
      for (double xx = 30e-3 + 2.5e-3; xx < xmax + 1e-9; xx += 2.5e-3) x.push_back(xx);
      const int nh = std::max(2, (int)std::ceil(39.0 * std::log(hmax / GMIN) / std::log(600.0)) + 1);   // as dense as the old 40 points over 0.1 .. 60 mm
      for (int i = 0; i < nh; i++) h.push_back(GMIN * std::pow(hmax / GMIN, i / (nh - 1.0)));
      // projection on a sphere of radius Rs: b_nm = -1/(n Rs^(n-1)) int H_r conj(Y_nm) dOmega, b_bar = b R^n
      const int nth = 8, nph = 16; const double Rs = 0.5 * R; std::vector<double> mu, w; LegGauss(nth, mu, w);
      struct Pt { double nx, ny, nz; }; std::vector<Pt> pts; std::vector<cd> proj; cd Y[64];
      for (int i = 0; i < nth; i++) for (int j = 0; j < nph; j++)
      {
         const double ct = mu[i], st = std::sqrt(std::max(0.0, 1 - ct * ct)), ph = 2 * PI * j / nph; pts.push_back({ st * std::cos(ph), st * std::sin(ph), ct });
         Ylm(L5, ct, st, ph, Y); const double wt = w[i] * 2 * PI / nph;
         for (int n = 1; n <= L5; n++) for (int m = -n; m <= n; m++) proj.push_back(-wt * std::conj(Y[lmi(n, m)]) / (n * std::pow(Rs, n - 1)) * std::pow(R, n));
      }
      const size_t np = pts.size(); b.assign(x.size() * h.size() * K5, cd(0, 0)); std::vector<double> Hr(np);
      for (size_t ix = 0; ix < x.size(); ix++) for (size_t ih = 0; ih < h.size(); ih++)
      {
         const double cx = x[ix], cz = R + h[ih]; if (!T.Inside(cx + R, cz + R)) continue;
         for (size_t p = 0; p < np; p++)
         {
            const double X = cx + Rs * pts[p].nx, Yy = Rs * pts[p].ny, Z = cz + Rs * pts[p].nz, rho = std::hypot(X, Yy); double v[2], vr[2], vz[2]; T.Eval(rho, Z, v, vr, vz);
            const double hx = rho > 1e-12 ? v[0] * X / rho : 0, hy = rho > 1e-12 ? v[0] * Yy / rho : 0; Hr[p] = hx * pts[p].nx + hy * pts[p].ny + v[1] * pts[p].nz;
         }
         cd *o = &b[(ix * h.size() + ih) * K5];
         for (size_t p = 0; p < np; p++) { const cd *pr = &proj[p * K5]; for (int k = 0; k < K5; k++) o[k] += Hr[p] * pr[k]; }
      }
   }
   bool DriveTable::Eval(double rho, double gap, cd *b5) const
   {
      if (b.empty() || rho > x.back() || gap > h.back()) return false;
      gap = std::max(gap, h.front());
      int i = 0; while (i < (int)x.size() - 2 && x[i + 1] < rho) i++; int k = 0; while (k < (int)h.size() - 2 && h[k + 1] < gap) k++;
      const double u = (rho - x[i]) / (x[i + 1] - x[i]), v = std::log(gap / h[k]) / std::log(h[k + 1] / h[k]); const size_t nh = h.size();
      const cd *b00 = &b[(i * nh + k) * K5], *b10 = &b[((i + 1) * nh + k) * K5], *b01 = b00 + K5, *b11 = b10 + K5;
      for (int q = 0; q < K5; q++) b5[q] = (1 - u) * (1 - v) * b00[q] + u * (1 - v) * b10[q] + (1 - u) * v * b01[q] + u * v * b11[q];
      return true;
   }
   bool MagnetDrive::Build(const MagnetCore::Model &m)
   {
      const auto t0 = std::chrono::steady_clock::now(); core = &m; bare.assign(m.levels.size(), {}); plug.assign(m.levels.size(), {});
      for (size_t i = 0; i < m.levels.size(); i++) { bare[i].Build(m.levels[i].bare, m.g.ballR); if (!m.levels[i].plug.r.empty()) plug[i].Build(m.levels[i].plug, m.g.ballR); else plug[i] = bare[i]; }
      Ops();
      std::ostringstream lg; lg << "ball-ball drive tables " << m.levels.size() << " level(s) x 2, pair operators " << Ops().gaps.size() << " gaps, "
                                << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() << " ms"; log = lg.str();
      return !m.levels.empty();
   }

   // ---- solver -------------------------------------------------------------------------------------------------------------------------------------------------
   struct Solver::Impl
   {
      Lean l4, l8, l9; const double R = BALL_R;   // l9: rotation of a magnet's drive (order 9, for its gradient to order 8) from its frame to the table frame
      // op: the pair operators at this gap, a table row itself when the gap is at the table's end (2026-11-02), else opBuf; pf: partner fields of this pass (at j, at i)
      struct Pair { int i, j; bool touch; double gap; Lean::Frame f4, f8; const double *op = nullptr; std::vector<double> opBuf; std::vector<cd> x, pf; double key = 0; };
      std::vector<Pair> P;                          // this Evaluate's pairs; kept between calls so their buffers are reused (build step B6)
      struct Warm { std::vector<cd> a4, ac; V3 pos; };   // a ball's last answer: global order-4 coefficients, corrected coefficients (order 8 layout), its place then
      std::unordered_map<uint64_t, Warm> warm;    // by ball id (build step B5); a ball not in it starts from its own response to the drive
      std::vector<V3> lastPos; std::vector<cd> lastA; std::vector<uint64_t> lastIds;   // the last Evaluate: positions, coefficients (order 8 layout), ids
      Impl() : l4(L4, BALL_R), l8(LP, BALL_R), l9(L9, BALL_R) {}
   };
   Solver::Solver() : p(new Impl) { Ops(); }
   Solver::~Solver() { delete p; }
   double Solver::BallRadius() const { return BALL_R; }
   void Solver::ClearLast() { p->lastPos.clear(); p->lastA.clear(); p->lastIds.clear(); p->warm.clear(); }

   bool Solver::FieldFromBalls(const V3 &pt, uint64_t exclude, V3 &B) const
   {
      const Impl &I = *p; B = V3();
      if (I.lastPos.empty())
         return false;
      std::vector<cd> aj(K4), acc(K4, 0.0), t1(K4), t2(K4);
      for (size_t j = 0; j < I.lastPos.size(); j++)
      {
         if (j < I.lastIds.size() && I.lastIds[j] == exclude)
            continue;
         const V3 D = sub(pt, I.lastPos[j]);
         if (len(D) < 2 * I.R * 0.999)
            continue;   // inside the ball's own volume: not a field point of this expansion
         for (int k = 0; k < K4; k++) aj[k] = I.lastA[j * KP + k];
         Lean::Frame F; I.l4.frame(D, F); I.l4.translate(F, +1, aj.data(), acc.data(), t1.data(), t2.data());
      }
      // uniform part of the incident field at pt: psi = sum_m b_1m r Y_1m (orthonormal, Condon-Shortley), b = b_bar / R; H = -grad psi
      const double c = std::sqrt(3.0 / (4 * PI)) / I.R, s2 = 1.0 / std::sqrt(2.0);
      const cd bm = acc[kidx(1, -1)], b0 = acc[kidx(1, 0)], bp = acc[kidx(1, 1)];
      const double Hx = c * ((bp - bm) * s2).real(), Hy = c * ((bp + bm) * cd(0, s2)).real(), Hz = -c * b0.real();
      B = { MU0 * Hx, MU0 * Hy, MU0 * Hz };
      return true;
   }

   void Solver::Evaluate(const std::vector<Source> &src, const std::vector<V3> &pos, std::vector<V3> &dF, Stats *st, const std::vector<uint64_t> *ids, std::vector<Part> *parts)
   {
      Impl &I = *p; const int N = (int)pos.size(); const double R = I.R; const PairOps &ops = Ops();
      dF.assign(N, V3());
      if (parts) parts->clear();
      if (st) *st = Stats();
      if (N < 2) return;
      // ---- drive per ball: b (n <= 8; n <= 4 except the on-axis plug) and its gradient
      std::vector<cd> b8((size_t)N * KP, 0.0), db8((size_t)N * 3 * KP, 0.0); std::vector<int> plugOf(N, -1); std::vector<double> facePhi(N, 1.0); std::vector<char> high(N, 0);
      std::vector<cd> b9(K9), b9r(K9), d9(3 * KP), b5(K5), b5p(K5);
      for (size_t s = 0; s < src.size(); s++)
      {
         const Source &S = src[s]; if (!S.drive || !S.drive->core || S.amps == 0) continue; const MagnetCore::Model &M = *S.drive->core;
         int li[2]; double lw[2]; const int nl = LevelWeights(M, S.amps, li, lw); if (!nl) continue;
         const double sgn = S.amps < 0 ? -1.0 : 1.0;
         // the magnet's frame: z' = its axis, x' y' from the Wigner frame of the axis (the tables are axisymmetric, so any x' serves); identity for an upright magnet
         const double an = len(S.axis); const V3 ax = an > 0 ? V3 { S.axis.x / an, S.axis.y / an, S.axis.z / an } : V3 { 0, 0, 1 };
         const bool tilted = ax.z < 1.0 - 1e-12;
         Lean::Frame Fm; V3 ex { 1, 0, 0 }, ey { 0, 1, 0 }, ez { 0, 0, 1 };
         if (tilted) { I.l9.frame(ax, Fm); ex = { Fm.Rm[0], Fm.Rm[3], Fm.Rm[6] }; ey = { Fm.Rm[1], Fm.Rm[4], Fm.Rm[7] }; ez = { Fm.Rm[2], Fm.Rm[5], Fm.Rm[8] }; }
         auto local = [&](int i, double &lx, double &ly, double &lz) { const V3 d = sub(pos[i], S.origin); lx = d.x * ex.x + d.y * ex.y + d.z * ex.z; ly = d.x * ey.x + d.y * ey.y + d.z * ey.z; lz = d.x * ez.x + d.y * ez.y + d.z * ez.z; };
         // plug: the ball closest to the face over the core
         int pl = -1; double pg = 1e9, prho = 0;
         for (int i = 0; i < N; i++) { double lx, ly, lz; local(i, lx, ly, lz); const double rho = std::hypot(lx, ly);
            if (lz > 0 && rho < M.g.coreRadius && lz - R < pg) { pg = lz - R; pl = i; prho = rho; } }
         if (pl >= 0 && pg > M.levels[0].hR.back()) pl = -1;
         double wsum = 0; for (int q = 0; q < nl; q++) wsum += lw[q];
         if (pl >= 0) { plugOf[pl] = (int)s; double f = 0; for (int q = 0; q < nl; q++) f += lw[q] / wsum * InterpLin(M.levels[li[q]].hR, M.levels[li[q]].sat, std::max(pg, GMIN)); facePhi[pl] = f;
            if (prho >= PLUG_AXIS && st) st->offAxisPlug++; }
         for (int i = 0; i < N; i++)
         {
            double lx, ly, lz; local(i, lx, ly, lz); const double rho = std::hypot(lx, ly), ph = std::atan2(ly, lx), gap = lz - R;
            if (lz <= 0) continue;
            for (int k = 0; k < K9; k++) b9[k] = 0;
            const bool zon = (i == pl && prho < PLUG_AXIS); int nmax = L4;
            for (int q = 0; q < nl; q++)
            {
               const MagnetCore::Level &Lv = M.levels[li[q]]; const double w = lw[q] * sgn;
               if (zon)
               {
                  // zonal drive from the on-axis solve, interpolated in log gap
                  const auto &h = Lv.hR; const double g = std::max(gap, h.front()); size_t a = 0; while (a + 2 < h.size() && g > h[a + 1]) a++;
                  const double u = g >= h.back() ? 1.0 : std::log(g / h[a]) / std::log(h[a + 1] / h[a]); const size_t a1 = g >= h.back() ? h.size() - 1 : a + 1, a0 = g >= h.back() ? a1 : a;
                  for (int n = 1; n <= L9 && n < (int)Lv.plugB[a0].size(); n++) b9[kidx(n, 0)] += w * ((1 - u) * Lv.plugB[a0][n] + u * Lv.plugB[a1][n]);
                  nmax = LP;
               }
               else
               {
                  const bool okb = S.drive->bare[li[q]].Eval(rho, gap, b5.data()); if (!okb) continue;
                  double wb = 0; if (pl >= 0) wb = std::min(1.0, std::max(0.0, InterpR(Lv.hR, Lv.R, std::max(pg, GMIN)) / std::max(Lv.R.front(), 1e-30)));
                  if (wb > 0 && S.drive->plug[li[q]].Eval(rho, gap, b5p.data())) for (int k = 0; k < K5; k++) b5[k] = (1 - wb) * b5[k] + wb * b5p[k];
                  for (int n = 1; n <= L5; n++) for (int m = -n; m <= n; m++) b9[kidx(n, m)] += w * b5[kidx(n, m)] * std::polar(1.0, -m * ph);
               }
            }
            if (tilted) { I.l9.rot(Fm, b9.data(), b9r.data(), true); b9.swap(b9r); }   // magnet frame -> table frame (the gradient below is then in the table frame)
            Ladder(b9.data(), nmax, R, d9.data());
            if (nmax > L4) high[i] = 1;
            for (int n = 1; n <= nmax; n++) for (int m = -n; m <= n; m++) { const int k = kidx(n, m); b8[(size_t)i * KP + k] += b9[k]; for (int c = 0; c < 3; c++) db8[(size_t)i * 3 * KP + c * KP + k] += d9[c * KP + k]; }
         }
      }
      // ---- each ball alone in its drive (no exchange): a = -beta b; force from the drive gradient
      auto extF = [&](int i, const cd *a, V3 &F) { double s3[3]; for (int c = 0; c < 3; c++) { cd e = 0; const cd *g = &db8[(size_t)i * 3 * KP + c * KP]; for (int k = 0; k < KP; k++) e += std::conj(g[k]) * (2.0 * I.l8.n_[k] + 1) * a[k]; s3[c] = -MU0 * R * e.real(); } F = { s3[0], s3[1], s3[2] }; };
      std::vector<V3> Flone(N);
      { std::vector<cd> al(KP); for (int i = 0; i < N; i++) { for (int k = 0; k < KP; k++) al[k] = -beta(I.l8.n_[k]) * b8[(size_t)i * KP + k]; extF(i, al.data(), Flone[i]); } }
      // ---- every pair, frames; pairs within the pair-fix range get order-8 frames and interpolated pair operators
      const double fixGap = std::min(set.pairFixGap, PAIR_FIX_MAX);
      using Pair = Impl::Pair;
      std::vector<Pair> &P = I.P; std::vector<int> pidx((size_t)N * N, -1); size_t np = 0;
      for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++)
      {
         const V3 D = sub(pos[i], pos[j]); const double d = len(D);
         if (np == P.size()) P.emplace_back();
         Pair &q = P[np]; q.i = i; q.j = j; q.key = 0; q.gap = d - 2 * R; q.touch = q.gap < fixGap; I.l4.frame(D, q.f4);   // touch = pair-fixed (within the range)
         if (q.touch)
         {
            I.l8.frame(D, q.f8); const double g = std::min(std::max(q.gap, ops.gaps.front()), ops.gaps.back()); size_t a = 0; while (a + 2 < ops.gaps.size() && g > ops.gaps[a + 1]) a++;
            const double u = std::log(g / ops.gaps[a]) / std::log(ops.gaps[a + 1] / ops.gaps[a]);
            const double *A0 = ops.data[a].data(), *A1 = ops.data[a + 1].data();
            if (u <= 0) q.op = A0;        // at the table's end (touching: the gap is under its first row): the row itself, the same numbers as the blend
            else if (u >= 1) q.op = A1;
            else { q.opBuf.resize(ops.size); for (size_t e = 0; e < ops.size; e++) q.opBuf[e] = (1 - u) * A0[e] + u * A1[e]; q.op = q.opBuf.data(); }
            q.x.assign(2 * KP, 0.0); q.pf.assign(2 * KP, 0.0);
         }
         pidx[(size_t)i * N + j] = pidx[(size_t)j * N + i] = (int)np; np++;
      }
      P.resize(np);
      if (st) { st->balls = N; st->pairs = (int)P.size(); for (auto &q : P) st->pairFix += q.touch; }
      // ---- global order-4 solve: Gauss-Seidel rounds from each ball's last answer (by id), a new ball from its own response (a = -beta b)
      std::vector<uint64_t> idv(N); for (int i = 0; i < N; i++) idv[i] = (ids && (int)ids->size() == N) ? (*ids)[i] : (uint64_t)i;
      std::vector<const Impl::Warm *> wp(N, nullptr);
      bool cold = false;   // a ball with no last answer: today's fixed rounds first, no shortcut
      for (int i = 0; i < N; i++) { const auto it = I.warm.find(idv[i]); if (it != I.warm.end()) wp[i] = &it->second; else cold = true; }
      // The last answer turned with its ball round the magnet's axis (2026-11-02). The drive is axisymmetric, so balls turned together about the axis have the same
      // answer turned: a_nm -> a_nm e^(-i m dphi). dphi: the turn about the axis that best carries the ball and its pair-fixed partners from their last places to
      // these, each weighted by its distance from the axis squared (a ball near the axis has no angle of its own; its partners give it one, a lone ball there gets
      // almost none). Only with one energised upright magnet: with two, or a tilted one, there is no one axis of symmetry and nothing is turned. It sets the start
      // only: the rounds below stop on the force they still move, so a poor turn costs rounds, not the answer.
      std::vector<cd> ph((size_t)N * (2 * LP + 1), cd(1, 0));
      {
         int on = -1, nOn = 0;
         for (size_t s = 0; s < src.size(); s++) if (src[s].drive && src[s].drive->core && src[s].amps != 0) { on = (int)s; nOn++; }
         const double an = nOn == 1 ? len(src[on].axis) : 0;
         if (an > 0 && src[on].axis.z / an >= 1.0 - 1e-12)
         {
            const V3 o = src[on].origin;
            for (int i = 0; i < N; i++)
            {
               if (!wp[i]) continue;
               double c = 0, s = 0;
               for (int j = 0; j < N; j++)
               {
                  if (!wp[j] || (j != i && (pidx[(size_t)i * N + j] < 0 || !P[pidx[(size_t)i * N + j]].touch))) continue;
                  const double x0 = wp[j]->pos.x - o.x, y0 = wp[j]->pos.y - o.y, x1 = pos[j].x - o.x, y1 = pos[j].y - o.y;
                  c += x0 * x1 + y0 * y1; s += x0 * y1 - y0 * x1;
               }
               const double r = std::hypot(c, s);
               if (r <= 0) continue;
               const cd e(c / r, -s / r); cd pw(1, 0); cd *pi = &ph[(size_t)i * (2 * LP + 1) + LP];
               for (int m = 1; m <= LP; m++) { pw *= e; pi[m] = pw; pi[-m] = std::conj(pw); }
               if (st) st->turned = true;
            }
         }
      }
      std::vector<cd> a4((size_t)N * K4, 0.0), t1(KP), t2(KP), t3(KP);
      for (int i = 0; i < N; i++)
      {
         if (wp[i]) for (int k = 0; k < K4; k++) a4[(size_t)i * K4 + k] = wp[i]->a4[k] * ph[(size_t)i * (2 * LP + 1) + LP + I.l4.m_[k]];
         else for (int k = 0; k < K4; k++) a4[(size_t)i * K4 + k] = -beta(I.l4.n_[k]) * b8[(size_t)i * KP + k];
      }
      // the force a round could still move: mu0 |a|max |da|max, 2n+1 weighted norms of each ball's coefficients (the force between two balls' moments is mu0 a a)
      const auto moved = [&](double a2, double d2) { return MU0 * std::sqrt(a2 * d2); };
      int sweeps = 0; double change = 0;
      while (sweeps < std::max(1, set.gsRounds))
      {
         double a2 = 0, d2 = 0;
         for (int i = 0; i < N; i++)
         {
            for (int k = 0; k < K4; k++) t3[k] = b8[(size_t)i * KP + k];
            for (int j = 0; j < N; j++) { const int q = pidx[(size_t)i * N + j]; if (j == i || q < 0) continue; I.l4.translate(P[q].f4, i < j ? +1 : -1, &a4[(size_t)j * K4], t3.data(), t1.data(), t2.data()); }
            double ai = 0, di = 0;
            for (int k = 0; k < K4; k++) { const cd an = -beta(I.l4.n_[k]) * t3[k]; di += I.l4.W_[k] * std::norm(an - a4[(size_t)i * K4 + k]); ai += I.l4.W_[k] * std::norm(an); a4[(size_t)i * K4 + k] = an; }
            a2 = std::max(a2, ai); d2 = std::max(d2, di);
         }
         sweeps++; change = moved(a2, d2);
         if (sweeps >= (cold ? set.gsCold : 1) && change <= 0.5 * set.stopForce) break;
      }
      if (st) { st->sweeps = sweeps; st->change = change; }
      // ---- coefficients at order 8: global + each ball's own response above order 4 (the plug), then feedback passes on the touching pairs (stop as the sweeps)
      std::vector<cd> ag((size_t)N * KP, 0.0);
      for (int i = 0; i < N; i++) { for (int k = 0; k < K4; k++) ag[(size_t)i * KP + k] = a4[(size_t)i * K4 + k]; if (high[i]) for (int k = K4; k < KP; k++) ag[(size_t)i * KP + k] = -beta(I.l8.n_[k]) * b8[(size_t)i * KP + k]; }
      std::vector<cd> ac = ag, tot((size_t)N * KP), next, din(2 * KP), dout(2 * KP), s4a(K4), s4b(K4);
      for (int i = 0; i < N; i++) if (wp[i]) for (int k = 0; k < KP; k++) ac[(size_t)i * KP + k] = wp[i]->ac[k] * ph[(size_t)i * (2 * LP + 1) + LP + I.l8.m_[k]];   // passes from the last corrected answer, turned
      bool anyTouch = false; for (auto &q : P) anyTouch |= q.touch;
      if (!anyTouch) ac = ag;
      int passes = 0;
      while (passes < set.passes && anyTouch)
      {
         for (int i = 0; i < N; i++)
         {
            cd *ti = &tot[(size_t)i * KP]; for (int k = 0; k < KP; k++) ti[k] = 0;
            for (int j = 0; j < N; j++)
            {
               const int q = pidx[(size_t)i * N + j]; if (j == i || q < 0) continue; const int dir = i < j ? +1 : -1;
               if (P[q].touch)   // the partner's field is kept for the pair step below (it is taken out of the field the pair receives)
               {
                  cd *pf = &P[q].pf[(i == P[q].i ? 1 : 0) * KP]; for (int k = 0; k < KP; k++) pf[k] = 0;
                  I.l8.translate(P[q].f8, dir, &ac[(size_t)j * KP], pf, t1.data(), t2.data()); for (int k = 0; k < KP; k++) ti[k] += pf[k];
               }
               else { for (int k = 0; k < K4; k++) { s4a[k] = ac[(size_t)j * KP + k]; s4b[k] = 0; } I.l4.translate(P[q].f4, dir, s4a.data(), s4b.data(), t1.data(), t2.data()); for (int k = 0; k < K4; k++) ti[k] += s4b[k]; }
            }
         }
         next = ag;
         for (auto &q : P)
         {
            if (!q.touch) continue;
            // frame of D = c_i - c_j: ball 2 (+z) = i, ball 1 (origin) = j
            for (int sball = 0; sball < 2; sball++)
            {
               const int me = sball ? q.i : q.j, ot = sball ? q.j : q.i; cd *dv = &din[sball * KP];
               const cd *pf = &q.pf[sball * KP];   // the partner's field at me, from the field step (build step B6: not translated again)
               for (int k = 0; k < KP; k++) t3[k] = b8[(size_t)me * KP + k] + tot[(size_t)me * KP + k] - pf[k];
               (void)ot;
               I.l8.rot(q.f8, t3.data(), dv, false);
            }
            for (int k = 0; k < 2 * KP; k++) q.x[k] = din[k];
            for (int k = 0; k < 2 * KP; k++) dout[k] = 0;
            for (int m = -LP; m <= LP; m++)
            {
               const int mi = m + LP, nP = ops.sP[mi], S = ops.s[mi]; const auto &rows = I.l8.rowsByM_[mi]; const double *Op = &q.op[ops.oOp[mi]];
               cd x[32]; for (int r = 0; r < nP; r++) { x[r] = din[rows[r]]; x[nP + r] = din[KP + rows[r]]; }
               for (int r = 0; r < S; r++) { cd acc = 0; for (int c = 0; c < S; c++) acc += Op[(size_t)r * S + c] * x[c]; if (r < nP) dout[rows[r]] = acc; else dout[KP + rows[r - nP]] = acc; }
            }
            I.l8.rot(q.f8, dout.data(), t1.data(), true); for (int k = 0; k < KP; k++) next[(size_t)q.j * KP + k] += t1[k];
            I.l8.rot(q.f8, dout.data() + KP, t1.data(), true); for (int k = 0; k < KP; k++) next[(size_t)q.i * KP + k] += t1[k];
         }
         double a2 = 0, d2 = 0;
         for (int i = 0; i < N; i++)
         {
            double ai = 0, di = 0;
            for (int k = 0; k < KP; k++) { const size_t e = (size_t)i * KP + k; di += I.l8.W_[k] * std::norm(next[e] - ac[e]); ai += I.l8.W_[k] * std::norm(next[e]); }
            a2 = std::max(a2, ai); d2 = std::max(d2, di);
         }
         ac = next; passes++; change = moved(a2, d2);
         if (passes >= (cold ? set.passesCold : 1) && change <= 0.5 * set.stopForce) break;
      }
      if (st) { st->passes = passes; st->change = change; }
      I.lastPos = pos; I.lastA = ac; if (ids) I.lastIds = *ids; else I.lastIds.clear();   // kept for FieldFromBalls (gate (b))
      {
         std::unordered_map<uint64_t, Impl::Warm> w; w.reserve(N);   // the warm start of the next Evaluate: this answer, for the balls in it only
         for (int i = 0; i < N; i++) { Impl::Warm &e = w[idv[i]]; e.pos = pos[i]; e.a4.assign(a4.begin() + (size_t)i * K4, a4.begin() + (size_t)(i + 1) * K4); e.ac.assign(ac.begin() + (size_t)i * KP, ac.begin() + (size_t)(i + 1) * KP); }
         I.warm.swap(w);
      }
      // ---- linear forces: drive + every pair (touching: the order-24 pair from its last drive; others: order 4)
      std::vector<V3> Flin(N);
      for (int i = 0; i < N; i++) extF(i, &ac[(size_t)i * KP], Flin[i]);
      const std::vector<V3> Fd = Flin;           // each ball's force from the drive (its magnetisation in the magnets' field)
      std::vector<V3> fpair(P.size());           // each pair's force on ball i (minus on j)
      for (auto &q : P)
      {
         double f3[3];
         if (q.touch)
         {
            cd sz = 0, sp = 0, sm = 0;
            auto xb = [&](int m, cd *x) { const int mi = m + LP, nP = ops.sP[mi]; const auto &rows = I.l8.rowsByM_[mi]; for (int r = 0; r < nP; r++) { x[r] = q.x[rows[r]]; x[nP + r] = q.x[KP + rows[r]]; } };
            cd xm[32], xo[32];
            for (int m = -LP; m <= LP; m++)
            {
               const int mi = m + LP, S = ops.s[mi]; xb(m, xm);
               auto quad = [&](const double *Q, const cd *xr, int S2) { cd acc = 0; for (int r = 0; r < S; r++) { cd v = 0; for (int c = 0; c < S2; c++) v += Q[(size_t)r * S2 + c] * xr[c]; acc += std::conj(xm[r]) * v; } return acc; };
               sz += quad(&q.op[ops.oQz[mi]], xm, S);
               if (m - 1 >= -LP) { xb(m - 1, xo); sp += quad(&q.op[ops.oQp[mi]], xo, ops.s[mi - 1]); }
               if (m + 1 <= LP) { xb(m + 1, xo); sm += quad(&q.op[ops.oQm[mi]], xo, ops.s[mi + 1]); }
            }
            const cd sx = 0.5 * (sp + sm), sy = (sp - sm) * cd(0, -0.5); const double fx = sx.real() / R, fy = sy.real() / R, fz = sz.real() / R; const double *M = q.f8.Rm;
            f3[0] = M[0] * fx + M[1] * fy + M[2] * fz; f3[1] = M[3] * fx + M[4] * fy + M[5] * fz; f3[2] = M[6] * fx + M[7] * fy + M[8] * fz;
            // contact key: axial field at the midpoint
            const int nP0 = ops.sP[LP]; const auto &rows0 = I.l8.rowsByM_[LP]; const double *kz = &q.op[ops.oKz]; double hz = 0;
            for (int r = 0; r < nP0; r++) hz += kz[r] * q.x[rows0[r]].real() + kz[nP0 + r] * q.x[KP + rows0[r]].real();
            q.key = MU0 * std::fabs(hz);
         }
         else
         {
            for (int k = 0; k < K4; k++) { s4a[k] = ac[(size_t)q.i * KP + k]; s4b[k] = ac[(size_t)q.j * KP + k]; }
            std::vector<cd> ri(K4), rl(K4); I.l4.rot(q.f4, s4a.data(), ri.data(), false); I.l4.rot(q.f4, s4b.data(), rl.data(), false); I.l4.sandwich(q.f4, +1, ri.data(), rl.data(), f3);
         }
         const V3 f = { -MU0 * R * f3[0], -MU0 * R * f3[1], -MU0 * R * f3[2] };
         Flin[q.i].x += f.x; Flin[q.i].y += f.y; Flin[q.i].z += f.z; Flin[q.j].x -= f.x; Flin[q.j].y -= f.y; Flin[q.j].z -= f.z;
         fpair[&q - &P[0]] = f;
      }
      // the source a ball's magnet parts turn round: its plug's magnet, else the nearest energised one
      std::vector<int> srcOf(N, 0);
      for (int i = 0; i < N; i++)
      {
         if (plugOf[i] >= 0) { srcOf[i] = plugOf[i]; continue; }
         double best = 1e30; for (size_t s = 0; s < src.size(); s++) { if (!src[s].drive || src[s].amps == 0) continue; const double dd = len(sub(pos[i], src[s].origin)); if (dd < best) { best = dd; srcOf[i] = (int)s; } }
      }
      // dF[b] = sum_k w_bk Flin_k - f_b Flone_b with Flin_k = Fd_k + sum_q (+-) fpair_q: split into the parts (w = identity without saturation)
      auto emit = [&](const std::vector<double> &w, const std::vector<double> &f)
      {
         for (int b = 0; b < N; b++)
         {
            V3 mb;   // the magnet part of ball b itself (its own drive term and its lone pull taken out) is merged into one part
            for (int k = 0; k < N; k++)
            {
               const double c = w[(size_t)b * N + k]; if (c == 0) continue;
               if (k == b) { mb = { c * Fd[k].x - f[b] * Flone[b].x, c * Fd[k].y - f[b] * Flone[b].y, c * Fd[k].z - f[b] * Flone[b].z }; continue; }
               Part pt; pt.ball = b; pt.a = k; pt.b = -1; pt.src = srcOf[k]; pt.F = { c * Fd[k].x, c * Fd[k].y, c * Fd[k].z }; parts->push_back(pt);
            }
            if (w[(size_t)b * N + b] == 0) mb = { -f[b] * Flone[b].x, -f[b] * Flone[b].y, -f[b] * Flone[b].z };
            { Part pt; pt.ball = b; pt.a = b; pt.b = -1; pt.src = srcOf[b]; pt.F = mb; parts->push_back(pt); }
            for (size_t q = 0; q < P.size(); q++)
            {
               const double c = w[(size_t)b * N + P[q].i] - w[(size_t)b * N + P[q].j]; if (c == 0) continue;
               Part pt; pt.ball = b; pt.a = P[q].i; pt.b = P[q].j; pt.src = srcOf[b]; pt.F = { c * fpair[q].x, c * fpair[q].y, c * fpair[q].z }; parts->push_back(pt);
            }
         }
      };
      if (!set.saturation)   // linear balls (offline tests): the change of each ball's force caused by the others
      {
         for (int i = 0; i < N; i++) dF[i] = set.linearTotal ? Flin[i] : V3 { Flin[i].x - Flone[i].x, Flin[i].y - Flone[i].y, Flin[i].z - Flone[i].z };
         if (parts && !set.linearTotal) { std::vector<double> w((size_t)N * N, 0.0), f(N, 1.0); for (int i = 0; i < N; i++) w[(size_t)i * N + i] = 1; emit(w, f); }
         return;
      }
      // ---- saturation, cut rule: clusters of pair-fixed balls (factor 1 beyond 3 mm by the FEM table), tree from the plug (or the most strongly pulled ball); each contact's factor scales the total beyond it
      std::vector<int> comp(N, -1); int nc = 0;
      for (int i = 0; i < N; i++) { if (comp[i] >= 0) continue; std::vector<int> stk = { i }; comp[i] = nc;
         while (!stk.empty()) { const int b = stk.back(); stk.pop_back(); for (int j = 0; j < N; j++) { const int q = pidx[(size_t)b * N + j]; if (q >= 0 && P[q].touch && comp[j] < 0) { comp[j] = nc; stk.push_back(j); } } } nc++; }
      std::vector<V3> R1(N);
      std::vector<int> parAll(N, -1); std::vector<double> facAll(N, 1.0);   // the tree of every cluster, for the parts
      for (int c = 0; c < nc; c++)
      {
         int root = -1; double best = -1;
         for (int i = 0; i < N; i++) if (comp[i] == c && plugOf[i] >= 0 && len(Flone[i]) > best) { best = len(Flone[i]); root = i; }
         if (root < 0) for (int i = 0; i < N; i++) if (comp[i] == c && len(Flone[i]) > best) { best = len(Flone[i]); root = i; }
         std::vector<int> par(N, -2), order = { root }; std::vector<double> fac(N, 1.0); par[root] = -1; fac[root] = plugOf[root] >= 0 ? facePhi[root] : 1.0;
         for (size_t o = 0; o < order.size(); o++)
         {
            const int b = order[o];
            for (int j = 0; j < N; j++) { const int q = pidx[(size_t)b * N + j]; if (q < 0 || !P[q].touch || par[j] != -2) continue; par[j] = b; fac[j] = PhiBB(std::max(P[q].gap, GMIN), P[q].key); order.push_back(j); }
         }
         std::vector<V3> tot(N);
         for (size_t o = order.size(); o-- > 0;) { const int b = order[o]; tot[b].x += Flin[b].x; tot[b].y += Flin[b].y; tot[b].z += Flin[b].z; if (par[b] >= 0) { tot[par[b]].x += tot[b].x; tot[par[b]].y += tot[b].y; tot[par[b]].z += tot[b].z; } }
         // tot now holds subtree sums (linear); scaled totals and R1 = scaled - sum of children's scaled
         std::vector<V3> sc(N); for (int b : order) sc[b] = { fac[b] * tot[b].x, fac[b] * tot[b].y, fac[b] * tot[b].z };
         for (int b : order) { R1[b] = sc[b]; }
         for (int b : order) if (par[b] >= 0) { R1[par[b]].x -= sc[b].x; R1[par[b]].y -= sc[b].y; R1[par[b]].z -= sc[b].z; }
         for (int b : order) { parAll[b] = par[b]; facAll[b] = fac[b]; }
      }
      for (int i = 0; i < N; i++)
      {
         const double f = plugOf[i] >= 0 ? facePhi[i] : 1.0;
         dF[i] = { R1[i].x - f * Flone[i].x, R1[i].y - f * Flone[i].y, R1[i].z - f * Flone[i].z };
      }
      if (parts)
      {
         // R1[b] = fac_b tot_b - sum over children c of fac_c tot_c, tot = subtree sums of Flin -> w_bb = fac_b; for k below b: fac_b - fac_(b's child toward k)
         std::vector<double> w((size_t)N * N, 0.0), f(N, 1.0);
         for (int k = 0; k < N; k++)
         {
            w[(size_t)k * N + k] = facAll[k];
            for (int prev = k, x = parAll[k]; x >= 0; prev = x, x = parAll[x]) w[(size_t)x * N + k] = facAll[x] - facAll[prev];
            f[k] = plugOf[k] >= 0 ? facePhi[k] : 1.0;
         }
         emit(w, f);
      }
   }
}
