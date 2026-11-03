// license:GPLv3+
// The magnets' ball-pair force as a stored energy of place (2026-11-03; the form of runs/2026-11-02-cluster-cost/notes.md section 6). No VPX dependencies.
// At each ball-ball solve: U(x) = U(xc) - sum_k g_k (q_k(x) - q_k(xc)), q_k = invariants of the balls' places: each pair's distance, each ball's half squared offset from
// a magnet's axis and its height along that axis (per energised magnet); g = the least-norm fit of the solve's dF (N). The force is -grad U = sum g_k grad q_k.
// U depends on no angle about an axis, so a turn of everything about one magnet does no work; a net torque in dF is left out of the fit (the residual).
#pragma once
#include <cmath>
#include <vector>

namespace PairForm
{
   // kind 0: |c_a - c_b| (m; b < 0: to the point o); 1: half the squared offset of c_a from the axis through o along n (m2); 2: the height of c_a along n above o (m)
   struct Term { int kind = 0, a = 0, b = -1; double g = 0, q0 = 0, o[3] = { 0, 0, 0 }, n[3] = { 0, 0, 1 }; };

   // q of a term at centres ca, cb (cb: kind 0 only); ga, gb: its slopes in ca, cb when wanted
   inline double Q(const Term &t, const double *ca, const double *cb, double *ga = nullptr, double *gb = nullptr)
   {
      if (t.kind == 0)
      {
         if (cb == nullptr) cb = t.o;   // (b < 0: the other ball held at o, outside the group)
         const double d[3] = { ca[0] - cb[0], ca[1] - cb[1], ca[2] - cb[2] }, l = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
         for (int c = 0; c < 3; c++)
         {
            const double u = l > 0 ? d[c] / l : 0.0;
            if (ga) ga[c] = u;
            if (gb) gb[c] = -u;
         }
         return l;
      }
      const double d[3] = { ca[0] - t.o[0], ca[1] - t.o[1], ca[2] - t.o[2] }, h = d[0] * t.n[0] + d[1] * t.n[1] + d[2] * t.n[2];
      if (t.kind == 2)
      {
         if (ga) for (int c = 0; c < 3; c++) ga[c] = t.n[c];
         return h;
      }
      const double r[3] = { d[0] - h * t.n[0], d[1] - h * t.n[1], d[2] - h * t.n[2] };
      if (ga) for (int c = 0; c < 3; c++) ga[c] = r[c];
      return 0.5 * (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
   }

   // The terms for N balls at pos (m, 3 per ball) and the axes (o, n: 3 per magnet, n unit), fitted to dF (N, 3 per ball): least squares, least norm over unit-scaled
   // columns (a tiny ridge, 1e-12 of the trace, picks the least norm among fits that give the same force). Returns the fit's residual |G g - dF| (N, all balls).
   inline double Fit(int N, const double *pos, int nMag, const double *o, const double *n, const double *dF, std::vector<Term> &terms)
   {
      terms.clear();
      for (int i = 0; i < N; i++)
         for (int j = i + 1; j < N; j++) { Term t; t.kind = 0; t.a = i; t.b = j; terms.push_back(t); }
      for (int s = 0; s < nMag; s++)
         for (int i = 0; i < N; i++)
            for (int k = 1; k <= 2; k++)
            {
               Term t; t.kind = k; t.a = i;
               for (int c = 0; c < 3; c++) { t.o[c] = o[3 * s + c]; t.n[c] = n[3 * s + c]; }
               terms.push_back(t);
            }
      const int K = (int)terms.size(), R = 3 * N;
      std::vector<double> G((size_t)R * K, 0.0), sc(K, 1.0);
      for (int k = 0; k < K; k++)
      {
         Term &t = terms[k]; double ga[3], gb[3];
         t.q0 = Q(t, &pos[3 * t.a], t.b >= 0 ? &pos[3 * t.b] : nullptr, ga, gb);
         double s2 = 0;
         for (int c = 0; c < 3; c++) { G[(size_t)(3 * t.a + c) * K + k] = ga[c]; s2 += ga[c] * ga[c]; }
         if (t.b >= 0) for (int c = 0; c < 3; c++) { G[(size_t)(3 * t.b + c) * K + k] = gb[c]; s2 += gb[c] * gb[c]; }
         sc[k] = s2 > 0 ? 1.0 / std::sqrt(s2) : 0.0;
      }
      // (A^T A + lam I) y = A^T dF with A = G diag(sc); g = sc y
      std::vector<double> M((size_t)K * K, 0.0), y(K, 0.0);
      for (int a = 0; a < K; a++)
      {
         for (int b = a; b < K; b++) { double s = 0; for (int r = 0; r < R; r++) s += G[(size_t)r * K + a] * G[(size_t)r * K + b]; M[(size_t)a * K + b] = M[(size_t)b * K + a] = s * sc[a] * sc[b]; }
         double s = 0; for (int r = 0; r < R; r++) s += G[(size_t)r * K + a] * dF[r]; y[a] = s * sc[a];
      }
      double tr = 0; for (int a = 0; a < K; a++) tr += M[(size_t)a * K + a];
      for (int a = 0; a < K; a++) M[(size_t)a * K + a] += 1e-12 * tr + 1e-300;
      for (int a = 0; a < K; a++)   // Cholesky, in place (lower)
      {
         for (int b = 0; b <= a; b++)
         {
            double s = M[(size_t)a * K + b]; for (int c = 0; c < b; c++) s -= M[(size_t)a * K + c] * M[(size_t)b * K + c];
            M[(size_t)a * K + b] = a == b ? std::sqrt(std::fmax(s, 1e-300)) : s / M[(size_t)b * K + b];
         }
      }
      for (int a = 0; a < K; a++) { double s = y[a]; for (int c = 0; c < a; c++) s -= M[(size_t)a * K + c] * y[c]; y[a] = s / M[(size_t)a * K + a]; }
      for (int a = K - 1; a >= 0; a--) { double s = y[a]; for (int c = a + 1; c < K; c++) s -= M[(size_t)c * K + a] * y[c]; y[a] = s / M[(size_t)a * K + a]; }
      for (int k = 0; k < K; k++) terms[k].g = sc[k] * y[k];
      double res = 0;
      for (int r = 0; r < R; r++) { double s = -dF[r]; for (int k = 0; k < K; k++) s += G[(size_t)r * K + k] * terms[k].g; res += s * s; }
      return std::sqrt(res);
   }
}
