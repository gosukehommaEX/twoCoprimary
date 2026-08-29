#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Conditional probability g(y2 | y1, N, p1, p2, gamma) of the bivariate
// binomial distribution, equation (3) of Homma and Yoshida (2025).
//
// For each pair the sum runs over the set M = {m: max(0, y2 - (N - y1)), ...,
// min(y1, y2)}. The R implementation padded that ragged set into a rectangular
// matrix, which allocates on the order of N ^ 3 doubles for a whole outcome
// grid.
//
// Every exponent and every binomial coefficient argument in the summand lies in
// 0, ..., N, so the powers and the coefficients are tabulated once and the
// inner loop reduces to six multiplications and an addition. The tables are
// filled with the same routines that a direct evaluation would call, and the
// factors are multiplied in the same order, so the result is bit for bit what
// the direct expression gives while the number of calls to pow and choose falls
// from order N ^ 3 to order N ^ 2.
//
// The accumulator is a long double, matching the extended precision that
// rowSums uses, so the result also agrees with the R reference implementation
// retained in .dbibinom_g_r().
//
// N      number of subjects in the group
// y1, y2 responder counts of the two endpoints, of equal length
// xi     xi = p2 + gamma (p2 - p1)
// gamma  dependence parameter obtained from the correlation
// [[Rcpp::export]]
NumericVector dbibinom_g(int N, IntegerVector y1, IntegerVector y2,
                         double xi, double gamma) {

  const R_xlen_t n = y1.size();
  NumericVector out(n);
  if (n == 0) return out;

  const double base = std::pow(1.0 + gamma, -static_cast<double>(N));
  const double a_xi = xi + gamma;
  const double b_xi = 1.0 - xi;
  const double c_xi = 1.0 - xi + gamma;

  const int M = N + 1;

  // Powers of each base for every exponent that can occur
  std::vector<double> pa(M), pb(M), px(M), pc(M);
  for (int k = 0; k < M; ++k) {
    const double e = static_cast<double>(k);
    pa[k] = std::pow(a_xi, e);
    pb[k] = std::pow(b_xi, e);
    px[k] = std::pow(xi, e);
    pc[k] = std::pow(c_xi, e);
  }

  // Binomial coefficients for every argument pair that can occur
  std::vector<double> cf(static_cast<size_t>(M) * static_cast<size_t>(M));
  for (int a = 0; a < M; ++a) {
    double* row = &cf[static_cast<size_t>(a) * static_cast<size_t>(M)];
    for (int b = 0; b < M; ++b) {
      row[b] = R::choose(a, b);
    }
  }

  for (R_xlen_t t = 0; t < n; ++t) {

    const int u = y1[t];
    const int v = y2[t];
    const int lo = std::max(0, v - (N - u));
    const int hi = std::min(u, v);

    const double* cu = &cf[static_cast<size_t>(u) * static_cast<size_t>(M)];
    const double* cw = &cf[static_cast<size_t>(N - u) * static_cast<size_t>(M)];

    long double s = 0.0L;
    for (int m = lo; m <= hi; ++m) {
      const int w = v - m;
      s += static_cast<long double>(
        cu[m] * cw[w] * pa[m] * pb[u - m] * px[w] * pc[N - u - w]
      );
    }

    out[t] = base * static_cast<double>(s);
  }

  return out;
}
