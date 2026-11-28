// Compiled versions of loops that were slow in R. Each one reproduces the R
// code it replaces: random numbers are drawn through R's own generators
// (R::runif, R::rbeta) in the same order, so seeded results are unchanged.

#include <RcppArmadillo.h>
#include <cmath>
#include <string>
#include <unordered_set>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Pearson correlation between a fixed vector `y` and the masked elements of
// `x` after permuting its rows and columns, for every permutation (one per
// COLUMN of `perms`, 1-based as in R). `rows` and `cols` are the 0-based
// positions of the masked elements, in the order R's x[mask] returns them.
// This is the inner loop of the Mantel tests: cor(x[p, p][mask], y) for each p.
// [[Rcpp::export(name = ".mantel_perm_cor_cpp")]]
NumericVector mantel_perm_cor(const NumericMatrix& x, const NumericVector& y,
                              const IntegerVector& rows, const IntegerVector& cols,
                              const IntegerMatrix& perms) {
  const R_xlen_t n_elem = rows.size();
  const int n = x.nrow();
  if (cols.size() != n_elem || y.size() != n_elem)
    Rcpp::stop("rows, cols and y must have the same length");
  if (perms.nrow() != n)
    Rcpp::stop("each permutation must have one entry per row of x");
  const int n_perm = perms.ncol();
  double y_mean = 0.0;
  for (R_xlen_t k = 0; k < n_elem; ++k) y_mean += y[k];
  y_mean /= n_elem;
  std::vector<double> yc(n_elem), v(n_elem);
  double syy = 0.0;
  for (R_xlen_t k = 0; k < n_elem; ++k) { yc[k] = y[k] - y_mean; syy += yc[k] * yc[k]; }
  const double y_norm = std::sqrt(syy);
  const double* xp = x.begin();
  NumericVector out(n_perm);
  for (int p = 0; p < n_perm; ++p) {
    const int* perm = perms.begin() + (R_xlen_t) p * n;
    double v_mean = 0.0;
    for (R_xlen_t k = 0; k < n_elem; ++k) {
      v[k] = xp[(R_xlen_t) (perm[cols[k]] - 1) * n + (perm[rows[k]] - 1)];
      v_mean += v[k];
    }
    v_mean /= n_elem;
    double sxy = 0.0, sxx = 0.0;
    for (R_xlen_t k = 0; k < n_elem; ++k) {
      const double d = v[k] - v_mean;
      sxy += d * yc[k];
      sxx += d * d;
    }
    out[p] = sxy / (std::sqrt(sxx) * y_norm);
  }
  return out;
}

// Random correlation matrix from uniformly drawn partial correlations: the
// loop of RandCorr(), before its final random permutation (done in R).
// [[Rcpp::export(name = ".rand_corr_partial_cpp")]]
NumericMatrix rand_corr_partial(int n, double ke) {
  NumericMatrix rc(n, n), b(n, n);
  for (int j = 0; j < n; ++j)
    for (int i = j; i < n; ++i) b(i, j) = 1.0;
  for (int i = 1; i < n; ++i)
    rc(i, 0) = -1.0 + 2.0 * std::nearbyint(R::runif(0.0, 1.0) * 1e8) / 1e8;
  for (int i = 1; i < n; ++i) {
    b(i, 0) = rc(i, 0);
    const double s = std::sqrt(1.0 - rc(i, 0) * rc(i, 0));
    for (int j = 1; j <= i; ++j) b(i, j) = s;
  }
  for (int i = 2; i < n; ++i) {
    for (int j = 1; j < i; ++j) {
      double b1 = 0.0;
      for (int l = 0; l < j; ++l) b1 += b(j, l) * b(i, l);
      const double b2 = b(j, j) * b(i, j);
      const double z = b1 + b2, y = b1 - b2;
      if (b2 < ke) rc(i, j) = b1;
      else rc(i, j) = y + (z - y) * std::nearbyint(R::runif(0.0, 1.0) * 1e8) / 1e8;
      const double cosinv = (rc(i, j) - b1) / b2;
      if (std::isfinite(cosinv)) {
        if (cosinv > 1) {
          for (int l = j + 1; l < n; ++l) b(i, l) = 0.0;
        } else if (cosinv < -1) {
          b(i, j) = -b(i, j);
          for (int l = j + 1; l < n; ++l) b(i, l) = 0.0;
        } else {
          b(i, j) = b(i, j) * cosinv;
          const double sin_theta = std::sqrt(1.0 - cosinv * cosinv);
          for (int l = j + 1; l < n; ++l) b(i, l) = b(i, l) * sin_theta;
        }
      }
    }
  }
  NumericMatrix out(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      out(i, j) = (i == j) ? rc(i, j) + rc(j, i) + 1.0 : rc(i, j) + rc(j, i);
  return out;
}

// LKJ(shape) correlation matrix by the vine method: the loop of RandLKJ().
// [[Rcpp::export(name = ".rand_lkj_cpp")]]
NumericMatrix rand_lkj(int k, double shape) {
  NumericMatrix P(k, k), rc(k, k);
  for (int i = 0; i < k; ++i) rc(i, i) = 1.0;
  for (int m = 1; m <= k - 1; ++m) {
    const double beta = shape + (k - 1 - m) / 2.0;
    for (int i = m + 1; i <= k; ++i) {
      P(m - 1, i - 1) = 2 * R::rbeta(beta, beta) - 1;
      double p = P(m - 1, i - 1);
      for (int l = m - 1; l >= 1; --l) {
        const double pli = P(l - 1, i - 1), plm = P(l - 1, m - 1);
        p = p * std::sqrt((1 - pli * pli) * (1 - plm * plm)) + pli * plm;
      }
      rc(m - 1, i - 1) = p;
      rc(i - 1, m - 1) = p;
    }
  }
  return rc;
}

// Combinations of modularity hypotheses for CombineHypot(). Row c of `masks`
// selects the hypothesis columns of combination c. Combining several columns
// joins traits that share any selected module (1 if they do, 0 otherwise); a
// single column x gives outer(x, x). Keeps the first occurrence of each
// distinct lower triangle that is not all zeros, as the R version did by
// comparing pasted lower triangles, and returns the kept matrices with a unit
// diagonal.
// [[Rcpp::export(name = ".combine_hypot_cpp")]]
List combine_hypot(const NumericMatrix& hypot, const LogicalMatrix& masks) {
  const int n = hypot.nrow(), h = hypot.ncol(), n_comb = masks.nrow();
  if (masks.ncol() != h) Rcpp::stop("masks must have one column per hypothesis");
  bool binary = h <= 30;
  for (R_xlen_t i = 0; binary && i < hypot.size(); ++i)
    binary = (hypot[i] == 0.0 || hypot[i] == 1.0);
  // for 0/1 hypotheses: bit k of member[i] says trait i is in module k
  std::vector<unsigned int> member(n, 0u);
  if (binary)
    for (int i = 0; i < n; ++i)
      for (int k = 0; k < h; ++k)
        if (hypot(i, k) == 1.0) member[i] |= (1u << k);
  const R_xlen_t n_low = (R_xlen_t) n * (n - 1) / 2;
  std::vector<double> value(n_low);
  std::vector<int> selected;
  std::string key;
  std::unordered_set<std::string> seen;
  LogicalVector keep(n_comb);
  List kept(n_comb); // filled from the front; List::push_back would copy every time
  R_xlen_t n_kept = 0;
  for (int c = 0; c < n_comb; ++c) {
    selected.clear();
    unsigned int comb = 0u;
    for (int k = 0; k < h; ++k)
      if (masks(c, k)) { selected.push_back(k); comb |= (1u << k); }
    bool all_zero = true;
    R_xlen_t pos = 0;
    for (int j = 0; j < n; ++j) {
      for (int i = j + 1; i < n; ++i, ++pos) { // column-major lower triangle
        double val;
        if (binary) {
          val = (member[i] & member[j] & comb) ? 1.0 : 0.0;
        } else if (selected.size() > 1) {
          double s = 0.0;
          for (int k : selected) s += hypot(i, k) * hypot(j, k);
          val = (s != 0.0) ? 1.0 : 0.0;
        } else {
          val = hypot(i, selected[0]) * hypot(j, selected[0]);
          if (val == 0.0) val = 0.0; // fold -0 into 0, as paste() did
        }
        value[pos] = val;
        if (val != 0.0) all_zero = false;
      }
    }
    if (all_zero) continue;
    if (binary) { // one bit per pair
      key.assign((n_low + 7) / 8, '\0');
      for (R_xlen_t q = 0; q < n_low; ++q)
        if (value[q] != 0.0) key[q / 8] |= (char) (1 << (q % 8));
    } else {
      key.assign(reinterpret_cast<const char*>(value.data()), n_low * sizeof(double));
    }
    if (!seen.insert(key).second) continue;
    keep[c] = true;
    NumericMatrix m(n, n);
    pos = 0;
    for (int j = 0; j < n; ++j) {
      m(j, j) = 1.0;
      for (int i = j + 1; i < n; ++i, ++pos) m(i, j) = m(j, i) = value[pos];
    }
    kept[n_kept++] = m;
  }
  List matrices(n_kept);
  for (R_xlen_t q = 0; q < n_kept; ++q) matrices[q] = kept[q];
  return List::create(Named("keep") = keep, Named("matrices") = matrices);
}
