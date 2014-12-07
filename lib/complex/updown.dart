/// CSparse: a Concise Sparse matrix package.
/// Copyright (c) 2006, Timothy A. Davis.
/// http://www.cise.ufl.edu/research/sparse/CSparse
///
/// CSparse is free software; you can redistribute it and/or
/// modify it under the terms of the GNU Lesser General Public
/// License as published by the Free Software Foundation; either
/// version 2.1 of the License, or (at your option) any later version.
///
/// CSparse is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
/// Lesser General Public License for more details.
///
/// You should have received a copy of the GNU Lesser General Public
/// License along with this Module; if not, write to the Free Software
/// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
part of edu.emory.mathcs.cxsparse;

/// Sparse Cholesky rank-1 update/downdate, L*L' + sigma*w*w' (sigma = +1 or -1).
/// Note that this CXSparse version is different than CSparse.
///
/// [L] factorization to update/downdate.
/// [sigma] +1 for update, -1 for downdate.
/// [parent] the elimination tree of L.
/// Returns true if successful, false on error.
bool cs_updown(DZcs L, int sigma, DZcs C, Int32List parent) {
  int n, p, f, j;
  Int32List Lp, Li, Cp, Ci;
  DZcsa Lx = new DZcsa(),
      Cx = new DZcsa(),
      w;
  Float64List alpha, gamma, w1, w2;
  double phase,
      beta = 1.0,
      delta,
      beta2 = 1.0;
  if (!CS_CSC(L) || !CS_CSC(C) || parent == null) {
    return false;
  }
  Lp = L.p;
  Li = L.i;
  Lx.x = L.x;
  n = L.n;
  Cp = C.p;
  Ci = C.i;
  Cx.x = C.x;
  if ((p = Cp[0]) >= Cp[1]) return true; // return if C empty
  w = new DZcsa.sized(n); // get workspace
  f = Ci[p];
  for ( ; p < Cp[1]; p++) {
    f = math.min(f, Ci[p]); // f = min (find (C))
  }
  for (j = f; j != -1; j = parent[j]) {
    w.set_list(j, cs_czero()); // clear workspace w
  }
  for (p = Cp[0]; p < Cp[1]; p++) {
    w.set_list(Ci[p], Cx.get(p)); // w = C
  }
  for (j = f; j != -1; j = parent[j]) // walk path f up to root
  {
    p = Lp[j];
    alpha = cs_cdiv_list(w.get(j), Lx.get(p)); // alpha = w(j) / L(j,j)
    /* CXSparse */
    beta2 = beta * beta + sigma * cs_creal(cs_cmult_list(alpha, cs_conj(alpha)));
    if (beta2 <= 0) break; // not positive definite
    beta2 = math.sqrt(beta2);
    delta = (sigma > 0) ? (beta / beta2) : (beta2 / beta);
    gamma = cs_cmult(cs_cdiv(cs_conj(alpha), beta2 * beta, 0.0), sigma.toDouble());
    Lx.set_list(p, cs_cplus(cs_cmult(Lx.get(p), delta), (sigma > 0) ? cs_cmult_list(gamma, w.get(j)) : cs_czero()));
    beta = beta2;
    /* CXSparse */
    phase = cs_cabs_list(cs_cdiv_list(Lx.get(p), Lx.get(p))); // phase = abs(L(j,j))/L(j,j)
    Lx.set_list(p, cs_cmult(Lx.get(p), phase)); // L(j,j) = L(j,j) * phase
    for (p++; p < Lp[j + 1]; p++) {
      w1 = w.get(Li[p]);
      w.set_list(Li[p], cs_cminus(w1, cs_cmult_list(alpha, Lx.get(p))));
      w2 = w.get(Li[p]);
      Lx.set_list(p, cs_cplus(cs_cmult(Lx.get(p), delta), cs_cmult_list(gamma, sigma > 0 ? w1 : w2)));
      /* CXSparse */
      Lx.set_list(p, cs_cmult(Lx.get(p), phase)); // L(i,j) = L(i,j) * phase
    }
  }
  w = null;
  return beta2 > 0;
}
