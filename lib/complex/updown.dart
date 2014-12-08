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
bool updown(Matrix L, int sigma, Matrix C, Int32List parent) {
  int n, p, f, j;
  Int32List Lp, Li, Cp, Ci;
  Vector Lx = new Vector(),
      Cx = new Vector(),
      w;
  Float64List alpha, gamma, w1, w2;
  double phase,
      beta = 1.0,
      delta,
      beta2 = 1.0;
  if (!csc(L) || !csc(C) || parent == null) {
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
  w = new Vector.sized(n); // get workspace
  f = Ci[p];
  for ( ; p < Cp[1]; p++) {
    f = math.min(f, Ci[p]); // f = min (find (C))
  }
  for (j = f; j != -1; j = parent[j]) {
    w.setList(j, czero()); // clear workspace w
  }
  for (p = Cp[0]; p < Cp[1]; p++) {
    w.setList(Ci[p], Cx.get(p)); // w = C
  }
  for (j = f; j != -1; j = parent[j]) // walk path f up to root
  {
    p = Lp[j];
    alpha = cdiv_list(w.get(j), Lx.get(p)); // alpha = w(j) / L(j,j)
    /* CXSparse */
    beta2 = beta * beta + sigma * _creal(cmult_list(alpha, conj(alpha)));
    if (beta2 <= 0) break; // not positive definite
    beta2 = math.sqrt(beta2);
    delta = (sigma > 0) ? (beta / beta2) : (beta2 / beta);
    gamma = cmult(cdiv(conj(alpha), beta2 * beta, 0.0), sigma.toDouble());
    Lx.setList(p, cplus(cmult(Lx.get(p), delta), (sigma > 0) ? cmult_list(gamma, w.get(j)) : czero()));
    beta = beta2;
    /* CXSparse */
    phase = cabs_list(cdiv_list(Lx.get(p), Lx.get(p))); // phase = abs(L(j,j))/L(j,j)
    Lx.setList(p, cmult(Lx.get(p), phase)); // L(j,j) = L(j,j) * phase
    for (p++; p < Lp[j + 1]; p++) {
      w1 = w.get(Li[p]);
      w.setList(Li[p], cminus(w1, cmult_list(alpha, Lx.get(p))));
      w2 = w.get(Li[p]);
      Lx.setList(p, cplus(cmult(Lx.get(p), delta), cmult_list(gamma, sigma > 0 ? w1 : w2)));
      /* CXSparse */
      Lx.setList(p, cmult(Lx.get(p), phase)); // L(i,j) = L(i,j) * phase
    }
  }
  w = null;
  return beta2 > 0;
}
