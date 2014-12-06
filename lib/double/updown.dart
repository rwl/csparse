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
part of edu.emory.mathcs.csparse;

/// Sparse Cholesky rank-1 update/downdate, L*L' + sigma*w*w'
/// (sigma = +1 or -1)
///
/// [L] factorization to update/downdate.
/// [sigma] +1 for update, -1 for downdate.
/// [parent] the elimination tree of L.
/// Returns true if successful, false on error.
bool cs_updown(Dcs L, int sigma, Dcs C, Int32List parent) {
  int n, p, f;
  Int32List Lp, Li, Cp, Ci;
  Float64List Lx, Cx, w;
  double alpha, delta, gamma, w1, w2;
  double beta = 1.0,
      beta2 = 1.0;
  if (!cs_csc(L) || !cs_csc(C) || parent == null) {
    return (false); // check inputs
  }
  Lp = L.p;
  Li = L.i;
  Lx = L.x;
  n = L.n;
  Cp = C.p;
  Ci = C.i;
  Cx = C.x;
  if ((p = Cp[0]) >= Cp[1]) {
    return true; // return if C empty
  }
  w = new Float64List(n); // get workspace
  f = Ci[p];
  for ( ; p < Cp[1]; p++) {
    f = math.min(f, Ci[p]); // f = min (find (C))
  }
  for (int j = f; j != -1; j = parent[j]) {
    w[j] = 0.0; // clear workspace w
  }
  for (p = Cp[0]; p < Cp[1]; p++) {
    w[Ci[p]] = Cx[p]; // w = C
  }
  for (int j = f; j != -1; j = parent[j]) // walk path f up to root
  {
    p = Lp[j];
    alpha = w[j] / Lx[p]; // alpha = w(j) / L(j,j)
    beta2 = beta * beta + sigma * alpha * alpha;
    if (beta2 <= 0) {
      break; // not positive definite
    }
    beta2 = math.sqrt(beta2);
    delta = (sigma > 0) ? (beta / beta2) : (beta2 / beta);
    gamma = sigma * alpha / (beta2 * beta);
    Lx[p] = delta * Lx[p] + ((sigma > 0) ? (gamma * w[j]) : 0);
    beta = beta2;
    for (p++; p < Lp[j + 1]; p++) {
      w1 = w[Li[p]];
      w[Li[p]] = w2 = w1 - alpha * Lx[p];
      Lx[p] = delta * Lx[p] + gamma * ((sigma > 0) ? w1 : w2);
    }
  }
  return beta2 > 0;
}
