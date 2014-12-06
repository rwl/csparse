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

/// Sparse Cholesky.
///
/// Numeric Cholesky factorization LL=PAP'.
///
/// Only the upper triangular part of [A] is used.
/// [S] symbolic Cholesky analysis, pinv is optional.
/// Returns numeric Cholesky factorization, null on error.
Dcsn cs_chol(Dcs A, Dcss S) {
  double d, lki;
  Float64List Lx, x, Cx;
  int n;
  Int32List Li, Lp, cp, pinv, s, c, parent, Cp, Ci;
  Dcs L, C;
  Dcsn N;
  if (!cs_csc(A) || S == null || S.cp == null || S.parent == null) {
    return (null);
  }
  n = A.n;
  N = new Dcsn(); // allocate result
  c = new Int32List(2 * n); // get int workspace
  x = new Float64List(n); // get double workspace
  cp = S.cp;
  pinv = S.pinv;
  parent = S.parent;
  C = pinv != null ? cs_symperm(A, pinv, true) : A;
  s = c;
  int s_offset = n;
  Cp = C.p;
  Ci = C.i;
  Cx = C.x;
  N.L = L = cs_spalloc(n, n, cp[n], true, false); // allocate result
  Lp = L.p;
  Li = L.i;
  Lx = L.x;
  for (int k = 0; k < n; k++) {
    Lp[k] = c[k] = cp[k];
  }
  for (int k = 0; k < n; k++) // compute L(k,:) for L*L' = C
  {
    /* Nonzero pattern of L(k,:) */
    int top = cs_ereach(C, k, parent, s, s_offset, c); // find pattern of L(k,:)
    x[k] = 0.0; // x (0:k) is now zero
    for (int p = Cp[k]; p < Cp[k + 1]; p++) // x = full(triu(C(:,k)))
    {
      if (Ci[p] <= k) {
        x[Ci[p]] = Cx[p];
      }
    }
    d = x[k]; // d = C(k,k)
    x[k] = 0.0; // clear x for k+1st iteration
    /* Triangular solve */
    for ( ; top < n; top++) // solve L(0:k-1,0:k-1) * x = C(:,k)
    {
      final i = s[s_offset + top]; // s [top..n-1] is pattern of L(k,:)
      lki = x[i] / Lx[Lp[i]]; // L(k,i) = x (i) / L(i,i)
      x[i] = 0.0; // clear x for k+1st iteration
      for (int p = Lp[i] + 1; p < c[i]; p++) {
        x[Li[p]] -= Lx[p] * lki;
      }
      d -= lki * lki; // d = d - L(k,i)*L(k,i)
      final p = c[i]++;
      Li[p] = k; // store L(k,i) in column i
      Lx[p] = lki;
    }
    /* Compute L(k,k) */
    if (d <= 0) {
      return null; // not pos def
    }
    final p = c[k]++;
    Li[p] = k; // store L(k,k) = sqrt (d) in column k
    Lx[p] = math.sqrt(d);
  }
  Lp[n] = cp[n]; // finalize L
  return N;
}
