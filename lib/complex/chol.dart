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

/// Numeric Cholesky factorization LL=PAP'.
///
/// [A] column-compressed matrix, only upper triangular part is used.
/// [S] symbolic Cholesky analysis, pinv is optional.
/// Returns numeric Cholesky factorization, null on error.
Numeric chol(Matrix A, Symbolic S) {
  Float64List d, lki;
  Vector Lx = new Vector(),
      x,
      Cx = new Vector();
  int top, i, p, k, n;
  Int32List Li, Lp, cp, pinv, s, c, parent, Cp, Ci;
  Matrix L, C, E;
  Numeric N;
  if (!csc(A) || S == null || S.cp == null || S.parent == null) {
    return (null);
  }
  n = A.n;
  N = new Numeric(); // allocate result
  c = new Int32List(2 * n); // get int workspace
  x = new Vector.sized(n); // get complex workspace
  cp = S.cp;
  pinv = S.pinv;
  parent = S.parent;
  C = pinv != null ? symperm(A, pinv, true) : A;
  E = pinv != null ? C : null; // E is alias for A, or a copy E=A(p,p)
  if (N == null || c == null || x == null || C == null) {
    return _ndone(N, E, c, x, false);
  }
  s = c;
  int s_offset = n;
  Cp = C.p;
  Ci = C.i;
  Cx.x = C.x;
  N.L = L = spalloc(n, n, cp[n], true, false); // allocate result
  if (L == null) {
    return _ndone(N, E, c, x, false);
  }
  Lp = L.p;
  Li = L.i;
  Lx.x = L.x;
  for (k = 0; k < n; k++) {
    Lp[k] = c[k] = cp[k];
  }
  for (k = 0; k < n; k++) // compute L(k,:) for L*L' = C
  {
    /* Nonzero pattern of L(k,:) */
    top = ereach(C, k, parent, s, s_offset, c); // find pattern of L(k,:)
    x.setList(k, czero()); // x (0:k) is now zero
    for (p = Cp[k]; p < Cp[k + 1]; p++) // x = full(triu(C(:,k)))
    {
      if (Ci[p] <= k) x.setList(Ci[p], Cx.get(p));
    }
    d = x.get(k); // d = C(k,k)
    x.setList(k, czero()); // clear x for k+1st iteration
    /* Triangular solve */
    for ( ; top < n; top++) // solve L(0:k-1,0:k-1) * x = C(:,k)
    {
      i = s[s_offset + top]; // s [top..n-1] is pattern of L(k,:)
      lki = cdiv_list(x.get(i), Lx.get(Lp[i])); // L(k,i) = x (i) / L(i,i)
      x.setList(i, czero()); // clear x for k+1st iteration
      for (p = Lp[i] + 1; p < c[i]; p++) {
        x.setList(Li[p], cminus(x.get(Li[p]), cmult_list(Lx.get(p), lki)));
      }
      d = cminus(d, cmult_list(lki, conj(lki))); // d = d - L(k,i)*L(k,i)
      p = c[i]++;
      Li[p] = k; // store L(k,i) in column i
      Lx.setList(p, conj(lki));
    }
    /* Compute L(k,k) */
    if (d[0] <= 0 || d[1] != 0) {
      return _ndone(N, E, c, x, false); // not pos def
    }
    p = c[k]++;
    Li[p] = k; // store L(k,k) = sqrt (d) in column k
    Lx.setList(p, csqrt(d));
  }
  Lp[n] = cp[n]; // finalize L
  return _ndone(N, E, c, x, true); // success: free E,s,x; return N
}
