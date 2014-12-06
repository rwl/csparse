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

/// Solve a least-squares problem (min ||Ax-b||_2, where A is m-by-n with
/// m >= n) or underdetermined system (Ax=b, where m < n).
///
/// [order] ordering method to use (0 to 3).
/// [b] size max(m,n), b (size m) on input, x(size n) on output.
/// Returns true if successful, false on error.
bool qrsol(int order, Matrix A, Float64List b) {
  Float64List x;
  Symbolic S;
  Numeric N;
  Matrix AT = null;
  int k, m, n;
  bool ok;
  if (!csc(A) || b == null) {
    return false; // check inputs
  }
  n = A.n;
  m = A.m;
  if (m >= n) {
    S = sqr(order, A, true); // ordering and symbolic analysis
    N = qr(A, S); // numeric QR factorization
    x = new Float64List(S != null ? S.m2 : 1); // get workspace
    ok = (S != null && N != null);
    if (ok) {
      ipvec(S.pinv, b, x, m); // x(0:m-1) = b(p(0:m-1)
      for (k = 0; k < n; k++) // apply Householder refl. to x
      {
        happly(N.L, k, N.B[k], x);
      }
      usolve(N.U, x); // x = R\x
      ipvec(S.q, x, b, n); // b(q(0:n-1)) = x(0:n-1)
    }
  } else {
    AT = transpose(A, true); // Ax=b is underdetermined
    S = sqr(order, AT, true); // ordering and symbolic analysis
    N = qr(AT, S); // numeric QR factorization of A'
    x = new Float64List(S != null ? S.m2 : 1); // get workspace
    ok = (AT != null && S != null && N != null);
    if (ok) {
      pvec(S.q, b, x, m); // x(q(0:m-1)) = b(0:m-1)
      utsolve(N.U, x); // x = R'\x
      for (k = m - 1; k >= 0; k--) // apply Householder refl. to x
      {
        happly(N.L, k, N.B[k], x);
      }
      pvec(S.pinv, x, b, n); // b(0:n-1) = x(p(0:n-1))
    }
  }
  return ok;
}
