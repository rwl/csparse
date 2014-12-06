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

/// Solves Ax=b using sparse LU factorization, where A is square and
/// nonsingular. b overwritten with solution. Partial pivoting if tol=1.
///
/// [order] ordering method to use (0 to 3).
/// [b] size n, b on input, x on output.
/// [tol] partial pivoting tolerance.
/// Returns true if successful, false on error.
bool lusol(int order, Matrix A, Float64List b, double tol) {
  Float64List x;
  Symbolic S;
  Numeric N;
  int n;
  bool ok;
  if (!csc(A) || b == null) {
    return false; // check inputs
  }
  n = A.n;
  S = sqr(order, A, false); // ordering and symbolic analysis
  N = lu(A, S, tol); // numeric LU factorization
  x = new Float64List(n); // get workspace
  ok = (S != null && N != null);
  if (ok) {
    ipvec(N.pinv, b, x, n); // x = b(p)
    lsolve(N.L, x); // x = L\x
    usolve(N.U, x); // x = U\x
    ipvec(S.q, x, b, n); // b(q) = x
  }
  return ok;
}
