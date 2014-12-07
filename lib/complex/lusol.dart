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

/// Solves Ax=b, where A is square and nonsingular. b overwritten with
/// solution. Partial pivoting if tol = 1.
///
/// [order] ordering method to use (0 to 3).
/// [A] column-compressed matrix.
/// [b] size n, b on input, x on output.
/// [tol] partial pivoting tolerance.
/// Returns true if successful, false on error.
bool cs_lusol(int order, DZcs A, DZcsa b, double tol) {
  DZcsa x;
  DZcss S;
  DZcsn N;
  int n;
  bool ok;
  if (!CS_CSC(A) || b == null) {
    return false;
  }
  n = A.n;
  S = cs_sqr(order, A, false); // ordering and symbolic analysis
  N = cs_lu(A, S, tol); // numeric LU factorization
  x = new DZcsa.sized(n); // get workspace
  ok = (S != null && N != null);
  if (ok) {
    cs_ipvec(N.pinv, b, x, n); // x = b(p)
    cs_lsolve(N.L, x); // x = L\x
    cs_usolve(N.U, x); // x = U\x
    cs_ipvec(S.q, x, b, n); // b(q) = x
  }
  x = null;
  S = null;
  N = null;
  return ok;
}
