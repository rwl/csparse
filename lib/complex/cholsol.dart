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

/// Solves Ax=b where A is symmetric positive definite; b is overwritten with
/// solution.
///
/// [order] ordering method to use (0 or 1)
/// [A] column-compressed matrix, symmetric positive definite, only
/// upper triangular part is used
/// [b] right hand side, b is overwritten with solution
/// Returns true if successful, false on error.
bool cs_cholsol(int order, DZcs A, DZcsa b) {
  DZcsa x;
  DZcss S;
  DZcsn N;
  int n;
  bool ok;
  if (!CS_CSC(A) || b == null) {
    return false; // check inputs
  }
  n = A.n;
  S = cs_schol(order, A); // ordering and symbolic analysis
  N = cs_chol(A, S); // numeric Cholesky factorization
  x = new DZcsa.sized(n); // get workspace
  ok = (S != null && N != null && x != null);
  if (ok) {
    cs_ipvec(S.pinv, b, x, n); // x = P*b
    cs_lsolve(N.L, x); // x = L\x
    cs_ltsolve(N.L, x); // x = L'\x
    cs_pvec(S.pinv, x, b, n); // b = P'*x
  }
  return ok;
}
