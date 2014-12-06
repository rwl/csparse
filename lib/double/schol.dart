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

/// Ordering and symbolic analysis for a Cholesky factorization.
///
/// [order] ordering option (0 or 1).
/// Returns symbolic analysis for Cholesky, null on error.
Symbolic schol(int order, Matrix A) {
  int n;
  Int32List c, _post, P;
  Matrix C;
  Symbolic S;
  if (!csc(A)) {
    return null; // check inputs
  }
  n = A.n;
  S = new Symbolic(); // allocate result S
  P = amd(order, A); // P = amd(A+A'), or natural
  S.pinv = pinv(P, n); // find inverse permutation
  if (order != 0 && S.pinv == null) {
    return null;
  }
  C = symperm(A, S.pinv, false); // C = spones(triu(A(P,P)))
  S.parent = etree(C, false); // find etree of C
  _post = post(S.parent, n); // postorder the etree
  c = counts(C, S.parent, _post, false); // find column counts of chol(C)
  S.cp = new Int32List(n + 1); // allocate result S.cp
  S.unz = S.lnz = cumsum(S.cp, c, n); // find column pointers for L
  return (S.lnz >= 0) ? S : null;
}
