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

/// Permutes a sparse matrix, C = PAQ.
///
/// [A] m-by-n, column-compressed matrix.
/// [pinv] a permutation vector of length m.
/// [q] a permutation vector of length n.
/// [values] allocate pattern only if false, values and pattern otherwise.
/// Returns C = PAQ, null on error.
Matrix permute(Matrix A, Int32List pinv, Int32List q, bool values) {
  int t, j, k, m, n;
  int nz = 0;
  Int32List Ap, Ai, Cp, Ci;
  Vector Cx = new Vector(),
      Ax = new Vector();
  Matrix C;
  if (!csc(A)) {
    return null;
  }
  m = A.m;
  n = A.n;
  Ap = A.p;
  Ai = A.i;
  Ax.x = A.x;
  C = spalloc(m, n, Ap[n], values && Ax.x != null, false); // alloc result
  Cp = C.p;
  Ci = C.i;
  Cx.x = C.x;
  for (k = 0; k < n; k++) {
    Cp[k] = nz; // column k of C is column q[k] of A
    j = q != null ? (q[k]) : k;
    for (t = Ap[j]; t < Ap[j + 1]; t++) {
      if (Cx.x != null) {
        Cx.setList(nz, Ax.get(t)); // row i of A is row pinv[i] of C
      }
      Ci[nz++] = pinv != null ? (pinv[Ai[t]]) : Ai[t];
    }
  }
  Cp[n] = nz; // finalize the last column of C
  return C;
}
