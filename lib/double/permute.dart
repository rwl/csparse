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

/// Permutes a sparse matrix, C = PAQ.
///
/// [A] m-by-n, column-compressed matrix.
/// [pinv] a permutation vector of length m.
/// [q] a permutation vector of length n.
/// [values] allocate pattern only if false, values and pattern otherwise.
/// Returns C = PAQ, null on error.
Dcs cs_permute(Dcs A, Int32List pinv, Int32List q, bool values) {
  int m, n;
  int nz = 0;
  Int32List Ap, Ai, Cp, Ci;
  Float64List Cx, Ax;
  Dcs C;
  if (!cs_csc(A)) {
    return null; // check inputs
  }
  m = A.m;
  n = A.n;
  Ap = A.p;
  Ai = A.i;
  Ax = A.x;
  C = cs_spalloc(m, n, Ap[n], values && Ax != null, false); // alloc result
  Cp = C.p;
  Ci = C.i;
  Cx = C.x;
  for (int k = 0; k < n; k++) {
    Cp[k] = nz; // column k of C is column q[k] of A
    final j = q != null ? (q[k]) : k;
    for (int t = Ap[j]; t < Ap[j + 1]; t++) {
      if (Cx != null) {
        Cx[nz] = Ax[t]; // row i of A is row pinv[i] of C
      }
      Ci[nz++] = pinv != null ? (pinv[Ai[t]]) : Ai[t];
    }
  }
  Cp[n] = nz; // finalize the last column of C
  return C;
}
