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

/// Removes and sums duplicate entries in a sparse matrix.
///
/// Returns true if successful, false on error.
bool cs_dupl(Dcs A) {
  int p, q, n, m;
  int nz = 0;
  Int32List Ap, Ai, w;
  Float64List Ax;
  if (!cs_csc(A)) {
    return false;
  }
  /* check inputs */
  m = A.m;
  n = A.n;
  Ap = A.p;
  Ai = A.i;
  Ax = A.x;
  w = new Int32List(m); // get workspace
  for (int i = 0; i < m; i++) {
    w[i] = -1; // row i not yet seen
  }
  for (int j = 0; j < n; j++) {
    q = nz; // column j will start at q
    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      final i = Ai[p]; // A(i,j) is nonzero
      if (w[i] >= q) {
        Ax[w[i]] += Ax[p]; // A(i,j) is a duplicate
      } else {
        w[i] = nz; // record where row i occurs
        Ai[nz] = i; // keep A(i,j)
        Ax[nz++] = Ax[p];
      }
    }
    Ap[j] = q; // record start of column j
  }
  Ap[n] = nz; // finalize A
  return cs_sprealloc(A, 0); // remove extra space from A
}
