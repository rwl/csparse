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

/// Drops entries from a sparse matrix;
///
/// [A] column-compressed matrix.
/// [fkeep] drop aij if fkeep.fkeep(i,j,aij,other) is false.
/// [other] optional parameter to fkeep.
/// Returns the new number of entries in A, -1 on error.
int fkeep(Matrix A, ifkeep fkeep, Object other) {
  int nz = 0,
      n;
  Int32List Ap, Ai;
  Vector Ax = new Vector();
  if (!csc(A)) {
    return -1;
  }
  n = A.n;
  Ap = A.p;
  Ai = A.i;
  Ax.x = A.x;
  for (int j = 0; j < n; j++) {
    int p = Ap[j]; // get current location of col j
    Ap[j] = nz; // record new location of col j
    for ( ; p < Ap[j + 1]; p++) {
      if (fkeep(Ai[p], j, Ax.x != null ? Ax.get(p) : cone(), other)) {
        if (Ax.x != null) Ax.setList(nz, Ax.get(p)); // keep A(i,j)
        Ai[nz++] = Ai[p];
      }
    }
  }
  Ap[n] = nz; // finalize A
  sprealloc(A, 0); // remove extra space from A
  return nz;
}
