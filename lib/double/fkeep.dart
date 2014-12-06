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

/// Drops entries from a sparse matrix;
///
/// [fkeep] drop aij if fkeep.fkeep(i,j,aij,other) is false
/// [other] optional parameter to fkeep
/// Returns the new number of entries in A, -1 on error.
int cs_fkeep(Dcs A, ifkeep fkeep, Object other) {
  int n;
  int nz = 0;
  Int32List Ap, Ai;
  Float64List Ax;
  if (!cs_csc(A)) {
    return -1; // check inputs
  }
  n = A.n;
  Ap = A.p;
  Ai = A.i;
  Ax = A.x;
  for (int j = 0; j < n; j++) {
    int p = Ap[j]; // get current location of col j
    Ap[j] = nz; // record new location of col j
    for ( ; p < Ap[j + 1]; p++) {
      if (fkeep(Ai[p], j, Ax != null ? Ax[p] : 1.0, other)) {
        if (Ax != null) {
          Ax[nz] = Ax[p]; // keep A(i,j)
        }
        Ai[nz++] = Ai[p];
      }
    }
  }
  Ap[n] = nz; // finalize A
  cs_sprealloc(A, 0); // remove extra space from A
  return nz;
}
