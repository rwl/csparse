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

/// Sparse matrix times dense column vector, y = A*x+y.
///
/// [A] column-compressed matrix.
/// [x] size n vector.
/// [y] size m vector.
/// Returns true if successful, false on error.
bool cs_gaxpy(DZcs A, DZcsa x, DZcsa y) {
  int n;
  Int32List Ap, Ai;
  DZcsa Ax = new DZcsa();
  if (!CS_CSC(A) || x == null || y == null) {
    return false;
  }
  n = A.n;
  Ap = A.p;
  Ai = A.i;
  Ax.x = A.x;
  for (int j = 0; j < n; j++) {
    for (int p = Ap[j]; p < Ap[j + 1]; p++) {
      y.set_list(Ai[p], cs_cplus(y.get(Ai[p]), cs_cmult_list(Ax.get(p), x.get(j))));
    }
  }
  return true;
}
