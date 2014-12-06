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

/// Computes the transpose of a sparse matrix, C =A';
///
/// [values] pattern only if false, both pattern and values otherwise.
/// Returns C=A', null on error.
Matrix transpose(Matrix A, bool values) {
  int q;
  Int32List Cp, Ci;
  int n, m;
  Int32List Ap, Ai, w;
  Float64List Cx, Ax;
  Matrix C;
  if (!csc(A)) {
    return (null); // check inputs
  }
  m = A.m;
  n = A.n;
  Ap = A.p;
  Ai = A.i;
  Ax = A.x;
  C = spalloc(n, m, Ap[n], values && (Ax != null), false); // allocate result
  w = new Int32List(m); // get workspace
  Cp = C.p;
  Ci = C.i;
  Cx = C.x;
  for (int p = 0; p < Ap[n]; p++) {
    w[Ai[p]]++; // row counts
  }
  cumsum(Cp, w, m); // row pointers
  for (int j = 0; j < n; j++) {
    for (int p = Ap[j]; p < Ap[j + 1]; p++) {
      Ci[q = w[Ai[p]]++] = j; // place A(i,j) as entry C(j,i)
      if (Cx != null) {
        Cx[q] = Ax[p];
      }
    }
  }
  return C;
}
