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

/// Convert a triplet form to compressed-column form.
///
/// C = compressed-column form of a triplet matrix T. The columns of C are
/// not sorted, and duplicate entries may be present in C.
/// Returns C if successful, null on error
Matrix compress(Matrix T) {
  int m, n, nz, p;
  Int32List Cp, Ci, w, Ti, Tj;
  Float64List Cx, Tx;
  Matrix C;
  if (!triplet(T)) {
    return null; // check inputs
  }
  m = T.m;
  n = T.n;
  Ti = T.i;
  Tj = T.p;
  Tx = T.x;
  nz = T.nz;
  C = spalloc(m, n, nz, Tx != null, false); // allocate result
  w = new Int32List(n); // get workspace
  Cp = C.p;
  Ci = C.i;
  Cx = C.x;
  for (int k = 0; k < nz; k++) {
    w[Tj[k]]++; // column counts
  }
  cumsum(Cp, w, n); // column pointers
  for (int k = 0; k < nz; k++) {
    Ci[p = w[Tj[k]]++] = Ti[k]; // A(i,j) is the pth entry in C
    if (Cx != null) {
      Cx[p] = Tx[k];
    }
  }
  return C;
}
