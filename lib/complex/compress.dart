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

/// Convert a triplet form to compressed-column form.
///
/// C = compressed-column form of a triplet matrix T. The columns of C are
/// not sorted, and duplicate entries may be present in C.
Matrix compress(Matrix T) {
  int m, n, nz, p, k;
  Int32List Cp, Ci, w, Ti, Tj;
  Vector Cx = new Vector(),
      Tx = new Vector();
  Matrix C;
  if (!triplet(T)) {
    return null;
  }
  m = T.m;
  n = T.n;
  Ti = T.i;
  Tj = T.p;
  Tx.x = T.x;
  nz = T.nz;
  C = spalloc(m, n, nz, Tx.x != null, false); // allocate result
  w = new Int32List(n); // get workspace
  if (C == null || w == null) {
    return _done(C, w, null, false); // out of memory
  }
  Cp = C.p;
  Ci = C.i;
  Cx.x = C.x;
  for (k = 0; k < nz; k++) w[Tj[k]]++; // column counts
  cumsum(Cp, w, n); // column pointers
  for (k = 0; k < nz; k++) {
    Ci[p = w[Tj[k]]++] = Ti[k]; // A(i,j) is the pth entry in C
    if (Cx.x != null) {
      Cx.setList(p, Tx.get(k));
    }
  }
  return _done(C, w, null, true); // success; free w and return C
}
