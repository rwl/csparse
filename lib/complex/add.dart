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

/// Add sparse matrices.
///
/// C = alpha*A + beta*B
Matrix add(Matrix A, Matrix B, Float64List alpha, Float64List beta) {
  int nz = 0,
      anz;
  Int32List Cp, Ci, Bp, w;
  int m, n, bnz;
  Vector x,
      Bx = new Vector(),
      Cx = new Vector();
  bool values;
  Matrix C;
  if (!csc(A) || !csc(B)) {
    return null;
  }
  if (A.m != B.m || A.n != B.n) {
    return null;
  }
  m = A.m;
  anz = A.p[A.n];
  n = B.n;
  Bp = B.p;
  Bx.x = B.x;
  bnz = Bp[n];
  w = new Int32List(m); // get workspace
  values = (A.x != null) && (Bx.x != null);
  x = values ? new Vector.sized(m) : null;
  /// get workspace
  C = spalloc(m, n, anz + bnz, values, false); // allocate result
  Cp = C.p;
  Ci = C.i;
  Cx.x = C.x;
  for (int j = 0; j < n; j++) {
    Cp[j] = nz; // column j of C starts here
    nz = scatter(A, j, alpha, w, x, j + 1, C, nz); // alpha*A(:,j)
    nz = scatter(B, j, beta, w, x, j + 1, C, nz); // beta*B(:,j)
    if (values) {
      for (int p = Cp[j]; p < nz; p++) {
        Cx.setList(p, x.get(Ci[p]));
      }
    }
  }
  Cp[n] = nz; // finalize the last column of C
  sprealloc(C, 0); // remove extra space from C
  return C; // success; free workspace, return C
}
