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

/// Add sparse matrices.
///
/// Returns C = alpha*A + beta*B, null on error
Dcs cs_add(Dcs A, Dcs B, double alpha, double beta) {
  int p, j, anz;
  int nz = 0;
  Int32List Cp, Ci, Bp;
  int m, n, bnz;
  Int32List w;
  Float64List x, Bx, Cx;
  bool values;
  Dcs C;
  if (!cs_csc(A) || !cs_csc(B)) {
    return null; // check inputs
  }
  if (A.m != B.m || A.n != B.n) {
    return null;
  }
  m = A.m;
  anz = A.p[A.n];
  n = B.n;
  Bp = B.p;
  Bx = B.x;
  bnz = Bp[n];
  w = new Int32List(m); // get workspace
  values = (A.x != null) && (Bx != null);
  x = values ? new Float64List(m) : null; // get workspace
  C = cs_spalloc(m, n, anz + bnz, values, false); // allocate result
  Cp = C.p;
  Ci = C.i;
  Cx = C.x;
  for (j = 0; j < n; j++) {
    Cp[j] = nz; // column j of C starts here
    nz = cs_scatter(A, j, alpha, w, x, j + 1, C, nz); // alpha*A(:,j)
    nz = cs_scatter(B, j, beta, w, x, j + 1, C, nz); // beta*B(:,j)
    if (values) {
      for (p = Cp[j]; p < nz; p++) {
        Cx[p] = x[Ci[p]];
      }
    }
  }
  Cp[n] = nz; // finalize the last column of C
  cs_sprealloc(C, 0); // remove extra space from C
  return C; // success; free workspace, return C
}
