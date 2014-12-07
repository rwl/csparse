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

/// Sparse matrix multiplication, C = A*B.
DZcs cs_multiply(DZcs A, DZcs B) {
  int p, j, anz, m, n, bnz;
  int nz = 0;
  Int32List Cp, Ci, Bp, w, Bi;
  DZcsa x,
      Bx = new DZcsa(),
      Cx = new DZcsa();
  bool values;
  DZcs C;
  if (!CS_CSC(A) || !CS_CSC(B)) {
    return null;
  }
  if (A.n != B.m) {
    return null;
  }
  m = A.m;
  anz = A.p[A.n];
  n = B.n;
  Bp = B.p;
  Bi = B.i;
  Bx.x = B.x;
  bnz = Bp[n];
  w = new Int32List(m); // get workspace
  values = (A.x != null) && (Bx.x != null);
  x = values ? new DZcsa.sized(m) : null; // get workspace
  C = cs_spalloc(m, n, anz + bnz, values, false); // allocate result
  if (C == null || w == null || (values && x == null)) {
    return cs_done(C, w, x, false);
  }
  Cp = C.p;
  for (j = 0; j < n; j++) {
    if (nz + m > C.nzmax && !cs_sprealloc(C, 2 * (C.nzmax) + m)) {
      return cs_done(C, w, x, false); // out of memory
    }
    Ci = C.i;
    Cx.x = C.x; // C.i and C.x may be reallocated
    Cp[j] = nz; // column j of C starts here
    for (p = Bp[j]; p < Bp[j + 1]; p++) {
      nz = cs_scatter(A, Bi[p], (Bx.x != null) ? Bx.get(p) : cs_cone(), w, x, j + 1, C, nz);
    }
    if (values) {
      for (p = Cp[j]; p < nz; p++) {
        Cx.set_list(p, x.get(Ci[p]));
      }
    }
  }
  Cp[n] = nz; // finalize the last column of C
  cs_sprealloc(C, 0); // remove extra space from C
  return (cs_done(C, w, x, true)); // success; free workspace, return C
}
