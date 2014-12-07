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

/// Prints a sparse matrix.
///
/// [A] sparse matrix (triplet ot column-compressed).
/// [brief] print all of A if false, a few entries otherwise.
/// Returns true if successful, false on error.
bool cs_print(DZcs A, bool brief, StringBuffer out) {
  int p, j, m, n, nzmax, nz;
  Int32List Ap, Ai;
  DZcsa Ax = new DZcsa();
  if (A == null) {
    out.write("(null)\n");
    return (false);
  }
  m = A.m;
  n = A.n;
  Ap = A.p;
  Ai = A.i;
  Ax.x = A.x;
  nzmax = A.nzmax;
  nz = A.nz;
  /*out.printf("CXSparseJ Version %d.%d.%d, %s.  %s\n",
		DZcs_common.CS_VER, DZcs_common.CS_SUBVER, DZcs_common.CS_SUBSUB,
		DZcs_common.CS_DATE, DZcs_common.CS_COPYRIGHT) ;*/
  if (nz < 0) {
    out.write("$m-by-$n, nzmax: $nzmax nnz: ${Ap[n]}, 1-norm: ${cs_norm (A)}\n");
    for (j = 0; j < n; j++) {
      out.write("    col $j : locations ${Ap [j]} to ${Ap [j+1] - 1}\n");
      for (p = Ap[j]; p < Ap[j + 1]; p++) {
        out.write("      ${Ai [p]} : (${Ax.x != null ? cs_creal (Ax.get(p)) : 1}, ${Ax.x != null ? cs_cimag (Ax.get(p)) : 0})\n");
        if (brief && p > 20) {
          out.write("  ...\n");
          return (true);
        }
      }
    }
  } else {
    out.write("triplet: $m-by-$n, nzmax: $nzmax nnz: $nz\n");
    for (p = 0; p < nz; p++) {
      out.write("    ${Ai [p]} ${Ap [p]} : (${Ax.x != null ? cs_creal (Ax.get(p)) : 1}, ${Ax.x != null ? cs_cimag (Ax.get(p)) : 0})\n");
      if (brief && p > 20) {
        out.write("  ...\n");
        return (true);
      }
    }
  }
  return (true);
}
