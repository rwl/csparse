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

/// Solves a lower triangular system Lx=b where x and b are dense. x=b on
/// input, solution on output.
///
/// [L] column-compressed, lower triangular matrix.
/// [x] size n, right hand side on input, solution on output.
/// Returns true if successful, false on error.
bool cs_lsolve(DZcs L, DZcsa x) {
  int n;
  Int32List Lp, Li;
  DZcsa Lx = new DZcsa();
  if (!CS_CSC(L) || x == null) {
    return false;
  }
  n = L.n;
  Lp = L.p;
  Li = L.i;
  Lx.x = L.x;
  for (int j = 0; j < n; j++) {
    x.set_list(j, cs_cdiv_list(x.get(j), Lx.get(Lp[j])));
    for (int p = Lp[j] + 1; p < Lp[j + 1]; p++) {
      x.set_list(Li[p], cs_cminus(x.get(Li[p]), cs_cmult_list(Lx.get(p), x.get(j))));
    }
  }
  return true;
}
