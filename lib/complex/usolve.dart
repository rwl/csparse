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

/// Solves an upper triangular system Ux=b, where x and b are dense vectors.
/// The diagonal of U must be the last entry of each column.
///
/// [U] upper triangular matrix in column-compressed form.
/// [x] size n, right hand side on input, solution on output.
/// Returns true if successful, false on error.
bool cs_usolve(DZcs U, DZcsa x) {
  int n;
  Int32List Up, Ui;
  DZcsa Ux = new DZcsa();
  if (!CS_CSC(U) || x == null) {
    return false;
  }
  n = U.n;
  Up = U.p;
  Ui = U.i;
  Ux.x = U.x;
  for (int j = n - 1; j >= 0; j--) {
    x.set_list(j, cs_cdiv_list(x.get(j), Ux.get(Up[j + 1] - 1)));
    for (int p = Up[j]; p < Up[j + 1] - 1; p++) {
      x.set_list(Ui[p], cs_cminus(x.get(Ui[p]), cs_cmult_list(Ux.get(p), x.get(j))));
    }
  }
  return true;
}
