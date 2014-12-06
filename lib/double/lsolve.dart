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

/// Solves a lower triangular system Lx=b where x and b are dense. x=b on
/// input, solution on output.
///
/// [L] column-compressed, lower triangular matrix.
/// [x] size n, right hand side on input, solution on output.
/// Returns true if successful, false on error.
bool cs_lsolve(Dcs L, Float64List x) {
  int n;
  Int32List Lp, Li;
  Float64List Lx;
  if (!cs_csc(L) || x == null) {
    return (false); // check inputs
  }
  n = L.n;
  Lp = L.p;
  Li = L.i;
  Lx = L.x;
  for (int j = 0; j < n; j++) {
    x[j] /= Lx[Lp[j]];
    for (int p = Lp[j] + 1; p < Lp[j + 1]; p++) {
      x[Li[p]] -= Lx[p] * x[j];
    }
  }
  return true;
}
