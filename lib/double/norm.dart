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

/// Sparse matrix 1-norm.
///
/// Computes the 1-norm of a sparse matrix = max (sum (abs (A))), largest
/// column sum.
///
/// Returns the 1-norm if successful, -1 on error.
double cs_norm(Dcs A) {
  int n;
  Int32List Ap;
  Float64List Ax;
  double norm = 0.0,
      s;
  if (!cs_csc(A) || A.x == null) {
    return -1.0; // check inputs
  }
  n = A.n;
  Ap = A.p;
  Ax = A.x;
  for (int j = 0; j < n; j++) {
    s = 0.0;
    for (int p = Ap[j]; p < Ap[j + 1]; p++) {
      s += Ax[p].abs();
    }
    norm = math.max(norm, s);
  }
  return norm;
}
