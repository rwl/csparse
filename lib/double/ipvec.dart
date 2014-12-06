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

/// Permutes a vector, x = P'b.
///
/// [p] permutation vector, p=null denotes identity.
/// [b] input vector. [x] output vector, x = P'b.
/// [n] length of p, b, and x.
/// Returns true if successful, false on error.
bool ipvec(Int32List p, Float64List b, Float64List x, int n) {
  if (x == null || b == null) {
    return false; // check inputs
  }
  for (int k = 0; k < n; k++) {
    x[p != null ? p[k] : k] = b[k];
  }
  return true;
}
