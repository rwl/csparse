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

/// Inverts a permutation vector. Returns pinv[i] = k if p[k] = i on input.
///
/// [p] a permutation vector if length n.
/// [n] length of p.
/// Returns pinv, null on error.
Int32List pinv(Int32List p, int n) {
  Int32List pinv;
  if (p == null) {
    return null; // p = NULL denotes identity
  }
  pinv = new Int32List(n); // allocate result
  for (int k = 0; k < n; k++) {
    pinv[p[k]] = k; // invert the permutation
  }
  return pinv; // return result
}
