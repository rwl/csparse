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

/// Cumulative sum.
///
/// p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c
///
/// [p] size n+1, cumulative sum of c.
/// [c] size n, overwritten with p [0..n-1] on output
/// [n] length of c.
/// Returns sum (c), null on error.
int cumsum(Int32List p, Int32List c, int n) {
  int nz = 0;
  double nz2 = 0.0;
  if (p == null || c == null) {
    return (-1); // check inputs
  }
  for (int i = 0; i < n; i++) {
    p[i] = nz;
    nz += c[i];
    nz2 += c[i]; // also in double to avoid int overflow
    c[i] = p[i]; // also copy p[0..n-1] back into c[0..n-1]
  }
  p[n] = nz;
  return nz2.toInt(); // return sum (c [0..n-1])
}
