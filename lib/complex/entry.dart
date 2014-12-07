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

/// Adds an entry to a triplet matrix. Memory-space and dimension of T are
/// increased if necessary.
///
/// [T] triplet matrix; new entry added on output.
/// [i] row index of new entry.
/// [j] column index of new entry.
/// [x] numerical value of new entry.
/// Returns true if successful, false otherwise.
bool cs_entry_list(DZcs T, int i, int j, Float64List x) {
  return cs_entry(T, i, j, x[0], x[1]);
}

/// Adds an entry to a triplet matrix. Memory-space and dimension of T are
/// increased if necessary.
///
/// [T] triplet matrix; new entry added on output.
/// [i] row index of new entry.
/// [j] column index of new entry.
/// [re] real value of new entry.
/// [im] imaginary value of new entry.
/// Returns true if successful, false otherwise.
bool cs_entry(DZcs T, int i, int j, double re, double im) {
  if (!CS_TRIPLET(T) || i < 0 || j < 0) {
    return false;
  }
  if ((T.nz >= T.nzmax) && !cs_sprealloc(T, 2 * (T.nzmax))) {
    return false;
  }
  if (T.x != null) {
    T.set(T.nz, re, im);
  }
  T.i[T.nz] = i;
  T.p[T.nz++] = j;
  T.m = math.max(T.m, i + 1);
  T.n = math.max(T.n, j + 1);
  return true;
}
