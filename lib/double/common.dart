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

const CS_VER = 1; // CSparse Version 1.0.0
const CS_SUBVER = 0;
const CS_SUBSUB = 0;
const CS_DATE = "June 13, 2009"; // CSparse release date
const CS_COPYRIGHT = "Copyright (c) Timothy A. Davis, 2006-2009";

/// Matrix in compressed-column or triplet form.
class Matrix {
  /// Maximum number of entries.
  int nzmax;

  /// number of rows
  int m;

  /// number of columns
  int n;

  /// column pointers (size n+1) or col indices (size nzmax)
  Int32List p;

  /// row indices, size nzmax
  Int32List i;

  /// numerical values, size nzmax
  Float64List x;

  /// # of entries in triplet matrix, -1 for compressed-col
  int nz;

  Matrix();
}

/// Output of symbolic Cholesky, LU, or QR analysis.
class Symbolic {
  /// inverse row perm. for QR, fill red. perm for Chol
  Int32List pinv;

  /// fill-reducing column permutation for LU and QR
  Int32List q;

  /// elimination tree for Cholesky and QR
  Int32List parent;

  /// column pointers for Cholesky, row counts for QR
  Int32List cp;

  /// leftmost[i] = min(find(A(i,:))), for QR
  Int32List leftmost;

  /// # of rows for QR, after adding fictitious rows
  int m2;

  /// # entries in L for LU or Cholesky; in V for QR
  int lnz;

  /// # entries in U for LU; in R for QR
  int unz;

  Symbolic();
}

/// Output of numeric Cholesky, LU, or QR factorization
class Numeric {
  /// L for LU and Cholesky, V for QR
  Matrix L;

  /// U for LU, R for QR, not used for Cholesky
  Matrix U;

  /// partial pivoting for LU
  Int32List pinv;

  /// beta [0..n-1] for QR
  Float64List B;

  Numeric();
}

/// Output of Dulmage-Mendelsohn decomposition.
class Decomposition {
  /// size m, row permutation
  Int32List p;

  /// size n, column permutation
  Int32List q;

  /// size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q)
  Int32List r;

  /// size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q)
  Int32List s;

  /// # of blocks in fine dmperm decomposition
  int nb;

  /// coarse row decomposition
  Int32List rr;

  /// coarse column decomposition
  Int32List cc;
}
