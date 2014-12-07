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

const int CS_VER = 2; // CXSparse Version 2.2.6
const int CS_SUBVER = 2;
const int CS_SUBSUB = 6;
const String CS_DATE = "Dec 15, 2011";
const String CS_COPYRIGHT = "Copyright (C) Timothy A. Davis, 2006-2011";

/// Complex array.
class DZcsa {
  /// Numerical values.
  Float64List x;

  DZcsa([this.x = null]);

  /// Constructs an array of the given length.
  factory DZcsa.sized(int len) {
    return new DZcsa(new Float64List(2 * len));
  }

  Float64List get(final int idx) {
    final offset = 2 * idx;
    return new Float64List.fromList([x[offset], x[offset + 1]]);
  }

  double real(final int idx) => x[2 * idx];

  double imag(final int idx) => x[(2 * idx) + 1];

  void set_list(final int idx, final Float64List val) {
    final offset = 2 * idx;

    x[offset] = val[0];
    x[offset + 1] = val[1];
  }

  void set(final int idx, final double re, final double im) {
    final offset = 2 * idx;

    x[offset] = re;
    x[offset + 1] = im;
  }

  String toString() {
    String s = "DZcsa [";
    for (int i = 0; i < x.length; i += 2) {
      if (i != 0) s += ", ";
      s += "${x[i]}+j${x[i + 1]}";
    }
    return s + "]";
  }
}

/// Complex matrix in compressed-column or triplet form.
class DZcs {
  /// Show a few entries in string representation.
  static bool BRIEF_PRINT = true;

  /// Maximum number of entries.
  int nzmax;

  /// Number of rows.
  int m;

  /// Number of columns.
  int n;

  /// Column pointers (size n+1) or col indices (size nzmax).
  Int32List p;

  /// Row indices, size nzmax.
  Int32List i;

  /// Numerical values, size 2 * nzmax.
  Float64List x;

  /// # of entries in triplet matrix, -1 for compressed-col.
  int nz;

  DZcs();

  Float64List get(final int idx) {
    final offset = 2 * idx;
    return new Float64List.fromList([x[offset], x[offset + 1]]);
  }

  void set_list(final int idx, final Float64List val) {
    set(idx, val[0], val[1]);
  }

  void set(final int idx, final double re, final double im) {
    final offset = 2 * idx;

    x[offset] = re;
    x[offset + 1] = im;
  }

  String toString() {
    StringBuffer out = new StringBuffer();
    cs_print(this, BRIEF_PRINT, out);
    return out.toString();
  }
}

/// Output of symbolic Cholesky, LU, or QR analysis.
class DZcss {
  /// Inverse row perm. for QR, fill red. perm for Chol.
  Int32List pinv;

  /// Fill-reducing column permutation for LU and QR.
  Int32List q;

  /// Elimination tree for Cholesky and QR.
  Int32List parent;

  /// Column pointers for Cholesky, row counts for QR.
  Int32List cp;

  /// leftmost[i] = min(find(A(i,:))), for QR.
  Int32List leftmost;

  /// # of rows for QR, after adding fictitious rows.
  int m2;

  /// # entries in L for LU or Cholesky; in V for QR.
  int lnz;

  /// # entries in U for LU; in R for QR.
  int unz;

  DZcss();
}

/// Output of numeric Cholesky, LU, or QR factorization.
class DZcsn {
  /// L for LU and Cholesky, V for QR.
  DZcs L;

  /// U for LU, R for QR, not used for Cholesky.
  DZcs U;

  /// Partial pivoting for LU.
  Int32List pinv;

  /// Beta [0..n-1] for QR.
  Float64List B;

  DZcsn();
}

/// Output of Dulmage-Mendelsohn decomposition.
class DZcsd {
  /// Size m, row permutation.
  Int32List p;

  /// Size n, column permutation.
  Int32List q;

  /// Size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q).
  Int32List r;

  /// Size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q).
  Int32List s;

  /// # of blocks in fine dmperm decomposition.
  int nb;

  /// Coarse row decomposition.
  Int32List rr;

  /// Coarse column decomposition.
  Int32List cc;
}
