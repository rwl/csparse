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

/// Scatters and sums a sparse vector A(:,j) into a dense vector,
/// x = x + beta * A(:,j).
///
/// [A] the sparse vector is A(:,j).
/// [j] the column of A to use.
/// [beta] scalar multiplied by A(:,j).
/// [w] size m, node i is marked if w[i] = mark.
/// [x] size m, ignored if null.
/// [mark] mark value of w.
/// [C] pattern of x accumulated in C.i.
/// [nz] pattern of x placed in C starting at C.i[nz].
/// Returns new value of nz, -1 on error.
int scatter(Matrix A, int j, double beta, Int32List w, Float64List x, int mark, Matrix C, int nz) {
  Int32List Ap, Ai, Ci;
  Float64List Ax;
  if (!csc(A) || w == null || !csc(C)) {
    return -1; // check inputs
  }
  Ap = A.p;
  Ai = A.i;
  Ax = A.x;
  Ci = C.i;
  for (int p = Ap[j]; p < Ap[j + 1]; p++) {
    final i = Ai[p];
    /* A(i,j) is nonzero */
    if (w[i] < mark) {
      w[i] = mark; // i is new entry in column j
      Ci[nz++] = i; // add i to pattern of C(:,j)
      if (x != null) {
        x[i] = beta * Ax[p]; // x(i) = beta*A(i,j)
      }
    } else if (x != null) {
      x[i] += beta * Ax[p]; // i exists in C(:,j) already
    }
  }
  return nz;
}
