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
library edu.emory.mathcs.csparse.cs_print;

import 'dart:typed_data';
import 'dart:io' show stdout;
import 'csparse.dart';

/// Prints a sparse matrix.
///
/// [A] sparse matrix (triplet or column-compressed).
/// [brief] print all of A if false, a few entries otherwise.
/// Returns true if successful, false on error.
bool print(Matrix A, bool brief) {
  int p, j, m, n, nzmax, nz;
  Int32List Ap, Ai;
  Float64List Ax;
  if (A == null) {
    stdout.write("(null)\n");
    return (false);
  }
  m = A.m;
  n = A.n;
  Ap = A.p;
  Ai = A.i;
  Ax = A.x;
  nzmax = A.nzmax;
  nz = A.nz;
  stdout.write("CSparse Version $CS_VER.$CS_SUBVER.$CS_SUBSUB, $CS_DATE.  $CS_COPYRIGHT\n");
  if (nz < 0) {
    stdout.write("$m-by-$n, nzmax: $nzmax nnz: ${Ap[n]}, 1-norm: ${norm(A)}\n");
    for (j = 0; j < n; j++) {
      stdout.write("    col $j : locations ${Ap[j]} to ${Ap[j + 1] - 1}\n");
      for (p = Ap[j]; p < Ap[j + 1]; p++) {
        stdout.write("      ${Ai[p]} : ${Ax != null ? Ax[p] : 1}\n");
        if (brief && p > 20) {
          stdout.write("  ...\n");
          return true;
        }
      }
    }
  } else {
    stdout.write("triplet: $m-by-$n, nzmax: $nzmax nnz: $nz\n");
    for (p = 0; p < nz; p++) {
      stdout.write("    ${Ai[p]} ${Ap[p]} : ${Ax != null ? Ax[p] : 1}\n");
      if (brief && p > 20) {
        stdout.write("  ...\n");
        return true;
      }
    }
  }
  return true;
}
