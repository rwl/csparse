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
library edu.emory.mathcs.csparse.cs_load;

import 'dart:io';
import 'csparse.dart' show Matrix, spalloc, entry;

/// Load a sparse matrix from [file].
///
/// Loads a triplet matrix T from [file]. Each line of the file contains
/// three values: a row index i, a column index j, and a numerical value aij.
///
/// If the file uses 1-based indexing simply set [base] to 1. Returns T if
/// successful, null on error.
Matrix load(File file, [int base = 0]) {
  final T = spalloc(0, 0, 1, true, true);
  try {
    for (String line in file.readAsLinesSync()) {
      List<String> tokens = line.trim().split(new RegExp(r"\s+"));
      if (tokens.length != 3) {
        return null;
      }
      final i = int.parse(tokens[0]) - base,
          j = int.parse(tokens[1]) - base,
          x = double.parse(tokens[2]);
      if (!entry(T, i, j, x)) {
        return null;
      }
    }
  } on FileSystemException catch (e) {
    print(e.message);
    return null;
  }
  return T;
}
