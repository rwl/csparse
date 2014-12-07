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
library edu.emory.mathcs.cxsparse.load;

import 'dart:io';

import 'cxsparse.dart' as cx;

/// Loads a triplet matrix T from a file. Each line of the file contains
/// four values: a row index i, a column index j, a real value aij, and an
/// imaginary value aij.
cx.DZcs cs_load(File file, [int base = 0]) {
  final T = cx.cs_spalloc(0, 0, 1, true, true); // allocate result

  try {
    for (String line in file.readAsLinesSync()) {
      List<String> tokens = line.trim().split(new RegExp(r"\s+"));
      if (tokens.length != 4) return null;
      final i = int.parse(tokens[0]) - base;
      final j = int.parse(tokens[1]) - base;
      final x_re = double.parse(tokens[2]);
      final x_im = double.parse(tokens[3]);
      if (!cx.cs_entry(T, i, j, x_re, x_im)) return (null);
    }
  } on FileSystemException catch (e) {
    print(e.message);
    return (null);
  }
  return (T);
}
