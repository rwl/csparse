/*
 * CXSparse: a Concise Sparse matrix package.
 * Copyright (C) 2006-2011, Timothy A. Davis.
 * Copyright (C) 2011-2012, Richard W. Lincoln.
 * http://www.cise.ufl.edu/research/sparse/CXSparse
 *
 * -------------------------------------------------------------------------
 *
 * CXSparseJ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * CXSparseJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 *
 */

library edu.emory.mathcs.cxsparse.load;

import 'dart:io';

import 'cxsparse.dart';

//import java.io.BufferedReader ;
//import java.io.IOException ;
//import java.io.InputStream;
//import java.io.InputStreamReader;
//import java.io.Reader;

//import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs ;

//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_util.cs_spalloc ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_entry.cs_entry ;

/**
 * Load a sparse matrix from a file.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
//public class DZcs_load {

/**
 * Loads a triplet matrix T from a file. Each line of the file contains
 * four values: a row index i, a column index j, a real value aij, and an
 * imaginary value aij.
 * The file is zero-based.
 *
 * @param fileName
 *            file name
 * @return T if successful, null on error
 */
DZcs cs_load(File file, [int base=0])
{
  final T = cs_spalloc(0, 0, 1, true, true) ;	/* allocate result */

	try
	{
    for (String line in file.readAsLinesSync()) {
      List<String> tokens = line.trim().split(new RegExp(r"\s+"));
			if (tokens.length != 4) return null ;
			final i = int.parse(tokens [0]) ;
			final j = int.parse(tokens [1]) ;
			final x_re = double.parse(tokens [2]) ;
			final x_im = double.parse(tokens [3]) ;
			if (!cs_entry(T, i, j, x_re, x_im)) return (null) ;
		}
	}
	on FileSystemException catch (e)
	{
    print(e.message);
		return (null) ;
	}
	return (T) ;
}

//}
