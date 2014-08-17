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
 * Lesser General License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 *
 */

part of edu.emory.mathcs.csparse.complex;

//import java.io.ByteArrayOutputStream;
//import java.io.UnsupportedEncodingException;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_print.cs_print;

/**
 * Common data structures.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
//class DZcs_common {

const int CS_VER = 2; /* CXSparseJ Version 2.2.6 */
const int CS_SUBVER = 2;
const int CS_SUBSUB = 6;
const String CS_DATE = "Dec 15, 2011"; /* CXSparseJ release date */
const String CS_COPYRIGHT = "Copyright (C) Timothy A. Davis, 2006-2011";

/**
 *
 * Complex array.
 *
 */
class DZcsa
{

	/**
	 * numerical values
	 */
	Float64List x;

	DZcsa()
	{

	}

	DZcsa(Float64List x)
	{
		this.x = x ;
	}

	/**
	 * Constructs an array of the given length.
	 */
	DZcsa(int len)
	{
		this.x = new Float64List(2*len) ;
	}

	/**
	 *
	 * @param idx
	 * @return
	 */
	Float64List get(final int idx)
	{
		int offset = 2 * idx ;
		return new Float64List.fromList([x [offset], x [offset + 1]]) ;
	}

	double real(final int idx) {
		return x [2 * idx] ;
	}

	double imag(final int idx) {
		return x [(2 * idx) + 1] ;
	}

	/**
	 *
	 * @param idx
	 * @param val
	 */
	void set(final int idx, final Float64List val)
	{
		int offset = 2 * idx ;

		x [offset] = val [0] ;
		x [offset + 1] = val [1] ;
	}

	void set(final int idx, final double re, final double im) {
		int offset = 2 * idx ;

		x [offset] = re ;
		x [offset + 1] = im ;
	}

	@Override
	String toString() {
		String s = "DZcsa [" ;
		for (int i = 0; i < x.length; i+=2) {
			if (i != 0) s += ", " ;
			s += String.format("%g+j%g", x[i], x[i + 1]) ;
		}
		return s + "]" ;
	}
}

/**
 *
 * Complex matrix in compressed-column or triplet form.
 *
 */
class DZcs
{

	/**
	 * show a few entries in string representation
	 */
	static bool BRIEF_PRINT = true;

	/**
	 * maximum number of entries
	 */
	int nzmax ;

	/**
	 * number of rows
	 */
	int m ;

	/**
	 * number of columns
	 */
	int n ;

	/**
	 * column pointers (size n+1) or col indices (size nzmax)
	 */
	Int32List p ;

	/**
	 * row indices, size nzmax
	 */
	Int32List i ;

	/**
	 * numerical values, size 2 * nzmax
	 */
	Float64List x ;

	/**
	 * # of entries in triplet matrix, -1 for compressed-col
	 */
	int nz ;

	DZcs()
	{

	}

	Float64List get(final int idx)
	{
		int offset = 2 * idx ;
		return new Float64List.fromList([x [offset], x [offset + 1]]) ;
	}

	void set(final int idx, final Float64List val) {
		set(idx, val [0], val [1]);
	}

	void set(final int idx, final double re, final double im)
	{
		int offset = 2 * idx ;

		x [offset] = re ;
		x [offset + 1] = im ;
	}

	String toString() {
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		cs_print(this, BRIEF_PRINT, out);
		try {
			return new String(out.toByteArray(), "UTF-8");
		} on UnsupportedEncodingException catch (e) {
			return out.toString();
		}
	}

}

/**
 *
 * Output of symbolic Cholesky, LU, or QR analysis.
 *
 */
class DZcss
{
	/**
	 * inverse row perm. for QR, fill red. perm for Chol
	 */
	Int32List pinv ;

	/**
	 * fill-reducing column permutation for LU and QR
	 */
	Int32List q ;

	/**
	 * elimination tree for Cholesky and QR
	 */
	Int32List parent ;

	/**
	 * column pointers for Cholesky, row counts for QR
	 */
	Int32List cp ;

	/**
	 * leftmost[i] = min(find(A(i,:))), for QR
	 */
	Int32List leftmost ;

	/**
	 * # of rows for QR, after adding fictitious rows
	 */
	int m2 ;

	/**
	 * # entries in L for LU or Cholesky; in V for QR
	 */
	int lnz ;

	/**
	 * # entries in U for LU; in R for QR
	 */
	int unz ;

	DZcss()
	{

	}
}

/**
 *
 * Output of numeric Cholesky, LU, or QR factorization
 *
 */
class DZcsn
{

	/**
	 * L for LU and Cholesky, V for QR
	 */
	DZcs L ;

	/**
	 * U for LU, R for QR, not used for Cholesky
	 */
	DZcs U ;

	/**
	 * partial pivoting for LU
	 */
	Int32List pinv ;

	/**
	 * beta [0..n-1] for QR
	 */
	Float64List B ;

	DZcsn()
	{

	}

}

/**
 *
 * Output of Dulmage-Mendelsohn decomposition.
 *
 */
class DZcsd
{

	/**
	 * size m, row permutation
	 */
	Int32List p ;

	/**
	 * size n, column permutation
	 */
	Int32List q ;

	/**
	 * size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q)
	 */
	Int32List r ;

	/**
	 * size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q)
	 */
	Int32List s ;

	/**
	 * # of blocks in fine dmperm decomposition
	 */
	int nb ;

	/**
	 * coarse row decomposition
	 */
	Int32List rr ;

	/**
	 * coarse column decomposition
	 */
	Int32List cc ;
}

//}
