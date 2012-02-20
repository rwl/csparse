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

package edu.emory.mathcs.csparsej.tdcomplex ;

/**
 * Cumulative sum.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
public class DZcs_cumsum {

	/**
	 * p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c
	 *
	 * @param p
	 *            size n+1, cumulative sum of c
	 * @param c
	 *            size n, overwritten with p [0..n-1] on output
	 * @param n
	 *            length of c
	 * @return sum (c), null on error
	 */
	public static double cs_cumsum(int [] p, int [] c, int n)
	{
		int i, nz = 0 ;
		double nz2 = 0 ;
		if (p == null || c == null)
			return (-1) ;		/* check inputs */
		for (i = 0 ; i < n ; i++)
		{
			p [i] = nz ;
			nz += c [i] ;
			nz2 += c [i] ;		/* also in double to avoid int overflow */
			c [i] = p [i] ;		/* also copy p[0..n-1] back into c[0..n-1]*/
		}
		p [n] = nz ;
		return (nz2) ;			/* return sum (c [0..n-1]) */
	}

}