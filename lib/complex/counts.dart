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

part of edu.emory.mathcs.csparse.complex;

//import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs ;

//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_util.CS_CSC ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_util.cs_idone ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_transpose.cs_transpose ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_leaf.cs_leaf ;

/**
 * Column counts for Cholesky and QR.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
//public class DZcs_counts {

int _HEAD (int k, int j, Int32List head, int head_offset, bool ata)
{
	return ata ? head [head_offset + k] : j ;
}

int _NEXT (int J, Int32List next, int next_offset, bool ata)
{
	return ata ? next [next_offset + J] : -1 ;
}

Int32List _init_ata (DZcs AT, Int32List post, Int32List w)
{
	int i, k, p, m = AT.n, n = AT.m;
	Int32List ATp = AT.p, ATi = AT.i ;
	Int32List head = w ;
	int head_offset = 4 * n ;
	Int32List next = w ;
	int next_offset = 5 * n + 1 ;
	for (k = 0 ; k < n ; k++) w [post [k]] = k ;  /* invert post */
	for (i = 0 ; i < m ; i++)
	{
	  k = n;
		for (p = ATp [i] ; p < ATp [i + 1] ; p++) k = Math.min(k, w [ATi [p]]) ;
		next [next_offset + i] = head [head_offset + k] ;  /* place row i in linked list k */
		head [head_offset + k] = i ;
	}
	return new Int32List.fromList([head_offset, next_offset]) ;
}

/**
 * Column counts of LL'=A or LL'=A'A, given parent & postordering
 *
 * @param A
 *            column-compressed matrix
 * @param parent
 *            elimination tree of A
 * @param post
 *            postordering of parent
 * @param ata
 *            analyze A if false, A'A otherwise
 * @return column counts of LL'=A or LL'=A'A, null on error
 */
Int32List cs_counts(DZcs A, Int32List parent, Int32List post, bool ata)
{
	int i, j, k, n, m, J, s, p, q;
	Int32List ATp, ATi, maxfirst, prevleaf, ancestor, colcount, w, first, delta ;
	Int32List head = null, next = null ;
	Int32List jleaf = new Int32List(1) ;
	int head_offset = 0, next_offset = 0 ;
	DZcs AT ;
	if (!CS_CSC(A) || parent == null || post == null)
		return (null) ;				/* check inputs */
	m = A.m ; n = A.n ;
	s = 4*n + (ata ? (n+m+1) : 0) ;
	delta = colcount = new Int32List(n) ;		/* allocate result */
	w = new Int32List(s) ;				/* get workspace */
	AT = cs_transpose (A, false);			/* AT = A' */
	if (AT == null || colcount == null || w == null)
		return (cs_idone (colcount, AT, w, false)) ;
	ancestor = w ;
	maxfirst = w ;
	int maxfirst_offset = n ;
	prevleaf = w ;
	int prevleaf_offset = 2 * n ;
	first = w ;
	int first_offset = 3 * n ;
	for (k = 0 ; k < s ; k++) w [k] = -1 ;		/* clear workspace w [0..s-1] */
	for (k = 0 ; k < n ; k++) 			/* find first [j] */
	{
		j = post [k] ;
		delta [j] = (first [first_offset + j] == -1) ? 1 : 0 ; /* delta[j]=1 if j is a leaf */
		for ( ; j != -1 && first [first_offset + j] == -1 ; j = parent [j])
			first [first_offset + j] = k ;
	}
	ATp = AT.p ; ATi = AT.i ;
	if (ata)
	{
		Int32List offsets = init_ata(AT, post, w) ;
		head = w ;
		head_offset = offsets [0] ;
		next = w ;
		next_offset = offsets[1] ;
	}
	for (i = 0 ; i < n ; i++) ancestor [i] = i ;	/* each node in its own set */
	for (k = 0 ; k < n ; k++)
	{
		j = post [k] ;				/* j is the kth node in postordered etree */
		if (parent [j] != -1) delta [parent [j]]-- ;  /* j is not a root */
		for (J = HEAD (k, j, head, head_offset, ata) ; J != -1 ; J = NEXT (J, next, next_offset, ata)) /* J=j for LL'=A case */
		{
			for (p = ATp [J]; p < ATp [J + 1]; p++)
			{
				i = ATi [p] ;
				q = cs_leaf(i, j, first, first_offset, maxfirst, maxfirst_offset, prevleaf,
					prevleaf_offset, ancestor, 0, jleaf);
				if (jleaf [0] >= 1) delta [j]++ ;  /* A(i,j) is in skeleton */
				if (jleaf [0] == 2) delta [q]-- ;  /* account for overlap in q */
			}
		}
		if (parent [j] != -1) ancestor [j] = parent [j] ;
	}
	for (j = 0 ; j < n ; j++)			/* sum up delta's of each child */
	{
		if (parent [j] != -1) colcount [parent [j]] += colcount [j] ;
	}
	return (cs_idone (colcount, AT, w, true)) ;	/* success: free workspace */
}

//}
