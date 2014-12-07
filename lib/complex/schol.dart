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

//import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs ;
//import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcss ;

//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_util.CS_CSC ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_amd.cs_amd ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_pinv.cs_pinv ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_symperm.cs_symperm ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_etree.cs_etree ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_post.cs_post ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_counts.cs_counts ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_cumsum.cs_cumsum ;

/**
 * Symbolic Cholesky ordering and analysis.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
//public class DZcs_schol {

/**
 * Ordering and symbolic analysis for a Cholesky factorization.
 *
 * @param order
 *            ordering option (0 or 1)
 * @param A
 *            column-compressed matrix
 * @return symbolic analysis for Cholesky, null on error
 */
DZcss cs_schol(int order, DZcs A)
{
	int n;
	Int32List c, post, P ;
	DZcs C ;
	DZcss S ;
	if (!CS_CSC (A)) return (null);		/* check inputs */
	n = A.n ;
	S = new DZcss() ; 			/* allocate result S */
	if (S == null) return (null) ;		/* out of memory */
	P = cs_amd (order, A) ;			/* P = amd(A+A'), or natural */
	S.pinv = cs_pinv (P, n) ;		/* find inverse permutation */
	P = null ;
	if (order != 0 && S.pinv == null) return null ;
	C = cs_symperm (A, S.pinv, false) ;	/* C = spones(triu(A(P,P))) */
	S.parent = cs_etree (C, false) ;	/* find etree of C */
	post = cs_post (S.parent, n) ;		/* postorder the etree */
	c = cs_counts (C, S.parent, post, false) ;  /* find column counts of chol(C) */
	post = null;
	C = null ;
	S.cp = new Int32List(n+1) ;			/* allocate result S.cp */
	S.unz = S.lnz = cs_cumsum (S.cp, c, n).toInt() ;  /* find column pointers for L */
	c = null ;
	return ((S.lnz >= 0) ? S : null) ;
}

//}
