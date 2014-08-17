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

part of edu.emory.mathcs.cxsparse;

//import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsa;

//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_czero ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cplus ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cmult ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cdiv ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_conj ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_csqrt ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cequal ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cabs ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_creal ;

/**
 * Compute Householder reflection.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
//public class DZcs_house {

/**
 * Compute a Householder reflection [v,beta,s]=house(x), overwrite x with v,
 * where (I-beta*v*v')*x = s*e1 and e1 = [1 0 ... 0]'.
 * Note that this CXSparseJ version is different than CSparseJ.  See Higham,
 * Accuracy & Stability of Num Algorithms, 2nd ed, 2002, page 357.
 *
 * @param x
 *            x on output, v on input
 * @param beta
 *            scalar beta
 * @param n
 *            the length of x
 * @return norm2(x), -1 on error
 */
Float64List cs_house(DZcsa x, int x_offset, Float64List beta, int n)
{
	Float64List s = cs_czero() ;
	int i ;
	if (x == null) return new Float64List.fromList([-1.0, 0.0]) ;	/* check inputs */
	/* s = norm(x) */
	for (i = 0 ; i < n ; i++)  // TODO: check i = 1
		s = cs_cplus(s, cs_cmult_list(x.get(x_offset + i), cs_conj(x.get(x_offset + i)))) ;
	s = cs_csqrt(s) ;
	if (cs_cequal(s, cs_czero()))
	{
		beta [0] = 0.0 ;
		x.set(x_offset + 0, 1.0, 0.0) ;
	}
	else
	{
		/* s = sign(x[0]) * norm (x) ; */
		if (!cs_cequal(x.get(x_offset + 0), cs_czero()))
		{
			s = cs_cmult_list(s, cs_cdiv_list(x.get(x_offset + 0), new Float64List.fromList([cs_cabs_list(x.get(x_offset + 0)), 0.0]))) ;
		}
		x.set_list(x_offset + 0, cs_cplus(x.get(x_offset + 0), s)) ;
		beta [0] = 1 / cs_creal( cs_cmult_list(cs_conj(s), x.get(x_offset + 0)) ) ;
	}
	return cs_cmult(s, -1.0) ;
}

//}
