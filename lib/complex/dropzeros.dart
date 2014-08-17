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

//import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs ;

//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_creal ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cimag ;
//import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_fkeep.cs_fkeep ;

/**
 * Drop zeros from a sparse matrix.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
//public class DZcs_dropzeros {

class _Cs_nonzero implements DZcs_ifkeep {

	bool fkeep(int i, int j, Float64List aij, Object other)
	{
		return ((cs_creal(aij) != 0) || (cs_cimag(aij) != 0)) ;
	}
}

/**
 * Removes numerically zero entries from a matrix.
 *
 * @param A
 *            column-compressed matrix
 * @return nz, new number of entries in A, -1 on error
 */
int cs_dropzeros(DZcs A)
{
	return (cs_fkeep(A, new _Cs_nonzero(), null));  /* keep all nonzero entries */
}

//}
