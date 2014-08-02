/* ***** BEGIN LICENSE BLOCK *****
 *
 * CSparse: a Concise Sparse matrix package.
 * Copyright (c) 2006, Timothy A. Davis.
 * http://www.cise.ufl.edu/research/sparse/CSparse
 *
 * -------------------------------------------------------------------------
 *
 * CSparse is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * CSparse is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * ***** END LICENSE BLOCK ***** */

part of edu.emory.mathcs.csparse;

/**
 * Common data structures.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
//class Dcs_common {

const CS_VER = 1; /* CSparse Version 1.0.0 */
const CS_SUBVER = 0;
const CS_SUBSUB = 0;
const CS_DATE = "June 13, 2009"; /* CSparse release date */
const CS_COPYRIGHT = "Copyright (c) Timothy A. Davis, 2006-2009";

/**
 *
 * Matrix in compressed-column or triplet form.
 *
 */
class Dcs {

    /**
     * maximum number of entries
     */
    int nzmax;

    /**
     * number of rows
     */
    int m;

    /**
     * number of columns
     */
    int n;

    /**
     * column pointers (size n+1) or col indices (size nzmax)
     */
    List<int> p;

    /**
     * row indices, size nzmax
     */
    List<int> i;

    /**
     * numerical values, size nzmax
     */
    List<double> x;

    /**
     * # of entries in triplet matrix, -1 for compressed-col
     */
    int nz;

    Dcs() {

    }

}

/**
 *
 * Output of symbolic Cholesky, LU, or QR analysis.
 *
 */
class Dcss {
    /**
     * inverse row perm. for QR, fill red. perm for Chol
     */
    List<int> pinv;

    /**
     * fill-reducing column permutation for LU and QR
     */
    List<int> q;

    /**
     * elimination tree for Cholesky and QR
     */
    List<int> parent;

    /**
     * column pointers for Cholesky, row counts for QR
     */
    List<int> cp;

    /**
     * leftmost[i] = min(find(A(i,:))), for QR
     */
    List<int> leftmost;

    /**
     * # of rows for QR, after adding fictitious rows
     */
    int m2;

    /**
     * # entries in L for LU or Cholesky; in V for QR
     */
    int lnz;

    /**
     * # entries in U for LU; in R for QR
     */
    int unz;

    Dcss() {
    }
}

/**
 *
 * Output of numeric Cholesky, LU, or QR factorization
 *
 */
class Dcsn {
    /**
     * L for LU and Cholesky, V for QR
     */
    Dcs L;

    /**
     * U for LU, R for QR, not used for Cholesky
     */
    Dcs U;

    /**
     * partial pivoting for LU
     */
    List<int> pinv;

    /**
     * beta [0..n-1] for QR
     */
    List<double> B;

    Dcsn() {
    }

}

/**
 *
 * Output of Dulmage-Mendelsohn decomposition.
 *
 */
class Dcsd {

    /**
     * size m, row permutation
     */
    List<int> p;

    /**
     * size n, column permutation
     */
    List<int> q;

    /**
     * size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q)
     */
    List<int> r;

    /**
     * size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q)
     */
    List<int> s;

    /**
     * # of blocks in fine dmperm decomposition
     */
    int nb;

    /**
     * coarse row decomposition
     */
    List<int> rr;

    /**
     * coarse column decomposition
     */
    List<int> cc;
}

