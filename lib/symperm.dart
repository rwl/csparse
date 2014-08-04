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
 * Symmetric permutation of a sparse matrix.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
//public class Dcs_symperm {

/**
 * Permutes a symmetric sparse matrix. C = PAP' where A and C are symmetric.
 *
 * @param A
 *            column-compressed matrix (only upper triangular part is used)
 * @param pinv
 *            size n, inverse permutation
 * @param values
 *            allocate pattern only if false, values and pattern otherwise
 * @return C = PAP', null on error
 */
Dcs cs_symperm(Dcs A, Int32List pinv, bool values) {
    int i, j, p, q, i2, j2, n;
    Int32List Ap, Ai, Cp, Ci, w;
    Float64List Cx, Ax;
    Dcs C;
    if (!CS_CSC(A))
        return (null); /* check inputs */
    n = A.n;
    Ap = A.p;
    Ai = A.i;
    Ax = A.x;
    C = cs_spalloc(n, n, Ap[n], values && (Ax != null), false); /* alloc result*/
    w = new Int32List(n); /* get workspace */
    Cp = C.p;
    Ci = C.i;
    Cx = C.x;
    for (j = 0; j < n; j++) /* count entries in each column of C */
    {
        j2 = pinv != null ? pinv[j] : j; /* column j of A is column j2 of C */
        for (p = Ap[j]; p < Ap[j + 1]; p++) {
            i = Ai[p];
            if (i > j)
                continue; /* skip lower triangular part of A */
            i2 = pinv != null ? pinv[i] : i; /* row i of A is row i2 of C */
            w[math.max(i2, j2)]++; /* column count of C */
        }
    }
    cs_cumsum(Cp, w, n); /* compute column pointers of C */
    for (j = 0; j < n; j++) {
        j2 = pinv != null ? pinv[j] : j; /* column j of A is column j2 of C */
        for (p = Ap[j]; p < Ap[j + 1]; p++) {
            i = Ai[p];
            if (i > j)
                continue; /* skip lower triangular part of A*/
            i2 = pinv != null ? pinv[i] : i; /* row i of A is row i2 of C */
            Ci[q = w[math.max(i2, j2)]++] = math.min(i2, j2);
            if (Cx != null)
                Cx[q] = Ax[p];
        }
    }
    return C;
}
