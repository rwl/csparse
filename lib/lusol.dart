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
 * Solve Ax=b using sparse LU factorization.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
//public class Dcs_lusol {

/**
 * Solves Ax=b, where A is square and nonsingular. b overwritten with
 * solution. Partial pivoting if tol = 1.
 *
 * @param order
 *            ordering method to use (0 to 3)
 * @param A
 *            column-compressed matrix
 * @param b
 *            size n, b on input, x on output
 * @param tol
 *            partial pivoting tolerance
 * @return true if successful, false on error
 */
bool cs_lusol(int order, Dcs A, Float64List b, double tol) {
    Float64List x;
    Dcss S;
    Dcsn N;
    int n;
    bool ok;
    if (!CS_CSC(A) || b == null)
        return (false); /* check inputs */
    n = A.n;
    S = cs_sqr(order, A, false); /* ordering and symbolic analysis */
    N = cs_lu(A, S, tol); /* numeric LU factorization */
    x = new Float64List(n); /* get workspace */
    ok = (S != null && N != null);
    if (ok) {
        cs_ipvec(N.pinv, b, x, n); /* x = b(p) */
        cs_lsolve(N.L, x); /* x = L\x */
        cs_usolve(N.U, x); /* x = U\x */
        cs_ipvec(S.q, x, b, n); /* b(q) = x */
    }
    return (ok);
}
