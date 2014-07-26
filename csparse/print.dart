/* ***** BEGIN LICENSE BLOCK *****
 *
 * CSparse: a Concise Sparse matrix package.
 * Copyright (c) 2006, Timothy A. Davis.
 * http://www.cise.ufl.edu/research/sparse/CSparse
 *
 * -------------------------------------------------------------------------
 *
 * CSparseJ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * CSparseJ is distributed in the hope that it will be useful,
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
 * Print a sparse matrix.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
//public class Dcs_print {

/**
 * Prints a sparse matrix.
 *
 * @param A
 *            sparse matrix (triplet ot column-compressed)
 * @param brief
 *            print all of A if false, a few entries otherwise
 * @return true if successful, false on error
 */
bool cs_print(Dcs A, bool brief) {
    int p, j, m, n, nzmax, nz;
    List<int> Ap, Ai;
    List<double> Ax;
    if (A == null) {
        print("(null)\n");
        return (false);
    }
    m = A.m;
    n = A.n;
    Ap = A.p;
    Ai = A.i;
    Ax = A.x;
    nzmax = A.nzmax;
    nz = A.nz;
    print("CSparseJ Version $CS_VER.$CS_SUBVER.$CS_SUBSUB, $CS_DATE.  $CS_COPYRIGHT\n");
    if (nz < 0) {
        print("$m-by-$n, nzmax: $nzmax nnz: ${Ap[n]}, 1-norm: ${cs_norm(A)}\n");
        for (j = 0; j < n; j++) {
            print("    col $j : locations ${Ap[j]} to ${Ap[j + 1] - 1}\n");
            for (p = Ap[j]; p < Ap[j + 1]; p++) {
                print("      ${Ai[p]} : ${Ax != null ? Ax[p] : 1}\n");
                if (brief && p > 20) {
                    print("  ...\n");
                    return (true);
                }
            }
        }
    } else {
        print("triplet: $m-by-$n, nzmax: $nzmax nnz: $nz\n");
        for (p = 0; p < nz; p++) {
            print("    ${Ai[p]} ${Ap[p]} : ${Ax != null ? Ax[p] : 1}\n");
            if (brief && p > 20) {
                print("  ...\n");
                return (true);
            }
        }
    }
    return (true);
}
