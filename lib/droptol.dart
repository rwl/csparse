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
 * Drop small entries from a sparse matrix.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
//public class Dcs_droptol {

class Cs_tol implements Dcs_ifkeep {
    bool fkeep(int i, int j, double aij, Object other) {
        return aij.abs() > (other as double);
    }
}

/**
 * Removes entries from a matrix with absolute value <= tol.
 *
 * @param A
 *            column-compressed matrix
 * @param tol
 *            drop tolerance
 * @return nz, new number of entries in A, -1 on error
 */
int cs_droptol(Dcs A, double tol) {
    return (cs_fkeep(A, new Cs_tol(), tol)); /* keep all large entries */
}
