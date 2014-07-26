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
 * Compute Householder reflection.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
//public class Dcs_house {

/**
 * Compute a Householder reflection, overwrite x with v, where
 * (I-beta*v*v')*x = s*e1. See Algo 5.1.1, Golub & Van Loan, 3rd ed.
 *
 * @param x
 *            x on output, v on input
 * @param x_offset
 *            the index of the first element in array x
 * @param beta
 *            scalar beta
 * @param n
 *            the length of x
 * @return norm2(x), -1 on error
 */
double cs_house(List<double> x, int x_offset, List<double> beta, int n) {
    double s, sigma = 0.0;
    int i;
    if (x == null || beta == null)
        return (-1.0); /* check inputs */
    for (i = 1; i < n; i++)
        sigma += x[x_offset + i] * x[x_offset + i];
    if (sigma == 0) {
        s = (x[x_offset + 0]).abs(); /* s = |x(0)| */
        beta[0] = (x[x_offset + 0] <= 0) ? 2.0 : 0.0;
        x[x_offset + 0] = 1.0;
    } else {
        s = math.sqrt(x[x_offset + 0] * x[x_offset + 0] + sigma); /* s = norm (x) */
        x[x_offset + 0] = (x[x_offset + 0] <= 0) ? (x[x_offset + 0] - s) : (-sigma / (x[x_offset + 0] + s));
        beta[0] = -1.0 / (s * x[x_offset + 0]);
    }
    return (s);
}
