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
 * Various utilities.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
//public class Dcs_util {

/**
 * Allocate a sparse matrix (triplet form or compressed-column form).
 *
 * @param m
 *            number of rows
 * @param n
 *            number of columns
 * @param nzmax
 *            maximum number of entries
 * @param values
 *            allocate pattern only if false, values and pattern otherwise
 * @param triplet
 *            compressed-column if false, triplet form otherwise
 * @return sparse matrix
 */
Dcs cs_spalloc(int m, int n, int nzmax, bool values, bool triplet) {
    Dcs A = new Dcs(); /* allocate the Dcs struct */
    A.m = m; /* define dimensions and nzmax */
    A.n = n;
    A.nzmax = nzmax = math.max(nzmax, 1);
    A.nz = triplet ? 0 : -1; /* allocate triplet or comp.col */
    A.p = triplet ? new List<int>(nzmax) : new List<int>(n + 1);
    A.i = new List<int>(nzmax);
    A.x = values ? new List<double>(nzmax) : null;
    return A;
}

/**
 * Change the max # of entries a sparse matrix can hold.
 *
 * @param A
 *            column-compressed matrix
 * @param nzmax
 *            new maximum number of entries
 * @return true if successful, false on error
 */
bool cs_sprealloc(Dcs A, int nzmax) {
    if (A == null)
        return (false);
    if (nzmax <= 0)
        nzmax = (CS_CSC(A)) ? (A.p[A.n]) : A.nz;
    List<int> Ainew = new List<int>(nzmax);
    int length = math.min(nzmax, A.i.length);
    System.arraycopy(A.i, 0, Ainew, 0, length);
    A.i = Ainew;
    if (CS_TRIPLET(A)) {
        List<int> Apnew = new List<int>(nzmax);
        length = math.min(nzmax, A.p.length);
        System.arraycopy(A.p, 0, Apnew, 0, length);
        A.p = Apnew;
    }
    if (A.x != null) {
        List<double> Axnew = new List<double>(nzmax);
        length = math.min(nzmax, A.x.length);
        System.arraycopy(A.x, 0, Axnew, 0, length);
        A.x = Axnew;
    }
    A.nzmax = nzmax;
    return (true);
}

/**
 * Allocate a Dcsd object (a Dulmage-Mendelsohn decomposition).
 *
 * @param m
 *            number of rows of the matrix A to be analyzed
 * @param n
 *            number of columns of the matrix A to be analyzed
 * @return Dulmage-Mendelsohn decomposition
 */
Dcsd cs_dalloc(int m, int n) {
    Dcsd D;
    D = new Dcsd();
    D.p = new List<int>(m);
    D.r = new List<int>(m + 6);
    D.q = new List<int>(n);
    D.s = new List<int>(n + 6);
    D.cc = new List<int>(5);
    D.rr = new List<int>(5);
    return D;
}

int CS_FLIP(int i) {
    return (-(i) - 2);
}

int CS_UNFLIP(int i) {
    return (((i) < 0) ? CS_FLIP(i) : (i));
}

bool CS_MARKED(List<int> w, int j) {
    return (w[j] < 0);
}

void CS_MARK(List<int> w, int j) {
    w[j] = CS_FLIP(w[j]);
}

/**
 * Returns true if A is in column-compressed form, false otherwise.
 *
 * @param A
 *            sparse matrix
 * @return true if A is in column-compressed form, false otherwise
 */
bool CS_CSC(Dcs A) {
    return (A != null && (A.nz == -1));
}

/**
 * Returns true if A is in triplet form, false otherwise.
 *
 * @param A
 *            sparse matrix
 * @return true if A is in triplet form, false otherwise
 */
bool CS_TRIPLET(Dcs A) {
    return (A != null && (A.nz >= 0));
}
