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
    A.p = triplet ? new Int32List(nzmax) : new Int32List(n + 1);
    A.i = new Int32List(nzmax);
    A.x = values ? new Float64List(nzmax) : null;
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
    Int32List Ainew = new Int32List(nzmax);
    int length = math.min(nzmax, A.i.length);
    //System.arraycopy(A.i, 0, Ainew, 0, length);
    for (int i = 0; i < length; i++) Ainew[i] = A.i[i];
    //Ainew.setAll(0, length, A.i);
    A.i = Ainew;
    if (CS_TRIPLET(A)) {
        Int32List Apnew = new Int32List(nzmax);
        length = math.min(nzmax, A.p.length);
        //System.arraycopy(A.p, 0, Apnew, 0, length);
        //Apnew.replaceRange(0, length, A.p);
        for (int i = 0; i < length; i++) Apnew[i] = A.p[i];
        A.p = Apnew;
    }
    if (A.x != null) {
        Float64List Axnew = new Float64List(nzmax);
        length = math.min(nzmax, A.x.length);
        //System.arraycopy(A.x, 0, Axnew, 0, length);
        //Axnew.replaceRange(0, length, A.x);
        for (int i = 0; i < length; i++) Axnew[i] = A.x[i];
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
    D.p = new Int32List(m);
    D.r = new Int32List(m + 6);
    D.q = new Int32List(n);
    D.s = new Int32List(n + 6);
    D.cc = new Int32List(5);
    D.rr = new Int32List(5);
    return D;
}

int CS_FLIP(int i) {
    return (-(i) - 2);
}

int CS_UNFLIP(int i) {
    return (((i) < 0) ? CS_FLIP(i) : (i));
}

bool CS_MARKED(Int32List w, int j) {
    return (w[j] < 0);
}

void CS_MARK(Int32List w, int j) {
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
