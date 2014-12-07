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
part of edu.emory.mathcs.cxsparse;

/// Allocate a sparse matrix (triplet form or compressed-column form).
///
/// [m] number of rows. [n] number of columns.
/// [nzmax] maximum number of entries.
/// [values] allocate pattern only if false, values and pattern otherwise.
/// [triplet] compressed-column if false, triplet form otherwise.
DZcs cs_spalloc(int m, int n, int nzmax, bool values, bool triplet)
{
	DZcs A = new DZcs() ;			// allocate the DZcs struct
	A.m = m ;				// define dimensions and nzmax
	A.n = n ;
	A.nzmax = nzmax = math.max (nzmax, 1) ;
	A.nz = triplet ? 0 : -1 ;		// allocate triplet or comp.col
	A.p = triplet ? new Int32List(nzmax) : new Int32List(n+1) ;
	A.i = new Int32List(nzmax) ;
	A.x = values ? new Float64List(2*nzmax) : null ;
	return A ;
}

/// Change the max # of entries a sparse matrix can hold.
bool cs_sprealloc(DZcs A, int nzmax)
{
	if (A == null) {
	  return false;
	}
	if (nzmax <= 0) {
	  nzmax = (CS_CSC (A)) ? (A.p [A.n]) : A.nz ;
	}
	Int32List Ainew = new Int32List(nzmax) ;
	int length = math.min (nzmax, A.i.length) ;
	//System.arraycopy (A.i, 0, Ainew, 0, length) ;
	Ainew.setAll(0, A.i);
	A.i = Ainew ;
	if (CS_TRIPLET (A))
	{
		Int32List Apnew = new Int32List(nzmax) ;
		length = math.min (nzmax, A.p.length) ;
		//System.arraycopy (A.p, 0, Apnew, 0, length) ;
		Apnew.setAll(0, A.p);
		A.p = Apnew ;
	}
	if (A.x != null)
	{
		Float64List Axnew = new Float64List(2*nzmax) ;
		length = math.min (2*nzmax, A.x.length) ;
		//System.arraycopy (A.x, 0, Axnew, 0, length) ;
		Axnew.setAll(0, A.x);
		A.x = Axnew ;
	}
	A.nzmax = nzmax ;
	return (true) ;
}

/// Allocate a Dulmage-Mendelsohn decomposition.
///
/// [m] number of rows of the matrix A to be analyzed.
/// [n] number of columns of the matrix A to be analyzed.
DZcsd cs_dalloc(int m, int n)
{
	DZcsd D ;
	D = new DZcsd() ;
	D.p = new Int32List(m) ;
	D.r = new Int32List(m+6) ;
	D.q = new Int32List(n) ;
	D.s = new Int32List(n+6) ;
	D.cc = new Int32List(5) ;
	D.rr = new Int32List(5) ;
	return D ;
}

int _CS_FLIP(int i) => -i - 2;

int _CS_UNFLIP(int i) => (((i) < 0) ? _CS_FLIP (i) : (i)) ;

bool _CS_MARKED(Int32List w, int j) => w [j] < 0;

void _CS_MARK(Int32List w, int j)
{
	w [j] = _CS_FLIP (w [j]) ;
}

/// Returns true if A is in column-compressed form, false otherwise.
bool CS_CSC(DZcs A) => A != null && A.nz == -1;

/// Returns true if A is in triplet form, false otherwise.
bool CS_TRIPLET(DZcs A) => A != null && A.nz >= 0;

/// free workspace and return a sparse matrix result
DZcs cs_done (DZcs C, Int32List w, DZcsa x, bool ok) => ok ? C : null;

/// free workspace and return CS_INT array result
Int32List cs_idone (Int32List p, DZcs C, Int32List w, bool ok) => ok ? p : null;

/// free workspace and return a numeric factorization (Cholesky, LU, or QR)
DZcsn cs_ndone (DZcsn N, DZcs C, Int32List w, DZcsa x, bool ok) => ok ? N : null;

/// free workspace and return a csd result
DZcsd cs_ddone (DZcsd D, DZcs C, Int32List w, bool ok) => ok ? D : null;
