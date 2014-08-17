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

/*package edu.emory.mathcs.csparsej.tdcomplex.test ;

import java.io.InputStream;
import java.util.Random;

import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_add.cs_add;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_chol.cs_chol;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cmult;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_complex.cs_cone;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_ipvec.cs_ipvec;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_lsolve.cs_lsolve;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_ltsolve.cs_ltsolve;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_multiply.cs_multiply;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_permute.cs_permute;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_pinv.cs_pinv;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_pvec.cs_pvec;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_schol.cs_schol;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_transpose.cs_transpose;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_updown.cs_updown;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_util.cs_spalloc;

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsa;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsn;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcss;*/

import 'dart:io';
import 'dart:math' as math;
import 'dart:typed_data';

import 'package:unittest/unittest.dart';
import 'package:csparse/complex/cxsparse.dart';
import 'package:csparse/complex/test_util.dart';

/**
 * Read a matrix, solve a linear system, update/downdate.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 *
 */
main() {

	/**
	 * Cholesky update/downdate
	 *
	 * @param prob
	 *            problem
	 * @return true if successful, false on error
	 */
	bool test3(DZproblem prob)
	{
		DZcs A, C, W = null, WW, WT, E = null, W2 ;
		int t, n, k, p1, p2;
		Int32List Li, Lp, Wi, Wp, p = null ;
		bool ok ;
		DZcsa b, x, resid, y = null, Lx, Wx ;
		Float64List s ;
		int t1 ;
		DZcss S = null ;
		DZcsn N = null ;
		if (prob == null || prob.sym == 0 || prob.A.n == 0) return (false) ;
		A = prob.A ; C = prob.C ; b = prob.b ; x = prob.x ; resid = prob.resid ;
		n = A.n ;
		if (prob.sym == 0 || n == 0) return (true) ;
		rhs (x, b, n) ;						/* compute right-hand side */
		stdout.write("\nchol then update/downdate ") ;
		print_order (1) ;
		y = new DZcsa.sized(n) ;
		t = tic() ;
		S = cs_schol (1, C) ;					/* symbolic Chol, amd(A+A') */
		stdout.write("\nsymbolic chol time ${toc (t)} ms\n") ;
		t = tic() ;
		N = cs_chol (C, S) ;					/* numeric Cholesky */
		stdout.write("numeric  chol time ${toc (t)} ms\n") ;
		if (S == null || N == null) return (false) ;
		t = tic() ;
		cs_ipvec (S.pinv, b, y, n) ;				/* y = P*b */
		cs_lsolve (N.L, y) ;					/* y = L\y */
		cs_ltsolve (N.L, y) ;					/* y = L'\y */
		cs_pvec (S.pinv, y, x, n) ;				/* x = P'*y */
		stdout.write("solve    chol time ${toc (t)} ms\n") ;
		stdout.write("original: ") ;
		print_resid (true, C, x, b, resid, prob) ;		/* print residual */
		k = n ~/ 2 ;						/* construct W  */
		W = cs_spalloc (n, 1, n, true, false) ;
		Lp = N.L.p ; Li = N.L.i ; Lx = new DZcsa(N.L.x) ;
		Wp = W.p ; Wi = W.i ; Wx = new DZcsa (W.x) ;
		Wp [0] = 0 ;
		p1 = Lp [k] ;
		Wp [1] = Lp [k+1] - p1 ;
		s = Lx.get(p1) ;
		math.Random r = new math.Random(1) ;
		for ( ; p1 < Lp [k+1] ; p1++)
		{
			p2 = p1 - Lp [k] ;
			Wi [p2] = Li [p1] ;
			Wx.set_list(p2, cs_cmult(s, r.nextDouble())) ;
		}
		t = tic() ;
		ok = cs_updown (N.L, 1, W, S.parent) ;			/* update: L*L'+W*W' */
		t1 = toc (t) ;
		stdout.write("update:   time: $t1 ms\n") ;
		if (!ok) return (false) ;
		t = tic() ;
		cs_ipvec (S.pinv, b, y, n) ;				/* y = P*b */
		cs_lsolve (N.L, y) ;					/* y = L\y */
		cs_ltsolve (N.L, y) ;					/* y = L'\y */
		cs_pvec (S.pinv, y, x, n) ;				/* x = P'*y */
		t = toc (t) ;
		p = cs_pinv (S.pinv, n) ;
		W2 = cs_permute (W, p, null, true) ;			/* E = C + (P'W)*(P'W)' */
		WT = cs_transpose (W2, true) ;
		WW = cs_multiply (W2, WT) ;
		WT = null ;
		W2 = null ;
		E = cs_add (C, WW, cs_cone(), cs_cone()) ;
		WW = null ;
		if (E == null || p == null) return (false) ;
		stdout.write("update:   time: ${t1 + t} ms(incl solve) ") ;
		print_resid (true, E, x, b, resid, prob) ;		/* print residual */
		N = null ;						/* clear N */
		t = tic() ;
		N = cs_chol (E, S) ;					/* numeric Cholesky */
		if (N == null) return (false) ;
		cs_ipvec (S.pinv, b, y, n) ;				/* y = P*b */
		cs_lsolve (N.L, y) ;					/* y = L\y */
		cs_ltsolve (N.L, y) ;					/* y = L'\y */
		cs_pvec (S.pinv, y, x, n) ;				/* x = P'*y */
		t = toc (t) ;
		stdout.write("rechol:   time: $t ms(incl solve) ") ;
		print_resid (true, E, x, b, resid, prob) ;		/* print residual */
		t = tic() ;
		ok = cs_updown (N.L, -1, W, S.parent) ;			/* downdate: L*L'-W*W' */
		t1 = toc (t) ;
		if (!ok) return (false) ;
		stdout.write("downdate: time: $t1\n") ;
		t = tic() ;
		cs_ipvec (S.pinv, b, y, n) ;				/* y = P*b */
		cs_lsolve (N.L, y) ;					/* y = L\y */
		cs_ltsolve (N.L, y) ;					/* y = L'\y */
		cs_pvec (S.pinv, y, x, n) ;				/* x = P'*y */
		t = toc (t) ;
		stdout.write("downdate: time: ${t1 + t} ms(incl solve) ") ;
		print_resid (true, C, x, b, resid, prob) ;		/* print residual */
		return (true) ;
	}

	test('c4', ()
	{
		final file = get_file (C4) ;
		DZproblem prob = get_problem (file, 0.0) ;

		assert_problem(prob, 4, 4, 10, -1, 16, 7.37e+01) ;

		test3(prob) ;

		double x_norm = 8.3862 ;
		expect(x_norm, closeTo(prob.norms[0], DELTA)) ;
		x_norm = 9.0064 ;
		expect(x_norm, closeTo(prob.norms[1], 0.3)) ;  // TODO: tighten accuracy
		expect(x_norm, closeTo(prob.norms[2], 0.3)) ;
		x_norm = 8.3862 ;
		expect(x_norm, closeTo(prob.norms[3], DELTA)) ;
	});

	test('mhd1280b', ()
	{
		final file = get_file (MHD1280B) ;
		DZproblem prob = get_problem (file, 0.0) ;

		assert_problem(prob, 1280, 1280, 12029, -1, 22778, 79.9740) ;

		test3(prob) ;

		double x_norm = 76270143066.4197 ;
		expect(x_norm, closeTo(prob.norms[0], DELTA)) ;
		expect(x_norm, closeTo(prob.norms[1], DELTA)) ;
		expect(x_norm, closeTo(prob.norms[2], DELTA)) ;
		expect(x_norm, closeTo(prob.norms[3], DELTA)) ;
	});

}
