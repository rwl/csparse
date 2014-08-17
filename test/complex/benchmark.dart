/*package edu.emory.mathcs.csparsej.tdcomplex.test;

import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_qrsol.cs_qrsol;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_lusol.cs_lusol;

import java.io.InputStream;

import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsa;*/

//class Benchmark extends DZcs_test {

import 'dart:io';
import 'dart:typed_data';

import 'package:unittest/unittest.dart';
import 'package:csparse/complex/cxsparse.dart';
import 'package:csparse/complex/load.dart';
import 'package:csparse/complex/test_util.dart';

final int N = 100 ;

final int ORDER = 0 ;

void main() {
	int m ;
	double tol ;
	DZcs A, C ;
	DZcsa b, x ;

	DZproblem prob ;

	final file = get_file (C_IBM32B) ;
	prob = get_problem (file, DROP_TOL) ;

	A = prob.A ; C = prob.C ; b = prob.b ; x = prob.x ;
	m = A.m ;

	/* partial pivoting tolerance */
	tol = prob.sym != 0 ? 0.001 : 1 ;

	rhs (x, b, m) ;

	stdout.writeln("CSparse: c_ibm32a") ;

	benchmark(C, x, b, m, tol, ORDER) ;
}

void benchmark(DZcs C, DZcsa x, DZcsa b, int m, double tol, int order) {

	int t ;
	final t_lu = new Float64List(N) ;
	final t_qr = new Float64List(N) ;

	for (int i = 0; i < N; i++) {
		t = tic() ;
		cs_lusol (order, C, x, tol) ;
		t_lu [i] = toc (t) ;

		//System.arraycopy(b.x, 0, x.x, 0, m) ;
		x.x.setAll(0, b.x);

		t = tic() ;
		cs_qrsol (order, C, x) ;
		t_qr [i] = toc (t) ;

		//System.arraycopy(b.x, 0, x.x, 0, m) ;
		x.x.setAll(0, b.x);
	}

	stdout.write("LU - min: ${min (t_lu)}, max: ${max (t_lu)}, avg: ${avg (t_lu)}") ;
	stdout.writeln() ;

	stdout.write("QR - min: ${min (t_qr)}, max: ${max (t_qr)}, avg: ${avg (t_qr)}") ;
	stdout.writeln() ;
}

double min(Float64List tt) {
	int l = tt.length ;
	assert(l > 0) ;

	double tmin = tt [0] ;

	for (int i = 1 ; i < l ; i++) {
		if (tt [i] < tmin) tmin = tt [i] ;
	}

	return tmin ;
}

double max(Float64List tt) {
	int l = tt.length ;
	assert(l > 0) ;

	double tmax = tt [0] ;

	for (int i = 1 ; i < l ; i++) {
		if (tt [i] > tmax) tmax = tt [i] ;
	}

	return tmax ;
}

double avg(Float64List tt) {
	int l = tt.length ;
	assert(l > 0) ;

	double sum = 0.0 ;
	for (int i = 0; i < l; i++) {
		sum = sum + tt [i] ;
	}

	return sum / l ;
}

//}
