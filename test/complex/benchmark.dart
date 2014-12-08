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
library cxsparse.test.benchmark;

import 'dart:io';
import 'dart:typed_data';

import 'package:csparse/complex/cxsparse.dart';
import 'package:csparse/complex/test_util.dart';

final int N = 100;

final int ORDER = 0;

void main() {
  int m;
  double tol;
  Matrix A, C;
  Vector b, x;

  Problem prob;

  final file = getFile(C_IBM32B);
  prob = getProblem(file, DROP_TOL);

  A = prob.A;
  C = prob.C;
  b = prob.b;
  x = prob.x;
  m = A.m;

  /* partial pivoting tolerance */
  tol = prob.sym != 0 ? 0.001 : 1;

  rhs(x, b, m);

  stdout.writeln("CSparse: c_ibm32a");

  benchmark(C, x, b, m, tol, ORDER);
}

void benchmark(Matrix C, Vector x, Vector b, int m, double tol, int order) {
  int t;
  final t_lu = new Float64List(N);
  final t_qr = new Float64List(N);

  for (int i = 0; i < N; i++) {
    t = tic();
    lusol(order, C, x, tol);
    t_lu[i] = toc(t);

    //System.arraycopy(b.x, 0, x.x, 0, m) ;
    x.x.setAll(0, b.x);

    t = tic();
    qrsol(order, C, x);
    t_qr[i] = toc(t);

    //System.arraycopy(b.x, 0, x.x, 0, m) ;
    x.x.setAll(0, b.x);
  }

  stdout.write("LU - min: ${min (t_lu)}, max: ${max (t_lu)}, avg: ${avg (t_lu)}");
  stdout.writeln();

  stdout.write("QR - min: ${min (t_qr)}, max: ${max (t_qr)}, avg: ${avg (t_qr)}");
  stdout.writeln();
}

double min(Float64List tt) {
  int l = tt.length;
  assert(l > 0);

  double tmin = tt[0];

  for (int i = 1; i < l; i++) {
    if (tt[i] < tmin) tmin = tt[i];
  }

  return tmin;
}

double max(Float64List tt) {
  int l = tt.length;
  assert(l > 0);

  double tmax = tt[0];

  for (int i = 1; i < l; i++) {
    if (tt[i] > tmax) tmax = tt[i];
  }

  return tmax;
}

double avg(Float64List tt) {
  int l = tt.length;
  assert(l > 0);

  double sum = 0.0;
  for (int i = 0; i < l; i++) {
    sum = sum + tt[i];
  }

  return sum / l;
}
