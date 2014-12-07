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

import 'dart:io';
import 'dart:typed_data';

import 'package:unittest/unittest.dart';
import 'package:csparse/complex/cxsparse.dart';
import 'package:csparse/complex/test_util.dart';

/// Solves a linear system using Cholesky, LU, and QR, with various
/// orderings.
bool test2(DZproblem prob) {
  DZcs A, C;
  DZcsa b, x, resid;
  double tol;
  int t, k, m, n, order, nb, ns, sprank;
  Int32List r, s, rr;
  bool ok;
  DZcsd D;
  if (prob == null) {
    return false;
  }
  A = prob.A;
  C = prob.C;
  b = prob.b;
  x = prob.x;
  resid = prob.resid;
  m = A.m;
  n = A.n;
  tol = prob.sym != 0 ? 0.001 : 1; // partial pivoting tolerance
  D = cs_dmperm(C, 1); // randomized dmperm analysis
  if (D == null) {
    return false;
  }
  prob.nb = nb = D.nb;
  r = D.r;
  s = D.s;
  rr = D.rr;
  prob.sprank = sprank = rr[3];
  ns = 0;
  for (k = 0; k < nb; k++) {
    if ((r[k + 1] == r[k] + 1) && (s[k + 1] == s[k] + 1)) {
      ns++;
    }
  }
  prob.ns = ns;
  stdout.write("blocks: $nb singletons: $ns structural rank: $sprank\n");
  D = null;
  for (order = 0; order <= 3; order += 3) // natural and amd(A'*A)
  {
    if (order == 0 && m > 1000) continue;
    stdout.write("QR   ");
    print_order(order);
    rhs(x, b, m); // compute right-hand side
    t = tic();
    ok = cs_qrsol(order, C, x); // min norm(Ax-b) with QR
    stdout.write("time: ${toc (t)} ms ");
    print_resid(ok, C, x, b, resid, prob); // print residual
  }
  if (m != n || sprank < n) {
    return true; // return if rect. or singular
  }
  for (order = 0; order <= 3; order++) // try all orderings
  {
    if (order == 0 && m > 1000) continue;
    stdout.write("LU   ");
    print_order(order);
    rhs(x, b, m); // compute right-hand side
    t = tic();
    ok = cs_lusol(order, C, x, tol); // solve Ax=b with LU
    stdout.write("time: ${toc (t)} ms ");
    print_resid(ok, C, x, b, resid, prob); // print residual
  }
  if (prob.sym == 0) {
    return true;
  }
  for (order = 0; order <= 1; order++) // natural and amd(A+A')
  {
    if (order == 0 && m > 1000) continue;
    stdout.write("Chol ");
    print_order(order);
    rhs(x, b, m); // compute right-hand side
    t = tic();
    ok = cs_cholsol(order, C, x); // solve Ax=b with Cholesky
    stdout.write("time: ${toc (t)} ms ");
    print_resid(ok, C, x, b, resid, prob); // print residual
  }
  return true;
}

/// Read a matrix from a file and solve a linear system.
main() {
  test('c_ibm32a', () {
    final file = get_file(C_IBM32A);
    DZproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 32, 31, 123, 0, 0, 9.90);
    assert_structure(prob, 1, 0, 31);

    double x_norm = 3.9456;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
  });

  test('c_ibm32b', () {
    final file = get_file(C_IBM32B);
    DZproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 31, 32, 123, 0, 0, 11.313);
    assert_structure(prob, 1, 0, 31);

    double x_norm = 3.7723;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
  });

  test('c_mbeacxc', () {
    final file = get_file(C_MBEACXC);
    DZproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 492, 490, 49920, 0, 0, 9.29e-01);
    assert_structure(prob, 10, 8, 448);

    expect(double.NAN, equals(prob.norms[0]));
    expect(double.NAN, equals(prob.norms[1]));
  });

  test('c_west0067', () {
    final file = get_file(C_WEST0067);
    DZproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 67, 67, 294, 0, 0, 6.17);
    assert_structure(prob, 2, 1, 67);

    double x_norm = 21.4903;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
  });

  test('c4', () {
    final file = get_file(C4);
    DZproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 4, 4, 10, -1, 16, 7.37e+01);
    assert_structure(prob, 1, 0, 4);

    double x_norm = 8.3862;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));

    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));

    expect(x_norm, closeTo(prob.norms[6], DELTA));
    expect(x_norm, closeTo(prob.norms[7], DELTA));
  });

  test('czero', () {
    final file = get_file(CZERO);
    DZproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 1, 1, 0, 1, 0, 0.0);
    assert_dropped(prob, 1, 0);
    assert_structure(prob, 2, 0, 0);

    expect(double.NAN, equals(prob.norms[0]));
    expect(double.NAN, equals(prob.norms[1]));
  });

  test('mhd1280b', () {
    final file = get_file(MHD1280B);
    DZproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 1280, 1280, 11963, -1, 22646, 79.9740);
    assert_dropped(prob, 0, 66);
    assert_structure(prob, 20, 14, 1280);

    double delta = 1e-02;
    double x_norm = 76270143066.4161;
    expect(x_norm, closeTo(prob.norms[0], delta));
    x_norm = 76270143066.4197;
    expect(x_norm, closeTo(prob.norms[1], delta));
    expect(x_norm, closeTo(prob.norms[2], delta));
    expect(x_norm, closeTo(prob.norms[3], delta));
    expect(x_norm, closeTo(prob.norms[4], delta));
  });

  test('qc324', () {
    final file = get_file(QC324);
    DZproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 324, 324, 26730, 0, 0, 1.71);
    assert_structure(prob, 1, 0, 324);

    double x_norm = 6355.8643;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
  });

  test('t2', () {
    final file = get_file(T2);
    DZproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 4, 4, 10, 0, 0, 106.0752);
    assert_structure(prob, 1, 0, 4);

    double x_norm = 0.6623;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
  });

  test('t3', () {
    final file = get_file(T3);
    DZproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 3, 4, 12, 0, 0, 3.06);
    assert_structure(prob, 1, 0, 3);

    double x_norm = 1.5357;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
  });

  test('t4', () {
    final file = get_file(T4);
    DZproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 2, 2, 3, 1, 4, 2.83);
    assert_structure(prob, 1, 0, 2);

    double x_norm = 0.9014;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
  });

  test('young1c', () {
    final file = get_file(YOUNG1C);
    DZproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 841, 841, 4089, 0, 0, 730.46);
    assert_structure(prob, 1, 0, 841);

    double x_norm = 0.0509;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
  });
}
