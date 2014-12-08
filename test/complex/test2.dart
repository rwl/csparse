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
library cxsparse.test.test2;

import 'dart:io';
import 'dart:typed_data';

import 'package:unittest/unittest.dart';
import 'package:csparse/complex/cxsparse.dart';
import 'package:csparse/complex/test_util.dart';

/// Solves a linear system using Cholesky, LU, and QR, with various
/// orderings.
bool test2(Problem prob) {
  Matrix A, C;
  Vector b, x, resid;
  double tol;
  int t, k, m, n, order, nb, ns, sprank;
  Int32List r, s, rr;
  bool ok;
  Decomposition D;
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
  tol = prob.sym != 0 ? 0.001 : 1.0; // partial pivoting tolerance
  D = dmperm(C, 1); // randomized dmperm analysis
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
    printOrder(order);
    rhs(x, b, m); // compute right-hand side
    t = tic();
    ok = qrsol(order, C, x); // min norm(Ax-b) with QR
    stdout.write("time: ${toc (t)} ms ");
    printResid(ok, C, x, b, resid, prob); // print residual
  }
  if (m != n || sprank < n) {
    return true; // return if rect. or singular
  }
  for (order = 0; order <= 3; order++) // try all orderings
  {
    if (order == 0 && m > 1000) continue;
    stdout.write("LU   ");
    printOrder(order);
    rhs(x, b, m); // compute right-hand side
    t = tic();
    ok = lusol(order, C, x, tol); // solve Ax=b with LU
    stdout.write("time: ${toc (t)} ms ");
    printResid(ok, C, x, b, resid, prob); // print residual
  }
  if (prob.sym == 0) {
    return true;
  }
  for (order = 0; order <= 1; order++) // natural and amd(A+A')
  {
    if (order == 0 && m > 1000) continue;
    stdout.write("Chol ");
    printOrder(order);
    rhs(x, b, m); // compute right-hand side
    t = tic();
    ok = cholsol(order, C, x); // solve Ax=b with Cholesky
    stdout.write("time: ${toc (t)} ms ");
    printResid(ok, C, x, b, resid, prob); // print residual
  }
  return true;
}

/// Read a matrix from a file and solve a linear system.
main() {
  test('c_ibm32a', () {
    final file = getFile(C_IBM32A);
    Problem prob = getProblem(file, DROP_TOL);

    test2(prob);

    assertProblem(prob, 32, 31, 123, 0, 0, 9.90);
    assertStructure(prob, 1, 0, 31);

    double x_norm = 3.9456;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
  });

  test('c_ibm32b', () {
    final file = getFile(C_IBM32B);
    Problem prob = getProblem(file, DROP_TOL);

    test2(prob);

    assertProblem(prob, 31, 32, 123, 0, 0, 11.313);
    assertStructure(prob, 1, 0, 31);

    double x_norm = 3.7723;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
  });

  test('c_mbeacxc', () {
    final file = getFile(C_MBEACXC);
    Problem prob = getProblem(file, DROP_TOL);

    test2(prob);

    assertProblem(prob, 492, 490, 49920, 0, 0, 9.29e-01);
    assertStructure(prob, 10, 8, 448);

    expect(double.NAN, equals(prob.norms[0]));
    expect(double.NAN, equals(prob.norms[1]));
  });

  test('c_west0067', () {
    final file = getFile(C_WEST0067);
    Problem prob = getProblem(file, DROP_TOL);

    test2(prob);

    assertProblem(prob, 67, 67, 294, 0, 0, 6.17);
    assertStructure(prob, 2, 1, 67);

    double x_norm = 21.4903;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
  });

  test('c4', () {
    final file = getFile(C4);
    Problem prob = getProblem(file, DROP_TOL);

    test2(prob);

    assertProblem(prob, 4, 4, 10, -1, 16, 7.37e+01);
    assertStructure(prob, 1, 0, 4);

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
    final file = getFile(CZERO);
    Problem prob = getProblem(file, DROP_TOL);

    test2(prob);

    assertProblem(prob, 1, 1, 0, 1, 0, 0.0);
    assertDropped(prob, 1, 0);
    assertStructure(prob, 2, 0, 0);

    expect(prob.norms[0].isNaN, isTrue);
    expect(prob.norms[1].isNaN, isTrue);
  });

  test('mhd1280b', () {
    final file = getFile(MHD1280B);
    Problem prob = getProblem(file, DROP_TOL);

    test2(prob);

    assertProblem(prob, 1280, 1280, 11963, -1, 22646, 79.9740);
    assertDropped(prob, 0, 66);
    assertStructure(prob, 20, 14, 1280);

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
    final file = getFile(QC324);
    Problem prob = getProblem(file, DROP_TOL);

    test2(prob);

    assertProblem(prob, 324, 324, 26730, 0, 0, 1.71);
    assertStructure(prob, 1, 0, 324);

    double x_norm = 6355.8643;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
  });

  test('t2', () {
    final file = getFile(T2);
    Problem prob = getProblem(file, DROP_TOL);

    test2(prob);

    assertProblem(prob, 4, 4, 10, 0, 0, 106.0752);
    assertStructure(prob, 1, 0, 4);

    double x_norm = 0.6623;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
  });

  test('t3', () {
    final file = getFile(T3);
    Problem prob = getProblem(file, DROP_TOL);

    test2(prob);

    assertProblem(prob, 3, 4, 12, 0, 0, 3.06);
    assertStructure(prob, 1, 0, 3);

    double x_norm = 1.5357;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
  });

  test('t4', () {
    final file = getFile(T4);
    Problem prob = getProblem(file, DROP_TOL);

    test2(prob);

    assertProblem(prob, 2, 2, 3, 1, 4, 2.83);
    assertStructure(prob, 1, 0, 2);

    double x_norm = 0.9014;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
  });

  test('young1c', () {
    final file = getFile(YOUNG1C);
    Problem prob = getProblem(file, DROP_TOL);

    test2(prob);

    assertProblem(prob, 841, 841, 4089, 0, 0, 730.46);
    assertStructure(prob, 1, 0, 841);

    double x_norm = 0.0509;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
  });
}
