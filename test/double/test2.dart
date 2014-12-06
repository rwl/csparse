/*
 * CSparse: a Concise Sparse matrix package.
 * Copyright (C) 2006-2011, Timothy A. Davis.
 * Copyright (C) 2011-2012, Richard W. Lincoln.
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 *
 */

import 'dart:typed_data';
import 'dart:io' show stdout;
import 'package:unittest/unittest.dart';
import 'package:csparse/double/csparse.dart';
import 'package:csparse/double/test_util.dart';

/// Read a matrix from a file and solve a linear system.
main() {

  /// Solves a linear system using Cholesky, LU, and QR, with various
  /// orderings.
  bool test2(Dproblem prob) {
    Matrix A, C;
    Float64List b, x, resid;
    double tol;
    int t, k, m, n, order, nb, ns, sprank;
    Int32List r, s, rr;
    bool ok;
    Decomposition D;
    if (prob == null) return (false);
    A = prob.A;
    C = prob.C;
    b = prob.b;
    x = prob.x;
    resid = prob.resid;
    m = A.m;
    n = A.n;
    /* partial pivoting tolerance */
    tol = prob.sym != 0 ? 0.001 : 1.0;
    /* randomized dmperm analysis */
    D = dmperm(C, 1);
    if (D == null) return (false);
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
    /* natural and amd(A'*A) */
    for (order = 0; order <= 3; order += 3)
    {
      if (order == 0 && m > 1000) continue;
      stdout.write("QR   ");
      print_order(order);
      /* compute right-hand side */
      rhs(x, b, m);
      t = tic();
      /* min norm(Ax-b) with QR */
      ok = qrsol(order, C, x);
      stdout.write("time: ${toc (t)} ms ");
      /* print residual */
      print_resid(ok, C, x, b, resid, prob);
    }
    /* return if rect. or singular*/
    if (m != n || sprank < n) return (true);
    /* try all orderings */
    for (order = 0; order <= 3; order++)
    {
      if (order == 0 && m > 1000) continue;
      stdout.write("LU   ");
      print_order(order);
      /* compute right-hand side */
      rhs(x, b, m);
      t = tic();
      /* solve Ax=b with LU */
      ok = lusol(order, C, x, tol);
      stdout.write("time: ${toc (t)} ms ");
      /* print residual */
      print_resid(ok, C, x, b, resid, prob);
    }
    if (prob.sym == 0) return (true);
    /* natural and amd(A+A') */
    for (order = 0; order <= 1; order++)
    {
      if (order == 0 && m > 1000) continue;
      stdout.write("Chol ");
      print_order(order);
      /* compute right-hand side */
      rhs(x, b, m);
      t = tic();
      /* solve Ax=b with Cholesky */
      ok = cholsol(order, C, x);
      stdout.write("time: ${toc (t)} ms ");
      /* print residual */
      print_resid(ok, C, x, b, resid, prob);
    }
    return (true);
  }

  test('ash219', () {
    final file = get_file(ASH219);
    Dproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 219, 85, 438, 0, 0, 9.0);
    assert_structure(prob, 1, 0, 85);

    expect(1.0052, closeTo(prob.norms[0], DELTA));
    expect(1.0052, closeTo(prob.norms[1], DELTA));
  });

  test('bcsstk01', () {
    final file = get_file(BCSSTK01);
    Dproblem prob = get_problem(file, DROP_TOL);

    test2(prob);
    assert_problem(prob, 48, 48, 224, -1, 400, 3.57094807469e+09);
    assert_structure(prob, 1, 0, 48);

    double x_norm = 0.0005;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
    expect(x_norm, closeTo(prob.norms[6], DELTA));
    expect(x_norm, closeTo(prob.norms[7], DELTA));
  });

  test('bcsstk16', () {
    final file = get_file(BCSSTK16);
    Dproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 4884, 4884, 147631, -1, 290378, 7.008379365769155e+09);
    assert_structure(prob, 75, 74, 4884);

    double x_norm = 1.9998;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
  });

  test('fs_183_1', () {
    final file = get_file(FS_183_1);
    Dproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 183, 183, 988, 0, 0, 1.7031774210073e+09);
    assert_dropped(prob, 71, 10);
    assert_structure(prob, 38, 37, 183);

    double x_norm = 212022.2099;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
  });

  test('ibm32a', () {
    final file = get_file(IBM32A);
    Dproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 32, 31, 123, 0, 0, 7.0);
    assert_structure(prob, 1, 0, 31);

    double x_norm = 5.5800;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
  });

  test('ibm32b', () {
    final file = get_file(IBM32B);
    Dproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 31, 32, 123, 0, 0, 8.0);
    assert_structure(prob, 1, 0, 31);

    double x_norm = 5.3348;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
  });

  test('lp_afiro', () {
    final file = get_file(LP_AFIRO);
    Dproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 27, 51, 102, 0, 0, 3.43);
    assert_structure(prob, 1, 0, 27);

    expect(2.4534, closeTo(prob.norms[0], DELTA));
    expect(2.4534, closeTo(prob.norms[1], DELTA));
  });

  test('mbeacxc', () {
    final file = get_file(MBEACXC);
    Dproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 492, 490, 49920, 0, 0, 9.29e-01);
    assert_structure(prob, 10, 8, 448);

    expect(prob.norms[0].isNaN, isTrue);
    expect(prob.norms[1].isNaN, isTrue);
  });

  test('t1', () {
    final file = get_file(T1);
    Dproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 4, 4, 10, 0, 0, 1.11e+01);
    assert_structure(prob, 1, 0, 4);

    double x_norm = 2.4550;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
  });

  test('west0067', () {
    final file = get_file(WEST0067);
    Dproblem prob = get_problem(file, DROP_TOL);

    test2(prob);

    assert_problem(prob, 67, 67, 294, 0, 0, 6.14);
    assert_structure(prob, 2, 1, 67);

    double x_norm = 21.9478;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
    expect(x_norm, closeTo(prob.norms[4], DELTA));
    expect(x_norm, closeTo(prob.norms[5], DELTA));
  });
}
