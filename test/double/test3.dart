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
import 'dart:math' as math;
import 'package:unittest/unittest.dart';
import 'package:csparse/double/csparse.dart';
import 'package:csparse/double/test_util.dart';

/// Read a matrix, solve a linear system, update/downdate.
main() {
  /// Cholesky update/downdate.
  bool test3(Dproblem prob) {
    Matrix A,
        C,
        W = null,
        WW,
        WT,
        E = null,
        W2;
    int n, k, p1, p2;
    Int32List Li,
        Lp,
        Wi,
        Wp,
        p = null;
    bool ok;
    Float64List b,
        x,
        resid,
        y = null,
        Lx,
        Wx;
    double s;
    int t, t1;
    Symbolic S = null;
    Numeric N = null;
    if (prob == null || prob.sym == 0 || prob.A.n == 0) return (false);
    A = prob.A;
    C = prob.C;
    b = prob.b;
    x = prob.x;
    resid = prob.resid;
    n = A.n;
    if (prob.sym == 0 || n == 0) return (true);
    /* compute right-hand side */
    rhs(x, b, n);
    print("\nchol then update/downdate ");
    print_order(1);
    y = new Float64List(n);
    t = tic();
    /* symbolic Chol, amd(A+A') */
    S = schol(1, C);
    print("\nsymbolic chol time ${toc (t)} ms\n");
    t = tic();
    /* numeric Cholesky */
    N = chol(C, S);
    print("numeric  chol time ${toc (t)} ms\n");
    if (S == null || N == null) return (false);
    t = tic();
    /* y = P*b */
    ipvec(S.pinv, b, y, n);
    /* y = L\y */
    lsolve(N.L, y);
    /* y = L'\y */
    ltsolve(N.L, y);
    /* x = P'*y */
    pvec(S.pinv, y, x, n);
    print("solve    chol time ${toc (t)} ms\n");
    print("original: ");
    /* print residual */
    print_resid(true, C, x, b, resid, prob);
    k = n ~/ 2;
    /* construct W  */
    W = spalloc(n, 1, n, true, false);
    Lp = N.L.p;
    Li = N.L.i;
    Lx = N.L.x;
    Wp = W.p;
    Wi = W.i;
    Wx = W.x;
    Wp[0] = 0;
    p1 = Lp[k];
    Wp[1] = Lp[k + 1] - p1;
    s = Lx[p1];
    final r = new math.Random(1);
    for ( ; p1 < Lp[k + 1]; p1++) {
      p2 = p1 - Lp[k];
      Wi[p2] = Li[p1];
      Wx[p2] = s * r.nextDouble();
    }
    t = tic();
    /* update: L*L'+W*W' */
    ok = updown(N.L, /*+*/1, W, S.parent);
    t1 = toc(t);
    print("update:   time: $t1 ms\n");
    if (!ok) return (false);
    t = tic();
    /* y = P*b */
    ipvec(S.pinv, b, y, n);
    /* y = L\y */
    lsolve(N.L, y);
    /* y = L'\y */
    ltsolve(N.L, y);
    /* x = P'*y */
    pvec(S.pinv, y, x, n);
    t = toc(t);
    p = pinv(S.pinv, n);
    /* E = C + (P'W)*(P'W)' */
    W2 = permute(W, p, null, true);
    WT = transpose(W2, true);
    WW = multiply(W2, WT);
    WT = null;
    W2 = null;
    E = add(C, WW, 1.0, 1.0);
    WW = null;
    if (E == null || p == null) return (false);
    print("update:   time: ${t1 + t} ms(incl solve) ");
    /* print residual */
    print_resid(true, E, x, b, resid, prob);
    /* clear N */
    N = null;
    t = tic();
    /* numeric Cholesky */
    N = chol(E, S);
    if (N == null) return (false);
    /* y = P*b */
    ipvec(S.pinv, b, y, n);
    /* y = L\y */
    lsolve(N.L, y);
    /* y = L'\y */
    ltsolve(N.L, y);
    /* x = P'*y */
    pvec(S.pinv, y, x, n);
    t = toc(t);
    print("rechol:   time: $t ms(incl solve) ");
    /* print residual */
    print_resid(true, E, x, b, resid, prob);
    t = tic();
    /* downdate: L*L'-W*W' */
    ok = updown(N.L, -1, W, S.parent);
    t1 = toc(t);
    if (!ok) return (false);
    print("downdate: time: $t1\n");
    t = tic();
    /* y = P*b */
    ipvec(S.pinv, b, y, n);
    /* y = L\y */
    lsolve(N.L, y);
    /* y = L'\y */
    ltsolve(N.L, y);
    /* x = P'*y */
    pvec(S.pinv, y, x, n);
    t = toc(t);
    print("downdate: time: ${t1 + t} ms(incl solve) ");
    /* print residual */
    print_resid(true, C, x, b, resid, prob);
    return (true);
  }

  test('bcsstk01', () {
    final file = get_file(BCSSTK01);
    Dproblem prob = get_problem(file);

    assert_problem(prob, 48, 48, 224, -1, 400, 3.5709480746974373e+09);

    test3(prob);

    double x_norm = 0.0005;
    expect(x_norm, closeTo(prob.norms[0], 1e-04));
    expect(x_norm, closeTo(prob.norms[1], 1e-04));
    expect(x_norm, closeTo(prob.norms[2], 1e-04));
    expect(x_norm, closeTo(prob.norms[3], 1e-04));
  });

  test('bcsstk16', () {
    final file = get_file(BCSSTK16);
    Dproblem prob = get_problem(file);

    assert_problem(prob, 4884, 4884, 147631, -1, 290378, 7.008379365769155e+09);

    test3(prob);

    double x_norm = 1.9998;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
  });
}
