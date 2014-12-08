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
library cxsparse.test.test3;

import 'dart:io';
import 'dart:math' as math;
import 'dart:typed_data';

import 'package:unittest/unittest.dart';
import 'package:csparse/complex/cxsparse.dart';
import 'package:csparse/complex/test_util.dart';

/// Cholesky update/downdate.
bool test3(Problem prob) {
  Matrix A, C, W, WW, WT, E, W2;
  int t, n, k, p1, p2;
  Int32List Li, Lp, Wi, Wp, p;
  bool ok;
  Vector b,
      x,
      resid,
      y = null,
      Lx,
      Wx;
  Float64List s;
  int t1;
  Symbolic S = null;
  Numeric N = null;
  if (prob == null || prob.sym == 0 || prob.A.n == 0) {
    return false;
  }
  A = prob.A;
  C = prob.C;
  b = prob.b;
  x = prob.x;
  resid = prob.resid;
  n = A.n;
  if (prob.sym == 0 || n == 0) {
    return true;
  }
  rhs(x, b, n); // compute right-hand side
  stdout.write("\nchol then update/downdate ");
  printOrder(1);
  y = new Vector.sized(n);
  t = tic();
  S = schol(1, C); // symbolic Chol, amd(A+A')
  stdout.write("\nsymbolic chol time ${toc (t)} ms\n");
  t = tic();
  N = chol(C, S); // numeric Cholesky
  stdout.write("numeric  chol time ${toc (t)} ms\n");
  if (S == null || N == null) {
    return false;
  }
  t = tic();
  ipvec(S.pinv, b, y, n); // y = P*b
  lsolve(N.L, y); // y = L\y
  ltsolve(N.L, y); // y = L'\y
  pvec(S.pinv, y, x, n); // x = P'*y
  stdout.write("solve    chol time ${toc (t)} ms\n");
  stdout.write("original: ");
  printResid(true, C, x, b, resid, prob); // print residual
  k = n ~/ 2; // construct W
  W = spalloc(n, 1, n, true, false);
  Lp = N.L.p;
  Li = N.L.i;
  Lx = new Vector(N.L.x);
  Wp = W.p;
  Wi = W.i;
  Wx = new Vector(W.x);
  Wp[0] = 0;
  p1 = Lp[k];
  Wp[1] = Lp[k + 1] - p1;
  s = Lx.get(p1);
  math.Random r = new math.Random(1);
  for ( ; p1 < Lp[k + 1]; p1++) {
    p2 = p1 - Lp[k];
    Wi[p2] = Li[p1];
    Wx.setList(p2, cmult(s, r.nextDouble()));
  }
  t = tic();
  ok = updown(N.L, 1, W, S.parent); // update: L*L'+W*W'
  t1 = toc(t);
  stdout.write("update:   time: $t1 ms\n");
  if (!ok) return (false);
  t = tic();
  ipvec(S.pinv, b, y, n); // y = P*b
  lsolve(N.L, y); // y = L\y
  ltsolve(N.L, y); // y = L'\y
  pvec(S.pinv, y, x, n); // x = P'*y
  t = toc(t);
  p = pinv(S.pinv, n);
  W2 = permute(W, p, null, true); // E = C + (P'W)*(P'W)'
  WT = transpose(W2, true);
  WW = multiply(W2, WT);
  WT = null;
  W2 = null;
  E = add(C, WW, cone(), cone());
  WW = null;
  if (E == null || p == null) {
    return false;
  }
  stdout.write("update:   time: ${t1 + t} ms(incl solve) ");
  printResid(true, E, x, b, resid, prob); // print residual
  N = null; // clear N
  t = tic();
  N = chol(E, S); // numeric Cholesky
  if (N == null) {
    return false;
  }
  ipvec(S.pinv, b, y, n); // y = P*b
  lsolve(N.L, y); // y = L\y
  ltsolve(N.L, y); // y = L'\y
  pvec(S.pinv, y, x, n); // x = P'*y
  t = toc(t);
  stdout.write("rechol:   time: $t ms(incl solve) ");
  printResid(true, E, x, b, resid, prob); // print residual
  t = tic();
  ok = updown(N.L, -1, W, S.parent); // downdate: L*L'-W*W'
  t1 = toc(t);
  if (!ok) {
    return false;
  }
  stdout.write("downdate: time: $t1\n");
  t = tic();
  ipvec(S.pinv, b, y, n); // y = P*b
  lsolve(N.L, y); // y = L\y
  ltsolve(N.L, y); // y = L'\y
  pvec(S.pinv, y, x, n); // x = P'*y
  t = toc(t);
  stdout.write("downdate: time: ${t1 + t} ms(incl solve) ");
  printResid(true, C, x, b, resid, prob); // print residual
  return true;
}

/// Read a matrix, solve a linear system, update/downdate.
main() {
  test('c4', () {
    final file = getFile(C4);
    Problem prob = getProblem(file, 0.0);

    assertProblem(prob, 4, 4, 10, -1, 16, 7.37e+01);

    test3(prob);

    double x_norm = 8.3862;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    x_norm = 9.0064;
    expect(x_norm, closeTo(prob.norms[1], 0.3)); // TODO: tighten accuracy
    expect(x_norm, closeTo(prob.norms[2], 0.3));
    x_norm = 8.3862;
    expect(x_norm, closeTo(prob.norms[3], DELTA));
  });

  test('mhd1280b', () {
    final file = getFile(MHD1280B);
    Problem prob = getProblem(file, 0.0);

    assertProblem(prob, 1280, 1280, 12029, -1, 22778, 79.9740);

    test3(prob);

    double x_norm = 76270143066.4197;
    expect(x_norm, closeTo(prob.norms[0], DELTA));
    expect(x_norm, closeTo(prob.norms[1], DELTA));
    expect(x_norm, closeTo(prob.norms[2], DELTA));
    expect(x_norm, closeTo(prob.norms[3], DELTA));
  });
}
