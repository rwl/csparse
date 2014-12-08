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
library edu.emory.mathcs.cxsparse.test_util;

import 'dart:io';
import 'dart:typed_data';
import 'dart:math' as math;

import 'package:unittest/unittest.dart';

import 'cxsparse.dart';
import 'load.dart';

final double DELTA = 1e-4;
final double DROP_TOL = 1e-14;

const String DIR = "matrix";

final String CZERO = "czero";
final String C4 = "c4";
final String T2 = "t2";
final String T3 = "t3";
final String T4 = "t4";

/// US economy, 1972.  Dan Szyld, while at NYU.
/// Numerically and structurally singular.
final String C_MBEACXC = "c_mbeacxc";

/// Cavett problem with 5 components (chemical eng., Westerberg).
final String C_WEST0067 = "c_west0067";

/// Alfven spectra in magnetohydrodynamics.
final String MHD1280B = "mhd1280b";

/// Model of H+ in an electromagnetic field.
final String QC324 = "qc324";

/// The Harwell/Boeing matrix ibm32, but with the last column removed.
final String C_IBM32A = "c_ibm32a";

/// The transpose of [C_IBM32A].
final String C_IBM32B = "c_ibm32b";

/// Aeronautical problem.
final String YOUNG1C = "young1c";


File getFile(String name, [String dir = DIR]) {
  return new File([Uri.base.toFilePath() + dir, name].join('/'));
}

void assertDimensions(Matrix A, int m, int n, int nzmax, int nnz, [double norm1 = null]) {
  expect(m, equals(A.m));
  expect(n, equals(A.n));
  expect(nzmax, equals(A.nzmax));

  int nz = (A.nz < 0) ? A.p[A.n] : A.nz;
  expect(nnz, equals(nz));

  if (norm1 != null) {
    expect(norm1, closeTo(norm(A), DELTA));
  }
}

void assertProblem(Problem prob, int m, int n, int nnz, int sym, int sym_nnz, double _norm) {
  expect(m, equals(prob.A.m));
  expect(n, equals(prob.A.n));
  expect(nnz, equals(prob.A.p[n]));
  expect(sym, equals(prob.sym));
  expect(sym_nnz, equals(sym != 0 ? prob.C.p[n] : 0));
  expect(_norm, closeTo(norm(prob.C), 1e-2));
}

void assertStructure(Problem prob, int blocks, int singletons, int rank) {
  expect(blocks, equals(prob.nb));
  expect(singletons, equals(prob.ns));
  expect(rank, equals(prob.sprank));
}

void assertDropped(Problem prob, int dropped_zeros, int dropped_tiny) {
  expect(dropped_zeros, equals(prob.dropped_zeros));
  expect(dropped_tiny, equals(prob.dropped_tiny));
}

/// A structure for a demo problem.
class Problem {
  Matrix A;
  Matrix C;
  int sym;
  Vector x;
  Vector b;
  Vector resid;

  List<double> norms = new List<double>();

  int nb;
  int ns;
  int sprank;

  int dropped_zeros;
  int dropped_tiny;

  Problem();
}

/// 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise
int isSym(Matrix A) {
  int n = A.n,
      m = A.m;
  Int32List Ap = A.p,
      Ai = A.i;
  bool is_upper, is_lower;
  if (m != n) return (0);
  is_upper = true;
  is_lower = true;
  for (int j = 0; j < n; j++) {
    for (int p = Ap[j]; p < Ap[j + 1]; p++) {
      if (Ai[p] > j) is_upper = false;
      if (Ai[p] < j) is_lower = false;
    }
  }
  return (is_upper ? 1 : (is_lower ? -1 : 0));
}

/// true for off-diagonal entries
bool dropdiag(int i, int j, aij, other) => i != j;

/// C = A + triu(A,1)'
Matrix makeSym(Matrix A) {
  Matrix AT, C;
  AT = transpose(A, true); // AT = A'
  fkeep(AT, dropdiag, null); // drop diagonal entries from AT
  C = add(A, AT, cone(), cone()); // C = A+AT
  AT = null;
  return C;
}

/// Create a right-hand side.
void rhs(Vector x, Vector b, int m) {
  for (int i = 0; i < m; i++) {
    b.setList(i, new Float64List.fromList([1 + (i.toDouble()) / m, 0.0]));
  }
  for (int i = 0; i < m; i++) {
    x.setList(i, b.get(i));
  }
}

/// Infinity-norm of x.
double normInf(Vector x, int n) {
  double normx = 0.0;
  for (int i = 0; i < n; i++) {
    normx = math.max(normx, cabs_list(x.get(i)));
  }
  return normx;
}

/// Compute residual, norm(A*x-b,inf) / (norm(A,1)*norm(x,inf) + norm(b,inf)).
void printResid(bool ok, Matrix A, Vector x, Vector b, Vector resid, Problem prob) {
  int m = A.m,
      n = A.n;
  if (!ok) {
    stdout.write("    (failed)\n");
    return;
  }
  for (int i = 0; i < m; i++) {
    resid.setList(i, cneg(b.get(i))); // resid = -b
  }
  gaxpy(A, x, resid); // resid = resid + A*x

  double r = normInf(resid, m) / ((n == 0) ? 1 : (norm(A) * normInf(x, n) + normInf(b, m)));
  stdout.write("resid: $r");

  double nrm = normInf(x, n);
  stdout.write(" (norm: $nrm, ${normInf(b, m)})\n");
  prob.norms.add(nrm);
}

int tic() => new DateTime.now().millisecondsSinceEpoch;

int toc(int t) => math.max(0, tic() - t);// / 1000000.0 ;

void printOrder(int order) {
  switch (order) {
    case 0:
      stdout.write("natural    ");
      break;
    case 1:
      stdout.write("amd(A+A')  ");
      break;
    case 2:
      stdout.write("amd(S'*S)  ");
      break;
    case 3:
      stdout.write("amd(A'*A)  ");
      break;
  }
}

/// Reads a problem from a file.
Problem getProblem(File file, double tol) {
  Matrix T, A, C;
  int sym, m, n, mn, nz1, nz2;
  Problem prob;
  prob = new Problem();
  T = load(file); // load triplet matrix T from a file
  prob.A = A = compress(T); // A = compressed-column form of T
  T = null; // clear T
  if (!dupl(A)) {
    return null; // sum up duplicates
  }
  prob.sym = sym = isSym(A); // determine if A is symmetric
  m = A.m;
  n = A.n;
  mn = math.max(m, n);
  nz1 = A.p[n];
  if (tol > 0) {
    dropzeros(A); // drop zero entries
  }
  nz2 = A.p[n];
  if (tol > 0) {
    droptol(A, tol); // drop tiny entries (just to test)
  }
  prob.C = C = sym != 0 ? makeSym(A) : A; // C = A + triu(A,1)', or C=A
  if (C == null) {
    return null;
  }
  stdout.write("\n--- Matrix: $m-by-$n, nnz: ${A.p [n]} (sym: $sym: nnz ${sym != 0 ? C.p [n] : 0}), norm: ${norm (C)}\n");
  prob.dropped_zeros = nz1 - nz2;
  if (nz1 != nz2) {
    stdout.write("zero entries dropped: ${nz1 - nz2}\n");
  }
  prob.dropped_tiny = nz2 - A.p[n];
  if (nz2 != A.p[n]) {
    stdout.write("tiny entries dropped: ${nz2 - A.p [n]}\n");
  }
  prob.b = new Vector.sized(mn);
  prob.x = new Vector.sized(mn);
  prob.resid = new Vector.sized(mn);
  return prob;
}
