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
library edu.emory.mathcs.csparse.test_util;

import 'dart:typed_data';
import 'dart:math' as math;
import 'dart:io' show stdout, File, Platform;

import 'package:unittest/unittest.dart';

import 'csparse.dart';
import 'load.dart';

const double DELTA = 1e-3;
const double DROP_TOL = 1e-14;

final String DIR = "matrix";

final String T1 = "t1";

/// Unsymmetric overdetermined pattern of Holland survey. Ashkenazi, 1974
///
/// http://www.cise.ufl.edu/research/sparse/matrices/HB/ash219.html
final String ASH219 = "ash219";

/// Symmetric stiffness matrix small generalized eigenvalue problem.
///
/// http://www.cise.ufl.edu/research/sparse/matrices/HB/bcsstk01.html
final String BCSSTK01 = "bcsstk01";

/// S stiffness matrix - Corp. of Engineers Dam
///
/// http://www.cise.ufl.edu/research/sparse/matrices/HB/bcsstk16.html
final String BCSSTK16 = "bcsstk16";

/// Unsymmetric facsimile convergence matrix.
///
/// http://www.cise.ufl.edu/research/sparse/matrices/HB/fs_183_1.html
final String FS_183_1 = "fs_183_1";

/// Unsymmetric pattern on leaflet advertising ibm 1971 conference,
/// but with the last column removed.
///
/// http://www.cise.ufl.edu/research/sparse/matrices/HB/ibm32.html
final String IBM32A = "ibm32a";

/// The transpose of [IBM32A].
final String IBM32B = "ibm32b";

/// Netlib LP problem afiro: minimize c'*x, where Ax=b, lo<=x<=hi
///
/// http://www.cise.ufl.edu/research/sparse/matrices/LPnetlib/lp_afiro.html
final String LP_AFIRO = "lp_afiro";

/// U Nonsymmetric matrix U.S. Economy 1972 -SZYLD-I.E.A.-NYU-
///
/// http://www.cise.ufl.edu/research/sparse/matrices/HB/mbeacxc.html
final String MBEACXC = "mbeacxc";

/// Cavett problem with 5 components (chemical eng., Westerberg)
///
/// http://www.cise.ufl.edu/research/sparse/matrices/HB/west0067.html
final String WEST0067 = "west0067";


File get_file(String name) {
  return new File([Uri.base.toFilePath() + DIR, name].join('/'));
}

void assert_dimensions(Dcs A, int m, int n, int nzmax, int nnz, [double norm1 = null, double delta = DELTA]) {
  expect(m, equals(A.m));
  expect(n, equals(A.n));
  expect(nzmax, equals(A.nzmax));

  int nz = (A.nz < 0) ? A.p[A.n] : A.nz;
  expect(nnz, equals(nz));

  if (norm1 != null) {
    expect(norm1, closeTo(cs_norm(A), delta));
  }
}

void assert_problem(Dproblem prob, int m, int n, int nnz, int sym, int sym_nnz, double norm) {
  expect(m, equals(prob.A.m));
  expect(n, equals(prob.A.n));
  expect(nnz, equals(prob.A.p[n]));
  expect(sym, equals(prob.sym));
  expect(sym_nnz, equals(sym != 0 ? prob.C.p[n] : 0));
  expect(norm, closeTo(cs_norm(prob.C), 1e-2));
}

void assert_structure(Dproblem prob, int blocks, int singletons, int rank) {
  expect(blocks, equals(prob.nb));
  expect(singletons, equals(prob.ns));
  expect(rank, equals(prob.sprank));
}

void assert_dropped(Dproblem prob, int dropped_zeros, int dropped_tiny) {
  expect(dropped_zeros, equals(prob.dropped_zeros));
  expect(dropped_tiny, equals(prob.dropped_tiny));
}

/// A structure for a demo problem.
class Dproblem {
  Dcs A;
  Dcs C;
  int sym;
  Float64List x;
  Float64List b;
  Float64List resid;

  final norms = new List<double>();

  int nb;
  int ns;
  int sprank;

  int dropped_zeros;
  int dropped_tiny;

  Dproblem();
}

/// 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise
int is_sym(Dcs A) {
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
  return is_upper ? 1 : (is_lower ? -1 : 0);
}

/// true for off-diagonal entries
bool dropdiag(int i, int j, double aij, Object other) => i != j;

/// C = A + triu(A,1)'
Dcs make_sym(Dcs A) {
  Dcs AT, C;
  AT = cs_transpose(A, true);
  // AT = A'
  cs_fkeep(AT, dropdiag, null);
  // drop diagonal entries from AT
  C = cs_add(A, AT, 1.0, 1.0);
  // C = A+AT
  AT = null;
  return C;
}

/// create a right-hand side
void rhs(Float64List x, Float64List b, int m) {
  for (int i = 0; i < m; i++) b[i] = 1 + i.toDouble() / m;
  for (int i = 0; i < m; i++) x[i] = b[i];
}

/// infinity-norm of x
double norm(Float64List x, int n) {
  double normx = 0.0;
  for (int i = 0; i < n; i++) {
    normx = math.max(normx, x[i].abs());
  }
  return normx;
}

/// compute residual, norm(A*x-b,inf) / (norm(A,1)*norm(x,inf) + norm(b,inf))
void print_resid(bool ok, Dcs A, Float64List x, Float64List b, Float64List resid, Dproblem prob) {
  int m, n;
  if (!ok) {
    stdout.write("    (failed)\n");
    return;
  }
  m = A.m;
  n = A.n;
  for (int i = 0; i < m; i++) resid[i] = -b[i];
  // resid = -b
  cs_gaxpy(A, x, resid);
  // resid = resid + A*x

  double r = norm(resid, m) / ((n == 0) ? 1 : (cs_norm(A) * norm(x, n) + norm(b, m)));
  stdout.write("resid: $r");

  double nrm = norm(x, n);
  stdout.write(" (norm: $nrm, $norm (b, m))\n");
  prob.norms.add(nrm);
}

int tic() => new DateTime.now().millisecondsSinceEpoch;

int toc(int t) => math.max(0, tic() - t);// / 1000000.0 ;

void print_order(int order) {
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
///
/// [tol] drop tolerance.
/// [base] file index base.
Dproblem get_problem(File file, [double tol = 0.0, int base = 0]) {
  Dcs T, A, C;
  int sym, m, n, mn, nz1, nz2;
  Dproblem prob;
  prob = new Dproblem();
  T = cs_load(file, base);
  // load triplet matrix T from a file
  prob.A = A = cs_compress(T);
  // A = compressed-column form of T
  T = null;
  /* clear T */
  if (!cs_dupl(A)) return (null);
  // sum up duplicates
  prob.sym = sym = is_sym(A);
  // determine if A is symmetric
  m = A.m;
  n = A.n;
  mn = math.max(m, n);
  nz1 = A.p[n];
  if (tol > 0) cs_dropzeros(A);
  // drop zero entries
  nz2 = A.p[n];
  if (tol > 0) cs_droptol(A, tol);
  // drop tiny entries (just to test)
  prob.C = C = sym != 0 ? make_sym(A) : A;
  // C = A + triu(A,1)', or C=A
  if (C == null) return null;
  stdout.write("\n--- Matrix: $m-by-$n, nnz: ${A.p[n]} (sym: $sym: nnz ${sym != 0 ? C.p[n] : 0}), norm: ${cs_norm(C)}\n");
  prob.dropped_zeros = nz1 - nz2;
  if (nz1 != nz2) {
    stdout.write("zero entries dropped: ${nz1 - nz2}\n");
  }
  prob.dropped_tiny = nz2 - A.p[n];
  if (nz2 != A.p[n]) {
    stdout.write("tiny entries dropped: ${nz2 - A.p[n]}\n");
}
  prob.b = new Float64List(mn);
  prob.x = new Float64List(mn);
  prob.resid = new Float64List(mn);
  return prob;
}
