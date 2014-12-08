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
library cxsparse.test.test1;

import 'dart:io';
import 'dart:typed_data';

import 'package:unittest/unittest.dart';
import 'package:csparse/complex/cxsparse.dart';
import 'package:csparse/complex/load.dart';
import 'package:csparse/complex/test_util.dart';

Matrix _load(File file) {
  Matrix T = load(file); // load triplet matrix T from file
  //System.out.print("T:\n") ;
  //cs_print (T, false) ;
  return T;
}

Matrix _compress(Matrix T) {
  Matrix A = compress(T); // A = compressed-column form of T
  //System.out.print("A:\n") ;
  //cs_print (A, false) ;
  return A;
}

Matrix _transpose(Matrix A) {
  Matrix AT = transpose(A, true); // AT = A'
  //System.out.print("AT:\n") ;
  //cs_print (AT, false) ;
  return AT;
}

Matrix _multiplyAdd(Matrix A, Matrix AT) {
  Matrix T, Eye, C, D;
  int m = A != null ? A.m : 0; // m = # of rows of A
  T = spalloc(m, m, m, true, true); // create triplet identity matrix
  for (int i = 0; i < m; i++) {
    entry(T, i, i, cone());
  }
  Eye = compress(T); // Eye = speye (m)
  T = null;
  C = multiply(A, AT); // C = A*A'
  D = add(C, Eye, cone(), new Float64List.fromList([norm(C), 0.0])); // D = C + Eye*norm (C,1)
  //System.out.print("D:\n") ;
  //cs_print(D, false) ;
  return D;
}

/// Read a matrix from a file and perform basic matrix operations.
main() {
  test('demo1_c_ibm32a', () {
    Matrix T, A, AT, D;

    final file = getFile(C_IBM32A);

    T = _load(file);
    assertDimensions(T, 32, 31, 128, 123);

    A = _compress(T);
    assertDimensions(A, 32, 31, 123, 123, 9.89949);

    AT = _transpose(A);
    assertDimensions(AT, 31, 32, 123, 123, 11.3137);

    D = _multiplyAdd(A, AT);
    assertDimensions(D, 32, 32, 386, 386, 140.0);
  });

  test('demo1_c_ibm32b', () {
    Matrix T, A, AT, D;

    final file = getFile(C_IBM32B);

    T = _load(file);
    assertDimensions(T, 31, 32, 128, 123);

    A = _compress(T);
    assertDimensions(A, 31, 32, 123, 123, 11.3137);

    AT = _transpose(A);
    assertDimensions(AT, 32, 31, 123, 123, 9.89949);

    D = _multiplyAdd(A, AT);
    assertDimensions(D, 31, 31, 373, 373, 128.0);
  });

  test('demo1_c_mbeacxc', () {
    Matrix T, A, AT, D;

    final file = getFile(C_MBEACXC);

    T = _load(file);
    assertDimensions(T, 492, 490, 65536, 49920);

    A = _compress(T);
    assertDimensions(A, 492, 490, 49920, 49920, 0.92863);

    AT = _transpose(A);
    assertDimensions(AT, 490, 492, 49920, 49920, 16.5516);

    D = _multiplyAdd(A, AT);
    assertDimensions(D, 492, 492, 157350, 157350, 19.6068);
  });

  test('demo1_c_west0067', () {
    Matrix T, A, AT, D;

    final file = getFile(C_WEST0067);

    T = _load(file);
    assertDimensions(T, 67, 67, 512, 299);

    A = _compress(T);
    assertDimensions(A, 67, 67, 299, 299, 6.16948);

    AT = _transpose(A);
    assertDimensions(AT, 67, 67, 299, 299, 6.62541);

    D = _multiplyAdd(A, AT);
    assertDimensions(D, 67, 67, 1041, 1041, 61.4518);
  });

  test('demo1_c4', () {
    Matrix T, A, AT, D;

    final file = getFile(C4);

    T = _load(file);
    assertDimensions(T, 4, 4, 16, 10);

    A = _compress(T);
    assertDimensions(A, 4, 4, 10, 10, 58.8904);

    AT = _transpose(A);
    assertDimensions(AT, 4, 4, 10, 10, 66.8771);

    D = _multiplyAdd(A, AT);
    assertDimensions(D, 4, 4, 16, 16, 5029.4944);
  });

  test('demo1_czero', () {
    Matrix T, A, AT, D;

    final file = getFile(CZERO);

    T = _load(file);
    assertDimensions(T, 1, 1, 1, 1);

    A = _compress(T);
    assertDimensions(A, 1, 1, 1, 1, 0.0);

    AT = _transpose(A);
    assertDimensions(AT, 1, 1, 1, 1, 0.0);

    D = _multiplyAdd(A, AT);
    assertDimensions(D, 1, 1, 1, 1, 0.0);
  });

  test('demo1_mhd1280b', () {
    Matrix T, A, AT, D;

    final file = getFile(MHD1280B);

    T = _load(file);
    assertDimensions(T, 1280, 1280, 16384, 12029);

    A = _compress(T);
    assertDimensions(A, 1280, 1280, 12029, 12029, 79.974);

    AT = _transpose(A);
    assertDimensions(AT, 1280, 1280, 12029, 12029, 64.1995);

    D = _multiplyAdd(A, AT);
    assertDimensions(D, 1280, 1280, 27642, 27642, 8519.6629);
  });

  test('demo1_qc324', () {
    Matrix T, A, AT, D;

    final file = getFile(QC324);

    T = _load(file);
    assertDimensions(T, 324, 324, 32768, 26730);

    A = _compress(T);
    assertDimensions(A, 324, 324, 26730, 26730, 1.70664);

    AT = _transpose(A);
    assertDimensions(AT, 324, 324, 26730, 26730, 1.70664);

    D = _multiplyAdd(A, AT);
    assertDimensions(D, 324, 324, 65934, 65934, 5.42006);
  });

  test('demo1_t2', () {
    Matrix T, A, AT, D;

    final file = getFile(T2);

    T = _load(file);
    assertDimensions(T, 4, 4, 16, 10);

    A = _compress(T);
    assertDimensions(A, 4, 4, 10, 10, 106.07516);

    AT = _transpose(A);
    assertDimensions(AT, 4, 4, 10, 10, 144.29640);

    D = _multiplyAdd(A, AT);
    assertDimensions(D, 4, 4, 16, 16, 25308.3283);
  });

  test('demo1_t3', () {
    Matrix T, A, AT, D;

    final file = getFile(T3);

    T = _load(file);
    assertDimensions(T, 3, 4, 16, 12);

    A = _compress(T);
    assertDimensions(A, 3, 4, 12, 12, 3.05601);

    AT = _transpose(A);
    assertDimensions(AT, 4, 3, 12, 12, 3.93809);

    D = _multiplyAdd(A, AT);
    assertDimensions(D, 3, 3, 9, 9, 21.3318);
  });

  test('demo1_t4', () {
    Matrix T, A, AT, D;

    final file = getFile(T4);

    T = _load(file);
    assertDimensions(T, 2, 2, 4, 3);

    A = _compress(T);
    assertDimensions(A, 2, 2, 3, 3, 2.82843);

    AT = _transpose(A);
    assertDimensions(AT, 2, 2, 3, 3, 2.82843);

    D = _multiplyAdd(A, AT);
    assertDimensions(D, 2, 2, 4, 4, 12.0);
  });

  test('demo1_young1c', () {
    Matrix T, A, AT, D;

    final file = getFile(YOUNG1C);

    T = _load(file);
    assertDimensions(T, 841, 841, 4096, 4089);

    A = _compress(T);
    assertDimensions(A, 841, 841, 4089, 4089, 730.46);

    AT = _transpose(A);
    assertDimensions(AT, 841, 841, 4089, 4089, 730.46);

    D = _multiplyAdd(A, AT);
    assertDimensions(D, 841, 841, 10357, 10357, 1.0671436232e+06);
  });
}
