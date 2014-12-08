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
part of edu.emory.mathcs.cxsparse;

/// Compute a Householder reflection [v,beta,s]=house(x), overwrite x with v,
/// where (I-beta*v*v')*x = s*e1 and e1 = [1 0 ... 0]'.
///
/// Note that this CXSparse version is different than CSparse. See Higham,
/// Accuracy & Stability of Num Algorithms, 2nd ed, 2002, page 357.
///
/// [x] x on output, v on input.
/// [n] the length of x.
/// Returns norm2(x), -1 on error.
Float64List house(Vector x, int x_offset, Float64List beta, int n) {
  Float64List s = czero();
  int i;
  if (x == null) {
    return new Float64List.fromList([-1.0, 0.0]);
  }
  /* s = norm(x) */
  for (i = 0; i < n; i++) { // TODO: check i = 1
    s = cplus(s, cmult_list(x.get(x_offset + i), conj(x.get(x_offset + i))));
  }
  s = csqrt(s);
  if (cequal(s, czero())) {
    beta[0] = 0.0;
    x.setParts(x_offset + 0, 1.0, 0.0);
  } else {
    /* s = sign(x[0]) * norm (x) ; */
    if (!cequal(x.get(x_offset + 0), czero())) {
      s = cmult_list(s, cdiv_list(x.get(x_offset + 0),
          new Float64List.fromList([cabs_list(x.get(x_offset + 0)), 0.0])));
    }
    x.setList(x_offset + 0, cplus(x.get(x_offset + 0), s));
    beta[0] = 1 / _creal(cmult_list(conj(s), x.get(x_offset + 0)));
  }
  return cmult(s, -1.0);
}
