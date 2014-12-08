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

/// Solve Gx=b(:,k), where G is either upper ([lo]=false) or lower ([lo]=true)
/// triangular.
///
/// [G] lower or upper triangular matrix in column-compressed form.
/// [B] right hand side, b=B(:,k).
/// [k] use kth column of B as right hand side.
/// [xi] size 2*n, nonzero pattern of x in xi[top..n-1].
/// [x] size n, x in x[xi[top..n-1]].
/// [pinv] mapping of rows to columns of G, ignored if null.
/// [lo] true if lower triangular, false if upper.
/// Returns top, -1 in error.
int spsolve(Matrix G, Matrix B, int k, Int32List xi, Vector x, Int32List pinv, bool lo) {
  int j, J, p, q, px, top, n;
  Int32List Gp, Gi, Bp, Bi;
  Vector Gx = new Vector(),
      Bx = new Vector();
  if (!csc(G) || !csc(B) || xi == null || x == null) {
    return -1;
  }
  Gp = G.p;
  Gi = G.i;
  Gx.x = G.x;
  n = G.n;
  Bp = B.p;
  Bi = B.i;
  Bx.x = B.x;
  top = reach(G, B, k, xi, pinv); // xi[top..n-1]=Reach(B(:,k))
  for (p = top; p < n; p++) {
    x.setList(xi[p], czero()); // clear x
  }
  for (p = Bp[k]; p < Bp[k + 1]; p++) {
    x.setList(Bi[p], Bx.get(p)); // scatter B
  }
  for (px = top; px < n; px++) {
    j = xi[px]; // x(j) is nonzero
    J = pinv != null ? (pinv[j]) : j; // j maps to col J of G
    if (J < 0) {
      continue; // column J is empty
    }
    x.setList(j, cdiv_list(x.get(j), Gx.get(lo ? (Gp[J]) : (Gp[J + 1] - 1)))); // x(j) /= G(j,j)
    p = lo ? (Gp[J] + 1) : (Gp[J]); // lo: L(j,j) 1st entry
    q = lo ? (Gp[J + 1]) : (Gp[J + 1] - 1); // up: U(j,j) last entry
    for ( ; p < q; p++) {
      x.setList(Gi[p], cminus(x.get(Gi[p]), cmult_list(Gx.get(p), x.get(j)))); // x(i) -= G(i,j) * x(j)
    }
  }
  return top; // return top of stack
}
