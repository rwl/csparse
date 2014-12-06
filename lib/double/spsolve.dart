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
part of edu.emory.mathcs.csparse;

/// Sparse lower or upper triangular solve. x=G\b where G, x, and b are sparse,
/// and G upper/lower triangular.
///
/// Solve Gx=b(:,k), where G is either upper (lo=false) or lower (lo=true)
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
int spsolve(Matrix G, Matrix B, int k, Int32List xi, Float64List x, Int32List pinv, bool lo) {
  int top, n;
  Int32List Gp, Gi, Bp, Bi;
  Float64List Gx, Bx;
  if (!csc(G) || !csc(B) || xi == null || x == null) {
    return -1;
  }
  Gp = G.p;
  Gi = G.i;
  Gx = G.x;
  n = G.n;
  Bp = B.p;
  Bi = B.i;
  Bx = B.x;
  top = reach(G, B, k, xi, pinv); // xi[top..n-1]=Reach(B(:,k))
  for (int p = top; p < n; p++) {
    x[xi[p]] = 0.0; // clear x
  }
  for (int p = Bp[k]; p < Bp[k + 1]; p++) {
    x[Bi[p]] = Bx[p]; // scatter B
  }
  for (int px = top; px < n; px++) {
    final j = xi[px]; // x(j) is nonzero
    final J = pinv != null ? (pinv[j]) : j; // j maps to col J of G
    if (J < 0) {
      continue; // column J is empty
    }
    x[j] /= Gx[lo ? (Gp[J]) : (Gp[J + 1] - 1)]; // x(j) /= G(j,j)
    int p = lo ? (Gp[J] + 1) : (Gp[J]); // lo: L(j,j) 1st entry
    final q = lo ? (Gp[J + 1]) : (Gp[J + 1] - 1); // up: U(j,j) last entry
    for ( ; p < q; p++) {
      x[Gi[p]] -= Gx[p] * x[j]; // x(i) -= G(i,j) * x(j)
    }
  }
  return top; // return top of stack
}
