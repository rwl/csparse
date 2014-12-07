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

/// Finds a nonzero pattern of x=L\b for sparse L and b.
///
/// [G] graph to search (G.p modified, then restored).
/// [B] right hand side, b = B(:,k).
/// [k] use kth column of B.
/// [xi] size 2*n, output in xi[top..n-1].
/// [pinv] mapping of rows to columns of G, ignored if null.
/// Returns top, -1 on error.
int cs_reach(DZcs G, DZcs B, int k, Int32List xi, Int32List pinv) {
  int p, n, top;
  Int32List Bp, Bi, Gp;
  if (!CS_CSC(G) || !CS_CSC(B) || xi == null) {
    return -1;
  }
  n = G.n;
  Bp = B.p;
  Bi = B.i;
  Gp = G.p;
  top = n;
  for (p = Bp[k]; p < Bp[k + 1]; p++) {
    if (!_CS_MARKED(Gp, Bi[p])) // start a dfs at unmarked node i
    {
      top = cs_dfs(Bi[p], G, top, xi, 0, xi, n, pinv, 0);
    }
  }
  for (p = top; p < n; p++) {
    _CS_MARK(Gp, xi[p]); // restore G
  }
  return top;
}
