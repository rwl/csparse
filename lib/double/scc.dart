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

/// Finds the strongly connected components of a square matrix.
///
/// [A] column-compressed matrix (A.p modified then restored).
/// Returns strongly connected components, null on error.
Dcsd cs_scc(Dcs A) {
  int n, top;
  int nb = 0;
  Int32List xi, pstack, p, r, Ap, ATp, rcopy, Blk;
  Dcs AT;
  Dcsd D;
  if (!CS_CSC(A)) {
    return null; // check inputs
  }
  n = A.n;
  Ap = A.p;
  D = cs_dalloc(n, 0); // allocate result
  AT = cs_transpose(A, false); // AT = A'
  xi = new Int32List(2 * n + 1); // get workspace
  if (D == null || AT == null) {
    return null;
  }
  Blk = xi;
  rcopy = xi;
  int rcopy_offset = n;
  pstack = xi;
  int pstack_offset = n;
  p = D.p;
  r = D.r;
  ATp = AT.p;
  top = n;
  for (int i = 0; i < n; i++) // first dfs(A) to find finish times (xi)
  {
    if (!CS_MARKED(Ap, i)) {
      top = cs_dfs(i, A, top, xi, 0, pstack, pstack_offset, null, 0);
    }
  }
  for (int i = 0; i < n; i++) {
    CS_MARK(Ap, i); // restore A; unmark all nodes
  }
  top = n;
  nb = n;
  for (int k = 0; k < n; k++) // dfs(A') to find strongly connnected comp
  {
    final i = xi[k]; // get i in reverse order of finish times
    if (CS_MARKED(ATp, i)) {
      continue; // skip node i if already ordered
    }
    r[nb--] = top; // node i is the start of a component in p
    top = cs_dfs(i, AT, top, p, 0, pstack, pstack_offset, null, 0);
  }
  r[nb] = 0; // first block starts at zero; shift r up
  for (int k = nb; k <= n; k++) {
    r[k - nb] = r[k];
  }
  D.nb = nb = n - nb; // nb = # of strongly connected components
  for (int b = 0; b < nb; b++) // sort each block in natural order
  {
    for (int k = r[b]; k < r[b + 1]; k++) {
      Blk[p[k]] = b;
    }
  }
  for (int b = 0; b <= nb; b++) {
    rcopy[rcopy_offset + b] = r[b];
  }
  for (int i = 0; i < n; i++) {
    p[rcopy[rcopy_offset + Blk[i]]++] = i;
  }
  return D;
}
