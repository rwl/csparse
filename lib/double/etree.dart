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

/// Compute the elimination tree of A or A'A (without forming A'A).
///
/// Analyze A if [ata] is false, A'A otherwise.
/// Returns elimination tree, null on error.
Int32List cs_etree(Dcs A, bool ata) {
  int p, m, n, inext;
  Int32List Ap, Ai, w, parent, ancestor, prev;
  if (!cs_csc(A)) {
    return null; // check inputs
  }
  m = A.m;
  n = A.n;
  Ap = A.p;
  Ai = A.i;
  parent = new Int32List(n); // allocate result
  w = new Int32List(n + (ata ? m : 0)); // get workspace
  ancestor = w;
  prev = w;
  int prev_offset = n;
  if (ata) {
    for (int i = 0; i < m; i++) {
      prev[prev_offset + i] = -1;
    }
  }
  for (int k = 0; k < n; k++) {
    parent[k] = -1; // node k has no parent yet
    ancestor[k] = -1; // nor does k have an ancestor
    for (p = Ap[k]; p < Ap[k + 1]; p++) {
      int i = ata ? (prev[prev_offset + Ai[p]]) : (Ai[p]);
      for ( ; i != -1 && i < k; i = inext) // traverse from i to k
      {
        inext = ancestor[i]; // inext = ancestor of i
        ancestor[i] = k; // path compression
        if (inext == -1) {
          parent[i] = k; // no anc., parent is k
        }
      }
      if (ata) {
        prev[prev_offset + Ai[p]] = k;
      }
    }
  }
  return parent;
}
