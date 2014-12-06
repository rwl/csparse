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

/// Find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k)).
///
/// If ok, s[top..n-1] contains pattern of L(k,:).
///
/// [A] column-compressed matrix; L is the Cholesky factor of A.
/// [k] find kth row of L. [parent] elimination tree of A.
/// [s] size n, s[top..n-1] is nonzero pattern of L(k,1:k-1).
/// [s_offset] the index of the first element in array s.
/// [w] size n, work array, w[0..n-1]>=0 on input, unchanged on output.
/// Returns top in successful, -1 on error.
int cs_ereach(Dcs A, int k, Int32List parent, Int32List s, int s_offset, Int32List w) {
  int i, p, n, len, top;
  Int32List Ap, Ai;
  if (!cs_csc(A) || parent == null || s == null || w == null) {
    return (-1); // check inputs
  }
  top = n = A.n;
  Ap = A.p;
  Ai = A.i;
  cs_mark(w, k); // mark node k as visited
  for (p = Ap[k]; p < Ap[k + 1]; p++) {
    i = Ai[p]; // A(i,k) is nonzero
    if (i > k) {
      continue; // only use upper triangular part of A
    }
    for (len = 0; !cs_marked(w, i); i = parent[i]) // traverse up etree
    {
      s[s_offset + len++] = i; // L(k,i) is nonzero
      cs_mark(w, i); // mark i as visited
    }
    while (len > 0) {
      s[s_offset + --top] = s[s_offset + --len]; // push path onto stack
    }
  }
  for (p = top; p < n; p++) {
    cs_mark(w, s[s_offset + p]); // unmark all nodes
  }
  cs_mark(w, k); // unmark node k
  return (top); // s [top..n-1] contains pattern of L(k,:)
}
