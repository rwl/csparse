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

/// find an augmenting path starting at column k and extend the match if found
void _augment(int k, Matrix A, Int32List jmatch, int jmatch_offset, Int32List cheap, int cheap_offset,
                Int32List w, int w_offset, Int32List js, int js_offset, Int32List _is, int is_offset,
                Int32List ps, int ps_offset) {
  int p, j;
  Int32List Ap = A.p,
      Ai = A.i;
  int head = 0,
      i = -1;
  bool found = false;
  js[js_offset + 0] = k; // start with just node k in jstack
  while (head >= 0) {
    /* Start (or continue) depth-first-search at node j */
    j = js[js_offset + head]; // get j from top of jstack
    if (w[w_offset + j] != k) // 1st time j visited for kth path
    {
      w[w_offset + j] = k; // mark j as visited for kth path
      for (p = cheap[cheap_offset + j]; p < Ap[j + 1] && !found; p++) {
        i = Ai[p]; // try a cheap assignment (i,j)
        found = (jmatch[jmatch_offset + i] == -1);
      }
      cheap[cheap_offset + j] = p; // start here next time j is traversed
      if (found) {
        _is[is_offset + head] = i; // column j matched with row i
        break; // end of augmenting path
      }
      ps[ps_offset + head] = Ap[j]; // no cheap match: start dfs for j
    }
    /* Depth-first-search of neighbors of j */
    for (p = ps[ps_offset + head]; p < Ap[j + 1]; p++) {
      i = Ai[p]; // consider row i
      if (w[w_offset + jmatch[jmatch_offset + i]] == k) {
        continue; // skip jmatch [i] if marked
      }
      ps[ps_offset + head] = p + 1; // pause dfs of node j
      _is[is_offset + head] = i; // i will be matched with j if found
      js[js_offset + (++head)] = jmatch[jmatch_offset + i]; // start dfs at column jmatch [i]
      break;
    }
    if (p == Ap[j + 1]) {
      head--; // node j is done; pop from stack
    }
  } // augment the match if path found:
  if (found) {
    for (p = head; p >= 0; p--) {
      jmatch[jmatch_offset + _is[is_offset + p]] = js[js_offset + p];
    }
  }
}

/// Find a maximum transveral (permutation for zero-free diagonal). Seed optionally selects a
/// randomized algorithm.
///
/// [seed] 0: natural, -1: reverse, randomized otherwise.
/// Returns row and column matching, size m+n.
Int32List maxtrans(Matrix A, int seed) // [jmatch [0..m-1]; imatch [0..n-1]]
{
  int i, j, k, n, m, p;
  int n2 = 0,
      m2 = 0;
  Int32List Ap, jimatch, w, cheap, js, _is, ps, Ai, Cp, jmatch, imatch, q;
  Matrix C;
  if (!csc(A)) {
    return null; // check inputs
  }
  n = A.n;
  m = A.m;
  Ap = A.p;
  Ai = A.i;
  w = jimatch = new Int32List(m + n); // allocate result

  k = 0;
  for (j = 0; j < n; j++) // count nonempty rows and columns
  {
    if (Ap[j] < Ap[j + 1]) {
      n2++;
    }
    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      w[Ai[p]] = 1;
      if (j == Ai[p]) {
        k++; // count entries already on diagonal
      }
    }
  }
  if (k == math.min(m, n)) // quick return if diagonal zero-free
  {
    jmatch = jimatch;
    imatch = jimatch;
    int imatch_offset = m;
    for (i = 0; i < k; i++) {
      jmatch[i] = i;
    }
    for ( ; i < m; i++) {
      jmatch[i] = -1;
    }
    for (j = 0; j < k; j++) {
      imatch[imatch_offset + j] = j;
    }
    for ( ; j < n; j++) {
      imatch[imatch_offset + j] = -1;
    }
    return jimatch;
  }
  for (i = 0; i < m; i++) {
    m2 += w[i];
  }
  C = (m2 < n2) ? transpose(A, false) : A; // transpose if needed
  if (C == null) {
    return null;
  }
  n = C.n;
  m = C.m;
  Cp = C.p;
  jmatch = jimatch;
  imatch = jimatch;
  int jmatch_offset = 0;
  int imatch_offset = 0;
  if (m2 < n2) {
    jmatch_offset = n;
  } else {
    imatch_offset = m;
  }
  w = new Int32List(5 * n); // get workspace
  cheap = w;
  int cheap_offset = n;
  js = w;
  int js_offset = 2 * n;
  _is = w;
  int is_offset = 3 * n;
  ps = w;
  int ps_offset = 4 * n;
  for (j = 0; j < n; j++) {
    cheap[cheap_offset + j] = Cp[j]; // for cheap assignment
  }
  for (j = 0; j < n; j++) {
    w[j] = -1; // all columns unflagged
  }
  for (i = 0; i < m; i++) {
    jmatch[jmatch_offset + i] = -1; // nothing matched yet
  }
  q = randperm(n, seed); // q = random permutation
  for (k = 0; k < n; k++) // augment, starting at column q[k]
  {
    _augment(q != null ? q[k] : k, C, jmatch, jmatch_offset, cheap, cheap_offset, w, 0, js, js_offset, _is, is_offset, ps, ps_offset);
  }
  q = null;
  for (j = 0; j < n; j++) {
    imatch[imatch_offset + j] = -1; // find row match
  }
  for (i = 0; i < m; i++) {
    if (jmatch[jmatch_offset + i] >= 0) {
      imatch[imatch_offset + jmatch[jmatch_offset + i]] = i;
    }
  }
  return jimatch;
}
