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

/// Breadth-first search for coarse decomposition (C0,C1,R1 or R0,R3,C3)
bool _bfs(Matrix A, int n, Int32List wi, Int32List wj, Int32List queue,
             Int32List imatch, int imatch_offset, Int32List jmatch, int jmatch_offset, int mark) {
  Int32List Ap, Ai;
  int head = 0,
      tail = 0;
  Matrix C;
  for (int j = 0; j < n; j++) // place all unmatched nodes in queue
  {
    if (imatch[imatch_offset + j] >= 0) {
      continue; // skip j if matched
    }
    wj[j] = 0; // j in set C0 (R0 if transpose)
    queue[tail++] = j; // place unmatched col j in queue
  }
  if (tail == 0) {
    return true; // quick return if no unmatched nodes
  }
  C = (mark == 1) ? A : transpose(A, false);
  if (C == null) {
    return false; // bfs of C=A' to find R3,C3 from R0
  }
  Ap = C.p;
  Ai = C.i;
  while (head < tail) // while queue is not empty
  {
    final j = queue[head++]; // get the head of the queue
    for (int p = Ap[j]; p < Ap[j + 1]; p++) {
      final i = Ai[p];
      if (wi[i] >= 0) {
        continue; // skip if i is marked
      }
      wi[i] = mark; // i in set R1 (C3 if transpose)
      final j2 = jmatch[jmatch_offset + i]; // traverse alternating path to j2
      if (wj[j2] >= 0) {
        continue; // skip j2 if it is marked
      }
      wj[j2] = mark; // j2 in set C1 (R3 if transpose)
      queue[tail++] = j2; // add j2 to queue
    }
  }
  if (mark != 1) {
    C = null; // free A' if it was created
  }
  return true;
}

/// Collect matched rows and columns into p and q.
void _matched(int n, Int32List wj, Int32List imatch, int imatch_offset,
                 Int32List p, Int32List q, Int32List cc, Int32List rr, int set, int mark) {
  int kc = cc[set],
      j;
  int kr = rr[set - 1];
  for (j = 0; j < n; j++) {
    if (wj[j] != mark) {
      continue; // skip if j is not in C set
    }
    p[kr++] = imatch[imatch_offset + j];
    q[kc++] = j;
  }
  cc[set + 1] = kc;
  rr[set] = kr;
}

/// Collect unmatched rows into the permutation vector p.
void _unmatched(int m, Int32List wi, Int32List p, Int32List rr, int set) {
  int i,
      kr = rr[set];
  for (i = 0; i < m; i++) if (wi[i] == 0) p[kr++] = i;
  rr[set + 1] = kr;
}

/// Returns true if row i is in R2.

bool _rprune(int i, int j, Float64List aij, Object other) {
  final rr = other as List<int>;
  return i >= rr[1] && i < rr[2];
}

/// Dulmage-Mendelsohn decomposition.
///
/// Compute coarse and then fine Dulmage-Mendelsohn decomposition. seed
/// optionally selects a randomized algorithm.
///
/// [seed] 0: natural, -1: reverse, random order otherwise
/// Returns Dulmage-Mendelsohn analysis, null on error.
Decomposition dmperm(Matrix A, int seed) {
  int m, n, i, j, k, cnz, nc, nb1, nb2;
  Int32List jmatch, imatch, wi, wj, _pinv, Cp, Ci, ps, rs, p, q, cc, rr, r, s;
  bool ok;
  Matrix C;
  Decomposition D, _scc;
  /* Maximum matching */
  if (!csc(A)) {
    return null;
  }
  m = A.m;
  n = A.n;
  D = dalloc(m, n); // allocate result
  if (D == null) {
    return null;
  }
  p = D.p;
  q = D.q;
  r = D.r;
  s = D.s;
  cc = D.cc;
  rr = D.rr;
  jmatch = maxtrans(A, seed); // max transversal
  imatch = jmatch; // imatch = inverse of jmatch
  int imatch_offset = m;
  if (jmatch == null) {
    return _ddone(D, null, jmatch, false);
  }
  /* Coarse decomposition */
  wi = r;
  wj = s; // use r and s as workspace
  for (j = 0; j < n; j++) wj[j] = -1; // unmark all cols for bfs
  for (i = 0; i < m; i++) wi[i] = -1; // unmark all rows for bfs
  _bfs(A, n, wi, wj, q, imatch, imatch_offset, jmatch, 0, 1); // find C1, R1 from C0
  ok = _bfs(A, m, wj, wi, p, jmatch, 0, imatch, imatch_offset, 3); // find R3, C3 from R0
  if (!ok) {
    return _ddone(D, null, jmatch, false);
  }
  _unmatched(n, wj, q, cc, 0); // unmatched set C0
  _matched(n, wj, imatch, imatch_offset, p, q, cc, rr, 1, 1); // set R1 and C1
  _matched(n, wj, imatch, imatch_offset, p, q, cc, rr, 2, -1); // set R2 and C2
  _matched(n, wj, imatch, imatch_offset, p, q, cc, rr, 3, 3); // set R3 and C3
  _unmatched(m, wi, p, rr, 3); // unmatched set R0
  jmatch = null;
  /* Fine decomposition */
  _pinv = pinv(p, m); // pinv=p'
  if (pinv == null) {
    return _ddone(D, null, null, false);
  }
  C = permute(A, _pinv, q, false); // C=A(p,q) (it will hold A(R2,C2))
  _pinv = null;
  if (C == null) {
    return _ddone(D, null, null, false);
  }
  Cp = C.p;
  nc = cc[3] - cc[2]; // delete cols C0, C1, and C3 from C
  if (cc[2] > 0) {
    for (j = cc[2]; j <= cc[3]; j++) {
      Cp[j - cc[2]] = Cp[j];
    }
  }
  C.n = nc;
  if (rr[2] - rr[1] < m) // delete rows R0, R1, and R3 from C
  {
    fkeep(C, _rprune, rr);
    cnz = Cp[nc];
    Ci = C.i;
    if (rr[1] > 0) {
      for (k = 0; k < cnz; k++) {
        Ci[k] -= rr[1];
      }
    }
  }
  C.m = nc;
  _scc = scc(C); // find strongly connected components of C
  if (_scc == null) {
    return _ddone(D, C, null, false);
  }
  /* Combine coarse and fine decompositions */
  ps = _scc.p; // C(ps,ps) is the permuted matrix
  rs = _scc.r; // kth block is rs[k]..rs[k+1]-1
  nb1 = _scc.nb; // # of blocks of A(R2,C2)
  for (k = 0; k < nc; k++) wj[k] = q[ps[k] + cc[2]];
  for (k = 0; k < nc; k++) q[k + cc[2]] = wj[k];
  for (k = 0; k < nc; k++) wi[k] = p[ps[k] + rr[1]];
  for (k = 0; k < nc; k++) p[k + rr[1]] = wi[k];
  nb2 = 0; // create the fine block partitions
  r[0] = s[0] = 0;
  if (cc[2] > 0) {
    nb2++; // leading coarse block A (R1, [C0 C1])
  }
  for (k = 0; k < nb1; k++) // coarse block A (R2,C2)
  {
    r[nb2] = rs[k] + rr[1]; // A (R2,C2) splits into nb1 fine blocks
    s[nb2] = rs[k] + cc[2];
    nb2++;
  }
  if (rr[2] < m) {
    r[nb2] = rr[2]; // trailing coarse block A ([R3 R0], C3)
    s[nb2] = cc[3];
    nb2++;
  }
  r[nb2] = m;
  s[nb2] = n;
  D.nb = nb2;
  _scc = null;
  return _ddone(D, C, null, true);
}
