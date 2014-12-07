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

/// Depth-first-search of the graph of a matrix, starting at node j.
///
/// [j] starting node.
/// [G] graph to search (G.p modified, then restored).
/// [top] stack[top..n-1] is used on input.
/// [xi] size n, stack containing nodes traversed.
/// [xi_offset] the index of the first element in array xi.
/// [pstack] size n, work array.
/// [pstack_offset] the index of the first element in array pstack.
/// [pinv] mapping of rows to columns of G, ignored if null.
/// [pinv_offset] the index of the first element in array pinv.
/// Returns new value of top, -1 on error.
int cs_dfs(int j, DZcs G, int top, Int32List xi, int xi_offset, Int32List pstack, int pstack_offset,
           Int32List pinv, int pinv_offset) {
  int i,
      p,
      p2,
      jnew,
      head = 0;
  Int32List Gp, Gi;
  bool done;
  if (!CS_CSC(G) || xi == null || pstack == null) {
    return -1;
  }
  Gp = G.p;
  Gi = G.i;
  xi[xi_offset + 0] = j; // initialize the recursion stack
  while (head >= 0) {
    j = xi[xi_offset + head]; // get j from the top of the recursion stack
    jnew = pinv != null ? (pinv[pinv_offset + j]) : j;
    if (!_CS_MARKED(Gp, j)) {
      _CS_MARK(Gp, j); // mark node j as visited
      pstack[pstack_offset + head] = (jnew < 0) ? 0 : _CS_UNFLIP(Gp[jnew]);
    }
    done = true; // node j done if no unvisited neighbors
    p2 = (jnew < 0) ? 0 : _CS_UNFLIP(Gp[jnew + 1]);
    for (p = pstack[pstack_offset + head]; p < p2; p++) // examine all neighbors of j
    {
      i = Gi[p]; // consider neighbor node i
      if (_CS_MARKED(Gp, i)) {
        continue; // skip visited node i
      }
      pstack[pstack_offset + head] = p; // pause depth-first search of node j
      xi[xi_offset + ++head] = i; // start dfs at node i
      done = false; // node j is not done
      break; // break, to start dfs (i)
    }
    if (done) // depth-first search at node j is done
    {
      head--; // remove j from the recursion stack
      xi[xi_offset + --top] = j; // and place in the output stack
    }
  }
  return top;
}
