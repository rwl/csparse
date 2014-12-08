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

/// Depth-first search and postorder of a tree rooted at node j
///
/// [j] postorder of a tree rooted at node j.
/// [k] number of nodes ordered so far.
/// [head] head[i] is first child of node i; -1 on output.
/// [head_offset] the index of the first element in array head.
/// [next] next[i] is next sibling of i or -1 if none.
/// [next_offset] the index of the first element in array next.
/// [post] postordering.
/// [post_offset] the index of the first element in array post.
/// [stack] size n, work array.
/// [stack_offset] the index of the first element in array stack.
/// Returns new value of k, -1 on error.
int tdfs(int j, int k, Int32List head, int head_offset, Int32List next, int next_offset,
            Int32List post, int post_offset, Int32List stack, int stack_offset) {
  int top = 0;
  if (head == null || next == null || post == null || stack == null) {
    return -1;
  }
  stack[stack_offset + 0] = j; // place j on the stack
  while (top >= 0) // while (stack is not empty)
  {
    final p = stack[stack_offset + top]; // p = top of stack
    final i = head[head_offset + p]; // i = youngest child of p
    if (i == -1) {
      top--; // p has no unordered children left
      post[post_offset + (k++)] = p; // node p is the kth postordered node
    } else {
      head[head_offset + p] = next[next_offset + i]; // remove i from children of p
      stack[stack_offset + (++top)] = i; // start dfs on child node i
    }
  }
  return k;
}
