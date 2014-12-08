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

/// Determines if j is a leaf of the skeleton matrix and find lowest common
/// ancestor (lca).
int leaf(int i, int j, Int32List first, int first_offset, Int32List maxfirst, int maxfirst_offset,
            Int32List prevleaf, int prevleaf_offset, Int32List ancestor, int ancestor_offset, Int32List jleaf) {
  int q, s, sparent, jprev;
  if (first == null || maxfirst == null || prevleaf == null || ancestor == null || jleaf == null) {
    return -1;
  }
  jleaf[0] = 0;
  if (i <= j || first[first_offset + j] <= maxfirst[maxfirst_offset + i]) {
    return -1; // j not a leaf
  }
  maxfirst[maxfirst_offset + i] = first[first_offset + j]; // update max first[j] seen so far
  jprev = prevleaf[prevleaf_offset + i]; // jprev = previous leaf of ith subtree
  prevleaf[prevleaf_offset + i] = j;
  jleaf[0] = (jprev == -1) ? 1 : 2; // j is first or subsequent leaf
  if (jleaf[0] == 1) {
    return i; // if 1st leaf, q = root of ith subtree
  }
  for (q = jprev; q != ancestor[ancestor_offset + q]; q = ancestor[ancestor_offset + q]);
  for (s = jprev; s != q; s = sparent) {
    sparent = ancestor[ancestor_offset + s]; // path compression
    ancestor[ancestor_offset + s] = q;
  }
  return q; // q = least common ancestor (jprev,j)
}
