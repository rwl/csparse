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

/// Sparse LU factorization of a square matrix, PAQ = LU.
///
/// [tol] partial pivoting threshold (1 for partial pivoting).
/// Returns numeric LU factorization, null on error.
Numeric lu(Matrix A, Symbolic S, double tol) {
  Matrix L, U;
  Numeric N;
  double pivot, a, t;
  Float64List Lx, Ux, x;
  Int32List Lp, Li, Up, Ui, pinv, xi, q;
  int n, ipiv, k, top, p, i, col, lnz, unz;
  if (!csc(A) || S == null) {
    return null; // check inputs
  }
  n = A.n;
  q = S.q;
  lnz = S.lnz;
  unz = S.unz;
  x = new Float64List(n); // get double workspace
  xi = new Int32List(2 * n); // get int workspace
  N = new Numeric(); // allocate result
  N.L = L = spalloc(n, n, lnz, true, false); // allocate result L
  N.U = U = spalloc(n, n, unz, true, false); // allocate result U
  N.pinv = pinv = new Int32List(n); // allocate result pinv
  Lp = L.p;
  Up = U.p;
  for (i = 0; i < n; i++) {
    x[i] = 0.0; // clear workspace
  }
  for (i = 0; i < n; i++) {
    pinv[i] = -1; // no rows pivotal yet
  }
  for (k = 0; k <= n; k++) {
    Lp[k] = 0; // no cols of L yet
  }
  lnz = unz = 0;
  for (k = 0; k < n; k++) // compute L(:,k) and U(:,k)
  {
    /* Triangular solve */
    Lp[k] = lnz; // L(:,k) starts here
    Up[k] = unz; // U(:,k) starts here
    if (lnz + n > L.nzmax) {
      sprealloc(L, 2 * L.nzmax + n);
    }
    if (unz + n > U.nzmax) {
      sprealloc(U, 2 * U.nzmax + n);
    }
    Li = L.i;
    Lx = L.x;
    Ui = U.i;
    Ux = U.x;
    col = q != null ? (q[k]) : k;
    top = spsolve(L, A, col, xi, x, pinv, true); // x = L\A(:,col)
    /* Find pivot */
    ipiv = -1;
    a = -1.0;
    for (p = top; p < n; p++) {
      i = xi[p]; // x(i) is nonzero
      if (pinv[i] < 0) // row i is not yet pivotal
      {
        if ((t = x[i].abs()) > a) {
          a = t; // largest pivot candidate so far
          ipiv = i;
        }
      } else // x(i) is the entry U(pinv[i],k)
      {
        Ui[unz] = pinv[i];
        Ux[unz++] = x[i];
      }
    }
    if (ipiv == -1 || a <= 0) {
      return null;
    }
    if (pinv[col] < 0 && x[col].abs() >= a * tol) {
      ipiv = col;
    }
    /* Divide by pivot */
    pivot = x[ipiv]; // the chosen pivot
    Ui[unz] = k; // last entry in U(:,k) is U(k,k)
    Ux[unz++] = pivot;
    pinv[ipiv] = k; // ipiv is the kth pivot row
    Li[lnz] = ipiv; // first entry in L(:,k) is L(k,k) = 1
    Lx[lnz++] = 1.0;
    for (p = top; p < n; p++) // L(k+1:n,k) = x / pivot
    {
      i = xi[p];
      if (pinv[i] < 0) // x(i) is an entry in L(:,k)
      {
        Li[lnz] = i; // save unpermuted row in L
        Lx[lnz++] = x[i] / pivot; // scale pivot column
      }
      x[i] = 0.0; // x [0..n-1] = 0 for next k
    }
  }
  /* Finalize L and U */
  Lp[n] = lnz;
  Up[n] = unz;
  Li = L.i; // fix row indices of L for final pinv
  for (p = 0; p < lnz; p++) {
    Li[p] = pinv[Li[p]];
  }
  sprealloc(L, 0); // remove extra space from L and U
  sprealloc(U, 0);
  return N;
}
