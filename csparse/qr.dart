/* ***** BEGIN LICENSE BLOCK *****
 *
 * CSparse: a Concise Sparse matrix package.
 * Copyright (c) 2006, Timothy A. Davis.
 * http://www.cise.ufl.edu/research/sparse/CSparse
 *
 * -------------------------------------------------------------------------
 *
 * CSparseJ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * CSparseJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * ***** END LICENSE BLOCK ***** */

part of edu.emory.mathcs.csparse;

/**
 * Sparse QR factorization.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
//public class Dcs_qr {

/**
 * Sparse QR factorization of an m-by-n matrix A, A= Q*R
 *
 * @param A
 *            column-compressed matrix
 * @param S
 *            symbolic QR analysis
 * @return numeric QR factorization, null on error
 */
Dcsn cs_qr(Dcs A, Dcss S) {
    List<double> Rx, Vx, Ax, x, Beta;
    int i, k, p, n, vnz, p1, top, m2, len, col, rnz;
    List<int> s, leftmost, Ap, Ai, parent, Rp, Ri, Vp, Vi, w, pinv, q;
    Dcs R, V;
    Dcsn N;
    if (!CS_CSC(A) || S == null)
        return (null);
    n = A.n;
    Ap = A.p;
    Ai = A.i;
    Ax = A.x;
    q = S.q;
    parent = S.parent;
    pinv = S.pinv;
    m2 = S.m2;
    vnz = S.lnz;
    rnz = S.unz;
    leftmost = S.leftmost;
    w = new List<int>(m2 + n); /* get int workspace */
    x = new List<double>(m2); /* get double workspace */
    N = new Dcsn(); /* allocate result */
    s = w;
    int s_offset = m2; /* s is size n */
    for (k = 0; k < m2; k++)
        x[k] = 0.0; /* clear workspace x */
    N.L = V = cs_spalloc(m2, n, vnz, true, false); /* allocate result V */
    N.U = R = cs_spalloc(m2, n, rnz, true, false); /* allocate result R */
    N.B = Beta = new List<double>(n); /* allocate result Beta */
    Rp = R.p;
    Ri = R.i;
    Rx = R.x;
    Vp = V.p;
    Vi = V.i;
    Vx = V.x;
    for (i = 0; i < m2; i++)
        w[i] = -1; /* clear w, to mark nodes */
    rnz = 0;
    vnz = 0;
    for (k = 0; k < n; k++) /* compute V and R */
    {
        Rp[k] = rnz; /* R(:,k) starts here */
        Vp[k] = p1 = vnz; /* V(:,k) starts here */
        w[k] = k; /* add V(k,k) to pattern of V */
        Vi[vnz++] = k;
        top = n;
        col = q != null ? q[k] : k;
        for (p = Ap[col]; p < Ap[col + 1]; p++) /* find R(:,k) pattern */
        {
            i = leftmost[Ai[p]]; /* i = min(find(A(i,q))) */
            for (len = 0; w[i] != k; i = parent[i]) /* traverse up to k */
            {
                s[s_offset + (len++)] = i;
                w[i] = k;
            }
            while (len > 0)
                s[s_offset + (--top)] = s[s_offset + (--len)]; /* push path on stack */
            i = pinv[Ai[p]]; /* i = permuted row of A(:,col) */
            x[i] = Ax[p]; /* x (i) = A(:,col) */
            if (i > k && w[i] < k) /* pattern of V(:,k) = x (k+1:m) */
            {
                Vi[vnz++] = i; /* add i to pattern of V(:,k) */
                w[i] = k;
            }
        }
        for (p = top; p < n; p++) /* for each i in pattern of R(:,k) */
        {
            i = s[s_offset + p]; /* R(i,k) is nonzero */
            cs_happly(V, i, Beta[i], x); /* apply (V(i),Beta(i)) to x */
            Ri[rnz] = i; /* R(i,k) = x(i) */
            Rx[rnz++] = x[i];
            x[i] = 0.0;
            if (parent[i] == k)
                vnz = cs_scatter(V, i, 0.0, w, null, k, V, vnz);
        }
        for (p = p1; p < vnz; p++) /* gather V(:,k) = x */
        {
            Vx[p] = x[Vi[p]];
            x[Vi[p]] = 0.0;
        }
        Ri[rnz] = k; /* R(k,k) = norm (x) */
        List<double> beta = new List<double>(1);
        beta[0] = Beta[k];
        Rx[rnz++] = cs_house(Vx, p1, beta, vnz - p1); /* [v,beta]=house(x) */
        Beta[k] = beta[0];
    }
    Rp[n] = rnz; /* finalize R */
    Vp[n] = vnz; /* finalize V */
    return N;
}
