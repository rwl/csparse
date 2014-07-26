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
/*
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
*/

/**
 * Load a sparse matrix from a file.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
//public class Dcs_load {

/**
 * Loads a triplet matrix T from a file. Each line of the file contains
 * three values: a row index i, a column index j, and a numerical value aij.
 *
 * @param fileName
 *            file name
 * @param base
 *            index base
 * @return T if successful, null on error
 */
Dcs cs_load(InputStream _in, [int base=0]) {
    int i, j;
    double x;
    Dcs T;
    Reader r = new InputStreamReader(_in);
    BufferedReader br = new BufferedReader(r);

    T = cs_spalloc(0, 0, 1, true, true); /* allocate result */
    String line;
    try {
        while ((line = br.readLine()) != null) {
            List<String> tokens = line.trim().split("\\s+");
            if (tokens.length != 3) {
                return null;
            }
            i = int.parse(tokens[0]) - base;
            j = int.parse(tokens[1]) - base;
            x = double.parse(tokens[2]);
            if (!cs_entry(T, i, j, x))
                return (null);
        }
        r.close();
        br.close();
    } on IOException catch (e) {
        return (null);
    }
    return (T);
}
