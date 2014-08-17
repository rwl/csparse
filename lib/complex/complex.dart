/*
 * CXSparse: a Concise Sparse matrix package.
 * Copyright (C) 2006-2011, Timothy A. Davis.
 * Copyright (C) 2011-2012, Richard W. Lincoln.
 * http://www.cise.ufl.edu/research/sparse/CXSparse
 *
 * -------------------------------------------------------------------------
 *
 * CXSparseJ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * CXSparseJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 *
 */

part of edu.emory.mathcs.csparse.complex;

//public class DZcs_complex {

Float64List cs_czero()
{
	return new Float64List.fromList([0.0, 0.0]) ;
}

Float64List cs_cone()
{
	return new Float64List.fromList([1.0, 0.0]) ;
}

double cs_creal(Float64List x)
{
	return x [0] ;
}

double cs_cimag(Float64List x)
{
	return x [1] ;
}

Float64List cs_cget(Float64List x, int idx)
{
	return new Float64List.fromList([x [idx], x [idx + 1]]) ;
}

void cs_cset(Float64List x, int idx, Float64List val)
{
	x [idx] = val [0] ;
	x [idx + 1] = val [1] ;
}

bool cs_cequal(Float64List x, Float64List y)
{
	return cs_cequal (x, y, 1e-14) ;
}

bool cs_cequal(Float64List x, Float64List y, double tol)
{
	if (cs_cabs (x [0] - y [0], x [1] - y [1]) <= Math.abs(tol))
	{
		return true ;
	}
	else
	{
		return false ;
	}
}

double cs_cabs(Float64List x)
{
	double absX = Math.abs(x [0]) ;
	double absY = Math.abs(x [1]) ;

	if (absX == 0.0 && absY == 0.0)
	{
		return 0.0 ;
	}
	else if (absX >= absY)
	{
		double d = x [1] / x [0] ;
		return absX * Math.sqrt(1.0 + d * d) ;
	}
	else
	{
	    double d = x[0] / x[1] ;
	    return absY * Math.sqrt(1.0 + d * d) ;
	}
}

double cs_cabs(double re, double im)
{
	double absX = Math.abs(re) ;
	double absY = Math.abs(im) ;

	if (absX == 0.0 && absY == 0.0)
	{
		return 0.0 ;
	}
	else if (absX >= absY)
	{
		double d = im / re ;
		return absX * Math.sqrt(1.0 + d * d) ;
	}
	else
	{
		double d = re / im ;
		return absY * Math.sqrt(1.0 + d * d) ;
	}
}

Float64List cs_conj(Float64List x)
{
	return new Float64List.fromList([x[0], -x[1]]) ;
}

Float64List cs_cdiv(Float64List x, double re, double im)
{
	Float64List z = new Float64List(2) ;
	double scalar ;

	if (Math.abs(re) >= Math.abs(im))
	{
		scalar = 1.0 / (re + im * (im / re)) ;

		z [0] = scalar * (x [0] + x [1] * (im / re)) ;
		z [1] = scalar * (x [1] - x [0] * (im / re)) ;
	}
	else
	{
		scalar = 1.0 / (re * (re / im) + im) ;

		z [0] = scalar * (x [0] * (re / im) + x [1]) ;
		z [1] = scalar * (x [1] * (re / im) - x [0]) ;
	}

	return z ;
}

Float64List cs_cdiv(Float64List x, Float64List y)
{
	return cs_cdiv (x, y [0], y [1]) ;
}

Float64List cs_cplus(Float64List x, Float64List y)
{
	return new Float64List.fromList([x [0] + y [0], x [1] + y [1]]) ;
}

Float64List cs_cminus(final Float64List x, final Float64List y)
{
	return new Float64List.fromList([x [0] - y [0], x [1] - y [1]]) ;
}

Float64List cs_cmult(Float64List x, double y)
{
	return new Float64List.fromList([x [0] * y, x [1] * y]); ;
}

Float64List cs_cmult(Float64List x, Float64List y)
{
	return new Float64List.fromList([
		x [0] * y [0] - x [1] * y [1],
		x [1] * y [0] + x [0] * y [1]
	]) ;
}

Float64List cs_cneg(Float64List x)
{
	Float64List neg_one = new Float64List.fromList([-1.0, 0.0]) ;
	return cs_cmult(neg_one, x) ;
}

Float64List cs_csqrt(Float64List x)
{
	Float64List z = new Float64List(2) ;
	double absx = cs_cabs (x) ;
	double tmp ;
	if (absx > 0.0)
	{
		if (x [0] > 0.0)
		{
			tmp = Math.sqrt(0.5 * (absx + x [0])) ;
			z [0] = tmp ;
			z [1] = 0.5 * (x [1] / tmp) ;
		}
		else
		{
			tmp = Math.sqrt(0.5 * (absx - x [0])) ;
			if (x [1] < 0.0)
			{
				tmp = -tmp ;
			}
			z [0] = 0.5 * (x [1] / tmp) ;
			z [1] = tmp ;
		}
	}
	else
	{
		z [0] = 0.0 ;
		z [1] = 0.0 ;
	}
	return z ;
}

Float64List cs_csquare(Float64List x)
{
	return cs_cmult (x, x) ;
}

//}
