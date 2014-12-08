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

Float64List czero() => new Float64List.fromList([0.0, 0.0]);

Float64List cone() => new Float64List.fromList([1.0, 0.0]);

double _creal(Float64List x) => x[0];

double _cimag(Float64List x) => x[1];

Float64List _cget(Float64List x, int idx) {
  return new Float64List.fromList([x[idx], x[idx + 1]]);
}

void _cset(Float64List x, int idx, Float64List val) {
  x[idx] = val[0];
  x[idx + 1] = val[1];
}

bool cequal(Float64List x, Float64List y, [double tol = 1e-14]) {
  return cabs(x[0] - y[0], x[1] - y[1]) <= tol.abs();
}

double cabs_list(Float64List x) {
  double absX = x[0].abs();
  double absY = x[1].abs();

  if (absX == 0.0 && absY == 0.0) {
    return 0.0;
  } else if (absX >= absY) {
    double d = x[1] / x[0];
    return absX * math.sqrt(1.0 + d * d);
  } else {
    double d = x[0] / x[1];
    return absY * math.sqrt(1.0 + d * d);
  }
}

double cabs(double re, double im) {
  double absX = re.abs();
  double absY = im.abs();

  if (absX == 0.0 && absY == 0.0) {
    return 0.0;
  } else if (absX >= absY) {
    double d = im / re;
    return absX * math.sqrt(1.0 + d * d);
  } else {
    double d = re / im;
    return absY * math.sqrt(1.0 + d * d);
  }
}

Float64List conj(Float64List x) {
  return new Float64List.fromList([x[0], -x[1]]);
}

Float64List cdiv(Float64List x, double re, double im) {
  Float64List z = new Float64List(2);
  double scalar;

  if (re.abs() >= im.abs()) {
    scalar = 1.0 / (re + im * (im / re));

    z[0] = scalar * (x[0] + x[1] * (im / re));
    z[1] = scalar * (x[1] - x[0] * (im / re));
  } else {
    scalar = 1.0 / (re * (re / im) + im);

    z[0] = scalar * (x[0] * (re / im) + x[1]);
    z[1] = scalar * (x[1] * (re / im) - x[0]);
  }

  return z;
}

Float64List cdiv_list(Float64List x, Float64List y) {
  return cdiv(x, y[0], y[1]);
}

Float64List cplus(Float64List x, Float64List y) {
  return new Float64List.fromList([x[0] + y[0], x[1] + y[1]]);
}

Float64List cminus(final Float64List x, final Float64List y) {
  return new Float64List.fromList([x[0] - y[0], x[1] - y[1]]);
}

Float64List cmult(Float64List x, double y) {
  return new Float64List.fromList([x[0] * y, x[1] * y]);
}

Float64List cmult_list(Float64List x, Float64List y) {
  return new Float64List.fromList([x[0] * y[0] - x[1] * y[1], x[1] * y[0] + x[0] * y[1]]);
}

Float64List cneg(Float64List x) {
  Float64List neg_one = new Float64List.fromList([-1.0, 0.0]);
  return cmult_list(neg_one, x);
}

Float64List csqrt(Float64List x) {
  Float64List z = new Float64List(2);
  double absx = cabs_list(x);
  double tmp;
  if (absx > 0.0) {
    if (x[0] > 0.0) {
      tmp = math.sqrt(0.5 * (absx + x[0]));
      z[0] = tmp;
      z[1] = 0.5 * (x[1] / tmp);
    } else {
      tmp = math.sqrt(0.5 * (absx - x[0]));
      if (x[1] < 0.0) {
        tmp = -tmp;
      }
      z[0] = 0.5 * (x[1] / tmp);
      z[1] = tmp;
    }
  } else {
    z[0] = 0.0;
    z[1] = 0.0;
  }
  return z;
}

Float64List csquare(Float64List x) {
  return cmult_list(x, x);
}
