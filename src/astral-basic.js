var τ = 2 * Math.PI,
  degrees = 360 / τ,
  radians = τ / 360;


// Return the angular data structure.
function degrees_minutes_seconds(d, m, s) {
  return [d, m, s];
}

// Return an angle in degrees:minutes:seconds from angle,
// 'alpha' in degrees.

function angle_from_degrees(alpha) {
  var d = ifloor(alpha),
    m = ifloor(60 * mod(alpha, 1)),
    s = mod(alpha * 60 * 60, 60);
  return degrees_minutes_seconds(d, m, s);
}
aa.angle_from_degrees = angle_from_degrees;

// Return decimal degrees

function decimal_degrees(d, m, s) {
  return d + (m + s / 60) / 60;
}
aa.decimal_degrees = decimal_degrees;

// define sign of number
function signum(x) {
  return x ? x < 0 ? -1 : 1 : 0;
};
aa.signum = signum;
Math.signum = signum;

// round a-la Common Lisp
function round(n) {
  return Math.signum(n) * Math.round(Math.abs(n));
}
aa.round = round;

// function integer(n) { return n - (mod(n,1));}
// aa.int = integer;

function mod(m, n) {
  return m - (n * ifloor(m / n));
}
aa.mod = mod;

function quotient(m, n) {
  return ifloor(m / n);
}
aa.quotient = quotient;

function ifloor(n) {
  return Math.floor(n);
}
aa.ifloor = ifloor;

function iround(n) {
  return round(n);
}
aa.iround = iround;

// Return same as mod() with y instead of 0

function amod(x, y) {
  return y + (mod(x, -y));
}
aa.amod = amod;

// Return first integer greater or equal to initial index, i,
// such that condition, p, holds.

function next(i, p) {
  return (p(i) ? i : next(i + 1, p));
}
aa.next = next;

// Return last integer greater or equal to initial index, i,
// such that condition, p, holds.

function final(i, p) {
  return (!p(i)) ? i - 1 : final(i + 1, p);
}
aa.final = final;


// Return the sum of f(i) from i=k, k+1, ... till p(i) holds true or 0.
// This is a tail recursive implementation."""

function summa(f, k, p) {
  return ((!p(k)) ? 0 : f(k) + summa(f, k + 1, p));
}
aa.summa = summa;

// Bisection search for x in [lo, hi] such that condition 'e' holds.
// p determines when to go left.

function binary_search(lo, hi, p, e) {
  var x = (lo + hi) / 2;
  if (p(lo, hi)) return x;
  if (e(x)) return binary_search(lo, x, p, e);
  return binary_search(x, hi, p, e);
}
aa.binary_search = binary_search;

// Find inverse of angular function 'f' at 'y' within interval [a,b].
// Default precision is 0.00001.

function invert_angular(f, y, a, b, prec) {
  prec = prec || 1e-5;
  return binary_search(a, b,

  function(l, h) {
    return ((h - l) <= prec);
  },

  function(x) {
    return mod((f(x) - y), 360) < 180;
  });
}
aa.invert_angular = invert_angular;

// homemade zip (could have used Underscore.js)

function zip(arrays) {
  return arrays.length == 0 ? [] : arrays[0].map(function(_, i) {
    return arrays.map(function(array) {
      return array[i];
    });
  });
}
aa.zip = zip;

function sigma(l, b) {
  return zip(l).map(function(v) {
    return b.apply(null, v);
  })
    .reduce(function(memo, n) {
    return memo + n;
  }, 0);
}
aa.sigma = sigma;

// Calculate polynomial with coefficients 'a' at point 'x'.
// The polynomial is
//  a[0] + a[1] * x + a[2] * x^2 + ... + a[n] * x^n
// The code below implements Horner's Rule

function poly(x, a) {
  var l = a.length,
    p = a[l - 1],
    i = l - 2;
  while (i >= 0) {
    p = p * x + a[i];
    --i;
  }
  return p;
}
aa.poly = poly;

function polyAlt(x, a) {
  var l = a.length,
    p = a[l - 1],
    y = x,
    i = l - 2;
  while (i >= 0) {
    p = p.multiply(y).add(a[i]);
    --i;
  }
  return p;
}
aa.polyAlt = polyAlt;


function sin_degrees(α) {
  return Math.sin(α * radians);
}
aa.sin_degrees = sin_degrees;

function cos_degrees(α) {
  return Math.cos(α * radians);
}
aa.cos_degrees = cos_degrees;

// Return a normalize angle α to range [0,360) degrees.

function normalized_degrees(α) {
  return mod(α, 360);
}
