import "astral-calendrica";

var J2000 = days_from_hours(12) + gregorian_new_year(2000);
aa.J2000 = J2000;

// Return Julian centuries since 2000 at moment 'tee' (in days and
// fractions since EPOCH).

function julian_centuries(tee) {
  return (dynamical_from_universal(tee) - J2000) / 36525;
}
aa.julian_centuries = julian_centuries;

// Return Dynamical time at Universal moment, tee.

function dynamical_from_universal(tee) {
  return tee + ephemeris_correction(tee);
}
aa.dynamical_from_universal = dynamical_from_universal;

// Return Universal moment from Dynamical time, tee.

function universal_from_dynamical(tee) {
  return tee - ephemeris_correction(tee);
}
aa.universal_from_dynamical = universal_from_dynamical;

// Return Dynamical Time minus Universal Time (in days) for
// moment, 'tee'.  Adapted from "Astronomical Algorithms"
// by Jean Meeus, Willmann_Bell, Inc., 1991.

function ephemeris_correction(tee) {
  var year = gregorian_year_from_fixed(ifloor(tee)),
    c = gregorian_date_difference(gregorian_date(1900, JANUARY, 1),
    gregorian_date(year, JULY, 1)) / 36525;

  if ((1988 <= year) && (year <= 2019)) {
    return (year - 1933) / 86400;
  } else if ((1900 <= year) && (year <= 1987)) {
    return poly(c, [-0.00002, 0.000297, 0.025184, -0.181133, 0.553040, -0.861938,
    0.677066, -0.212591]);
  } else if ((1800 <= year) && (year <= 1899)) {
    return poly(c, [-0.000009, 0.003844, 0.083563,
    0.865736, 4.867575, 15.845535,
    31.332267, 38.291999, 28.316289,
    11.636204, 2.043794]);
  } else if ((1700 <= year) && (year <= 1799)) {
    return (1 / 86400 * poly(year - 1700, [8.118780842, -0.005092142,
    0.003336121, -0.0000266484]));
  } else if ((1620 <= year) && (year <= 1699)) {
    return (1 / 86400 * poly(year - 1600, [196.58333, -4.0675, 0.0219167]));
  } else {
    var x = (12 / 24 + gregorian_date_difference(gregorian_date(1810, JANUARY, 1),
    gregorian_date(year, JANUARY, 1)));
    return (1 / 86400 * (((x * x) / 41048480) - 15));
  }
}
aa.ephemeris_correction = ephemeris_correction;


// Return mean longitude of moon (in degrees) at moment
// given in Julian centuries c (including the constant term of the
// effect of the light-time (-0".70).
// Adapted from eq. 47.1 in "Astronomical Algorithms" by Jean Meeus,
// Willmann_Bell, Inc., 2nd ed. with corrections, 2005.
function mean_lunar_longitude(c) {
  return normalized_degrees(poly(c, [218.3164477,
                                     481267.88123421,
                                     -0.0015786,
                                     1 / 538841,
                                     -1 / 65194000]));
}
aa.mean_lunar_longitude = mean_lunar_longitude;

// Return elongation of moon (in degrees) at moment
// given in Julian centuries c.
// Adapted from eq. 47.2 in "Astronomical Algorithms" by Jean Meeus,
// Willmann_Bell, Inc., 2nd ed. with corrections, 2005.

function lunar_elongation(c) {
  return normalized_degrees(poly(c, [297.8501921,
  445267.1114034, -0.0018819,
  1 / 545868, -1 / 113065000]));
}
aa.lunar_elongation = lunar_elongation;


// Return mean anomaly of sun (in degrees) at moment
// given in Julian centuries c.
// Adapted from eq. 47.3 in "Astronomical Algorithms" by Jean Meeus,
// Willmann_Bell, Inc., 2nd ed. with corrections, 2005.

function solar_anomaly(c) {
  return normalized_degrees(poly(c, [357.5291092,
  35999.0502909, -0.0001536,
  1 / 24490000]));
}
aa.solar_anomaly = solar_anomaly;


// Return mean anomaly of moon (in degrees) at moment
// given in Julian centuries c.
// Adapted from eq. 47.4 in "Astronomical Algorithms" by Jean Meeus,
// Willmann_Bell, Inc., 2nd ed. with corrections, 2005."""

function lunar_anomaly(c) {
  return normalized_degrees(poly(c, [134.9633964,
  477198.8675055,
  0.0087414,
  1 / 69699, -1 / 14712000]));
}
aa.lunar_anomaly = lunar_anomaly;


// Return Moon's argument of latitude (in degrees) at moment
// given in Julian centuries 'c'.
// Adapted from eq. 47.5 in "Astronomical Algorithms" by Jean Meeus,
// Willmann_Bell, Inc., 2nd ed. with corrections, 2005."""

function moon_node(c) {
  return normalized_degrees(poly(c, [93.2720950,
  483202.0175233, -0.0036539, -1 / 3526000,
  1 / 863310000]));
}
aa.moon_node = moon_node;


// Return the longitudinal nutation at moment, tee."""

function nutation(tee) {
  var c = julian_centuries(tee),
    cap_A = poly(c, [124.90, -1934.134, 0.002063]),
    cap_B = poly(c, [201.11, 72001.5377, 0.00057]);
  return (-0.004778 * sin_degrees(cap_A) + -0.0003667 * sin_degrees(cap_B));
}
aa.nutation = nutation;


// Return longitude of moon (in degrees) at moment tee.
// Adapted from "Astronomical Algorithms" by Jean Meeus,
// Willmann_Bell, Inc., 2nd ed., 1998.

function lunar_longitude(tee) {
  var c = julian_centuries(tee),
    cap_L_prime = mean_lunar_longitude(c),
    cap_D = lunar_elongation(c),
    cap_M = solar_anomaly(c),
    cap_M_prime = lunar_anomaly(c),
    cap_F = moon_node(c),
    cap_E = poly(c, [1, -0.002516, -0.0000074]),
    args_lunar_elongation = [0, 2, 2, 0, 0, 0, 2, 2, 2, 2, 0, 1, 0, 2, 0, 0, 4, 0, 4, 2, 2, 1,
    1, 2, 2, 4, 2, 0, 2, 2, 1, 2, 0, 0, 2, 2, 2, 4, 0, 3, 2, 4, 0, 2,
    2, 2, 4, 0, 4, 1, 2, 0, 1, 3, 4, 2, 0, 1, 2],
    args_solar_anomaly = [0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1,
    0, 1, -1, 0, 0, 0, 1, 0, -1, 0, -2, 1, 2, -2, 0, 0, -1, 0, 0, 1, -1, 2, 2, 1, -1, 0, 0,
    -1, 0, 1, 0, 1, 0, 0, -1, 2, 1, 0],
    args_lunar_anomaly = [1, -1, 0, 2, 0, 0, -2, -1, 1, 0, -1, 0, 1, 0, 1, 1, -1, 3, -2,
    -1, 0, -1, 0, 1, 2, 0, -3, -2, -1, -2, 1, 0, 2, 0, -1, 1, 0, -1, 2, -1, 1, -2, -1, -1,
    -2, 0, 1, 4, 0, -2, 0, 2, 1, -2, -3, 2, 1, -1, 3],
    args_moon_node = [0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, -2, 2, -2, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, -2, 2, 0, 2, 0, 0, 0, 0,
    0, 0, -2, 0, 0, 0, 0, -2, -2, 0, 0, 0, 0, 0, 0, 0],
    sine_coefficients = [6288774, 1274027, 658314, 213618, -185116, -114332,
    58793, 57066, 53322, 45758, -40923, -34720, -30383,
    15327, -12528, 10980, 10675, 10034, 8548, -7888, -6766, -5163, 4987, 4036, 3994, 3861,
    3665, -2689, -2602, 2390, -2348, 2236, -2120, -2069, 2048, -1773, -1595, 1215, -1110,
    -892, -810, 759, -713, -700, 691, 596, 549, 537, 520, -487, -399, -381, 351, -340, 330,
    327, -323, 299, 294],
    correction = (1 / 1000000 * sigma([sine_coefficients, args_lunar_elongation,
    args_solar_anomaly, args_lunar_anomaly,
    args_moon_node],
    function(v, w, x, y, z) {
      return (v * Math.pow(cap_E, Math.abs(x)) * sin_degrees((w * cap_D) + (x * cap_M) + (y * cap_M_prime) + (z * cap_F)));
    })),
    A1 = 119.75 + (c * 131.849),
    venus = 3958 / 1000000 * sin_degrees(A1),
    A2 = 53.09 + (c * 479264.29),
    jupiter = 318 / 1000000 * sin_degrees(A2),
    flat_earth = 1962 / 1000000 * sin_degrees(cap_L_prime - cap_F);

  return mod(cap_L_prime + correction + venus + jupiter + flat_earth + nutation(tee), 360);
}
aa.lunar_longitude = lunar_longitude;


// Return the latitude of moon (in degrees) at moment, tee.
// Adapted from "Astronomical Algorithms" by Jean Meeus,
// Willmann_Bell, Inc., 1998.
function lunar_latitude(tee) {
  var c = julian_centuries(tee),
      cap_L_prime = mean_lunar_longitude(c),
      cap_D = lunar_elongation(c),
      cap_M = solar_anomaly(c),
      cap_M_prime = lunar_anomaly(c),
      cap_F = moon_node(c),
      cap_E = poly(c, [1, -0.002516, -0.0000074]),
      args_lunar_elongation = [
        0, 0, 0, 2, 2, 2, 2, 0, 2, 0, 2, 2, 2, 2, 2, 2, 2, 0, 4, 0, 0, 0,
        1, 0, 0, 0, 1, 0, 4, 4, 0, 4, 2, 2, 2, 2, 0, 2, 2, 2, 2, 4, 2, 2,
        0, 2, 1, 1, 0, 2, 1, 2, 0, 4, 4, 1, 4, 1, 4, 2],
      args_solar_anomaly = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, -1, -1, -1, 1, 0, 1,
        0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 1,
        0, -1, -2, 0, 1, 1, 1, 1, 1, 0, -1, 1, 0, -1, 0, 0, 0, -1, -2],
      args_lunar_anomaly = [
        0, 1, 1, 0, -1, -1, 0, 2, 1, 2, 0, -2, 1, 0, -1, 0, -1, -1, -1,
        0, 0, -1, 0, 1, 1, 0, 0, 3, 0, -1, 1, -2, 0, 2, 1, -2, 3, 2, -3,
        -1, 0, 0, 1, 0, 1, 1, 0, 0, -2, -1, 1, -2, 2, -2, -1, 1, 1, -2,
        0, 0],
      args_moon_node = [
        1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1, -1, 1, -1, 1, 1, -1, -1,
        -1, 1, 3, 1, 1, 1, -1, -1, -1, 1, -1, 1, -3, 1, -3, -1, -1, 1,
        -1, 1, -1, 1, 1, 1, 1, -1, 3, -1, -1, 1, -1, -1, 1, -1, 1, -1,
        -1, -1, -1, -1, -1, 1],
      sine_coefficients = [
        5128122, 280602, 277693, 173237, 55413, 46271, 32573,
        17198, 9266, 8822, 8216, 4324, 4200, -3359, 2463, 2211,
        2065, -1870, 1828, -1794, -1749, -1565, -1491, -1475,
        -1410, -1344, -1335, 1107, 1021, 833, 777, 671, 607,
        596, 491, -451, 439, 422, 421, -366, -351, 331, 315,
        302, -283, -229, 223, 223, -220, -220, -185, 181,
        -177, 176, 166, -164, 132, -119, 115, 107],
      beta = ((1/1000000) *
              sigma([sine_coefficients, 
                     args_lunar_elongation,
                     args_solar_anomaly,
                     args_lunar_anomaly,
                     args_moon_node],
                    function(v, w, x, y, z) {
                      return (v *
                              Math.pow(cap_E, Math.abs(x)) *
                              sin_degrees((w * cap_D) +
                                          (x * cap_M) +
                                          (y * cap_M_prime) +
                                          (z * cap_F)));
                    })),
      venus = ((175/1000000) *
                  (sin_degrees(119.75 + c * 131.849 + cap_F) +
                   sin_degrees(119.75 + c * 131.849 - cap_F))),
      flat_earth = ((-2235/1000000) * sin_degrees(cap_L_prime) +
                      (127/1000000) * sin_degrees(cap_L_prime - cap_M_prime) +
                     (-115/1000000) * sin_degrees(cap_L_prime + cap_M_prime)),
      extra = ((382/1000000) * sin_degrees(313.45 + c * 481266.484));

  return beta + venus + flat_earth + extra;
}
aa.lunar_latitude = lunar_latitude;



// Return the distance to moon (in meters) at moment, tee.
// Adapted from "Astronomical Algorithms" by Jean Meeus,
// Willmann_Bell, Inc., 2nd ed.
function lunar_distance(tee) {
  var c = julian_centuries(tee),
      cap_D = lunar_elongation(c),
      cap_M = solar_anomaly(c),
      cap_M_prime = lunar_anomaly(c),
      cap_F = moon_node(c),
      cap_E = poly(c, [1, -0.002516, -0.0000074])
      args_lunar_elongation = [
        0, 2, 2, 0, 0, 0, 2, 2, 2, 2, 0, 1, 0, 2, 0, 0, 4, 0, 4, 2, 2, 1,
        1, 2, 2, 4, 2, 0, 2, 2, 1, 2, 0, 0, 2, 2, 2, 4, 0, 3, 2, 4, 0, 2,
        2, 2, 4, 0, 4, 1, 2, 0, 1, 3, 4, 2, 0, 1, 2, 2],
      args_solar_anomaly = [
        0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 1, -1, 0, 0, 0, 1, 0, -1, 0, -2, 1, 2, -2, 0, 0, -1, 0, 0, 1,
        -1, 2, 2, 1, -1, 0, 0, -1, 0, 1, 0, 1, 0, 0, -1, 2, 1, 0, 0],
      args_lunar_anomaly = [
        1, -1, 0, 2, 0, 0, -2, -1, 1, 0, -1, 0, 1, 0, 1, 1, -1, 3, -2,
        -1, 0, -1, 0, 1, 2, 0, -3, -2, -1, -2, 1, 0, 2, 0, -1, 1, 0,
        -1, 2, -1, 1, -2, -1, -1, -2, 0, 1, 4, 0, -2, 0, 2, 1, -2, -3,
        2, 1, -1, 3, -1],
      args_moon_node = [
        0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, -2, 2, -2, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, -2, 2, 0, 2, 0, 0, 0, 0,
        0, 0, -2, 0, 0, 0, 0, -2, -2, 0, 0, 0, 0, 0, 0, 0, -2],
      cosine_coefficients = [
        -20905355, -3699111, -2955968, -569925, 48888, -3149,
        246158, -152138, -170733, -204586, -129620, 108743,
        104755, 10321, 0, 79661, -34782, -23210, -21636, 24208,
        30824, -8379, -16675, -12831, -10445, -11650, 14403,
        -7003, 0, 10056, 6322, -9884, 5751, 0, -4950, 4130, 0,
        -3958, 0, 3258, 2616, -1897, -2117, 2354, 0, 0, -1423,
        -1117, -1571, -1739, 0, -4421, 0, 0, 0, 0, 1165, 0, 0,
        8752],
      correction = sigma ([cosine_coefficients,
                           args_lunar_elongation,
                           args_solar_anomaly,
                           args_lunar_anomaly,
                           args_moon_node],
                          function(v, w, x, y, z) {
                            return (v *
                                    Math.pow(cap_E, Math.abs(x)) * 
                                    cos_degrees((w * cap_D) +
                                                   (x * cap_M) +
                                                   (y * cap_M_prime) +
                                                   (z * cap_F)));});
    return 385000560 + correction;
}
aa.lunar_distance = lunar_distance;


// Return the geocentric apparent lunar diameter of the moon (in
// degrees) at moment, tee.  Adapted from 'Astronomical
// Algorithms' by Jean Meeus, Willmann_Bell, Inc., 2nd ed.
function lunar_diameter(tee) {
  return (1792367000/9) / lunar_distance(tee);
}
aa.lunar_diameter = lunar_diameter;


// Return the moon position (geocentric latitude and longitude [in degrees]
// and distance [in meters]) at moment, tee.
// Adapted from "Astronomical Algorithms" by Jean Meeus,
// Willmann_Bell, Inc., 2nd ed.
function lunar_position(tee) {
  return [lunar_latitude(tee), lunar_longitude(tee), lunar_distance(tee)];
}
aa.lunar_position = lunar_position;


var MEAN_SYNODIC_MONTH = 29.530588861;
aa.MEAN_SYNODIC_MONTH = MEAN_SYNODIC_MONTH;


// Return the moment of n-th new moon after (or before) the new moon
// of January 11, 1.  Adapted from "Astronomical Algorithms"
// by Jean Meeus, Willmann_Bell, Inc., 2nd ed., 1998.
function nth_new_moon(n) {
  var n0 = 24724,
      k = n - n0,
      c = k / 1236.85,
      approx = (J2000 +
                poly(c, [5.09766,
                         MEAN_SYNODIC_MONTH * 1236.85,
                         0.0001437,
                         -0.000000150,
                         0.00000000073])),
      cap_E = poly(c, [1, -0.002516, -0.0000074]),
      solar_anomaly = poly(c, [2.5534,
                               1236.85 * 29.10535669,
                               -0.0000014, -0.00000011]),
      lunar_anomaly = poly(c, [201.5643,
                               (385.81693528 * 1236.85),
                               0.0107582, 0.00001238,
                              -0.000000058]),
      moon_argument = poly(c, [160.7108,
                               (390.67050284 * 1236.85),
                               -0.0016118, -0.00000227,
                               0.000000011]),
      cap_omega = poly(c, [124.7746,
                           (-1.56375588 * 1236.85),
                           0.0020672, 0.00000215]),
      E_factor = [0, 1, 0, 0, 1, 1, 2, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0],
      solar_coeff = [0, 1, 0, 0, -1, 1, 2, 0, 0, 1, 0, 1, 1, -1, 2,
                   0, 3, 1, 0, 1, -1, -1, 1, 0],
      lunar_coeff = [1, 0, 2, 0, 1, 1, 0, 1, 1, 2, 3, 0, 0, 2, 1, 2,
                   0, 1, 2, 1, 1, 1, 3, 4],
      moon_coeff = [0, 0, 0, 2, 0, 0, 0, -2, 2, 0, 0, 2, -2, 0, 0,
                  -2, 0, -2, 2, 2, 2, -2, 0, 0],
      sine_coeff = [-0.40720, 0.17241, 0.01608,
                    0.01039, 0.00739, -0.00514,
                    0.00208, -0.00111, -0.00057,
                    0.00056, -0.00042, 0.00042,
                    0.00038, -0.00024, -0.00007,
                    0.00004, 0.00004, 0.00003,
                    0.00003, -0.00003, 0.00003,
                    -0.00002, -0.00002, 0.00002],
    correction = (-0.00017 * sin_degrees(cap_omega) +
                  sigma([sine_coeff, E_factor, solar_coeff,
                         lunar_coeff, moon_coeff],
                        function(v, w, x, y, z) {
                          return (v *
                                  Math.pow(cap_E, w) *
                                  sin_degrees((x * solar_anomaly) + 
                                              (y * lunar_anomaly) +
                                              (z * moon_argument)));})),
    add_const = [251.88, 251.83, 349.42, 84.66,
                 141.74, 207.14, 154.84, 34.52,
                 207.19, 291.34, 161.72, 239.56,
                 331.55],
    add_coeff = [0.016321, 26.651886, 36.412478,
                 18.206239, 53.303771, 2.453732,
                 7.306860, 27.261239, 0.121824,
                 1.844379, 24.198154, 25.513099,
                 3.592518],
    add_factor = [0.000165, 0.000164, 0.000126,
                  0.000110, 0.000062, 0.000060,
                  0.000056, 0.000047, 0.000042,
                  0.000040, 0.000037, 0.000035,
                  0.000023],
    extra = (0.000325 *
             sin_degrees(poly(c, [299.77, 132.8475848,
                                      -0.009173]))),
    additional = sigma([add_const, add_coeff, add_factor],
                       function(i, j, l) {return l * sin_degrees(i + j * k)});

    return universal_from_dynamical(approx + correction + extra + additional);
}
aa.nth_new_moon = nth_new_moon;

// compute altitude and azimuth
function compute()
{
    var dec_degrees = parseInt(calc.dec_degrees.value, 10);
    var dec_minutes = parseInt(calc.dec_minutes.value, 10);
    var ra_hours    = parseInt(calc.ra_hours.value, 10);
    var ra_minutes  = parseInt(calc.ra_minutes.value, 10);
    var lat_degrees = parseInt(calc.lat_degrees.value, 10);
    var lat_minutes = parseInt(calc.lat_minutes.value, 10);
    var lon_degrees = parseInt(calc.lon_degrees.value, 10);
    var lon_minutes = parseInt(calc.lon_minutes.value, 10);
    var dtg_year    = parseInt(calc.dtg_year.value, 10);
    var dtg_month   = parseInt(calc.dtg_month.value, 10);
    var dtg_day     = parseInt(calc.dtg_day.value, 10);
    var dtg_hour    = parseInt(calc.dtg_hour.value, 10);
    var dtg_minute  = parseInt(calc.dtg_minute.value, 10);
    var dtg_second  = parseInt(calc.dtg_second.value, 10); 

    if (isNaN(dec_degrees)||(Math.abs(dec_degrees) >=  90)||
        isNaN(dec_minutes)||(dec_minutes < 0)||(dec_minutes >=  60)||
        isNaN(ra_hours)   ||(ra_hours    < 0)||(ra_hours    >=  24)||
        isNaN(ra_minutes) ||(ra_minutes  < 0)||(ra_minutes  >=  60))
    {
        window.alert("Invalid Object data!");
        return;
    }

    if (isNaN(lat_degrees)||(Math.abs(lat_degrees) >=  90)||
        isNaN(lat_minutes)||(lat_minutes < 0)||(lat_minutes >=  60)||
        isNaN(lon_degrees)||(Math.abs(lon_degrees) >= 180)||
        isNaN(lon_minutes)||(lon_minutes < 0)||(lon_minutes >=  60))
    {
        window.alert("Invalid Observer data!");
        return;
    }

    var now = new Date(dtg_year, dtg_month - 1, dtg_day, dtg_hour, dtg_minute, dtg_second);

    if (isNaN(now.getTime())||
        (dtg_month  < 1)||(dtg_month  >=  13)||
        (dtg_day    < 1)||(dtg_day    >=  32)||
        (dtg_hour   < 0)||(dtg_hour   >=  24)||
        (dtg_minute < 0)||(dtg_minute >=  60)||
        (dtg_second < 0)||(dtg_second >=  60))
    {
        window.alert("Invalid Date/Time data!");
        return;
    }

    var ra  = ra2real( ra_hours, ra_minutes );
    var dec = dms2real( dec_degrees, dec_minutes );
    var lat = dms2real( lat_degrees, lat_minutes );
    var lon = dms2real( lon_degrees, lon_minutes );

    coord_to_horizon(now, ra, dec, lat, lon);
    calc.alt_degrees.value = hrz_altitude;
    calc.azm_degrees.value = hrz_azimuth;
    calc.date_time.value   = now.toUTCString();
}

// compute horizon coordinates from ra, dec, lat, lon, and utc
// ra, dec, lat, lon in  degrees
// utc is a Date object
// results returned in hrz_altitude, hrz_azimuth
function coord_to_horizon( utc, ra, dec, lat, lon )
{
    // compute hour angle in degrees
    var ha = mean_sidereal_time( utc, lon ) - ra;
    if (ha < 0) ha = ha + 360;

    // convert degrees to radians
    ha  = ha*Math.PI/180
    dec = dec*Math.PI/180
    lat = lat*Math.PI/180

    // compute altitude in radians
    var sin_alt = Math.sin(dec)*Math.sin(lat) + Math.cos(dec)*Math.cos(lat)*Math.cos(ha);
    var alt = Math.asin(sin_alt);
    
    // compute azimuth in radians
    // divide by zero error at poles or if alt = 90 deg
    var cos_az = (Math.sin(dec) - Math.sin(alt)*Math.sin(lat))/(Math.cos(alt)*Math.cos(lat));
    var az  = Math.acos(cos_az);

    // convert radians to degrees
    hrz_altitude = alt*180/Math.PI;
    hrz_azimuth  = az*180/Math.PI;

    // choose hemisphere
    if (Math.sin(ha) > 0) hrz_azimuth = 360 - hrz_azimuth;
}

// Compute the Mean Sidereal Time in units of degrees. 
// Use lon := 0 to get the Greenwich MST. 
// East longitudes are positive; West longitudes are negative
// returns: time in degrees
function mean_sidereal_time(now, lon)
{
    var year   = now.getUTCFullYear();
    var month  = now.getUTCMonth() + 1;
    var day    = now.getUTCDate();
    var hour   = now.getUTCHours();
    var minute = now.getUTCMinutes();
    var second = now.getUTCSeconds();

    if ((month == 1)||(month == 2))
    {
        year  = year - 1;
        month = month + 12;
    }

    var a = Math.floor(year/100);
    var b = 2 - a + Math.floor(a/4);
    var c = Math.floor(365.25*year);
    var d = Math.floor(30.6001*(month + 1));

    // days since J2000.0
    var jd = b + c + d - 730550.5 + day + (hour + minute/60.0 + second/3600.0)/24.0;
    
    // julian centuries since J2000.0
    var jt = jd/36525.0;

    // the mean sidereal time in degrees
    var mst = 280.46061837 + 360.98564736629*jd + 0.000387933*jt*jt - jt*jt*jt/38710000 + lon;

    // in degrees modulo 360.0
    if (mst > 0.0) 
        while (mst > 360.0) mst = mst - 360.0;
    else
        while (mst < 0.0)   mst = mst + 360.0;
        
    return mst;
}

function getDateTime() 
{
    // date and time at Greenwich (UTC)
    var now = new Date();
    calc.date_time.value = now.toUTCString();
}

// format two digits with leading zero if needed
function d2( n )
{
    var s = "";

    if (n < 10) s = "0" + n;
    else        s = n;

    return s;
}

// convert right ascension (hours, minutes) to degrees as real
function ra2real( hr, min )
{
    return 15*(hr + min/60);
}

// convert angle (deg, min) to degrees as real
function dms2real( deg, min )
{
    var rv;

    if (deg < 0)
        rv = deg - min/60;
    else
        rv = deg + min/60;

    return rv;
}




// Return the lunar phase, as an angle in degrees, at moment tee.
// An angle of 0 means a new moon, 90 degrees means the
// first quarter, 180 means a full moon, and 270 degrees
// means the last quarter.
function lunar_phase(tee) {
  var phi = mod(lunar_longitude(tee) - solar_longitude(tee), 360),
      t0 = nth_new_moon(0),
      n = iround((tee - t0) / MEAN_SYNODIC_MONTH),
      phi_prime = (deg(360) * mod((tee - nth_new_moon(n)) / MEAN_SYNODIC_MONTH, 1));

  if (Math.abs(phi - phi_prime) > 180) return phi_prime;

  return phi;
}
aa.lunar_phase = lunar_phase;
