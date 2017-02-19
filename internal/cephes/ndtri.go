// Copyright ©2016 The gonum Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*
 * Cephes Math Library Release 2.1:  January, 1989
 * Copyright 1984, 1987, 1989 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

package cephes

import "math"

// TODO(btracey): There is currently an implementation of this functionality
// in gonum/stat/distuv. Find out which implementation is better, and rectify
// by having distuv call this, or moving this implementation into
// gonum/mathext/internal/gonum.

// math.Sqrt(2*pi)
const s2pi = 2.50662827463100050242E0

// approximation for 0 <= |y - 0.5| <= 3/8
var P0 = [5]float64{
	-5.99633501014107895267E1,
	9.80010754185999661536E1,
	-5.66762857469070293439E1,
	1.39312609387279679503E1,
	-1.23916583867381258016E0,
}

var Q0 = [8]float64{
	/* 1.00000000000000000000E0, */
	1.95448858338141759834E0,
	4.67627912898881538453E0,
	8.63602421390890590575E1,
	-2.25462687854119370527E2,
	2.00260212380060660359E2,
	-8.20372256168333339912E1,
	1.59056225126211695515E1,
	-1.18331621121330003142E0,
}

// Approximation for interval z = math.Sqrt(-2 log y ) between 2 and 8
// i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
var P1 = [9]float64{
	4.05544892305962419923E0,
	3.15251094599893866154E1,
	5.71628192246421288162E1,
	4.40805073893200834700E1,
	1.46849561928858024014E1,
	2.18663306850790267539E0,
	-1.40256079171354495875E-1,
	-3.50424626827848203418E-2,
	-8.57456785154685413611E-4,
}

var Q1 = [8]float64{
	/*  1.00000000000000000000E0, */
	1.57799883256466749731E1,
	4.53907635128879210584E1,
	4.13172038254672030440E1,
	1.50425385692907503408E1,
	2.50464946208309415979E0,
	-1.42182922854787788574E-1,
	-3.80806407691578277194E-2,
	-9.33259480895457427372E-4,
}

// Approximation for interval z = math.Sqrt(-2 log y ) between 8 and 64
// i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
var P2 = [9]float64{
	3.23774891776946035970E0,
	6.91522889068984211695E0,
	3.93881025292474443415E0,
	1.33303460815807542389E0,
	2.01485389549179081538E-1,
	1.23716634817820021358E-2,
	3.01581553508235416007E-4,
	2.65806974686737550832E-6,
	6.23974539184983293730E-9,
}

var Q2 = [8]float64{
	/*  1.00000000000000000000E0, */
	6.02427039364742014255E0,
	3.67983563856160859403E0,
	1.37702099489081330271E0,
	2.16236993594496635890E-1,
	1.34204006088543189037E-2,
	3.28014464682127739104E-4,
	2.89247864745380683936E-6,
	6.79019408009981274425E-9,
}

// Ndtri returns the argument, x, for which the area under the
// Gaussian probability density function (integrated from
// minus infinity to x) is equal to y.
func Ndtri(y float64) float64 {
	// For small arguments 0 < y < exp(-2), the program computes
	// z = math.Sqrt( -2.0 * math.Log(y) );  then the approximation is
	// x = z - math.Log(z)/z  - (1/z) P(1/z) / Q(1/z).
	// There are two rational functions P/Q, one for 0 < y < exp(-32)
	// and the other for y up to exp(-2).  For larger arguments,
	// w = y - 0.5, and  x/math.Sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).
	if y <= 0.0 {
		if y < 0 {
			panic(badParamOutOfBounds)
		}
		return math.Inf(-1)
	}
	if y >= 1.0 {
		if y > 1 {
			panic(badParamOutOfBounds)
		}
		return math.Inf(1)
	}
	code := 1
	if y > 1.0-0.13533528323661269189 { /* 0.135... = exp(-2) */
		y = 1.0 - y
		code = 0
	}

	if y > 0.13533528323661269189 {
		y -= 0.5
		y2 := y * y

		ap0 := (((P0[0]*y2+P0[1])*y2+P0[2])*y2+P0[3])*y2 + P0[4]
		aq0 := (((((((y2+Q0[0])*y2+Q0[1])*y2+Q0[2])*y2+Q0[3])*y2+Q0[4])*y2+Q0[5])*y2+Q0[6])*y2 + Q0[7]
		return s2pi * (y + y*(y2*ap0/aq0))
	}

	x := math.Sqrt(-2.0 * math.Log(y))
	x0 := x - math.Log(x)/x

	z := 1.0 / x
	var x1 float64
	if x < 8.0 { /* y > exp(-32) = 1.2664165549e-14 */
		ap1 := (((((((P1[0]*z+P1[1])*z+P1[2])*z+P1[3])*z+P1[4])*z+P1[5])*z+P1[6])*z+P1[7])*z + P1[8]
		aq1 := (((((((z+Q1[0])*z+Q1[1])*z+Q1[2])*z+Q1[3])*z+Q1[4])*z+Q1[5])*z+Q1[6])*z + Q1[7]
		x1 = z * ap1 / aq1
	} else {
		ap2 := (((((((P2[0]*z+P2[1])*z+P2[2])*z+P2[3])*z+P2[4])*z+P2[5])*z+P2[6])*z+P2[7])*z + P2[8]
		aq2 := (((((((z+Q2[0])*z+Q2[1])*z+Q2[2])*z+Q2[3])*z+Q2[4])*z+Q2[5])*z+Q2[6])*z + Q2[7]
		x1 = z * ap2 / aq2
	}
	x = x0 - x1
	if code != 0 {
		x *= -1
	}
	return x
}
