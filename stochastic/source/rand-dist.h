/*
Stochastic simulator for zebrafish segmentation
Copyright (C) 2012 Ahmet Ay, Jack Holland, Adriana Sperlea

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef RAND_DIST_H
#define RAND_DIST_H

#include <math.h>

// uniform distribution (returns a double from 0.0-1.0)
inline double unif_dist () {
	return rand() / (((double)RAND_MAX) + 1);
}

// binomial distribution (returns an integer from 0-n)
inline int bino_dist (int n, double p) {
	int x = 0;
	for (int i = 0; i < n; i++) {
		if (unif_dist() < p) {
			x++;
		}
	}
	return x;
}

// exponential distribution (returns a double based on mean)
inline double expo_dist (double mean) {
	return -log(1.0 - unif_dist()) / mean;
}

// Poisson distribution (returns an integer based on mean)
inline int pois_dist (double mean) {
	double L = pow(M_E, -mean);
	int k = 0;
	double p = 1;
	do {
		k++;
		p *= unif_dist();
	} while (p > L);
	return k - 1;
}

// Pk distribution, which is just a logarithmic uniform distribution
inline double pk_dist () {
	return log(1 / unif_dist());
}

#endif

