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

#ifndef UPDATES_H
#define UPDATES_H

#include "macros.h"
#include "parameters.h"

inline void update_a0_27 (double* ca, parameters* p, int* cx, int xcell, double* a0) { // If her1 mRNA is updated
	double old = ca[0] + ca[27];
	ca[0] = p->psh1 * cx[xcell + 0]; 		// Her1 protein synthesis
	ca[27] = p->mdh1 * cx[xcell + 0]; 		// her1 mRNA degradation
	*a0 += ca[0] + ca[27] - old;
}

inline void update_a1_2_4_6 (double* ca, parameters* p, int* cx, int xcell, double* a0) { // If Her1 protein is updated
	double old = ca[1] + ca[2] + ca[4] + ca[6];
	ca[1] = p->pdh1 * cx[xcell + 4]; 								// Her1 protein degradation
	ca[2] = p->dah1h1 * cx[xcell + 4] * (cx[xcell + 4] - 1) / 2; 	// Her1-Her1 dimer association
	ca[4] = p->dah1h7 * cx[xcell + 4] * cx[xcell + 5]; 				// Her1-Her7 dimer association
	ca[6] = p->dah1h13 * cx[xcell + 4] * cx[xcell + 6]; 			// Her1-Her13 dimer association
	*a0 += ca[1] + ca[2] + ca[4] + ca[6] - old;
}

inline void update_a3_18 (double* ca, parameters* p, int* cx, int xcell, double* a0) { // If Her1-Her1 dimer is updated
	double old = ca[3] + ca[18];
	ca[3] = p->ddh1h1 * cx[xcell + 8]; 		// Her1-Her1 dimer dissociation
	ca[18] = p->ddgh1h1 * cx[xcell + 8]; 	// Her1-Her1 dimer degradation
	*a0 += ca[3] + ca[18] - old;
}

inline void update_a4_9_10_12 (double* ca, parameters* p, int* cx, int xcell, double* a0) { // If Her7 protein is updated
	double old = ca[4] + ca[9] + ca[10] + ca[12];
	ca[4] = p->dah1h7 * cx[xcell + 4] * cx[xcell + 5]; 				// Her1-Her7 dimer association
	ca[9] = p->pdh7 * cx[xcell + 5]; 								// Her7 protein degradation
	ca[10] = p->dah7h7 * cx[xcell + 5] * (cx[xcell + 5] - 1) / 2; 	// Her7-Her7 dimer association
	ca[12] = p->dah7h13 * cx[xcell + 5] * cx[xcell + 6]; 			// Her7-Her13 dimer association
	*a0 += ca[4] + ca[9] + ca[10] + ca[12] - old;
}

inline void update_a5_19 (double* ca, parameters* p, int* cx, int xcell, double* a0) { // If Her1-Her7 dimer is updated
	double old = ca[5] + ca[19];
	ca[5] = p->ddh1h7 * cx[xcell + 9];		// Her1-Her7 dimer dissociation
	ca[19] = p->ddgh1h7 * cx[xcell + 9];	// Her1-Her7 dimer degradation
	*a0 += ca[5] + ca[19] - old;
}

inline void update_a6_12_15_16 (double* ca, parameters* p, int* cx, int xcell, double* a0) { // If Her13 protein is updated
	double old = ca[6] + ca[12] + ca[15] + ca[16];
	ca[6] = p->dah1h13 * cx[xcell + 4] * cx[xcell + 6];					// Her1-Her13 dimer association
	ca[12] = p->dah7h13 * cx[xcell + 5] * cx[xcell + 6];				// Her7-Her13 dimer association
	ca[15] = p->pdh13 * cx[xcell + 6];									// Her13 protein degradation
	ca[16] = p->dah13h13 * cx[xcell + 6] * (cx[xcell + 6] - 1) / 2;		// Her13-Her13 dimer association
	*a0 += ca[6] + ca[12] + ca[15] + ca[16] - old;
}

inline void update_a7_20 (double* ca, parameters* p, int* cx, int xcell, double* a0) { // If Her1-Her13 dimer is updated
	double old = ca[7] + ca[20];
	ca[7] = p->ddh1h13 * cx[xcell + 10];		// Her1-Her13 dimer dissociation
	ca[20] = p->ddgh1h13 * cx[xcell + 10];		// Her1-Her13 dimer degradation
	*a0 += ca[7] + ca[20] - old;
}

inline void update_a8_29 (double* ca, parameters* p, int* cx, int xcell, double* a0) { // If her7 mRNA is updated
	double old = ca[8] + ca[29];
	ca[8] = p->psh7 * cx[xcell + 1];			// Her7 protein synthesis
	ca[29] = p->mdh7 * cx[xcell + 1];			// her7 mRNA degradation
	*a0 += ca[8] + ca[29] - old;
}

inline void update_a11_21 (double* ca, parameters* p, int* cx, int xcell, double* a0) { // If Her7-Her7 dimer is updated
	double old = ca[11] + ca[21];
	ca[11] = p->ddh7h7 * cx[xcell + 11];		// Her7-Her7 dimer dissociation
	ca[21] = p->ddgh7h7 * cx[xcell + 11];		// Her7-Her7 dimer degradation
	*a0 += ca[11] + ca[21] - old;				
}

inline void update_a13_22 (double* ca, parameters* p, int* cx, int xcell, double* a0) { // If Her7-Her13 dimer is updated
	double old = ca[13] + ca[22];
	ca[13] = p->ddh7h13 * cx[xcell + 12];			// Her7-Her13 dimer dissociation
	ca[22] = p->ddgh7h13 * cx[xcell + 12];			// Her7-Her13 dimer degradation
	*a0 += ca[13] + ca[22] - old;
}

inline void update_a14_31 (double* ca, parameters* p, int* cx, int xcell, double* a0) { // If her13 mRNA is updated
	double old = ca[14] + ca[31];
	ca[14] = p->psh13 * cx[xcell + 2];			// Her13 protein synthesis
	ca[31] = p->mdh13 * cx[xcell + 2];			// her13 mRNA degradation
	*a0 += ca[14] + ca[31] - old;
}

inline void update_a17_23 (double* ca, parameters* p, int* cx, int xcell, double* a0) {	// If Her13-Her13 dimer is updated
	double old = ca[17] + ca[23];
	ca[17] = p->ddh13h13 * cx[xcell + 13];		// Her13-Her13 dimer dissociation
	ca[23] = p->ddgh13h13 * cx[xcell + 13];		// Her13-Her13 dimer degradation
	*a0 += ca[17] + ca[23] - old;
}

inline void update_a24_33 (double* ca, parameters* p, int* cx, int xcell, double* a0) { // If delta mRNA is updated
	double old = ca[24] + ca[33];
	ca[24] = p->psd * cx[xcell + 3];			// Delta protein synthesis
	ca[33] = p->mdd * cx[xcell + 3];			// delta mRNA degradation
	*a0 += ca[24] + ca[33] - old;
}

inline void update_a25 (double* ca, parameters* p, int* cx, int xcell, double* a0) { // If delta protein is updated
	double old = ca[25];
	ca[25] = p->pdd * cx[xcell + 7];			// Delta protein degradation
	*a0 += ca[25] - old;
}

inline void update_a26_28_32 (double* ca, parameters* p, int* cx, int xcell, int neighbors, int* nc, double* a0) {	// If Her1-Her1 dimer, Her7-Her13 dimer and Delta protein is updated
	double old = ca[26] + ca[28] + ca[32];
	double x11 = cx[xcell + 8] / p->critph1h1;
	double x713 = cx[xcell + 12] / p->critph7h13;
	double reaction_sum = 0;
	for (int nn = 1; nn < neighbors; nn++) {
		reaction_sum += cx[nc[nn] * species + 7];
	}
	double y = (reaction_sum / (neighbors - 1)) / p->critpd;
	double result1 = 1 + x11 * x11 + x713 * x713;
	double result2 = (1 + y) / (y + result1);
	ca[26] = p->msh1 * result2;
	ca[28] = p->msh7 * result2;
	ca[32] = p->msd / result1;
	*a0 += ca[26] + ca[28] + ca[32] - old;
}

#endif

