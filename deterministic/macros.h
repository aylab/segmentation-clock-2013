/*
 Deterministic simulator for the zebrafish segmentation clock.
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

#ifndef MACROS_H
#define MACROS_H

// mRNA synthesis rates
#define RMSH1			0
#define RMSH7			1
#define RMSH13			2
#define RMSDELTA		3

// mRNA degradation rates
#define RMDH1			4
#define RMDH7			5
#define RMDH13			6
#define RMDDELTA		7

// Protein synthesis rates
#define RPSH1			8
#define RPSH7			9
#define RPSH13			10
#define RPSDELTA		11

// Protein degradation rates
#define RPDH1			12
#define RPDH7			13
#define RPDH13			14
#define RPDDELTA		15

// Dimer association rates
#define RDAH1H1			16
#define RDAH1H7			17
#define RDAH1H13		18
#define RDAH7H7			19
#define RDAH7H13		20
#define RDAH13H13		21

// Dimer dissociation rates
#define RDDIH1H1		22
#define RDDIH1H7		23
#define RDDIH1H13		24
#define RDDIH7H7		25
#define RDDIH7H13		26
#define RDDIH13H13		27

// Dimer degradation rates
#define RDDGH1H1		28
#define RDDGH1H7		29
#define RDDGH1H13		30
#define RDDGH7H7		31
#define RDDGH7H13		32
#define RDDGH13H13		33

// mRNA transcription delays
#define RDELAYMH1		34
#define RDELAYMH7		35
#define RDELAYMH13		36
#define RDELAYMDELTA	37

// Protein translation delays
#define RDELAYPH1		38
#define RDELAYPH7		39
#define RDELAYPH13		40
#define RDELAYPDELTA	41

// Critical number of molecules of proteins per cell for inhibition of transcription
#define RCRITPH1H1 		42
#define RCRITPH7H13		43

#define RCRITPDELTA		44

#define NUM_RATES		45 // How big an array holding rates must be
#define MIN_DELAY		34 // The smallest index referring to a delay
#define MAX_DELAY		41 // The largest index referring to a delay

#endif