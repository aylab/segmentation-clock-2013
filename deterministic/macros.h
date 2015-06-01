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

// Protein synthesis rates
#define RPSH1			0
#define RPSH7			1
#define RPSH13			2
#define RPSDELTA		3

// Protein degradation rates
#define RPDH1			4
#define RPDH7			5
#define RPDH13			6
#define RPDDELTA		7

// mRNA synthesis rates
#define RMSH1			8
#define RMSH7			9
#define RMSH13			10
#define RMSDELTA		11

// mRNA degradation rates
#define RMDH1			12
#define RMDH7			13
#define RMDH13			14
#define RMDDELTA		15

// Dimer degradation rates
#define RDDGH1H1		16
#define RDDGH1H7		17
#define RDDGH1H13		18
#define RDDGH7H7		19
#define RDDGH7H13		20
#define RDDGH13H13		21

// mRNA transcription delays
#define RDELAYMH1		22
#define RDELAYMH7		23
#define RDELAYMH13		24
#define RDELAYMDELTA	25

// Protein translation delays
#define RDELAYPH1		26
#define RDELAYPH7		27
#define RDELAYPH13		28
#define RDELAYPDELTA	29

// Dimer association rates
#define RDAH1H1			30
#define RDAH1H7			32
#define RDAH1H13		34
#define RDAH7H7			36
#define RDAH7H13		38
#define RDAH13H13		40

// Dimer dissociation rates
#define RDDIH1H1		31
#define RDDIH1H7		33
#define RDDIH1H13		35
#define RDDIH7H7		37
#define RDDIH7H13		39
#define RDDIH13H13		41


// Critical number of molecules of proteins per cell for inhibition of transcription
#define RCRITPH1H1 		42
#define RCRITPH7H13		43

#define RCRITPDELTA		44

#define NUM_RATES		45 // How big an array holding rates must be
#define MIN_DELAY		34 // The smallest index referring to a delay
#define MAX_DELAY		41 // The largest index referring to a delay

#endif