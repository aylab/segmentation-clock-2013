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

#ifndef MACROS_H
#define MACROS_H

#define beta 0.05 // increase this to merge more delayed queues for id-leaping at the cost of less accuracy (should be 0.05-0.1)
#define delta_factor 0.05 // increase this to make partial equilibrium a more stringent condition, allowing more reactions to be considered for implicit tau (should be around 0.05)
#define epsilon 0.01 // increase this to increase the timesteps of tau-leaping (should be 0.03-0.05)
#define MB pow(2, 20) // the number of bytes in a megabyte
#define ncrit 10 // the maximum number of molecules / update value of a species for it to be considered critical (should be around 10)
#define structure_twocell 0 // indicates a two-cell simulation
#define structure_chain 1 // indicates a chain simulation
#define structure_tissue 2 // indicates a tissue simulation
#define neighbors_for_twocell 2 // how many neighbors a cell in a two-cell systems has (including itself)
#define neighbors_for_chain 3 // how many neighbors a cell in a chain has (including itself)
#define neighbors_for_tissue 7 // how many neighbors a cell in a haxagonal tissue grid has (including itself)
#define num_of_delayed_reactions 7 // the number of reactions that are delayed in the zebrafish segmentation system
#define num_of_parameters 45 // the number of parameters each line in the input file should have
#define nstiff 100 // increase this to increase the rate of implicit tau-leaping over explicit tau-leaping (should be around 100)
#define reactions 34 // the number of reactions in the zebrafish segmentation system
#define skip_steps_ex 100 // how many steps of the next-reaction-method to perform before resuming tau-leaping if the last tau-leap was explicit
#define skip_steps_im 10 // how many steps of the next-reaction-method to perform before resuming tau-leaping if the last tau-leap was implicit
#define species 14 // the number of species in the zebrafish segmentation system
#define tau1_mult 10 // the factor used (in combination with a0) to determine if tau-leaping is efficient enough to use over the next-reaction-method

// these are convenient macros to calculate the maximum and minimum of two given values and the absolute value of the given value
// don't give compound expressions as arguments to these because macros do stupid things with them (e.g. don't put max(i++, j) because then i will be incremented more times than expected)
/*
#define min(x, y) x < y ? x : y
#define max(x, y) x > y ? x : y
#define abs(x) x < 0 ? -x : x
*/
// escape codes to color the terminal output and shortcuts for common outputs (set -c or --no-color to disable these)
#define terminal_blue_d "\x1b[34m"
#define terminal_red_d "\x1b[31m"
#define terminal_reset_d "\x1b[0m"
#define terminal_done terminal_blue << "Done" << terminal_reset
#define terminal_no_memory terminal_red << "Not enough memory!" << terminal_reset

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

