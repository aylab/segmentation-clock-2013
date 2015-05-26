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
#define min(x, y) x < y ? x : y
#define max(x, y) x > y ? x : y
#define abs(x) x < 0 ? -x : x

// escape codes to color the terminal output and shortcuts for common outputs (set -c or --no-color to disable these)
#define terminal_blue_d "\x1b[34m"
#define terminal_red_d "\x1b[31m"
#define terminal_reset_d "\x1b[0m"
#define terminal_done terminal_blue << "Done" << terminal_reset
#define terminal_no_memory terminal_red << "Not enough memory!" << terminal_reset

#endif

