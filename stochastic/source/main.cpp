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

/*
We used the following articles and their associated algorithms to construct this simulation:
"A modified next reaction method for simulating chemical systems with time dependent propensities and delays" by Anderson, J. of Chemical Physics, 2007
"Adaptive explicit-implicit tau-leaping method with automatic tau selection" by Cao et al., J. of Chemical Physics, 2007
"D-leaping: Accelerating stochastic simulation algorithms for reactions with delays" by Bayati et al., J. of Computational Physics, 2009
"Improved delay-leaping simulation algorithm for biochemical reaction systems with delays" by Ni et al., J. of Chemical Physics, 2012
*/

#include "main.h"

#include <errno.h>
#include <fstream>
#include <list>
#include <map>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <sys/stat.h>

#include "file-io.h"
#include "macros.h"
#include "parameters.h"
#include "rand-dist.h"
#include "rq-node.h"
#include "updates.h"
#include "utility.h"

using namespace std;

// global variables used in main.cpp and file-io.cpp
char* terminal_blue;
char* terminal_red;
char* terminal_reset;


int main (int argc, char** argv) {
	/*
	Initialize constant variables:
	These arrays and matrices do not change over the course of the simulation and contain information specific to the zebrafish segmentation notch pathway.
	*/
	
	const int delayed_reactions[num_of_delayed_reactions] = {0, 8, 14, 24, 26, 28, 32}; // list of delayed reaction indices
	
	const int initial_concentrations[species] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // initial (time=0) concentration values for each species
	
	const int species_update_indices[reactions][4] = {{0, 0, 0, 0}, {1, 4, 0, 0}, {2, 4, 8, 0}, {2, 4, 8, 0}, {3, 4, 5, 9}, {3, 4, 5, 9}, {3, 4, 6, 10}, {3, 4, 6, 10}, {0, 0, 0, 0}, {1, 5, 0, 0}, {2, 5, 11, 0}, {2, 5, 11, 0}, {3, 5, 6, 12}, {3, 5, 6, 12}, {0, 0, 0, 0}, {1, 6, 0, 0}, {2, 6, 13, 0}, {2, 6, 13, 0}, {1, 8, 0, 0}, {1, 9, 0, 0}, {1, 10, 0, 0}, {1, 11, 0, 0}, {1, 12, 0, 0}, {1, 13, 0, 0}, {0, 0, 0, 0}, {1, 7, 0, 0}, {0, 0, 0, 0}, {1, 0, 0, 0}, {0, 0, 0, 0}, {1, 1, 0, 0}, {1, 2, 0, 0}, {1, 2, 0, 0}, {0, 0, 0, 0}, {1, 3, 0, 0}}; // the species indices to update for each nChains of cells look like this:on-delayed reaction (the first value of each array is the number of species to update, so if the first value is 2 then only the next 2 indices are looked at to determine which species to update)
	
	const int species_update_indices_delayed[num_of_delayed_reactions] = {4, 5, 6, 7, 0, 1, 3}; // the species indices to update for each delayed reaction (in this system delayed reactions update only 1 species each so the structure is less complicated than the one above for non-delayed reactions)
	
	const int species_update_values[species][reactions] = // the values (i.e. number of molecules) to change each species by for each reaction
	//	  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33
	       {{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, +1, -1,  0,  0,  0,  0,  0,  0},
		{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, +1, -1,  0,  0,  0,  0},
		{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, +1, -1,  0,  0},
		{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, +1, -1},
		{+1, -1, -2, +2, -1, +1, -1, +1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
		{ 0,  0,  0,  0, -1, +1,  0,  0, +1, -1, -2, +2, -1, +1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
		{ 0,  0,  0,  0,  0,  0, -1, +1,  0,  0,  0,  0, -1, +1, +1, -1, -2, +2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
		{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, +1, -1,  0,  0,  0,  0,  0,  0,  0,  0},
		{ 0,  0, +1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
		{ 0,  0,  0,  0, +1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
		{ 0,  0,  0,  0,  0,  0, +1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
		{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, +1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
		{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, +1, -1,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
		{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, +1, -1,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}};
	
	const int par_eq_pairs[reactions] = {-1, -1, 3, 2, 5, 4, 7, 6, -1, -1, 11, 10, 13, 12, -1, -1, 17, 16, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}; // if a reaction has a pair then the value at its index contains its pair's index, and if a reaction doesn't have a pair then the value at its index is -1
	
	map<string, int> levels; // a map between strings describing concentration levels and their indices used in the program
	levels["her1"] = 1;
	levels["her7"] = 2;
	levels["her13"] = 3;
	levels["delta"] = 4;
	levels["Her1"] = 5;
	levels["Her7"] = 6;
	levels["Her13"] = 7;
	levels["Delta"] = 8;
	levels["Her1Her1"] = 9;
	levels["Her1Her7"] = 10;
	levels["Her1Her13"] = 11;
	levels["Her7Her7"] = 12;
	levels["Her7Her13"] = 13;
	levels["Her13Her13"] = 14;
	
	/*
	Check and accept user-given arguments:
	This program can accept a number of different arguments from the command-line.
	Each argument follows the form "-abbreviated_option value" or "--option value" where "-abbreviated" is a single letter preceded by a hyphen, "--option" is a full word preceded by two hyphens, and value is the value to give to the option.
	The abbreviated and full option strings are entirely interchangeable; one is for brevity and the other for readability. Every argument that is omitted from the command-line is given a default value below.
	Note that if you enter "-h" or "-l" (help and licensing information, respectively) no simulations will run. Also note that for both "-h" and "-l", no other arguments should be given afterwards.
	*/
	
	unsigned int xcells = 2; // the tissue width in cells (-x or --width changes this)
	unsigned int ycells = 1; // the tissue height in cells (-y or --height changes this)
	unsigned int max_minutes = 1200; // the maximum number of minutes the simulation should run before being terminated (-m or --minutes changes this)
	unsigned long max_timesteps = 1000000000000L; // the maximum number of timesteps the simulation should run before being terminated (-t or --timesteps changes this)
	unsigned int runs = 1; // the number of times to run the simulation (-r or --runs changes this)
	unsigned int seed = time(NULL); // the seed used to generate random numbers (using the same seed will produce the same results each time) (-s or --seed changes this)
	char* input_file = NULL; // the path and filename containing the simulation parameters (absolute and relative paths work, this defaults to "input.txt") (-i or --input changes this)
	char* output_path = NULL; // the directory path to print the results to (absolute and relative paths work, this defaults to "output") (-o or --output changes this)
	char* gradients_file = NULL; // the path and filename containing the gradient (absolute and relative paths work, this defaults to NULL) (-i or --input changes this)
	unsigned int con_level = 0; // the index of the concentration level to print as output
	double granularity = 0.1; // the amount of time to skip between each timestep when printing the output file (-g or --granularity changes this)
	unsigned int print_interval = 1200; // how often (in minutes) the output file should be printed to (this does not change the simulation results, but is useful to see progress in very slow simulations) (-p or --print changes this)
	char* seed_file = NULL; // the filename containing the seed used to generate random numbers (relative to the output path, this defaults to "seed.txt") (-k --keepseed changes this)
	bool appx = false; // if this is set to true then the simulation will use approximation algorithms to create faster but potentially less accurate results (-a or --algorithm followed by "exact" or "appx" changes this)
	
	terminal_color();

	checkArgs(argc, argv, xcells, ycells, max_minutes, max_timesteps, runs, seed, &input_file, &output_path, con_level, levels, granularity, print_interval, &seed_file, appx);
	
	unsigned int cells = xcells * ycells; // the total number of cells
	int structure; // two-cell, chain, or tissue
	int neighbors; // the number of neighbors each cell has
	checkSize(xcells, ycells, input_file, output_path, structure, neighbors);
	
	int nc[cells][neighbors];
	cells_neighbors(structure, neighbors, cells, xcells, &nc);
	
	/*
	Initialize the seed used to generate random numbers:
	1) Store the seed in the seed file
	2) initialize the random generator with the seed
	*/
	
	if (seed_file != NULL) {
		ofstream ofile_seed;
		ofile_seed.open(seed_file, fstream::out);
		ofile_seed << seed << endl;
		ofile_seed.close();
	}
	srand(seed);
	
	/*
	File input:
	1) Read the input file
	2) Parse the first line of the file (mechanisms to parse more are available, but not used, 
	   since our system uses a bash script to split parameter files into single lines 
	   so a cluster list system can assign different parameter sets to different processors with more control/granularity)
	3) Store each comma-separated value in the globally-used parameter struct
	4) Create an array of delay times, each index corresponding to each delayed reaction
	*/
	
	cout << terminal_blue << "Reading " << terminal_reset << input_file << " ... ";
	char* all_parameters = NULL; // the buffer to contain the file as a string
	read_file(input_file, &all_parameters); // read the file into the all_parameters buffer
	double par_array[num_of_parameters]; // create an array for each parameter
	int if_index = 0; // this keeps track of the most recently read index in the file so multiple lines can be read
	parse_line(all_parameters, par_array, &if_index); // read the first line into the par_array array
	parameters pars(par_array); // create the parameters struct using par_array
	double delay_times[num_of_delayed_reactions] = {pars.delayph1, pars.delayph7, pars.delayph13, pars.delaypd, pars.delaymh1, pars.delaymh7, pars.delaymd}; // create the delay times array
	free(all_parameters);
	cout << terminal_done << endl;
	
	/*
	File output:
	1) Created the output directory if necessasry
	2) Create an array of file streams, one file per run
	3) Print the cell width and height in the first line of each file (this is used later by data-processing scripts)
	*/
	
	int path_length = strlen(output_path); // get the path length and remove the trailing slash from the path if it was given with one
	if (output_path[path_length - 1] == '/') {
		output_path[--path_length] = '\0';
	}
	cout << terminal_blue << "Creating " << terminal_reset << output_path << " directory if necessary ... ";
	if (mkdir(output_path, 0755) != 0 && errno != EEXIST) { // make the output directory if necessary or exit the program if there was an error trying to do so
		cout << terminal_red << "Couldn't create " << output_path << " directory!" << terminal_reset << endl;
		exit(1);
	}
	cout << terminal_done << endl;
	
	ofstream ofile[runs]; // the array of file streams
	for (unsigned int r = 0; r < runs; r++) { // for every run that will be simulated
		int rsize = r != 0 ? log10(r) + 1 : 1; // run numbers 0-9 take one byte to write, run numbers 10-99 take two, etc.
		char output_file[path_length + rsize + 9]; // the buffer containing the path and file name for the run output
		memcpy(output_file, output_path, path_length);
		output_file[path_length] = '/';
		memcpy(output_file + path_length + 1, "run" , 3);
		output_file[path_length + rsize + 3] = r % 10 + '0';
		int rdig = r / 10;
		for (int i = rsize - 2; i >= 0; i--) { // write each digit of the run number
			output_file[path_length + i + 4] = rdig + '0';
			rdig /= 10;
		}
		strcpy(output_file + path_length + rsize + 4, ".txt");
		cout << terminal_blue << "Creating " << terminal_reset << output_file << " ... ";
		try { // try to create the file or exit the program if there was an error trying to do so
			ofile[r].open(output_file, fstream::out);
			ofile[r] << xcells << " " << ycells << endl; // write the tab-separated cell width and height to the file on the first line
		} catch (ofstream::failure) {
			cout << terminal_red << "Couldn't create " << output_file << "!" << terminal_reset << endl;
			exit(1);
		}
		cout << terminal_done << endl;
	}	
	
	int chunk = max_minutes / granularity + 1;
	int** x; // concentrations
	double* T; // timesteps
	memory_alloc(chunk, cells, &x, &T);

	/*
	Run the simulations:
	1) For every run (by default 1, but can be changed with -r or --runs), do the following:
	2) Initialize the initial concentration and timestep values
	3) Initialize delayed reaction lists, propensity values, Pk and Tk arrays (for calculating delta using next-reaction-method), etc.
	4) Iterate through the allowed number of iterations doing the following:
		a) End the simulation if the number of simulation minutes exceeds the limit
		b) Otherwise, if tau-leaping is not temporarily disabled, try to leap
		c) If tau-leaping is disabled, run next-reaction-method until it isn't anymore
		d) If the simulation has progressed enough to print (by default every 60 minutes, but can be changed with -p or --print), print the results
	5) When the simulation is done, clear any remaining items in the delayed reaction lists
	6) Print any unprinted results
	7) If there is another run then go back to step 4
	8) Otherwise, if every run is finished, clear the memory allocated for the concentrations and timestep structures
	*/
	
	for (unsigned int r = 0; r < runs; r++) {
		// indicate to the user that a new run is starting
		cout << terminal_blue << "Simulating " << terminal_reset << "run #" << r << " ... ";
		cout.flush();

		// initialize concentration and timestep values (the 0th index isn't used because T[chunk_index - 1] must always exist, so the results start at 1)
		for (unsigned int i = 0; i < cells; i++) {
			for (int j = 0; j < species; j++) {
				x[1][i * species + j] = initial_concentrations[j];
			}
		}
		T[0] = 0; // placeholder value
		T[1] = 0; // actual first time entry
		int last_print = T[1]; // the last time that was printed

		// create an array of lists to store the delayed reactions
		list<double> rq[cells][num_of_delayed_reactions]; // used for next-reaction-method and id-leaping
		list<rq_node> rq_idl[cells][num_of_delayed_reactions]; // used just for id-leaping
		/* create Tk and Pk arrays, used to calculate delta for the next-reaction-method 
		   (Tk = current internal time of reaction k, Pk = first internal time after Tk at which reaction k fires)
		*/
		double Tk[cells][reactions];
		double Pk[cells][reactions];
		// create the propensity values array and the sum of them, a0
		double a[cells][reactions];
		double a0 = 0;
		// create an array to keep track of how many times each reaction fires for tau-leaping
		int firings[cells][reactions];
		
		// initialize everything that was just created
		for (unsigned int i = 0; i < cells; i++) {
			for (int k = 0; k < reactions; k++) {
				Tk[i][k] = 0;
				Pk[i][k] = pk_dist(); // log(1 / unif_dist())
				a[i][k] = 0;
				firings[i][k] = 0;
			}
			
			// initialize delayed transcription reactions properly
			a[i][26] = pars.msh1;
			a[i][28] = pars.msh7;
			a[i][30] = pars.msh13; // remains constant throughout the simulation
			a[i][32] = pars.msd;
			a0 += a[i][26] + a[i][28] + a[i][30] + a[i][32]; // keep track of the sum of the propensities as they are changed
		}
		
		int HOR[3]; // the higher order reactions array stores a dividing factor for Hill-type reactions, used for tau-leaping
		for (int j = 0; j < 3; j++) {
			HOR[j] = 0;
		}
		
		bool last_step_ex = true; // whether or not the last tau-leaping step was explicit (as opposed to implicit)
		int skip_steps = !appx * -1; // if the simulation shouldn't approximate then set skip_steps to -1, preventing it from ever becoming 0

		// iterate through each time step
		int cell_index = -1; // the cell index of the most recently fired reaction for the next-reaction-method
		int reaction_index = -1; // the reaction index of the most recently fired reaction
		bool is_delayed = false; // whether or not the most recently fired reaction was a delayed reaction
		int chunk_index = 1; // the current chunk index in which to store data
		int* cx = x[chunk_index]; // the current concentrations array in which to store data
		
		// iterate through the maximum allowed number of timesteps
		for (unsigned long iter = 0; iter < max_timesteps; iter++) {
			// if the current simulation time is more than the maximum simulation time then end the simulation
			if (T[chunk_index] >= max_minutes) {
				break;
			}
			
			if (skip_steps == 0) { // if tau-leaping isn't temporarily disabled
				
				/*
				Adaptive tau-leaping with improved delay leaping (id-leaping):
				1) Calculate non-critical tau candidate (tau1)
					a) Calculate explicit tau1 candidate (tau_ex) based on all non-critical reactions
					b) Calculate implicit tau1 candidate (tau_im) based on non-critical non-partial-equilibrium reactions
					c) Take tau1 as the minimum of tau_ex and tau_im
				2) If tau1 is so small that the next-reaction-method would be more efficient, stop trying to tau-leap
				3) Otherwise, if tau-leaping is useful, calculate tau2
					a) Store the sum of the propensities of the critical reactions in a0_crit
					b) Calculate tau2 based on an exponential distribution of a0_crit
				4) Take tau as the minimum of tau1 and tau2
				5) Fire any delayed reactions that should be fired and remove them from their list
				6) Calculate the number of firings for each reaction
				7) If the new concentration levels would be negative, don't update the concentrations and instead halve tau1 and try from step 3 again
				8) Add any appropriate delayed reactions (or merge them with existing ones via id-leaping method)
				9) Update concentrations, propensities, and the simulation timestep
				*/
				
				double tau1; // non-critical tau candidate
				double tau_ex = INFINITY; // explicit tau1 candidate
				double tau_im = INFINITY; // implicit tau1 candidate
				bool critical[cells][reactions]; // a matrix that stores whether each reaction of each cell is critical or not at this point in the simulation
				bool ncrit_exists = false; // true if at least one non-critical reaction exists, false otherwise
				for (unsigned int i = 0; i < cells; i++) {
					int co = i * species; // cell offset for indexing the concentrations array
					for (int k = 0; k < reactions; k++) {
						if (a[i][k] > 0) { // if the reaction has a chance of occuring
							int min = (int)-INFINITY; // since this considers only negative update reactions, the least negative value is chosen
							for (int j = 0; j < species; j++) {
								if (species_update_values[j][k] < 0) { // consider only negative update reactions (i.e. degradations)
									int min_temp = cx[co + j] / species_update_values[j][k];
									if (min_temp > min) { // if min_temp is less negative than min then overwrite min with min_temp
										min = min_temp;
									}
								}
							}
							critical[i][k] = min != (int)-INFINITY ? min * -1 < ncrit : false; // each reaction is considered critical if the absolute value of min is less than ncrit
						} else { // reactions with no chance of occuring can't be critical
							critical[i][k] = false;
						}
						if (!ncrit_exists && !critical[i][k]) { // if the first non-critical reaction is encountered then indicate that it exists
							ncrit_exists = true;
						}
					}
					
					if (ncrit_exists) { // if at least one non-critical reaction exists
						for (int j = 4; j <= 6; j++) { // species 4, 5, and 6 are higher order reactions, so their division factors must be calculated based on their current concentration values
							HOR[j - 4] = 2 + 1 / (cx[co + j] + 1);
						}
						
						for (int j = 0; j < species; j++) {
							double exg = epsilon * cx[co + j]; // epsilon is usually 0.03-0.05
							if (j >= 4 && j <= 6) { // higher order reactions divide by a factor more than 1, lower order reactions divide by 1, so no division is executed
								exg /= HOR[j - 4];
							}
							double max_ep = max(exg, 1); // make sure exg is at least one
							
							double mu_ex = 0; // mu for explicit tau
							double mu_im = 0; // mu for implicit tau
							double sigma_ex = 0; // sigma for explicit tau
							double sigma_im = 0; // sigma for implicit tau
							for (int k = 0; k < reactions; k++) {
								if (!critical[i][k]) { // mu and sigma values only consider non-critical reactions
									int r_up = species_update_values[j][k]; // reaction update value
									double r_a = a[i][k]; // reaction propensity value
									double mu_change = r_up * r_a; // add to explicit tau
									mu_ex += mu_change;
									sigma_ex += mu_change * r_up;
									int pair = par_eq_pairs[k]; // the index of the reaction's pair, or -1 if it has none
									if (pair == -1) { // a reaction without a pair is never in partial equilibrium
										mu_im += mu_change;
										sigma_im += mu_change * r_up;
									} else if (pair > k) { // if the reaction has a pair and the pair has not already been considered
										int pair_up = species_update_values[j][pair]; // pair update value
										double pair_a = a[i][pair]; // pair propensity value
										if ((r_a < pair_a && pair_a - r_a <= delta_factor * r_a) || (pair_a <= r_a && r_a - pair_a <= delta_factor * pair_a)) { // if the pair isn't in partial equilibrium then add to implicit tau
											mu_im += mu_change + pair_up * pair_a;
											sigma_im += mu_change * r_up + pair_up * pair_up * pair_a;
										}
									}
								}
							}
							mu_ex = abs(mu_ex);
							mu_im = abs(mu_im);
							
							// calculate new candidates for tau explicit and tau implicit and update them if a new minimum value has been found
							double min1 = mu_ex != 0 ? max_ep / mu_ex : INFINITY;
							double min2 = sigma_ex != 0 ? (max_ep * max_ep) / sigma_ex : INFINITY;
							min1 = min(min1, min2);
							if (min1 < tau_ex) {
								tau_ex = min1;
							}
							min1 = mu_im != 0 ? max_ep / mu_im : INFINITY;
							min2 = sigma_im != 0 ? (max_ep * max_ep) / sigma_im : INFINITY;
							min1 = min(min1, min2);
							if (min1 < tau_im) {
								tau_im = min1;
							}
						}
					}
				}
				
				bool temp_lse = last_step_ex; // store whether the last step was explicit
				if (tau_im > nstiff * tau_ex) { // if the system is considered stiff then set tau1 to tau implicit and mark the current step as implicit
					tau1 = tau_im;
					last_step_ex = false;
				} else { // if the system isn't considered stiff then set tau1 to tau explicit and mark the current step as explicit
					tau1 = tau_ex;
					last_step_ex = true;
				}
				
				bool repeat;
				do { // under some conditions these steps must be repeated multiple times per iteration but by default they aren't (this only repeats when repeat=true)
					repeat = false;
					if (tau1 < tau1_mult / a0) { // if the stepsize of tau1 is so small that it would be more efficient to use the next-reaction-method then skip tau-leaping for a number of steps depending on the stiffness of the system
						if (temp_lse) {
							skip_steps = skip_steps_ex;
						} else {
							skip_steps = skip_steps_im;
						}
					} else { // if tau-leaping is more efficient than the next-reaction-method then continue
						// calculate the sum of the propensities for all critical reactions
						double a0_crit = 0;
						for (unsigned int i = 0; i < cells; i++) {
							for (int k = 0; k < reactions; k++) {
								if (critical[i][k]) {
									a0_crit += a[i][k];
								}
							}
						}
						
						double tau2 = a0_crit != 0 ? expo_dist(a0_crit) : INFINITY; // critical tau candidate
						double tau = min(tau1, tau2); // pick tau based on the smaller of the non-critical and critical candidates
						
						double nT = T[chunk_index] + tau; // the next simulation timestep
						for (unsigned int i = 0; i < cells; i++) { // iterate through every delayed reaction for every cell
							for (int d = 0; d < num_of_delayed_reactions; d++) {
								int ri = delayed_reactions[d];
								list<double>::iterator earliest = rq[i][d].begin(); // iterate through this reaction's list, split into the earlist firing time list and the number of firings and span length list
								list<rq_node>::iterator fs = rq_idl[i][d].begin();
								while (earliest != rq[i][d].end()) {
									if (*earliest < nT) { // if the delayed reaction's delay is over
										double qd = nT - *earliest; // the difference between the delayed reaction's earliest firing time and the new time
										int kd = bino_dist(fs->firings, min(qd, fs->span) / fs->span); // a binomial distribution based on the number of times the reaction should fire and the minimum of qd and the reaction's span
										fs->firings -= kd; // update firings, span, and earliest
										fs->span -= qd;
										*earliest = nT;
										for (int j = 0; j < species; j++) { // update the concentrations based on kd
											cx[i * species + j] += kd * species_update_values[j][ri];
										}
										if (fs->firings == 0) { // if the delayed reaction is done then remove it from its list
											rq[i][d].erase(earliest++);
											rq_idl[i][d].erase(fs++);
										}
									} else {
										earliest++;
										fs++;
									}
								}
							}
						}
						
						if (tau2 > tau1) { // if tau=tau1
							for (unsigned int i = 0; i < cells; i++) {
								for (int k = 0; k < reactions; k++) {
									if (critical[i][k]) {
										firings[i][k] = 0; // no critical reactions should fire
									} else {
										firings[i][k] = pois_dist(a[i][k] * tau); // non-critical reactions should fire based on a Poisson distribution of their propensities and the time difference, tau
									}
								}
							}
						} else { // if tau=tau2
							unsigned int jc = 0; // the index of the critical reaction that should fire this iteration
							double random_start = 0; // use a point probability distribution to pick jc
							double random_end = unif_dist();
							bool breakout = false;
							for (unsigned int i = 0; i < cells; i++) {
								if (breakout) {
									break;
								}
								for (int k = 0; k < reactions; k++) {
									if (critical[i][k]) {
										random_start += a[i][k] / a0_crit; // incrementally work up to a0_crit
										if (random_start >= random_end) { // if the intermediate sum is more than the uniformly distributed random variable chosen above then set jc to the kth reaction of the ith cell
											jc = i * reactions + k;
											breakout = true; // stop iterating because jc has already been found
											break;
										}
									}
								}
							}
							
							for (unsigned int i = 0; i < cells; i++) {
								for (int k = 0; k < reactions; k++) {
									if (critical[i][k]) {
										firings[i][k] = (i * reactions + k == jc); // no critical reaction but jc should fire
									} else {
										firings[i][k] = pois_dist(a[i][k] * tau); // non-critical reactions should fire based on a Poisson distribution of their propensities and the time difference, tau
									}
									
									for (int j = 0; j < species; j++) {
										if (cx[i * species + j] + firings[i][k] * species_update_values[j][k] < 0) { // if the new concentration levels for any species would be negative then halve tau1 and try again
											repeat = true;
											tau1 /= 2;
											break;
										}
									}
									if (repeat) {
										break;
									}
								}
							}
						}
						
						if (!repeat) { // if the concentrations won't be negative (i.e. everything is fine)
							for (unsigned int i = 0; i < cells; i++) {
								for (int d = 0; d < num_of_delayed_reactions; d++) {
									int ri = delayed_reactions[d]; // reaction index
									if (firings[i][ri] > 0) { // if the delayed reaction should fire
										bool add = true;
										if (!rq_idl[i][d].empty()) { // if the list for this reaction isn't empty
											rq_node* fs = &rq_idl[i][d].back();
											double fs_ratio = fs->firings / fs->span;
											double ratio_diff = firings[i][ri] / tau - fs_ratio;
											if (abs(ratio_diff) < beta * fs_ratio) { // if the last list element and this one can be merged then merge them
												fs->firings += firings[i][ri];
												fs->span += tau;
												add = false;
											}
										}
										if (add) { // if the last list element couldn't be merged with this one then add it to the reaction's list
											try {
												rq[i][d].push_back(T[chunk_index] + delay_times[d]);
												rq_idl[i][d].push_back(rq_node(firings[i][ri], tau));
											} catch (bad_alloc) { // if there isn't enough memory to allocate a new list element then exit the program
												cout << terminal_no_memory << endl;
												exit(1);
											}
										}
									}
								}
							}
							
							// update the appropriate concentrations
							for (unsigned int i = 0; i < cells; i++) {
								for (int j = 0; j < species; j++) {
									for (int k = 0; k < reactions; k++) {
										bool not_delayed = true;
										for (int d = 0; d < num_of_delayed_reactions; d++) { // check if the kth reaction is delayed
											if (k == delayed_reactions[d]) {
												not_delayed = false;
												break;
											}
										}
										if (not_delayed) { // if the reaction isn't delayed (and therefore hasn't updated its appropriate concentrations) update the concentrations for it
											cx[i * species + j] += firings[i][k] * species_update_values[j][k];
										}
									}
								}
								
								// update every propensity value since any could have changed
								double* acell = a[i];
								int xcell = i * species;
								update_a0_27(acell, &pars, cx, xcell, &a0);
								update_a1_2_4_6(acell, &pars, cx, xcell, &a0);
								update_a3_18(acell, &pars, cx, xcell, &a0);
								update_a4_9_10_12(acell, &pars, cx, xcell, &a0);
								update_a5_19(acell, &pars, cx, xcell, &a0);
								update_a6_12_15_16(acell, &pars, cx, xcell, &a0);
								update_a7_20(acell, &pars, cx, xcell, &a0);
								update_a8_29(acell, &pars, cx, xcell, &a0);
								update_a11_21(acell, &pars, cx, xcell, &a0);
								update_a13_22(acell, &pars, cx, xcell, &a0);
								update_a14_31(acell, &pars, cx, xcell, &a0);
								update_a17_23(acell, &pars, cx, xcell, &a0);
								update_a24_33(acell, &pars, cx, xcell, &a0);
								update_a25(acell, &pars, cx, xcell, &a0);
								update_a26_28_32(acell, &pars, cx, xcell, neighbors, nc[i], &a0);
							}
							
							// update the simulation timestep
							T[chunk_index] = nT;
						}
					}
				} while (repeat); // repeat if the concentrations would have become negative
			} else { // if tau-leaping has been temporarily disabled for the sake of efficiency
				if (skip_steps > 0) { 
					/* 
					   if tau-leaping should eventually be resumed 
					   then decrement the number of times the next-reaction-method should be run before doing so
					*/
					skip_steps--;
				}
				
				/*
				Next reaction method with delayed reactions:
				1) Update the propensity values related to the most recently fired reaction (skip this step the first iteration)
				2) Calculate delta
					a) Find the minimum of (Pk - Tk) / ak for each reaction to find delta
					b) Find the minimum delay time - current time for each delayed reaction to find delta
					c) Take delta as the minimum of steps a and b
				3) Update the simulation timestep
				4) Update the concentrations for the active reaction and 
				   either pop it off its list if it's a finished delayed reaction or add it to the list if it's a new one
				5) Update Pk and Tk with the appropriate values
				*/
				
				// update the propensity functions
				if (cell_index != -1) { // if there is a reaction to update the propensities with
					int* cnc = nc[cell_index]; // get the cell's neighbors
					int maxcell;
					if ((is_delayed && reaction_index == 24) || reaction_index == 25) { 
					    // reactions 24 and 25 require each of the cell's neighbors to update as well
						maxcell = neighbors;
					} else {
						maxcell = 1;
					}
					for (int cn = 0; cn < maxcell; cn++) { 
					    // for every cell that should update (usually just the active cell but see above about reactions 24 and 25)
						double* acell = a[cnc[cn]]; // the cell index for the propensities array
						int xcell = cnc[cn] * species; // the cell index for the concentrations array
						bool enter = true; // update the propensities for all non-delayed reactions and delayed reactions that are finishing, not starting
						for (int j = 0; j < num_of_delayed_reactions; j++) {
							if (reaction_index == delayed_reactions[j]) {
								enter = is_delayed;
								break;
							}
						}
						if (enter) {
							switch (reaction_index) { // update the propensities based on the reaction index (these functions are inlined and combined for efficiency)
								case 0:
								case 1:
									update_a1_2_4_6(acell, &pars, cx, xcell, &a0);
									break;
								case 2:
								case 3:
									update_a1_2_4_6(acell, &pars, cx, xcell, &a0);
									update_a3_18(acell, &pars, cx, xcell, &a0);
									update_a26_28_32(acell, &pars, cx, xcell, neighbors, nc[cnc[cn]], &a0);
									break;
								case 4:
								case 5:
									update_a1_2_4_6(acell, &pars, cx, xcell, &a0);
									update_a5_19(acell, &pars, cx, xcell, &a0);
									update_a4_9_10_12(acell, &pars, cx, xcell, &a0);
									break;
								case 6:
								case 7:
									update_a1_2_4_6(acell, &pars, cx, xcell, &a0);
									update_a7_20(acell, &pars, cx, xcell, &a0);
									update_a6_12_15_16(acell, &pars, cx, xcell, &a0);
									break;
								case 8:
								case 9:
									update_a4_9_10_12(acell, &pars, cx, xcell, &a0);
									break;
								case 10:
								case 11:
									update_a4_9_10_12(acell, &pars, cx, xcell, &a0);
									update_a11_21(acell, &pars, cx, xcell, &a0);
									break;
								case 12:
								case 13:
									update_a4_9_10_12(acell, &pars, cx, xcell, &a0);
									update_a13_22(acell, &pars, cx, xcell, &a0);
									update_a6_12_15_16(acell, &pars, cx, xcell, &a0);
									update_a26_28_32(acell, &pars, cx, xcell, neighbors, nc[cnc[cn]], &a0);
									break;
								case 14:
								case 15:
									update_a6_12_15_16(acell, &pars, cx, xcell, &a0);
									break;
								case 16:
								case 17:
									update_a6_12_15_16(acell, &pars, cx, xcell, &a0);
									update_a17_23(acell, &pars, cx, xcell, &a0);
									break;
								case 18:
									update_a3_18(acell, &pars, cx, xcell, &a0);
									update_a26_28_32(acell, &pars, cx, xcell, neighbors, nc[cnc[cn]], &a0);
									break;
								case 19:
									update_a5_19(acell, &pars, cx, xcell, &a0);
									break;
								case 20:
									update_a7_20(acell, &pars, cx, xcell, &a0);
									break;
								case 21:
									update_a11_21(acell, &pars, cx, xcell, &a0);
									break;
								case 22:
									update_a13_22(acell, &pars, cx, xcell, &a0);
									update_a26_28_32(acell, &pars, cx, xcell, neighbors, nc[cnc[cn]], &a0);
									break;
								case 23:
									update_a17_23(acell, &pars, cx, xcell, &a0);
									break;
								case 24:
								case 25:
									update_a25(acell, &pars, cx, xcell, &a0);
									update_a26_28_32(acell, &pars, cx, xcell, neighbors, nc[cnc[cn]], &a0);
									break;
								case 26:
								case 27:
									update_a0_27(acell, &pars, cx, xcell, &a0);
									break;
								case 28:
								case 29:
									update_a8_29(acell, &pars, cx, xcell, &a0);
									break;
								case 30:
								case 31:
									update_a14_31(acell, &pars, cx, xcell, &a0);
									break;
								case 32:
								case 33:
									update_a24_33(acell, &pars, cx, xcell, &a0);
									break;
							}
						}
					}
				}
				
				double delta = INFINITY; // the time change
				cell_index = -1; // the cell index of the active reaction
				reaction_index = -1; // the reaction index of the active reaction
				int dr_index = -1; // the delayed reaction index of the active reaction
				is_delayed = false; // whether or not the reaction is delayed
				
				// calculate the minimum delta using non-delayed reactions
				for (unsigned int i = 0; i < cells; i++) {
					for (int j = 0; j < reactions; j++) {
						if (a[i][j] != 0) {
							double delta_temp = (Pk[i][j] - Tk[i][j]) / a[i][j];
							if (delta_temp < delta) {
								delta = delta_temp;
								cell_index = i;
								reaction_index = j;
							}
						}
					}
				}
				
				// calculate the minimum delta using delayed reactions
				for (unsigned int i = 0; i < cells; i++) {
					for (int d = 0; d < num_of_delayed_reactions; d++) {
						if (!rq[i][d].empty()) {
							double delta_temp = rq[i][d].front() - T[chunk_index];
							if (delta_temp < delta) {
								delta = delta_temp;
								cell_index = i;
								reaction_index = delayed_reactions[d];
								dr_index = d;
								is_delayed = true;
							}
						}
					}
				}
				
				// update the simulation timestep
				T[chunk_index] += delta;
				
				// update the concentrations and appropriate delayed list if necessary
				int cell_offset = cell_index * species;
				if (is_delayed) { 
				    // if the current reaction is delayed then update the concentrations according to delayed update values
					cx[cell_offset + species_update_indices_delayed[dr_index]]++;
					
					// remove the current reaction from its delayed list
					if (!rq[cell_index][dr_index].empty()) {
						rq[cell_index][dr_index].pop_front();
						if (appx) { 
						    /* 
						       if tau-leaping is also activated (if -a or --algorithm is set to "appx") 
						       then update the firings and span list as well
							*/
							rq_idl[cell_index][dr_index].pop_front();
						}
					}
				} else { 
				    /* 
				       if the current reaction isn't a delayed one finishing 
				       then update the concentrations according to non-delayed update values
					   if the reaction is a delayed one starting then add it to the corresponding delayed reactions list
					*/
					bool do_update = true;
					for (int d = 0; d < num_of_delayed_reactions; d++) {
						if (reaction_index == delayed_reactions[d]) {
							try {
								rq[cell_index][d].push_back(T[chunk_index] + delay_times[d]);
								if (appx) {
									rq_idl[cell_index][d].push_back(rq_node(1, delta));
								}
							} catch (bad_alloc) { // if there isn't enough memory to allocate a new list index then exit the program
								cout << terminal_no_memory << endl;
								exit(1);
							}
							do_update = false; // don't update the concentrations if the reaction is a delayed one starting
							break;
						}
					}
					
					// if the reaction is not a delayed one starting then update the concentrations according to non-delayed update values
					if (do_update) {
						int end = species_update_indices[reaction_index][0];
						for (int i = 1; i <= end; i++) {
							int species_index = species_update_indices[reaction_index][i];
							cx[cell_offset + species_index] += species_update_values[species_index][reaction_index];
						}
					}
					
					// update the appropriate Pk with a new random value
					Pk[cell_index][reaction_index] += pk_dist(); // log(1 / unif_dist())
				}
				
				// update every Tk value according to delta and the propensity values
				for (unsigned int i = 0; i < cells; i++) {
					for (int j = 0; j < reactions; j++) {
						Tk[i][j] += a[i][j] * delta;
					}
				}
			}
			
			/*
			Print the results:
			1) If the difference in simulation time exceeds the granularity 
			   (0.1 minutes by default, but can be changed with -g or --granularity) then do the following:
			2) If the sum of the concentrations is negative then something went wrong and the simulation ends prematurely
			3) Increment the chunk index so the current one can be stored without being overwritten
			4) If the the final chunk index has been reached or the difference in simulation time warrants printing 
			   (60 minutes by default, but can be changed with -p or --print):
				a) Print every chunk element from 1 to chunk_index
				b) Reset the chunk index to 1
				c) Wrap the most recent concentration levels and simulation timestep back to index 1 
				   so they can be considered the previous ones
			*/
			
			// if the difference in timesteps exceeds the granularity then move to the next index in the chunk
			if (T[chunk_index] - T[chunk_index - 1] >= granularity) {
				// if the sum of the concentration values is less than 0 then end the simulation
				int sum = 0;
				for (unsigned int i = 0; i < cells; i++) {
					for (int j = 0; j < species; j++) {
						sum += x[chunk_index][i * species + j];
					}
				}
				if (sum < 0) {
					break;
				}
				
				// move to the next chunk index
				int prev_ci = chunk_index;
				chunk_index++;
				
				// at the end of each chunk
				if (chunk_index == chunk || (T[chunk_index - 1] - last_print >= print_interval)) {
					// print the current chunk's results
					store_results(&ofile[r], x, T, cells, chunk_index, r, con_level);
					
					// reset the chunk
					last_print = T[prev_ci];
					T[0] = T[chunk_index - 1];
					chunk_index = 1;
				}
				
				// wrap around the concentration levels and simulation timestep so they can be considered the previous ones
				for (unsigned int i = 0; i < cells; i++) {
					for (int j = 0; j < species; j++) {
						x[chunk_index][i * species + j] = x[prev_ci][i * species + j];
					}
				}
				cx = x[chunk_index];
				T[chunk_index] = T[prev_ci];
			}
		}
		
		// free any delayed reactions from their list
		for (unsigned int i = 0; i < cells; i++) {
			for (int d = 0; d < num_of_delayed_reactions; d++) {
				while (!rq[i][d].empty()) {
					rq[i][d].pop_front();
					if (appx) {
						rq_idl[i][d].pop_front();
					}
				}
			}
		}

		// close the output file and indicate the end of the run
		store_results(&ofile[r], x, T, cells, chunk_index + 1, r, con_level);
		ofile[r].close();
		cout << terminal_done << endl;
		cout.flush();
	}
	
	memory_dealloc(x, T, input_file, output_path, chunk);	
	return 0;
}



