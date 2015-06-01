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

/*
 The program simulates the behavior of the zebrafish segmentation clock for two cells,
 through a matematical model using delay differential equations.
 To solve the equations, the program uses Euler's method, with an adjustable stepsize.
 The program tries to find parameter sets which replicate the behavior of the system in wild type
 and several mutants.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>

#include "functions.h"
#include "io_function.h"
using namespace std;

const int CHUNK_SIZE = 10000; // size of chunk which is read from the input buffer at a time
// global variables set in functions.cpp
extern char* terminal_blue;
extern char* terminal_red;
extern char* terminal_reset;

int main(int argc, char** argv)
{
    /*
     Process user input and handle output files.
     */
    // Check for user inputed arguments and initialize default values
    char *input_file = NULL, *output_path = NULL, *gradients_file = NULL, *ofeat_file = NULL; // input and output file paths (either absolute or relative)
    int PARS = 1; // number of parameters to simulate, default value is 1
    int x = 2, y = 1; // width and height of the tissue being simulated, default is 2x1
    int minutes = 1200; // number of minutes to run each simulation for, default is 1200
    int seed = time(0), pid = getpid(); // seed to be used for generating random numbers, default is time
	seed = abs(((seed*181)*((pid-83)*359))%805306457);    
    double eps = 0.01; // time step to be used for Euler's method, default is 0.01
    double max_prop = INFINITY; // maximum threshold for propensity functions, default is INFINITY
    bool toPrint = false, ofeat = false; // boolean marking whether or not concentrations should be printed to a text file
    checkArgs(argc, argv, &input_file, &output_path, &gradients_file, &ofeat_file, ofeat, PARS, seed, minutes, eps, max_prop, toPrint, x, y);
    rates rateValues[CHUNK_SIZE];

    // Read the entire input file into a char array buffer to speed up I/O
    char *buffer;
    char *gradients = NULL;
    int index = 0;
    create_buffer(buffer, input_file);
    create_buffer(gradients, gradients_file);  


    //Create output files
    ofstream allpassed, oft;
    create_output(output_path, toPrint, ofeat, ofeat_file, &allpassed, &oft);


    // Iterate through every paramater set
    for (int p = 0; p < PARS; p += CHUNK_SIZE) {
        int STEP = (PARS - p > CHUNK_SIZE ? CHUNK_SIZE : PARS - p);
        srand(seed);
        store_values(input_file, buffer, index, rateValues, STEP, seed); // Read the parameter sets from the buffer
        glevels gene(minutes / eps, x * y); // Create the structure that contains the 2D arrays which hold the gene levels
		for (int i = 0; i < STEP; i++) {
            string res;
	    	ostringstream convert;
	    	convert << i;
	    	res = convert.str();
            cerr << "Simulating set " << p + i << endl; // Used for creating output file names specific to the paramater
	    	int t_steps = int(minutes / eps); // Set the amount of time steps to be used in the simulation
		 
	    	rates temp_rate = rateValues[i]; // Temporary rate structure used to alter the protein synthesis rates in order to create mutants
	    	data of_wt, of_her1, of_her7, of_her13, of_her713, of_delta; // Data structures for storing oscillation features
	    	bool wt = false, her1 = false, her7 = false, her13 = false, her713 = false, delta = false;  // Booleans to see if we met mutant conditions on each iteration
		 
	    	// Clear data from previous iterations
	    	clear_data(of_wt);
	    	clear_data(of_her1);
	    	clear_data(of_her7);
	    	clear_data(of_her13);
	    	clear_data(of_her713);
	    	clear_data(of_delta);
			
            /* 
             For the wild type and every mutant, perform the following steps:
             1) Adjust the appropriate protein synthesis rates to create mutants if necessary
             2) Run the simulation
             3) Go to the next parameter set if the propensities for the current one have gone above the set threshold
             4) Otherwise, test oscillation features
             5) Go to the next parameter set if the current one did not produce oscillations or did not satisfy the mutant conditions
             */
            wt = run_mutant(&gene, t_steps, eps, temp_rate, of_wt, true, max_prop, x, y);
            if (toPrint) {
                printForPlotting(mutants[0] + "/run0.txt", &gene, t_steps, eps);
            }
	    	if (!wt) continue; 
            wt = fwildtype(of_wt.peaktotrough1, of_wt.peaktotrough2); 
            if (!wt) continue; 

            temp_rate.psd = 0.0; 
            delta = run_mutant(&gene, t_steps, eps, temp_rate, of_delta, false, max_prop, x, y);
            if (toPrint) {
                printForPlotting(mutants[1] + "/run0.txt", &gene, t_steps, eps);
            }
	    	if (!delta) continue; 
            delta = fd_mutant(of_delta.period, of_delta.amplitude, of_wt.period); 
            if (!delta) continue; 
            temp_rate.psd = rateValues[i].psd; 
            
            temp_rate.psh13 = 0.0; 
            her13 = run_mutant(&gene, t_steps, eps, temp_rate, of_her13, false, max_prop, x, y);
            if (toPrint) {
                printForPlotting(mutants[2] + "/run0.txt", &gene, t_steps, eps);
            }
	    	if (!her13) continue; 
            her13 = f13_mutant(of_her13.period, of_her13.amplitude, of_wt.period);
            if (!her13) continue; 
            temp_rate.psh13 = rateValues[i].psh13; 
            
            temp_rate.psh1 = 0.0;
            her1 = run_mutant(&gene, t_steps, eps, temp_rate, of_her1, false, max_prop, x, y);
            if (toPrint) {
                printForPlotting(mutants[3] + "/run0.txt", &gene, t_steps, eps);
            }
	    	if (!her1) continue;
            her1 = f1_mutant(of_her1.period, of_her1.amplitude, of_wt.period);
            if (!her1) continue;
            temp_rate.psh1 = rateValues[i].psh1;
            
            temp_rate.psh7 = 0.0;
            her7 = run_mutant(&gene, t_steps, eps, temp_rate, of_her7, true, max_prop, x, y);
            if (toPrint) {
                printForPlotting(mutants[4] + "/run0.txt", &gene, t_steps, eps);
            }
            if (!her7) continue;
            her7 = f7_mutant(of_her7.period, of_her7.amplitude, of_wt.period);
            if (!her7) continue;
            temp_rate.psh7 = rateValues[i].psh7;
            
            temp_rate.psh7 = 0.0, temp_rate.psh13 = 0.0;
            her713 = run_mutant(&gene, t_steps, eps, temp_rate, of_her713, true, max_prop, x, y);
            if (toPrint) {
                printForPlotting(mutants[5] + "/run0.txt", &gene, t_steps, eps);
            }
	    	if (!her713) continue;
            her713 = f713_mutant(of_her713.period, of_her713.amplitude, of_wt.period);
            if (!her713) continue;
            temp_rate.psh7 = rateValues[i].psh7;
            temp_rate.psh13 = rateValues[i].psh13;

            /*
             If the paramater set created oscillatory behavior in wild type and all the mutant conditions were satisfied:
             1) Print the appropriate message
             2) Print the oscillation features into the appropriate files
             2) Print the parameter set into the output file containing parameter that passed the conditions.
            */
            
            cerr << terminal_blue << "Parameter set " << i << " passed." << terminal_reset << endl;
            if (ofeat) {
                oft << res << "," << of_wt.period << "," << of_wt.amplitude << "," << of_wt.peaktotrough1 << ",";
                oft << of_delta.period << "," << of_delta.amplitude << "," << of_delta.peaktotrough1 << ",";
                oft << of_her1.period << "," << of_her1.amplitude << "," << of_her1.peaktotrough1 << ",";
                oft << of_her7.period << "," << of_her7.amplitude << "," << of_her7.peaktotrough1 << ",";
                oft << of_her13.period << "," << of_her13.amplitude << "," << of_her13.peaktotrough1 << ",";
                oft << of_her713.period << "," << of_her713.amplitude << "," << of_her713.peaktotrough1 << ",";
            }
            
            allpassed<<temp_rate.psh1<<","<<temp_rate.psh7<<","<<temp_rate.psh13<<","<<temp_rate.psd<<","<<temp_rate.pdh1<<","<<temp_rate.pdh7<<",";
	    	allpassed<<temp_rate.pdh13<<","<<temp_rate.pdd<<","<<temp_rate.msh1<<","<<temp_rate.msh7<<","<<temp_rate.msh13<<","<<temp_rate.msd<<",";
	    	allpassed<<temp_rate.mdh1<<","<<temp_rate.mdh7<<","<<temp_rate.mdh13<<","<<temp_rate.mdd<<"," <<temp_rate.ddgh1h1 << "," << temp_rate.ddgh1h7 << ",";
	    	allpassed<<temp_rate.ddgh1h13<<","<<temp_rate.ddgh7h7<<","<<temp_rate.ddgh7h13<<","<<temp_rate.ddgh13h13<<",";
	    	allpassed<<temp_rate.delaymh1<<","<<temp_rate.delaymh7<<","<<temp_rate.delaymh13<<","<<temp_rate.delaymd<<","<<temp_rate.delayph1<<",";
	    	allpassed<<temp_rate.delayph7<<","<<temp_rate.delayph13<<","<<temp_rate.delaypd<<","<<temp_rate.dah1h1<<","<<temp_rate.ddh1h1<<","<<temp_rate.dah1h7<<",";
	    	allpassed<<temp_rate.ddh1h7<<","<<temp_rate.dah1h13<<","<<temp_rate.ddh1h13<<","<<temp_rate.dah7h7<<","<<temp_rate.ddh7h7<<","<<temp_rate.dah7h13<<",";
	    	allpassed<<temp_rate.ddh7h13<<","<<temp_rate.dah13h13<<","<<temp_rate.ddh13h13<<","<<temp_rate.critph1h1 << "," << temp_rate.critph7h13 <<","<<temp_rate.critpd<<endl;
	    }
        
        cerr << terminal_blue << "Done with " << terminal_reset << STEP << " parameter sets." << endl;
	}

	if (input_file != NULL) {
        free(buffer);
	}
    allpassed.close();
    oft.close();
    return 0;
}
