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

#include <errno.h>
#include <fstream>
#include <list>
#include <map>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <sys/stat.h>

#include "macros.h"
#include "utility.h"

using namespace std;

// global variables used in main.cpp and file-io.cpp
extern char* terminal_blue;
extern char* terminal_red;
extern char* terminal_reset;

void terminal_color(){
	// allocate memory for the terminal color code strings
	terminal_blue = (char*)malloc(sizeof(terminal_blue_d));
	terminal_red = (char*)malloc(sizeof(terminal_red_d));
	terminal_reset = (char*)malloc(sizeof(terminal_reset_d));
	if (terminal_blue == NULL || terminal_red == NULL || terminal_reset == NULL) {
		cout << terminal_no_memory << endl;
		exit(1);
	}
	strcpy(terminal_blue, terminal_blue_d);
	strcpy(terminal_red, terminal_red_d);
	strcpy(terminal_reset, terminal_reset_d);
}

void checkArgs(int argc, char **argv, int& xcells, int& ycells, int& max_minutes, long& max_timesteps, int& runs, int& seed, char **input_file, char **output_path, int& con_level, map<string, int> levels, double& granularity, int& print_interval, char **seed_file, bool& appx){
	if (argc > 1) { // if arguments were given
		for (int i = 1; i < argc; i += 2) { // iterate through each argument
			char* option = argv[i];
			char* value;
			if (i < argc - 1) {
				value = argv[i + 1];
			} else {
				value = NULL;
			}
			
			/*
			Check for each possible argument option and overrides the default value for each specified option. 
			If the option isn't recognized or the value given for an option doesn't appear valid 
			then the usage information for the program is printed with an error message and no simulations are run. 
			The code should be fairly self-explanatory with a few exceptions:
			1) atoi converts a string to an integer, atof converts a string to a floating point number (i.e. rational)
			2) strings should always be compared using strcmp, not ==, and strcmp returns 0 if the two strings match
			3) usage(true) prints the usage information with an error message while usage(false) prints it without one
			*/
			
			if (strcmp(option, "-x") == 0 || strcmp(option, "--width") == 0) {
				xcells = atoi(value);
				if (xcells < 1) {
					usage("The tissue width must be at least two cells. Set -x or --width to at least 2.");
				}
			} else if (strcmp(option, "-y") == 0 || strcmp(option, "--height") == 0) {
				ycells = atoi(value);
				if (ycells < 1) {
					usage("The tissue height must be at least one cell. Set -y or --height to at least 1.");
				}
			} else if (strcmp(option, "-m") == 0 || strcmp(option, "--minutes") == 0) {
				max_minutes = atoi(value);
				if (max_minutes < 1) {
					usage("The simulation must run for at least one minute. Set -m or --minutes to at least 1.");
				}
			} else if (strcmp(option, "-t") == 0 || strcmp(option, "--time-steps") == 0) {
				max_timesteps = atol(value);
				if (max_timesteps < 1) {
					usage("The simulation must run for at least one timestep. Set -t or --timesteps to at least 1.");
				}
			} else if (strcmp(option, "-r") == 0 || strcmp(option, "--runs") == 0) {
				runs = atoi(value);
				if (runs < 1) {
					usage("The simulation must run at least once. Set -r or --runs to at least 1.");
				}
			} else if (strcmp(option, "-s") == 0 || strcmp(option, "--seed") == 0) {
				seed = atoi(value);
				if (seed < 1) {
					usage("The seed to generate random numbers must be a positive integer. Set -s or --seed to at least 1.");
				}
			} else if (strcmp(option, "-i") == 0 || strcmp(option, "--input") == 0) {
				store_filename(input_file, value);
			} else if (strcmp(option, "-o") == 0 || strcmp(option, "--output") == 0) {
				store_filename(output_path, value);
			} else if (strcmp(option, "-l") == 0 || strcmp(option, "--con-level") == 0) {
				con_level = levels[value];
				if (con_level == 0) {
					usage("The concentration level to print must be her1, her7, her13, or delta mRNA, Her1, Her7, Her13, or Delta protein, or a dimer between any two proteins (e.g. Her1Her13)");
				} else {
					con_level--;
				}
			} else if (strcmp(option, "-g") == 0 || strcmp(option, "--granularity") == 0) {
				granularity = atof(value);
				if (granularity < 0.0) {
					usage("The output granularity must be a nonnegative real. Set -g or --granularity to at least 0.0.");
				}
			} else if (strcmp(option, "-p") == 0 || strcmp(option, "--print") == 0) {
				print_interval = atoi(value);
				if (print_interval < 1) {
					usage("The simulation cannot print its output more than once a simulation-minute. Set -p or --print to at least 1.");
				}
			} else if (strcmp(option, "-k") == 0 || strcmp(option, "--keep-seed") == 0) {
				store_filename(seed_file, value);
			} else if (strcmp(option, "-a") == 0 || strcmp(option, "--approximate") == 0) {
				appx = true;
				i--;
			} else if (strcmp(option, "-c") == 0 || strcmp(option, "--no-color") == 0) {
				strcpy(terminal_blue, "");
				strcpy(terminal_red, "");
				strcpy(terminal_reset, "");
				i--;
			} else if (strcmp(option, "-q") == 0 || strcmp(option, "--quiet") == 0) {
				ofstream nullout("/dev/null");
				cout.rdbuf(nullout.rdbuf());
				i--;
			} else if (strcmp(option, "-h") == 0 || strcmp(option, "--help") == 0) {
				usage("");
				i--;
			} else if (strcmp(option, "-l") == 0 || strcmp(option, "--licensing") == 0) {
				licensing();
				i--;
			} else {
				const char* temp_message = " is not a valid option.";
				char message[strlen(option) + strlen(temp_message) + 1];
				strcat(message, option);
				strcat(message, temp_message);
				usage(message);
			}
		}
	}

}

void checkSize(int xcells, int ycells, char *input_file, char *output_path, int& structure, int& neighbors){
	/*
	Check if the tissue size is valid:
	1) If it's a horizontal chain then it must be a minimum of 3x1
	2) If it's a hexagonal tissue grid then it must be a minimum of 4x4 and always even by even, not odd by even or even by odd
	3) If the tissue size isn't valid then print the usage information and don't run any simulations
	*/
	
	if ((ycells == 1 && xcells < 2) || (ycells == 2 || ycells == 3) || (ycells > 3 && (xcells < 4 || ycells % 2 == 1 || xcells % 2 == 1))) {
		usage("Invalid tissue size. For two cell systems, x=2, y=1. For chains, x>=3, y=1. For tissues, x>=4 and even, y>=4 and even.");
	}
	
    if (input_file == NULL) { // if the input file wasn't given by command-line then initialize it to "input.txt"
        store_filename(&input_file, "input.txt");
    }
    if (output_path == NULL) { // if the output directory path wasn't given then initialize it to "output"
        store_filename(&output_path, "output");
    }	
	
	if (ycells == 1) {
		if (xcells == 2) {
			structure = structure_twocell;
			neighbors = neighbors_for_twocell; // if there are only two cells then they are neighbors of each other, so they have two neighbors including themselves
		} else {
			structure = structure_chain;
			neighbors = neighbors_for_chain; // if the cells are in a chain then they have 3 neighbors including themselves
		}
	} else {
		structure = structure_tissue;
		neighbors = neighbors_for_tissue; // if the cells are in a tissue then they have 7 neighbors including themselves
	}
}

void cells_neighbors(int structure, int neighbors, int cells, int xcells, int ***nc){
	/*
	Pre-compute the neighbors for each cell. In two-cell systems, both cells are neighbors of each other. Chains wrap horizontally and hexagonal tissue grids wrap horizontally and vertically like a honeycomb.
	
										 ___  ___
	Two-cell systems look like this:	/   \/   \				where 1 and 2 are neighbors of each other
										| 1 || 2 |
										\___/\___/
	
										 ___  ___  ___  ___
	Chains of cells look like this:		/   \/   \/   \/   \	where x has neighbors n
										| n || x || n ||   |
										\___/\___/\___/\___/
	
										 ___  ___  ___  ___
	Tissues of cells look like this:	/   \/   \/   \/   \	where x has neighbors n
										|   || n || n ||   |
										\___/\___/\___/\___/_
										  /   \/   \/   \/   \
										  | n || x || n ||   |
										  \___/\___/\___/\___/
										/   \/   \/   \/   \
										|   || n || n ||   |
										\___/\___/\___/\___/
										  /   \/   \/   \/   \
										  |   ||   ||   ||   |
										  \___/\___/\___/\___/
	*/
	
	
	if (structure == structure_twocell) { // for two-cell systems
		for (unsigned int i = 0; i < cells; i++) {
			nc[i][0] = i;
			nc[i][1] = 1 - i;
		}
	} else if (structure == structure_chain) { // for chains
		for (unsigned int i = 0; i < cells; i++) {
			nc[i][0] = i;
			unsigned int end = xcells - 1;
			if (i == 0) {
				nc[i][1] = end;
				nc[i][2] = 1;
			} else if (i == end) {
				nc[i][1] = end - 1;
				nc[i][2] = 0;
			} else {
				nc[i][1] = i - 1;
				nc[i][2] = i + 1;
			}
		}
	} else { // for tissues
		for (unsigned int i = 0; i < cells; i++) {
			nc[i][0] = i;
			unsigned int not_fc = xcells * (i == 0);
			unsigned int not_lc = xcells * (i == xcells - 1);
			if ((i / xcells) % 2 == 0) {
				nc[i][1] = i - xcells + not_fc - 1;
				nc[i][2] = i - xcells;
				nc[i][3] = i + not_fc - 1;
				nc[i][4] = i + not_lc + 1;
				nc[i][5] = i + xcells + not_fc - 1;
				nc[i][6] = i + xcells;
			} else {
				nc[i][1] = i - xcells;
				nc[i][2] = i + not_fc - 1;
				nc[i][3] = i + not_lc + 1;
				nc[i][4] = i + xcells + not_fc - 1;
				nc[i][5] = i + xcells;
				nc[i][6] = i + xcells + not_lc + 1;
			}
			if (i <= xcells || cells - i <= xcells) {
				for (int j = 1; j < 7; j++) {
					nc[i][j] = (nc[i][j] + cells) % cells;
				}
			}
		}
	}
}

void memory_alloc(int chunk, int cells, int ***x, double **T){
	/*
	Allocate memory:
	1) Indicate to the user how much memory is initially allocated 
	   (more will be allocated later when delayed reactions must be added to lists)
	2) Establish a chunk size proportional to the maximum number of simulation minutes divided by the output granularity 
	   (+ 1 for wrapping purposes)
		Each element in the chunk stores the concentrations for a discrete timestep that will eventually be printed. 
		By storing many concentration snapshots (specifically, the maximum number required for one run), 
		printing can be done either at the end or at various intervals (depending on the printing interval command-line argument). 
		This way the program is not slowed down by repeated printing every timestep.
	3) Allocate memory for the concentrations matrix and the simulation timestep array
	*/
	
	
	cout << terminal_blue << "Allocating memory " << terminal_reset << "(" << (chunk * (cells * species * sizeof(int) + sizeof(double))) / MB << " MB) ... ";
	cout.flush();	
	try {
		/* 
		   create an array of size chunk that will contain enough entries for at least a whole run 
		   (it may not be fully utilized by one run if timesteps jump quickly enough)
		*/
		*x = new int*[chunk]; 
		for (int i = 0; i < chunk; i++) {
			/* 
			   the cells x species concentrations matrix is represented as a flat array, 
			   where each cell gets an element for each species, followed by the next cell's species elements
			*/
			(*x)[i] = new int[cells * species];
		}

		*T = new double[chunk];
	} catch (bad_alloc) { // if there isn't enough memory to allocate the structures then exit the program
		cout << terminal_no_memory << endl;
		exit(1);
	}
	cout << terminal_done << endl;
}

void memory_dealloc(int **x, double *T, char *input_file, char *output_path, int chunk){
	/*
	Deallocate heap memory:
	1) Delete the concentrations matrix and simulation timestep array
	2) Free the memory used to store the input and output paths
	*/
	
	for (int i = 0; i < chunk; i++) {
		delete[] x[i];
	}
	delete[] x;
	delete[] T;
	
	free(input_file);
    free(output_path);
	free(terminal_blue);
	free(terminal_red);
	free(terminal_reset);
}

// usage message shown when "-h" or "--help" is given or when an an argument option or option value is invalid
void usage (const char* message) {
	bool has_message = strcmp(message, "") != 0;
	if (has_message) { // if there is an error message to print then print it
		cout << terminal_red << message << terminal_reset << endl << endl;
	}
	cout << "Usage: [-option [value]]... [--option [value]]..." << endl;
	cout << "-x, --width       : the tissue width (in cells), min=3 for chain, min=4 and even for tissue, default=2" << endl;
	cout << "-y, --height      : the tissue height (in cells), min=1 for chain, min=4 and even for tissue, default=1" << endl;
	cout << "-m, --minutes     : the maximum number of minutes to simulate before ending, min=1, default=1200" << endl;
	cout << "-t, --time-steps  : the maximum number of timesteps to simulate before ending, min=1, default=10^12" << endl;
	cout << "-r, --runs        : the number of runs, min=1, default=1" << endl;
	cout << "-s, --seed        : the seed to generate random numbers, min=1, default=time" << endl;
	cout << "-i, --input       : the input path and file to accept parameters from, default=input.txt" << endl;
	cout << "-o, --output      : the path to print the output (i.e. results) to, default=output" << endl;
	cout << "-l, --con-level   : the concentration level to print (her1, her7, her13, delta, Her1, Her7, Her13, Delta, or a dimer (e.g. Her1Her13)), default=her1" << endl;
	cout << "-g, --granularity : the granularity of the output, i.e. print values for every x number of minutes simulated, min=0, default=0.1" << endl;
	cout << "-p, --print       : printing interval in minutes, min=1, default=1200" << endl;
	cout << "-k, --keep-seed   : store the seed in the specified file relative to the output directory, default=seed.txt" << endl;
	cout << "-a, --approximate : approximate the simulation for faster results, default=unused" << endl;
	cout << "-c, --no-color    : disable coloring the terminal output, default=unused" << endl;
	cout << "-q, --quiet       : hide the terminal output, default=unused" << endl;
	cout << "-l, --licensing   : view licensing information (no simulations will be run)" << endl;
	cout << "-h, --help        : view usage information (i.e. this) (no simulations will be run)" << endl;
	cout << endl << terminal_blue << "Example: ./stochastic -x 10 -y 6 --runs 5 --minutes 1200 -a --input set.csv --output results" << terminal_reset << endl << endl;
	exit(has_message);
}

// licensing information
void licensing () {
	cout << "Stochastic simulator for zebrafish segmentation" << endl;
	cout << "Copyright (C) 2012 Ahmet Ay (aay@colgate.edu), Jack Holland (jholland@colgate.edu), Adriana Sperlea (asperlea@colgate.edu)" << endl;
    cout << "This program comes with ABSOLUTELY NO WARRANTY" << endl;
    cout << "This is free software, and you are welcome to redistribute it under certain conditions;" << endl;
	cout << "You can use this code and modify it as you wish under the condition that you refer to the article: ???" << endl;
	exit(0);
}

void fill_gradients (rates& rs, char* gradients) {
	if (gradients != NULL) {
		static const char* usage_message = "There was an error reading the given gradients file.";
		int con; // The index of the concentration
		int column; // The column in the cell tissue
		double factor; // The factor to apply
		int last_column = 0; // The last column in the cell tissue given a gradient factor
		int start_column; // The first column to apply the gradient from
		int i = 0; // The index in the buffer
		while (gradients[i] != '\0') {
			// Read the concentration value
			if (sscanf(gradients + i, "%d", &con) != 1) {
				usage(usage_message);
			}
			if (con < 0 || con > NUM_RATES) {
				usage("The given gradients file includes rate indices outside of the valid range. Please adjust the gradients file or add the appropriate rates by editing the macros file and recompiling.");
			}
			rs.using_gradients = true; // Mark that at least one concentration has a gradient
			rs.has_gradient[con] = true; // Mark that this concentration has a gradient
			
			// Read every (position factor) pair
			while (gradients[i++] != ' ') {} // Skip past the concentration index
			while (not_EOL(gradients[i])) {
				// Read a position factor pair
				if (sscanf(gradients + i, "(%d %lf)", &column, &factor) != 2) {
					usage(usage_message);
				}
				if (column < 0 || column >= rs.width) {
					usage("The given gradients file includes positions outside of the given simulation width. Please adjust the gradients file or increase the width of the simulation using -x or --total-width.");
				}
				if (factor < 0) {
					usage("The given gradients file includes factors less than 0. Please adjusted the gradients file.");
				}
				
				// Apply the gradient factor
				factor /= 100;
				start_column = last_column;
				last_column = column;
				for (int j = start_column + 1; j < column; j++) {
					rs.factors_gradient[con][j] = interpolate(j, start_column, column, rs.factors_gradient[con][start_column], factor);
				}
				rs.factors_gradient[con][column] = factor;
				while (gradients[i++] != ')') {} // Skip past the end of the pair
				while (gradients[i] == ' ') {i++;} // Skip any whitespace before the next pair
			}
			
			// Apply the last gradient factor to the rest of the columns
			for (int j = column + 1; j < rs.width; j++) {
				rs.factors_gradient[con][j] = rs.factors_gradient[con][column];
			}
			i++;
		}
	}
}