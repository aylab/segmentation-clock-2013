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

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <math.h>

#include "functions.h"
#include "io_functions.h"


using namespace std;

void readFile(char **buffer, char* input_file)
{
    /*
     Reads the specified input file into a char* buffer.
     */    
    
    FILE *pFile;
    long lSize;
    size_t result;
    static const int MB = 1048576;

    // open file for reading and print header	
    pFile = fopen(input_file, "rb");
    if (pFile == NULL) {
	cerr << terminal_red << "Couldn't open file " << input_file << "! Exit status 1." << terminal_reset << endl;
        exit(1);
    }
	
    // obtain file size and check that it is under 400MB
    fseek(pFile, 0, SEEK_END);
    lSize = ftell(pFile);
    if (lSize > 400 * MB) {
        cerr << terminal_red << "Input file is too large (400MB limit)! Exit status 1." << terminal_reset;
        exit(1);
    }    
    rewind(pFile);
	
    //allocate memory to contain the whole file:
    *buffer = (char*) malloc (sizeof(char)*lSize);
    if (*buffer == NULL) {
	cerr << terminal_no_memory << " Exit status 2." << endl;
        exit(2);
    }
	
    //copy the file into the buffer:
    result = fread(*buffer, 1, lSize, pFile);
    if ((int)result != lSize) {
	cerr << terminal_red << "Reading error. Exit status 3." << terminal_reset << endl;
        exit(3);
    }
	
    //the whole file is now loaded in the memory buffer
    fclose (pFile);
}


void create_buffer (char *buffer, char *input_file){
    // Read the entire input file into a char array buffer to speed up I/O
    if (input_file != NULL) {
        cout << terminal_blue << "Reading file " << terminal_reset << input_file << " . . . ";
        readFile(&buffer, input_file);
        cout << terminal_done << endl;
    }
}

void terminal_color(){
    /* 
     Allocate memory for the terminal color code strings.
     */
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

void checkArgs(int argc, char** argv, char** input_file, char** output_path, char** ofeat_file, bool& ofeat, int& pars, int& seed, int& minutes, double& eps, double& max_prop, bool& toPrint, int &x, int &y) {
    terminal_color();    
    
    /*
     Parses the arguments given by the user and sets default options for the arguments which weren't specified.
     */
    if (argc > 1) { // if arguments were given and each argument option is followed by a value
        for (int i = 1; i < argc; i += 2) { // iterate through each argument pair
            char* option = argv[i];
            char* value;
            if (i < argc - 1) {
                value = argv[i + 1];
            } else {
                value = NULL;
            }
            
            /*
             Check for each possible argument option and overrides the default value for each specified option. If the option isn't recognized or the value given for an option doesn't appear valid then the usage information for the program is printed with an error message and no simulations are run. The code should be fairly self-explanatory with a few exceptions:
             1) atoi converts a string to an integer, atof converts a string to a floating point number (i.e. rational)
             2) strings should always be compared using strcmp, not ==, and strcmp returns 0 if the two strings match
             3) usage(true) prints the usage information with an error message while usage(false) prints it without one
             */
            
            if (strcmp(option, "-i") == 0 || strcmp(option, "--input") == 0) {
                store_filename(input_file, value);
            } else if (strcmp(option, "-o") == 0 || strcmp(option, "--output") == 0) {
                store_filename(output_path, value);
            } else if (strcmp(option, "-s") == 0 || strcmp(option, "--seed") == 0) {
                seed = atof(value);
                if (seed <= 0) {
                    usage("The seed to generate random numbers must be a positive integer. Set -s or --seed to at least 1.");
                }
            } else if (strcmp(option, "-p") == 0 || strcmp(option, "--parameters") == 0) {
                pars = atoi(value);
                if (pars < 1) {
                    usage("The number of parameteres to generate must be a positive integer. Set -p or --parameters to at least 1.");
                }
            } else if (strcmp(option, "-m") == 0 || strcmp(option, "--minutes") == 0) {
                minutes = atoi(value);
                if (minutes < 1) {
                    usage("The simulation must be run for at least one minute. Set -m or --minutes to at least 1.");
                }
            } else if (strcmp(option, "-e") == 0 || strcmp(option, "--epsilon") == 0) {
                eps = atof(value);
                if (eps <= 0) {
                    usage("The time step for Euler's method must be a positive real number. Set -e or --epsilon to be greater than 0.");
                }
            } else if (strcmp(option, "-f") == 0 || strcmp(option, "--ofeatures") == 0) {
                ofeat = true;
                store_filename(ofeat_file, value);
            } else if (strcmp(option, "-a") == 0 || strcmp(option, "--propensities") == 0) {
                max_prop = atof(value);
                if (max_prop == 0) {
                    usage("The propensities threshold must be a positive real number. Set -a or --propensities to be greater than 0.");
                }
            } else if (strcmp(option, "-x") == 0 || strcmp(option, "--width") == 0) {
                x = atoi(value);
                if (x < 2) {
                    usage("The tissue width must be at least two cells. Set -x or --width to at least 2.");
                }
            } else if (strcmp(option, "-y") == 0 || strcmp(option, "--height") == 0) {
                y = atoi(value);
                if (y < 1) {
                    usage("The tissue heigiht must be at least one cell. Set -y or --height to at least 1.");
                }
            } else if (strcmp(option, "-w") == 0 || strcmp(option, "--write") == 0) {
                toPrint = true;
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
            }
        }
        if (*output_path == NULL) {
            store_filename(output_path, "output");
        }
        if ((y == 1 && x < 2) || (y == 2 || y == 3) || (y > 3 && (x < 4 || y % 2 == 1 || x % 2 == 1))) {
            usage("Invalid simulation size. For two cell systems, x=2, y=1. For chains, x>=3, y=1. For tissues, x>=4 and even, y>=4 and even.");
        }
    }
}

