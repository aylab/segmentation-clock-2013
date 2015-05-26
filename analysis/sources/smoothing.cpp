/*
Oscillation smoother for stochastic runs
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

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <string.h>

// escape codes to color the terminal output and shortcuts for common outputs (set -c or --color to no to disable these)
#define terminal_blue_d "\x1b[34m"
#define terminal_red_d "\x1b[31m"
#define terminal_reset_d "\x1b[0m"
#define terminal_done terminal_blue << "Done" << terminal_reset
#define terminal_no_memory terminal_red << "Not enough memory!" << terminal_reset

using namespace std;

char* terminal_blue;
char* terminal_red;
char* terminal_reset;

void usage (const char*);

int getfsize(char* fname) {
	int count = 0;
	ifstream ifile;
	string line;	
	ifile.open(fname, fstream::in);
	while (!ifile.eof()) {
		getline(ifile, line);
		count++;
	}
	ifile.close();
	return count;
}

int main(int argc, char** argv) {
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
	
    if (argc < 4) {
		usage("Smoothing requires input and output file names and the size of the interval to average on.");
	}
	
	const int FSIZE = getfsize(argv[1]);	
    ifstream ifile;
    ifile.open(argv[1], fstream::in);
    ofstream ofile;
    ofile.open(argv[2], fstream::out);
	const int MSIZE = atoi(argv[3]);
	int r, c;
	ifile >> r >> c;
	const int CELLS = r * c;
	double cell[CELLS][FSIZE];
	double tstep[FSIZE];
	double sum[CELLS];
	memset(sum, 0, sizeof(sum));
	ofile << r << " " << c << endl;
	
    //calculating the first value requires the first MSIZE / 2 values 
    for (int i = 0; i <= MSIZE / 2; i++) {
        ifile >> tstep[i];
		for (int n = 0; n < CELLS; n++) {
			ifile >> cell[n][i];
			sum[n] += cell[n][i];
		}
    }
	ofile << tstep[0] << "\t";
	for (int n = 0; n < CELLS; n++) {
		ofile << sum[n] / (MSIZE / 2) << "\t";
    }    
    ofile << endl;
    
	int cellindex = 1;
    //the next MSIZE / 2 values require incrementally more values
    for (int i = (MSIZE / 2) + 1; i < MSIZE; i++) {
        ifile >> tstep[i];
		for (int n = 0; n < CELLS; n++) {
			ifile >> cell[n][i];
			sum[n] += cell[n][i];
		}
        ofile << tstep[cellindex] << "\t";
		for (int n = 0; n < CELLS; n++) {
			ofile << sum[n] / (i + 1) << "\t";
		}
        ofile << endl;
        cellindex++;
    }
	
    // the next values require all MSIZE values 
    int addindex = MSIZE, removeindex = 0;
    while (!ifile.eof() && addindex < FSIZE) {
        ifile >> tstep[addindex];
		for (int n = 0; n < CELLS; n++) {
			ifile >> cell[n][addindex];
			sum[n] += cell[n][addindex];
		}
        addindex++;
		
        for (int n = 0; n < CELLS; n++) {
			sum[n] -= cell[n][removeindex];
		}
        removeindex++;
		
        ofile << tstep[cellindex] << "\t";
		for (int n = 0; n < CELLS; n ++) {
			ofile << sum[n] / MSIZE << "\t";
		}
		ofile << endl;
        cellindex++;
    }
    
    ofile.close();
    
    return 0;
}

// usage message shown when an argument option or option value is invalid
void usage (const char* message) {
	if (strcmp(message, "") != 0) { // if there is an error message to print then print it
		cout << terminal_red << message << terminal_reset << endl << endl;
	}
	cout << "Usage: [-c|--no-color] <input file> <output file> <size of averaging interval>" << endl;
	exit(0);
}

