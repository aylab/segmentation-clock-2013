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

#include "t-test.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

char* terminal_blue;
char* terminal_green;
char* terminal_red;
char* terminal_reset;

int main (int argc, char** argv) {
	char* table = NULL;
	char* run1 = NULL;
	char* run2 = NULL;
	terminal_blue = (char*)malloc(sizeof(terminal_blue_d));
	terminal_green = (char*)malloc(sizeof(terminal_green_d));
	terminal_red = (char*)malloc(sizeof(terminal_red_d));
	terminal_reset = (char*)malloc(sizeof(terminal_reset_d));
	if (terminal_blue == NULL || terminal_red == NULL || terminal_reset == NULL) {
		cout << terminal_no_memory << endl;
		exit(1);
	}
	strcpy(terminal_blue, terminal_blue_d);
	strcpy(terminal_green, terminal_green_d);
	strcpy(terminal_red, terminal_red_d);
	strcpy(terminal_reset, terminal_reset_d);
	
	if (argc < 10) {
		if (argc == 2) {
			if (strcmp(argv[1], "-h") == 0) {
				usage();
			} else if (strcmp(argv[1], "-l") == 0) {
				licensing();
			} else {
				usage();
			}
		} else if (argc > 2) {
			if (argc % 2 == 1) {
				for (int a = 1; a < argc; a += 2) {
					char* option = argv[a];
					char* value = argv[a + 1];
					if (strcmp(option, "-t") == 0 || strcmp(option, "--table") == 0) {
						table = value;
					} else if (strcmp(argv[a], "-1") == 0 || strcmp(option, "--run1") == 0) {
						run1 = value;
					} else if (strcmp(option, "-2") == 0 || strcmp(option, "--run2") == 0) {
						run2 = value;
					} else if (strcmp(option, "-c") == 0 || strcmp(option, "--color") == 0) {
						if (strcmp(value, "no") == 0) {
							strcpy(terminal_blue, "");
							strcpy(terminal_red, "");
							strcpy(terminal_reset, "");
						} else if (strcmp(argv[a + 1], "yes") != 0) {
							usage();
						}
					} else {
						usage();
					}
				}
			} else {
				usage();
			}
		} else {
			usage();
		}
	} else {
		usage();
	}
	
	cout << terminal_blue << "Opening " << terminal_reset << table << ", " << run1 << ", and " << run2 << " ... ";
	cout.flush();
	int ch;
	
	FILE* table_file = open_file(table);
	int table_columns = get_columns(table_file) + 1;
	rewind(table_file);
	int table_rows = 0;
	do {
		ch = getc(table_file);
		if (ch == '\n') {
			table_rows++;
		}
	} while (ch != EOF);
	
	FILE* run1_file = open_file(run1);
	int run1_columns = get_columns(run1_file);
	
	FILE* run2_file = open_file(run2);
	int run2_columns = get_columns(run2_file);
	
	cout << terminal_done << endl;
	
	int table_bytes = table_rows * table_columns * sizeof(double);
	int run1_bytes = run1_columns * sizeof(double);
	int run2_bytes = run2_columns * sizeof(double);
	int total_bytes = table_bytes + run1_bytes + run2_bytes;
	cout << terminal_blue << "Allocating memory " << terminal_reset << "(" << (total_bytes) << "B) " << "... ";
	cout.flush();
	double* table_data = (double*)allocate(table_bytes);
	double* run1_data = (double*)allocate(run1_bytes);
	double* run2_data = (double*)allocate(run2_bytes);
	
	cout << terminal_done << endl;
	
	cout << terminal_blue << "Reading data " << terminal_reset << "... ";
	cout.flush();
	
	rewind(table_file);
	table_data[0] = 0;
	for (int c = 1; c < table_columns; c++) {
		table_data[c] = read_number(table_file);
	}
	for (int r = 1; r < table_rows; r++) {
		for (int c = 0; c < table_columns; c++) {
			table_data[r * table_columns + c] = read_number(table_file);
		}
	}
	
	rewind(run1_file);
	for (int c = 0; c < run1_columns; c++) {
		run1_data[c] = read_number(run1_file);
	}
	
	rewind(run2_file);
	for (int c = 0; c < run2_columns; c++) {
		run2_data[c] = read_number(run2_file);
	}
	
	cout << terminal_done << endl;
	
	cout << terminal_blue << "Calculating results " << terminal_reset << "... ";
	cout.flush();
	
	double mean1 = mean(run1_data, run1_columns);
	double mean2 = mean(run2_data, run2_columns);
	double sample_var1 = sample_var(run1_data, run1_columns, mean1);
	double sample_var2 = sample_var(run2_data, run2_columns, mean2);
	double t = (mean1 - mean2) / sqrt(sample_var1 / run1_columns + sample_var2 / run2_columns);
	int df = deg_free(sample_var1, run1_columns, sample_var2, run2_columns);
	if (df >= table_rows) {
		cout << terminal_red << "Your table data does not extend to " << df << " degrees of freedom!" << terminal_reset << endl;
		exit(1);
	}
	
	int pvalue_index = table_columns;
	for (int c = table_columns - 1; c > 0; c--) {
		if (t < table_data[df * table_columns + c]) {
			pvalue_index = c;
		}
	}
	
	double pvalue;
	if (pvalue_index == table_columns) {
		pvalue = table_data[table_columns - 1];
	} else if (pvalue_index == 1) {
		pvalue = table_data[1];
	} else {
		int index = df * table_columns + pvalue_index;
		pvalue = (table_data[pvalue_index] - table_data[pvalue_index - 1]) / (table_data[index] - table_data[index - 1]) * (t - table_data[index - 1]) + table_data[pvalue_index - 1];
	}
	
	cout << terminal_done << endl;
	cout << terminal_green << "t=" << t << ", df=" << df << ", p-value=" << pvalue << terminal_reset << endl;
	
	free(terminal_blue);
	free(terminal_green);
	free(terminal_red);
	free(terminal_reset);
	
	return 0;
}

inline void* allocate (int bytes) {
	void* block = malloc(bytes);
	if (block == NULL) {
		cout << terminal_no_memory << endl;
		exit(1);
	}
	return block;
}

inline FILE* open_file (char* filename) {
	FILE* f = fopen(filename, "r");
	if (f == NULL) {
		cout << terminal_red << "Couldn't open " << filename << "!" << endl;
		exit(1);
	}
	return f;
}

int get_columns (FILE* f) {
	int columns = 1;
	int ch;
	do {
		ch = getc(f);
		if (ch == ',') {
			columns++;
		}
	} while (ch != '\n' && ch != EOF);
	return columns;
}

double read_number (FILE* f) {
	char buffer[50];
	int ch;
	int i = 0;
	do {
		ch = getc(f);
		buffer[i++] = ch;
	} while (!(ch == ',' || ch == '\n' || ch == EOF));
	buffer[i] = '\0';
	return atof(buffer);
}

double mean (double* data, int size) {
	double sum = 0.0;
	for (int i = 0; i < size; i++) {
		sum += data[i];
	}
	return sum / size;
}

double sample_var (double* data, int size, double mean) {
	double sum = 0.0;
	for (int i = 0; i < size; i++) {
		double diff = mean - data[i];
		sum += diff * diff;
	}
	return sum / (size - 1);
}

int deg_free (double sample_var1, int deg_free1, double sample_var2, int deg_free2) {
	double a = sample_var1 / deg_free1;
	double b = sample_var2 / deg_free2;
	double c = a + b;
	c = c * c;
	double d = a * a / (deg_free1 - 1) + b * b / (deg_free2 - 1);
	return c / d;
}

// usage message shown when "-h" or "--help" is given as the first argument or when an an argument option or option value is invalid
void usage () {
	cout << "Usage: [-h] | [-l] | [ [-option1 option_value] [-option2 option_value] ... ]" << endl;
	cout << "-t, --table       : the filename for your critical t values table" << endl;
	cout << "-1, --run1        : the filename for your first run's dataset" << endl;
	cout << "-2, --run2        : the filename for your second run's dataset" << endl;
	cout << "-c, --color       : whether or not to color the terminal output (yes/no), default=yes" << endl;
	cout << "-l, --licensing   : view licensing information (no simulations will be run)" << endl;
	cout << "-h, --help        : view usage information (i.e. this)" << endl;
	cout << endl << terminal_blue << "Example: ./t-test --table t-test-table.csv -1 run1.csv -2 run2.csv" << terminal_reset << endl << endl;
	exit(0);
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

