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

#include <stdio.h>

void* allocate(int);
FILE* open_file(char*);
int get_columns(FILE*);
double read_number(FILE*);
double mean(double*, int);
double sample_var(double*, int, double);
int deg_free(double, int, double, int);
void usage();
void licensing();

// escape codes to color the terminal output and shortcuts for common outputs (set -c or --color to no to disable these)
#define terminal_blue_d "\x1b[34m"
#define terminal_green_d "\x1b[32m"
#define terminal_red_d "\x1b[31m"
#define terminal_reset_d "\x1b[0m"
#define terminal_done terminal_blue << "Done" << terminal_reset
#define terminal_no_memory terminal_red << "Not enough memory!" << terminal_reset

