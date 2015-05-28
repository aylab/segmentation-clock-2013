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

using namespace std;

void terminal_color();
void checkArgs(int argc, char **argv, int& xcells, int& ycells, int& max_minutes, long& max_timesteps, int& runs, int& seed, char **input_file, char **output_path, int& con_level, map<string, int> levels, double& granularity, int& print_interval, char **seed_file, bool& appx);
void checkSize(int xcells, int ycells, char *input_file, char *output_path, int& structure, int& neighbors);
void cells_neighbors(int structure, int neighbors, int cells, int xcells, int ***nc);
void memory_alloc(int chunk, int cells, int ***x, double **T);
void memory_dealloc(int **x, double *T, char *input_file, char *output_path, int chunk);
void usage (const char* message);
void licensing ();

/* interpolate linearly interpolates the value at the given location between two given points
	parameters:
		x: the location at which to interpolate the value
		x0: the left location
		x1: the right location
		y0: the left location's value
		y1: the right location's value
	returns: the interpolated value
	notes:
		x0 and x1 are integers because this function is used to interpolate in a cell tissue.
	todo:
*/
inline double interpolate (double x, int x0, int x1, double y0, double y1) {
	return y0 + (y1 - y0) * ((x - x0) / (x1 - x0));
}