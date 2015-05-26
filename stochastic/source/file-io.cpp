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

#include "file-io.h"

#include <fstream>
#include <stdlib.h>
#include <string.h>

#include "macros.h"
#include "main.h"

// global variables set in main.cpp
extern char* terminal_blue;
extern char* terminal_red;
extern char* terminal_reset;

// store the her1 values of the simulation up to the chunk index specified by chunk
void store_results (ofstream* file, int** x, double* T, int cells, int chunk, int run, int con_level) {
	try {
		for (int sn = 1; sn < chunk; sn++) {
			*file << T[sn] << "\t";
			for (int i = 0; i < cells; i++) {
				*file << x[sn][i * species + con_level] << "\t";
			}
			*file << endl;
		}
	} catch (ofstream::failure) { // if there was a problem writing the results then exit the program
		cout << terminal_red << "Couldn't write to run #" << run << "'s output file!" << terminal_reset << endl;
		exit(1);
	}
}

// store the filename for the input or the file path for the output
void store_filename (char** field, const char* value) {
	*field = (char*)malloc(strlen(value) + 1);
	if (*field == 0) { // if the name couldn't be allocated then exit the program
		cout << terminal_no_memory << endl;
		exit(1);
	}
	strcpy(*field, value);
}

// read the given file into a buffer
void read_file (char* filename, char** buffer) {
	FILE* file;
	long size;
	size_t result;
	
	file = fopen(filename, "rb");
	if (file == NULL) {
		cout << terminal_red << "Couldn't open " << filename << "!" << terminal_reset << endl;
		exit(1);
	}
	
	fseek(file, 0, SEEK_END);
	size = ftell(file);
	rewind(file);
	
	*buffer = (char*)malloc(sizeof(char)*size);
	if (*buffer == NULL) {
		cout << terminal_no_memory << endl;
		exit(1);
	}
	
	result = fread(*buffer, 1, size, file);
	if ((int)result != size) {
		cout << terminal_red << "Couldn't read " << file << "!" << terminal_reset << endl;
	}
	
	fclose(file);
}

// parse a line in a string buffer and update the current buffer index
void parse_line(char* buffer, double items[], int* index) {
	char db[50]; // maximum size of a number in printed bytes
	db[0] = '\0';
	int itindex = 0;
	int i;
	for (i = *index; buffer[i] != '\n' && buffer[i] != '\0'; i++) {
		if (buffer[i] == ',') {
			items[itindex++] = atof(db);
			db[0] = '\0';
		} else {
			int orsize = strlen(db);
			db[orsize] = buffer[i];
			db[orsize + 1] = '\0';
		}
	}
	*index = i + 1;
	items[itindex] = atof(db);
}

