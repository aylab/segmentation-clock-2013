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

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
// escape codes to color the terminal output and shortcuts for common outputs (set -c or --no-color to disable these)
#ifndef INPUT_H
#define INPUT_H

#define terminal_blue_d "\x1b[34m"
#define terminal_red_d "\x1b[31m"
#define terminal_reset_d "\x1b[0m"
#define terminal_done terminal_blue << "Done" << terminal_reset
#define terminal_no_memory terminal_red << "Not enough memory!" << terminal_reset

void readFile(char **buffer, char* input_file);
void create_buffer (char *buffer, char *input_file);
void terminal_color();
void store_filename (char** field, const char* value);
void checkArgs(int, char**, char**, char**, char**, char**, bool&, int&, int&, int&, double&, double&, bool&, int&, int&);

#endif