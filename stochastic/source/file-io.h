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

#ifndef FILE_IO_H
#define FILE_IO_H

#include <iostream>

using namespace std;

#include "main.h"

void store_results(ofstream*, int**, double*, int, int, int, int);
void store_filename(char**, const char*);
void read_file(char*, char**);
void parse_line(char*, double[], int*);

#endif

