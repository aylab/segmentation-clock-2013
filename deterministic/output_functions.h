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

void create_directory(char *output_path);
void create_file(char output_file);
void create_mutant_dir(char output_file);
void create_ofeat(char *str);
void create_output(char *output_path, bool toPrint, bool ofeat, char *ofeat_name, ofstream *allpassed, ofstream *oft);