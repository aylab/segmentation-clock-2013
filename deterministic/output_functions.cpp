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

/*
 The program simulates the behavior of the zebrafish segmentation clock for two cells,
 through a matematical model using delay differential equations.
 To solve the equations, the program uses Euler's method, with an adjustable stepsize.
 The program tries to find parameter sets which replicate the behavior of the system in wild type
 and several mutants.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>

#include "output_functions.h"
#include "functions.h" 
using namespace std;

// global variables
extern char* terminal_blue;
extern char* terminal_red;
extern char* terminal_reset;

void create_directory(char *output_path){
    cout << terminal_blue << "Creating " << terminal_reset << output_path << " directory if necessary ... ";
    if (mkdir(output_path, 0755) != 0 && errno != EEXIST) { 
        /* 
         *  make the output directory if necessary 
         *  or exit the program if there was an error trying to do so
        */
        cout << terminal_red << "Couldn't create " << output_path << " directory!" << terminal_reset << endl;
        exit(1);
    }
    cout << terminal_done << endl;
}

void create_file(char *output_file, ofstream *allpassed){  
    try {
        cout << terminal_blue << "Creating output file " << terminal_reset << output_file << " . . . ";        
        (*allpassed).open(output_file, fstream::out);
    } catch (ofstream::failure) {
        cout << terminal_red << "Couldn't create " << output_file << "!" << " Exit status 1." << terminal_reset << endl;
        exit(1);
    }
    cout << terminal_done << endl;
}

void create_mutant_dir(char *output_path, string mutants[]){    
    char output_dir[strlen(output_path)];
    for (int mut = 0; mut < 6; mut++) {
        strcpy(output_dir, output_path);
        strcat(output_dir, mutants[mut].c_str());
            
        string mutant;
        ostringstream convert;
        convert << output_dir;
        mutant = convert.str();
        mutants[mut] = mutant;
           
        cout << terminal_blue << "Creating " << terminal_reset << output_dir << " directory if necessary ... ";
        if (mkdir(output_dir, 0755) != 0 && errno != EEXIST) { // make the output directory if necessary or exit the program if there was an error trying to do so
            cout << terminal_red << "Couldn't create " << output_dir << " directory!" << terminal_reset << endl;
            exit(1);
        }
        cout << terminal_done << endl;
    }
}

void create_ofeat(char *str, char *ofeat_file, ofstream *oft){
    try {
        cout << terminal_blue << "Creating oscillation features file " << terminal_reset << ofeat_file << " . . . " << endl;
        (*oft).open(ofeat_file, fstream::out);
        (*oft)<< str << endl;
        cout << terminal_done << endl;
    } catch (ofstream::failure) {
        cout << terminal_red << "Couldn't create " << ofeat_file << "!" << " Exit status 1." << terminal_reset << endl;
        exit(1);
    }
}

void create_file_name(char *buff, char *output_path, int path_length, char *file){
    int file_name_length = strlen(file);
    memcpy(buff, output_path, path_length);
    buff[path_length] = '/';
    memcpy(buff + path_length + 1, file , file_name_length); 
    buff[path_length + file_name_length + 1] = '\0';
}

void create_output(char *output_path, bool toPrint, bool ofeat, char *ofeat_name, ofstream *allpassed, ofstream *oft, string mutants[]){
    // Create output files    
    
    int path_length = strlen(output_path); // get the path length and remove the trailing slash from the path if it was given with one
    if (output_path[path_length - 1] == '/') {
        output_path[--path_length] = '\0';
    }

    create_directory(output_path);    

    char output_file[path_length + 30]; // the buffer containing the output file names
    char ofeat_file[path_length + 30]; // the buffer containing the ofeat file names
    char output_name[] = "det-passed.csv";
    create_file_name(output_file, output_path, path_length, output_name);
    create_file_name(ofeat_file, output_path, path_length, ofeat_name);      

    create_file(output_file, allpassed);
    
    if (toPrint) {
        create_mutant_dir(output_path, mutants);
    }
    
    if (ofeat) {
        char ofeat[] = "set,per wt,amp wt,peak to trough wt,per delta,amp delta,peak to trough delta,per her1,amp her1,peak to trough her1,per her7,amp her7,peak to trough her7,per her13,amp her13,peak to trough her13,per her713,amp her713,peak to trough her713";
        create_ofeat(ofeat, ofeat_file, oft);
    }
}
