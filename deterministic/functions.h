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

#include "macros.h"

using namespace std;

#ifndef FUNC_H
#define FUNC_H

// escape codes to color the terminal output and shortcuts for common outputs (set -c or --no-color to disable these)

#define terminal_blue_d "\x1b[34m"
#define terminal_red_d "\x1b[31m"
#define terminal_reset_d "\x1b[0m"
#define terminal_done terminal_blue << "Done" << terminal_reset
#define terminal_no_memory terminal_red << "Not enough memory!" << terminal_reset


struct glevels {
    /*
     Structure for storing concentration levels.
     */
    int nfinal, cells;
 
    // matrices which will contain concentration levels for each protein and mRNA
    double** mh1;
    double** mh7;
    double** mh13;
    double** md;
    double** ph1;
    double** ph7;
    double** ph13;
    double** pd;
    double** ph11;
    double** ph17;
    double** ph113;
    double** ph77;
    double** ph713;
    double** ph1313;
    
    glevels(int, int);
    ~glevels();
};

struct rates {
    double rates_base[NUM_RATES]; // Base rates taken from the current parameter set
    double curr_rates[NUM_RATES]; // Current rates calculated using base rates and gradient factors
    bool using_gradients; // Whether or not any rates have specified perturbations
    double* factors_gradient[NUM_RATES]; // Gradients (as arrays of (step, percentage with 1=100%) pairs) taken from the gradients input file
    bool has_gradient[NUM_RATES]; // Whether each rate has a specified gradient
    int steps;
    //double* rates_active[NUM_RATES]; // Rates per cell position that factor in the base rates, each cell's perburations, and the gradients at each position
    
    explicit rates (int steps) {
        memset(this->rates_base, 0, sizeof(this->rates_base));
        memset(this->curr_rates, 0, sizeof(this->curr_rates));
        this->using_gradients = false;
        this->steps = steps;
        for (int i = 0; i < NUM_RATES; i++) {
            this->factors_gradient[i] = new double[steps];
            for (int j = 0; j < steps; j++) {
                this->factors_gradient[i][j] = 1;
            }
            this->has_gradient[i] = false;
            //this->rates_active[i] = new double[cells];
        }
    }
    
    ~rates () {
        for (int i = 0; i < NUM_RATES; i++) {
            delete[] this->factors_gradient[i];
            //delete[] this->rates_active[i];
        }
    }
    
};

struct data{
    /*
     Structure for storing oscillation features.
     */
    double period, amplitude, peaktotrough1, peaktotrough2;
    bool w;
};

bool checkPropensities(glevels*, rates, int, double);
void printForPlotting(string, glevels*, int, double);
void test_print(rates);
void store_values(char*, char*, int&, rates*, const int, int);
void clear_levels(glevels*, int, int);
bool model(double, int, glevels*, rates, double);
void ofeatures(glevels*, double, int, bool, data&);
void parseLine(char*, int[], int&);
void skipFirstLine(char*, int&);
void usage(const char*);
void licensing();
bool run_mutant(glevels*, int, double, rates, data&, bool, double, int, int);
void fill_rates(rates& rs, char *buffer, int *index);
void fill_gradients (rates& rs, char* gradients);
void print_rate(rates *rs);
void update_rate(rates& rs, int step);
void reset_rate(rates& rs);

inline void clear_data(data &d){
    d.period = 0.0;
    d.amplitude = 0.0;
    d.peaktotrough1 = 0.0;
    d.peaktotrough2 = 0.0;
    d.w = false;
}

/*
 Inline functions for computing mRNA transcription rates
 */
inline double fh1(double xh11, double xh713, double yd, double msh1, double critph1h1, double critph7h13, double critpd)
{
    // her1 mRNA
    double x11 = xh11 / critph1h1, x713 = xh713 / critph7h13, y = yd / critpd;
    return msh1 * ((1 + y) / (1 + y + x11 * x11 + x713 * x713));
}

inline double fh7(double xh11, double xh713, double yd, double msh7, double critph1h1, double critph7h13, double critpd)
{
    // her7 mRNA
    double x11 = xh11 / critph1h1, x713 = xh713 / critph7h13, y = yd/critpd;
    return msh7 * ((1+y) / (1 + y + x11 * x11 + x713 * x713));
}

inline double fd(double xh11, double xh713, double yd, double msd, double critph1h1, double critph7h13, double critpd)
{
    // delta mRNA
    double x11 = xh11 / critph1h1, x713 = xh713 / critph7h13;
    return msd / (1 + x11 * x11 + x713 * x713);
}

/*
 Inline functions for testing mutant conditions -- this is where you may change the condition ranges.
 */
inline bool fwildtype(double peaktotrough, double peaktotrough2){
    // Wild type
    return (peaktotrough2 >= 1.5 && peaktotrough >= 1.5 && (peaktotrough2 / peaktotrough)<=1.5);
}

inline bool f1_mutant(double h1period, double h1amplitude, double wperiod){
    // Her1 mutant
    return (((h1period / wperiod) > 0.97) && ((h1period / wperiod) < 1.03));
}

inline bool f7_mutant(double h7period, double h7amplitude, double wperiod){
    // Her7 mutant
    return (((h7period / wperiod) > 0.97) && ((h7period / wperiod) < 1.03));
}

inline bool f13_mutant(double h13period, double h13amplitude, double wperiod){
    // Her13 mutant
    return (((h13period / wperiod) > 1.03) && ((h13period / wperiod) < 1.09));
}

// Her7 and Her6 (13) mutant 
inline bool f713_mutant(double h713period, double h713amplitude, double wperiod){
    return (((h713period / wperiod) > 1.03) && ((h713period / wperiod) < 1.09));
}

inline bool fd_mutant(double dperiod, double damplitude, double wperiod){
    return (((dperiod / wperiod) > 1.04) && ((dperiod / wperiod) < 1.30));
}

#endif 
