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

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <math.h>

#include "functions.h"
#include "macros.h"

using namespace std;

// global variables used in main.cpp and functions.cpp
extern char* terminal_blue;
extern char* terminal_red;
extern char* terminal_reset;

void fill_rates(rates& rs, double items[]){
    for (int i = 0; i < NUM_RATES; i++){
        rs.rates_base[i] = items[i];
        rs.curr_rates[i] = items[i];
    }
}

void fill_gradients (rates& rs, char* gradients) {
    if (gradients != NULL) {
        static const char* usage_message = "There was an error reading the given gradients file.";
        int con; // The index of the concentration
        int step; // The column in the cell tissue
        double factor; // The factor to apply
        int last_step = 0; // The last column in the cell tissue given a gradient factor
        int start_step; // The first column to apply the gradient from
        int i = 0; // The index in the buffer
        while (gradients[i] != '\0') {
            // Read the concentration value
            if (sscanf(gradients + i, "%d", &con) != 1) {
                usage(usage_message);
            }
            if (con < 0 || con > NUM_RATES) {
                usage("The given gradients file includes rate indices outside of the valid range. Please adjust the gradients file or add the appropriate rates by editing the macros file and recompiling.");
            }
            rs.using_gradients = true; // Mark that at least one concentration has a gradient
            rs.has_gradient[con] = true; // Mark that this concentration has a gradient
            
            // Read every (position factor) pair
            while (gradients[i++] != ' ') {} // Skip past the concentration index
            while (not_EOL(gradients[i])) {
                // Read a position factor pair
                if (sscanf(gradients + i, "(%d %lf)", &step, &factor) != 2) {
                    usage(usage_message);
                }
                if (step < 0 || step >= rs.steps) {
                    usage("The given gradients file includes positions outside of the given simulation width. Please adjust the gradients file or increase the width of the simulation using -x or --total-width.");
                }
                if (factor < 0) {
                    usage("The given gradients file includes factors less than 0. Please adjusted the gradients file.");
                }
                
                // Apply the gradient factor
                factor /= 100;
                start_step = last_step;
                last_step = step;
                for (int j = start_step + 1; j < step; j++) {
                    rs.factors_gradient[con][j] = interpolate(j, start_step, step, rs.factors_gradient[con][start_step], factor);
                }
                rs.factors_gradient[con][step] = factor;
                while (gradients[i++] != ')') {} // Skip past the end of the pair
                while (gradients[i] == ' ') {i++;} // Skip any whitespace before the next pair
            }
            
            // Apply the last gradient factor to the rest of the steps
            for (int j = step + 1; j < rs.steps; j++) {
                rs.factors_gradient[con][j] = rs.factors_gradient[con][step];
            }
            i++;
        }
    }
}

//Print rates structure for testing purpose
void print_rate(rates *rs){ 
    for (int i = 0; i < NUM_RATES; i++){
        cout << "Rate " << i << ": " << rs->rates_base[i] << endl;      
        if (rs->has_gradient[i]){
            for (int j = 0; j < rs->steps; j++){
                cout << "Step " << j << ": " << rs->factors_gradient[i][j] << endl;
            }               
        } 
    }
}

void update_rate(rates& rs, int step){
    for (int i = 0; i < NUM_RATES; i++){
        if (rs.has_gradient[i]){
            rs.curr_rates[i] = rs.rates_base[i]*rs.factors_gradient[i][step];
        }
    }
}

void printForPlotting(string file, glevels* x, int nfinal, double eps) {
    /*
     Prints the concentrations of her1 mRNA into a specified output file.
     The concentrations are printed every 0.1 minutes of the simulation.
     */
    ofstream plot;
    cerr << file.c_str() << endl;
    try {
        plot.open(file.c_str(), fstream::out);
    } catch (ofstream::failure) {
        cerr << terminal_red << "Couldn't open file " << file << " for plotting! Exit status 1." << terminal_reset << endl;
        exit(1);
    }
    int step = (int) (0.1 / eps);
    for (int n = 0; n < nfinal; n += step) {
        plot << n << " " << x->mh1[0][n] << endl;
    }
    plot.close();
}


void usage (const char* message) {
    if (strcmp(message, "") != 0) { // if there is an error message to print then print it
        cout << terminal_red << message << terminal_reset << endl << endl;
    }
    cout << "Usage: [-option [value]]... [--option [value]]..." << endl;
    cout << "-x, --width        : the tissue width (in cells), 2 for two-cell system, min=3 for chain, min=4 and even for tissue, default=3" << endl;
    cout << "-y, --height       : the tissue height (in cells), 1 for two-cell system and chain, min=4 and even for tissue, default=1" << endl;
    cout << "-e, --epsilon      : the size of the timestep to be used for solving the DDEs using Euler's method" << endl;
    cout << "-m, --minutes      : the maximum number of minutes to simulate before ending, min=1, default=1200" << endl;
    cout << "-f, --ofeatures    : the path and file in which to print oscillation features" << endl;
    cout << "-p, --parameters   : the number of parameters for which to simulate the model, min=1, default=1" << endl;
    cout << "-s, --seed         : the seed to generate random numbers, min=1, default=time" << endl;
    cout << "-i, --input        : the input path and file to accept parameters from, default=input.txt" << endl;
    cout << "-o, --output       : the path and file to print the output (i.e. parameters which passed conditions) to, default=det-allpassed.csv" << endl;
    cout << "-a, --propensities : the threshold for the propensity functions which could be used in the stochastic simulation, min=1, default=none" << endl;
    cout << "-w, --write        : print the concentrations of the simulations to file, default=unused" << endl;
    cout << "-c, --no-color     : disable coloring the terminal output, default=unused" << endl;
    cout << "-q, --quiet        : hide the terminal output, default=unused" << endl;
    cout << "-l, --licensing    : view licensing information (no simulations will be run)" << endl;
    cout << "-h, --help         : view usage information (i.e. this)" << endl;
    cout << endl << terminal_blue << "Example: ./" << terminal_reset << endl << endl;
    exit(0);
}

void licensing () {
    /*
     Licensing information.
     */
    cout << "Stochastic simulator for zebrafish segmentation" << endl;
    cout << "Copyright (C) 2012 Ahmet Ay (aay@colgate.edu), Jack Holland (jholland@colgate.edu), Adriana Sperlea (asperlea@colgate.edu)" << endl;
    cout << "This program comes with ABSOLUTELY NO WARRANTY" << endl;
    cout << "This is free software, and you are welcome to redistribute it under certain conditions;" << endl;
    cout << "You can use this code and modify it as you wish under the condition that you refer to the article: ???" << endl;
    exit(0);
}


bool checkPropensities(glevels *x, rates pars, int sn, double CUTOFF) {
    /*
     Checks that the propensity functions which could be used in a stochastic simulation do not go over the set CUTOFF.
     Returns true if all propensities are < CUTOFF and false otherwise.
     */
    if (pars.curr_rates[RPSH1]       * x->mh1[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RPDH1]       * x->ph1[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RDAH1H1]     * x->ph1[0][sn] * (x->ph1[0][sn] - 1) / 2 > CUTOFF) return false;
    if (pars.curr_rates[RDDIH1H1]    * x->ph11[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RDAH1H7]     * x->ph1[0][sn] * x->ph7[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RDDIH1H7]    * x->ph17[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RDAH1H13]    * x->ph1[0][sn] * x->ph13[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RDDIH1H13]   * x->ph113[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RPSH7]       * x->mh7[0][sn] > CUTOFF) return false;                                 
    if (pars.curr_rates[RPDH7]       * x->ph7[0][sn] > CUTOFF) return false;                                 
    if (pars.curr_rates[RDAH7H7]     * x->ph7[0][sn] * (x->ph7[0][sn] - 1) / 2 > CUTOFF) return false;      
    if (pars.curr_rates[RDDIH7H7]    * x->ph77[0][sn] > CUTOFF) return false;                               
    if (pars.curr_rates[RDAH7H13]    * x->ph7[0][sn]  * x->ph13[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RDDIH7H13]   * x->ph713[0][sn] > CUTOFF) return false;                              
    if (pars.curr_rates[RPSH13]      * x->mh13[0][sn] > CUTOFF) return false;                               
    if (pars.curr_rates[RPDH13]      * x->ph13[0][sn] > CUTOFF) return false;                               
    if (pars.curr_rates[RDAH13H13]   * x->ph13[0][sn] * (x->ph13[0][sn] - 1) / 2  > CUTOFF) return false;   
    if (pars.curr_rates[RDDIH13H13]  * x->ph1313[0][sn] > CUTOFF) return false;                             
    if (pars.curr_rates[RDDGH1H1]    * x->ph11[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RDDGH1H7]    * x->ph17[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RDDGH1H13]   * x->ph113[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RDDGH7H7]    * x->ph77[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RDDGH7H13]   * x->ph713[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RDDGH13H13]  * x->ph1313[0][sn] > CUTOFF) return false;
    if (pars.curr_rates[RPSDELTA]    * x->md[0][sn] > CUTOFF) return false;                                                         
    if (pars.curr_rates[RPDDELTA]    * x->pd[0][sn] > CUTOFF) return false;                                                         
    if (fh1(x->ph11[0][sn], x->ph713[0][sn], x->pd[1][sn], pars.curr_rates[RMSH1], pars.curr_rates[RCRITPH1H1], pars.curr_rates[RCRITPH7H13], pars.curr_rates[RCRITPDELTA]) > CUTOFF) return 0; 
    if (pars.curr_rates[RMDH1]       * x->mh1[0][sn] > CUTOFF) return false;                                                            
    if (fh7(x->ph11[0][sn], x->ph713[0][sn], x->pd[1][sn], pars.curr_rates[RMSH7], pars.curr_rates[RCRITPH1H1], pars.curr_rates[RCRITPH7H13], pars.curr_rates[RCRITPDELTA]) > CUTOFF) return 0;
    if (pars.curr_rates[RMDH7]       * x->mh7[0][sn] > CUTOFF) return false;                                                        
    if (pars.curr_rates[RPSH13] > CUTOFF) return 0;                                                                                 
    if (pars.curr_rates[RMDH13]      * x->mh13[0][sn] > CUTOFF) return false;                                                           
    if (fd(x->ph11[0][sn], x->ph713[0][sn], x->pd[1][sn], pars.curr_rates[RMSDELTA], pars.curr_rates[RCRITPH1H1], pars.curr_rates[RCRITPH7H13], pars.curr_rates[RCRITPDELTA]) > CUTOFF) return 0;   
    if (pars.curr_rates[RMDDELTA]    * x->md[0][sn] > CUTOFF) return false;
    return true;
}

void parseLine(char* buffer, double items[], int &index)
{
    /*
     Parses a line of the file buffer into the items array, to then be transfered in the rateValues structure.
     */
    char db[50];
    db[0] = '\0';
    int itindex = 0;
    int i;
    
    // Check every character in the buffer
    for (i = index; buffer[i] != '\n' && buffer[i] != '\0'; i++) {
        // When a comma is found, it marks the end of one number
        if (buffer[i] == ',') {
            // Convert it into a double and add it to index
            items[itindex ++] = atof(db);
            db[0] = '\0';
        } else {
            // If the current character is not a comma then continue adding digits to the number currently being created
            int orsize = strlen(db);
            db[orsize] = buffer[i];
            db[orsize + 1] = '\0';
        }
    }
    index = i + 1;
    items[itindex] = atof(db);
}

double make_random(double range[])
{
    /*
     Creates and returns a random real number in the interval [range[0], range[1]]
     */
    return double((range[1] - range[0]) * rand() / (RAND_MAX + 1.0) + range[0]);
}

void generate_set(double items[])
{
    /*
     Generate a random set of parameters to be used with the simulation.
     The ranges for every parameter can be adjusted in their specific 2 element arrays.
     */
    // ranges of protein synthesis rates
    double psh1[2] = {30, 60}, psh7[2] = {10, 57}, psh13[2] = {27, 57}, psd[2] = {22, 59};
    items[0] = make_random(psh1); 
    items[1] = make_random(psh7);
    items[2] = make_random(psh13); 
    items[3] = make_random(psd);  
    
    // ranges of protein degradation rates
    double pdh1[2] = {0.12, 0.37}, pdh7[2] = {0.11, 0.4}, pdh13[2] = {0.11, 0.39}, pdd[2] = {0.15, 0.38};
    items[4] = make_random(pdh1); 
    items[5] = make_random(pdh7);
    items[6] = make_random(pdh13);
    items[7] = make_random(pdd); 
    
    // ranges of mRNA synthesis rates
    double msh1[2] = {32, 63}, msh7[2] = {34, 62}, msh13[2] = {31, 62}, msd[2] = {31, 65};
    items[8] = make_random(msh1); 
    items[9] = make_random(msh7); 
    items[10] = make_random(msh13); 
    items[11] = make_random(msd); 

    // ranges of mRNA degradation rates
    double mdh1[2] = {0.2, 0.38}, mdh7[2] = {0.28, 0.4}, mdh13[2] = {0.13, 0.39}, mdd[2] = {0.12, 0.39};
    items[12] = make_random(mdh1); 
    items[13] = make_random(mdh7);  
    items[14] = make_random(mdh13);
    items[15] = make_random(mdd);
    
    // ranges of dimer degradation rates
    double ddgh1h1[2] = {0.25, 0.4}, ddgh1h7[2] = {0.16, 0.34}, ddgh1h13[2] = {0.1, 0.36}, ddgh7h7[2] = {0.12, 0.4}, ddgh7h13[2] = {0.26, 0.4}, ddgh13h13[2] = {0.11, 0.34};
    items[16] = make_random(ddgh1h1); 
    items[17] = make_random(ddgh1h7); 
    items[18] = make_random(ddgh1h13); 
    items[19] = make_random(ddgh7h7); 
    items[20] = make_random(ddgh7h13); 
    items[21] = make_random(ddgh13h13); 
    
    // ranges of mRNA transcription delays
    double delaymh1[2] = {8.8, 12.0}, delaymh7[2] = {8.6, 11.6}, delaymd[2] = {6.1, 12.0};
    items[22] = make_random(delaymh1);
    items[23] = make_random(delaymh7);
    items[24] = -1;
    items[25] = make_random(delaymd);
    
    // ranges mRNA translation delays
    double delayph1[2] = {0.8, 2.0}, delayph7[2] = {0.4, 1.8}, delayph13[2] = {0.6, 1.8}, delaypd[2] = {10, 18};
    items[26] = make_random(delayph1);
    items[27] = make_random(delayph7);
    items[28] = make_random(delayph13);
    items[29] = make_random(delaypd);
    
    // ranges of dimer association rates
    double dah1h1[2] = {0.005, 0.03}, dah1h7[2] = {0.0006, 0.009}, dah1h13[2] = {0.006, 0.029}, dah7h7[2] = {0.002, 0.024}, dah7h13[2] = {0.007, 0.03}, dah13h13[2] = {0.001, 0.016};
    items[30] = make_random(dah1h1);
    items[32] = make_random(dah1h7);
    items[34] = make_random(dah1h13);
    items[36] = make_random(dah7h7);
    items[38] = make_random(dah7h13);
    items[40] = make_random(dah13h13);
    
    double ddh1h1[2] = {0.06, 0.3}, ddh1h7[2] = {0.03, 0.28}, ddh1h13[2] = {0.004, 0.18}, ddh7h7[2] = {0.07, 0.3}, ddh7h13[2] = {0.03, 0.3}, ddh13h13[2] = {0.05, 0.29};
    items[31] = make_random(ddh1h1);
    items[33] = make_random(ddh1h7);
    items[35] = make_random(ddh1h13);
    items[37] = make_random(ddh7h7);
    items[39] = make_random(ddh7h13);
    items[41] = make_random(ddh13h13);
    
    // critical number of Her1-Her1, Her7-Her13 dimer proteins for inhibition of transcription
    double critph1h1[2] = {160, 720}, critph7h13[2] = {200, 920};
    items[42] = make_random(critph1h1);
    items[43] = make_random(critph7h13);
    
    // critical Delta proteins for inhibition 
    double critpd[2] = {240, 720};
    items[44] = make_random(critpd);
}

void store_values(char* input_file, char* buffer, int &index, rates *rateValues, const int STEP, int seed, char *gradients){
    /*
     Stores sets of parameters into the rateValues structure.
     The sets are either taken from the input file given by the user, or generated according to a random seed.
     */
    double items[45];
    for (int i = 0; i < STEP; i++) {
        if (input_file != NULL) {
            // If an input file was specified, parse a line from it.
            parseLine(buffer, items, index);
        } else {
            // If there is no input file, generate a random set and print the seed to a file so that results can be replicated later.
            generate_set(items);
            ofstream seedfile;
            seedfile.open("seed.txt", fstream::out);
            seedfile << seed;
            seedfile.close();
        }
        fill_rates(rateValues[i], items);
        if (gradients != NULL){
            fill_gradients(rateValues[i], gradients);
        }
    }
}

void clear_levels(glevels *old, int nfinal, int cells){
    /*
     Clears concentrations from previous simulations.
     */
    for (int i = 0; i < cells; i++) {
        for(int j = 0; j < nfinal; j++) {
            old->mh1[i][j] = 0, old->ph1[i][j] = 0;
            old->mh7[i][j] = 0, old->ph7[i][j] = 0;
            old->mh13[i][j] = 0, old->ph13[i][j] = 0;
            old->md[i][j] = 0, old->pd[i][j] = 0;
            old->ph11[i][j] = 0, old->ph17[i][j] = 0;
            old->ph113[i][j] = 0, old->ph77[i][j] = 0;
            old->ph713[i][j] = 0, old->ph1313[i][j] = 0;
        }
    }   
}

int fix(int x, int end)
{
    if (x < 0) return (end - 1);
    if (x >= end) return 0;
    return x;
}

bool model(double eps, int nfinal, glevels *g, rates r, double max_prop, int columns, int rows){
    /*
     Runs the deterministic simulation of the model.
     For each time step:
     1) Iterate through every cell (2 cells in this case) and update the concentrations of proteins and mRNA.
        These concentration values are obtained by solving the differential equations for that time step, using Euler's method.
     2) Check that the concentrations do not become negative -- a negative amount of protein is not biologically sensible
     3) Check that the propensity functions do not go above the set threshold -- if one was specified.
     */
    
    // Convert the time delay values to integers, because the deterministic simulation uses discrete time points.
    int ndelaymd, ndelayph1, ndelayph7, ndelayph13, ndelaypd, ndelaymh1, ndelaymh7; 
    ndelayph1 = int(r.curr_rates[RDELAYPH1]/eps);
    ndelayph7 = int(r.curr_rates[RDELAYPH7]/eps);
    ndelayph13 = int(r.curr_rates[RDELAYPH13]/eps);
    ndelaypd = int(r.curr_rates[RDELAYPDELTA]/eps);
    ndelaymd = int(r.curr_rates[RDELAYMDELTA]/eps);
    ndelaymh1 = int(r.curr_rates[RDELAYMH1]/eps);
    ndelaymh7 = int(r.curr_rates[RDELAYMH7]/eps);

    int nmh1, nph1; 
    int nmh7, nph7;
    int nph13;
    int nmd, npd;

    int cells = rows * columns;
    int last_step = 1; //last step when we recalculate the rates based on gradient factors
    int step_size = nfinal/50; //the distance between 2 points we recalculate the rates
    for (int n = 1; n < nfinal; n++) {
        if ((n-last_step) >= step_size){
            update_rate(rs, last_step/step_size + 1);
            last_step = n;
        }
        for (int i = 0; i < cells; i++) {
            nmh1 = ndelaymh1, nph1 = ndelayph1;
            nmh7 = ndelaymh7, nph7 = ndelayph7;
            nph13 = ndelayph13;
            nmd = ndelaymd, npd = ndelaypd;
    
            //Protein synthesis
            g->ph1[i][n] = g->ph1[i][n - 1] + eps * ((n > nph1 ? r.curr_rates[RPSH1] * g->mh1[i][n - nph1]:0)-r.curr_rates[RPDH1]*g->ph1[i][n - 1]-2*r.curr_rates[RDAH1H1]*g->ph1[i][n - 1]*g->ph1[i][n - 1]+2*r.curr_rates[RDDIH1H1]*g->ph11[i][n - 1]-r.curr_rates[RDAH1H7]*g->ph1[i][n - 1]*g->ph7[i][n - 1]+r.curr_rates[RDDIH1H7]*g->ph17[i][n - 1]-r.curr_rates[RDAH1H13]*g->ph1[i][n - 1]*g->ph13[i][n - 1]+r.curr_rates[RDDIH1H13]*g->ph113[i][n - 1]);
            g->ph7[i][n] = g->ph7[i][n - 1] + eps * ((n > nph7 ? r.curr_rates[RPSH7]*g->mh7[i][n - nph7]:0)-r.curr_rates[RPDH7]*g->ph7[i][n - 1]-2*r.curr_rates[RDAH7H7]*g->ph7[i][n - 1]*g->ph7[i][n - 1]+2*r.curr_rates[RDDIH7H7]*g->ph77[i][n - 1]-r.curr_rates[RDAH1H7]*g->ph1[i][n - 1]*g->ph7[i][n - 1]+r.curr_rates[RDDIH1H7]*g->ph17[i][n - 1]-r.curr_rates[RDAH7H13]*g->ph7[i][n - 1]*g->ph13[i][n - 1]+r.curr_rates[RDDIH7H13]*g->ph713[i][n - 1]);
            g->ph13[i][n] = g->ph13[i][n - 1] + eps * ((n > nph13 ? r.curr_rates[RPSH13]*g->mh13[i][n - nph13]:0)-r.curr_rates[RPDH13]*g->ph13[i][n - 1]-2*r.curr_rates[RDAH13H13]*g->ph13[i][n - 1]*g->ph13[i][n - 1]+2*r.curr_rates[RDDIH13H13]*g->ph1313[i][n - 1]-r.curr_rates[RDAH1H13]*g->ph1[i][n - 1]*g->ph13[i][n - 1]+r.curr_rates[RDDIH1H13]*g->ph113[i][n - 1]-r.curr_rates[RDAH7H13]*g->ph7[i][n - 1]*g->ph13[i][n - 1]+r.curr_rates[RDDIH7H13]*g->ph713[i][n - 1]);
            if (g->ph1[i][n] < 0 || g->ph7[i][n] < 0 || g->ph13[i][n] < 0) {
                return false;
            }

            //Dimer proteins
            g->ph11[i][n] = g->ph11[i][n - 1] + eps * (r.curr_rates[RDAH1H1]*g->ph1[i][n - 1]*g->ph1[i][n - 1]-r.curr_rates[RDDIH1H1]*g->ph11[i][n - 1]-r.curr_rates[RDDGH1H1]*g->ph11[i][n - 1]);
            g->ph17[i][n] = g->ph17[i][n - 1] + eps * (r.curr_rates[RDAH1H7]*g->ph1[i][n - 1]*g->ph7[i][n - 1]-r.curr_rates[RDDIH1H7]*g->ph17[i][n - 1]-r.curr_rates[RDDGH1H7]*g->ph17[i][n - 1]);
            g->ph113[i][n] = g->ph113[i][n - 1] + eps * (r.curr_rates[RDAH1H13]*g->ph1[i][n - 1]*g->ph13[i][n - 1]-r.curr_rates[RDDIH1H13]*g->ph113[i][n - 1]-r.curr_rates[RDDGH1H13]*g->ph113[i][n - 1]);
            g->ph77[i][n] = g->ph77[i][n - 1] + eps * (r.curr_rates[RDAH7H7]*g->ph7[i][n - 1]*g->ph7[i][n - 1]-r.curr_rates[RDDIH7H7]*g->ph77[i][n - 1]-r.curr_rates[RDDGH7H7]*g->ph77[i][n - 1]);
            g->ph713[i][n] = g->ph713[i][n - 1] + eps * (r.curr_rates[RDAH7H13]*g->ph7[i][n - 1]*g->ph13[i][n - 1]-r.curr_rates[RDDIH7H13]*g->ph713[i][n - 1]-r.curr_rates[RDDGH7H13]*g->ph713[i][n - 1]);
            g->ph1313[i][n] = g->ph1313[i][n - 1] + eps * (r.curr_rates[RDAH13H13]*g->ph13[i][n - 1]*g->ph13[i][n - 1]-r.curr_rates[RDDIH13H13]*g->ph1313[i][n - 1]-r.curr_rates[RDDGH13H13]*g->ph1313[i][n - 1]);
                        
            // Delta Protein
            g->pd[i][n] = g->pd[i][n - 1] + eps*((n > npd ? r.curr_rates[RPSDELTA]*g->md[i][n-npd]:0)-r.curr_rates[RPDDELTA]*g->pd[i][n - 1]);
            
            if (g->ph11[i][n] < 0 || g->ph17[i][n] < 0 || g->ph113[i][n] < 0 || g->ph77[i][n] < 0 || g->ph713[i][n] < 0 || g->ph1313[i][n] < 0 || g->pd[i][n] < 0) {
                return false;
            }
            
            /*
             Compute the value of the delta protein coming from the neighbors for each cell. In two-cell systems, both cells are neighbors of each other. Chains wrap horizontally and hexagonal
             tissue grids wrap horizontally and vertically like a honeycomb.
             Two-cell systems look like this:
             ___  ___
             /   \/   \             where 1 and 2 are neighbors of each other
             | 1 || 2 |
             \___/\___/
             
             Chains of cells look like this:
             ___  ___  ___  ___
             /   \/   \/   \/   \   where x has neighbors n
             | n || x || n ||   |
             \___/\___/\___/\___/
             
             Tissues of cells look like this:
              ___  ___  ___  ___
             /   \/   \/   \/   \   where x has neighbors n
             |   || n || n ||   |
             \___/\___/\___/\___/_
             /   \/   \/   \/   \
             | n || x || n ||   |
             \___/\___/\___/\___/
             /   \/   \/   \/   \
             |   || n || n ||   |
             \___/\___/\___/\___/
             /   \/   \/   \/   \
             |   ||   ||   ||   |
             \___/\___/\___/\___/
             */
            double avgpdh1 = 0, avgpdh7 = 0, avgpdd = 0;
            if (rows == 1) {
                if (columns == 2) {
                    // If there are only two cells, the only delta input is coming from the other cell
                    avgpdh1 = g->pd[1 - i][n - nmh1];
                    avgpdh7 = g->pd[1 - i][n - nmh7];
                    avgpdd = g->pd[1 - i][n - nmd];
                } else {
                    // If there is a chain of cells, the delta input is coming from the two cells to the left and right
                    int left = fix(i - 1, columns), right = fix(i + 1, columns);
                    avgpdh1 = (g->pd[left][n - nmh1] + g->pd[right][n - nmh1]) / 2;
                    avgpdh7 = (g->pd[left][n - nmh7] + g->pd[right][n - nmh7]) / 2;
                    avgpdd = (g->pd[left][n - nmd] + g->pd[right][n - nmd]) / 2;
                }
            } else {
                int curi = i / columns, curj = i % columns, xx, yy;
                int dx[] = {0, 0, -1, -1, 1, 1};
                if (curi % 2 == 0) {
                    int dy[] = {1, -1, -1, 0, -1, 0};
                    for (int k = 0; k < 6; k++) {
                        xx = curi + dx[k];
                        xx = fix(xx, rows); // check for wrapping
                        yy = curj + dy[k];
                        yy = fix(yy, columns); // check for wrapping;
                        int newi = xx * columns + yy;
                    
                        // add the delta of all neighbours
                        if (n > nmh1) avgpdh1 += g -> pd[newi][n - nmh1];
                        if (n > nmh7) avgpdh7 += g -> pd[newi][n - nmh7];
                        if (n > nmd) avgpdd += g -> pd[newi][n - nmd];
                    }
                } else {
                    int dy[] = {1, -1, 0, 1, 0, 1};
                    for (int k = 0; k < 6; k++) {
                        xx = curi + dx[k];
                        xx = fix(xx, rows);
                        yy = curj + dy[k];
                        yy = fix(yy, columns);
                        int newi = xx * columns + yy;
                    
                        // add the delta of all neighbours
                        if (n > nmh1) avgpdh1 += g -> pd[newi][n - nmh1];
                        if (n > nmh7) avgpdh7 += g -> pd[newi][n - nmh7];
                        if (n > nmd) avgpdd += g -> pd[newi][n - nmd];
                    }
                }
                avgpdh1 /= 6;
                avgpdh7 /= 6;
                avgpdd /= 6;
            }
            
            // mRNA Synthesis
            g->mh1[i][n] = g->mh1[i][n - 1] + eps * ((n > nmh1 ? fh1(g->ph11[i][n - nmh1], g->ph713[i][n - nmh1], avgpdh1, r.curr_rates[RMSH1], r.curr_rates[RCRITPH1H1], r.curr_rates[RCRITPH7H13], r.curr_rates[RCRITPDELTA]):fh1(0, 0, 0, r.curr_rates[RMSH1], r.curr_rates[RCRITPH1H1], r.curr_rates[RCRITPH7H13], r.curr_rates[RCRITPDELTA]))-r.curr_rates[RMDH1]*g->mh1[i][n - 1]);
            g->mh7[i][n] = g->mh7[i][n - 1] + eps * ((n > nmh7 ? fh7(g->ph11[i][n - nmh7], g->ph713[i][n - nmh7], avgpdh7, r.curr_rates[RMSH7], r.curr_rates[RCRITPH1H1], r.curr_rates[RCRITPH7H13], r.curr_rates[RCRITPDELTA]):fh7(0, 0, 0, r.curr_rates[RMSH7], r.curr_rates[RCRITPH1H1], r.curr_rates[RCRITPH7H13], r.curr_rates[RCRITPDELTA]))-r.curr_rates[RMDH7]*g->mh7[i][n - 1]);
            g->mh13[i][n] = g->mh13[i][n - 1] + eps * (r.curr_rates[RMSH13]-r.curr_rates[RMDH13]*g->mh13[i][n - 1]); 
            g->md[i][n] = g->md[i][n - 1] + eps * ((n > nmd ? fd(g->ph11[i][n - nmd], g->ph713[i][n - nmd], avgpdd, r.curr_rates[RMSDELTA], r.curr_rates[RCRITPH1H1], r.curr_rates[RCRITPH7H13], r.curr_rates[RCRITPDELTA]):fd(0, 0, 0, r.curr_rates[RMSDELTA], r.curr_rates[RCRITPH1H1], r.curr_rates[RCRITPH7H13], r.curr_rates[RCRITPDELTA]))-r.curr_rates[RMDDELTA]*g->md[i][n - 1]);
            if (g->mh1[i][n] < 0 || g->mh7[i][n] < 0 || g-> mh13[i][n] < 0 || g-> md[i][n] < 0) {
                return false;
            }
            if (max_prop != INFINITY && !checkPropensities(g, r, n, max_prop)) {
                return false;
            }
        }
    }
    return true;
}

bool run_mutant(glevels *g, int t_steps, double eps, rates temp_rate, data &of, bool wild, double max_prop, int x, int y)
{
    /*
     Performs the steps necessary in the simulation and analysis of the wild type or a certain mutant.
     1) Clear the levels from the previous simulation.
     2) Run the model for the specified duration.
     3) Create the oscillation features of the simulation.
     4) Return whether concentrations in the model where positive values below the propensity threshold
     */
    clear_levels(g, t_steps, x * y);
    bool pass = model(eps, t_steps, g, temp_rate, max_prop, x, y);
    ofeatures(g, eps, t_steps, wild, of); 
    return pass;
}
/
void ofeatures(glevels *g, double eps, int nfinal, bool wild, data &d) {
    /*
     Calculates the oscillation features -- period, amplitude and peak to trough
     ratio for a set of concentration levels.
     The values are calculated using the last peak and trough of the oscillations,
     since the amplitude of the first few oscillations can be slightly unstable.
     For the wild type, the peak and trough at the middle of the graph are also calculated
     in order to ensure that the oscillations are sustained.
     */
    double tmaxlast = 0, tmaxpenult = 0, mmaxlast = 0, mmaxpenult = 0;
    double tminlast = 0, tminpenult = 0, mminlast = 0, mminpenult = 0;
        
    for (int n = 1; n < nfinal - 1; n++) {
        // check if the current point is a peak
        if (g -> mh1[0][n + 1] < g -> mh1[0][n] && g -> mh1[0][n] > g -> mh1[0][n - 1]) {
            tmaxpenult = tmaxlast;
            tmaxlast = n*eps;
            mmaxpenult = mmaxlast;
            mmaxlast = g -> mh1[0][n];
        }
        
        // check if the current point is a trough
        if (g -> mh1[0][n + 1] > g -> mh1[0][n] && g -> mh1[0][n] < g -> mh1[0][n - 1]){
            tminpenult = tminlast;
            tminlast = n * eps;
            mminpenult = mminlast;
            mminlast = g -> mh1[0][n];
        }
    }
    
    if (wild) {
        // calculate the peak and trough at the middle of the graph
        double mmaxlast2 = 0.0, mminlast2 = 0.0;
        for(int m = 2; m < nfinal/2; m++){
            if(g -> mh1[0][m + 1] < g -> mh1[0][m] && g -> mh1[0][m] > g -> mh1[0][m - 1]){
                mmaxlast2 = g -> mh1[0][m];
            }
            if(g -> mh1[0][m + 1] > g -> mh1[0][m] && g -> mh1[0][m] < g -> mh1[0][m - 1]){
                mminlast2 = g -> mh1[0][m];
            }
        
            // in order to avoid dividing by zero in case a trough is 0, set it to 1
            if(mminlast2 == 0.0 || mminlast == 0.0) {
                mminlast2 = 1.0;
                mminlast = 1.0;
            }
            d.peaktotrough2 = mmaxlast2/mminlast2;
        }
    }
    d.period = tmaxlast-tmaxpenult;
    d.amplitude = mmaxlast-mminlast;
    d.peaktotrough1 = mmaxlast/mminlast;
}

glevels::glevels (int nfinal, int cells) {
    /*
     Constructor for the structure in which the concentration levels will be stored.
     Each array has two elements because the simulation will be run for two cells.
     */
    this->nfinal = nfinal;
    this->cells = cells;
    
    mh1 = new double*[cells];
    mh7 = new double*[cells];
    mh13 = new double*[cells];
    md = new double*[cells];
    ph1 = new double*[cells];
    ph7 = new double*[cells];
    ph13 = new double*[cells];
    pd = new double*[cells];
    ph11 = new double*[cells];
    ph17 = new double*[cells];
    ph113 = new double*[cells];
    ph77 = new double*[cells];
    ph713 = new double*[cells];
    ph1313 = new double*[cells];
    
    for (int i = 0; i < cells; i++) {
        mh1[i] = new double[nfinal];
        mh7[i] = new double[nfinal];
        mh13[i] = new double[nfinal];
        md[i] = new double[nfinal];
        ph1[i] = new double[nfinal];
        ph7[i] = new double[nfinal];
        ph13[i] = new double[nfinal];
        pd[i] = new double[nfinal];
        ph11[i] = new double[nfinal];
        ph17[i] = new double[nfinal];
        ph113[i] = new double[nfinal];
        ph77[i] = new double[nfinal];
        ph713[i] = new double[nfinal];
        ph1313[i] = new double[nfinal];
    }
}

glevels::~glevels () {
    /*
     Destructor for the structure containing the concentration levels.
     */
    for (int i = 0; i < this->cells; i++) {
        delete[] mh1[i];
        delete[] mh7[i];
        delete[] mh13[i];
        delete[] md[i];
        delete[] ph1[i];
        delete[] ph7[i];
        delete[] ph13[i];
        delete[] pd[i];
        delete[] ph11[i];
        delete[] ph17[i];
        delete[] ph113[i];
        delete[] ph77[i];
        delete[] ph713[i];
        delete[] ph1313[i];
    }
    
    delete[] mh1;
    delete[] mh7;
    delete[] mh13;
    delete[] md;
    delete[] ph1;
    delete[] ph7;
    delete[] ph13;
    delete[] pd;
    delete[] ph11;
    delete[] ph17;
    delete[] ph113;
    delete[] ph77;
    delete[] ph713;
    delete[] ph1313;
}


