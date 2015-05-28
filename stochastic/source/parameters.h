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

using namespace std;

#ifndef PARAMETERS_H
#define PARAMETERS_H

struct rates {
	double rates_base[NUM_RATES]; // Base rates taken from the current parameter set
	bool using_gradients; // Whether or not any rates have specified perturbations
	double* factors_gradient[NUM_RATES]; // Gradients (as arrays of (position, percentage with 1=100%) pairs) taken from the gradients input file
	bool has_gradient[NUM_RATES]; // Whether each rate has a specified gradient
	double* rates_active[NUM_RATES]; // Rates per cell position that factor in the base rates, each cell's perburations, and the gradients at each position
	
	explicit rates (int steps) {
		memset(this->rates_base, 0, sizeof(this->rates_base));
		this->using_gradients = false;
		for (int i = 0; i < NUM_RATES; i++) {
			this->factors_gradient[i] = new double[steps];
			for (int j = 0; j < steps; j++) {
				this->factors_gradient[i][j] = 1;
			}
			this->has_gradient[i] = false;
			this->rates_active[i] = new double[cells];
		}
	}
	
	~rates () {
		for (int i = 0; i < NUM_RATES; i++) {
			delete[] this->factors_gradient[i];
			delete[] this->rates_cell[i];
			delete[] this->rates_active[i];
		}
	}

	ostream& operator<<(ostream &strm){
		for (int i = 0; i < NUM_RATES; i++){
			strm << "Rate 1: " << this->rates_base[i] << endl;
			if (this->using_gradients){
				strm << "Gradient 1: " << this->factors_gradient[i] << endl;
			} 
		}
	}
};

#endif

