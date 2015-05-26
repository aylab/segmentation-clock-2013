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

#ifndef PARAMETERS_H
#define PARAMETERS_H

struct parameters {
	// protein synthesis rates per mRNA molecule for HerX & Delta for the two cells
	double psh1;
	double psh7;
	double psh13;
	double psd;
	
	// protein degradation rates i.e. inverse lifetimes for HerX & Delta for the two cells
	double pdh1;
	double pdh7;
	double pdh13;
	double pdd;
	
	// mRNA synthesis rates for HerX & Delta for the two cells
	double msh1;
	double msh7;
	double msh13;
	double msd;
	
	// mRNA degradation rates i.e. inverse lifetimes for HerX & Delta for the two cells
	double mdh1;
	double mdh7;
	double mdh13;
	double mdd;

	// dimer degradation rates
	double ddgh1h1;
	double ddgh1h7;
	double ddgh1h13;
	double ddgh7h7;
	double ddgh7h13;
	double ddgh13h13;
	
	double delaymh1; // time to make a molecule of her1 mRNA
	double delaymh7; // time to make a molecule of her7 mRNA
	double delaymh13; // time to make a molecule of her13 mRNA
	double delaymd; // time to make a molecule of delta mRNA
	double delayph1; // time to make a molecule of her1 protein
	double delayph7; // time to make a molecule of her7 protein
	double delayph13; // time to make a molecule of her13 protein
	double delaypd; // time to make a molecule of delta protein (mature delivered to the cell surface, activating Notch) 
	
	// dimer association rates
	double dah1h1;
	double dah1h7;
	double dah1h13;
	double dah7h7;
	double dah7h13;
	double dah13h13;
	
	// dimers dissociation rates
	double ddh1h1;
	double ddh1h7;
	double ddh1h13;
	double ddh7h7;
	double ddh7h13;
	double ddh13h13;
	
	double critph1h1;
	double critph7h13;
	double critpd;
	
	parameters (double pa[]) {
					
		this->psh1 = pa[0];
		this->psh7 = pa[1];
		this->psh13 = pa[2];
		this->psd = pa[3];
		this->pdh1 = pa[4];
		this->pdh7 = pa[5];
		this->pdh13 = pa[6];
		this->pdd = pa[7];
		this->msh1 = pa[8];
		this->msh7 = pa[9];
		this->msh13 = pa[10];
		this->msd = pa[11];	
		this->mdh1 = pa[12];
		this->mdh7 = pa[13];
		this->mdh13 = pa[14];
		this->mdd = pa[15];
		this->ddgh1h1 = pa[16];
		this->ddgh1h7 = pa[17];
		this->ddgh1h13 = pa[18];
		this->ddgh7h7 = pa[19];
		this->ddgh7h13 = pa[20];
		this->ddgh13h13 = pa[21];
		this->delaymh1 = pa[22];
		this->delaymh7 = pa[23];
		this->delaymh13 = pa[24];
		this->delaymd = pa[25];
		this->delayph1 = pa[26];
		this->delayph7 = pa[27];
		this->delayph13 = pa[28];
		this->delaypd = pa[29];
		this->dah1h1 = pa[30];
		this->ddh1h1 = pa[31];
		this->dah1h7 = pa[32];
		this->ddh1h7 = pa[33];
		this->dah1h13 = pa[34];
		this->ddh1h13 = pa[35];
		this->dah7h7 = pa[36];
		this->ddh7h7 = pa[37];
		this->dah7h13 = pa[38];
		this->ddh7h13 = pa[39];
		this->dah13h13 = pa[40];
		this->ddh13h13 = pa[41];
		this->critph1h1 = pa[42];
		this->critph7h13 = pa[43];
		this->critpd = pa[44];
	};
	
	~parameters () {};
};

#endif

