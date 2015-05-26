/*
 *  randnum.cpp
 *  
 *
 *  Created by Adam Lock on 3/27/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <math.h>
#include <ctime>
#include <cstring>
#include <string>
#include <sstream>

using namespace std;

int main(){
	
	ofstream ran;
	ran.open("randnums.csv", fstream::out);
	/* Names of Random Numbers, and dXXX = delayXXX so dmh1 is really delaymh1 */
	
	/* Set the first line to Numbers to values that we know actually work */
	
	/* Create the variables for the multiple randoms */
	double psh1, psh7, psh13, psd;
	double pdh1, pdh7, pdh13, pdd;
	double msh1, msh7, msh13, msd;
	double mdh1, mdh7, mdh13, mdd;
	double ddgh1h1, ddgh1h7, ddgh1h13, ddgh7h7, ddgh7h13, ddgh13h13;
	double delaymh1, delaymh7, delaymd, delaymh13; 
	double delayph1, delayph7, delayph13, delaypd; 
	double dah1h1, ddh1h1, dah1h7, ddh1h7, dah1h13, ddh1h13, dah7h7, ddh7h7, dah7h13, ddh7h13, dah13h13, ddh13h13;
	double critph1h1, critph7h13, critpd; 
	
	/* Number of Random Numbers */
	int count = 200000;
	
	srand((unsigned)time(0));
	double psh1a = 30, psh1b = 60;
	double psh7a = 10, psh7b = 57;
	double psh13a = 27, psh13b = 57;
	double psda = 22, psdb = 59;
	double pdh1a = 0.12, pdh1b = 0.37;
	double pdh7a = 0.11, pdh7b = 0.4;
	double pdh13a = 0.11, pdh13b = 0.39;
	double pdda = 0.15, pddb = 0.38;
	double msh1a = 32, msh1b = 63;
	double msh7a = 34, msh7b = 62;
	double msh13a = 31, msh13b = 62;
	double msda = 31, msdb = 65;
	double mdh1a = 0.2, mdh1b = 0.38;
	double mdh7a = 0.28, mdh7b = 0.4;
	double mdh13a = 0.13, mdh13b = 0.39;
	double mdda = 0.12, mddb = 0.39;
	double ddgh1h1a = 0.25, ddgh1h1b = 0.4;
	double ddgh1h7a = 0.16, ddgh1h7b = 0.34;
	double ddgh1h13a = 0.1, ddgh1h13b = 0.36;
	double ddgh7h7a = 0.12, ddgh7h7b = 0.4;
	double ddgh7h13a = 0.26, ddgh7h13b = 0.4;
	double ddgh13h13a = 0.11, ddgh13h13b = 0.34;
	double delaymh1a = 8.8, delaymh1b = 12.0;
	double delaymh7a = 8.6, delaymh7b = 11.6;
	double delaymda = 6.1, delaymdb = 12.0;
	double delayph1a = 0.8, delayph1b = 2.0;
	double delayph7a = 0.4, delayph7b = 1.8;
	double delayph13a = 0.6, delayph13b = 1.8;
	double delaypda = 10, delaypdb = 18;
	double dah1h1a = 0.005, dah1h1b = 0.03;
	double ddh1h1a = 0.06, ddh1h1b = 0.3;
	double dah1h7a = 0.0006, dah1h7b = 0.009;
	double ddh1h7a = 0.03, ddh1h7b = 0.28;
	double dah1h13a = 0.006, dah1h13b = 0.029;
	double ddh1h13a = 0.004, ddh1h13b = 0.18;
	double dah7h7a = 0.002, dah7h7b = 0.024;
	double ddh7h7a = 0.07, ddh7h7b = 0.3;
	double dah7h13a = 0.007, dah7h13b = 0.03;
	double ddh7h13a = 0.03, ddh7h13b = 0.3;
	double dah13h13a = 0.001, dah13h13b = 0.016;
	double ddh13h13a = 0.05, ddh13h13b = 0.29;
	double critph1h1a = 160, critph1h1b = 720;
        double critph7h13a = 200, critph7h13b = 920;
	double critpda = 240, critpdb = 720;
	
	while(count!=0){
		// Random doubles for every object 
		psh1 = double(((psh1b - psh1a)*rand()/(RAND_MAX + 1.0))+psh1a); // 3.75-60 
		psh7 = double(((psh7b - psh7a)*rand()/(RAND_MAX + 1.0))+psh7a); // 3.75-60 
		psh13 = double(((psh13b - psh13a)*rand()/(RAND_MAX + 1.0))+psh13a); // 10-20
		psd = double(((psdb - psda)*rand()/(RAND_MAX + 1.0))+psda); // 3.75-60 
		pdh1 = double(((pdh1b - pdh1a)*rand()/(RAND_MAX + 1.0))+pdh1a); // 0.2-0.8 
		pdh7 = double(((pdh7b - pdh7a)*rand()/(RAND_MAX + 1.0))+pdh7a); // 0.2-0.8 
		pdh13 = double(((pdh13b - pdh13a)*rand()/(RAND_MAX + 1.0))+pdh13a); // 0.2-0.6 
		pdd = double(((pddb - pdda)*rand()/(RAND_MAX + 1.0))+pdda); // 0.2-0.8 
		msh1 = double(((msh1b - msh1a)*rand()/(RAND_MAX + 1.0))+msh1a); // 8.25 - 132
		msh7 = double(((msh7b - msh7a)*rand()/(RAND_MAX + 1.0))+msh7a); // 8.25 - 132
		msh13 = double(((msh13b - msh13a)*rand()/(RAND_MAX + 1.0))+msh13a);// 10-20
		msd = double(((msdb - msda)*rand()/(RAND_MAX + 1.0))+msda); // 8.25 - 132
		mdh1 = double(((mdh1b - mdh1a)*rand()/(RAND_MAX + 1.0))+mdh1a); // 0.2-0.8 
		mdh7 = double(((mdh7b - mdh7a)*rand()/(RAND_MAX + 1.0))+mdh7a); // 0.2-0.8
		mdh13 = double(((mdh13b - mdh13a)*rand()/(RAND_MAX + 1.0))+mdh13a); // 0.2-0.6
		mdd = double(((mddb - mdda)*rand()/(RAND_MAX + 1.0))+mdda); // 0.2-0.8 
		ddgh1h1 = double (((ddgh1h1b - ddgh1h1a)*rand()/(RAND_MAX + 1.0))+ddgh1h1a); //0.2-0.8
		ddgh1h7 = double (((ddgh1h7b - ddgh1h7a)*rand()/(RAND_MAX + 1.0))+ddgh1h7a); //0.2-0.8
		ddgh1h13 = double (((ddgh1h13b - ddgh1h13a)*rand()/(RAND_MAX + 1.0))+ddgh1h13a); //0.2-0.8
		ddgh7h7 = double (((ddgh7h7b - ddgh7h7a)*rand()/(RAND_MAX + 1.0))+ddgh7h7a); //0.2-0.8
		ddgh7h13 = double (((ddgh7h13b - ddgh7h13a)*rand()/(RAND_MAX + 1.0))+ddgh7h13a); //0.2-0.8
		ddgh13h13 = double (((ddgh1h1b - ddgh13h13a)*rand()/(RAND_MAX + 1.0))+ddgh13h13a); //0.2-0.8
		delaymh1 = double(((delaymh1b - delaymh1a)*rand()/(RAND_MAX + 1.0))+delaymh1a); // 3-12 
		delaymh7 = double(((delaymh7b - delaymh7a)*rand()/(RAND_MAX + 1.0))+delaymh7a); // 3-12  
		delaymd = double(((delaymdb - delaymda)*rand()/(RAND_MAX + 1.0))+delaymda); // 3-12 
		delaymh13 = -1;
		delayph1 = double(((delayph1b - delayph1a)*rand()/(RAND_MAX + 1.0))+delayph1a); // 0.3-2.4 
		delayph7 = double(((delayph7b - delayph7a)*rand()/(RAND_MAX + 1.0))+delayph7a); // 0.3-2.4 
		delayph13 = double(((delayph13b - delayph13a)*rand()/(RAND_MAX + 1.0))+delayph13a); // 0.3-2.4 
		delaypd = double(((delaypdb - delaypda)*rand()/(RAND_MAX + 1.0))+delaypda); // 9-36 
		dah1h1 = double(((dah1h1b - dah1h1a)*rand()/(RAND_MAX + 1.0))+dah1h1a); // 0.00003-0.003 
		ddh1h1 = double(((ddh1h1b - ddh1h1a)*rand()/(RAND_MAX + 1.0))+ddh1h1a); // 0.003-0.3 
		dah1h7 = double(((dah1h7b - dah1h7a)*rand()/(RAND_MAX + 1.0))+dah1h7a); // 0.00003-0.003 
		ddh1h7 = double(((ddh1h7b - ddh1h7a)*rand()/(RAND_MAX + 1.0))+ddh1h7a); // 0.003-0.3 
		dah1h13 = double(((dah1h13b - dah1h13a)*rand()/(RAND_MAX + 1.0))+dah1h13a); // 0.00003-0.003 
		ddh1h13 = double(((ddh1h13b - ddh1h13a)*rand()/(RAND_MAX + 1.0))+ddh1h13a); // 0.003-0.3 
		dah7h7 = double(((dah7h7b - dah7h7a)*rand()/(RAND_MAX + 1.0))+dah7h7a); // 0.00003-0.003 
		ddh7h7 = double(((ddh7h7b - ddh7h7a)*rand()/(RAND_MAX + 1.0))+ddh7h7a); // 0.003-0.3 
		dah7h13 = double(((dah7h13b - dah7h13a)*rand()/(RAND_MAX + 1.0))+dah7h13a); // 0.00003-0.003 
		ddh7h13 = double(((ddh7h13b - ddh7h13a)*rand()/(RAND_MAX + 1.0))+ddh7h13a); // 0.003-0.3 
		dah13h13 = double(((dah13h13b - dah13h13a)*rand()/(RAND_MAX + 1.0))+dah13h13a); // 0.00003-0.003 
		ddh13h13 = double(((ddh13h13b - ddh13h13a)*rand()/(RAND_MAX + 1.0))+ddh13h13a); // 0.003-0.3 
		critph1h1 = double(((critph1h1b - critph1h1a)*rand()/(RAND_MAX + 1.0))+critph1h1a); // 30-480 
        critph7h13 = double(((critph7h13b - critph7h13a)*rand()/(RAND_MAX + 1.0))+critph7h13a); // 30-480 
		critpd = double(((critpdb - critpda)*rand()/(RAND_MAX + 1.0))+critpda); // 150-2400 
		
		/* Place them all in the file separated by commas */
		ran<<psh1<<","<<psh7<<","<<psh13<<","<<psd<<","<<pdh1<<","<<pdh7<<","<<pdh13<<","<<pdd<<","<<msh1<<","<<msh7<<","<<msh13<<","<<msd<<",";
		ran<<mdh1<<","<<mdh7<<","<<mdh13<<","<<mdd<<","<<ddgh1h1<<","<<ddgh1h7<<","<<ddgh1h13<<","<<ddgh7h7<<","<<ddgh7h13<<","<<ddgh13h13<<",";
		ran<<delaymh1<<","<<delaymh7<<","<<delaymh13<<","<<delaymd<<","<<delayph1<<","<<delayph7<<","<<delayph13<<","<<delaypd<<","<<dah1h1<<",";
		ran<<ddh1h1<<","<<dah1h7<<","<<ddh1h7<<","<<dah1h13<<","<<ddh1h13<<","<<dah7h7<<","<<ddh7h7<<","<<dah7h13<<","<<ddh7h13<<","<<dah13h13<<","<<ddh13h13<<","<<critph1h1 << "," << critph7h13 <<","<<critpd<<endl;
		
		count--;
		
		if (count % 100000 == 0)
			cout << count << " numbers left" << endl;
	}
	ran.close(); /* Close the file when done */
}
