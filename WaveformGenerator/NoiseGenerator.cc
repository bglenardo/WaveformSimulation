#include <cstdlib>
#include <ctime> 
#include <vector>
#include <iostream>
#include <cmath>
#include "TRandom3.h"
#define PI 3.14159265358979323846
/**********************************************************
* This function includes the digitized data set from the surface lab noise study PDF. 
* The data points correspond to plot B5 (Preamp + Postamp to STRUCK).
* callMyFunction.cc provides an example of how one would use this function to generate noise at some given time based on the
* noise spectrum from plot B5.
***********************************************************/

using namespace std;

double unifRand()
{
    	return (2*PI)*(rand() / double(RAND_MAX));
}

vector <double> randNoiseGeneratorNoRoot(double time_input){

	srand(time(NULL));

	int N = (sizeof(binCenter)/sizeof(*binCenter)); // number of bins of log histogram. 
	int num_samples = time_input*1e8; // time * 1e8 [samples/sec]
	double randomPhase;
	vector <double> waveform;

        double phe_scale_factor = 4.6/1000; // 4.6 mV/phe / 1000 mV/V

	for(int i=0;i<N;i++){

		randomPhase = unifRand();

		//the waveform at some given time is the sum of the noise over all the bins
		for(int s=0; s <= num_samples; s++) {
			if( i==0 ) {
				waveform.push_back(sqrt(2*binWidth[i])*binHeight[i]*sin(2*PI*binCenter[i]*s*1e-8+randomPhase) / phe_scale_factor );
			}
			else {
				waveform[s] += sqrt(2*binWidth[i])*binHeight[i]*sin(2*PI*binCenter[i]*s*1e-8+randomPhase / phe_scale_factor);
			}
		}

	}   

	return waveform;

}
