// Fits.h
#ifndef __FITS_H__
#define __FITS_H__ 1
#include <vector>
#include <cmath>
#include <iostream>
#include <TMinuit.h>


////////////////////////////////////////////////////////////////////////////////
// Global variables 
////////////////////////////////////////////////////////////////////////////////
extern std::vector<double>::iterator pulse_start;
extern std::vector<double>::iterator pulse_end;
extern TMinuit *minuit;
extern double chisq;
extern double gerror;
extern double minimizer_small;


void Stats(double *stats);
void FCN1(int &npar, double *grad, double &fval, double *par, int flag);
void FCN2(int &npar, double *grad, double &fval, double *par, int flag);
void FCN3(int &npar, double *grad, double &fval, double *par, int flag);
void FCN4(int &npar, double *grad, double &fval, double *par, int flag);
void FCN5(int &npar, double *grad, double &fval, double *par, int flag);
//void FNC2(int &npar, double *grad, double &fval, double *par, int flag);
void Fit(std::vector<double> wave, double trace_start, double perror, double *results, double *results_err);
void FitFullS1(std::vector<double> *wave, size_t start, size_t end, double perror, double *results, double *results_err);


#endif
