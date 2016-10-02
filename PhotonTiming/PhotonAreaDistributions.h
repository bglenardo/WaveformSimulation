#ifndef _PHOTONAREADISTRIBUTIONS_H_
#define _PHOTONAREADISTRIBUTIONS_H_

double Baseline_area_PDF(double peak_area_phe);
double SPE_area_PDF(double peak_area_phe);
double DPE_area_PDF(double peak_area_phe);
double TPE_area_PDF(double peak_area_phe);
double QPE_area_PDF(double peak_area_phe);
double QIPE_area_PDF(double peak_area_phe);
double SPE_area_CDF_int(double a, double b);

#endif
