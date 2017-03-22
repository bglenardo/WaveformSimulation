#ifndef __SMEARING_PARAMETERS_H__
#define __SMEARING_PARAMETERS_H__

// This file contains the first iteration of the smearing parameters.
// // The array is organized by bin, and each bin contains:
// // { low-edge-energy, mean_offset, smearing_width }


 double DD_smearing[7][3] = { 
 {2., 1.070118, 0.630031},
 {4., 1.044759, 0.511007},
 {6., 1.037924, 0.439378},
 {8., 1.033532, 0.407538},
 {10.,1.026731, 0.385569},
 {12.,1.013810, 0.369530},
 {14.,1.036969, 0.348577}         };  

 double TT_smearing[7][3] = { 
 {2., 1.087896, 0.743315},
 {4., 1.052010, 0.571158},
 {6., 1.046328, 0.478643},
 {8., 1.040452, 0.441292},
 {10.,1.031122, 0.408190},
 {12.,1.018680, 0.397237},
 {14.,1.030091, 0.379531}         };  

 #endif
