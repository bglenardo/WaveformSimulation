#ifndef __SMEARING_PARAMETERS_H__
#define __SMEARING_PARAMETERS_H__

// This file contains the first iteration of the smearing parameters.
// The array is organized by bin, and each bin contains:
// { low-edge-energy, mean_offset, smearing_width }


double DD_smearing[7][3] = {
{2., 1.182, 0.651},
{4., 1.148, 0.510},
{6., 1.133, 0.447},
{8., 1.123, 0.409},
{10.,1.117, 0.387},
{12.,1.112, 0.371},
{14.,1.126, 0.353}         };

double TT_smearing[7][3] = {
{2., 1.184, 0.775},
{4., 1.144, 0.576},
{6., 1.154, 0.492},
{8., 1.141, 0.453},
{10.,1.117, 0.419},
{12.,1.124, 0.398},
{14.,1.117, 0.383}         };

#endif
