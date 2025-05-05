/***************************************************************************/
/*                                                                         */
/*  Program: sor06.h                                                       */
/*  Version: 2.0                                                           */
/*  Date: 6 July 2021                                                     */
/*                                                                         */
/*  Primary Authors: M.S. Olufsen                                          */
/*  Key Contributers: M.U. Qureshi & M.J. Colebank                         */
/*  Institute: North Carolina State University, Raleigh, NC, USA           */
/*  Funded by: NSF-DMS # 1615820 & AHA #19PRE34380459                      */
/*                                                                         */
/*  This header file defines the global parameters                         */
/*                                                                         */
/***************************************************************************/


#ifndef _SOR06_H
#define _SOR06_H

#include <cmath>
#include "tools.h"

int    nbrves;                            // Number of vessels in the tree.

int    tmstps = 8192,                     // The number of timesteps per period.
       plts   = 512;                     // Number of plots per period.

const char *CO_filename = "Qin.dat";      // Input flow file at the heart.

int max_D = 4; // MJC: For passing bifurcation pointer
     
const int VEL_MODEL = 2 ; // Define which velocity profile to use. (1: stokes, 2: power law)
double vel_power = 9;    // Power for the power-law profile (2 - poiseuille, 9 - nearly plug flow)
const int WALL_MODEL = 2;// Define which wall model to use.

double conv   = 1333.220,               // Conversion from mmHg to SI-units.
       rho    = 1.055,                  // Density of blood [g/cm^3].
       mu     = 0.049,                  // Viscosity of blood [g/cm/s].
       nu     = mu/rho,                 // Dynamic viscosity of blood [cm^2/s].
       Lr     = 1.0,                    // Characteristic radius of the
                                        // vessels in the tree [cm].
       Lr2    = sq(Lr),                 // The squared radius [cm2].
       Lr3    = cu(Lr),                 // The radius to the third power [cm^3].
       g      = 981.0,                  // The gravitational force [cm/s^2].
       q      = 10.0*Lr2,               // The characteristic flow [cm^3/s].
       Tper   = 0.90,//0.85,                   // Dimensional period [s].
       Fr2    = sq(q)/g/pow(Lr,5),      // The squared Froudes number.
       Re     = q*rho/mu/Lr,            // Reynolds number.
       Period = Tper*q/Lr3,             // The dimension-less period.
       k      = Period/tmstps,         // Length of a timestep.
       Deltat = Period/plts,           // Interval between each point plotted.
       bound_thick = sqrt(Tper*nu/(2.0*M_PI))/Lr, // Boundary layer thickness
        Fcst   = 20.0,//sqrt(2*M_PI/nu/Tper)/Lr, // FROM OLD CODE: assumes R/delta const.
       p0     = 0.0/rho/g/Lr*conv,      // Ensures a certain diastolic pressure.
       *fjac[24],                       // Work space used by bound_bif.
       xr, f, df;                       // Work space used by bound_right.

int max_cycles = 100,
    cycles     = 1;

#endif
