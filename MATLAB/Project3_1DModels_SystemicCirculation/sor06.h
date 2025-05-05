/***************************************************************************/
/*                                                                         */
/*  Program: sor06.h                                                       */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This header file defines the global parameters                         */
/*                                                                         */
/***************************************************************************/

// $Id: sor06.h,v 1.8 2010-10-20 15:38:28 mette Exp $

#ifndef _SOR06_H
#define _SOR06_H

#include <cmath>
#include "tools.h"

int    nbrves;               // Number of vessels in the tree.
// MJC
int size_conn;

int    tmstps = 16384*2,   //8192              // The number of timesteps per period. Match Qin, power of 2
       plts   = 512;                 // Number of plots per period.

// for a control
const char *CO_filename = "newDORV2q.dat";
// for a diseased
// const char *CO_filename = "newHLHS7q.dat";

int max_D = 4; // MJC: For passing bifurcation pointer
        
const int WALL_MODEL = 2; // Define which wall model to use.

double conv   = 1333.220,              // Conversion from mmHg to SI-units.
       rho    = 1.057,                // Density of blood [g/cm^3].
       mu     = 0.032,                // Viscosity of blood [g/cm/s].
       mu_pl  = mu,                   // Viscosity of blood [g/cm/s].
       nu     = mu/rho,

        // for a control human
       Tper   = 0.7,        //1.1         // The period of one heart beat [s]. for exercise decrease

        // for a diseased human
       // Tper   = 0.6,        //1.1         // The period of one heart beat [s]. for exercise decrease

//       Fcst = 10.0,                 // Determines the damping coeff.
//       Fcst   = 2.0*sqrt(2.0*M_PI/nu/Tper),//DID HAVE A 1.33 MULTIPLIED HERE        // Determines the damping coeff.
                                      // for the friction.
       Lr     = 1.0,                  // Characteristic radius of the
                                      // vessels in the tree [cm].
       Lr2    = sq(Lr),               // The squared radius [cm2].
       Lr3    = cu(Lr),               // The radius to the third power [cm^3].
       g      = 981.0,                // The gravitational force [cm/s^2].
       q      = 10.0*Lr2,             // The characteristic flow [cm^3/s].
       Fr2    = sq(q)/g/pow(Lr,5),    // The squared Froudes number.
       Re     = q*rho/mu/Lr,          // Reynolds number.
       bound_thick = sqrt(nu*Tper/(2*M_PI))/Lr, // Boundary layer thickness (MJC)
       Period = Tper*q/Lr3,           // The dimension-less period.
       k      = Period/tmstps,        // Length of a timestep.
       Deltat = Period/plts,          // Interval between each point plottet.
        
        // for a control human
       p0     = 65.0*rho/g/Lr*conv,   // Ensures a certain diastolic pressure.
        
        // for a diseased human
       // p0     = 20.0*rho/g/Lr*conv,   // Ensures a certain diastolic pressure.


       *fjac[24],   // Work space used by bound_bif.
       xr, f, df;                     // Work space used by bound_right.

#endif
