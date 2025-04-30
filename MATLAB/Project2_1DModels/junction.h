
/***************************************************************************/
/*                                                                         */
/* The junction.h header program                                            */
/*  Version: 1.0                                                           */
/*  Date: 24 March 2020                                                    */
/*                                                                         */
/*  Primary Authors: M.J. Colebank & J. Mackenzie                          */
/*  Key Contributers: M.U. Qureshi & M.S. Olufsen                          */
/*  Institute: North Carolina State University, Raleigh, NC, USA           */
/*  Funded by: NSF-DMS # 1615820 & AHA #19PRE34380459                      */
/*                                                                         */
/***************************************************************************/


#ifndef _JUNCTION_H
#define _JUNCTION_H

#include <cmath>
#include "tools.h"

extern double vel_power;
extern double Fcst;

// Updates bifurcation conditions. Uses daughter vessels, and should
// only be called when such exist.

void junction   (double theta, double gamma, Tube *Arteries[], int parent, int WM, int VP);
void bound_monf (double theta, double gamma, Tube *Arteries[], int parent, int D1, int WM, int VP);
void bound_screen_lnsearch (double theta, double gamma, Tube *Arteries[], int parent, int D1, int WM, int VP);
void bound_sten (double theta, double gamma, Tube *Arteries[], int parent, int D1, int WM, int VP);
void bound_bif  (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int WM, int VP);
void bound_sten_bif  (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int WM, int VP);
void bound_sten_daughter  (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int D3, int WM, int VP);

void bound_trif (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int D3, int WM, int VP);


//double *fjac[24]; // Work space used by bound_junc.


// In order to ensure a more efficient execution of the program the following
// functions is made as in-line functions.

// A function returning the Friction of the system. The definition of this
// function is given according to the derivation in the mathematical model.
// The constant cst, determines the amount of damping in the system.
inline double F (double Q, double A, int VP)
  {
      double tmp1,tmp2,tmp3;
      
      if (VP==1) { // Stokes boundary layer
          tmp1 = -2.0*sqrt(M_PI)*Q; //MJC
          tmp2 = bound_thick*Re*sqrt(A);         //MJC
          tmp3 = tmp1/tmp2;
      }
      else if (VP==-1) { // Stokes boundary layer, fixed
          tmp1 = (-2.0*Fcst*M_PI*Q); //MJC
          tmp2 = (Re*A);         //MJC
          tmp3 = tmp1/tmp2;
      }
  else
  { // Power law
        tmp1 = -2.0*M_PI*(vel_power+2.0)*Q; //MJC
        tmp2 = Re*A;         //MJC
        tmp3 = tmp1/tmp2;
  }
    return(tmp3);
  }

inline double dFdQ (double A, int VP)
{
    double tmp1,tmp2,tmp3;
    if (VP==1) { // Stokes boundary layer
            tmp1 = (-2.0*sqrt(M_PI)); //MJC
            tmp2 = (bound_thick*Re*sqrt(A));         //MJC
            tmp3 = tmp1/tmp2;
        }
    else if (VP==-1) { // Stokes boundary layer, fixed
        tmp1 = (-2.0*Fcst*M_PI); //MJC
        tmp2 = (Re*sqrt(A));         //MJC
        tmp3 = tmp1/tmp2;
    }
    else
    { // Power law
          tmp1 = -2.0*M_PI*(vel_power+2.0); //MJC
          tmp2 = Re*A;         //MJC
          tmp3 = tmp1/tmp2;
    }
      return(tmp3);
}
    
inline double dFdA (double Q, double A, int VP)
{
    double tmp1,tmp2,tmp3;
    if (VP==1) { // Stokes boundary layer
            tmp1 = Q*sqrt(M_PI); //MJC
            tmp2 = bound_thick*Re*sqrt(cu(A));         //MJC
            tmp3 = tmp1/tmp2;
        }
    else if (VP==-1) { // Stokes boundary layer, fixed
        tmp1 = (2.0*Fcst*M_PI*Q); //MJC
        tmp2 = Re*sq(A);         //MJC
        tmp3 = tmp1/tmp2;
        
        tmp1 = (-2.0*M_PI*Q); //MJC
        tmp2 = (bound_thick*Re*A);         //MJC
        tmp3 = tmp1/tmp2;
    }
    else
    { // Power law
          tmp1 = 4.0*M_PI*(vel_power+2.0)*Q; //MJC
          tmp2 = Re*sq(A);         //MJC
          tmp3 = tmp1/tmp2;
    }
    return(tmp3); //MJC
}

#endif
