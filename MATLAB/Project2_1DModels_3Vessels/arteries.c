/***************************************************************************/
/*                                                                         */
/*  Program: arteries.C                                                    */
/*  Version: 1.0                                                           */
/*  Date: 30 Dec. 2019                                                     */
/*                                                                         */
/*  Primary Authors: M.S. Olufsen                                          */
/*  Key Contributers: M.U. Qureshi & M.J. Colebank                         */
/*  Institute: North Carolina State University, Raleigh, NC, USA           */
/*  Funded by: NSF-DMS # 1615820 & AHA #19PRE34380459                      */
/*                                                                         */
/*  This module can predict the flow and pressure in an tree of elastic    */
/*  vessels as described in REFERENCE.pdf. The dependencies of the vessels */
/*  in the tree must be specified in the main module according to the tree */
/*  in question.                                                           */
/*                                                                         */
/*  This module includes all the functions needed to solve the system      */
/*  of equations. That is the description of all functions in the class    */
/*  containing the vessel (for further details see arteries.h), and in     */
/*  particular the functions needed to solve the system of equations nu-   */
/*  merically.                                                             */
/*                                                                         */
/*  The module is dependent on the utilities in tools.C, and               */
/*  their corresponding h-files, and also arteries.h that includes the     */
/*  declaration of the vessel-object.                                      */
/*                                                                         */
/***************************************************************************/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>

#include "tools.h"
#include "arteries.h"
#include "junction.h"

using namespace std;


extern int nbrves;
extern int tmstps;
extern char* CO_filename;

/* Methods of class Tube, see arteries.h for description of this. */

// The constructor. When an object is made this function will initialize
// all the attributes of the specific tube. The parameters for the length
// of the specific vessel, the top and bottom radii, and if applicable
// the pointers to the daughter arteries will be initialized according
// to the actual parameters passed in the call.
// If the tube is terminal then the peripheral resistance must be set,
// and the daughter vessels should be NIL. Otherwise the pointers to
// the daughter vessels must be given.
// Further all the work arrays are declared and initialized, and the
// initial condition for the system equations is applied.
Tube :: Tube (double Length,
              double topradius, double botradius,
              int *daughters,
              double points, int init,
              double f1, double f2, double f3, double Res1_in, double Res2_in,
              double CT_in, double sten_factor, double sten_length):
    L(Length),
    rtop(topradius),
    rbot(botradius),
    daughters(daughters),
    pts(points),
    ff1(f1),
    ff2(f2),
    ff3(f3),
    Res1(Res1_in),
    Res2(Res2_in),
    CT(CT_in),
    sten_factor(sten_factor),
    sten_length(sten_length)
{
  // Initialization of the basic parameters
  N      = int(pts*L);
  h      = 1.0/pts/Lr;

  // Declaration and Initialization of the needed intermediate arrays.
  Qnew      = new double[N+1];
  Anew      = new double[N+1];
  Qold      = new double[N+1];
  Aold      = new double[N+1];
  Qprv      = new double[N+1];
  Aprv      = new double[N+1];
  R1      = new double[N+1];
  R2      = new double[N+1];
  S1      = new double[N+1];
  S2      = new double[N+1];
  r0      = new double[N+1];
  r0h      = new double[N+2];
  dr0dx   = new double[N+1];
  dr0dxh  = new double[N+2];
  A0      = new double[N+1];
  A0h     = new double[N+2];
  fr      = new double[N+1];
  frh     = new double[N+2];
  dfrdr0  = new double[N+1];
  dfrdr0h = new double[N+2];
  Ah      = new double[N];
  Qh      = new double[N];
  R1h      = new double[N];
  R2h      = new double[N];
  S1h      = new double[N];
  S2h      = new double[N];

  // MJC: Define nondimensional stiffness ahead of time
  double f1ND = ff1/rho/g/Lr;
  double f2ND = ff2*Lr;
  double f3ND = ff3/rho/g/Lr;

  // Vessel geometry is tabulated and initial conditions are applied
  for (int i=0; i<=N; i++)
  {
    r0 [i]     = rtop*exp(i*log(rbot/rtop)/N)/Lr; // ND
    r0h[i]     = rtop*exp((i-0.5)*log(rbot/rtop)/N)/Lr; // ND
    dr0dx [i]  = log(rbot/rtop)/h/N*r0 [i]; // ND
    dr0dxh[i]  = log(rbot/rtop)/h/N*r0h[i]; // ND
    A0 [i]     = M_PI*sq(r0 [i]); // ND
    A0h[i]     = M_PI*sq(r0h[i]); // ND
    fr [i]     = (4.0/3.0)*(f1ND*exp(f2ND*r0 [i])+f3ND); // ND
    frh[i]     = (4.0/3.0)*(f1ND*exp(f2ND*r0h[i])+f3ND); // ND
    dfrdr0 [i] = (4.0/3.0)*f1ND*f2ND*exp(f2ND*r0 [i]); // ND
    dfrdr0h[i] = (4.0/3.0)*f1ND*f2ND*exp(f2ND*r0h[i]); // ND
    Qnew[i]    = 1e-8;//0.0;
    Anew[i]    = A0[i]; // ND
  }
  r0h[N+1]     = rtop*exp((N+0.5)*log(rbot/rtop)/N)/Lr; // ND
  dr0dxh[N+1]  = log(rbot/rtop)/h/N*r0h[N+1]; // ND
  A0h[N+1]     = M_PI*sq(r0h[N+1]); // ND
  frh[N+1]     = (4.0/3.0)*(f1ND*exp(f2ND*r0h[N+1])+f3ND); // ND
  dfrdr0h[N+1] = (4.0/3.0)*f1ND*f2ND*exp(f2ND*r0h[N+1]); // ND

  // Read from file data for the inflow profile.
  if (init == 1)
  {
      
    Q0 = new double[tmstps+1];

    FILE *fi = fopen (CO_filename, "r");

    for (int i=0; i<=tmstps; i++)
    {
      fscanf(fi,"%lf",&Q0[i]); // ND
      Q0[i] = Q0[i]/q; // If the indata have dimensions they should be made
                       // non-dimensional.
    }
  }

}

// The destructor. When the tube-objects terminates, all arrays are deleted,
// in order to free the memory occupied by the object.
Tube :: ~Tube ()
{
  delete[] Anew;
  delete[] Qnew;
  delete[] Aold;
  delete[] Qold;
  delete[] Aprv;
  delete[] Qprv;
  delete[] Ah;
  delete[] Qh;
  delete[] y;
  delete[] pL;
  delete[] R1h;
  delete[] R2h;
  delete[] S1h;
  delete[] S2h;
  delete[] R1;
  delete[] R2;
  delete[] S1;
  delete[] S2;
  delete[] r0;
  delete[] r0h;
  delete[] dr0dx;
  delete[] dr0dxh;
  delete[] A0;
  delete[] A0h;
  delete[] fr;
  delete[] frh;
  delete[] dfrdr0;
  delete[] dfrdr0h;
}

// ----------------------PLOTTING ROUTINES WITH DIMENSIONS ------------

void Tube :: printQ0 (FILE *fd)
{
  for (int i=0; i<=tmstps; i++)
  {
    fprintf (fd, "%15.10f\n", Q0[i]*q);
  }
}

// The following functions prints p, q(x,t) in terms of the re-dimensionalized
// variables. The parameters for the function are the  position (x),
// and the time (t).
void Tube :: printPt (FILE *fd, double t, int i, int WM)
{
  fprintf (fd, "%13.10f %15.10f\n", t*Lr3/q, (P(i,Anew[i],WM)+p0)*rho*g*Lr/conv);
}

void Tube :: printQt (FILE *fd, double t, int i)
{
  fprintf (fd, "%13.10f %15.10f\n", t*Lr3/q, Qnew[i]*q);
}

void Tube :: printAt (FILE *fd, double t, int i)
{
  fprintf (fd, "%13.10f %15.10f\n", t*Lr3/q, Anew[i]*Lr2);
}

void Tube :: printFt (FILE *fd, double t, int i, int VP)
{
  fprintf (fd, "%13.10f %15.10f\n", t*Lr3/q, (sq(q)/Lr3)*F(Qnew[i],Anew[i],VP));
    //.(Re*mu*q/rho/Lr2)*
}

void Tube :: printAllt (FILE *fd,double t, int i, int WM)
{
    fflush(stdout);
    fprintf (fd, "%13.10f %13.10f %15.10f %15.10f %15.10f %17.10f\n", t*Lr3/q, i*h*Lr, (P(i,Anew[i],WM)+p0)*rho*g*Lr/conv, Qnew[i]*q, Anew[i]*Lr2, c(i, Anew[i],WM));// NEED TO REDIMENSIONALIZE WAVESPEED
}

// The following functions prints P, Q, A, and F as functions of
// (x, t). This is done in terms of the re-dimensionalized variables.
// In this case the functions is plotted for a
// fixed time, but for all x along the vessel in question. Since the
// doesn't have to be the first vessel in the tree, it would have
// some offset from the heart. Which determines the position for x.
// Therefore there are two arguments passed to this function the time
// and the offset.
void Tube :: printPxt (FILE *fd, double t, int offset, int WM)
{
  if (offset == 0) fprintf (fd, "\n");
  for (int i=0; i<=N; i++)
  {
      fprintf (fd, "%13.10f %13.10f %15.10f %15.10f %15.10f %17.10f\n",
               t*Lr3/q, (i+offset)*h*Lr, (P(i,Anew[i],WM)+p0)*rho*g*Lr/conv, Qnew[i]*q, Anew[i]*Lr2, c(i, Anew[i],WM)); // modifeid to print more variables
  }
}

void Tube :: printQxt (FILE *fd, double t, int offset)
{
  if (offset == 0) fprintf (fd, "\n");
  for (int i=0; i<N; i++)
  {
    fprintf (fd, "%13.10f %13.10f %15.10f\n",
             t*Lr3/q, (i+offset)*h*Lr, Qnew[i]*q);
  }
}

void Tube :: printAxt (FILE *fd, double t, int offset)
{
  for (int i=0; i<=N; i++)
  {
    fprintf (fd, "%13.10f %13.10f %15.10f\n",
             t*Lr3/q, (i+offset)*h*Lr, Anew[i]*Lr2);
  }
}

void Tube :: printFxt (FILE *fd, double t, int offset, int VP)
{
  for (int i=0; i<=N; i++)
  {
    fprintf (fd, "%13.10f %13.10f %15.10f\n",
             t*Lr3/q, (i+offset)*h*Lr, F(Qnew[i],Anew[i],VP)*sq(q)/Lr3);
  }
}

// A function that prints p(Q) for all t. This is done in terms of
// the re-dimensionalized variables. In this case the plot is made for a
// fixed point in space, but for all t along the vessel in question.
void Tube :: printPQ (FILE *fd, int i, int WM)
{
  fprintf (fd,"%15.10f %15.10f\n",
           Qnew[i]*q, (P(i,Anew[i],WM)+p0)*rho*g*Lr/conv);
}

void Tube :: printPA (FILE *fd, int i, int WM)
{
  fprintf (fd,"%15.10f %15.10f\n",
           (P(i,Anew[i],WM)+p0)*rho*g*Lr/conv, Anew[i]*Lr2);
}

// Plotting the terms in the continuity equation on dimension-less form.
void Tube :: printdQdx (FILE *fd, double t, int i)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           (Qnew[i+1]-Qnew[i-1])/2.0/h);

}

void Tube :: printdAdt (FILE *fd, double t, int i, double Aprev, double tmst)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           (Anew[i]-Aprev)/tmst);
}

void Tube :: printTotConRes (FILE *fd, double t, int i, double Aprev, double tmst)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
          (Qnew[i+1] - Qnew[i-1])/2.0/h +
          (Anew[i]-Aprev)/tmst);
}

// Plotting the terms in the momentum equation on dimension-less form.
void Tube :: printdQdt (FILE *fd, double t, int i, double Qprev, double tmst)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           (Qnew[i]-Qprev)/tmst);
}

void Tube :: printddxQ2divA (FILE *fd, double t, int i)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           (sq(Qnew[i+1])/Anew[i+1] - sq(Qnew[i-1])/Anew[i-1])/2.0/h);

}

void Tube :: printdPdx (FILE *fd, double t, int i, int WM)
{
    fprintf (fd, "%13.10f %15.10f\n", t,
             Anew[i]*(P(i+1,Anew[i+1],WM)-P(i-1,Anew[i-1],WM))/Fr2/2.0/h);
}

void Tube :: printFric (FILE *fd, double t, int i, int VP)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           F(Qnew[i],Anew[i],VP));
}

void Tube :: printTotMomRes (FILE *fd, double t, int i, double Qprev, double tmst, int WM, int VP)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           Anew[i]*(P(i+1,Anew[i+1],WM)-P(i-1,Anew[i-1],WM))/Fr2/2.0/h +
           (sq(Qnew[i+1])/Anew[i+1] - sq(Qnew[i-1])/Anew[i-1])/2.0/h +
           (Qnew[i]-Qprev)/tmst + F(Qnew[i],Anew[i],VP));
}

/* Mathematical definitions*/

// Pressure area relation
double Tube :: P (int i, double A, int WM)
{
      double pold;
      if (WM==1){
          pold = fr[i]*(1.0-sqrt(A0[i]/A));
      }
      else {
          pold = fr[i]*(sqrt(A/A0[i])-1.0);
      }
   
    return pold;
}

// Derivative of the tube-law with respect to area
double Tube :: dPdA (int i, double A, int WM)
{
    double pold;
      if (WM==1){
          pold = 0.5*fr[i]*sqrt(A0[i]/cu(A));
      }
      else {
          pold = 0.5*fr[i]/sqrt(A0[i]*A);
      }
    return pold;
}

double Tube :: dPdx1(int i, double A, int WM)
{
 
    double pold;
    if (WM==1){
        pold = (dfrdr0[i]*(1.0-sqrt(A0[i]/A))-fr[i]*sqrt(M_PI/A))*dr0dx[i];
    }
    else {
        pold = (dfrdr0[i]*(sqrt(A/A0[i])-1.0)-fr[i]*sqrt(M_PI*A)/A0[i])*dr0dx[i];
    }
    return pold;
}

double Tube :: B (int i, double A, int WM)
{
    double pold;
    if (WM==1){
       pold = fr[i]*(sqrt(A0[i]*A)-A0[i])/Fr2;
    }
    else {
       pold = fr[i]*(sqrt(cu(A)/A0[i])-A0[i])/Fr2/3.0;  // Linear B
    }
    return pold;
}

double Tube :: Bh (int i, double A, int WM)
{
   int ip1 = i+1;
    double pold;
    if (WM==1){
        pold = frh[ip1]*(sqrt(A0h[ip1]*A)-A0h[ip1])/Fr2;
    }
    else {
        pold = frh[ip1]*(sqrt(cu(A)/A0h[ip1])-A0h[ip1])/Fr2/3.0;  // Linear B
    }
   return pold;

}

double Tube :: B_RHS (int i, double A, int WM)
{
    double dfr = dfrdr0[i];
    double pold, inteval;
    if (WM==1){
           pold = dr0dx[i]/Fr2*(
                  sqrt(A)*(2.0*sqrt(M_PI)*fr[i]+2.0*sqrt(A0[i])*dfr)-
                  A*dfr-2.0*M_PI*r0[i]*fr[i]-A0[i]*dfr);
        
        }
       else {

           inteval = 2.0*M_PI*r0[i]*fr[i]+A0[i]*dfr;
           pold = dr0dx[i]*(2.0*sqrt(cu(A))*(sqrt(M_PI)*fr[i]-sqrt(A0[i])*dfr)/A0[i]/3.0
                                   + A*dfr - inteval/3.0)/Fr2; // Linear B
           
       }

    return pold;
}

double Tube :: B_RHSh (int i, double A, int WM)
{
    int ip1 = i+1;

      double dfr = dfrdr0h[ip1];
      double pold,inteval;
    
      if (WM==1){
          pold = dr0dxh[ip1]/Fr2*(
               sqrt(A)*(2.0*sqrt(M_PI)*frh[ip1]+2.0*sqrt(A0h[ip1])*dfr)-
               A*dfr-2.0*M_PI*r0h[ip1]*frh[ip1]-A0h[ip1]*dfr);
       }
      else {
          inteval = 2.0*M_PI*r0h[ip1]*frh[ip1]+A0h[ip1]*dfr;
          pold = dr0dxh[ip1]*(2.0*sqrt(cu(A))*(sqrt(M_PI)*frh[ip1]-sqrt(A0h[ip1])*dfr)/A0h[ip1]/3.0
                                     + A*dfr - inteval/3.0)/Fr2; // Linear B
      }

  return pold;
}

double Tube :: dBdAh (int i, double A, int WM)
{
  int ip1      = i+1;
    double pold;
    if (WM==1){
          pold = 0.5*frh[ip1]*sqrt(A0h[ip1]/A)/Fr2;
     }
    else {
          pold = 0.5*frh[ip1]*sqrt(A/A0h[ip1])/Fr2;  // Linear B
    }
  return pold;
}

double Tube :: dB_RHSdAh (int i, double A, int WM)
{
   int ip1 = i+1;

    double dfr = dfrdr0h[ip1];
     double pold;
   if (WM==1){
         pold = (-dfr+1.0/sqrt(A)*(sqrt(M_PI)*frh[ip1]+
                sqrt(A0h[ip1])*dfr))*dr0dxh[ip1]/Fr2;
    }
   else {
           pold = (dfr+sqrt(A)*(sqrt(M_PI)*frh[ip1] - sqrt(A0h[ip1])*dfr)/A0h[ip1])*dr0dxh[ip1]/Fr2;  // Linear B
   }
    return pold;
}

// When determining or checking the step-size (k) the CFL-condition is applied.
// This is determined according to the result reached from the analysis
// made using the method of characteristics (See IMFUFATEKST no 297).
// In this function the minimal step-size fulfilling this condition for this
// tube is returned.
double Tube :: CFL (int WM) // The CFL-condition
{
  double minimum = 64000000.0;
  for (int i=0; i<=N; i++)
  {
    double c_tmp = c(i, Anew[i],WM);
    double Vnew  = Qnew[i]/Anew[i];
    double temp = min (h / fabs (Vnew - c_tmp),
                h / fabs (Vnew + c_tmp));
    if (temp < minimum) minimum = temp;
  }
  return (minimum);
}

// When taking a Lax-Wendroff step, the flux of the system must be determined.
// This is evaluated at i + j/2, and the prediction is given as described
// in IMFUFATEKST no 297 and D2.1-4. The integer k determines whether we deal
// with the first or the second component of the vector.
double Tube :: Rvec (int k, int i, int j, double Q, double A, int WM, int VP)
{
    double qA2_term=0.0;
    if (VP==2) {
        qA2_term=((vel_power+2.0)/(vel_power+1.0)); // Power law
    } else if (VP<2) {

        qA2_term = 1.0; // Stokes boundary layer
    }
    else error ("arteries.cxx","Velocity profile doesn't exist.");
        
    
    if(k==1) {
        return(Q);
    }
    else if(k==2) {
        return(qA2_term*sq(Q)/A + ((j==0)?B(i,A,WM):Bh(i,A,WM)));
    }
  else error ("arteries.cxx","Call of non-existing vector-component of R");
  return(0);
}

// Similarly the right hand side of the system of equations must be determined
// at i + j/2. Also in this case the function is given as stated in
// the mathematical model, and also in this case k states the needed component
// of the vector.
double Tube :: Svec (int k, int i, int j, double Q, double A, int WM, int VP)
{
  if(k==1) return(0.0); else
  if(k==2) return(F(Q,A,VP) + ((j==0)?B_RHS(i,A,WM):B_RHSh(i,A,WM)));
  else error ("arteries.cxx","Call of non-existing vector-component of S");
  return(0);
}

// The solutions of Anew and Qnew are found for all interior points
// of the vessel at (t+k), where k is the length of the current
// time-step. This function saves the results in the arrays Anew and
// Qnew, and the function is made according to Lax-Wendroff's method
// as described in Olufsen, et al., Ann Biomed Eng 28, 1281-1299, 2000.


void Tube :: step (double k, int WM, int VP)
{
  double theta = k/h;    // Theta is determined.
  double gamma = 0.5*k;  // Gamma is determined.

  for (int i=0; i<=N; i++)  // Remember the values at this time level.
  {
    Qold[i] = Qnew[i];
    Aold[i] = Anew[i];
  }

  // Anew and Qnew are predicted at the new time level (t+k).
  for (int i=0; i<=N; i++)
  {
    R1[i] = Rvec(1,i,0,Qold[i],Aold[i],WM,VP);
    R2[i] = Rvec(2,i,0,Qold[i],Aold[i],WM,VP);
    S1[i] = Svec(1,i,0,Qold[i],Aold[i],WM,VP);
    S2[i] = Svec(2,i,0,Qold[i],Aold[i],WM,VP);
  }

  for (int i=0; i<N; i++)
  {
    Ah[i]  = 0.5*(Aold[i+1]+Aold[i]) - 0.5*theta*(R1[i+1]-R1[i]) +
         0.5*gamma*(S1[i+1]+S1[i]);
    Qh[i]  = 0.5*(Qold[i+1]+Qold[i]) - 0.5*theta*(R2[i+1]-R2[i]) +
         0.5*gamma*(S2[i+1]+S2[i]);
    R1h[i] = Rvec(1,i,1,Qh[i],Ah[i],WM,VP);
    R2h[i] = Rvec(2,i,1,Qh[i],Ah[i],WM,VP);
    S1h[i] = Svec(1,i,1,Qh[i],Ah[i],WM,VP);
    S2h[i] = Svec(2,i,1,Qh[i],Ah[i],WM,VP);
  }
  for (int i=1; i<N; i++)
  {
    Anew[i] = Aold[i] - theta*(R1h[i]-R1h[i-1]) + gamma*(S1h[i]+S1h[i-1]);
    Qnew[i] = Qold[i] - theta*(R2h[i]-R2h[i-1]) + gamma*(S2h[i]+S2h[i-1]);
  }
}

// The left boundary (x=0) uses this function to model an inflow into
// the system. The actual parameter given to the function is the model time.
// As stated in the mathematical model the constants of the function are
// chosen in order to ensure a certain CO (specified in main.h). Hence we have
// the specified value of b. Further the period (dimension-less) is assumed
// to be Period.
double Tube :: Q0_init (double t, double k, double Period)
{
  if (t <= Period) return (Q0[int(t/k)]); else
  if (t >  Period) return (Q0_init((t-Period),k,Period));
  else return (0);
}


// Update of the left boundary at time t. This function uses Q0 to determine
// the flow rate at the next time-step. From this the value of A is predicted
// using Lax-Wendroff's numerical scheme. This function is only relevant
// when the tube is an inlet vessel.
void Tube :: bound_left (double t, double k, double Period, int WM, int VP)
{
    
    Qnew[0]   = Q0_init(t,k,Period);
    
    if (int(t/k) < 0)
        
        printf("t/k negative in bound_left\n");
    
    // Use Negative characteristic
    double qS, aS, cS, HnS, uS;
    qS = aS = cS = HnS = 0.0;
    negchar(k/h, qS, aS, cS, HnS, WM, VP);
    uS = qS/aS;
    Anew[0]   = aS + (Qnew[0] - qS)/(uS + cS) + k*HnS;
    
    // Or use Lax-Wendroff
//    double Qhm05 = Qnew[0]+Q0_init(t-k,k,Period) - Qh[0];
//    double R1hm05    = Qhm05;
//    Anew[0]   = Aold[0] - k*(R1h[0] - R1hm05)/h;
//    fprintf(stdout,"Ah/Qh in step: %lf %lf",Ah[i],Qh[i]);
}

// The value at the right boundary at time t is predicted. NB: This should
// only be used with terminal vessels, i.e. for vessels that don't bifurcate
// into further branches.
// In that situation the bifurcation boundary function should be called
// instead. Again the procedure specified is given according to the mathemati-
// cal theory presented in Olufsen, et al., Ann Biomed Eng 28, 1281?1299, 2000.

double Tube :: c (int i, double A, int WM) // The wave speed through aorta.
{
    double cnst;
    // NOTE: Need to come back and fix this; the wave speed
    // only needs to be nondimensionalized by sqrt(g*Lr)
      if (WM==1){
          cnst =  0.5*fr[i]*sqrt(A0[i]/A)/Fr2;
      }
      else {
          cnst =  0.5*fr[i]*sqrt(A/A0[i])/Fr2;  //Linear B
      }
    return sqrt (cnst);
}

double Tube :: Hp (int i, double Q, double A, int WM, int VP)
{
  return (F(Q,A,VP) - A*dPdx1(i,A,WM))/(-Q/A + c(i,A,WM));
}

double Tube :: Hn (int i, double Q, double A, int WM, int VP)
{
    return (F(Q,A,VP) - A*dPdx1(i,A,WM))/(-Q/A - c(i,A,WM));
}

void Tube :: poschar (double theta, double &qR, double &aR, double &cR, double &HpR, int WM, int VP)
{
  double ctm1  = c  (N, Aold[N],WM);
  double Hptm1 = Hp (N, Qold[N], Aold[N],WM,VP);
  double uR    = Qold[N] / Aold[N];
  double ch    = (uR + ctm1) * theta;

  if (uR + ctm1 < 0)
  {
    printf("uR + ctm1 < 0, CFL condition violated\n");
      exit(1);
  }

  qR  = Qold[N] - (Qold[N] - Qold[N-1])*ch;
  aR  = Aold[N] - (Aold[N] - Aold[N-1])*ch;
  cR  = ctm1    - (ctm1  - c (N-1,Aold[N-1],WM))*ch;
  HpR = Hptm1   - (Hptm1 - Hp(N-1,Qold[N-1],Aold[N-1],WM,VP))*ch;
}

void Tube :: negchar (double theta, double &qS, double &aS, double &cS, double &HnS, int WM, int VP)
{
    double ctm1  = c(0, Aold[0],WM);
    double Hntm1 = Hn(0, Qold[0], Aold[0],WM,VP);
    double uS    = Qold[0]/Aold[0];
    double ch    = (uS - ctm1) * theta;
    
    if ( ctm1 - uS < 0)
    {
        printf("ctm1 - uS < 0, CFL condition violated\n");
        exit(1);
    }
    
    qS  = Qold[0] + (Qold[0] - Qold[1])*ch;
    aS  = Aold[0] + (Aold[0] - Aold[1])*ch;
    cS  = ctm1    + (ctm1  - c (1,Aold[1],WM))*ch;
    HnS = Hntm1   + (Hntm1 - Hn(1,Qold[1],Aold[1],WM,VP))*ch;
}


// FROM UMAR
void Tube :: bound_right (double k, double theta, double t, int WM, int VP)
{
    double qR, aR, cR, HpR, uR, cst, k1, k2;
    int j = 1, ok = false, ch, ntrial = 50;

    qR = aR = cR = HpR = 0.0;
    poschar(theta, qR, aR, cR, HpR,WM,VP);

    uR = qR/aR;
    k1 = 1.0/(1.0 + k*(Res1+Res2)/(Res1*Res2*CT));
    k2 = k/(Res1*Res2*CT);
    cst = (k1*(Qold[N]-P(N,Aold[N],WM)/Res1) - qR)/(cR-uR) - aR - HpR*k;

    // Initial guesses

    xr = Anew[N-1];
    f  = 0;
    df = 0;

    while (j <= ntrial && ok==false)
    {
        f  = xr + cst + k1*(1.0/Res1 +k2)*P(N,xr,WM)/(cR-uR);
        df = 1.0 + k1*(1.0/Res1+k2)*dPdA(N,xr,WM)/(cR-uR);
        ch   = zero_1d (&xr, f, df, 1.0e-4);
        if (xr <= 0.0)
        {
            printf("WARNING (arteries.C): Bound_right: x was negative xr = %f t = %f L =%f\n", xr, t, L);
            xr = Anew[N-1]; // Bound xr[1] away from zero.
        }
        if (ch == 1) ok = true;
        j = j+1;
    }
    // Solutions are applied, and right boundary and the intermediate array QL
    // are updated.
    Anew[N] = xr;
    Qnew[N] = k1*(Qold[N] + (P(N,Anew[N],WM)-P(N,Aold[N],WM))/Res1 + k2*P(N,Anew[N],WM));

    // If the solution is not found print an error message. We don't use
    // subroutine error,
    // since it can't take the function values as arguments.
    if (j >= ntrial)
    {
        printf ("WARNING (arteries.C): Root not found in the right boundary, ");
        printf ("x=%f, f=%f, df=%f, j=%d, t=%f\n",xr,f,df,j,t);

        Anew[N]    = Ah[N-1];
        Qnew[N]    = Qh[N-1];
    }
}

void Tube :: call_junc (double theta, double gamma, Tube *Arteries[], int parent, int WM, int VP)
{
    junction(theta,gamma,Arteries,parent,WM,VP);
}

// Solves the non-linear PDE's (momentum and continuity eqn's.
// from t = tstart to t= tend.
//
// This function checks the maximal possible size of the next time-step,
// reduces it to make sure that we don't walk to far, and takes the
// next step. This is done by executing the step routine, then updating
// the left boundary and finally updating bifurcation points and the
// right boundaries. This is carried out as long as the time hasn't passed
// the desired ending time (tend) which is passed to the function as a
// parameter.
void solver (Tube *Arteries[], double tstart, double tend, double k, double Period, int WM, int VP)
{
  // The following definitions only used when a variable time-stepping is
  // used.

  double t    = tstart;
  int qLnb = (int) fmod(t/k,tmstps);

  // As long as we haven't passed the desired ending time do:
  while (t < tend)
  {
    // Check that the step we take is valid. If this error occurs when forcing
    // a constant step-size the program must be terminated.
    if (t+k > tend)
    {
      double kold = k;
      k = tend - t;
      printf("ERROR (arteries.C): Step-size changed, t+k=%10.15f, tend=%10.15f k=%10.15f kold=%10.15f\n",t+kold,tend,k,kold);
    }

    // Check that the CFL-condition applies.
    for (int i=0; i<nbrves; i++)
    {
      if (k > Arteries[i] -> CFL(WM))
      {
        error("arteries.C","Step-size too large CFL-condition violated\n");
//          exit(1);
      }
      
    }
    // solve for interior points, by calling step.
    for (int i=0; i<nbrves; i++)
    {
      Arteries[i] -> step (k,WM,VP);

    }
    // Update left and right boundaries, and the bifurcation points.
    Arteries[0] -> bound_left(t+k, k, Period, WM,VP);
//      fprintf(stdout,"BOUND LEFT\n");
     
    for (int i=0; i<nbrves; i++)
    {
//       fprintf(stdout,"Vessel: %d\n",i);
        if (Arteries[i]->Res1 > 0)
      {
//         fprintf(stdout,"BOUND RIGHT\n");
          Arteries[i] -> bound_right ( k, k/Arteries[i]->h, t,WM,VP);
      }
      else
      {
//         fprintf(stdout,"INTERIOR\n");
//         fprintf(stdout,"Q:%lf\n",Arteries[0]->Qold[0]);
        double theta = k/Arteries[i]->h;
        double gamma = k/2;
        Arteries[i] -> call_junc (theta, gamma, Arteries,i,WM,VP);
      }
    }
    // Update the time and position within one period.
    t = t + k;
    qLnb = (qLnb + 1) % tmstps;

  }
}
