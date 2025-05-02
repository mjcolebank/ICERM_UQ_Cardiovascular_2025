/***************************************************************************/
/*                                                                         */
/*  Program: arteries.C                                                    */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This module can predict the flow and pressure in an tree of elastic    */
/*  vessels as described in IMFUFATEKST NR 297, and D2.1-4. The dependen-  */
/*  cies of the vessels in the tree must be specified in the main module   */
/*  according to the tree in question (for further details see documenta-  */
/*  tion in problem.default.cxx).                                          */
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

// $Id: arteries.C,v 1.17 2010-10-20 15:38:28 mette Exp $

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
extern int size_conn;
extern int tmstps;
extern char* CO_filename;

// A f90 subroutine that determines the impedance at the root of a
// structured tree. Used as BC for the tree consisting of the larger
// arteries.

extern "C" void impedance_driver_(int *tmstps, double *Period,
				  double *ffs1, double *ffs2, double *ffs3,
				  double *rho, double *mu, 
                                  double *r_root, double *rmin,
                                  double *y, double *Lr, double *Fr2, double *q, double *g,
                                  double *par1, double *par2, double *lrr, int *term_ID, int *Z_flag);

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
              double k1r, double k2r, double k3r,
              int *daughters,
              double rmin_in, double points, int init,
              double f1, double f2, double f3, double fs1, double fs2, double fs3, double par1, double par2, double lrr_in, int terminal_id, int Z_flag, double ves_angle):
	L(Length),
	rtop(topradius),
	rbot(botradius),
	daughters(daughters),
    rmin(rmin_in),
	pts(points),
    init(init),
	ff1(f1),
	ff2(f2),
	ff3(f3),
    ffs1(fs1),
    ffs2(fs2),
    ffs3(fs3),
    par1(par1),
    par2(par2),
    lrr(lrr_in),
    term_ID(terminal_id),
    Z_flag(Z_flag),
    ves_angle(ves_angle)
{
  // Initialization of the basic parameters
  N	  = int(pts*L);
  h	  = 1.0/pts/Lr;

//  fprintf(stdout, "K_loss=%8.4f\n",K_loss);

  // Declaration and Initialization of the needed intermediate arrays.
  Qnew	  = new double[N+1];
  Anew	  = new double[N+1];
  Qold	  = new double[N+1];
  Aold	  = new double[N+1];
  Qprv	  = new double[N+1];
  Aprv	  = new double[N+1];
  R1	  = new double[N+1];
  R2	  = new double[N+1];
  S1	  = new double[N+1];
  S2	  = new double[N+1];
  r0	  = new double[N+1];
  r0h	  = new double[N+2];
  dr0dx   = new double[N+1];
  dr0dxh  = new double[N+2];
  A0      = new double[N+1];
  A0h     = new double[N+2];
  fr      = new double[N+1];
  frh     = new double[N+2];
  dfrdr0  = new double[N+1];
  dfrdr0h = new double[N+2];
  p1      = new double[N+1];
  p1h     = new double[N+2];
  dp1dr0  = new double[N+1];
  dp1dr0h = new double[N+2];
  Ah	  = new double[N];
  Qh	  = new double[N];
  R1h	  = new double[N];
  R2h	  = new double[N];
  S1h	  = new double[N];
  S2h	  = new double[N];
  pL	  = new double[tmstps];
  y	      = new double[tmstps];

  double rgLr  = 4.0/3.0/rho/g/Lr;
  double rgLr2 = 4.0/3.0/rho/g/Lr2;

  // Vessel geometry is tabulated and initial conditions are applied
  for (int i=0; i<=N; i++)
  {
    
    r0 [i]     = (k1r*exp(i*(-k2r*Lr)*L/N)+k3r/Lr);  //if x[0,1] do not need L
    r0h[i]     = (k1r*exp((i-0.5)*(-k2r*Lr)*L/N)+k3r/Lr);   //deltax = r/N 
    dr0dx [i]  = L*-k2r*Lr*(k1r*exp(i*(-k2r*Lr)*L/N))/h/N;  //do we need L here? 
    dr0dxh[i]  = L*-k2r*Lr*(k1r*exp((i-0.5)*(-k2r*Lr)*L/N))/h/N;

    A0 [i]     = M_PI*sq(r0 [i]);
    A0h[i]     = M_PI*sq(r0h[i]);
    fr [i]     = (ff1*exp(ff2*r0 [i])+ff3)*rgLr;
    frh[i]     = (ff1*exp(ff2*r0h[i])+ff3)*rgLr;
    dfrdr0 [i] = ff1*ff2*exp(ff2*r0 [i])*rgLr2;
    dfrdr0h[i] = ff1*ff2*exp(ff2*r0h[i])*rgLr2;
    p1 [i]     = fr [i]/M_PI;
    p1h[i]     = frh[i]/M_PI;
    dp1dr0 [i] = dfrdr0 [i]/M_PI;
    dp1dr0h[i] = dfrdr0h[i]/M_PI;
    Qnew[i]    = 1e-12;//0.0;
    Anew[i]    = A0[i];
  }

  r0h[N+1]     = (k1r*exp((N+0.5)*(-k2r*Lr)*L/N)+k3r/Lr);       
  dr0dxh[N+1]  = L*-k2r*Lr*(k1r*exp((N+0.5)*(-k2r*Lr)*L/N))/h/N;

  
  A0h[N+1]     = M_PI*sq(r0h[N+1]);
  frh[N+1]     = (ff1*exp(ff2*r0h[N+1])+ff3)*rgLr;
  dfrdr0h[N+1] = ff1*ff2*exp(ff2*r0h[N+1])*rgLr2;
  p1h[N+1]     = frh[N+1]/M_PI;
  dp1dr0h[N+1] = dfrdr0h[N+1]/M_PI;
  
//  fprintf(stdout,"k1r %lf, k2r %lf, k3r %lf, r0[0] %lf r0[N] %lf, L %lf\n",k1r, k2r, k3r, r0[0],r0[N], L);
  
  // Read from file data for the inflow profile.
  if (init == 1)
  {
    Q0 = new double[tmstps+1];

    FILE *fi = fopen (CO_filename, "r");
      if (fi); //fprintf(stdout, "Q0 opened\n");
    else error ("arteries.C"," Q0 NOT OK");

    for (int i=0; i<=tmstps; i++)
    {
      fscanf(fi,"%lf",&Q0[i]);
      Q0[i] = Q0[i]/q; // If the indata have dimensions they should be made
                       // non-dimensional.
    }
  }

  // In case of an end-tube evaluate the impedances for the boundary condition.
  // This is done by calling the f90 routine root_imp which calculates the
  // impedance at the root of a structured tree. The underscores is sensitive
  // to the compiler but can be seen at the bottom of the file root_imp.o.
  if (rmin > 0)
  {
//    fprintf(stdout,"Calling f90 subroutines\n");
//          fprintf(stdout,"\n");
//      printf("r_root: %lf, r_min: %lf, alpha: %lf, beta: %lf, lrr: %lf\n",rbot,rmin,par1,par2,lrr);
    impedance_driver_(&tmstps,&Period,&ffs1,&ffs2,&ffs3,&rho,&mu_pl,&rbot,&rmin,y,&Lr,&Fr2,&q,&g,&par1,&par2,&lrr,&term_ID,&Z_flag);
//    printf("Finished with f90 subroutines.\n\n");

    // Initialize the array pL used when determining the convolution
    // in the right boundary condition (see the subroutine bound_right).
    for (int j=0; j<tmstps; j++)
    {
      pL[j] = 0.0;
    };
  };
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
  delete[] p1;
  delete[] p1h;
  delete[] dp1dr0;
  delete[] dp1dr0h;
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
void Tube :: printPt (FILE *fd, double t, int i,int WM)
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

void Tube :: printFt (FILE *fd, double t, int i)
{
  fprintf (fd, "%13.10f %15.10f\n", t*Lr3/q, F(Qnew[i],Anew[i])*sq(q)/Lr3);
}

// ADDED BY MJC

void Tube :: printAllt (FILE *fd,double t, int i,int WM)
{
    fprintf (fd, "%13.10f %13.10f %15.10f %15.10f %15.10f %17.10f\n", t*Lr3/q, i*h*Lr, (P(i,Anew[i],WM)+p0)*rho*g*Lr/conv, Qnew[i]*q, Anew[i]*Lr2, c(i, Anew[i]*10,WM));
}

// The following functions prints P, Q, A, and F as functions of
// (x, t). This is done in terms of the re-dimensionalized variables.
// In this case the functions is plotted for a
// fixed time, but for all x along the vessel in question. Since the
// doesn't have to be the first vessel in the tree, it would have
// some offset from the heart. Which determines the position for x.
// Therefore there are two arguments passed to this function the time
// and the offset.
void Tube :: printPxt (FILE *fd, double t, int offset,int WM)
{
  if (offset == 0) fprintf (fd, "\n");
  for (int i=0; i<N; i++)
  {
    //fprintf (fd, "%13.10f %13.10f %15.10f\n",
      //       t*Lr3/q, (i+offset)*h*Lr, (P(i,Anew[i])+p0)*rho*g*Lr/conv);
      fprintf (fd, "%13.10f %13.10f %15.10f %15.10f %15.10f %17.10f\n",
               //t*Lr3/q, (i+offset)*h*Lr, (P(i,Anew[i])+p0)*rho*g*Lr/conv, (Qnew[i]*q)/(Anew[i]*Lr2), 2*sqrt((Anew[i]*Lr2)/M_PI), c(i, Anew[i])*10);
               t*Lr3/q, (i+offset)*h*Lr, (P(i,Anew[i],WM)+p0)*rho*g*Lr/conv, Qnew[i]*q, Anew[i]*Lr2, c(i, Anew[i]*10,WM));
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

void Tube :: printFxt (FILE *fd, double t, int offset)
{
  for (int i=0; i<=N; i++)
  {
    fprintf (fd, "%13.10f %13.10f %15.10f\n",
             t*Lr3/q, (i+offset)*h*Lr, F(Qnew[i],Anew[i])*sq(q)/Lr3);
  }
}

// A function that prints p(Q) for all t. This is done in terms of
// the re-dimensionalized variables. In this case the plot is made for a
// fixed point in space, but for all t along the vessel in question.
void Tube :: printPQ (FILE *fd, int i,int WM)
{
  fprintf (fd,"%15.10f %15.10f\n",
           Qnew[i]*q, (P(i,Anew[i],WM)+p0)*rho*g*Lr/conv);
}

void Tube :: printPA (FILE *fd, int i,int WM)
{
  fprintf (fd,"%15.10f %15.10f\n",
           (P(i,Anew[i],WM)+p0)*rho*g*Lr/conv, Anew[i]*Lr2);
}

// Plotting the terms in the continuity equation on dimension-less form.
void Tube :: printdQdx (FILE *fd, double t, int i)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           (Qnew[i+1]-Qnew[i-1])/2/h);

}

void Tube :: printdAdt (FILE *fd, double t, int i, double Aprev, double tmst)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           (Anew[i]-Aprev)/tmst);
}

void Tube :: printTotConRes (FILE *fd, double t, int i, double Aprev, double tmst)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
          (Qnew[i+1] - Qnew[i-1])/2/h +
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
           (sq(Qnew[i+1])/Anew[i+1] - sq(Qnew[i-1])/Anew[i-1])/2/h);

}

void Tube :: printdPdx (FILE *fd, double t, int i,int WM)
{
    fprintf (fd, "%13.10f %15.10f\n", t,
             Anew[i]*(P(i+1,Anew[i+1],WM)-P(i-1,Anew[i-1],WM))/Fr2/2/h);
}

void Tube :: printFric (FILE *fd, double t, int i,int WM)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           F(Qnew[i],Anew[i]));
}

void Tube :: printTotMomRes (FILE *fd, double t, int i, double Qprev, double tmst, int WM)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           Anew[i]*(P(i+1,Anew[i+1],WM)-P(i-1,Anew[i-1],WM))/Fr2/2/h +
           (sq(Qnew[i+1])/Anew[i+1] - sq(Qnew[i-1])/Anew[i-1])/2/h +
           (Qnew[i]-Qprev)/tmst + F(Qnew[i],Anew[i]));
}

// Further print functions can be added, and they would look similar to
// the two functions above!

// The next function returns the pressure p as a function of a fixed x,
// and the corresponding cross-sectional area A. The pressure is defined
// according to the  mathematical model, described in IMFUFATEKST no 297,
// and D2.1-4.
double Tube :: P (int i, double A, int WM)
{
    double pold;
    if (WM==1){
        pold = fr[i]*(1.0-sqrt(A0[i]/A));
    }
    else {
        pold = fr[i]*(sqrt(A/A0[i])-1.0);  // Linear B
    }

  return pold;
}

double Tube :: dPdA (int i, double A, int WM)
{
    double pold;
    if (WM==1){
        pold = 0.5*fr[i]*sqrt(A0[i]/cu(A));
    }
    else {
        pold = 0.5*fr[i]/sqrt(A0[i]*A);  // Linear B
    }
  return pold;
}

double Tube :: dPdx1(int i, double A, int WM)
{
    double pold;
  if (WM==1){
      pold = (dfrdr0[i]*(1-sqrt(A0[i]/A))-fr[i]*sqrt(M_PI/A))*dr0dx[i];
  }
  else {
      pold = (dfrdr0[i]*(sqrt(A/A0[i])-1)-fr[i]*sqrt(M_PI*A)/A0[i])*dr0dx[i];  // Linear B
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

double Tube :: dBdx1 (int i, double A, int WM)
{

    double dfr = dfrdr0[i];
    double pold, inteval;
    if (WM==1){
           pold = dr0dx[i]/Fr2*(
                  sqrt(A)*(2*sqrt(M_PI)*fr[i]+2*sqrt(A0[i])*dfr)-
                  A*dfr-2.0*M_PI*r0[i]*fr[i]-A0[i]*dfr);
        
        }
       else {

           inteval = 2.0*M_PI*r0[i]*fr[i]+A0[i]*dfr;
           pold = dr0dx[i]*(2.0*sqrt(cu(A))*(sqrt(M_PI)*fr[i]-sqrt(A0[i])*dfr)/A0[i]/3.0
                                   + A*dfr - inteval/3.0)/Fr2; // Linear B
           
       }


  return pold;
}

double Tube :: dBdx1h (int i, double A, int WM)
{
    int ip1 = i+1;

    double dfr = dfrdr0h[ip1];
    double pold,inteval;
  
    if (WM==1){
        pold = dr0dxh[ip1]/Fr2*(
             sqrt(A)*(2*sqrt(M_PI)*frh[ip1]+2*sqrt(A0h[ip1])*dfr)-
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

double Tube :: d2BdAdxh (int i, double A, int WM)
{
   int ip1 = i+1;

   double dfr = dfrdr0h[ip1];
    double pold;
  if (WM==1){
        pold = (-dfr+1/sqrt(A)*(sqrt(M_PI)*frh[ip1]+
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
//      if (i==0)//(Vnew>c_tmp)
//      fprintf(stdout,"Rtop: %lf C:%lf U:%lf\n",rtop,c_tmp,Vnew);
  }
  return (minimum);
}

// When taking a Lax-Wendroff step, the flux of the system must be determined.
// This is evaluated at i + j/2, and the prediction is given as described
// in IMFUFATEKST no 297 and D2.1-4. The integer k determines whether we deal
// with the first or the second component of the vector.
double Tube :: Rvec (int k, int i, int j, double Q, double A, int WM)
{
  if(k==1) return(Q); else
  if(k==2) return(sq(Q)/A + ((j==0)?B(i,A,WM):Bh(i,A,WM)));
  else error ("arteries.cxx","Call of non-existing vector-component of R");
  return(0);
}

// Similarly the right hand side of the system of equations must be determined
// at i + j/2. Also in this case the function is given as stated in
// the mathematical model, and also in this case k states the needed component
// of the vector.
double Tube :: Svec (int k, int i, int j, double Q, double A, int WM)
{   
    
//  fprintf(stdout, "G: %lf\n", G(A,ves_angle));
  if(k==1) return(0.0); else
  if(k==2) return(F(Q,A) + G(ves_angle, A) + ((j==0)?dBdx1(i,A,WM):dBdx1h(i,A,WM))); // -gAcostheta
  else error ("arteries.cxx","Call of non-existing vector-component of S");
  return(0);

}

// The solutions of Anew and Qnew are found for all interior points
// of the vessel at (t+k), where k is the length of the current
// time-step. This function saves the results in the arrays Anew and
// Qnew, and the function is made according to Lax-Wendroff's method
// as described in IMFUFATEKST no 297 and D2.1-4.
void Tube :: step (double k, int WM)
{
  double theta = k/h;  // Theta is determined.
  double gamma = 0.5*k;  // Gamma is determined.

  for (int i=0; i<=N; i++)  // Remember the values at this time level.
  {
    Qold[i] = Qnew[i];
    Aold[i] = Anew[i];
  }

  // Anew and Qnew are predicted at the new time level (t+k).
  for (int i=0; i<=N; i++)
  {
    R1[i] = Rvec(1,i,0,Qold[i],Aold[i],WM);
    R2[i] = Rvec(2,i,0,Qold[i],Aold[i],WM);
    S1[i] = Svec(1,i,0,Qold[i],Aold[i],WM);
    S2[i] = Svec(2,i,0,Qold[i],Aold[i],WM);
  }

  for (int i=0; i<N; i++)
  {
    Ah[i]  = 0.5*(Aold[i+1]+Aold[i]) - 0.5*theta*(R1[i+1]-R1[i]) +
	     0.5*gamma*(S1[i+1]+S1[i]);
      Qh[i]  = 0.5*(Qold[i+1]+Qold[i]) - 0.5*theta*(R2[i+1]-R2[i]) +
	     0.5*gamma*(S2[i+1]+S2[i]);

    R1h[i] = Rvec(1,i,1,Qh[i],Ah[i],WM);
    R2h[i] = Rvec(2,i,1,Qh[i],Ah[i],WM);
    S1h[i] = Svec(1,i,1,Qh[i],Ah[i],WM);
    S2h[i] = Svec(2,i,1,Qh[i],Ah[i],WM);
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
void Tube :: bound_left (double t, double k, double Period)
{
  Qnew[0]   = Q0_init(t,k,Period);

  if (int(t/k) < 0)
    printf("t/k negative in bound_left\n");
  double Qhm05 = Qnew[0]+Q0_init(t-k,k,Period) - Qh[0];
  double R1hm05    = Qhm05;
  Anew[0]   = Aold[0] - k*(R1h[0] - R1hm05)/h;
}

// The value at the right boundary at time t is predicted. NB: This should
// only be used with terminal vessels, i.e. for vessels that don't bifurcate
// into further branches.
// In that situation the bifurcation boundary function should be called
// instead. Again the procedure specified is given according to the mathemati-
// cal theory presented in IMFUFATEKST no 297 and D2.1-4.

double Tube :: c (int i, double A, int WM) // The wave speed through aorta.
{
    double cnst;
    if (WM==1){
        cnst =  0.5*fr[i]*sqrt(A0[i]/A)/Fr2;
    }
    else {
        cnst =  0.5*fr[i]*sqrt(A/A0[i])/Fr2;  //Linear B
    }

  return sqrt (cnst);
}

double Tube :: Hp (int i, double Q, double A, int WM)
{
  return (F(Q,A) - A*dPdx1(i,A,WM)/Fr2)/(-Q/A + c(i,A,WM));
}

void Tube :: poschar (double theta, double &qR, double &aR, double &cR, double &HpR, int WM)
{
  double ctm1  = c  (N, Aold[N], WM);
  double Hptm1 = Hp (N, Qold[N], Aold[N], WM);
  double uR    = Qold[N] / Aold[N];
  double ch    = (uR + ctm1) * theta;

  if (uR + ctm1 < 0)
  {
    printf("uR + ctm1 < 0, CFL condition violated\n");
  }

  qR  = Qold[N] - (Qold[N] - Qold[N-1])*ch;
  aR  = Aold[N] - (Aold[N] - Aold[N-1])*ch;
  cR  = ctm1    - (ctm1  - c (N-1,Aold[N-1],WM))*ch;
  HpR = Hptm1   - (Hptm1 - Hp(N-1,Qold[N-1],Aold[N-1],WM))*ch;
}

void Tube :: bound_right (int qLnb, double k, double theta, double t, int WM)
{
  int j = 1;
  int ntrial = 500;

  int qLnb_1 = qLnb + 1;

  // Make sure that qLnb_1 runs in the interval [0:tmstps-1].
  if (qLnb_1 == (int) tmstps)
  {
    qLnb_1 = 0;
  }

  // In order to make a relation between P(x_L, t+dt) and Q(x_L, t+dt), and
  // P(x_L, t+dt/2) and Q(x_L,t+dt/2) we need to extract the term involving
  // y[0] (see mathematical derivation). 
  // This term corresponds to the peripheral
  // The remaining terms in the convolution present at the boundary,
  // see mathematical derivation.
  double pterms = 0.0;

  if (t > Period)
  {
    for (int m=1; m<tmstps; m++)
    {
      int pindex  = (qLnb_1 + tmstps - m) % tmstps;
      pterms  = pterms  + (pL[pindex])*y[m];
    }
    pterms  = k*pterms;
  }

  double qR = 0.0;
  double aR = 0.0;
  double cR = 0.0;
  double HpR = 0.0;

  poschar(theta, qR, aR, cR, HpR, WM);

  double uR = qR/aR;

  double cst = (pterms - qR) / (cR -uR) - aR - HpR*k;

  xr = Anew[N-1];    // Initial guess (A_L^n+0.5)
  //finit = xr + cst + y[0]*k*(P(N,xr))/(cR-uR);
  //dfinit =1 + y[0]*k*dPdA(N,xr)/(cR-uR);

  // The value of the function and the derivative is initialized to 0.
  f  = 0;
  df = 0;

  // Solve the equation using Newton Raphson's method
  // (using subroutine zero_1d in tools.C.
  bool ok = false;
  while (j <= ntrial && ok==false)
  {
    f  = xr + cst + y[0]*k*(P(N,xr,WM))/(cR-uR);
    df = 1.0 + y[0]*k*dPdA(N,xr,WM)/(cR-uR);
    ok  = zero_1d (&xr, f, df, 1.0e-4);
    if (xr <= 0.0 && j==1)
    {
      printf("WARNING (arteries.C): Bound_right: x was negative xr = %f t = %f L =%f\n", xr, t, L);
      xr = Anew[N-1]; // Bound xr away from zero.
    }
    j = j+1;
  }
  // Solutions are applied, and right boundary and the intermediate array QL
  // are updated.
  Anew[N]    = xr;
  pL[qLnb_1] = P(N,Anew[N],WM); // ADD P0???
  Qnew[N]    = k*pL[qLnb_1]*y[0] + pterms;

  // If the solution is not found print an error message. We don't use
  // subroutine error,
  // since it can't take the function values as arguments.
  if (j >= ntrial)
  {
    printf ("WARNING (arteries.C): Root not found in the right boundary, ");
    printf ("x=%f, f=%f, df=%f, j=%d, t=%f\n",xr,f,df,j,t);

    Anew[N]    = Ah[N-1];
    Qnew[N]    = Qh[N-1];
    pL[qLnb_1] = P(N,Ah[N-1],WM);
  }
}

// The value at the bifurcation point at time t is predicted. NB: This should
// only be done for tubes that do bifurcate into further branches. If
// this is not the case we have a terminal vessel and bound_right should be
// called instead. The procedure operates according to the specifications
// in the mathematical model as a link between this tube and its daughters.
// Therefore there will be three tubes involved in this function.
// One problem is however, that the rather complicated system of equations does
// not converge for all choices of parameters (the peripheral resistance, the
// top radius, and the bottom radius).

void Tube :: call_junc (double theta, double gamma, Tube *Arteries[], int parent, int conn_id, int WM)
{
    junction(theta,gamma,Arteries,parent,conn_id,WM);
}

///////////////////////////////////////////////////////////////
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

void solver (Tube *Arteries[], double tstart, double tend, double k, int WM)
{
  // The following definitions only used when a variable time-stepping is
  // used.
  double t    = tstart;
  int qLnb = (int) fmod(t/k,tmstps);
    // MJC: junction index updated when we have converging flow
    int curr_ves=0;

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
          fprintf(stdout,"dt: %lf CFL: %lf\n",k,Arteries[i] -> CFL(WM));
      }
    }

    // solve for interior points, by calling step.

    for (int i=0; i<nbrves; i++)
    {
      Arteries[i] -> step (k,WM);
    }
    // Update left and right boundaries, and the bifurcation points.
      
    Arteries[0] -> bound_left(t+k, k, Period);
    //fprintf(stdout,"r0[0] %lf r0[N] %lf\n", Arteries[0]->r0[0], Arteries[0]->r0[Arteries[0]->N]);
      
    /* MJC: Instead of running from i=1 -> nbrves, walk through
            the connectivity matrix to determine which vessel
            is next in the array.
     */
    for (int i=0; i<=size_conn; i++)
   {
       // MJC increment based on connectivity matrix
       curr_ves = Arteries[0]->daughters[max_D*i];
//       fprintf(stdout,"ves: %d\n",curr_ves);
        // Here, we test if there is a converging vessel that is terminal
        if (Arteries[curr_ves] ->rmin>0 && Arteries[curr_ves]->daughters[max_D*i+1]!=0) //daughters[max_D*i+1]
        {
                double theta = k/Arteries[curr_ves]->h;
            double gamma = k/2.0;
            
//                          fprintf(stdout,"IN BOUND BIF_CONV %d\n", curr_ves);
            Arteries[i] -> call_junc (theta, gamma, Arteries, curr_ves,i, WM); //increment or add another flag

//              fprintf(stdout,"IN BOUND RIGHT %d\n", curr_ves);
            Arteries[i] -> bound_right (qLnb, k, k/Arteries[curr_ves]->h, t, WM);
           
                 
        }
        // Here, we test if a vessel is terminal
        else if (Arteries[curr_ves] ->rmin > 0)
          {
//               fprintf(stdout,"IN BOUND RIGHT %d\n", curr_ves);
              Arteries[curr_ves] -> bound_right (qLnb, k, k/Arteries[curr_ves]->h, t, WM);
          }
            
        // ELSE IF INIT==0 DO BOUND_BIF, IF INIT==-1 DO NOTHING
        else if (Arteries[curr_ves] -> init == 0 || Arteries[curr_ves] -> init == 1)
          {
            double theta = k/Arteries[curr_ves]->h;
            double gamma = k/2;
//               fprintf(stdout,"IN BOUND BIF %d\n", curr_ves);
            Arteries[curr_ves] -> call_junc (theta, gamma, Arteries, curr_ves,i, WM);
          }
        }
    // Update the time and position within one period.
    t = t + k;
    qLnb = (qLnb + 1) % tmstps;
//     fprintf(stdout,"t:%lf\n",t);
//     exit(1);
//     if (t>0.1)
//     {
//         exit(1);
//     }
    // since the outer loop is time, restart junc_flag at zero
  }
}
