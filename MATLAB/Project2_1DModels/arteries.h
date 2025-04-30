/***************************************************************************/
/*                                                                         */
/*  Program: arteries.h                                                    */
/*  Version: 1.1                                                           */
/*  Date: 31 March 2021                                                    */
/*                                                                         */
/*  Primary Authors: M.S. Olufsen                                          */
/*  Key Contributers: M.U. Qureshi & M.J. Colebank                         */
/*  Institute: North Carolina State University, Raleigh, NC, USA           */
/*  Funded by: NSF-DMS # 1615820 & AHA #19PRE34380459                      */
/*                                                                         */
/*  This header file includes the definition of the global constants, and  */
/*  the definition of the class representing one blood vessel.             */
/*                                                                         */
/***************************************************************************/

#ifndef _ARTERIES_H
#define _ARTERIES_H

#include <cstdio>
#include <cmath>

// Global parameters imported from main.h

extern double   conv, rho, mu, nu, Lr, Lr2, Lr3, g, q, Fr2,
                Re, p0, pmean, tmst, Period, Fcst, CO, COm,
                Deltat, tau,
            *fjac[24], xr, f, df;

extern int max_D;
extern double bound_thick;
extern double vel_power, Fcst;

// The class structure.
class Tube {
public:
  double L;                    // The length of the vessel
  double rtop, rbot;           // The top and bottom radii of the vessel
    int *daughters;          //MJC: try to assign the tube daughters without passing Tube
    
    
  double pts;                  // The number of grid points per cm
  double ff1, ff2, ff3;
  double Res1, Res2, CT;
  double sten_factor, sten_length;

  int N;                       // The number of grid points along the vessel
  double h;                    // The interval length of delta x
  double RLrb;                 // The peripheral resistance of the vessel

  double *Qnew, *Qold, *Qh,    // The arrays needed to store data during
         *Anew, *Aold, *Ah,    // the numerical solution of the system.
         *R1, *R2, *R1h, *R2h,
         *S1, *S2, *S1h, *S2h,
         *Qprv, *Aprv,
         *pL, *y, *QL, *Q0, *Z,
     *r0, *r0h,
     *dr0dx, *dr0dxh,
     *A0, *A0h, *wom,
     *fr, *frh,
    *dfrdr0, *dfrdr0h;

    //MJC new version
    Tube (double Length,
          double topradius, double botradius,
          int *daughters,
          double points, int init,
          double f1, double f2, double f3,
          double Res1_in, double Res2_in, double CT_in, double sten_factor, double sten_length);
                                                         // Constructor.
  ~Tube ();                                              // Destructor.


 // Prints P(x_fixed,t), A(x_fixed,t), F(x_fixed,t), or Q(x_fixed,t) for all
  // t's along the tube.
  void printQ0  (FILE *fd);
  void printPt  (FILE *fd, double t, int i, int WM);
  void printQt  (FILE *fd, double t, int i);
  void printAt  (FILE *fd, double t, int i);
  void printFt  (FILE *fd, double t, int i, int VP);


  // Prints P(x,t), A(x,t), F(x,t), or Q(x,t) for all x's and t's
  // along the tube. The argument offset makes sure that the vessel
  // is located with the right offset from the inlet.
  void printPxt (FILE *fd, double t, int offset, int WM);
  void printAxt (FILE *fd, double t, int offset);
  void printFxt (FILE *fd, double t, int offset, int VP);
  void printQxt (FILE *fd, double t, int offset);
    
  void printAllt (FILE *fd, double t, int i, int WM);

  // Prints P as a function of Q and A, respectively,  at a given location x.
  void printPQ (FILE *fd, int i, int WM);
  void printPA (FILE *fd, int i, int WM);

  // Prints dP/dx(x_fixed,t), dQ/dx(x_fixed,t), dA/dt(x_fixed,t),
  // dQ/dt(x_fixed,t), Fric(x_fixed, t), TotConRes(x_fixed,t),
  // TotMomRes(x_fixed,t) for any fixed x and for all t.
  // I.e. all terms including the sums of the momentum and
  // continuity equations.
  void printdAdt      (FILE *fd, double t, int i, double Aprev, double tmst);
  void printdQdx      (FILE *fd, double t, int i);
  void printTotConRes (FILE *fd, double t, int i, double Qprev, double tmst);

  void printdQdt      (FILE *fd, double t, int i, double Qprev, double tmst);
  void printddxQ2divA (FILE *fd, double t, int i);
  void printdPdx      (FILE *fd, double t, int i, int WM);
  void printFric      (FILE *fd, double t, int i, int VP);
  void printTotMomRes (FILE *fd, double t, int i, double Qprev, double tmst, int WM, int VP);


  // Defines P(x,A(x,t)).
  double P (int i, double A, int WM);

  // Defines dPdA(x,A(x,t)).
  double dPdA (int i, double A, int WM);

  // Defines dPdx1(x,A(x,t)).
  double dPdx1 (int i, double A, int WM);

  // Defines B(x,A(x,t)).
  double B (int i, double A, int WM);

  // Defines Bh(x,A(x,t)).
  double Bh (int i, double A, int WM);

  // Defines dBdx1(x,A(x,t)).
  double B_RHS (int i, double A, int WM);

  // Defines dBdx1h(x,A(x,t)).
  double B_RHSh (int i, double A, int WM);

  // Defines dBdAh (x,A(x,t)).
  double dBdAh (int i, double A, int WM);

  // Defines d2BdAdxh (x, A(x,t));
  double dB_RHSdAh (int i, double A, int WM);


  // Tests that the CFL-condition is valid throughout the tube.
  double CFL (int WM);


  // Finds the flux acc. to sys. eq.
  double Rvec (int k, int i, int j, double Q, double A, int WM, int VP);

  // Finds the rhs. of system eq.
  double Svec (int k, int i, int j, double Q, double A, int WM, int VP);


  // Steps through interior points.
  void step (double k, int WM, int VP);
    
  // Use characteristics
    double c  (int i, double A, int WM); // The wave speed through aorta.
    double Hp (int i, double Q, double A, int WM, int VP);
    double Hn (int i, double Q, double A, int WM, int VP);
    void poschar (double theta, double &qR, double &aR, double &cR, double &HpR, int WM, int VP);
    void negchar (double theta, double &qS, double &aS, double &cS, double &HnS, int WM, int VP);


  // Updates left bndry. This should only be done for the inlet tube.
  void bound_left (double t, double k, double Period, int WM, int VP);


  // Updates right bndry. This should only be done for terminal vessels.
  void bound_right ( double k, double theta, double t, int WM, int VP);


  // Updates bifurcation conditions. Uses daughter vessels, and should
  // only be called when such exist.
    // CALLED IN JUNCTION.C

    void call_junc (double theta, double gamma, Tube *Arteries[], int parent, int WM, int VP);
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
    
private:
  // The private function Q0 may only be accessed from the left boundary
  // function. It ensures a certain and given CO (defined in main.h).
  double Q0_init (double t, double k, double Period);
};

void solver (Tube *Arteries[], double tstart, double tend, double k, double Period, int WM, int VP);

#endif
