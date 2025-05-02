/***************************************************************************/
/*                                                                         */
/*  Program: arteries.h                                                    */
/*  Version: 1.0                                                           */
/*  Date: 30 Dec. 2019                                                     */
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

extern double   conv, rho, mu, mu_pl, nu, Lr, Lr2, Lr3, g, q, Fr2,
                Re, p0, pmean, tmst, Period, Fcst, CO, COm,
                Deltat, tau,// m1, m2, m3, m4,
	        *fjac[24], xr, f, df;
extern int max_D;
extern double bound_thick; // MJC


// The class structure.
class Tube {
public:
  double L;                    // The length of the vessel
  double rtop, rbot;           // The top and bottom radii of the vessel
  double k1r, k2r, k3r;        //
  int *daughters;
  double rmin;                      // should be set to zero.
  int pts;                  // The number of grid points per cm
  int init;
  double ff1, ff2, ff3;
  double ffs1, ffs2, ffs3;
  double par1, par2, lrr;
  int term_ID;
  int Z_flag;
  double ves_angle;
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
	 *A0, *A0h, *wom, *Cu,
	 *fr, *frh,
	 *dfrdr0, *dfrdr0h,
	 *p1, *p1h,
	 *dp1dr0, *dp1dr0h;

  Tube (double Length,
        double topradius, double botradius,
        double k1r, double k2r, double k3r,
        int *daughters,
        double rmin, double points, int init,
        double f1, double f2, double f3,
        double fs1, double fs2, double fs3,
        double par1, double par2, double lrr,
        int term_ID, int Z_flag, double ves_angle);
                                                         // Constructor.
  ~Tube ();                                              // Destructor.


  // Prints P(x_fixed,t), A(x_fixed,t), F(x_fixed,t), or Q(x_fixed,t) for all
  // t's along the tube.
  void printQ0  (FILE *fd);
  void printPt  (FILE *fd, double t, int i, int WM);
  void printQt  (FILE *fd, double t, int i);
  void printAt  (FILE *fd, double t, int i);
  void printFt  (FILE *fd, double t, int i);
    // ADDED BY MJC
    void printAllt (FILE *fd, double t, int i, int WM);

  // Prints P(x,t), A(x,t), F(x,t), or Q(x,t) for all x's and t's
  // along the tube. The argument offset makes sure that the vessel
  // is located with the right offset from the inlet.
  void printPxt (FILE *fd, double t, int offset ,int WM);
  void printAxt (FILE *fd, double t, int offset);
  void printFxt (FILE *fd, double t, int offset);
  void printQxt (FILE *fd, double t, int offset);

  // Prints P as a function of Q and A, respectively,  at a given location x.
  void printPQ (FILE *fd, int i ,int WM);
  void printPA (FILE *fd, int i ,int WM);

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
  void printFric      (FILE *fd, double t, int i, int WM);
  void printTotMomRes (FILE *fd, double t, int i, double Qprev, double tmst, int WM);


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
  double dBdx1 (int i, double A, int WM);

  // Defines dBdx1h(x,A(x,t)).
  double dBdx1h (int i, double A, int WM);

  // Defines dBdAh (x,A(x,t)).
  double dBdAh (int i, double A, int WM);

  // Defines d2BdAdxh (x, A(x,t));
  double d2BdAdxh (int i, double A, int WM);


  // Tests that the CFL-condition is valid throughout the tube.
  double CFL (int WM);


  // Finds the flux acc. to sys. eq.
  double Rvec (int k, int i, int j, double Q, double A, int WM);

  // Finds the rhs. of system eq.
  double Svec (int k, int i, int j, double Q, double A, int WM);


  // Steps through interior points.
  void step (double k, int WM);


  // Updates left bndry. This should only be done for the inlet tube.
  void bound_left (double t, double k, double Period);


  // Updates right bndry. This should only be done for terminal vessels.
  double c  (int i, double A, int WM); // The wave speed.
  double Hp (int i, double Q, double A, int WM);
  void poschar (double theta, double &qR, double &aR, double &cR, double &HpR, int WM);
  void bound_right (int qLnb, double k, double theta, double t, int WM);

  // Updates bifurcation conditions. Uses daughter vessels, and should
  // only be called when such exist.
    void call_junc (double theta, double gamma, Tube *Arteries[], int parent, int conn_id, int WM);
    

// A function returning the Friction of the system. The definition of this
// function is given according to the derivation in the mathematical model.
// The constant cst, determines the amount of damping in the system.
inline double F (double Q, double A)
  {
  //  double tmp1 = Fcst*M_PI*Q;
  //  double tmp2 = A*Re;
      double tmp1 = -2.0*sqrt(M_PI)*Q; //MJC
      double tmp2 = bound_thick*Re*sqrt(A);         //MJC
    double tmp3 = tmp1/tmp2;
//      fprintf(stdout,"Re: %lf F:%lf\n",Re,tmp3);
    return(tmp3);
  }

inline double dFdQ (double A)
{
//  return(-Fcst*M_PI/(A*Re));
    return((-2.0*sqrt(M_PI))/(bound_thick*Re*sqrt(A))); //MJC
}
    
inline double dFdA (double Q, double A)
{
//  return(Fcst*M_PI*Q/(sq(A)*Re));
    double tmp1 = Q*sqrt(M_PI);
    double tmp2 = bound_thick*Re*sqrt(cu(A));
    return(tmp1/tmp2); //MJC
}
    
    // MJC: When modeling the pulmonary circulation, it is important
    // to account for gravitational effects
inline double G (double angle, double A)
    {
        double theta_rad = angle*M_PI/180.0;
        
//         fprintf(stdout, "G: %lf A: %lf g: %lf cos: %lf F: %lf T: %lf\n", -1.0*A*g*cos(theta_rad)/Fr2, A, g, cos(theta_rad), Fr2, theta_rad);
        return(1.0*A*cos(theta_rad)/Fr2);
    }
    
    inline double dGdA (double angle,double A)
    {
        double theta_rad = angle*M_PI/180.0;
        return(1.0*cos(theta_rad)/Fr2);
    }

private:
  // The private function Q0 may only be accessed from the left boundary
  // function. It ensures a certain and given CO (defined in main.h).
  double Q0_init (double t, double k, double Period);
};

void solver (Tube *Arteries[], double tstart, double tend, double k, int WM);

#endif
