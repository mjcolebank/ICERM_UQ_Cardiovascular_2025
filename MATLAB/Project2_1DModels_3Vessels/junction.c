/***************************************************************************/
/*                                                                         */
/* The junction.C main program                                             */
/*  Version: 1.0                                                           */
/*  Date: 24 March 2020                                                    */
/*                                                                         */
/*  Primary Authors: M.J. Colebank & J. Mackenzie                          */
/*  Key Contributers: M.U. Qureshi & M.S. Olufsen                          */
/*  Institute: North Carolina State University, Raleigh, NC, USA           */
/*  Funded by: NSF-DMS # 1615820 & AHA #19PRE34380459                      */
/*                                                                         */
/*  In contrast to previous versions, this code calls a generalized        */
/*  junction condition, which can be a single vessel, bifurcation, or n-   */
/*  furcation, and consructs the appropriate Jacobian.                     */
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

extern int max_D;

void junction (double theta, double gamma, Tube *Arteries[], int parent, int WM, int VP)
{
    int D1,D2,D3;//,D4;
    D1 = Arteries[parent]->daughters[max_D*parent + 1];
    D2 = Arteries[parent]->daughters[max_D*parent + 2];
    D3 = Arteries[parent]->daughters[max_D*parent + 3];
    

    if (D2==0) {
       bound_monf(theta, gamma, Arteries, parent,D1,WM,VP); // Only call mono if wanting a stenosis
    }
    else if (D3==0){
//        fprintf(stdout,"Call bif with D1:%d D2:%d D3:%d\n",D1,D2,D3);
        bound_bif(theta, gamma, Arteries, parent,D1,D2,WM,VP);
    }
    else if (D3==-1){
//           fprintf(stdout,"Call sten bif with D1:%d D2:%d D3:%d\n",D1,D2,D3);
            bound_sten_bif(theta, gamma, Arteries, parent,D1,D2,WM,VP);
    }
    else if (D3<-1){
          // fprintf(stdout,"Call sten daughter with D1:%d D2:%d D3:%d\n",D1,D2,D3);
            bound_sten_daughter(theta, gamma, Arteries, parent,D1,D2,D3,WM,VP);
        }
    else
    {
//        fprintf(stdout,"Call trif with D1:%d D2:%d D3:%d\n",D1,D2,D3);
        bound_trif(theta, gamma, Arteries, parent,D1,D2,D3,WM,VP);
    }
    // Can add quadfurcation here if necessary.
}



void bound_monf (double theta, double gamma, Tube *Arteries[], int parent, int D1, int WM, int VP)
{
    // MJC: Try to define the daughters here
    Tube* PV = Arteries[ parent]; // Parent vessel
    Tube* B1 = Arteries[ D1];      // Branch 1
    int N = PV->N;
    int j = 1;
    int ok = false;
    const int ntrial = 100;
    
    double qA2_term = 0.0;
    if (VP==2) {
        qA2_term=((vel_power+2.0)/(vel_power+1.0)); // Power law
    } else if (VP<2) {

        qA2_term = 1.0; // Stokes boundary layer
    }
    else error ("arteries.cxx","Velocity profile doesn't exist.");
    
    double k1[2];
    double k2[2];
    double k3[2];
    double k4[2];
    
    // These are the flows at the half time step.
      k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
    
      k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);
    
    k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
    k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);
    
    k3[0]   = PV->Qh[N-1]/2.0;
    k3[1]   = B1->Qh[0]/2.0;
    
    k4[0]   = PV->Ah[N-1]/2.0;
    k4[1]   = B1->Ah[0]/2.0;
    
    double xb[12];
    
    // The approximative initial guesses are applied.
    xb[ 0] =  PV->Qh[N-1];                      //Initial guess for Q1_xb n+1
    xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;       //Initial guess for Q1_xb^n+0.5
    xb[ 2] =  PV->Qold[N];                      //Initial guess for Q1_xb+0.5 n+0.5
    xb[ 3] =  B1->Qh[0];                    //Initial guess for Q2_xb n+1
    xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0; //Initial guess for Q2_xb n+0.5
    xb[ 5] =  B1->Qold[0];                  //Initial guess for Q2_xb+0.5 n+0.5
    xb[ 6] =  PV->Ah[N-1];                      //Initial guess for A1_xb n+1
    xb[ 7] = (PV->Aold[N-1] + PV->Aold[N])/2.0;       //Initial guess for A1_xb^n+0.5
    xb[ 8] =  PV->Aold[N];                      //Initial guess for A1_xb+0.5 n+0.5
    xb[ 9] =  B1->Ah[0];                    //Initial guess for A2_xb n+1
    xb[10] = (B1->Aold[0] + B1->Aold[1])/2.0; //Initial guess for A2_xb n+0.5
    xb[11] =  B1->Aold[0];                  //Initial guess for A2_xb+0.5 n+0.5

    
    
    double k7nh  = 0.0;// Bernoulli loss
    double k7n   = 0.0;//

    
    // The residuals (fvec), and the Jacobian is determined, and if possible
    // the system of equations is solved.
    while (j <= ntrial && ok==false) // Find the zero
    {
        double fvec[12];
        
        // The residuals
        
        // Characteristic Q residual at n+1
        fvec[0]  = k1[0]  - xb[0] -
        theta*(qA2_term*sq(xb[2])/xb[8] + PV->Bh(N,xb[8],WM)) +
        gamma*(F(xb[2],xb[8],VP)  + PV->B_RHSh(N,xb[8],WM));
        
        fvec[1]  = k1[1]  - xb[3] +
        theta*(qA2_term*sq(xb[5])/xb[11] + B1->Bh(-1,xb[11],WM)) +
        gamma*(F(xb[5],xb[11],VP)  + B1->B_RHSh(-1,xb[11],WM));
        
        // Characteristic A residual at n+1
        fvec[2]  = - theta*xb[2] - xb[6]  + k2[0];
        fvec[3]  =   theta*xb[5] - xb[9]  + k2[1];
        
        // Flow residuals at n+1/2 (ghost points)
        fvec[4]  = - xb[ 1] + xb[ 2]/2.0  + k3[0];
        fvec[5]  = - xb[ 4] + xb[ 5]/2.0  + k3[1];
        
        // Area residuals at n+1/2 (ghost points)
        fvec[6]  = - xb[ 7] + xb[ 8]/2.0  + k4[0];
        fvec[7]  = - xb[10] + xb[11]/2.0  + k4[1];
        
        // Flow conservation residuals (n+1/2 and n+1)
        fvec[8]  = - xb[ 1] + xb[ 4];
        fvec[9]  = - xb[ 0] + xb[ 3];
        
        // Use these terms if you want to have Benoulli loss: u^2=(q/A)^2
        double sq211 = sq(xb[1]/xb[7]);
        double sq110 = sq(xb[0]/xb[6]);

            fvec[10] =  - PV->P(N,xb[7],WM) + B1->P(0,xb[10],WM) + ab(k7nh)*sq211;
        
            fvec[11] = - PV->P(N,xb[6],WM) + B1->P(0,xb[9],WM) + ab(k7n)*sq110;
        
        double chi[8];
        
        // Here are the residuals for the characteristic matching for flow
        chi[0] = -2.0*theta*qA2_term*xb[ 2]/xb[8] + gamma*dFdQ(xb[8],VP);
        chi[2] =  2.0*theta*qA2_term*xb[ 5]/xb[11] + gamma*dFdQ(xb[11],VP);
        
        // Here are the residuals for the area characteristic matching
        chi[1] = theta*(qA2_term*sq(xb[2]/xb[8]) - PV->dBdAh(N,xb[8],WM)) +
                   gamma*(dFdA(xb[2],xb[8],VP) + PV->dB_RHSdAh(N,xb[8],WM));
        
        chi[3] = theta*( -qA2_term*sq(xb[5]/xb[11]) + B1->dBdAh(-1,xb[11],WM)) +
                  gamma*(dFdA(xb[5],xb[11],VP) + B1->dB_RHSdAh(-1,xb[11],WM));
        
        // Here is pressure conservation (n+1/2)
        chi[4]  = -PV->dPdA(N,xb[7],WM) + sq(xb[1])/cu(xb[7])*(-2.0*ab(k7n)); //Loss term
        chi[5]  = B1->dPdA(0,xb[10],WM);
        
        // Here is pressure conservation (n+1)
        chi[6] = -PV->dPdA(N,xb[6],WM) + sq(xb[0])/cu(xb[6])*(-2.0*ab(k7nh)); //Loss term
        chi[7] = B1->dPdA(0,xb[9],WM);
                                              
        
        for (int row = 0; row < 12; row++)
            for (int col = 0; col < 12; col++)
                fjac[row][col] = 0.0;
        
        // The Jacobian.
        // Order is [row][column]
        
        fjac[ 0][ 0] = -1.0;
        fjac[ 0][ 2] = chi[0];
        fjac[ 0][ 8] = chi[1];
        
        fjac[ 1][ 3] = -1.0;
        fjac[ 1][ 5] = chi[2];
        fjac[ 1][11] = chi[3];
        
        fjac[ 2][ 2] = -theta;
        fjac[ 2][ 6] = -1.0;
        
        fjac[ 3][ 5] =  theta;
        fjac[ 3][ 9] = -1.0;
        
        fjac[ 4][ 1] = -1.0;
        fjac[ 4][ 2] =  0.5;
        
        fjac[ 5][ 4] = -1.0;
        fjac[ 5][ 5] =  0.5;
        
        fjac[ 6][ 7] = -1.0;
        fjac[ 6][ 8] = 0.5;
        
        fjac[ 7][10] = -1.0;
        fjac[ 7][11] = 0.5;
        
        fjac[ 8][ 1] = -1.0;
        fjac[ 8][ 4] =  1.0;
        
        fjac[ 9][ 0] = -1.0;
        fjac[ 9][ 3] =  1.0;
        
        fjac[10][ 7] = chi[4];
        fjac[10][10] = chi[5];
        
        fjac[11][ 6] = chi[6];
        fjac[11][ 9] = chi[7];
        
        
        
        
        // Check whether solution is close enough. If not run the loop again.
        int ch = zero (xb, 12, 1.0e-8, 1.0e-8, fvec, fjac);
        if (ch == 1) ok = true;
//        fprintf(stdout,"AREA DIFFERENCE 2: %lf AT t=%lf  \n",xb[6]-xb[9],t);

        j = j+1;
    }
    
    // Solutions is applied, and right boundary is updated.
    PV->Anew[N] = xb[ 6];
    PV->Qnew[N] = xb[ 0];
    B1->Anew[0] = xb[ 9];
    B1->Qnew[0] = xb[ 3];
    

      if (j >=ntrial) error ("arteries.C","Root not found in the bifurcation");
}

///////////////////////////////////////////////////////////////
/* ADDED BY MJC: Create a solver for a stenosis that uses nearly the same
   details as bound bif, but only solves for a connection between a parent
   and one daughter. For a stenosis, we employ this and add a loss term at
   the inlet of the daughter branch
   11/8/2019
 */
void bound_sten(double theta, double gamma,Tube *Arteries[], int parent, int D1, int WM, int VP)
{
    // MJC: Try to define the daughters here
    Tube* PV = Arteries[ parent]; // Parent vessel
    Tube* B1 = Arteries[ D1];      // Branch 1
    int N = PV->N;
    int j = 1;
    int ok = false;
    const int ntrial = 200;
    
    double qA2_term = 0.0;
    if (VP==2) {
        qA2_term=((vel_power+2.0)/(vel_power+1.0)); // Power law
    } else if (VP<2) {

        qA2_term = 1.0; // Stokes boundary layer
    }
    else error ("arteries.cxx","Velocity profile doesn't exist.");

        
    // Define the stenosis parameters based on derivations
    // For more details, see Young and Tsai, 1973
        double C,Kt,Ku,Ls,sten_factor;
        Kt = 1.52;             // Dissapation due to turbulent forces
        Ku = 1.2;              // Dissapation due to inertial forces
        sten_factor = 0.1;     // Severity of stenoses
        Ls = (PV->L)*0.25;     // Define this outside in future versions
        C = (1.0-sten_factor); // Define degree of stenosis as As/Ap, or 1-sten_factor

        double dp1,dp2,dp3;
        if (C==1.0) {
            dp1=0.0; dp2=0.0; dp3=0.0;
        }
        else
        {
            dp1 = (8.0*Ls*mu*M_PI)/(C); // From Karniadakis group 2019
            dp2 = rho*Kt*sq(1.0/C - 1.0)/2.0;
            dp3 = rho*Ku*Ls;
        }
    
    double k1[2];
    double k2[2];
    double k3[2];
    double k4[2];
    
    k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
    k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);
    
    k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
    k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);
    
    k3[0]   = PV->Qh[N-1]/2.0;
    k3[1]   = B1->Qh[0]/2.0;
    
    k4[0]   = PV->Ah[N-1]/2.0;
    k4[1]   = B1->Ah[0]/2.0;
    
    double xb[12];
    
    // The approximative initial guesses are applied.
    xb[ 0] =  PV->Qh[N-1];                      //Initial guess for Q1_xb n+1
    xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;       //Initial guess for Q1_xb^n+0.5
    xb[ 2] =  PV->Qold[N];                      //Initial guess for Q1_xb+0.5 n+0.5
    xb[ 3] =  B1->Qh[0];                    //Initial guess for Q2_xb n+1
    xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0; //Initial guess for Q2_xb n+0.5
    xb[ 5] =  B1->Qold[0];                  //Initial guess for Q2_xb+0.5 n+0.5
    xb[ 6] =  PV->Ah[N-1];                      //Initial guess for A1_xb n+1
    xb[ 7] = (PV->Aold[N-1] + PV->Aold[N])/2.0;       //Initial guess for A1_xb^n+0.5
    xb[ 8] =  PV->Aold[N];                      //Initial guess for A1_xb+0.5 n+0.5
    xb[ 9] =  B1->Ah[0];                    //Initial guess for A2_xb n+1
    xb[10] = (B1->Aold[0] + B1->Aold[1])/2.0; //Initial guess for A2_xb n+0.5
    xb[11] =  B1->Aold[0];                  //Initial guess for A2_xb+0.5 n+0.5

    double k7nh  = 0.0;//Bernoulli Loss
    double k7n   = 0.0;//

    
    // The residuals (fvec), and the Jacobian is determined, and if possible
    // the system of equations is solved.
    while (j <= ntrial && ok==false) // Find the zero
    {
        double fvec[12];
        
        // Characteristic Q residual at n+1
        fvec[0]  = k1[0]  - xb[0] -
        theta*(qA2_term*sq(xb[2])/xb[8] + PV->Bh(N,xb[8],WM)) +
        gamma*(F(xb[2],xb[8],VP)+PV->B_RHSh(N,xb[8],WM));
        
        fvec[1]  = k1[1]  - xb[3] +
        theta*(qA2_term*sq(xb[5])/xb[11] + B1->Bh(-1,xb[11],WM)) +
        gamma*(F(xb[5],xb[11],VP)  + B1->B_RHSh(-1,xb[11],WM));
        
        // Characteristic A residual at n+1
        fvec[2]  = - theta*xb[2] - xb[6]  + k2[0];
        fvec[3]  =   theta*xb[5] - xb[9]  + k2[1];
        
        // Flow residuals at n+1/2 (ghost points)
        fvec[4]  = - xb[ 1] + xb[ 2]/2.0  + k3[0];
        fvec[5]  = - xb[ 4] + xb[ 5]/2.0  + k3[1];
        
        // Area residuals at n+1/2 (ghost points)
        fvec[6]  = - xb[ 7] + xb[ 8]/2.0  + k4[0];
        fvec[7]  = - xb[10] + xb[11]/2.0  + k4[1];
        
        // Flow conservation residuals (n+1/2 and n+1)
        fvec[8]  = - xb[ 1] + xb[ 4];
        fvec[9]  = - xb[ 0] + xb[ 3];
        
        // Use these terms if you want to have Benoulli loss: u^2=(q/A)^2
        double sq211 = sq(xb[1]/xb[7]);
        double sq110 = sq(xb[0]/xb[6]);
        
        double delta_PA,delta_PB,term1A,term2A,term3A,term1B,term2B,term3B;
        // Define the pressure loss terms. Order is viscous, turbulent, and inertial
            term1A = dp1 * xb[1]/sq(xb[10]);
          
            term2A = (dp2 / sq(xb[10])) * ab(xb[1])*xb[1]; // Using absolute value
          
            term3A = (dp3/xb[10]) * (xb[1]-(PV->Qold[N]))/(2.0*gamma);
          
            delta_PA =  term1A + term2A + term3A;
            
        
            term1B = dp1 * xb[0]/sq(xb[9]);
          
            term2B = (dp2 / sq(xb[9])) * ab(xb[0])*xb[0]; // Using absolute value
          
            term3B = (dp3/xb[9]) * (xb[0]-(PV->Qold[N]))/(2.0*gamma);
          
            delta_PB =  term1B + term2B + term3B;
        
        fvec[10] = - PV->P(N,xb[7],WM) + B1->P(0,xb[10],WM) + ab(k7nh)*sq211 + delta_PA;
        fvec[11] = - PV->P(N,xb[6],WM) + B1->P(0,xb[9],WM)  + ab(k7n)*sq110  + delta_PB;

     
        double chi[10];
        
        // Here are the residuals for the characteristic matching for flow
        chi[0] = -2.0*theta*qA2_term*xb[ 2]/xb[8] + gamma*dFdQ(xb[8],VP);
        chi[2] =  2.0*theta*qA2_term*xb[ 5]/xb[11] + gamma*dFdQ(xb[11],VP);
        
        // Here are the residuals for the area characteristic matching
        chi[1] = theta*(qA2_term*sq(xb[2]/xb[8]) - PV->dBdAh(N,xb[8],WM)) +
                   gamma*(dFdA(xb[2],xb[8],VP) + PV->dB_RHSdAh(N,xb[8],WM));
        
        chi[3] = theta*( -qA2_term*sq(xb[5]/xb[11]) + B1->dBdAh(-1,xb[11],WM)) +
                  gamma*(dFdA(xb[5],xb[11],VP) + B1->dB_RHSdAh(-1,xb[11],WM));
        
        // Here is pressure conservation (n+1/2)
        chi[4]  = -PV->dPdA(N,xb[7],WM) - 2.0*dp1*xb[1]/cu(xb[7])
        - 2.0*dp2*ab(xb[1])*xb[1]/cu(xb[7]) - (dp3/sq(xb[7])) * (xb[1]-(PV->Qold[N]))/(2.0*gamma);
        
        chi[5]  =  B1->dPdA(0,xb[10],WM);
        
        // Here is pressure conservation (n+1)
        chi[6]  = -PV->dPdA(N,xb[6],WM) - 2.0*dp1*xb[0]/cu(xb[6])
        - 2.0*dp2*ab(xb[0])*xb[0]/cu(xb[6]) - (dp3/sq(xb[6])) * (xb[0]-(PV->Qold[N]))/(2.0*gamma);
        
        chi[7] =  B1->dPdA(0,xb[9],WM);
        
        // Additional terms that arrive because of the pressure loss (depends on parent flow)
        chi[8]  = dp1/sq(xb[7]) + (dp2/sq(xb[7]))*2.0*xb[1]/ab(xb[1])
        + dp3/(2.0*gamma*xb[7]);
        
        chi[9] = dp1/sq(xb[6]) + (dp2/sq(xb[6]))*2.0*xb[0]/ab(xb[0])
        + dp3/(2.0*gamma*xb[6]);
        
                                                   
        
        for (int row = 0; row < 12; row++)
            for (int col = 0; col < 12; col++)
                fjac[row][col] = 0.0;
        
        // The Jacobian.
        // MJC: Ordered by column
        // Column 0:
        fjac[ 0][ 0] = -1.0;
        fjac[ 0][ 2] = chi[0];
        fjac[ 0][ 8] = chi[1];
        
        
        fjac[ 1][ 3] = -1.0;
        fjac[ 1][ 5] = chi[2];
        fjac[ 1][11] = chi[3];
        
        
        fjac[ 2][2] = -theta;
        fjac[ 2][6] = -1.0;
        
        fjac[ 3][5] =  theta;
        fjac[ 3][9] = -1.0;
        
        
        fjac[ 4][1] = -1.0;
        fjac[ 4][2] =  0.5;
        
        
        fjac[ 5][4] = -1.0;
        fjac[ 5][5] =  0.5;
        
        
        fjac[ 6][ 7] = -1.0;
        fjac[ 6][ 8] = 0.5;
        
        
        fjac[7][10] = -1.0;
        fjac[7][11] = 0.5;
        
        
        fjac[8][1] = -1.0;
        fjac[8][4] =  1.0;
        
        
        fjac[9][0]  = -1.0;
        fjac[ 9][3] =  1.0;
        
        fjac[10][ 1] = chi[8]; //Loss
        fjac[10][ 7] = chi[4];
        fjac[10][10] = chi[5];
        
        fjac[11][ 0] = chi[9]; //Loss
        fjac[11][ 6] = chi[6]; //
        fjac[11][ 9] = chi[7]; //

        // Check whether solution is close enough. If not run the loop again.
        // int ch = zero (xb, 18, 1.0e-4, 1.0e-4, fvec, fjac);
        int ch = zero (xb, 12, 1.0e-8, 1.0e-8, fvec, fjac);
        if (ch == 1) ok = true;
//        fprintf(stdout,"AREA DIFFERENCE 2: %lf AT t=%lf  \n",xb[6]-xb[9],t);

        j = j+1;
    }
    
    // Solutions is applied, and right boundary is updated.
    PV->Anew[N] = xb[ 6];
    PV->Qnew[N] = xb[ 0];
    B1->Anew[0] = xb[ 9];
    B1->Qnew[0] = xb[ 3];

      if (j >=ntrial) error ("arteries.C","Root not found in the bifurcation");
}





void bound_bif (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int WM, int VP)
{
    Tube* PV  = Arteries[ parent];
    Tube* B1  = Arteries[ D1];
    Tube* B2  = Arteries[ D2];
    
    double qA2_term = 0.0;
    if (VP==2) {
        qA2_term=((vel_power+2.0)/(vel_power+1.0)); // Power law
    } else if (VP<2) {

        qA2_term = 1.0; // Stokes boundary layer
    }
    else error ("arteries.cxx","Velocity profile doesn't exist.");
    
    
  int N = PV->N;
  double PN;
  int j = 1;
  int ok = false;
  const int ntrial = 500;
    

    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];
    
  k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
  k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);
  k1[2]= B2->Qold[0] - theta*(B2->R2h[0]) + gamma*(B2->S2h[0]);
    

  k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
  k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);
  k2[2] = B2->Aold[0] - theta*(B2->R1h[0]);

  k3[0] = PV->Qh[N-1]/2.0;
  k3[1] = B1->Qh[0]/2.0;
  k3[2] = B2->Qh[0]/2.0;

  k4[0]  = PV->Ah[N-1]/2.0;
  k4[1]  = B1->Ah[0]/2.0;
  k4[2]  = B2->Ah[0]/2.0;

  double xb[18];

  // The approximative initial guesses are applied.
  xb[ 0] =  PV->Qh[N-1];                        //Initial guess for Q1_xb n+1
  xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;       //Initial guess for Q1_xb^n+0.5
  xb[ 2] =  PV->Qold[N];                        //Initial guess for Q1_xb+0.5 n+0.5
  xb[ 3] =  B1->Qh[0];                      //Initial guess for Q2_xb n+1
  xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0; //Initial guess for Q2_xb n+0.5
  xb[ 5] =  B1->Qold[0];                    //Initial guess for Q2_xb+0.5 n+0.5
  xb[ 6] =  B2->Qh[0];                      //Initial guess for Q3_xb n+1
  xb[ 7] = (B2->Qold[0] + B2->Qold[1])/2.0; //Initial guess for Q3_xb n+0.5
  xb[ 8] =  B2->Qold[0];                    //Initial guess for Q3_xb+0.5 n+0.5
  xb[ 9] =  PV->Ah[N-1];                        //Initial guess for A1_xb n+1
  xb[10] = (PV->Aold[N-1] + PV->Aold[N])/2.0;       //Initial guess for A1_xb^n+0.5
  xb[11] =  PV->Aold[N];                        //Initial guess for A1_xb+0.5 n+0.5
  xb[12] =  B1->Ah[0];                      //Initial guess for A2_xb n+1
  xb[13] = (B1->Aold[0] + B1->Aold[1])/2.0; //Initial guess for A2_xb n+0.5
  xb[14] =  B1->Aold[0];                    //Initial guess for A2_xb+0.5 n+0.5
  xb[15] =  B2->Ah[0];                      //Initial guess for A3_xb n+1
  xb[16] = (B2->Aold[0] + B2->Aold[1])/2.0; //Initial guess for A3_xb n+0.5
  xb[17] =  B2->Aold[0];                    //Initial guess for A3_xb+0.5 n+0.5

  double k7nh  = 0; //Bernoulli loss term
  double k7n   = 0;
  double k7anh = 0;
  double k7an  = 0;

  // The residuals (fvec), and the Jacobian is determined, and if possible
  // the system of equations is solved.
  while (j <= ntrial && ok==false) // Find the zero
  {
    double fvec[18];
    // The residuals.
      
      // Characteristic Q residual at n+1
    fvec[0]  = k1[0]  - xb[0] -
    theta*(qA2_term*sq(xb[2])/xb[11] + PV->Bh(N,xb[11],WM)) +
    gamma*(F(xb[2],xb[11],VP)+PV->B_RHSh(N,xb[11],WM));

    fvec[1]  = k1[1]  - xb[3] +
                 theta*(qA2_term*sq(xb[5])/xb[14] + B1->Bh(-1,xb[14],WM)) +
                 gamma*(F(xb[5],xb[14],VP)  + B1->B_RHSh(-1,xb[14],WM));

    fvec[2]  = k1[2] - xb[6] +
    theta*(qA2_term*sq(xb[8])/xb[17] + B2->Bh(-1,xb[17],WM)) +
    gamma*(F(xb[8],xb[17],VP)  + B2->B_RHSh(-1,xb[17],WM));

    // Characteristic A residual at n+1
    fvec[3]  = - theta*xb[2] - xb[ 9] + k2[0];
    fvec[4]  =   theta*xb[5] - xb[12] + k2[1];
    fvec[5]  =   theta*xb[8] - xb[15] + k2[2];
      
    // Flow residuals at n+1/2 (ghost points)
    fvec[6]  = - xb[ 1] + xb[ 2]/2.0 + k3[0];
    fvec[7]  = - xb[ 4] + xb[ 5]/2.0 + k3[1];
    fvec[8]  = - xb[ 7] + xb[ 8]/2.0 + k3[2];
      
    // Area residuals at n+1/2 (ghost points)
    fvec[9]  = - xb[10] + xb[11]/2.0 + k4[0];
    fvec[10] = - xb[13] + xb[14]/2.0 + k4[1];
    fvec[11] = - xb[16] + xb[17]/2.0 + k4[2];
      
    // Flow conservation residuals (n+1/2 and n+1)
    fvec[12] = - xb[ 1] + xb[ 4] + xb[7];
    fvec[13] = - xb[ 0] + xb[ 3] + xb[6];

    // Pressure continuity at the n+1/2 time step
    PN    = PV->P(N,xb[10],WM);
    double u_n_half = sq(xb[1]/xb[10]);

    // The if statements here only matter if a minor loss
    // is included
      fvec[14] =  - PN + B1->P(0,xb[13],WM) + ab(k7nh)*u_n_half;
      fvec[15] =  - PN + B2->P(0,xb[16],WM) + ab(k7anh)*u_n_half;

    // Pressure continuity at the n+1 time step
    PN    = PV->P(N,xb[9],WM);
    double u_n_1 = sq(xb[0]/xb[9]);
      fvec[16] = - PN + B1->P(0,xb[12],WM) + ab(k7n)*u_n_1;
      fvec[17] = - PN + B2->P(0,xb[15],WM) + ab(k7an)*u_n_1;

    for (int row = 0; row < 18; row++)
      for (int col = 0; col < 18; col++)
        fjac[row][col] = 0.0;

      
      
    // The Jacobian.
      
      double chi[12];
      
       // Here are the residuals for the characteristic matching for flow
      chi[0] = -2.0*theta*qA2_term*xb[ 2]/xb[11] + gamma*dFdQ(xb[11],VP);
      chi[2] =  2.0*theta*qA2_term*xb[ 5]/xb[14] + gamma*dFdQ(xb[14],VP);
      chi[4] =  2.0*theta*qA2_term*xb[ 8]/xb[17] + gamma*dFdQ(xb[17],VP);
      
      // Here are the residuals for the area characteristic matching
      chi[1] = theta*(qA2_term*sq(xb[2]/xb[11]) - PV->dBdAh(N,xb[11],WM)) +
                 gamma*(dFdA(xb[2],xb[11],VP) + PV->dB_RHSdAh(N,xb[11],WM));
      
      chi[3] = theta*( -qA2_term*sq(xb[5]/xb[14]) + B1->dBdAh(-1,xb[14],WM)) +
                gamma*(dFdA(xb[5],xb[14],VP) + B1->dB_RHSdAh(-1,xb[14],WM));
      
      chi[5] = theta*( -qA2_term*sq(xb[8]/xb[17]) + B2->dBdAh(-1,xb[17],WM)) +
                gamma*(dFdA(xb[8],xb[17],VP) + B2->dB_RHSdAh(-1,xb[17],WM));
      
      // Here is pressure conservation
      chi[6]  = -PV->dPdA(N,xb[10],WM) + sq(xb[1])/cu(xb[10])*(-2.0*ab(k7nh)); //Loss term
      chi[7]  = B1->dPdA(0,xb[13],WM);
      chi[8] = B2->dPdA(0,xb[16],WM);
      
      chi[9] = -PV->dPdA(N,xb[9],WM) + sq(xb[1])/cu(xb[10])*(-2.0*ab(k7nh)); //Loss term
      chi[10] = B1->dPdA(0,xb[12],WM);
      chi[11] = B2->dPdA(0,xb[15],WM);

      
      //NEW JACOBIAN
            // Order is [row][column]
            fjac[ 0][ 0]  = -1.0;
            fjac[ 0][ 2] = chi[0];
            fjac[ 0][11] = chi[1];
            
            fjac[ 1][ 3] = -1.0;
            fjac[ 1][ 5] = chi[2];
            fjac[ 1][14] = chi[3];

            fjac[ 2][ 6] = -1.0;
            fjac[ 2][ 8] = chi[4];
            fjac[ 2][17] = chi[5];

            fjac[ 3][ 2] = -theta;
            fjac[ 3][ 9] = -1.0;
            
            fjac[ 4][ 5] = theta;
            fjac[ 4][12] = -1.0;
            
            fjac[ 5][ 8] = theta;
            fjac[ 5][15] = -1.0;
            
            fjac[ 6][ 1] = -1.0;
            fjac[ 6][ 2] = 0.5;
            
            fjac[ 7][ 4] = -1.0;
            fjac[ 7][ 5] = 0.5;
            
            fjac[ 8][ 7] = -1.0;
            fjac[ 8][ 8] = 0.5;
            
            fjac[ 9][ 10] = -1.0;
            fjac[ 9][ 11] = 0.5;
            
            fjac[10][13] = -1.0;
            fjac[10][14] = 0.5;
            
            fjac[11][16] = -1.0;
            fjac[11][17] = 0.5;
            
            fjac[12][ 1] = -1.0;
            fjac[12][ 4] = 1.0;
            fjac[12][ 7] = 1.0;
            
            fjac[13][ 0] = -1.0;
            fjac[13][ 3] = 1.0;
            fjac[13][ 6] = 1.0;
            
      //      fjac[14][ 1] = ab(k7n)*(xb[1]/xb[10]); // Bernoulli Loss
            fjac[14][10] = chi[6];
            fjac[14][13] = chi[7];
                                    
                                    
      //      fjac[15][ 1] = ab(k7n)*(xb[1]/xb[10]); // Bernoulli Loss
            fjac[15][10] = chi[6];
            fjac[15][16] = chi[8];
            
      //      fjac[16][ 0] = ab(k7n)*(xb[0]/xb[9]); // Bernoulli Loss
            fjac[16][ 9] = chi[9];
            fjac[16][12] = chi[10];
            
            
      //      fjac[17][ 0] = ab(k7n)*(xb[0]/xb[9]); // Bernoulli Loss
            fjac[17][ 9] = chi[9];
            fjac[17][15] = chi[11];
      
    // Check whether solution is close enough. If not run the loop again.
    int ch = zero (xb, 18, 1.0e-8, 1.0e-8, fvec, fjac);
    if (ch == 1) ok = true;

    j = j+1;
  }

  // Solutions is applied, and right boundary is updated.
  PV->Anew[N] = xb[ 9];
  PV->Qnew[N] = xb[ 0];
  B1->Anew[0] = xb[12];
  B1->Qnew[0] = xb[ 3];
  B2->Anew[0] = xb[15];
  B2->Qnew[0] = xb[ 6];

    if (j >=ntrial) {error ("arteries.C","Root not found in the bifurcation");
        exit(1);}
}


void bound_sten_bif (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int WM, int VP)
{
    Tube* PV   = Arteries[ parent];
    Tube* B1  = Arteries[ D1];
    Tube* B2  = Arteries[ D2];
    
    
  int N = PV->N;
  int j = 1;
  int ok = false;
  const int ntrial = 500000;
    
    double qA2_term = 0.0;
    if (VP==2) {
        qA2_term=((vel_power+2.0)/(vel_power+1.0)); // Power law
    } else if (VP<2) {

        qA2_term = 1.0; // Stokes boundary layer
    }
    else error ("arteries.cxx","Velocity profile doesn't exist.");
    
    // Define the stenosis parameters based on derivations
    // For more details, see Young and Tsai, 1973
        double C,Kt,Ku,Ls,sten_factor;
        Kt = 1.52;             // Dissapation due to turbulent forces
        Ku = 1.2;              // Dissapation due to inertial forces
        sten_factor = (PV->sten_factor);     // Severity of stenoses
        Ls = PV->sten_length;     // Define this outside in future versions
        C = (1.0-sten_factor); // Define degree of stenosis as As/Ap, or 1-sten_factor

        double dp1,dp2,dp3;
        if (C==1.0) {
            dp1=0.0; dp2=0.0; dp3=0.0;
        }
        else
        {
            dp1 = (8.0*Ls*mu*M_PI)/(C); // From Karniadakis group 2019
            dp2 = rho*Kt*sq(1.0/C - 1.0)/2.0;
            dp3 = rho*Ku*Ls;
        }
    

    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];
    
  k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
  k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);
  k1[2] = B2->Qold[0] - theta*(B2->R2h[0]) + gamma*(B2->S2h[0]);
    

  k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
  k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);
  k2[2] = B2->Aold[0] - theta*(B2->R1h[0]);

  k3[0] = PV->Qh[N-1]/2.0;
  k3[1] = B1->Qh[0]/2.0;
  k3[2] = B2->Qh[0]/2.0;

  k4[0]  = PV->Ah[N-1]/2.0;
  k4[1]  = B1->Ah[0]/2.0;
  k4[2]  = B2->Ah[0]/2.0;

  double xb[18];

  // The approximative initial guesses are applied.
  xb[ 0] =  PV->Qh[N-1];                        //Initial guess for Q1_xb n+1
  xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;       //Initial guess for Q1_xb^n+0.5
  xb[ 2] =  PV->Qold[N];                        //Initial guess for Q1_xb+0.5 n+0.5
  xb[ 3] =  B1->Qh[0];                      //Initial guess for Q2_xb n+1
  xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0; //Initial guess for Q2_xb n+0.5
  xb[ 5] =  B1->Qold[0];                    //Initial guess for Q2_xb+0.5 n+0.5
  xb[ 6] =  B2->Qh[0];                      //Initial guess for Q3_xb n+1
  xb[ 7] = (B2->Qold[0] + B2->Qold[1])/2.0; //Initial guess for Q3_xb n+0.5
  xb[ 8] =  B2->Qold[0];                    //Initial guess for Q3_xb+0.5 n+0.5
  xb[ 9] =  PV->Ah[N-1];                        //Initial guess for A1_xb n+1
  xb[10] = (PV->Aold[N-1] + PV->Aold[N])/2.0;       //Initial guess for A1_xb^n+0.5
  xb[11] =  PV->Aold[N];                        //Initial guess for A1_xb+0.5 n+0.5
  xb[12] =  B1->Ah[0];                      //Initial guess for A2_xb n+1
  xb[13] = (B1->Aold[0] + B1->Aold[1])/2.0; //Initial guess for A2_xb n+0.5
  xb[14] =  B1->Aold[0];                    //Initial guess for A2_xb+0.5 n+0.5
  xb[15] =  B2->Ah[0];                      //Initial guess for A3_xb n+1
  xb[16] = (B2->Aold[0] + B2->Aold[1])/2.0; //Initial guess for A3_xb n+0.5
  xb[17] =  B2->Aold[0];                    //Initial guess for A3_xb+0.5 n+0.5


  // The residuals (fvec), and the Jacobian is determined, and if possible
  // the system of equations is solved.
  while (j <= ntrial && ok==false) // Find the zero
  {
    double fvec[18];
    // The residuals.
      
      // Characteristic Q residual at n+1
    fvec[0]  = k1[0]  - xb[0] -
    theta*(qA2_term*sq(xb[2])/xb[11] + PV->Bh(N,xb[11],WM)) +
    gamma*(F(xb[2],xb[11],VP)+PV->B_RHSh(N,xb[11],WM));

    fvec[1]  = k1[1]  - xb[3] +
                 theta*(qA2_term*sq(xb[5])/xb[14] + B1->Bh(-1,xb[14],WM)) +
                 gamma*(F(xb[5],xb[14],VP)  + B1->B_RHSh(-1,xb[14],WM));

    fvec[2]  = k1[2] - xb[6] +
    theta*(qA2_term*sq(xb[8])/xb[17] + B2->Bh(-1,xb[17],WM)) +
    gamma*(F(xb[8],xb[17],VP)  + B2->B_RHSh(-1,xb[17],WM));

    // Characteristic A residual at n+1
    fvec[3]  = - theta*xb[2] - xb[ 9] + k2[0];
    fvec[4]  =   theta*xb[5] - xb[12] + k2[1];
    fvec[5]  =   theta*xb[8] - xb[15] + k2[2];
      
    // Flow residuals at n+1/2 (ghost points)
    fvec[6]  = - xb[ 1] + xb[ 2]/2.0 + k3[0];
    fvec[7]  = - xb[ 4] + xb[ 5]/2.0 + k3[1];
    fvec[8]  = - xb[ 7] + xb[ 8]/2.0 + k3[2];
      
    // Area residuals at n+1/2 (ghost points)
    fvec[9]  = - xb[10] + xb[11]/2.0 + k4[0];
    fvec[10] = - xb[13] + xb[14]/2.0 + k4[1];
    fvec[11] = - xb[16] + xb[17]/2.0 + k4[2];
      
    // Flow conservation residuals (n+1/2 and n+1)
    fvec[12] = - xb[ 1] + xb[ 4] + xb[7];
    fvec[13] = - xb[ 0] + xb[ 3] + xb[6];

      double delta_PA,delta_PB,term1A,term2A,term3A,term1B,term2B,term3B;
    // Define the pressure loss terms. Order is viscous, turbulent, and inertial
        term1A = dp1 * xb[1]/sq(xb[10]);
      
        term2A = (dp2 / sq(xb[10])) * ab(xb[1])*xb[1]; // Using absolute value
      
        term3A = (dp3/xb[10]) * (xb[1]-(PV->Qold[N]))/(2.0*gamma);
      
        delta_PA =  term1A + term2A + term3A;
        
    
        term1B = dp1 * xb[0]/sq(xb[9]);
      
        term2B = (dp2 / sq(xb[9])) * ab(xb[0])*xb[0]; // Using absolute value
      
        term3B = (dp3/xb[9]) * (xb[0]-(PV->Qold[N]))/(2.0*gamma);
      
        delta_PB =  term1B + term2B + term3B;

    // Pressure loss at the n+1/2 time step
      fvec[14] = - PV->P(N,xb[ 10],WM) + B1->P(0,xb[13],WM) + delta_PA;
      fvec[15] = - PV->P(N,xb[ 10],WM) + B2->P(0,xb[16],WM) + delta_PA;

    // Pressure continuity at the n+1 time step
    fvec[16] = - PV->P(N,xb[ 9],WM) + B1->P(0,xb[12],WM) + delta_PB;
    fvec[17] = - PV->P(N,xb[ 9],WM) + B2->P(0,xb[15],WM) + delta_PB;

    for (int row = 0; row < 18; row++)
      for (int col = 0; col < 18; col++)
        fjac[row][col] = 0.0;

      
      
    // The Jacobian.
      
      double chi[14];
      
       // Here are the residuals for the characteristic matching for flow
      chi[0] = -2.0*theta*qA2_term*xb[ 2]/xb[11] + gamma*dFdQ(xb[11],VP);
      chi[2] =  2.0*theta*qA2_term*xb[ 5]/xb[14] + gamma*dFdQ(xb[14],VP);
      chi[4] =  2.0*theta*qA2_term*xb[ 8]/xb[17] + gamma*dFdQ(xb[17],VP);
      
      // Here are the residuals for the area characteristic matching
      chi[1] = theta*(qA2_term*sq(xb[2]/xb[11]) - PV->dBdAh(N,xb[11],WM)) +
                 gamma*(dFdA(xb[2],xb[11],VP) + PV->dB_RHSdAh(N,xb[11],WM));
      
      chi[3] = theta*( -qA2_term*sq(xb[5]/xb[14]) + B1->dBdAh(-1,xb[14],WM)) +
                gamma*(dFdA(xb[5],xb[14],VP) + B1->dB_RHSdAh(-1,xb[14],WM));
      
      chi[5] = theta*( -qA2_term*sq(xb[8]/xb[17]) + B2->dBdAh(-1,xb[17],WM)) +
                gamma*(dFdA(xb[8],xb[17],VP) + B2->dB_RHSdAh(-1,xb[17],WM));
      
      // Here are the pressure loss terms at n+1/2
      // Derivative with respect to flow
      chi[6]  = dp1/sq(xb[10]) + (dp2/sq(xb[10]))*2.0*xb[1]/ab(xb[1])
                           + dp3/(2.0*gamma*xb[10]);
      
      // Derivative with respect to area
      chi[7]  = -PV->dPdA(N,xb[10],WM) - 2.0*dp1*xb[1]/cu(xb[10])
      - 2.0*dp2*ab(xb[1])*xb[1]/cu(xb[10]) - (dp3/sq(xb[10])) * (xb[1]-(PV->Qold[N]))/(2.0*gamma);
      
      chi[8]  = B1->dPdA(0,xb[13],WM);
      
      chi[9]  = B2->dPdA(0,xb[16],WM);
      
      
      // Here are the pressure loss terms at n+1
      // Derivative with respect to flow
      chi[10]  = dp1/sq(xb[9]) + (dp2/sq(xb[9]))*2.0*xb[0]/ab(xb[0])
                           + dp3/(2*gamma*xb[9]);
      
      // Derivative with respect to area
      chi[11]  = -PV->dPdA(N,xb[9],WM) - 2.0*dp1*xb[0]/cu(xb[9])
      - 2.0*dp2*ab(xb[0])*xb[0]/cu(xb[9]) - (dp3/sq(xb[9])) * (xb[0]-(PV->Qold[N]))/(2.0*gamma);
      
      chi[12] = B1->dPdA(0,xb[12],WM);
      
      chi[13] = B2->dPdA(0,xb[15],WM);

      
      //NEW JACOBIAN
            // Order is [row][column]
            fjac[ 0][ 0]  = -1.0;
            fjac[ 0][ 2] = chi[0];
            fjac[ 0][11] = chi[1];
            
            fjac[ 1][ 3] = -1.0;
            fjac[ 1][ 5] = chi[2];
            fjac[ 1][14] = chi[3];

            fjac[ 2][ 6] = -1.0;
            fjac[ 2][ 8] = chi[4];
            fjac[ 2][17] = chi[5];

            fjac[ 3][ 2] = -theta;
            fjac[ 3][ 9] = -1.0;
            
            fjac[ 4][ 5] = theta;
            fjac[ 4][12] = -1.0;
            
            fjac[ 5][ 8] = theta;
            fjac[ 5][15] = -1.0;
            
            fjac[ 6][ 1] = -1.0;
            fjac[ 6][ 2] = 0.5;
            
            fjac[ 7][ 4] = -1.0;
            fjac[ 7][ 5] = 0.5;
            
            fjac[ 8][ 7] = -1.0;
            fjac[ 8][ 8] = 0.5;
            
            fjac[ 9][ 10] = -1.0;
            fjac[ 9][ 11] = 0.5;
            
            fjac[10][13] = -1.0;
            fjac[10][14] = 0.5;
            
            fjac[11][16] = -1.0;
            fjac[11][17] = 0.5;
            
            fjac[12][ 1] = -1.0;
            fjac[12][ 4] = 1.0;
            fjac[12][ 7] = 1.0;
            
            fjac[13][ 0] = -1.0;
            fjac[13][ 3] = 1.0;
            fjac[13][ 6] = 1.0;
            
            fjac[14][ 1] = chi[6];
            fjac[14][10] = chi[7];
            fjac[14][13] = chi[8];
                                    
                                    
            fjac[15][ 1] = chi[6];
            fjac[15][10] = chi[7];
            fjac[15][16] = chi[9];
            
            fjac[16][ 0] = chi[10];
            fjac[16][ 9] = chi[11];
            fjac[16][12] = chi[12];
            
            
            fjac[17][ 0] = chi[10];
            fjac[17][ 9] = chi[11];
            fjac[17][15] = chi[13];
      

    // Check whether solution is close enough. If not run the loop again.
    int ch = zero (xb, 18, 1.0e-8, 1.0e-8, fvec, fjac);
    if (ch == 1) ok = true;

    j = j+1;
  }

  // Solutions is applied, and right boundary is updated.
  PV->Anew[N] = xb[ 9];
  PV->Qnew[N] = xb[ 0];
  B1->Anew[0] = xb[12];
  B1->Qnew[0] = xb[ 3];
  B2->Anew[0] = xb[15];
  B2->Qnew[0] = xb[ 6];
    
//    fprintf(stdout,"Daughter 1 Qold: %d Daughter 2: %d\n",*daughters,*(daughters+1));
//    fprintf(stdout,"IN BIF: Parent Qold: %lf Daughter 1 Qold: %lf Daughter 2 Qold: %lf\n",Qold[N],B1->Qold[0],B2->Qold[0]);

    if (j >=ntrial) {error ("arteries.C","Root not found in the bifurcation");
        exit(1);}
}
void bound_sten_daughter (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int D3, int WM, int VP)
{
    Tube* PV   = Arteries[ parent];
    Tube* B1  = Arteries[ D1];
    Tube* B2  = Arteries[ D2];

    int which_daughter = -D3;
    double d1_coeff, d2_coeff;

  int N = PV->N;
  int j = 1;
  int ok = false;
  const int ntrial = 500000;

    // fprintf(stdout,"Q at parent %d: %lf %lf\n",parent,PV->Qold[0],PV->Qold[N]);

    double qA2_term = 0.0;
    if (VP==2) {
        qA2_term=((vel_power+2.0)/(vel_power+1.0)); // Power law
    } else if (VP<2) {

        qA2_term = 1.0; // Stokes boundary layer
    }
    else error ("arteries.cxx","Velocity profile doesn't exist.");
    
    // Define the stenosis parameters based on derivations
    // For more details, see Young and Tsai, 1973
        double C,Kt,Ku,Ls,sten_factor;
        Kt = 1.52;             // Dissapation due to turbulent forces
        Ku = 1.2;              // Dissapation due to inertial forces
        if (which_daughter==2){
          d1_coeff=1.0; // Only apply stenosis conditions to first branch
          d2_coeff=0.0;
          sten_factor = (B1->sten_factor);     // Severity of stenoses
          Ls = B1->sten_length;     // Define this outside in future versions
        }
        else if (which_daughter==3){
          d1_coeff=0.0; // Only apply stenosis conditions to second branch
          d2_coeff=1.0;
          sten_factor = (B2->sten_factor);     // Severity of stenoses
          Ls = B2->sten_length;     // Define this outside in future versions
        }
        else{
          error("Junction.cxx","Incompatible stenosis");
          exit(1);
        }
        
        C = (1.0-sten_factor); // Define degree of stenosis as As/Ap, or 1-sten_factor

        double dp1,dp2,dp3;
        if (C==1.0) {//No stenosis
            dp1=0.0; dp2=0.0; dp3=0.0;
        }
        else
        {
            dp1 = (8.0*Ls*mu*M_PI)/(C); // From Karniadakis group 2019
            dp2 = rho*Kt*sq(1.0/C - 1.0)/2.0;
            dp3 = rho*Ku*Ls;
        }

    

    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];
    
  k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
  k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);
  k1[2] = B2->Qold[0] - theta*(B2->R2h[0]) + gamma*(B2->S2h[0]);
    

  k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
  k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);
  k2[2] = B2->Aold[0] - theta*(B2->R1h[0]);

  k3[0] = PV->Qh[N-1]/2.0;
  k3[1] = B1->Qh[0]/2.0;
  k3[2] = B2->Qh[0]/2.0;

  k4[0]  = PV->Ah[N-1]/2.0;
  k4[1]  = B1->Ah[0]/2.0;
  k4[2]  = B2->Ah[0]/2.0;

  double xb[18];

  // The approximative initial guesses are applied.
  xb[ 0] =  PV->Qh[N-1];                        //Initial guess for Q1_xb n+1
  xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;       //Initial guess for Q1_xb^n+0.5
  xb[ 2] =  PV->Qold[N];                        //Initial guess for Q1_xb+0.5 n+0.5
  xb[ 3] =  B1->Qh[0];                      //Initial guess for Q2_xb n+1
  xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0; //Initial guess for Q2_xb n+0.5
  xb[ 5] =  B1->Qold[0];                    //Initial guess for Q2_xb+0.5 n+0.5
  xb[ 6] =  B2->Qh[0];                      //Initial guess for Q3_xb n+1
  xb[ 7] = (B2->Qold[0] + B2->Qold[1])/2.0; //Initial guess for Q3_xb n+0.5
  xb[ 8] =  B2->Qold[0];                    //Initial guess for Q3_xb+0.5 n+0.5
  xb[ 9] =  PV->Ah[N-1];                        //Initial guess for A1_xb n+1
  xb[10] = (PV->Aold[N-1] + PV->Aold[N])/2.0;       //Initial guess for A1_xb^n+0.5
  xb[11] =  PV->Aold[N];                        //Initial guess for A1_xb+0.5 n+0.5
  xb[12] =  B1->Ah[0];                      //Initial guess for A2_xb n+1
  xb[13] = (B1->Aold[0] + B1->Aold[1])/2.0; //Initial guess for A2_xb n+0.5
  xb[14] =  B1->Aold[0];                    //Initial guess for A2_xb+0.5 n+0.5
  xb[15] =  B2->Ah[0];                      //Initial guess for A3_xb n+1
  xb[16] = (B2->Aold[0] + B2->Aold[1])/2.0; //Initial guess for A3_xb n+0.5
  xb[17] =  B2->Aold[0];                    //Initial guess for A3_xb+0.5 n+0.5


  // The residuals (fvec), and the Jacobian is determined, and if possible
  // the system of equations is solved.
  while (j <= ntrial && ok==false) // Find the zero
  {
    double fvec[18];
    // The residuals.
      
      // Characteristic Q residual at n+1
    fvec[0]  = k1[0]  - xb[0] -
    theta*(qA2_term*sq(xb[2])/xb[11] + PV->Bh(N,xb[11],WM)) +
    gamma*(F(xb[2],xb[11],VP)+PV->B_RHSh(N,xb[11],WM));

    fvec[1]  = k1[1]  - xb[3] +
                 theta*(qA2_term*sq(xb[5])/xb[14] + B1->Bh(-1,xb[14],WM)) +
                 gamma*(F(xb[5],xb[14],VP)  + B1->B_RHSh(-1,xb[14],WM));

    fvec[2]  = k1[2] - xb[6] +
    theta*(qA2_term*sq(xb[8])/xb[17] + B2->Bh(-1,xb[17],WM)) +
    gamma*(F(xb[8],xb[17],VP)  + B2->B_RHSh(-1,xb[17],WM));

    // Characteristic A residual at n+1
    fvec[3]  = - theta*xb[2] - xb[ 9] + k2[0];
    fvec[4]  =   theta*xb[5] - xb[12] + k2[1];
    fvec[5]  =   theta*xb[8] - xb[15] + k2[2];
      
    // Flow residuals at n+1/2 (ghost points)
    fvec[6]  = - xb[ 1] + xb[ 2]/2.0 + k3[0];
    fvec[7]  = - xb[ 4] + xb[ 5]/2.0 + k3[1];
    fvec[8]  = - xb[ 7] + xb[ 8]/2.0 + k3[2];
      
    // Area residuals at n+1/2 (ghost points)
    fvec[9]  = - xb[10] + xb[11]/2.0 + k4[0];
    fvec[10] = - xb[13] + xb[14]/2.0 + k4[1];
    fvec[11] = - xb[16] + xb[17]/2.0 + k4[2];
      
    // Flow conservation residuals (n+1/2 and n+1)
    fvec[12] = - xb[ 1] + xb[ 4] + xb[7];
    fvec[13] = - xb[ 0] + xb[ 3] + xb[6];

      double delta_P1_B1, delta_P2_B1, delta_P1_B2, delta_P2_B2;
      double term1A,term2A,term3A,term1B,term2B,term3B;

    // Define the pressure loss terms. Order is viscous, turbulent, and inertial

    // Branch 1
        term1A = dp1 * xb[4]/sq(xb[13]);
        term2A = (dp2 / sq(xb[13])) * ab(xb[4])*xb[4]; // Using absolute value
        term3A = (dp3/xb[13]) * (xb[4]-(B1->Qold[0]))/(2.0*gamma);
        delta_P1_B1 =  d1_coeff*(term1A + term2A + term3A);

        term1B = dp1 * xb[3]/sq(xb[12]);
        term2B = (dp2 / sq(xb[12])) * ab(xb[3])*xb[3]; // Using absolute value
        term3B = (dp3/xb[12]) * (xb[3]-(B1->Qold[0]))/(2.0*gamma);
        delta_P2_B1 =  d1_coeff*(term1B + term2B + term3B);

        // Branch 2
        term1A = dp1 * xb[7]/sq(xb[16]);
        term2A = (dp2 / sq(xb[16])) * ab(xb[7])*xb[7]; // Using absolute value
        term3A = (dp3/xb[16]) * (xb[7]-(B2->Qold[0]))/(2.0*gamma);
        delta_P1_B2 =  d2_coeff*(term1A + term2A + term3A);

        term1B = dp1 * xb[6]/sq(xb[15]);
        term2B = (dp2 / sq(xb[15])) * ab(xb[6])*xb[6]; // Using absolute value
        term3B = (dp3/xb[15]) * (xb[6]-(B2->Qold[0]))/(2.0*gamma);
        delta_P2_B2 =  d2_coeff*(term1B + term2B + term3B);

    // Pressure loss at the n+1/2 time step
      fvec[14] = - PV->P(N,xb[ 10],WM) + B1->P(0,xb[13],WM) + delta_P1_B1;
      fvec[15] = - PV->P(N,xb[ 10],WM) + B2->P(0,xb[16],WM) + delta_P1_B2;

    // Pressure continuity at the n+1 time step
    fvec[16] = - PV->P(N,xb[ 9],WM) + B1->P(0,xb[12],WM) + delta_P2_B1;
    fvec[17] = - PV->P(N,xb[ 9],WM) + B2->P(0,xb[15],WM) + delta_P2_B2;

    // fprintf(stdout,"DP11: %lf DP12: %lf DP21: %lf DP22: %lf \n",delta_P1_B1,delta_P2_B1,delta_P1_B2,delta_P2_B2);
    // fprintf(stdout,"d1: %lf d2: %lf\n",d1_coeff,d2_coeff);

    for (int row = 0; row < 18; row++)
      for (int col = 0; col < 18; col++)
        fjac[row][col] = 0.0;

      
      
    // The Jacobian.
      
      double chi[6];
      
       // Here are the residuals for the characteristic matching for flow
      chi[0] = -2.0*theta*qA2_term*xb[ 2]/xb[11] + gamma*dFdQ(xb[11],VP);
      chi[2] =  2.0*theta*qA2_term*xb[ 5]/xb[14] + gamma*dFdQ(xb[14],VP);
      chi[4] =  2.0*theta*qA2_term*xb[ 8]/xb[17] + gamma*dFdQ(xb[17],VP);
      
      // Here are the residuals for the area characteristic matching
      chi[1] = theta*(qA2_term*sq(xb[2]/xb[11]) - PV->dBdAh(N,xb[11],WM)) +
                 gamma*(dFdA(xb[2],xb[11],VP) + PV->dB_RHSdAh(N,xb[11],WM));
      
      chi[3] = theta*( -qA2_term*sq(xb[5]/xb[14]) + B1->dBdAh(-1,xb[14],WM)) +
                gamma*(dFdA(xb[5],xb[14],VP) + B1->dB_RHSdAh(-1,xb[14],WM));
      
      chi[5] = theta*( -qA2_term*sq(xb[8]/xb[17]) + B2->dBdAh(-1,xb[17],WM)) +
                gamma*(dFdA(xb[8],xb[17],VP) + B2->dB_RHSdAh(-1,xb[17],WM));
      
      // Here are the pressure loss terms at n+1/2
      // Derivative with respect to flow
      double df14dx4 =  d1_coeff*(dp1/sq(xb[13]) + (dp2/sq(xb[13]))*2.0*xb[4]/ab(xb[4]) + dp3/(2.0*gamma*xb[13]));
      double df14dx10 = -PV->dPdA(N,xb[10],WM);
      double df14dx13 = B1->dPdA(0,xb[13],WM) - d1_coeff*(2.0*dp1*xb[4]/cu(xb[13])
      - 2.0*dp2*ab(xb[4])*xb[4]/cu(xb[13]) - (dp3/sq(xb[13])) * (xb[4]-(B1->Qold[0]))/(2.0*gamma));

      double df15dx7  =  d2_coeff*(dp1/sq(xb[16]) + (dp2/sq(xb[16]))*2.0*xb[7]/ab(xb[7]) + dp3/(2.0*gamma*xb[16]));
      double df15dx10 = -PV->dPdA(N,xb[10],WM);
      double df15dx16 = B2->dPdA(0,xb[16],WM) - d2_coeff*(2.0*dp1*xb[7]/cu(xb[16])
      - 2.0*dp2*ab(xb[7])*xb[7]/cu(xb[16]) - (dp3/sq(xb[16])) * (xb[7]-(B2->Qold[0]))/(2.0*gamma));

      double df16dx3  =  d1_coeff*(dp1/sq(xb[12]) + (dp2/sq(xb[12]))*2.0*xb[7]/ab(xb[3]) + dp3/(2.0*gamma*xb[12]));
      double df16dx9  = -PV->dPdA(N,xb[9],WM);
      double df16dx12 = B1->dPdA(0,xb[12],WM) - d1_coeff*(2.0*dp1*xb[3]/cu(xb[12])
      - 2.0*dp2*ab(xb[3])*xb[3]/cu(xb[12]) - (dp3/sq(xb[12])) * (xb[3]-(B1->Qold[0]))/(2.0*gamma));

      double df17dx6  =  d2_coeff*(dp1/sq(xb[15]) + (dp2/sq(xb[15]))*2.0*xb[6]/ab(xb[6]) + dp3/(2.0*gamma*xb[15]));
      double df17dx9 = -PV->dPdA(N,xb[9],WM);
      double df17dx15 = B2->dPdA(0,xb[15],WM) - d2_coeff*(2.0*dp1*xb[6]/cu(xb[15])
      - 2.0*dp2*ab(xb[6])*xb[6]/cu(xb[15]) - (dp3/sq(xb[15])) * (xb[6]-(B2->Qold[0]))/(2.0*gamma));


      
      //NEW JACOBIAN
            // Order is [row][column]
            fjac[ 0][ 0]  = -1.0;
            fjac[ 0][ 2] = chi[0];
            fjac[ 0][11] = chi[1];
            
            fjac[ 1][ 3] = -1.0;
            fjac[ 1][ 5] = chi[2];
            fjac[ 1][14] = chi[3];

            fjac[ 2][ 6] = -1.0;
            fjac[ 2][ 8] = chi[4];
            fjac[ 2][17] = chi[5];

            fjac[ 3][ 2] = -theta;
            fjac[ 3][ 9] = -1.0;
            
            fjac[ 4][ 5] = theta;
            fjac[ 4][12] = -1.0;
            
            fjac[ 5][ 8] = theta;
            fjac[ 5][15] = -1.0;
            
            fjac[ 6][ 1] = -1.0;
            fjac[ 6][ 2] = 0.5;
            
            fjac[ 7][ 4] = -1.0;
            fjac[ 7][ 5] = 0.5;
            
            fjac[ 8][ 7] = -1.0;
            fjac[ 8][ 8] = 0.5;
            
            fjac[ 9][ 10] = -1.0;
            fjac[ 9][ 11] = 0.5;
            
            fjac[10][13] = -1.0;
            fjac[10][14] = 0.5;
            
            fjac[11][16] = -1.0;
            fjac[11][17] = 0.5;
            
            fjac[12][ 1] = -1.0;
            fjac[12][ 4] = 1.0;
            fjac[12][ 7] = 1.0;
            
            fjac[13][ 0] = -1.0;
            fjac[13][ 3] = 1.0;
            fjac[13][ 6] = 1.0;
            
            fjac[14][ 4] = df14dx4;
            fjac[14][10] = df14dx10;
            fjac[14][13] = df14dx13;
                                    
                                    
            fjac[15][ 7] = df15dx7;
            fjac[15][10] = df15dx10;
            fjac[15][16] = df15dx16;
            
            fjac[16][ 3] = df16dx3;
            fjac[16][ 9] = df16dx9;
            fjac[16][12] = df16dx12;
            
            
            fjac[17][ 6] = df17dx6;
            fjac[17][ 9] = df17dx9;
            fjac[17][15] = df17dx15;
      
//       for (int kk=0; kk<18; kk++){
//       for (int jj=0; jj<18; jj++){
//         fprintf(stdout,"%lf  ",fjac[kk][jj]);
//       }
//       fprintf(stdout,"\n");
//     }

// fprintf(stdout,"\n\n\n");

    // Check whether solution is close enough. If not run the loop again.
    int ch = zero (xb, 18, 1.0e-8, 1.0e-8, fvec, fjac);
    if (ch == 1) ok = true;

    j = j+1;
  }
  
    // for (int kk=0; kk<18; kk++){
    //   for (int jj=0; jj<18; jj++){
    //     fprintf(stdout,"%lf  ",fjac[kk][jj]);
    //   }
    //   fprintf(stdout,"\n");
    // }

  // Solutions is applied, and right boundary is updated.
  PV->Anew[N] = xb[ 9];
  PV->Qnew[N] = xb[ 0];
  B1->Anew[0] = xb[12];
  B1->Qnew[0] = xb[ 3];
  B2->Anew[0] = xb[15];
  B2->Qnew[0] = xb[ 6];

  //  fprintf(stdout,"Daughter 1 Qold: %d Daughter 2: %d\n",*daughters,*(daughters+1);
  //  fprintf(stdout,"IN BIF: Parent Qold: %lf Daughter 1 Qold: %lf Daughter 2 Qold: %lf\n",PV->Qnew[N],B1->Qnew[0],B2->Qnew[0]);

    if (j >=ntrial) {error ("arteries.C","Root not found in the bifurcation");
        exit(1);}
}

void bound_trif (double theta, double gamma, Tube *Arteries[], int parent, int D1,int D2, int D3, int WM, int VP)
{
    // MJC: Try to define the daughters here
    Tube* PV  = Arteries[ parent];
    Tube* B1  = Arteries[ D1];
    Tube* B2  = Arteries[ D2];
    Tube* B3  = Arteries[ D3];
//    int ND = 3; // Number daughters: COME BACK AND MAKE XB AND FJAC FUNCTIONS OF THIS NUMBER
      int N = PV->N;
      double PN;
      int j = 1;
      int ok = false;
      const int ntrial = 100;
    
    double qA2_term = 0.0;
    if (VP==2) {
        qA2_term=((vel_power+2.0)/(vel_power+1.0)); // Power law
    } else if (VP<2) {

        qA2_term = 1.0; // Stokes boundary layer
    }
    else error ("arteries.cxx","Velocity profile doesn't exist.");
    
    double k1[4];
    double k2[4];
    double k3[4];
    double k4[4];
    
      // These are the flows at the half time step.
            k1[0] = PV->Qold[N]  + theta*(PV->R2h[N-1]) + gamma*(PV->S2h[N-1]);
          
            k1[1] = B1->Qold[0] - theta*(B1->R2h[0]) + gamma*(B1->S2h[0]);
          
            k1[2] = B2->Qold[0] - theta*(B2->R2h[0]) + gamma*(B2->S2h[0]);
          
            k1[3] = B3->Qold[0] - theta*(B3->R2h[0]) + gamma*(B3->S2h[0]);
          
            // These are areas at the half time step
            k2[0] = PV->Aold[N] + theta*(PV->R1h[N-1]);
            k2[1] = B1->Aold[0] - theta*(B1->R1h[0]);
            k2[2] = B2->Aold[0] - theta*(B2->R1h[0]);
            k2[3] = B3->Aold[0] - theta*(B3->R1h[0]);
          
            // These are flows at the half time stpe
            k3[0] = PV->Qh[N-1]/2.0;
            k3[1] = B1->Qh[0]/2.0;
            k3[2] = B2->Qh[0]/2.0;
            k3[3] = B3->Qh[0]/2.0;
          
            
            k4[0] = PV->Ah[N-1]/2.0;
            k4[1] = B1->Ah[0]/2.0;
            k4[2] = B2->Ah[0]/2.0;
            k4[3] = B3->Ah[0]/2.0;
          
            double xb[24];

            // The approximative initial guesses are applied.
              // Initial flows
          xb[ 0] =  PV->Qh[N-1];                        //Initial guess for Qp_xb n+1
          xb[ 1] = (PV->Qold[N-1] + PV->Qold[N])/2.0;   //Initial guess for Qp_xb^n+0.5
          xb[ 2] =  PV->Qold[N];                        //Initial guess for Qp_xb+0.5 n+0.5
          
          xb[ 3] =  B1->Qh[0];                          //Initial guess for Qd1_xb n+1
          xb[ 4] = (B1->Qold[0] + B1->Qold[1])/2.0;     //Initial guess for Qd1_xb n+0.5
          xb[ 5] =  B1->Qold[0];                        //Initial guess for Qd1_xb+0.5 n+0.5
          
          xb[ 6] =  B2->Qh[0];                          //Initial guess for Qd2_xb n+1
          xb[ 7] = (B2->Qold[0] + B2->Qold[1])/2.0;     //Initial guess for Qd2_xb n+0.5
          xb[ 8] =  B2->Qold[0];                        //Initial guess for Qd2_xb+0.5 n+0.5
          
          xb[ 9] =  B3->Qh[0];                          //Initial guess for Qd3_xb n+1
          xb[10] = (B3->Qold[0] + B3->Qold[1])/2.0;     //Initial guess for Qd3_xb n+0.5
          xb[11] =  B3->Qold[0];                        //Initial guess for Qd3_xb+0.5 n+0.5
          
          xb[12] =  PV->Ah[N-1];                        //Initial guess for Ap_xb n+1
          xb[13] = (PV->Aold[N-1] + PV->Aold[N])/2.0;   //Initial guess for Ap_xb^n+0.5
          xb[14] =  PV->Aold[N];                        //Initial guess for Ap_xb+0.5 n+0.5
          
          xb[15] =  B1->Ah[0];                          //Initial guess for Ad1_xb n+1
          xb[16] = (B1->Aold[0] + B1->Aold[1])/2.0;     //Initial guess for Ad1_xb n+0.5
          xb[17] =  B1->Aold[0];                        //Initial guess for Ad1_xb+0.5 n+0.5
          
          xb[18] =  B2->Ah[0];                          //Initial guess for Ad2_xb n+1
          xb[19] = (B2->Aold[0] + B2->Aold[1])/2.0;     //Initial guess for Ad2_xb n+0.5
          xb[20] =  B2->Aold[0];                        //Initial guess for Ad2_xb+0.5 n+0.5

          xb[21] =  B3->Ah[0];                          //Initial guess for Ad3_xb n+1
          xb[22] = (B3->Aold[0] + B3->Aold[1])/2.0;     //Initial guess for Ad3_xb n+0.5
          xb[23] =  B3->Aold[0];                        //Initial guess for Ad3_xb+0.5 n+0.5

          // This is where a Bernoulli term can be perscribed; set to zero otherwise
          double k7[4];
          double k7h[4];
          k7[ 0] = 0.0;
          k7[ 1] = 0.0;
          k7[ 2] = 0.0;
          k7h[0] = 0.0;
          k7h[1] = 0.0;
          k7h[2] = 0.0;

            // The residuals (fvec), and the Jacobian is determined, and if possible
            // the system of equations is solved.
            while (j <= ntrial && ok==false) // Find the zero
            {
              double fvec[24];
                
                
              // The residuals.
                
                // Characteristic Q residual at n+1
              fvec[0]  = k1[0]  - xb[0] -
                      theta*(qA2_term*sq(xb[2])/xb[14] + PV->Bh(N,xb[14],WM)) +
                      gamma*(F(xb[2],xb[14],VP)+PV->B_RHSh(N,xb[14],WM));

              fvec[1]  = k1[1]  - xb[3] +
                      theta*(qA2_term*sq(xb[5])/xb[17] + B1->Bh(-1,xb[17],WM)) +
                      gamma*(F(xb[5],xb[17],VP)  + B1->B_RHSh(-1,xb[17],WM));

              fvec[2]  = k1[2] - xb[6] +
                      theta*(qA2_term*sq(xb[8])/xb[20] + B2->Bh(-1,xb[20],WM)) +
                      gamma*(F(xb[8],xb[20],VP)  + B2->B_RHSh(-1,xb[20],WM));
                
              fvec[3]  = k1[3] - xb[9] +
                      theta*(qA2_term*sq(xb[11])/xb[23] + B3->Bh(-1,xb[23],WM)) +
                      gamma*(F(xb[11],xb[23],VP)  + B3->B_RHSh(-1,xb[23],WM));
                
                
                // Characteristic A residual at n+1
              fvec[4]  = - theta*xb[ 2] - xb[12]  + k2[0];
              fvec[5]  =   theta*xb[ 5] - xb[15]  + k2[1];
              fvec[6]  =   theta*xb[ 8] - xb[18]  + k2[2];
              fvec[7]  =   theta*xb[11] - xb[21]  + k2[3];
                
              // Flow residuals at n+1/2 (ghost points)
              fvec[ 8]  = - xb[ 1] + xb[ 2]/2.0 + k3[0];
              fvec[ 9]  = - xb[ 4] + xb[ 5]/2.0 + k3[1];
              fvec[10]  = - xb[ 7] + xb[ 8]/2.0 + k3[2];
              fvec[11]  = - xb[10] + xb[11]/2.0 + k3[3];
                
              // Area residuals at n+1/2 (ghost points)
              fvec[12] = - xb[13] + xb[14]/2.0 + k4[0];
              fvec[13] = - xb[16] + xb[17]/2.0 + k4[1];
              fvec[14] = - xb[19] + xb[20]/2.0 + k4[2];
              fvec[15] = - xb[22] + xb[23]/2.0 + k4[3];
                
                
              // Flow conservation residuals (n+1/2 and n+1)
              fvec[16] = - xb[ 1] + xb[ 4] + xb[7] + xb[10];
              fvec[17] = - xb[ 0] + xb[ 3] + xb[6] + xb[9];
                
                
              // Pressure continuity at the n+1/2 time step
              PN    = PV->P(N,xb[13],WM);
              double u_n_half = sq(xb[1]/xb[13]);

              // The if statements here only matter if a minor loss
              // is included.
                fvec[18] =  - PN + B1->P(0,xb[16],WM) + ab(k7h[0])*u_n_half;
                fvec[19] =  - PN + B2->P(0,xb[19],WM) + ab(k7h[1])*u_n_half;
                fvec[20] =  - PN + B3->P(0,xb[22],WM) + ab(k7h[2])*u_n_half;
              // Pressure continuity at the n+1 time step
              PN    = PV->P(N,xb[12],WM);
              double u_n_1 = sq(xb[0]/xb[12]);
                fvec[21] = - PN + B1->P(0,xb[15],WM) + ab(k7[0])*u_n_1;
                fvec[22] = - PN + B2->P(0,xb[18],WM) + ab(k7[1])*u_n_1;
                fvec[23] = - PN + B3->P(0,xb[21],WM) + ab(k7[2])*u_n_1;


              for (int row = 0; row < 4*6; row++)
                for (int col = 0; col < 4*6; col++)
                  fjac[row][col] = 0.0;

                
                
      //        // The Jacobian.
                double chi[16];
                
                // Here are the residuals for the characteristic matching for flow
                chi[0] = -2.0*theta*qA2_term*xb[ 2]/xb[14] + gamma*dFdQ(xb[14],VP);
                chi[2] =  2.0*theta*qA2_term*xb[ 5]/xb[17] + gamma*dFdQ(xb[17],VP);
                chi[4] =  2.0*theta*qA2_term*xb[ 8]/xb[20] + gamma*dFdQ(xb[20],VP);
                chi[6] =  2.0*theta*qA2_term*xb[11]/xb[23] + gamma*dFdQ(xb[23],VP);
                
                
                // Here are the residuals for the area characteristic matching
                chi[1] = theta*(qA2_term*sq(xb[2]/xb[14]) - PV->dBdAh(N,xb[14],WM)) +
                           gamma*(dFdA(xb[2],xb[14],VP) + PV->dB_RHSdAh(N,xb[14],WM));
                
                chi[3] = theta*( -qA2_term*sq(xb[5]/xb[17]) + B1->dBdAh(-1,xb[17],WM)) +
                          gamma*(dFdA(xb[5],xb[17],VP) + B1->dB_RHSdAh(-1,xb[17],WM));
                
                chi[5] = theta*( -qA2_term*sq(xb[8]/xb[20]) + B2->dBdAh(-1,xb[20],WM)) +
                          gamma*(dFdA(xb[8],xb[20],VP) + B2->dB_RHSdAh(-1,xb[20],WM));
                
                chi[7] = theta*( -qA2_term*sq(xb[11]/xb[23]) + B3->dBdAh(-1,xb[23],WM)) +
                          gamma*(dFdA(xb[11],xb[23],VP) + B3->dB_RHSdAh(-1,xb[23],WM));
                
                // Here is pressure conservation (n+1/2)
                chi[8]  = -PV->dPdA(N,xb[13],WM) + sq(xb[1])/cu(xb[13])*(-2.0*ab(k7h[0])); //Loss term
                chi[9]  = B1->dPdA(0,xb[16],WM);
                chi[10] = B2->dPdA(0,xb[19],WM);
                chi[11] = B3->dPdA(0,xb[22],WM);
                
                // Here is pressure conservation (n+1)
                chi[12] = -PV->dPdA(N,xb[12],WM) + sq(xb[1])/cu(xb[12])*(-2.0*ab(k7h[0])); //Loss term
                chi[13] = B1->dPdA(0,xb[15],WM);
                chi[14] = B2->dPdA(0,xb[18],WM);
                chi[15] = B3->dPdA(0,xb[21],WM);
                
                // The jacobian
                // Order is [row][column]
                
                fjac[0][ 0] = -1.0;
                fjac[0][ 2] = chi[0];
                fjac[0][14] = chi[1];
                
                fjac[1][ 3] = -1.0;
                fjac[1][ 5] = chi[2];
                fjac[1][17] = chi[3];
                
                fjac[2][ 6] = -1.0;
                fjac[2][ 8] = chi[4];
                fjac[2][20] = chi[5];
                
                fjac[3][ 9] = -1.0;
                fjac[3][11] = chi[6];
                fjac[3][23] = chi[7];
                
                fjac[4][ 2] = -theta;
                fjac[4][12] = -1.0;
                
                fjac[5][ 5] = theta;
                fjac[5][15] = -1.0;
                
                fjac[6][ 8] = theta;
                fjac[6][18] = -1.0;
                
                fjac[7][11] = theta;
                fjac[7][21] = -1.0;
                 
                fjac[8][ 1] = -1.0;
                fjac[8][ 2] = 0.5;
                
                fjac[9][ 4] = -1.0;
                fjac[9][ 5] = 0.5;
                
                fjac[10][ 7] = -1.0;
                fjac[10][ 8] = 0.5;
                
                fjac[11][10] = -1.0;
                fjac[11][11] = 0.5;
                
                fjac[12][13] = -1.0;
                fjac[12][14] = 0.5;
                
                fjac[13][16] = -1.0;
                fjac[13][17] = 0.5;
                
                fjac[14][19] = -1.0;
                fjac[14][20] = 0.5;
                
                fjac[15][22] = -1.0;
                fjac[15][23] = 0.5;
                
                fjac[16][ 1] = -1.0;
                fjac[16][ 4] = 1.0;
                fjac[16][ 7] = 1.0;
                fjac[16][10] = 1.0;
                
                fjac[17][ 0] = -1.0;
                fjac[17][ 3] = 1.0;
                fjac[17][ 6] = 1.0;
                fjac[17][ 9] = 1.0;
                
                fjac[18][13] = chi[8];
                fjac[18][16] = chi[9];
                
                fjac[19][13] = chi[8];
                fjac[19][19] = chi[10];
                
                fjac[20][13] = chi[8];
                fjac[20][22] = chi[11];
                
                fjac[21][12] = chi[12];
                fjac[21][15] = chi[13];
                
                fjac[22][12] = chi[12];
                fjac[22][18] = chi[14];
                
                fjac[23][12] = chi[12];
                fjac[23][21] = chi[15];
          

        
        // Check whether solution is close enough. If not run the loop again.
        int ch = zero (xb, 24, 1.0e-12, 1.0e-12, fvec, fjac);
        if (ch == 1) ok = true;

        j = j+1;
      }

      // Solutions is applied, and right boundary is updated.
      PV->Anew[N] = xb[ 12];
      PV->Qnew[N] = xb[ 0];
      B1->Anew[0] = xb[15];
      B1->Qnew[0] = xb[ 3];
      B2->Anew[0] = xb[18];
      B2->Qnew[0] = xb[ 6];
      B3->Anew[0] = xb[21];
      B3->Qnew[0] = xb[ 9];
//
        if (j >=ntrial) {error ("arteries.C","Root not found in the bifurcation");
            exit(1);}
    }






