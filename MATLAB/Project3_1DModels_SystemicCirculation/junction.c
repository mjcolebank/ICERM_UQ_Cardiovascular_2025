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


void junction (double theta, double gamma, Tube *Arteries[], int parent, int conn_id, int WM)
{
    int D1,D2,D3;//,D4;
    
        D1 = Arteries[parent]->daughters[max_D*conn_id + 1]; //change indexing of second parent argument
        D2 = Arteries[parent]->daughters[max_D*conn_id + 2];
        D3 = Arteries[parent]->daughters[max_D*conn_id + 3];

//    fprintf(stdout,"p d1 d2: %d %d %d\n",parent,D1,D2);

    if (D2==0) {
        bound_monf(theta, gamma, Arteries, parent,D1,WM); // Only call mono if wanting a stenosis
    }
    else if (D3==0 && D2>0){
        bound_bif(theta, gamma, Arteries, parent,D1,D2,WM);
    }
    else if (D1<0 && D2<0){
        bound_bif_converge(theta, gamma, Arteries, parent,-D1,-D2,WM);
    }
    else
    {
        bound_trif(theta, gamma, Arteries, parent,D1,D2,D3,WM);
    }
    // Can add quadfurcation here if necessary.
}


//begin diverging monofurcation
void bound_monf (double theta, double gamma, Tube *Arteries[], int parent, int D1, int WM)
{
    // MJC: Try to define the daughters here
    Tube* PV = Arteries[ parent]; // Parent vessel
    Tube* B1 = Arteries[ D1];      // Branch 1
    int N = PV->N;
    int j = 1;
    int ok = false;
    const int ntrial = 100;
    
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

    
    
    double k7nh  = 0.0;//B1->K_loss/2.0; //32*mu/(2*B1->rtop*rho*q); Bernoulli loss
    double k7n   = 0.0;//B1->K_loss/2.0; 

    
    // The residuals (fvec), and the Jacobian is determined, and if possible
    // the system of equations is solved.
    while (j <= ntrial && ok==false) // Find the zero
    {
        double fvec[12];
        
        // The residuals, 
        
        // Characteristic Q residual at n+1
        fvec[0]  = k1[0]  - xb[0] -
        theta*(sq(xb[2])/xb[8] + PV->Bh(N,xb[8],WM)) +
        gamma*(F(xb[2],xb[8])  + PV->dBdx1h(N,xb[8],WM));
        
        fvec[1]  = k1[1]  - xb[3] +
        theta*(sq(xb[5])/xb[11] + B1->Bh(-1,xb[11],WM)) +
        gamma*(F(xb[5],xb[11])  + B1->dBdx1h(-1,xb[11],WM));
        
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
        chi[0] = -2.0*theta*xb[ 2]/xb[8] + gamma*dFdQ(xb[8]);
        chi[2] =  2.0*theta*xb[ 5]/xb[11] + gamma*dFdQ(xb[11]);
        
        // Here are the residuals for the area characteristic matching
        chi[1] = theta*(sq(xb[2]/xb[8]) - PV->dBdAh(N,xb[8],WM)) +
                   gamma*(dFdA(xb[2],xb[8]) + PV->d2BdAdxh(N,xb[8],WM));
        
        chi[3] = theta*( -sq(xb[5]/xb[11]) + B1->dBdAh(-1,xb[11],WM)) +
                  gamma*(dFdA(xb[5],xb[11]) + B1->d2BdAdxh(-1,xb[11],WM));
        
        // Here is pressure conservation (n+1/2)
        chi[4]  = -PV->dPdA(N,xb[7],WM) + sq(xb[1])/cu(xb[7])*(-2.0*ab(k7n)); //Loss term
        chi[5]  =  B1->dPdA(0,xb[10],WM);
        
        // Here is pressure conservation (n+1)
        chi[6] = -PV->dPdA(N,xb[6],WM) + sq(xb[0])/cu(xb[6])*(-2.0*ab(k7nh)); //Loss term
        chi[7] =  B1->dPdA(0,xb[9],WM);
                                              
        
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
    
//    fprintf(stdout,"DISCREP:%lf\n",(xb[6]-xb[9])/xb[6]);

      if (j >=ntrial) error ("arteries.C","Root not found in the bifurcation");
}
//end diverging monofurcation

//begin diverging bifurcation
void bound_bif (double theta, double gamma, Tube *Arteries[], int parent, int D1, int D2, int WM)
{
    // MJC: Try to define the daughters here
    Tube* PV  = Arteries[ parent];
    Tube* B1  = Arteries[ D1];
    Tube* B2  = Arteries[ D2];
    
    
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

  double xb[6*3];

  // The approximative initial guesses are applied.
  //anywhere F need gravity term, and dFdA
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

  double k7nh  = 0;
  double k7n   = 0;
  double k7anh = 0;
  double k7an  = 0;

  // The residuals (fvec), and the Jacobian is determined, and if possible
  // the system of equations is solved.
  while (j <= ntrial && ok==false) // Find the zero
  {
    double fvec[18];
    // The residuals.
      
      // Characteristic Q residual at n+1, gravity term here
    fvec[0]  = k1[0]  - xb[0] -
    theta*(sq(xb[2])/xb[11] + PV->Bh(N,xb[11],WM)) +
    gamma*(F(xb[2],xb[11])+PV->dBdx1h(N,xb[11],WM));

    fvec[1]  = k1[1]  - xb[3] +
                 theta*(sq(xb[5])/xb[14] + B1->Bh(-1,xb[14],WM)) +
                 gamma*(F(xb[5],xb[14])  + B1->dBdx1h(-1,xb[14],WM));

    fvec[2]  = k1[2] - xb[6] +
    theta*(sq(xb[8])/xb[17] + B2->Bh(-1,xb[17],WM)) +
    gamma*(F(xb[8],xb[17])  + B2->dBdx1h(-1,xb[17],WM));

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
      chi[0] = -2.0*theta*xb[ 2]/xb[11] + gamma*dFdQ(xb[11]);
      chi[2] =  2.0*theta*xb[ 5]/xb[14] + gamma*dFdQ(xb[14]);
      chi[4] =  2.0*theta*xb[ 8]/xb[17] + gamma*dFdQ(xb[17]);
      
      // Here are the residuals for the area characteristic matching
      chi[1] = theta*(sq(xb[2]/xb[11]) - PV->dBdAh(N,xb[11],WM)) +
                 gamma*(dFdA(xb[2],xb[11]) + PV->d2BdAdxh(N,xb[11],WM) + dGdA(PV->ves_angle));
      
      chi[3] = theta*( -sq(xb[5]/xb[14]) + B1->dBdAh(-1,xb[14],WM)) +
                gamma*(dFdA(xb[5],xb[14]) + B1->d2BdAdxh(-1,xb[14],WM) + dGdA(B1->ves_angle));
      
      chi[5] = theta*( -sq(xb[8]/xb[17]) + B2->dBdAh(-1,xb[17],WM)) +
                gamma*(dFdA(xb[8],xb[17]) + B2->d2BdAdxh(-1,xb[17],WM) + dGdA(B2->ves_angle));
      
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
            
      //      fjac[14][ 1] = J_lossP1;
            fjac[14][10] = chi[6];
            fjac[14][13] = chi[7];
                                    
                                    
      //      fjac[15][ 1] = J_lossP2;
            fjac[15][10] = chi[6];
            fjac[15][16] = chi[8];
            
      //      fjac[16][ 0] = J_lossP3;
            fjac[16][ 9] = chi[9];
            fjac[16][12] = chi[10];
            
            
      //      fjac[17][ 0] = J_lossP4;
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
//end diverging bifurcation

// begin diverging trifurcation
void bound_trif (double theta, double gamma, Tube *Arteries[], int parent, int D1,int D2, int D3, int WM)
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
          
            double xb[6*4];

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
                      theta*(sq(xb[2])/xb[14] + PV->Bh(N,xb[14],WM)) +
                      gamma*(F(xb[2],xb[14])+PV->dBdx1h(N,xb[14],WM));

              fvec[1]  = k1[1]  - xb[3] +
                      theta*(sq(xb[5])/xb[17] + B1->Bh(-1,xb[17],WM)) +
                      gamma*(F(xb[5],xb[17])  + B1->dBdx1h(-1,xb[17],WM));

              fvec[2]  = k1[2] - xb[6] +
                      theta*(sq(xb[8])/xb[20] + B2->Bh(-1,xb[20],WM)) +
                      gamma*(F(xb[8],xb[20])  + B2->dBdx1h(-1,xb[20],WM));
                
              fvec[3]  = k1[3] - xb[9] +
                      theta*(sq(xb[11])/xb[23] + B3->Bh(-1,xb[23],WM)) +
                      gamma*(F(xb[11],xb[23])  + B3->dBdx1h(-1,xb[23],WM));
                
                
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
                chi[0] = -2.0*theta*xb[ 2]/xb[14] + gamma*dFdQ(xb[14]);
                chi[2] =  2.0*theta*xb[ 5]/xb[17] + gamma*dFdQ(xb[17]);
                chi[4] =  2.0*theta*xb[ 8]/xb[20] + gamma*dFdQ(xb[20]);
                chi[6] =  2.0*theta*xb[11]/xb[23] + gamma*dFdQ(xb[23]);
                
                
                // Here are the residuals for the area characteristic matching
                chi[1] = theta*(sq(xb[2]/xb[14]) - PV->dBdAh(N,xb[14],WM)) +
                           gamma*(dFdA(xb[2],xb[14]) + PV->d2BdAdxh(N,xb[14],WM));
                
                chi[3] = theta*( -sq(xb[5]/xb[17]) + B1->dBdAh(-1,xb[17],WM)) +
                          gamma*(dFdA(xb[5],xb[17]) + B1->d2BdAdxh(-1,xb[17],WM));
                
                chi[5] = theta*( -sq(xb[8]/xb[20]) + B2->dBdAh(-1,xb[20],WM)) +
                          gamma*(dFdA(xb[8],xb[20]) + B2->d2BdAdxh(-1,xb[20],WM));
                
                chi[7] = theta*( -sq(xb[11]/xb[23]) + B3->dBdAh(-1,xb[23],WM)) +
                          gamma*(dFdA(xb[11],xb[23]) + B3->d2BdAdxh(-1,xb[23],WM));
                
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
          
          // For debugging, you can print the jacobian
//                  if(j==1){
//                      printf("X\n");
//                      for (int i=0; i<24; i++)
//                      {
//                       printf("%lf ",xb[i]);
//                      }
//                      printf("\n\n");
//                      printf("JACOBIAN\n");
//                      for (int i=0; i<24; i++)
//                      {
//                          for (int j=0; j<24; j++) {
//                              printf("%lf ",fjac[i][j]);
//                          }
//                          printf("\n");
//                      }
//                      printf("\n\n");
//                  }
        
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

// end diverging trifurcation

// begin converging flow conditions
void bound_bif_converge (double theta, double gamma, Tube *Arteries[], int daughter, int P1, int P2, int WM)
{
//       fprintf(stdout,"d p1 p2: %d %d %d\n",daughter,P1,P2);
    Tube* PV  = Arteries[ daughter];
    Tube* B1  = Arteries[ P1];
    Tube* B2  = Arteries[ P2];
    

    
  int B1N = B1->N;
  int B2N = B2->N;
  double PN;
  int j = 1;
  int ok = false;
  const int ntrial = 500;
    

    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];
    
  k1[0] = PV->Qold[0]  - theta*(PV->R2h[0]) + gamma*(PV->S2h[0]);
  k1[1] = B1->Qold[B1N] + theta*(B1->R2h[B1N-1]) + gamma*(B1->S2h[B1N-1]);
  k1[2]= B2->Qold[B2N] + theta*(B2->R2h[B2N-1]) + gamma*(B2->S2h[B2N-1]);
    

  k2[0] = PV->Aold[0] - theta*(PV->R1h[0]);
  k2[1] = B1->Aold[B1N] + theta*(B1->R1h[B1N-1]);
  k2[2] = B2->Aold[B2N] + theta*(B2->R1h[B2N-1]);

  k3[0] = PV->Qh[0]/2.0;
  k3[1] = B1->Qh[B1N-1]/2.0;
  k3[2] = B2->Qh[B2N-1]/2.0;

  k4[0]  = PV->Ah[0]/2.0;
  k4[1]  = B1->Ah[B1N-1]/2.0;
  k4[2]  = B2->Ah[B2N-1]/2.0;

  double xb[18];

  // The approximative initial guesses are applied.
  xb[ 0] =  PV->Qh[0];                        //Initial guess for Q1_xb n+1
  xb[ 1] = (PV->Qold[0] + PV->Qold[1])/2.0;       //Initial guess for Q1_xb^n+0.5
  xb[ 2] =  PV->Qold[0];                        //Initial guess for Q1_xb+0.5 n+0.5
  xb[ 3] =  B1->Qh[B1N-1];                      //Initial guess for Q2_xb n+1
  xb[ 4] = (B1->Qold[B1N-1] + B1->Qold[B1N])/2.0; //Initial guess for Q2_xb n+0.5
  xb[ 5] =  B1->Qold[B1N];                    //Initial guess for Q2_xb+0.5 n+0.5
  xb[ 6] =  B2->Qh[B2N-1];                      //Initial guess for Q3_xb n+1
  xb[ 7] = (B2->Qold[B2N-1] + B2->Qold[B2N])/2.0; //Initial guess for Q3_xb n+0.5
  xb[ 8] =  B2->Qold[B2N];                    //Initial guess for Q3_xb+0.5 n+0.5
  xb[ 9] =  PV->Ah[0];                        //Initial guess for A1_xb n+1
  xb[10] = (PV->Aold[0] + PV->Aold[1])/2.0;       //Initial guess for A1_xb^n+0.5
  xb[11] =  PV->Aold[0];                        //Initial guess for A1_xb+0.5 n+0.5
  xb[12] =  B1->Ah[B1N-1];                      //Initial guess for A2_xb n+1
  xb[13] = (B1->Aold[B1N-1] + B1->Aold[B1N])/2.0; //Initial guess for A2_xb n+0.5
  xb[14] =  B1->Aold[B1N];                    //Initial guess for A2_xb+0.5 n+0.5
  xb[15] =  B2->Ah[B2N-1];                      //Initial guess for A3_xb n+1
  xb[16] = (B2->Aold[B2N-1] + B2->Aold[B2N])/2.0; //Initial guess for A3_xb n+0.5
  xb[17] =  B2->Aold[B2N];                    //Initial guess for A3_xb+0.5 n+0.5

  double k7nh  = 0.0; //Bernoulli loss term
  double k7n   = 0.0;
  double k7anh = 0.0;
  double k7an  = 0.0;

  // The residuals (fvec), and the Jacobian is determined, and if possible
  // the system of equations is solved.
  while (j <= ntrial && ok==false) // Find the zero
  {
    double fvec[18];
    // The residuals.
      
      // Characteristic Q residual at n+1
//    fvec[0]  = k1[0]  - xb[0] -
//    theta*(sq(xb[2])/xb[11] + PV->Bh(-1,xb[11],WM)) +
//    gamma*(F(xb[2],xb[11])+PV->dBdx1h(-1,xb[11],WM));
      
      // MJC
      fvec[0]  = k1[0]  - xb[0] +
         theta*(sq(xb[2])/xb[11] + PV->Bh(-1,xb[11],WM)) +
         gamma*(F(xb[2],xb[11])+PV->dBdx1h(-1,xb[11],WM));

    fvec[1]  = k1[1]  - xb[3] -
                 theta*(sq(xb[5])/xb[14] + B1->Bh(B1N,xb[14],WM)) +
                 gamma*(F(xb[5],xb[14])  + B1->dBdx1h(B1N,xb[14],WM));

    fvec[2]  = k1[2] - xb[6] -
    theta*(sq(xb[8])/xb[17] + B2->Bh(B2N,xb[17],WM)) +
    gamma*(F(xb[8],xb[17])  + B2->dBdx1h(B2N,xb[17],WM));

    // Characteristic A residual at n+1
    fvec[3]  =   theta*xb[2] - xb[ 9] + k2[0];
    fvec[4]  =  -theta*xb[5] - xb[12] + k2[1];
    fvec[5]  =  -theta*xb[8] - xb[15] + k2[2];
      
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
    PN    = PV->P(0,xb[10],WM);
    double u_n_half = sq(xb[1]/xb[10]);

    // The if statements here only matter if a minor loss
    // is included
      fvec[14] =  - PN + B1->P(B1N,xb[13],WM) + ab(k7nh)*u_n_half;
      fvec[15] =  - PN + B2->P(B2N,xb[16],WM) + ab(k7anh)*u_n_half;

    // Pressure continuity at the n+1 time step
    PN    = PV->P(0,xb[9],WM);
    double u_n_1 = sq(xb[0]/xb[9]);
      fvec[16] = - PN + B1->P(B1N,xb[12],WM) + ab(k7n)*u_n_1;
      fvec[17] = - PN + B2->P(B2N,xb[15],WM) + ab(k7an)*u_n_1;

    for (int row = 0; row < 18; row++)
      for (int col = 0; col < 18; col++)
        fjac[row][col] = 0.0;

      
      
    // The Jacobian.
      
      double chi[12];
      
       // Here are the residuals for the characteristic matching for flow
      chi[0] =  2.0*theta*xb[ 2]/xb[11] + gamma*dFdQ(xb[11]);
      chi[2] = -2.0*theta*xb[ 5]/xb[14] + gamma*dFdQ(xb[14]);
      chi[4] = -2.0*theta*xb[ 8]/xb[17] + gamma*dFdQ(xb[17]);
      
      // Here are the residuals for the area characteristic matching
      chi[1] = theta*(-sq(xb[2]/xb[11]) + PV->dBdAh(-1,xb[11],WM)) +
                 gamma*(dFdA(xb[2],xb[11]) + PV->d2BdAdxh(-1,xb[11],WM) + dGdA(PV->ves_angle));
      
      chi[3] = theta*( sq(xb[5]/xb[14]) - B1->dBdAh(B1N,xb[14],WM)) +
                gamma*(dFdA(xb[5],xb[14]) + B1->d2BdAdxh(B1N,xb[14],WM) + dGdA(B1->ves_angle));
      
      chi[5] = theta*( sq(xb[8]/xb[17]) - B2->dBdAh(B2N,xb[17],WM)) +
                gamma*(dFdA(xb[8],xb[17]) + B2->d2BdAdxh(B2N,xb[17],WM) + dGdA(B2->ves_angle));
      
      // Here is pressure conservation
      chi[6]  = -PV->dPdA(0,xb[10],WM) + sq(xb[1])/cu(xb[10])*(-2.0*ab(k7nh)); //Loss term
      chi[7]  = B1->dPdA(B1N,xb[13],WM);
      chi[8]  = B2->dPdA(B2N,xb[16],WM);
      
      chi[9]  = -PV->dPdA(0,xb[9],WM) + sq(xb[0])/cu(xb[9])*(-2.0*ab(k7nh)); //Loss term
      chi[10] = B1->dPdA(B1N,xb[12],WM);
      chi[11] = B2->dPdA(B2N,xb[15],WM);

      
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

//            fjac[ 3][ 2] = -theta;
            // MJC
            fjac[ 3][ 2] = theta;
            fjac[ 3][ 9] = -1.0;
            
//            fjac[ 4][ 5] = theta;
            // MJC
            fjac[ 4][ 5] = -theta;
            fjac[ 4][12] = -1.0;
            
//            fjac[ 5][ 8] = theta;
            fjac[ 5][ 8] = -theta;
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

      // MJC print statements
//       fprintf(stdout,"IN BIF CONV: Daughter Q: %lf Parent 1 Q: %lf Parent 2 Q: %lf\n",xb[0],xb[3],xb[6]);
//       fprintf(stdout,"IN BIF CONV: Daughter A: %lf Parent 1 A: %lf Parent 2 A: %lf\n",xb[9],xb[12],xb[15]);
    j = j+1;
  }

  // Solutions is applied, and right boundary is updated.
  PV->Anew[0]   = xb[ 9];
  PV->Qnew[0]   = xb[ 0];
  B1->Anew[B1N] = xb[12];
  B1->Qnew[B1N] = xb[ 3];
  B2->Anew[B2N] = xb[15];
  B2->Qnew[B2N] = xb[ 6];
    
//   fprintf(stdout,"IN BIF CONV: Daughter Qold: %lf Parent 1 Qold: %lf Parent 2 Qold: %lf\n",PV->Qold[0],B1->Qold[B1N],B2->Qold[B2N]);
//     fprintf(stdout,"IN BIF CONV: Daughter A: %lf Parent 1 A: %lf Parent 2 A: %lf\n",xb[9],xb[12],xb[15]);


    if (j >=ntrial) {error ("arteries.C","Root not found in the bifurcation");
        exit(1);}
}

//end converging flow conditions



