/***************************************************************************/
/*                                                                         */
/* The sor06.C main program                                                */
/*  Version: 1.0                                                           */
/*  Date: 30 Dec. 2019                                                     */
/*                                                                         */
/*  Primary Authors: M.S. Olufsen                                          */
/*  Key Contributers: M.U. Qureshi & M.J. Colebank                         */
/*  Institute: North Carolina State University, Raleigh, NC, USA           */
/*  Funded by: NSF-DMS # 1615820 & AHA #19PRE34380459                      */
/*                                                                         */
/***************************************************************************/

#include "sor06.h"
#include "tools.h"
#include "arteries.h"
#include "junction.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>

using namespace std;

// The vessel-network is specified and initialized, and the flow and
// pressures are to be determed. All global constants must be defined
// in the header file. That is the number of dots per cm, and the number
// of vessels.
int main(int argc, char *argv[])
{
    double tstart, tend, finaltime;
    // Declare stiffness parameters
    double f1, f2, f3;
    double Rp1,Rp2,Rd1,Rd2,CT1,CT2;
    // f1 = 0.0;
    // f2 = 0.0;
    // Declare the parameters describing the network size
    int total_vessels, total_terminal, total_conn,number_of_points;
    int total_sten, fileID; // Number of stenoses
    fflush(stdout);
    //==========adjustment in resistances (WK parameters) for the control========
    
    if (argc != 11) //argv[0] is the name of the program, here sor06
    {
        fprintf(stdout,"Not enough input arguments: only %d. Exiting. \n", argc);
        fflush(stdout);
        return 1;
    }
    
    f1     = atof(argv[1]);
    f2     = atof(argv[2]);
    f3     = atof(argv[3]);
    Rp1    = atof(argv[4]);
    Rp2    = atof(argv[5]);
    Rd1    = atof(argv[6]);
    Rd2    = atof(argv[7]);
    CT1    = atof(argv[8]);
    CT2    = atof(argv[9]);
    fileID = atoi(argv[10]);


    total_vessels    = 3;
    total_terminal   = 2;
    number_of_points = 10;

    total_sten = 0;
    
    total_conn         = total_vessels-total_terminal;
    nbrves             = total_vessels;
            
    /* Declare string vectors and files to hold 
     * the output from the model. This can 
     * be expanded to have more vessels
     */
    
    char name_all [20];
    sprintf(name_all, "output_%d.2d", fileID);
    FILE *art_ALL = fopen (name_all, "w");
    
    
    // Workspace used by bound_bif
    for(int i=0; i<24; i++) fjac[i] = new double[24];
    tstart    = 0.0;            // Starting time.
    
    // The number of vessels in the network is given when the governing array of
    // vessels is declared.
    
// =========================NETWORK =================================================================================
    
    Tube   *Arteries[nbrves];                     // Array of blood vessels.
    int connectivity_matrix[total_vessels][max_D];// Matrix with vessels
    int terminal_vessels[total_terminal];         // ID's for term. ves.
    double dimensions_matrix[total_vessels][3];   // Length and radius
    FILE *conn;
    conn = fopen("connectivity.txt","rt");
    int parent, daughter1, daughter2, daughter3, r_in;
    int conn_id = 0;
    
    
    // Check to see if we have the connectivity file
   if (conn == NULL)
   {
       fprintf(stdout,"Error: Connectivity File Does Not Exist \n");
       return 1;
   }
    while ((r_in = fscanf(conn, "%d %d %d %d", &parent, &daughter1, &daughter2, &daughter3)) != EOF)
    {
        connectivity_matrix[conn_id][0] = parent;
        connectivity_matrix[conn_id][1] = daughter1;
        connectivity_matrix[conn_id][2] = daughter2;
        connectivity_matrix[conn_id][3] = daughter3;
        conn_id++;
    }
    fclose(conn);
    

    dimensions_matrix[0][0] = 4.30;
    dimensions_matrix[1][0] = 2.50;
    dimensions_matrix[2][0] = 5.75;

    dimensions_matrix[0][1] = 1.35;
    dimensions_matrix[1][1] = 0.90;
    dimensions_matrix[2][1] = 1.10;

    dimensions_matrix[0][2] = 1.35;
    dimensions_matrix[1][2] = 0.90;
    dimensions_matrix[2][2] = 1.10;

    terminal_vessels[0] = 1;
    terminal_vessels[1] = 2;
    
    
   
   
    /* Initialization of the Arteries.
    // The tube class takes in the following values (in this order)
         Length: the length of the artery
     
         Top radius: the radius at the entry of the vessel
     
         Bottom radius: the radius at the end of the vessel
     
         Daughter Pointer: the pointer to the connectivity of all the blood
                           vessels, passed to junction.c
          
         Num Points: the number of pts (per non-dimensional length) to use in numerical solution
     
         Init: Set to 1 to make this vessel the inlet vessel (i.e., it gets a flow prescribed); set 0 else
     
         f1,f2,f3: Stiffness in the large blood vessels (Eh/r0 = f1*exp(f2*r0)+f3)
     
         R1,R2,CT: the proximal and distal resistance parameters and the peripheral compliance, used in the windkessel boundary condition
     
        */
   
    // MJC: For trying to pass junction condition around
    int *daughter_ptr = &connectivity_matrix[0][0];

    Arteries[2] = new Tube( dimensions_matrix[2][0], dimensions_matrix[2][1], dimensions_matrix[2][2],
                                    daughter_ptr, number_of_points, 0,f1,f2,f3, Rp2, Rd2, CT2, 0,0);
    Arteries[1] = new Tube( dimensions_matrix[1][0], dimensions_matrix[1][1], dimensions_matrix[1][2],
                                    daughter_ptr, number_of_points, 0,f1,f2,f3, Rp1, Rd1, CT1, 0,0);
    Arteries[0] = new Tube( dimensions_matrix[0][0], dimensions_matrix[0][1], dimensions_matrix[0][2],
                                    daughter_ptr, number_of_points, 1,f1,f2,f3, 0,0,0,0,0);

   
    
    // Solves the equations until time equals tend.///////////////////////////
    /* ADDED BY M. J. Colebank
     * Rather than specifying the number of cycles as an input to the function,
     * we want to test to see if the solution has converged. If so, we should exit.*/
    
    int period_counter = 1; // Count the number of periods you have solved for
    double norm_sol = 1e+6;
    double sol_tol  = 1e-5;
    double sol_p1[tmstps],sol_p2[tmstps];
    tend      = Deltat;
    
    
    // Solve the model once for initial convergence diagnostics
    int sol_ID = 0;
    while (tend<=period_counter*Period)
    {
    solver (Arteries, tstart, tend, k, Period, WALL_MODEL, VEL_MODEL);
    sol_p1[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0], WALL_MODEL);
    sol_p1[sol_ID] *= rho*g*Lr/conv;
    tstart = tend;
    tend   = tend + Deltat; // The current ending time is increased by Deltat.
    sol_ID++;
    }
    
    
    // Loop for convergence
    double sse;
    while (norm_sol>=sol_tol)
    {
        sol_ID = 0;
        sse    = 0;
        period_counter++;
        if (period_counter>max_cycles)
        {
            printf("ERROR: TOO MANY CYCLES (BAD PARAMETERS). EXITING. \n");
            fflush(stdout);
            return(1);
        }
        while (tend<=period_counter*Period)
        {
            solver (Arteries, tstart, tend, k, Period, WALL_MODEL, VEL_MODEL);
            sol_p2[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0],WALL_MODEL);
            sol_p2[sol_ID] *= rho*g*Lr/conv;
            sse = sse+ sq(sol_p1[sol_ID]-sol_p2[sol_ID]);
            tstart = tend;
            tend   = tend + Deltat; // The current ending time is increased by Deltat.
            sol_ID++;
        }
        norm_sol = sse;
        memcpy (sol_p1, sol_p2, sizeof(sol_p2));
    }

    //////////////////////////////////////////////
  // The loop is continued until the final time
  // is reached. If one wants to make a plot of
  // the solution versus x, tend is set to final-
  // time in the above declaration.

  period_counter++;
  finaltime = (period_counter+(cycles-1))*Period;
    while (tend <= finaltime)
    {
        for (int j=0; j<nbrves; j++)
        {
            int ArtjN = Arteries[j]->N;
            for (int i=0; i<ArtjN; i++)
            {
                Arteries[j]->Qprv[i+1] = Arteries[j]->Qnew[i+1];
                Arteries[j]->Aprv[i+1] = Arteries[j]->Anew[i+1];
            }
        }
        
        // Solves the equations until time equals tend.
        solver (Arteries, tstart, tend, k, Period, WALL_MODEL, VEL_MODEL);
        
        // A 2D plot of P(x_fixed,t) is made. The resulting 2D graph is due to
        // the while loop, since for each call of printPt only one point is set.
        // To print more vessels, use:
        //
        // Arteries[ connectivity_matrix[ROW][COL]] -> printPxt (artI, tend, 0, WALL_MODEL);
        
        
        // Prints the proximal, midpoint, and distal waveforms for
        // every vessel and combines it all to one file.
        
        for (int save_id=0; save_id<nbrves; save_id++)
         {

             Arteries[ save_id] -> printAllt (art_ALL, tend, 0, WALL_MODEL);
            
         }
        for (int save_id=0; save_id<nbrves; save_id++)
        {
            Arteries[ save_id] -> printAllt (art_ALL, tend, Arteries[save_id]->N/2, WALL_MODEL);
        }
        for (int save_id=0; save_id<nbrves; save_id++)
        {
            Arteries[ save_id] -> printAllt (art_ALL, tend, Arteries[save_id]->N, WALL_MODEL);
        }
        // Arteries[0] ->printFt (wssfile, tend, Arteries[0]->N/2, 2);
        fflush(art_ALL);

        // The time within each print is increased.
        tstart = tend;
        tend   = tend + Deltat; // The current ending time is increased by Deltat.
    }
    
    // In order to termate the program correctly the vessel network and hence
    // all the vessels and their workspace are deleted.
    for (int i=0; i<nbrves; i++) delete Arteries[i];
    
    // Matrices and arrays are deleted
    for (int i=0; i<24; i++) delete[] fjac[i];
    
    // Add additional close statements if you are writing to multiple files
    fclose (art_ALL);
    
    return 0;
}
