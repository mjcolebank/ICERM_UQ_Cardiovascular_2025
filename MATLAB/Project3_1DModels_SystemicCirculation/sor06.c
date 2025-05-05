/* The sor06.C main program */

// $Id: sor06.C,v 1.13 2010-10-20 15:38:28 mette Exp $

#include "sor06.h"
#include "tools.h"
#include "arteries.h"
#include "junction.h"

extern "C"  void impedance_init_driver_(int *tmstps);

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>

using namespace std;

// The vessel-network is specified and initialized, and the flow and
// pressures are to be determined. All global constants must be defined
// in the header file. That is the number of dots per cm, and the number
// of vessels.
int main(int argc, char *argv[])
{
    double tstart, tend, finaltime;
    // Declare parameters: AML 01/27/23
    double f1, f2, f3, fs1, fs2, fs3;
    double STpar_1,STpar_2, lrr, rmScaling;
    //double k1r, k2r, k3r;
    int total_vessels, total_terminal, num_pts;
    int ST_calc;
    int fileID;
    
    double ves_angle; // NOTE: Setting this to 1 will use the
    // ST parameters that have been previously generated
    
    if (argc < 7) {
        fprintf(stdout,"Not enough entries: Only %d inputs. Exiting now.\n",argc);
        return 0;
    }
    ////// Load in fluids parameters
    // Known values
     f1      = 200000000;
     fs1     = 200000000;
     STpar_2 = 0.60;
     lrr     = 50.0;
     //rm      = 0.001;    
     
     rmScaling = 1;
     
    // Unknown
     f2        = atof(argv[1]);
     f3        = atof(argv[2]);
     fs2       = atof(argv[3]);
     fs3       = atof(argv[4]);
     STpar_1   = atof(argv[5]);
     //rmScaling = atof(argv[6]); // global scaling factor for rmin
     fileID = atoi(argv[6]);
     
    int cycles = 1; // NOTE: Can make this a parameter to pass in
    
    // Load in network parameters
    total_vessels  = 9; // atoi(argv[1]);
    total_terminal = 5; // atoi(argv[2]);
    num_pts        = 10; // atoi(argv[3]);
    ves_angle      = 0;
    //total_conn     = total_vessels-total_terminal;
    nbrves         = total_vessels;
    
    // A parameter that can alleviate calculating the impedance for each iteration
    ST_calc = 1;//atoi(argv[14]);
  
   // mihaela added this
   char namepuALL[20];
    sprintf(namepuALL, "output_%d.2d", fileID);
	FILE *fpALL = fopen (namepuALL, "w");
    
    // Workspace used by bound_bif
    for(int i=0; i<24; i++) fjac[i] = new double[24];
    
    tstart    = 0.0;            // Starting time.
    
    // The number of vessels in the network is given when the governing array of
    // vessels is declared.
    
    impedance_init_driver_(&tmstps);
    Tube   *Arteries[nbrves];                    // Array of blood vessels.
        
    int conn_rows = 2*total_vessels; // Initialize with double the number of vessels
    int conn_cols = max_D; // This is defined so that we can have trifurcations if needed.
    int connectivity_matrix[conn_rows][conn_cols];
    int terminal_vessels[total_terminal];
    
        // Begin sequence of loading text files
    FILE *conn;
    conn = fopen("connectivity.txt","rt");
    int parent, daughter1, daughter2, daughter3, r_in;
    
    // MJC: Define boolean vector for converging flow
    //int conv_flag[total_vessels];
   // for (int i=0; i<total_vessels; i++)
    //    conv_flag[i]=0;
    
    // Check to see if we have the connectivity file
    int conn_id = 0;
    if (conn == NULL)
    {
        fprintf(stdout,"Error: Connectivity File Does Not Exist \n");
        return 0;
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
    size_conn = conn_id-1;

    
    // Load in terminal vessels
    
    FILE *terminal_file;
    terminal_file = fopen("terminal_vessels.txt","rt");
    for (int i=0; i<total_terminal; i++){
        fscanf(terminal_file, "%d", &terminal_vessels[i]);
    }
    fclose(terminal_file);
    
    
    /* Initialization of the Arteries.
     * // The tube class takes in the following values (in this order)
     * Length: the length of the artery
     *
     * Top radius: the radius at the entry of the vessel
     *
     * Bottom radius: the radius at the end of the vessel
     *
     * Daughter Pointer: the pointer to the connectivity of all the blood
     * vessels, passed to junction.c
     *
     * Minimum radius: minimum radius (in ST)
     *
     * Num Points: the number of pts (per non-dimensional length) to use in numerical solution
     *
     * Init: Set to 1 to make this vessel the inlet vessel (i.e., it gets a flow prescribed); set 0 else
     *
     * f1,f2,f3: Stiffness in the large blood vessels (Eh/r0 = f1*exp(f2*r0)+f3)
     *
     * fs1,fs2,fs3: Stiffness in the small blood vessels (Eh/r0 = fs1*exp(fs2*r0)+fs3) (in ST)
     *
     * STpar_1: the first parameter describing microvascular branching (could be alpha or asym) (in ST)
     *
     * STpar_2: the second parameter describing microvascular branching (could be beta or expo) (in ST)
     *
     * Length-radius ratio: the third parameter for microvasculature; relates length to radius (in ST)
     *
     * Terminal-ID: ID of the terminal vessel, which is utilized when simultations are run without
     * recalculating the structured tree parameters (in ST)
     *
     * ST_calc: A flag to tell the code whether to recalculate the ST parameters (set to 0 to not calculate ST value)
     *Definition of Class Tube: (Length, topradius, botradius, *LeftDaughter, *RightDaughter, points, init, K, f1, f2, f3, R1, R2,  CT);*/
       int *daughter_ptr = &connectivity_matrix[0][0];
    
        Arteries[8] = new Tube( 17.22, 0.25, 0.25, 0, 0, 0.25, daughter_ptr, rmScaling*0.001, num_pts, 0, f1,f2,f3,fs1,fs2,fs3, STpar_1, STpar_2, lrr,5,ST_calc,ves_angle);
        
        Arteries[7] = new Tube( 20.23, 0.35, 0.30, 0, 0, 0.30, daughter_ptr, rmScaling*0.05, num_pts, 0, f1,f2,f3,fs1,fs2,fs3, STpar_1, STpar_2, lrr,4,ST_calc,ves_angle);
        
        Arteries[6] = new Tube( 15.17, 0.88, 0.59, 0, 0, 0.59, daughter_ptr, rmScaling*0.08, num_pts, 0, f1,f2,f3,fs1,fs2,fs3, STpar_1, STpar_2, lrr,3,ST_calc,ves_angle);
        
        Arteries[5] = new Tube( 3.67, 0.59, 0.37, 0, 0, 0.37, daughter_ptr, 0, num_pts, 0, f1,f2,f3,0,0,0, 0, 0, lrr,0,ST_calc,ves_angle);
        
        Arteries[4] = new Tube( 20.23, 0.25, 0.25, 0, 0, 0.25, daughter_ptr, rmScaling*0.0005, num_pts, 0, f1,f2,f3,fs1,fs2,fs3, STpar_1, STpar_2, lrr,2,ST_calc,ves_angle);
        
        Arteries[3] = new Tube( 1.94, 0.95, 0.88, 0, 0, 0.88, daughter_ptr, 0, num_pts, 0, f1,f2,f3,0,0,0, 0, 0, 0,0,ST_calc,ves_angle);
        
        Arteries[2] = new Tube( 1.23, 0.25, 0.25, 0, 0, 0.25, daughter_ptr, rmScaling*0.08, num_pts, 0, f1,f2,f3,fs1,fs2,fs3, STpar_1, STpar_2, lrr,1,ST_calc,ves_angle);
       
        Arteries[1] = new Tube( 1.95, 1.10, 0.95, 0, 0, 0.95, daughter_ptr, 0, num_pts, 0, f1,f2,f3,0,0,0, 0, 0, 0,0,ST_calc,ves_angle);
        
        Arteries[0] = new Tube( 4.07, 1.20, 1.10, 0, 0, 1.10,
                daughter_ptr, 0, num_pts, 1,f1,f2,f3,0,0,0, 0, 0, 0,0,ST_calc,ves_angle);
    
    
    // Solves the equations until time equals tend.///////////////////////////
    /* ADDED BY MJ Colebank
     * Rather than specifying the number of cycles as an input to the function,
     * we want to test to see if the solution has converged. If so, we should exit.*/

    int period_counter = 1; // Count the number of periods you have solved for
    double norm_sol = 1e+6;
    double sol_tol  = 1e-5;
//    printf("NORM_SOL: %f\n",norm_sol);
    double sol_p1[tmstps],sol_p2[tmstps];
    tend      = Deltat;
    
    
    // SOLVE THE MODEL ONCE
    // Note: Only want to test the pressure at the inlet
    int sol_ID = 0;
    while (tend<=period_counter*Period)
    {
        solver (Arteries, tstart, tend, k, WALL_MODEL);
        sol_p1[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0],WALL_MODEL); // for printing
        sol_p1[sol_ID] *= rho*g*Lr/conv;
        tstart = tend;
        tend   = tend + Deltat; // The current ending time is increased by Deltat.
        sol_ID++;
    }
    
    
    // LOOP FOR CONVERGENCE
    double sse;
    while (norm_sol>=sol_tol)
    {
        sol_ID = 0;
        sse    = 0;
        period_counter++;
        while (tend<=period_counter*Period)
        {
            solver (Arteries, tstart, tend, k, WALL_MODEL);
            sol_p2[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0],WALL_MODEL); // for printing
            sol_p2[sol_ID] *= rho*g*Lr/conv;
            sse = sse+ sq(sol_p1[sol_ID]-sol_p2[sol_ID]);
            tstart = tend;
            tend   = tend + Deltat; // The current ending time is increased by Deltat.
            sol_ID++;
        }
        norm_sol = sse;
        memcpy (sol_p1, sol_p2, sizeof(sol_p2));
       //printf("NORM_SOL:%f\n",norm_sol);
    }
//    printf("num_cylces:%d\n",period_counter);
    
    
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
        solver (Arteries, tstart, tend, k, WALL_MODEL);
        
        
        for (int save_id=0; save_id<nbrves; save_id++)
        {
            Arteries[ save_id] -> printAllt (fpALL, tend, 0,WALL_MODEL);
        }
        for (int save_id=0; save_id<nbrves; save_id++)
        {
            Arteries[ save_id] -> printAllt (fpALL, tend, Arteries[save_id]->N/2,WALL_MODEL);
        }
        for (int save_id=0; save_id<nbrves; save_id++)
        {
            Arteries[ save_id] -> printAllt (fpALL, tend, Arteries[save_id]->N,WALL_MODEL);
        }
        // The time within each print is increased.
        tstart = tend;
        tend   = tend + Deltat; // The current ending time is increased by Deltat.
    }
    
    
    
    // In order to termate the program correctly the vessel network and hence
    // all the vessels and their workspace are deleted.
    for (int i=0; i<nbrves; i++) delete Arteries[i];
    
    // Matrices and arrays are deleted
    for (int i=0; i<18; i++) delete[] fjac[i];
    
    fclose(fpALL);
    
    return 1;
}
