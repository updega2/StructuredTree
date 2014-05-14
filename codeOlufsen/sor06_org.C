/* The sor06.C main program */

// $Id: sor06.C,v 1.10 2005/10/14 18:05:59 heine Exp $

#include "sor06.h"
#include "tools.h"
#include "arteries.h"

extern "C"  void impedance_init_driver_(int *tmstps);

extern "C" void g95_runtime_start(int argc, char *argv[]);
extern "C" void g95_runtime_stop();


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
  g95_runtime_start(argc,argv);
  double tstart, tend, finaltime;

  double rm1  = 0.04;
  double rm12 = 0.03;
  double rm2  = 0.02;
  double rm3  = 0.01;

  double f1   = 1.99925e07;
  double f2   = -25.5267;
  double f3AA = 365251;
  double f3A  = 465251;
  double f3AB = 565251;
  double f3B  = 665251;

  char *nameq1  = "q1.2d";
  char *nameq2  = "q2.2d";
  char *nameq3  = "q3.2d";
  char *nameq4  = "q4.2d";
  char *nameq6  = "q6.2d";
  char *nameq8  = "q8.2d";
  char *nameq10 = "q10.2d";
  char *nameq12 = "q12.2d";
  char *nameq14 = "q14.2d";
  char *nameq16 = "q16.2d";
  char *nameq18 = "q18.2d";
  char *nameq20 = "q20.2d";
  char *nameq22 = "q22.2d";
  char *nameq24 = "q24.2d";
  char *namep1  = "p1.2d";
  char *namep2  = "p2.2d";
  char *namep3  = "p3.2d";
  char *namep4  = "p4.2d";
  char *namep6  = "p6.2d";
  char *namep8  = "p8.2d";
  char *namep10 = "p10.2d";
  char *namep12 = "p12.2d";
  char *namep14 = "p14.2d";
  char *namep16 = "p16.2d";
  char *namep18 = "p18.2d";
  char *namep20 = "p20.2d";
  char *namep22 = "p22.2d";
  char *namep24 = "p24.2d";

  FILE *fp1 = fopen (namep1, "w");
  if (fp1) fprintf(stdout, "File 1p OK \n"); else error ("main.C","File 1p NOT OK");
  FILE *fp2 = fopen (namep2, "w");
  if (fp2) fprintf(stdout, "File 2p OK \n"); else error ("main.C","File 2p NOT OK");
  FILE *fp3 = fopen (namep3, "w");
  if (fp3) fprintf(stdout, "File 3p OK \n"); else error ("main.C","File 3p NOT OK");
  FILE *fp4 = fopen (namep4, "w");
  if (fp4) fprintf(stdout, "File 4p OK \n"); else error ("main.C","File 4p NOT OK");
  FILE *fp6 = fopen (namep6, "w");
  if (fp6) fprintf(stdout, "File 6p OK \n"); else error ("main.C","File 6p NOT OK");
  FILE *fp8 = fopen (namep8, "w");
  if (fp8) fprintf(stdout, "File 8p OK \n"); else error ("main.C","File 8p NOT OK");
  FILE *fp10 = fopen (namep10, "w");
  if (fp10) fprintf(stdout, "File 10p OK \n"); else error ("main.C","File 10p NOT OK");
  FILE *fp12 = fopen (namep12, "w");
  if (fp12) fprintf(stdout, "File 12p OK \n"); else error ("main.C","File 12p NOT OK");
  FILE *fp14 = fopen (namep14, "w");
  if (fp14) fprintf(stdout, "File 14p OK \n"); else error ("main.C","File 14p NOT OK");
  FILE *fp16 = fopen (namep16, "w");
  if (fp16) fprintf(stdout, "File 16p OK \n"); else error ("main.C","File 16p NOT OK");
  FILE *fp18 = fopen (namep18, "w");
  if (fp18) fprintf(stdout, "File 18p OK \n"); else error ("main.C","File 18p NOT OK");
  FILE *fp20 = fopen (namep20, "w");
  if (fp20) fprintf(stdout, "File 20p OK \n"); else error ("main.C","File 20p NOT OK");
  FILE *fp22 = fopen (namep22, "w");
  if (fp22) fprintf(stdout, "File 22p OK \n"); else error ("main.C","File 22p NOT OK");
  FILE *fp24 = fopen (namep24, "w");
  if (fp24) fprintf(stdout, "File 24p OK \n"); else error ("main.C","File 24p NOT OK");

  FILE *fq1 = fopen (nameq1, "w");
  if (fq1) fprintf(stdout, "File 1q OK \n"); else error ("main.C","File 1q NOT OK");
  FILE *fq2 = fopen (nameq2, "w");
  if (fq2) fprintf(stdout, "File 2q OK \n"); else error ("main.C","File 2q NOT OK");
  FILE *fq3 = fopen (nameq3, "w");
  if (fq3) fprintf(stdout, "File 3q OK \n"); else error ("main.C","File 3q NOT OK");
  FILE *fq4 = fopen (nameq4, "w");
  if (fq4) fprintf(stdout, "File 4q OK \n"); else error ("main.C","File 4q NOT OK");
  FILE *fq6 = fopen (nameq6, "w");
  if (fq6) fprintf(stdout, "File 6q OK \n"); else error ("main.C","File 6q NOT OK");
  FILE *fq8 = fopen (nameq8, "w");
  if (fq8) fprintf(stdout, "File 8q OK \n"); else error ("main.C","File 8q NOT OK");
  FILE *fq10 = fopen (nameq10, "w");
  if (fq10) fprintf(stdout, "File 10q OK \n"); else error ("main.C","File 10q NOT OK");
  FILE *fq12 = fopen (nameq12, "w");
  if (fq12) fprintf(stdout, "File 12q OK \n"); else error ("main.C","File 12q NOT OK");
  FILE *fq14 = fopen (nameq14, "w");
  if (fq14) fprintf(stdout, "File 14q OK \n"); else error ("main.C","File 14q NOT OK");
  FILE *fq16 = fopen (nameq16, "w");
  if (fq16) fprintf(stdout, "File 16q OK \n"); else error ("main.C","File 16q NOT OK");
  FILE *fq18 = fopen (nameq18, "w");
  if (fq18) fprintf(stdout, "File 18q OK \n"); else error ("main.C","File 18q NOT OK");
  FILE *fq20 = fopen (nameq20, "w");
  if (fq20) fprintf(stdout, "File 20q OK \n"); else error ("main.C","File 20q NOT OK");
  FILE *fq22 = fopen (nameq22, "w");
  if (fq22) fprintf(stdout, "File 22q OK \n"); else error ("main.C","File 22q NOT OK");
  FILE *fq24 = fopen (nameq24, "w");
  if (fq24) fprintf(stdout, "File 24q OK \n"); else error ("main.C","File 24q NOT OK");

  // Workspace used by bound_right

  // Workspace used by bound_bif
  for(int i=0; i<18; i++) fjac[i] = new double[18];

  clock_t c1 = clock();        // Only used when timing the program.
  nbrves    = 24;             // The number of vessels in the network.
  tstart    = 0.0;            // Starting time.
  finaltime = 4*Period;       // Final end-time during a simulation.
  tend      = 3*Period;       // Timestep before the first plot-point
                              // is reached.

  // The number of vessels in the network is given when the governing array of
  // vessels is declared.

  impedance_init_driver_(&tmstps);

  Tube   *Arteries[nbrves];                    // Array of blood vessels.
  // Arteries = new Tube*[nbrves];

  // Initialization of the vessel-network.
  // Right Radial: Estimate.

  // Initialization of the vessel-network.
  // DEBANJAN: WE HAVE ADDED A LAST TERM TO SIMULATE A TERMINAL RESISTANCE TO ALL OF THESE CONSTRUTORS
  Arteries[23] = new Tube(44.0,0.3973,0.300,0,0,rm3,4,0,0,f1,f2,f3AB, 0.0);
  // Femoral: Estimate.
  Arteries[22] = new Tube(11.0,0.200,0.200,0,0,rm3,4,0,0,f1,f2,f3A, 0.0);
  // Deep Femoral: Estimate.
  Arteries[21] = new Tube( 4.5,0.200,0.200,0,0,rm3,4,0,0,f1,f2,f3A, 0.0);
  // Internal Iliac: Pure estimate.
  Arteries[20] = new Tube(13.0,0.4317,0.3973,Arteries[22],Arteries[23],0,4,0,0,f1,f2,f3AB, 0.0);
  // Femoral artery:
  Arteries[19] = new Tube( 6.5,0.45,0.4317,Arteries[20],Arteries[21],0,4,0,0,f1,f2,f3A, 0.0);
  // External Iliac:
  Arteries[18] = new Tube( 3,0.7292,0.700, Arteries[19],Arteries[19],0,4,0,0,f1,f2,f3A, 0.0);
  // Abdominal Aorta:
  Arteries[17] = new Tube( 4.0,0.200,0.175,0,0,rm3,4,0,0,f1,f2,f3B, 0.0);
  // Inf Mesenteric: Pure estimate.
  Arteries[16] = new Tube( 6.0,0.7914,0.7292, Arteries[17],Arteries[18],0,4,0,0,f1,f2,f3A, 0.0);
  // Abdominal Aorta:
  Arteries[15] = new Tube( 3.0,0.275,0.25,0,0,rm2,4,0,0,f1,f2,f3B, 0.0);
  // Left Renal: Location and area should be estimated.
  Arteries[14] = new Tube( 1.0,0.8023,0.7914, Arteries[15],Arteries[16],0,4,0,0,f1,f2,f3A, 0.0);
  // Abdominal Aorta:
  Arteries[13] = new Tube( 3.0,0.275,0.25,0,0,rm2,4,0,0,f1,f2,f3B, 0.0);
  // Right Renal: Location and area should be estimated.
  Arteries[12] = new Tube( 2.0,0.8245,0.8023, Arteries[13],Arteries[14],0,4,0,0,f1,f2,f3A, 0.0);
  // Abdominal Aorta:
  Arteries[11] = new Tube( 5.0,0.325,0.325,0,0,rm1,4,0,0,f1,f2,f3B, 0.0);
  // Sup. Mes: Area of should be estimated.
  Arteries[10] = new Tube( 2.0,0.8473,0.8245, Arteries[11],Arteries[12],0,4,0,0,f1,f2,f3A, 0.0);
  // Descending Aorta: Cut to allow for leak from intercostal arteries.
  Arteries[ 9] = new Tube( 3.0,0.325,0.300,0,0,rm2,4,0,0,f1,f2,f3B, 0.0);
  // Celiac Axis: Location and area should be estimated.
  Arteries[ 8] = new Tube(18.75,1.1093,0.8473, Arteries[ 9],Arteries[10],0,4,0,0,f1,f2,f3A, 0.0);
  // Descending Aorta: Cut to allow for leak from intercostal arteries.
  Arteries[ 7] = new Tube(43.0,0.445,0.165,0,0,rm3,4,0,0,f1,f2,f3A, 0.0);
  // Left Sub: Length and bot radius estimated.
  Arteries[ 6] = new Tube( 1.0,1.1093,1.10943, Arteries[ 7],Arteries[ 8],0,4,0,0,f1,f2,f3A, 0.0);
  // Aorta des:
  Arteries[ 5] = new Tube(19.0,0.285,0.275, 0,0,rm12,4,0,0,f1,f2,f3AA, 0.0);
  // Left Carotid: Length and bot radius estimated.
  Arteries[ 4] = new Tube( 1.75,1.1361,1.1093, Arteries[ 5],Arteries[ 6],0,4,0,0.75,f1,f2,f3A, 0.0);
  // Aortic arc:
  Arteries[ 3] = new Tube(17.0,0.285,0.275, 0,0,rm12,4,0,0,f1,f2,f3AA, 0.0);
  // Right Carotid: Length and bot radius estimated.
  Arteries[ 2] = new Tube(43.0,0.445,0.165,0,0,rm3,4,0,0,f1,f2,f3A, 0.0);
  // Right Sub: Data taken from left sub, length and bot radius  estimated !
  Arteries[ 1] = new Tube( 3.5,0.700,0.700, Arteries[ 2],Arteries[ 3],0,4,0,0,f1,f2,f3AA, 0.0);
  // Anonymous: Length estimated - should be measured.
  Arteries[ 0] = new Tube( 7.0,1.250,1.1361, Arteries[ 1],Arteries[ 4],0,4,1,0,f1,f2,f3A, 0.0);
  // Ascending Aorta:

  // In the next three statements the simulations are performed until
  // tstart = tend. That is this is without making any output during this
  // first period of time. If one needs output during this period, these three
  // lines should be commented out, and the entire simulation done within the
  // forthcomming while-loop.

  // Solves the equations until time equals tend.
  solver (Arteries, tstart, tend, k);
  tstart = tend;
  tend = tend + Deltat;

  // fprintf (stdout,"saves Q0\n");
  // Arteries[ 0] -> printQ0(fq0);

  fprintf (stdout,"plots start\n");

  // The loop is continued until the final time
  // is reached. If one wants to make a plot of
  // the solution versus x, tend is set to final-
  // time in the above declaration.
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
    solver (Arteries, tstart, tend, k);
    fprintf (stdout,".");

    // A 2D plot of P(x_fixed,t) is made. The resulting 2D graph is due to
    // the while loop, since for each call of printPt only one point is set.
    Arteries[ 0] -> printPxt (fp1, tend, 0);
    Arteries[ 4] -> printPxt (fp1, tend, int( 7.00*Arteries[ 4]->pts));
    Arteries[ 6] -> printPxt (fp1, tend, int( 8.75*Arteries[ 6]->pts));
    Arteries[ 8] -> printPxt (fp1, tend, int( 9.75*Arteries[ 8]->pts));
    Arteries[10] -> printPxt (fp1, tend, int(28.50*Arteries[10]->pts));
    Arteries[12] -> printPxt (fp1, tend, int(30.50*Arteries[12]->pts));
    Arteries[14] -> printPxt (fp1, tend, int(31.50*Arteries[14]->pts));
    Arteries[16] -> printPxt (fp1, tend, int(32.50*Arteries[16]->pts));
    Arteries[18] -> printPxt (fp1, tend, int(38.50*Arteries[18]->pts));

    Arteries[ 0] -> printQxt (fq1, tend, 0);
    Arteries[ 4] -> printQxt (fq1, tend, int( 7.00*Arteries[ 4]->pts));
    Arteries[ 6] -> printQxt (fq1, tend, int( 8.75*Arteries[ 6]->pts));
    Arteries[ 8] -> printQxt (fq1, tend, int( 9.75*Arteries[ 8]->pts));
    Arteries[10] -> printQxt (fq1, tend, int(28.50*Arteries[10]->pts));
    Arteries[12] -> printQxt (fq1, tend, int(30.50*Arteries[12]->pts));
    Arteries[14] -> printQxt (fq1, tend, int(31.50*Arteries[14]->pts));
    Arteries[16] -> printQxt (fq1, tend, int(32.50*Arteries[16]->pts));
    Arteries[18] -> printQxt (fq1, tend, int(38.50*Arteries[18]->pts));

    Arteries[19] -> printPxt (fp20, tend, 0);
    Arteries[20] -> printPxt (fp20, tend, int( 6.5*Arteries[20]->pts));
    Arteries[22] -> printPxt (fp20, tend, int(19.5*Arteries[22]->pts));

    Arteries[19] -> printQxt (fq20, tend, 0);
    Arteries[20] -> printQxt (fq20, tend, int( 6.5*Arteries[20]->pts));
    Arteries[22] -> printQxt (fq20, tend, int(19.5*Arteries[22]->pts));

    Arteries[ 1] -> printPxt (fp2, tend, 0);
    Arteries[ 1] -> printQxt (fq2, tend, 0);

    Arteries[ 2] -> printPxt (fp3, tend, 0);
    Arteries[ 2] -> printQxt (fq3, tend, 0);

    Arteries[ 3] -> printPxt (fp4, tend, 0);
    Arteries[ 3] -> printQxt (fq4, tend, 0);

    Arteries[ 5] -> printPxt (fp6, tend, 0);
    Arteries[ 5] -> printQxt (fq6, tend, 0);

    Arteries[ 7] -> printPxt (fp8, tend, 0);
    Arteries[ 7] -> printQxt (fq8, tend, 0);

    Arteries[ 9] -> printPxt (fp10, tend, 0);
    Arteries[ 9] -> printQxt (fq10, tend, 0);

    Arteries[11] -> printPxt (fp12, tend, 0);
    Arteries[11] -> printQxt (fq12, tend, 0);

    Arteries[13] -> printPxt (fp14, tend, 0);
    Arteries[13] -> printQxt (fq14, tend, 0);

    Arteries[15] -> printPxt (fp16, tend, 0);
    Arteries[15] -> printQxt (fq16, tend, 0);

    Arteries[17] -> printPxt (fp18, tend, 0);
    Arteries[17] -> printQxt (fq18, tend, 0);

    Arteries[21] -> printPxt (fp22, tend, 0);
    Arteries[21] -> printQxt (fq22, tend, 0);

    Arteries[23] -> printPxt (fp24, tend, 0);
    Arteries[23] -> printQxt (fq24, tend, 0);


    // The time within each print is increased.
    tstart = tend;
    tend   = tend + Deltat; // The current ending time is increased by Deltat.
  }
  fprintf(stdout,"\n");

  // The following statements is used when timing the simulation.
  fprintf(stdout,"nbrves = %d, Lax, ", nbrves);
  clock_t c2 = clock(); // FIXME clock() may wrap after about 72 min.
  int tsec = (int) ((double) (c2-c1)/CLOCKS_PER_SEC + 0.5);
  fprintf(stdout,"cpu-time %d:%02d\n", tsec / 60, tsec % 60);
  fprintf(stdout,"\n");

  // In order to termate the program correctly the vessel network and hence
  // all the vessels and their workspace are deleted.
  for (int i=0; i<nbrves; i++) delete Arteries[i];

  // Matrices and arrays are deleted
  for (int i=0; i<18; i++) delete[] fjac[i];

  if (fclose (fp1)  != EOF) fprintf(stdout,"Close 1p OK\n");
    else error("main.C","Close 1p NOT OK");
  if (fclose (fp2)  != EOF) fprintf(stdout,"Close 2p OK\n");
    else error("main.C","Close 2p NOT OK");
  if (fclose (fp3)  != EOF) fprintf(stdout,"Close 3p OK\n");
    else error("main.C","Close 3p NOT OK");
  if (fclose (fp4)  != EOF) fprintf(stdout,"Close 4p OK\n");
    else error("main.C","Close 4p NOT OK");
  if (fclose (fp6)  != EOF) fprintf(stdout,"Close 6p OK\n");
    else error("main.C","Close 6p NOT OK");
  if (fclose (fp8)  != EOF) fprintf(stdout,"Close 8p OK\n");
    else error("main.C","Close 8p NOT OK");
  if (fclose (fp10)  != EOF) fprintf(stdout,"Close 10p OK\n");
    else error("main.C","Close 10p NOT OK");
  if (fclose (fp12)  != EOF) fprintf(stdout,"Close 12p OK\n");
    else error("main.C","Close 12p NOT OK");
  if (fclose (fp14)  != EOF) fprintf(stdout,"Close 14p OK\n");
    else error("main.C","Close 14p NOT OK");
  if (fclose (fp16)  != EOF) fprintf(stdout,"Close 16p OK\n");
    else error("main.C","Close 16p NOT OK");
  if (fclose (fp18)  != EOF) fprintf(stdout,"Close 18p OK\n");
    else error("main.C","Close 18p NOT OK");
  if (fclose (fp20)  != EOF) fprintf(stdout,"Close 20p OK\n");
    else error("main.C","Close 20p NOT OK");
  if (fclose (fp22)  != EOF) fprintf(stdout,"Close 22p OK\n");
    else error("main.C","Close 22p NOT OK");
  if (fclose (fp24)  != EOF) fprintf(stdout,"Close 24p OK\n");
    else error("main.C","Close 24p NOT OK");

  if (fclose (fq1)  != EOF) fprintf(stdout,"Close 1q OK\n");
    else error("main.C","Close 1q NOT OK");
  if (fclose (fq2)  != EOF) fprintf(stdout,"Close 2q OK\n");
    else error("main.C","Close 2q NOT OK");
  if (fclose (fq3)  != EOF) fprintf(stdout,"Close 3q OK\n");
    else error("main.C","Close 3q NOT OK");
  if (fclose (fq4)  != EOF) fprintf(stdout,"Close 4q OK\n");
    else error("main.C","Close 4q NOT OK");
  if (fclose (fq6)  != EOF) fprintf(stdout,"Close 6q OK\n");
    else error("main.C","Close 6q NOT OK");
  if (fclose (fq8)  != EOF) fprintf(stdout,"Close 8q OK\n");
    else error("main.C","Close 8q NOT OK");
  if (fclose (fq10)  != EOF) fprintf(stdout,"Close 10q OK\n");
    else error("main.C","Close 10q NOT OK");
  if (fclose (fq12)  != EOF) fprintf(stdout,"Close 12q OK\n");
    else error("main.C","Close 12q NOT OK");
  if (fclose (fq14)  != EOF) fprintf(stdout,"Close 14q OK\n");
    else error("main.C","Close 14q NOT OK");
  if (fclose (fq16)  != EOF) fprintf(stdout,"Close 16q OK\n");
    else error("main.C","Close 16q NOT OK");
  if (fclose (fq18)  != EOF) fprintf(stdout,"Close 18q OK\n");
    else error("main.C","Close 18q NOT OK");
  if (fclose (fq20)  != EOF) fprintf(stdout,"Close 20q OK\n");
    else error("main.C","Close 20q NOT OK");
  if (fclose (fq22)  != EOF) fprintf(stdout,"Close 22q OK\n");
    else error("main.C","Close 22q NOT OK");
  if (fclose (fq24)  != EOF) fprintf(stdout,"Close 24q OK\n");
    else error("main.C","Close 24q NOT OK");

  g95_runtime_stop();
  return 0;
}
