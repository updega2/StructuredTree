/***************************************************************************/
/*                                                                         */
/*  Program: sor06.h                                                       */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This header file defines the global parameters                         */
/*                                                                         */
/***************************************************************************/

// $Id: sor06.h,v 1.7 2005/07/07 22:19:51 heine Exp $

#ifndef _SOR06_H
#define _SOR06_H

#include <cmath>
#include "tools.h"

int    nbrves, N_aorta;               // Number of vessels in the tree.

int    tmstps = 8192,                 // The number of timesteps per period.
       plts   = 1024;                 // Number of plots per period.

char*  CO_filename = "sor06.dat";     // Input flow file at the heart.

double conv   = 1332.20,              // Conversion from mmHg to SI-units.
       rho    = 1.055,                // Density of blood [g/cm^3].
       mu     = 0.049,                // Viscosity of blood [g/cm/s].
       mu_pl  = mu,                   // Viscosity of blood [g/cm/s].
       nu     = mu/rho,

       Tper   = 1.10,                 // The period of one heart beat [s].

       ff1    = 0.999639e07,
       ff2    = -11.2635,
       ff3    = 432631,
       //ff1  = 1.99925e+07
       //ff2  = -22.5267,
       //ff3  =  865251,
       //ff1  = 1.99925e07,
       //ff2  = -25.5267,
       //ff3  = 465251,
       // GUESS FROM REVIEW PAPER
       //ff1 = 0.999639e07,
       //ff2 = -5.2635,
       //ff3 = 6432631,

       m1     = 7.3693,
       m2     = 1.2122,
       m3     = 5.6517,
       m4     = 0.21763,

       //Fcst = 10.0,                 // Determines the damping coeff.
       Fcst   = 50, //17.7778,        // Determines the damping coeff.
                                      // for the friction.
       Lr     = 1.0,                  // Characteristic radius of the
                                      // vessels in the tree [cm].
       Lr2    = sq(Lr),               // The squared radius [cm2].
       Lr3    = cu(Lr),               // The radius to the third power [cm^3].
       g      = 981.0,                // The gravitational force [cm/s^2].
       q      = 10.0*Lr2,             // The characteristic flow [cm^3/s].
       Fr2    = sq(q)/g/pow(Lr,5),    // The squared Froudes number.
       Re     = q*rho/mu/Lr,          // Reynolds number.
       Period = Tper*q/Lr3,           // The dimension-less period.
       tau    = 0.08*q/Lr3,           // End of pulse, dimension-less.
       k      = Period/tmstps,        // Length of a timestep.
       Deltat = Period/plts,          // Interval between each point plottet.
       p0     = 55.0/rho/g/Lr*conv,   // Ensures a certain diastolic pressure.

       *fjac[18],   // Work space used by bound_bif.
       xr, f, df;                     // Work space used by bound_right.

#endif
