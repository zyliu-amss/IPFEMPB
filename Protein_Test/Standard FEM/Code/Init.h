/* Initialize the constants and VAR. Read data from *.pqr */

#ifndef __INIT_H__
#define __INIT_H__

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#define r2(x, y, z, x0, y0, z0) (sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0)))

/************************************************************
* INITIALIZE THE CONSTANTS AND VARs
*************************************************************/
const FLOAT Na = 6.02214129e23;		/* 1mol = NA */
const FLOAT e = 1.6021765e-19;		/*elementary charge in unit of C*/
const FLOAT kB = 1.3806488e-23;		/* Boltzmann constant J/K */
const FLOAT e0 = 8.8541878e-12;		/* vacuum permittivity C^2/(J*K) */
FLOAT T;				/* temperature */
FLOAT e_m;				/* epsilon in \OMEGA_m */
FLOAT e_s;				/* epsilon in \OMEGA_s */
FLOAT e_e;              /* epsilon in medium */
FLOAT alp;				/*relaxation coefficient:alp=1.0 represent totally the new result*/
FLOAT tol;				/*tolerance*/
FLOAT K1, K2;				/* K1=e^2*\beta/\epsilon_0 in unit of A; K2= K1*mol/L in unit of A^-2 */ 
FLOAT kBT;				/* kB*T in unit of kcal/mol */
FLOAT z1,z2,energy;
FLOAT Kappa;            /*Debye-Huckel (2*c_b*ec*ec/(KB*T*e0*e_s))^(1/2) in unit of A^(-1) */
int MNATOMS = 40000;
int NATOMS;				/* The number of atoms */
int NION;				/* The number of ion species*/
//FLOAT p_bulk;
//FLOAT *a;				/* ion sizes*/
//int which_ion;			/* ion no. */
int BC_MAP[2];

FLOAT num, delta;                       /*num is the maximum radius of sv func, delta is the unit step*/

/* ION STRUCT */
typedef struct{
	FLOAT q; 			/* valent */
	FLOAT c; 			/* concentration */
	FLOAT D; 			/* diffusion coefficient */
	BOOLEAN reactive;
}ions;

ions *ion;
//ions  ion[]={{+1, 0.01, 0.196, FALSE}, {+2, 0.01, 0.203, FALSE}, {+3, 0.01, 0.203, FALSE}, {-1, 0.06, 0.203, FALSE}};

/* ATOM STRUCT */
typedef struct{
	FLOAT x, y, z; /* 3D coefficient */
	FLOAT q; /* value of atom = q * ec */
	FLOAT r; /* radius */
}atoms;

atoms *atom;
FLOAT energy;
INT level;
#endif
