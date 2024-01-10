#ifndef __PBE_COEF_H__
#define __PBE_COEF_H__

char *fn_mesh = "sphere.mesh";
char *fn_pqr = "1atm_02_08.pqr";
char *fn_vtk = "PBE.vtk";
extern FLOAT C_b = 0.0;
//char *fn_mesh = "diALAn.mesh";
//char *fn_pqr = "diALAn.pqr";

static void PBE_coef(){
	T = 298.15; 			         /* temperature */
	e_m = 2.0; 			         /* epsilon in \OMEGA_m */
	e_s = 80.0;			         /* epsilon in \OMEGA_s */
	e_e = 2.0;                  /* epsilon in medium */
	alp = 0.9;
    z1 = -15.0;
    z2 = 15.0;
	tol = 1e-6;   				 /*tolerance*/
	K1 = e * e / (kB * T * e0) * 1.0e10;     /* m to A */
	K2 = K1 * Na * 1.0e-27;			 /* A^-2*/
	kBT = kB * T * Na / 4.184e3;
	BC_MAP[0] = 5;				 /* mesh file boundary condition:DIRICHLET, if -1, no marking*/
	BC_MAP[1] = 9;				 /* mesh file boundary condition: BDRY_USER1, if -1, no marking*/

	NION = 2;
	/*initialize ion conditions*/
//	if(!(ion = (struct ions *)phgAlloc(NION * sizeof(*ion))))
	if(!(ion = phgAlloc(NION * sizeof(*ion))))
		phgError(1, "Error in memory allocation for ions.\n");
	ion[0].q = +1;
	ion[0].c = C_b;
	ion[0].D = 0.196;
	ion[0].reactive = FALSE;
	ion[1].q = -1;
	ion[1].c = C_b;
	ion[1].D = 0.203;
	ion[1].reactive = FALSE;


	int i;
	phgPrintf("Number of ion species              %d\n",NION);
	FLOAT Is=0.0;
	for (i = 0; i < NION; i++){
		Is += ion[i].c * ion[i].q * ion[i].q;
		phgPrintf("ion[%d].q  %+4.2f e    ion[%d].c  %f M\n",i,ion[i].q,i,ion[i].c);
	}
	Is *= 0.5;
	Kappa = Sqrt(2.0 * Is * K1 * 1.0e-10 * Na * 1000 / e_s);
	phgPrintf("Ionic strength Is                  %lf M\n\n",(double)Is); 
	phgPrintf("Avogadro constant                  %.4e M-1\n",(double)Na);
	phgPrintf("Elementary charge                  %.4e C\n",(double)e);
	phgPrintf("Boltzmann constant                 %.4e J/K\n",(double)kB);
	phgPrintf("vacuum permittivity                %.4e C^2/(J*K)\n",(double)e0);
	phgPrintf("K1=e*e/(kB*T*e0)*1.0e10            %.4e A\n", (double)K1);
	phgPrintf("K2=K1*Na*1.0e-27                   %.4e A-2\n", (double)K2);
	phgPrintf("kBT=kB*T*Na/4.184e3                %.4e kcal/mol\n", (double)kBT);
	phgPrintf("tolerance                          %.4e\n",(double)tol);
	phgPrintf("Temperature                        %f K\n",(double)T);
	phgPrintf("biomolecule relative permittivity  %f\n",(double)e_m);
	phgPrintf("solvent relative permittivity      %f\n",(double)e_s);
	phgPrintf("relaxation alp                     %f\n",(double)alp);
	phgPrintf("Debye-Huckel                       %.4e A-1\n", (double)Kappa);
//	phgPrintf("The maximum radius                 %f A\n", num);
//	phgPrintf("The step length                    %f A\n", delta);   
}
#endif
