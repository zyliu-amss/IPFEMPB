#include"phg.h"
#include"Init.h"
#include"PBE_coef.h"
#include"func.h"
#include"PBE.h"
//#include"tool.h"

int main(int argc, char **argv){

	INT mem_max = 20000;
	size_t mem, mem_peak;

	/**************************************************/
        /* define DOFs to be used in the solvation*/    
	/* G, H, u are the three parts of Poisson equation*/
	/* tu is the old value of u */ 
	/* deltau is used to store the increasement of u*/
	/* u2 = phgDofNormL2(u) */

	GRID *g;
	SOLVER *solver;
	DOF *G, *H, *Grad_G, *Grad_H;
	DOF *u, *tu, *deltau;
	FLOAT u2;
	/**************************************************/
    phgOptionsRegisterFilename("fn_mesh", "mesh_file", &fn_mesh);
    phgOptionsRegisterFilename("fn_pqr", "pqr_file", &fn_pqr);
    phgOptionsRegisterFilename("fn_vtk", "vtk_file", &fn_vtk);
	phgOptionsRegisterFloat("tol", "Tolerance", &tol);
	phgOptionsRegisterFloat("-C_b", "bulk C_b", &C_b);
    phgOptionsRegisterInt("level", "Level", &level);
	phgOptionsRegisterInt("mem_max", "Maximum memory (MB)", &mem_max);
	phgOptionsPreset("-mem_max 5120");
    
	phgInit(&argc,&argv);
    PBE_coef();
	Read_pqr(fn_pqr);
      

	g = phgNewGrid(-1);
	phgImportSetBdryMapFunc(bc_map);   
//        phgSetPeriodicity(g, X_MASK|Y_MASK);
	if (!phgImport(g, fn_mesh, FALSE))
		phgError(1, "can't read file \"%s\".\n", fn_mesh);
        phgRefineAllElements(g, level);      

	if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	    	phgPrintf("Repartition mesh, load imbalance: %lg\n", (double)g->lif);	
   phgPrintf("\n!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	/**************************************************/
	/*Initialize G, Grad_G, H*/

	G = phgDofNew(g, DOF_DEFAULT, 1, "G", func_G);
	phgPrintf("||G||_2 = %.4le\n", (double)phgDofNormL2(G));
	Grad_G = phgDofNew(g, DOF_ANALYTIC, 3, "Grad_G", func_gradG);
	H = phgDofNew(g, DOF_DEFAULT, 1, "H", DofInterpolation);
	/**************************************************/

	/**************************************************/
	/*Solve H and get Grad_H*/

	solver = phgSolverCreate(SOLVER_DEFAULT, H, NULL);
	Harmonic_stiff_rhs(solver, H);
	phgSolverSolve(solver, TRUE, H, NULL);
	phgPrintf("||H||_2 = %.4le, ",(double)phgDofNormL2(H));
	phgPrintf("nits = %d,residual = %.4le\n", solver->nits, (double)solver->residual);
	phgSolverDestroy(&solver);		
	Grad_H = phgDofGradient(H, NULL, NULL, NULL);
	phgPrintf("||Grad_H||_2 = %.4e\n\n",(double)phgDofNormL2(Grad_H));
	/**************************************************/
    
	/**************************************************/
	/* Initialize u, tmu, ttu, deltau*/
	//u=phgDofNew(g,DOF_DEFAULT,1,"u",DofInterpolation);
    u=phgDofNew(g,DOF_DEFAULT,1,"u", func_uc);
	tu = phgDofNew(g, DOF_DEFAULT, 1, "tu", DofInterpolation);
	deltau=phgDofNew(g,DOF_DEFAULT, 1, "deltau",DofInterpolation);
	
//	phgDofSetDataByValue(u, 0.0);
	phgDofCopy(u, &deltau, NULL, NULL);
	phgDofCopy(u, &tu, NULL, NULL);
	/**************************************************/

	/**************************************************/
	/* define dof Region to store region_mark */
//	DOF *Region = phgDofNew(g,DOF_DEFAULT,1,"Region",DofInterpolation);
//	phgDofSetLambdaFunction(Region, region);
	/**************************************************/

        FLOAT t0_total = phgGetTime(NULL), t0, t1, t1_total;
	int count = 0;

	while(TRUE){
		count ++;
		phgPrintf("STEP = %d\n",count);
		t0 = phgGetTime(NULL);		
		phgPrintf("%d DOF, %d elements, %d submeshes, load imbalance: %lg\n",
				DofGetDataCountGlobal(u), g->nleaf_global, g->nprocs,
				(double)g->lif);
			
		//get the solution of (un-un+1)
		solver = phgSolverCreate(SOLVER_DEFAULT, deltau, NULL);
		build_stiff_rhs(solver, u, Grad_G, Grad_H);
		phgSolverSolve(solver, 1, deltau, NULL);
		phgDofAXPY(-1.0, deltau, &u);	
		phgPrintf("Before Relaxation:\n||u||_2=%.4e, nits = %d, residual = %.4e\n",(double)phgDofNormL2(u), solver->nits, (double)solver->residual);
		phgSolverDestroy(&solver);

		/*relaxation to u*/
		phgDofAXPBY(1-alp, tu, alp, &u);		
		phgPrintf("After Relaxation:\n||u||_2=%.4e\n", (double)phgDofNormL2(u));

		phgDofAXPY(-1.0, u,&tu);
		u2=phgDofNormL2(tu);
		phgPrintf("||delta_u||_2=%12.8e\n",(double)u2);
		if(u2<tol) break;
		else phgDofCopy(u, &tu, NULL, NULL);

		t1 = phgGetTime(NULL);
		mem = phgMemoryUsage(g, &mem_peak);
		phgPrintf("Memory = %.2lf MB, Wall time = %.2f s\n",(double)mem_peak/(1024*1024),(double)(t1-t0));
	}

	t1_total = phgGetTime(NULL);
	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("\nMemory = %.2lf MB, Total Wall Time = %.2f s.\n", (double)mem_peak/(1024*1024), (double)(t1_total - t0_total));

	/* solvation energy calculation*/
	s_E(u, H);//solvation energy

	/* radical distributions (u, c[NION], ch) for spherical*/
	//Sphere(u, 0, num, delta);
	/**************************************************/	

	/*For DNA case, get concentration files*/
//	DNA(u, Region);
	/***************************************************/
#if 1
	/* creating vtk file*/
	INT i, dc = DofGetDataCount(u);
	for (i = 0; i < dc; i++)	u->data[i] *= (kB*T/e);	
	for (i = 0; i < dc; i++)	H->data[i] *= (kB*T/e);
	DOF* U_h;
	U_h = phgDofCopy(u, NULL, NULL, "U_h");
	phgDofAXPY(1.0, H, &U_h);
	phgDofAXPY((kB*T/e), G, &U_h);
	phgPrintf("\nCreating \"%s\".\n", phgExportVTK(g, fn_vtk, u, H,U_h , NULL));
	phgDofFree(&U_h);
#endif

#if 0
        char *delims = ".";
        char str[50];
        strcpy(str, fn_mesh);
        phgPrintf("fn_mesh = %s, str = %s\n", fn_mesh, str);
        char *result = NULL;
        result = strtok(str, delims);
        strcat(result, ".txt");
        FILE *fu;
        fu = fopen(result,"w");
        phgPrintf("The energy is %lf\n", energy);
        fprintf(fu, "%lf\n", (double)energy);
#endif       
	phgDofFree(&G);
	phgDofFree(&H);
	phgDofFree(&Grad_G);
	phgDofFree(&Grad_H);
	phgDofFree(&u);
	phgDofFree(&tu);
	phgDofFree(&deltau);
//	phgDofFree(&Region);
	phgFreeGrid(&g);
	phgFinalize();
    
	return 0;
}

