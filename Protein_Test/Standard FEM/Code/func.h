/* func_G, func_gradG, func_pb and Harmonic_stiff_rhs */

#ifndef __FUNC_H__
#define __FUNC_H__

/*Read atom from *.pqr*/
static void Read_pqr(const char *fn){
	FILE *f_in;
	char s[100];
	FLOAT charge=0.0;

	f_in = fopen(fn, "r");

	if (f_in == NULL){
		phgError(1, "Can not open file %s.\n", fn);
	}

	NATOMS = 0;
//	if(!(atom = (struct atoms *)phgAlloc(MNATOMS * sizeof(*atom))))
	if(!(atom = phgAlloc(MNATOMS * sizeof(*atom))))
		phgError(1, "Error in memory allocation for atoms!\n");

	double x = 0.0, y = 0.0, z = 0.0, zi = 0.0, r = 0.0;
	int number;
    while(fscanf(f_in, "%s", s) != EOF) {
		if(strcmp(s, "ATOM") != 0) continue;
		fscanf(f_in, "%s%s%s%s%lf%lf%lf%lf%lf",s,s,s,s,&x,&y,&z,&zi,&r);
//#if 1
//			if(r < 1e-3){
//				phgPrintf("warning %d atom: r = %lf\n", number, r);
//				r = 0.1;
//			}
//#else
//			if(r < 0.5){
//				r = 0.6;
//			}
//#endif
		atom[NATOMS].x = (FLOAT)x;
		atom[NATOMS].y = (FLOAT)y;
		atom[NATOMS].z = (FLOAT)z;
		atom[NATOMS].q = (FLOAT)zi;
		atom[NATOMS].r = (FLOAT)r;
		charge +=atom[NATOMS].q;
		NATOMS++;
	}
	phgPrintf("\n");
	phgPrintf("NATOMS                             %d\n",NATOMS);
	phgPrintf("net charge                         %lf e\n",charge);
	fclose(f_in);
}

/*set boundary map*/
static int bc_map(int bctype){
	// if(bctype == BC_MAP[0]) 
	// 	return DIRICHLET;
	// else if(bctype == BC_MAP[1])  
	// 	return BDRY_USER1;
	// else{ 
	// 	return -1;
	//  	phgError(1, "UNDEFINED BOUNDARY TYPR!\n");
	// }
	if(bctype == 9)
		return BDRY_USER1;
	else
		return DIRICHLET;
}

/*to determine the \epsilon */
static void epsilon(SIMPLEX * e, FLOAT * eps, FLOAT *lam){
	if (e->region_mark == 1){
		*lam = 0;
		*eps = e_m;
	}
    else if (e->region_mark == 3){
		*eps = e_e;
		*lam = 0;
	}
	else if (e->region_mark == 2){
		*eps = e_s;
		*lam = 1.0;
	}
	else	phgError(-1, "%s%d: UNKOWN REGION MARK!\n", __FILE__, __LINE__);
}

static void epsilon1(SIMPLEX *e, FLOAT *eps, FLOAT *lam, FLOAT z1, FLOAT z2){
	FLOAT z[4], Z;
	GRID *g;
	INT i;
	for(i=0; i<4; i++)
		z[i]=g->verts[e->verts[i]][2];
	Z=(z[0]+z[1]+z[2]+z[3])/4.0;
	if (e->region_mark == 1)
	{	*eps = e_m;  *lam = 0;}
	else if (e->region_mark == 2 && (Z > z2 || Z < z1))
         {   	*eps = e_s;  *lam = 1;}
	else if (e->region_mark == 2 && z1 < Z &&Z < z2)        
         {      *eps = e_e;  *lam = 0;}
	     else	phgError(-1, "%s%d: UNKOWN REGION MARK!\n", __FILE__, __LINE__);
}

/*define function G*/
static void func_G(FLOAT x, FLOAT y, FLOAT z, FLOAT *value){
	FLOAT dist;
	int i;
	*value = 0;
	for( i = 0; i < NATOMS; i++){
		dist = r2(x, y, z, atom[i].x, atom[i].y, atom[i].z);
        if(dist<1e-10)  dist = 1e-6;
		*value += (atom[i].q/dist);
	}
	*value *= K1/(4 * M_PI * e_m);
}

/*define function Grad_G*/
static void func_gradG(FLOAT x, FLOAT y, FLOAT z, FLOAT *value){
	FLOAT ret[3], r;
	int i;
	ret[0] = ret[1] = ret[2] = 0.0;
	for (i = 0; i < NATOMS; i++) {
		r = r2(x, y, z, atom[i].x, atom[i].y, atom[i].z);
        if(r<1e-10) r = 1e-6;
		ret[0] += -1.0 * atom[i].q * (x - atom[i].x) / pow(r,3);
		ret[1] += -1.0 * atom[i].q * (y - atom[i].y) / pow(r,3);
		ret[2] += -1.0 * atom[i].q * (z - atom[i].z) / pow(r,3);
        }
	*value++ = ret[0] * K1 /(4*M_PI* e_m);
	*value++ = ret[1] * K1 /(4*M_PI* e_m);       
	*value++ = ret[2] * K1 /(4*M_PI* e_m);  
}

static void func_uc(FLOAT x, FLOAT y, FLOAT z, FLOAT *value){
        FLOAT dist;
	int i;
	*value = 0;
	for( i = 0; i < NATOMS; i++){
		dist = r2(x, y, z, atom[i].x, atom[i].y, atom[i].z);
		if(Fabs(dist) < 1e-6){
			if(dist > 0){
				dist = 1e-6;
			}
			else{
				dist = -1e-6;
			}
		}
		*value += (atom[i].q * Exp(-Kappa * dist * 1.0e-10) / dist);
	}
	*value *= e * e /(4 * M_PI * e_s * e0 * kB * T) * 1.0e10;
}
/*set boundary value of u*/
static void func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value){
	*value = 0.0;
}

/*build stiffness matrix and RHS for H*/
static void Harmonic_stiff_rhs(SOLVER *solver,DOF *H){
	int i, j;
	DOF *HH;
	HH = phgDofCopy(H,NULL,NULL,NULL);
	GRID *g = H->g;
	SIMPLEX *e;

	ForAllElements(g,e){
		int N = DofGetNBas(HH,e);
		FLOAT A[N][N],rhs[N],buffer[N];
		INT I[N];
        
		bzero(rhs, N * sizeof(*rhs));
      
		for(i=0;i<N;i++)  I[i] = phgSolverMapE2L(solver, 0, e, i);
    
		if(e->region_mark == 1){
			for (i = 0; i < N; i++) {
				for (j = 0; j <= i; j++)
					A[j][i] = A[i][j] = phgQuadGradBasDotGradBas(e, H, i, H, j, QUAD_DEFAULT);
			}
		}
		else if(e->region_mark == 2){
			bzero(A, N * N * sizeof(**A));
			for (i = 0; i < N; i++) {
				if (phgDofGetElementBoundaryType(H, e, i) & (BDRY_USER1))
					A[i][i] = 0.0;
				else
					A[i][i] = 1.0;
			}
		     }
        
		phgDofSetDirichletBoundaryMask(HH, BDRY_USER1);
		for(i = 0; i < N; i++){
			BOOLEAN is_bdry = FALSE;
			if(e->region_mark == 1){
				is_bdry = phgDofDirichletBC(HH, e, i, func_G, buffer, rhs + i, DOF_PROJ_NONE);
				if(is_bdry)	*(rhs+i) *= -1.0;
			}
			if (is_bdry)
				phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer);
			else
				phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]);    
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}
	phgDofFree(&HH);
}

/*to calculate solvation energy*/
static void s_E(DOF *phi_r, DOF *phi_h){
	int i;
	COORD pos[NATOMS];
	FLOAT u_r[NATOMS], u_h[NATOMS], sum = 0.0;
 	DOF *u1_r, *u1_h;
	u1_r = phgDofCopy(phi_r, NULL, NULL, NULL);
	u1_h = phgDofCopy(phi_h, NULL, NULL, NULL);

	for(i = 0; i < NATOMS; i++){
		pos[i][0] = atom[i].x;
		pos[i][1] = atom[i].y;
		pos[i][2] = atom[i].z;
	}
	phgInterGridDofEval(u1_r, NATOMS, pos, u_r, -1);
	phgInterGridDofEval(u1_h, NATOMS, pos, u_h, -1);
	for(i = 0; i < NATOMS; i++){
		sum +=(u_r[i] + u_h[i]) * atom[i].q;
		//phgPrintf("u_r[%d] = %f\n", i, u_r[i]);
		//phgPrintf("u_h[%d] = %f\n", i, u_h[i]);
	} 
	sum *= 0.5;
    energy = sum * kBT;
	phgPrintf("\nSolvation energy = %lf\n", (double)(sum * kBT));
	phgDofFree(&u1_r);
	phgDofFree(&u1_h);
}
#endif
