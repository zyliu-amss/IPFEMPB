#ifndef __PBE_H__
#define __PBE_H__

static void quad_left(QUAD *quad, SIMPLEX *e, DOF *u, const FLOAT *DATA1, int i, int j, FLOAT h, FLOAT *value){
	if(e->region_mark==1) *value=0;
	else{
		int I, J, K;
		GRID *g = u->g;
		const FLOAT *DATA2, *DATA3;
		FLOAT sum = 0.0, res;
		h=phgGeomGetVolume(g, e);
		DATA1=phgQuadGetDofValues(e, u, quad);
		DATA2 = phgQuadGetBasisValues(e, u, i, quad);
		DATA3 = phgQuadGetBasisValues(e, u, j, quad);
		I = quad->npoints;
		for(J = 0; J < I; J++){
			res = 0.0;
			for(K = 0; K < NION; K++){
				res += ion[K].c * ion[K].q * ion[K].q * Exp(-1.0 * ion[K].q * DATA1[J]); //concentration is in unit of M
			}
			sum += DATA2[J] * DATA3[J] * res * quad->weights[J];//e^2*\beta/\epsilon_0*M=4.2434194 A_-2
		}
		*value = K2 * sum * h* e0;
	}
}

static void quad_right(QUAD *quad, SIMPLEX *e, DOF *u, const FLOAT *DATA1, int i, FLOAT h, FLOAT *value){
	if(e->region_mark==1) *value=0;
        else{
		int I, J, K;
		GRID *g=u->g;
		const FLOAT *DATA2;
		FLOAT sum = 0.0, res;
             	h=phgGeomGetVolume(g, e);
		DATA1=phgQuadGetDofValues(e, u, quad);
		DATA2 = phgQuadGetBasisValues(e, u, i, quad);
		I = quad->npoints;
		for(J = 0;J < I;J++){
			res = 0.0;
			for(K = 0; K < NION; K++){
				res += ion[K].c * ion[K].q * Exp(-1.0 * ion[K].q * DATA1[J]); 
			}
			sum -= res * DATA2[J] * quad->weights[J];
		}
		*value= K2 * sum * h * e0;
	}
}

/*solution to the PBE*/
static void build_stiff_rhs(SOLVER *solver, DOF *u, DOF *Grad_G, DOF *Grad_H){
        int i, j;
	GRID *g = u->g;
	SIMPLEX *e;
	QUAD *quad;
        quad = phgQuadGetQuad3D(3);
        DOF *grad_u;
        grad_u = phgDofGradient(u,NULL,NULL,NULL);
	
	ForAllElements(g, e){
		INT N = DofGetNBas(u, e);
		FLOAT A[N][N], rhs[N], buffer[N];
		FLOAT eps = 0.0,lam = 1.0, value;
		INT I[N];	
			
		memset(A, 0, N*N*sizeof(**A));
		memset(rhs, 0, N*sizeof(*rhs));
		
	    epsilon(e, &eps, &lam);
		//epsilon1(e, &eps, &lam,z1, z2);
		const FLOAT *DATA1;
		FLOAT h;
		h = phgGeomGetVolume(g, e);
		DATA1 = phgQuadGetDofValues(e, u, quad);

		for(i = 0; i < N; i++){
			I[i] = phgSolverMapE2L(solver, 0, e, i);
			for(j = 0; j <= i; j++){
                        	quad_left(quad, e, u, DATA1, i, j, h, &value);
				A[i][j] = A[j][i] =eps * phgQuadGradBasDotGradBas(e, u, j, u, i, 3) + lam * value; 
			}
		}
		
		for (i = 0; i < N; i++){
			if (phgDofDirichletBC(u, e, i, func_u, buffer, rhs + i, DOF_PROJ_NONE))
				phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer);
			else {
				quad_right(quad, e, u, DATA1, i, h, &value);
				rhs[i] = eps * phgQuadDofDotGradBas(e, grad_u, u, i, QUAD_DEFAULT) + lam* value;
				if(e->region_mark == 1){
					for(j = 0; j < NFace; j++){
						if(e->bound_type[j] & (BDRY_USER1 )){
                                			rhs[i] += e_m * phgQuadFaceDofDotBas(e, j, Grad_G, DOF_PROJ_DOT, u, i, u->type->order);
                                			rhs[i] += e_m * phgQuadFaceDofDotBas(e, j, Grad_H, DOF_PROJ_DOT, u, i, u->type->order);
                             			}
					}
				}	
				phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]);
			}
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}
	phgDofFree(&grad_u);
}
#endif
