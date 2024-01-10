/* calculate potential and concentration profiles for a sphere cavity*/
static void Sphere(DOF *u, FLOAT ra, FLOAT rb, FLOAT tr){
	FLOAT phi, theta, r;
	FLOAT tp = M_PI/10, tt = M_PI/5;
	INT cc;
	INT i, j, k;
	GRID *g = u->g;
	COORD x[100];
	FLOAT res1[100], *res;
	FLOAT sam[100][3];
	int n = (int) rb/tr;
	if(!(res = (FLOAT *)phgAlloc(n * sizeof(*res))))
		phgError(1, "Error in memory allocation for sv func.\n");
	r = ra;

	/*sample points on a unit sphere*/
	phi = 0;
	for(i = 0; i < 10; i++){
		theta = 0;
		for(j = 0; j < 10; j++){
			sam[10 * i + j][0] = Cos(phi) * Cos(theta);
			sam[10 * i + j][1] = Cos(phi) * Sin(theta);
			sam[10 * i + j][2] = Sin(phi);
			theta += tt;
		}
		phi += tp;
	}
	/**********************************/

	/*calculate reduced potential values e_c * \beta * \phi*/
	for(i = 0; r < rb; i++){
		cc = 0;
		/*sample points on a sphere of radius r*/
		for(j = 0; j < 10; j++){
			for(k = 0; k < 10; k++){
				x[cc][0] = r * sam[cc][0];
				x[cc][1] = r * sam[cc][1];
				x[cc][2] = r * sam[cc][2];
				cc++;
			}
		}
		phgInterGridDofEval(u, 100, x, res1, -1);

//		if(g->rank == 0){
			res[i] = 0;
			for(k = 0; k < 100; k++) 
			{
				res[i] += res1[k];
			}
			res[i] /= 100;
//		}
		r += tr;
	}
	/******************************************/

	/**************calculate charge*************/
	int start = (int)atom[0].r/tr;
	FLOAT *rho, *ch, rp, rq;
	if(!(rho = (FLOAT *)phgAlloc(n * sizeof(FLOAT))))
		phgError(1, "Error in memory allcoation for rho.\n");
	if(!(ch = (FLOAT *)phgAlloc(n * sizeof(FLOAT))))
		phgError(1, "Error in memory allcoation for ch.\n");
	for(i = 0; i < n; i++)
		rho[i] = 0.0;
	for(i = 0; i < n; i++){
		if(i >= start)
			for(j = 0; j < NION; j++)
				rho[i] += ion[j].q * ion[j].c * Exp(-ion[j].q * res[i]) * Na * 1e-27;
	}
	for(i = 0; i < n; i++){
		if(i <= start)	ch[i] = atom[0].q;
		else{
//			rp = i * delta;
//			rq = rp - delta;
//			ch[i] = ch[i-1] + 4/3 * M_PI * (pow(rp, 3)-pow(rq, 3)) * rho[i-1];
			ch[i] = ch[i-1] + 4 * M_PI * pow(delta, 3) * (i-1) * (i-1) * rho[i-1];
		}
	}
	/***********************************************/

	/*************print results******************/
	phgPrintf("   r          u            ");
	for(j = 0; j < NION; j++){
		phgPrintf("c[%d]          ", j);
	}
	phgPrintf("ch\n");

	for(i = start; i * delta < rb; i++){
		phgPrintf("%6.2f%14.4e", i * delta, res[i]);
		for(j = 0; j < NION; j++){
			phgPrintf("%14.4e", ion[j].c * Exp(-ion[j].q * res[i]));
		}
		phgPrintf("%14.4e\n", ch[i]);
	}

	if(g->rank==0){
		FILE *fu;
		fu = fopen("u","w");
		fprintf(fu,"   r          u            ");
		for(j = 0; j < NION; j++){
			fprintf(fu, "c[%d]          ", j);
		}
		fprintf(fu, "ch\n");
		for(i = start; i * delta < rb; i++){
			fprintf(fu, "%6.2f%14.4e", i * delta, res[i]);
			for(j = 0; j < NION; j++){
				fprintf(fu, "%14.4e", ion[j].c * Exp(-ion[j].q * res[i]));
			}
			fprintf(fu, "%14.4e\n", ch[i]);
		}
		for(j = 0; j < NION; j++){
			fprintf(fu, "c[%d]_bulk:%12.4f\n", j, ion[j].c);
		}
		fclose(fu);
	}
	/***********************************************/
	phgFree(res);
	res = NULL;
	phgFree(rho);
	rho = NULL;
	phgFree(ch);
	ch = NULL;
}

//define a new dof to store region_mark
static void region(DOF *dof, SIMPLEX *e, int bno, const FLOAT lambda[], FLOAT *value){
	*value = e->region_mark;
}

//for obtain the concentration distribution of dna
static void DNA(DOF *u, DOF *reg){
	//define box
	FLOAT xmin = -12, xmax = 12;
	FLOAT ymin = -12, ymax = 12;
	FLOAT zmin = -20, zmax = 20;

	//sample point freqency
	FLOAT r;//x,y are polar coordinates,z is Dicare coordinate
	FLOAT tp = 2*M_PI/10, td = 0.4, tr = 0.1;

	INT n;
	INT i, j, k;
	FLOAT **cv;
	if(!(cv = (FLOAT *)phgAlloc(NION * sizeof(*cv))))
		phgError(1, "Error in memory allcoation for cv.\n");
	for(i = 0; i < NION; i++)
		cv[i] = (FLOAT *)phgAlloc(300*sizeof(FLOAT));

	GRID *g = u->g;
	COORD x[100];
	FLOAT phi[100], res[100];
	FLOAT sam[100][3];

	//initialize sample points
	for(i=0;i<10;i++){
		for(j=0;j<10;j++){
			sam[10*i+j][0] = Cos(j*tp);
			sam[10*i+j][1] = Sin(j*tp);
			sam[10*i+j][2] = -2.0+i*td;
//			sam[10*i+j][2] = 0;
		}
	}
//for(i=0;i<100;i++){phgPrintf("%d %lf  ", i, sam[i][0]*sam[i][0]+sam[i][1]*sam[i][1]+sam[i][2]*sam[i][2]);}

	for(i=0;i*tr<30;i++){
		n = 0;
		r = i*tr;
		for(j=0;j<10;j++){
			for(k=0;k<10;k++){
				x[n][0] = r*sam[n][0];
				x[n][1] = r*sam[n][1];
				x[n][2] = sam[n][2];
//if(i<2) phgPrintf("%lf  %lf  %lf	",x[n][0],x[n][1],x[n][2]);
				n++;
			}
		}
		phgInterGridDofEval(u, 100, x, phi, 0);
		phgInterGridDofEval(reg, 100, x, res, 0);
		if(g->rank == 0){
			for(j = 0; j < NION; j++){
				cv[j][i] = 0.0;
				for(k = 0; k < 100; k++){
					if(res[k] == 1)	continue;
					cv[j][i] +=ion[j].c * Exp(-ion[j].q * phi[k]);
				}
				cv[j][i] /=100;
			}
		}
	}
//for(i=0;i*tr<30;i++){phgPrintf("%lf	",res[i]);}

	/**********output results*******************/
	phgPrintf("  Distance");
	for(i = 0; i < NION; i++)
		phgPrintf("     c[%d]     ", i);
	phgPrintf("\n");
	for(i = 0; i * 0.1 < 30; i++){
		phgPrintf("%8.2f", i * 0.1);
		for(j = 0; j < NION; j++)
			phgPrintf("%14.4e", cv[j][i]);
		phgPrintf("\n");
	}

	if(g->rank==0){
		FILE *fc;
		fc = fopen("c", "w");
		fprintf(fc, "  Distance");
		for(j = 0; j < NION; j++){
			fprintf(fc, "     c[%d]     ", j);
		}
		fprintf(fc, "\n");
		for(i = 0; i * 0.1 < 30; i++){
			fprintf(fc, "%6.2f", i * 0.1);
			for(j = 0; j < NION; j++){
				fprintf(fc, "%14.4e", cv[j][i]);
			}
			fprintf(fc, "\n");
		}
		for(j = 0; j < NION; j++){
			fprintf(fc, "c[%d]_bulk:%12.4f\n", j, ion[j].c);
		}
		fclose(fc);
	}
	/********************************************/
	phgFree(cv);
	cv = NULL;
}
