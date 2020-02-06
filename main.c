/* =============================================================================

	A finite volume code to solve the Saint-Venant system in an open channel.

   =============================================================================

	Copyright (C) 2013 Mathieu Besson

	This program is free software: you can redistribute
	it and/or modify it under the terms of the
	GNU General Public License as published by the Free Software Foundation, either 
	version 3 of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful, 
	but WITHOUT ANY WARRANTY; without even the implied warranty 
	of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
	See the GNU General Public License for more details.

	You should have received a copy of the GNU General Public License 
	along with this program. If not, see <http://www.gnu.org/licenses/>.*/

/* This program is the main program to solve the Saint-Venant system.
 * Before choosing your options in this file, make sure that you have choose the constants you need.
 * You must check : 
 * 						- bInf which is the lower bound of the domain.
 * 						- bSup which is the upper bound of the domain.
 * 						- nbCell which is the number of cells for the discretization into finite volumes.
 * 						- freq which is the frequence to save data; For example, each 5 iteration, you save your results in a file ".txt".
 * 						- T which is the duration (in seconds) of the experience.
 * 
 * Of course, you can also change : 
 * 
 * 						- epsilon which is a criterion usefull for testing if a quantity is equal to zero.
 * 						- N which is a parameter to to a quadrature to compute integrals.
 * 
 * but this is not necessary.
 * Then, when you have the constants you want, you have to check your options in the main program : 
 * 
 * 						- First, select the topography you want.
 * 						- Then, select the initial condition.
 * 						- If you have several solver, choose the one you want to use.
 * 						- Select the good function to compute the time step.
 * 						- Select the boundary condition.
 * 						- And run the program !!
 * 
 * In some case (dam break on wet/dry bed etc) you can compare your results with the exact solution.
 * If you want to do this :
 * 
 * 						- First, initialize the variable hExact (exact water height) ant huexact ( exact discharge).
 * 						- Then, choose the right exact solution (Ritter solution for a dam break on dry bed, Stoker solution for a dam break on wet bed etc).
 * 						- Don't forget to save the exact solution.
 */						

#include "kineticSolver.h"
#include "solverRusanov.h"
#include "initialConditions.h"
#include "boundaryCondition.h"
#include "usefull.h"
#include "constant.h"
#include "topography.h"
#include "timeSetting.h"
#include "testCase.h"
#include <time.h>
#include <stdio.h>

int main(void)
{
	float dx; /* space step */
	float dt; /* time step */
	float tOld,tNew;p
	float rat; /* ratio dt/dx */
	float hOld[nbCell+2];
	float hNew[nbCell+2];
	float hExact[nbCell+2];
	float huOld[nbCell+2];
	float huNew[nbCell+2];
	float huExact[nbCell+2];
	float Z[nbCell+2];
	int i,compt,tmp;
	float norm;
	FILE *fichier2;
	FILE *fichier3;
	float temps;
    clock_t t1, t2;
	
	/* We compute the space step. */
	/*-----------------------------------------------------------------*/
		dx = (bSup-bInf)/(nbCell+1); /* we compute the space step */
	/*-----------------------------------------------------------------*/
	
	/* Select the topography you want */
	/*-----------------------------------------------------------------*/
		bumpBed(Z,dx);
		/*stepBed(Z,dx);*/
		/*noBed(Z);*/
		/*constantBed(Z);*/
		
		/* We save the topograhy in the file topo.txt*/
		fichier3 = fopen("topo.txt","w");
		for(i=0;i<=nbCell+1;i++)
		{
			fprintf(fichier3,"%f ",Z[i]);
		}
		fclose(fichier3);
	/*-----------------------------------------------------------------*/
	
	/* If you want to compute exact solution, first initialise hExact and huExact */
	/*-----------------------------------------------------------------*/
		for(i=0;i<=nbCell+1;i++) 
		{
			hExact[i] = 2.0 - Z[i];
			huExact[i] = 0.0;
		}
	/*-----------------------------------------------------------------*/
	
	/* Choose the initial condition you want */
	/*-----------------------------------------------------------------*/
		/*damBreakWetBed(hOld,huOld,Z,dx);*/
		/*damBreakDryBed(hOld,huOld,Z,dx);*/
		lakeAtRest(hOld,huOld,Z);
	/*-----------------------------------------------------------------*/
	
	dt = 0.0;
	tOld = 0.0;
	tNew = 0.0;
	tmp = 1;
	compt = 0;
	norm = 0;
	
	save(hOld,huOld,tmp); /* saving initial data */
	
	fichier2 = fopen("time.txt","w");
	fprintf(fichier2,"%f ",tNew);
	t1 = clock();
	while(tNew<T) /* loop on time*/
	{
		compt = compt + 1;
		
		/* We compute the time step thanks to te CFL condition. */
		/*-----------------------------------------------------------------*/
			dt = timeStepKinetic(dx,huOld,hOld); /* Choose this one if you want to use the kinetic solver. */
			/*dt = timeStepRusanov(dx,huOld,hOld);*/ /* Choose this one if you want to use the Rusanov solver. */
		/*-----------------------------------------------------------------*/
		
		tNew = tOld + dt;
		if(tNew > T)
		{
			dt = T - tOld;
			tNew = T;
		}
		printf("T = %f		dt = %f\n",tNew,dt); /* a print to have a look on the program. */
		
		rat = dt/dx;
		/* Choose the solver you want */
		/*-----------------------------------------------------------------*/
			updateKinetic(hNew,huNew,hOld,huOld,Z,rat);
			/*updateRusanov(hNew,huNew,hOld,huOld,Z,rat);*/
		/*-----------------------------------------------------------------*/

		/* Choose the boundary conditions you want */
		/*-----------------------------------------------------------------*/
			/*updateBC(hNew,huNew);*/ /* use this for free boundary conditions. */
			/*updateBCPeriodic(hNew,huNew,hOld,huOld);*/ /* use this to periodic boundary conditions. */
			updateBCLakeAtRest(hNew,huNew); /* we need to use this method to compute a lake at rest. */
		/*-----------------------------------------------------------------*/
		
		
		/* If you want to compare your results with the exact solution, choose the right one. */
		/*-----------------------------------------------------------------*/
			/*Ritter(hExact,huExact,dx,tNew,1.0);*/ /* Use this function for a dam break on dry bed */
		/*-----------------------------------------------------------------*/
		
		/* We save our results in file ".txt" */
		/*-----------------------------------------------------------------*/
			if((compt%freq==0)||(tNew==T))
			{
				tmp = tmp + 1;
				save(hNew,huNew,tmp);
				saveExact(hExact,huExact,tmp);
				norm = getNorme(hNew,hExact,dx);
				printf("|h-hex| = %f\n",norm);
				norm = getNorme(huNew,huExact,dx);
				printf("|hu-huex| = %f\n",norm);
				fprintf(fichier2,"%f ",tNew);
			}
		/*-----------------------------------------------------------------*/
		
		for(i=0;i<=nbCell+1;i++)
		{
			hOld[i] = hNew[i];
			huOld[i] = huNew[i];
		}
		tOld = tNew;
	}
	fclose(fichier2);
	t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    printf("Computation time = %f\n", temps);
    printf("number of iterations = %d\n",compt);
	return 0;
}
