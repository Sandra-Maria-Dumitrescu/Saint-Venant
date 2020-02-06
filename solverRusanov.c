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

/* =============================================================================

		   Bodies of functions and procedures for the Rusanov solver

   ============================================================================= */
#include "solverRusanov.h"
void updateRusanov(float hNew[nbCell+2],float huNew[nbCell+2],float hOld[nbCell+2],float huOld[nbCell+2],float Z[nbCell+2],float rat)
{
	int i;
	float fh[nbCell+2];
	float fhu[nbCell+2];
	float c1[nbCell+2]; /* here c1 is define as c1 = max(|u[i+1]|+sqrt(gh[i+1]),|u[i]|+sqrt(gh[i])) */
	float c2[nbCell+2]; /* here c2 is define as c2 = max(|u[i-1]|+sqrt(gh[i+1]),|u[i-1]|+sqrt(gh[i])) */
	
	for(i=0;i<=nbCell+1;i++) /* we initialize the vectors fh and fhu before computing. */
	{
		fh[i] = 0.0;
		fhu[i] = 0.0;
		c1[i] = 0.0;
		c2[i] = 0.0;
	}

	for(i=1;i<=nbCell;i++) /* Then we compute c1 and c2 for each cell in the domain (so not for cells number 0 and n+1). */
	{
		c1[i] = max(vAbs(huOld[i+1]/hOld[i+1]) + sqrt(g*hOld[i+1]),vAbs(huOld[i]/hOld[i]) + sqrt(g*hOld[i]));
		c2[i] = max(vAbs(huOld[i-1]/hOld[i-1]) + sqrt(g*hOld[i-1]),vAbs(huOld[i]/hOld[i]) + sqrt(g*hOld[i]));
	}
	
	/* Before updating hNew and huNew we compute the differents flux. */
	flux_hRusanov(hOld,huOld,fh,c1,c2);
	flux_huRusanov(hOld,huOld,fhu,c1,c2);
	
	for(i=1;i<=nbCell;i++)
	{
		hNew[i] = hOld[i] - rat*fh[i];
		huNew[i] = huOld[i] - rat*fhu[i];
	}
}

void flux_hRusanov(float h[nbCell+2],float hu[nbCell+2],float fh[nbCell+2],float c1[nbCell+2],float c2[nbCell+2])
{
	int i;
	
	/* In this case we don't mind about the flux at cell number 0 and n+1. That is why we put the value 0.0. */
	for(i=1;i<=nbCell;i++)
	{
		fh[i] = 0.5*(hu[i+1] - hu[i-1]) - 0.5*c1[i]*(h[i+1] - h[i]) + 0.5*c2[i]*(h[i] - h[i-1]);
	}
}

void flux_huRusanov(float h[nbCell+2],float hu[nbCell+2],float fhu[nbCell+2],float c1[nbCell+2],float c2[nbCell+2])
{
	int i;
	
	/* In this case we don't mind about the flux at cell number 0 and n+1. That is why we put the value 0.0. */
	for(i=1;i<=nbCell;i++)
	{
		fhu[i] = 0.5*(((hu[i+1]*hu[i+1])/h[i+1]+0.5*g*h[i+1]*h[i+1])-((hu[i-1]*hu[i-1])/h[i-1]+0.5*g*h[i-1]*h[i-1]));
		fhu[i] = fhu[i] - 0.5*c1[i]*(hu[i+1] - hu[i]) + 0.5*c2[i]*(hu[i] - hu[i-1]);
	}	
}
