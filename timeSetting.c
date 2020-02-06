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

		Bodies of functions and procedures for the time step computation

   ============================================================================= */
#include "timeSetting.h"

float timeStepRusanov(float dx,float huOld[nbCell+2],float hOld[nbCell+2])	
{
	float dt = 0.0;
	float EGPlus[nbCell+2]; /* the eigen value u + c */
	float EGMinus[nbCell+2]; /* the eigen value u - c */
	float maxEGPlus = 0.0; /* max eigen value of the vector EGPlus */
	float maxEGMinus = 0.0; /* max eigen value of the vector EGMinus */
	float maxEG = 0.0; /* maximum of maxEGPlus and maxEGMinus */
	int i;
	
	for(i=0;i<=nbCell+1;i++)
	{
		EGPlus[i] = 0.0;
		EGMinus[i] = 0.0;
	}
	
	getEigenValues(huOld,hOld,EGPlus,EGMinus); /* first we compute eigen values */
	
	/* we compute the maximum eigen value for each vector */
	maxEGPlus = maxVal(EGPlus);
	maxEGMinus = maxVal(EGMinus);
	
	/* Then we take the maximum of the two different eigen values */
	maxEG = max(maxEGPlus,maxEGMinus);
	
	if(maxEG <= 0.001)
	{
		dt = dx*0.1;
	}
	else
	{
		dt = 0.5*(dx/maxEG);
	}
	return dt;
}

float timeStepKinetic(float dx,float huOld[nbCell+2],float hOld[nbCell+2])
{
	float dt = 0.0;
	float EG[nbCell+2]; /* the eigen value |u| + sqrt(2gh) */
	float maxEG = 0.0; /* maximum value of the vector EG */
	int i;
	
	for(i=0;i<=nbCell+1;i++)
	{
		if(hOld[i]<=epsilon)
		{
			EG[i] = 0.0;
		}
		else
		{
			EG[i] = vAbs(huOld[i]/hOld[i]) + sqrt(2.0*g*hOld[i]);
		}
	}
	
	/* we compute the maximum eigen value of EG */
	maxEG = maxVal(EG);
	
	if(maxEG <= epsilon)
	{
		dt = dx*0.1;
	}
	else
	{
		dt = dx/maxEG;
	}
	return dt;
}
