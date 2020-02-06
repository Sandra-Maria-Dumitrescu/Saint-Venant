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

		 Bodies of functions and procedures for the initial conditions

   ============================================================================= */
#include "initialConditions.h"
void damBreakWetBed(float h[nbCell+2],float hu[nbCell+2],float Z[nbCell+2],float dx)
{
	int i;
	float midDomain;
	float x;
	midDomain = (bSup+bInf)/2.0;
	
	for(i=0;i<=nbCell+1;i++)
	{
		x = (bInf + dx/2.0) + (i-1)*dx;
		if(x < midDomain)
		{
			h[i] = 2.0;
		}
		else
		{
			h[i] = 1.0;
		}
		hu[i] = 0.0;
	}
}

void damBreakDryBed(float h[nbCell+2],float hu[nbCell+2],float Z[nbCell+2],float dx)
{
	int i;
	float x;
	float midDomain;
	
	midDomain = (bSup+bInf)/2.0;
	
	for(i=0;i<=nbCell+1;i++)
	{
		x = (bInf + dx/2.0) + (i-1)*dx;
		if(x < midDomain)
		{
			h[i] = 1.0;
		}
		else
		{
			h[i] = 0.0;
		}
		hu[i] = 0.0;
	}
}

void lakeAtRest(float h[nbCell+2],float hu[nbCell+2],float Z[nbCell+2])
{
	int i;
	
	for(i=0;i<=nbCell+1;i++)
	{
		h[i] = 2.0 - Z[i];
		hu[i] = 0.0;
	}
}
