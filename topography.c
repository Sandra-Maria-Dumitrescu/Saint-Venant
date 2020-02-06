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

		 Bodies of functions and procedures for the topography

   ============================================================================= */
#include "topography.h"

void bumpBed(float Z[nbCell+2],float dx)
{
	int i;
	float x;
	
	for(i=0;i<=nbCell+1;i++)
	{
		x = (bInf + dx/2.0) + (i-1)*dx;
		if((x>=-2.0)&&(x<=2.0))
		{
			Z[i] = 0.2 - 0.05*x*x;
		}
		else
		{
			Z[i] = 0.0;
		}
	}
}

void bumpBedConstantPiecwise(float Z[nbCell+2],float dx)
{
	int i;
	float xPlus,xMinus;
	
	for(i=0;i<=nbCell;i++)
	{
		xMinus = bInf + i*dx;
		xPlus = bInf + (i+1)*dx;
		
		if((xPlus<=2.0)&&(xMinus>=-2.0))
		{
			Z[i] = 0.2 - (0.05/3.0)*(xPlus*xPlus + xPlus*xMinus + xMinus*xMinus);
		}
		else
		{
			if(((xPlus<=2.0)&&(xPlus>=-2.0))&&(xMinus<-2.0))
			{
				Z[i] = (1.0/dx)*(0.2*xPlus - (0.05/3.0)*xPlus*xPlus*xPlus + 0.8/3.0);
			}
			else
			{
				if(((xMinus>=-2.0)&&(xMinus<=2.0))&&(xPlus>2.0))
				{
					Z[i] = (1.0/dx)*((0.8/3.0) - 0.2*xMinus + (0.05/3.0)*xMinus*xMinus*xMinus);
				}
				else
				{
					Z[i] = 0.0;
				}
			}
		}
	}
	Z[nbCell+1] = 0.0;
}

void stepBed(float Z[nbCell+2],float dx)
{
	int i;
	float xCentre,xMinus,xPlus;
	float midDomain;
	
	midDomain = (bSup+bInf)/2.0;
	
	for(i=0;i<=nbCell;i++)
	{
		xCentre = (bInf + dx/2.0) + (i-1)*dx;
		xMinus = bInf + i*dx;
		xPlus = bInf + (i+1)*dx;
		
		if((xMinus<midDomain)&&(xPlus>midDomain))
		{
			Z[i] = (0.3/dx)*(xPlus - midDomain);
		}
		else
		{
			if(xCentre>=midDomain)
			{
				Z[i] = 0.3;
			}
			else
			{
				Z[i] = 0.0;
			}
		}	
	}
	Z[nbCell+1] = Z[nbCell];
}

void noBed(float Z[nbCell+2])
{
	int i;
	
	for(i=0;i<=nbCell+1;i++)
	{
		Z[i] = 0.0;
	}
}

void constantBed(float Z[nbCell+2])
{
	int i;
	
	for(i=0;i<=nbCell+1;i++)
	{
		Z[i] = 0.5;
	}
}


