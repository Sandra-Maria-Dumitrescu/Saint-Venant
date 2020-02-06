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

			Bodies of functions and procedures for boundary conditions

   ============================================================================= */
#include "boundaryCondition.h"
void updateBC(float hNew[nbCell+2],float huNew[nbCell+2])
{
	hNew[0] = hNew[1];
	huNew[0] = huNew[1];
	hNew[nbCell+1] = hNew[nbCell];
	huNew[nbCell+1] = huNew[nbCell];
}

void updateBCPeriodic(float hNew[nbCell+2],float huNew[nbCell+2],float hOld[nbCell+2],float huOld[nbCell+2])
{
	hNew[0] = hOld[nbCell];
	huNew[0] = huOld[nbCell];
	hNew[nbCell+1] = hOld[1];
	huNew[nbCell+1] = huOld[1];
}

void updateBCLakeAtRest(float hNew[nbCell+2],float huNew[nbCell+2])
{
	huNew[0] = 0.0;
	huNew[nbCell+1] = 0.0;
	
	hNew[0] = 2.0;
	hNew[nbCell+1] = 2.0 ;
}
