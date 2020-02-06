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

			Prototypes of functions and procedures for boundary conditions

   ============================================================================= */
#ifndef boundaryCondition_h
	#define boundaryCondition_h

	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "constant.h"
	#include <math.h>

	void updateBC(float hNew[nbCell+2],float huNew[nbCell+2]); /* This function updates the boundary conditions (free boundary conditions). */
	void updateBCPeriodic(float hNew[nbCell+2],float huNew[nbCell+2],float hOld[nbCell+2],float huOld[nbCell+2]); /* This function updates the boundary conditions (periodic conditions). */
	void updateBCLakeAtRest(float hNew[nbCell+2],float huNew[nbCell+2]); /* we need to use this method to compute a steady state. */
	
#endif
