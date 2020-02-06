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

		Prototypes of functions and procedures for the time step computation

   ============================================================================= */
#ifndef timeSetting_h
	#define timeSetting_h

	#include <stdio.h>
	#include "constant.h"
	#include "usefull.h"
	#include <math.h>
	
	float timeStepRusanov(float dx,float huOld[nbCell+2],float hOld[nbCell+2]); /* to compute the right time for the Rusanov solver. */
	float timeStepKinetic(float dx,float huOld[nbCell+2],float hOld[nbCell+2]); /* computes the right time step for each time iteration for the kinetic flux */
#endif

