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

		   Prototypes of functions and procedures for the Rusanov solver

   ============================================================================= */
#ifndef solverRusanov_h
	#define solverRusanov_h

	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "usefull.h"
	#include "constant.h"
	#include <math.h>

	void updateRusanov(float hNew[nbCell+2],float huNew[nbCell+2],float hOld[nbCell+2],float huOld[nbCell+2],float Z[nbCell+2],float rat); /* This function updates the water height (hNew) and the discharge (huNew) at time n+1 with respect to th water height (hOld) and the discharge (huOld) at time n. The argument rat is the ration between dt and dx. */
	void flux_hRusanov(float h[nbCell+2],float hu[nbCell+2],float fh[nbCell+2],float c1[nbCell+2],float c2[nbCell+2]); /* This function computes the flux for the equation of mass conservation. In this case the flux is the Rusanov flux  */
	void flux_huRusanov(float h[nbCell+2],float hu[nbCell+2],float fhu[nbCell+2],float c1[nbCell+2],float c2[nbCell+2]); /* This function computes the flux for the equation of momentum conservation. In this case the flux is the Rusanov flux. */
	
#endif
