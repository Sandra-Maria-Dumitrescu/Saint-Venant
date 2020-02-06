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

		 Prototypes of functions and procedures for the initial conditions

   ============================================================================= */
#ifndef initialConditions_h
	#define initialConditions_h

	#include <stdio.h>
	#include "constant.h"
	
	void damBreakWetBed(float h[nbCell+2],float hu[nbCell+2],float Z[nbCell+2],float dx); /* with this function you can modelise a dam break on a wet bed. */
	void damBreakDryBed(float h[nbCell+2],float hu[nbCell+2],float Z[nbCell+2],float dx); /* with this function you can modelise a dam break on a dry bed. */
	void lakeAtRest(float h[nbCell+2],float hu[nbCell+2],float Z[nbCell+2]); /* to compute a lake at rest. */
#endif
