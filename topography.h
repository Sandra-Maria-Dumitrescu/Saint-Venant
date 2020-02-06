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

		 Prototypes of functions and procedures for the topography

   ============================================================================= */
#ifndef topography_h
	#define topography_h

	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "constant.h"
	#include <math.h>

	void noBed(float Z[nbCell+2]); /* Here there is no topography */
	void bumpBed(float Z[nbCell+2],float dx); /* Here the function Z is given like Z(x) = (0.2-0.005x²) x in [-2;2]*/
	void bumpBedConstantPiecwise(float Z[nbCell+2],float dx); /* This is the piecwise constant representaion of the bump function. */
	void stepBed(float Z[nbCell+2],float dx); /* Here the topography is discontinous with  a step at the mid-domain. */
	void constantBed(float Z[nbCell+2]); /* Here the topography is constant. */
	
#endif
