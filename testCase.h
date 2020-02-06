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

			 Prototypes of functions and procedures for test cases

   ============================================================================= */
#ifndef testCase_h
	#define testCase_h

	#include <stdio.h>
	#include "constant.h"
	#include "usefull.h"
	#include <math.h>
	
	void Ritter(float h[nbCell+2],float hu[nbCell+2],float dx, float t,float hg); /* This procedure we compute the Ritter solution for a dam break on dry bed. */
#endif
