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

		 Prototypes of several functions and procedures needed for the code

   ============================================================================= */
#ifndef init_h
	#define init_h

	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include <math.h>
	#include "constant.h"

	float maxVal(float u[nbCell+2]); /* this function takes a vector and return the maximum value of this vector. */
	float max(float a,float b); /* this function takes two values and return the greatest. */
	float min(float a,float b); /* this function takes two values and return min value. */
	float vAbs(float x); /* returns the absolute value of a given parameter. */
	void save(float tab1[nbCell+2],float tab2[nbCell+2],int t); /* save the computation of a vector tab at iteration t. */
	void saveExact(float tab1[nbCell+2],float tab2[nbCell+2],int t); /* We will use this procedure to save exact values.*/
	void getEigenValues(float hu[nbCell+2],float h[nbCell+2],float EGPlus[nbCell+2],float EGMoins[nbCell+2]); /* computes all the eigen values (EGPlus = u + c and EGMinus = u-c) for the two vectore hu (dicharge) and h( water height) */
	float power3Over2(float a); /* computes a^(3/2)*/
	float getNorme(float h[nbCell+2],float hex[nbCell+2],float dx); /* computes the intinite norme*/
	
#endif
