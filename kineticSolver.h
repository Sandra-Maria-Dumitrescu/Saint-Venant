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

		   Prototypes of functions and procedures for the kinetic solver

   ============================================================================= */
#ifndef kineticSolver_h
	#define kineticSolver_h

	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "usefull.h"
	#include "constant.h"
	#include <math.h>

	void updateKinetic(float hNew[nbCell+2],float huNew[nbCell+2],float hOld[nbCell+2],float huOld[nbCell+2],float Z[nbCell+2],float rat); /* This function updates the water height (hNew) and the discharge (huNew) at time n+1 with respect to th water height (hOld) and the discharge (huOld) at time n. The argument rat is the ration between dt and dx. */
	void flux_hKinetic(float h[nbCell+2],float hu[nbCell+2],float fh[nbCell+2],float Z[nbCell+2]); /* This function computes the flux for the equation of mass conservation. In this case the flux is the Rusanov flux  */
	void flux_huKinetic(float h[nbCell+2],float hu[nbCell+2],float fhu[nbCell+2],float Z[nbCell+2]); /* This function computes the flux for the equation of momentum conservation. In this case the flux is the Rusanov flux. */
	void computeIntegralsFPlus_h(float a[nbCell+2],float Fr[nbCell+2], float I[3][nbCell+2],float KPlus[nbCell+2],float KMinus[nbCell+2]);/* computes the 3 integrals for the numercial flux F(i+0.5) for the mass conservation equation for the kinetic solver. */
	void computeIntegralsFMinus_h(float a[nbCell+2],float Fr[nbCell+2], float I[3][nbCell+2],float KPlus[nbCell+2],float KMinus[nbCell+2]);/* computes the 3 integrals for the numercial flux F(i-0.5) for the mass conservation equation for the kinetic solver. */
	void computeIntegralsFPlus_hu(float a[nbCell+2],float Fr[nbCell+2], float I[3][nbCell+2],float KPlus[nbCell+2],float KMinus[nbCell+2],float K[nbCell+2]);/* computes the 3 integrals for the numercial flux F(i+0.5) for the momentum conservation equation for the kinetic solver. */
	void computeIntegralsFMinus_hu(float a[nbCell+2],float Fr[nbCell+2], float I[3][nbCell+2],float KPlus[nbCell+2],float KMinus[nbCell+2],float K[nbCell+2]);/* computes the 3 integrals for the numercial flux F(i-0.5) for the momentum conservation equation for the kinetic solver. */
	float getFroudeNumber(float hu,float h);/* computes the Froude number. */
	/* The following procedure computes
	 * if face = 0
	 * 		if coef = 1
	 * 			kPlus[i] = D(Z(i+0.5))_+/h[i+1] where D(Z(i+0.5))_+ = max(0,Z(i+1)-Z(i))
	 * 		else, if coef = -1
	 * 			kPlus[i] = D(Z(i+0.5))_-/h[i+1] where D(Z(i+0.5))_- = max(0,-(Z(i+1)-Z(i)))
	 * 		end if
	 * else, if face = -1
	 * 		if coef = 1
	 * 			kPlus[i] = D(Z(i-0.5))_+/h[i] where D(Z(i-0.5))_+ = max(0,Z(i)-Z(i-1))
	 * 		else, if coef = -1
	 * 			kPlus[i] = D(Z(i-0.5))_-/h[i] where D(Z(i-0.5))_- = max(0,-(Z(i)-Z(i-1)))
	 * 		end if
	 * end if
	 * Note that face = 0 or face = -1 and coef = -1 or coef = 1.
	 * If you put an other value for "face" or "coef" you will have mistakes!!! */
	void getKPlus(float KPlus[nbCell+2],float h[nbCell+2],float Z[nbCell+2],int face,float coef);
	/* The following procedure computes
	 * if face = 0
	 * 		if coef = 1
	 * 			kMinus[i] = D(Z(i+0.5))_+/h[i] where D(Z(i+0.5))_+ = max(0,Z(i+1)-Z(i))
	 * 		else, if coef = -1
	 * 			kMinus[i] = D(Z(i+0.5))_-/h[i] where D(Z(i+0.5))_- = max(0,-(Z(i+1)-Z(i)))
	 * 		end if
	 * else, if face = -1
	 * 		if coef = 1
	 * 			kMinus[i] = D(Z(i-0.5))_+/h[i-1] where D(Z(i-0.5))_+ = max(0,Z(i)-Z(i-1))
	 * 		else, if coef = -1
	 *  		kMinus[i] = D(Z(i-0.5))_-/h[i-1] where D(Z(i-0.5))_- = max(0,-(Z(i)-Z(i-1)))
	 * 		end if
	 * end if
	 * Note that face = 0 or face = -1 and coef = -1 or coef = 1.
	 * If you put an other value for "face" or "coef" you will have mistakes!!!*/
	void getKMinus(float KMinus[nbCell+2],float h[nbCell+2],float Z[nbCell+2],int face,float coef);
	float f(float omega,float Fr, float K); /* This function is needed to compute the integrals I_3 for each flux */
	void getK(float K[nbCell+2],float h[nbCell+2],float Z[nbCell+2],int face); /* face = 1 or face = -1 */
#endif
