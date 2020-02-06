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

			 Bodies of functions and procedures for test cases

   ============================================================================= */
#include "testCase.h"

void Ritter(float h[nbCell+2],float hu[nbCell+2],float dx, float t,float hg)
{
	/* we assume that the dam is at the point x = 0;
	 * hg is the water height at the left side of the point x0.
	 * notice than t > 0 */
	 
	 float cg,xA,xB,x; /* some variables needed for the computation. */
	 int i;
	 
	 cg = sqrt(g*hg);
		 
	 xA = -cg*t;
	 xB = 2.0*cg*t;
	 
	 for(i=1;i<=nbCell;i++)
	 { 
		 x = (bInf + dx/2.0) + i*dx;
		 
		 if(x<=xA)
		 {
			 h[i] = hg;
			 hu[i] = 0.0;
		 }
		 else
		 {
			 if((x>=xA)&&(x<=xB))
			 {
				 h[i] = (4.0/(9.0*g))*(cg - x/(2.0*t))*(cg - x/(2.0*t));
				 hu[i] = h[i]*(2.0/3.0)*(x/t + cg);
			 }
			 else
			 {
				 if(x>xB)
				 {
					 h[i] = 0.0;
					 hu[i] = 0.0;
				 }
			 }
		 }
	}
	
	h[0] = h[1];
	h[nbCell+1] = h[nbCell];
	
	hu[0] = hu[1];
	hu[nbCell+1] = hu[nbCell];
}
