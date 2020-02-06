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

		   Bodies of functions and procedures for the kinetic solver

   ============================================================================= */
#include "kineticSolver.h"

void updateKinetic(float hNew[nbCell+2],float huNew[nbCell+2],float hOld[nbCell+2],float huOld[nbCell+2],float Z[nbCell+2], float rat)
{
	int i;
	float fh[nbCell+2];
	float fhu[nbCell+2];
	
	for(i=0;i<=nbCell+1;i++) /* we initialize the vectors fh and fhu before computing. */
	{
		fh[i] = 0.0;
		fhu[i] = 0.0;
	}
	
	/* Before updating hNew and huNew we compute the differents flux. */
	flux_hKinetic(hOld,huOld,fh,Z);
	flux_huKinetic(hOld,huOld,fhu,Z);

	/*for(i=1;i<=nbCell;i++)
	{
		printf("fh[%d] = %f		fhu[%d] = %f\n",i,fh[i],i,fhu[i]);
	}
	printf("\n\n");*/
	for(i=1;i<=nbCell;i++)
	{
		hNew[i] = hOld[i] - rat*fh[i];
		huNew[i] = huOld[i] - rat*fhu[i];
	}
}

void flux_hKinetic(float h[nbCell+2],float hu[nbCell+2],float fh[nbCell+2],float Z[nbCell+2])
{
	int i,j;
	float FPlus[nbCell+2]; /* it represents the flux at time tat the point i + 0.5 */
	float FMinus[nbCell+2]; /* it represents the flux at time tat the point i - 0.5 */
	float FrNumber[nbCell+2]; /* this vector will contain all Froude numbers (for each cell). */
	float IPlus[3][nbCell+2]; /* This vectors reresents the 3 integrals for the numercial flux FPlus. */
	float IMinus[3][nbCell+2]; /* This vectors reresents the 3 integrals for the numercial flux FMinus. */
	float alpha[nbCell+2]; /* alpha[i] = min(1,max(-1,-FrNumber[i])) */
	float hPower[nbCell+2]; /* This vector contains h^(3/2) */
	float KPlus[nbCell+2];
	float KMinus[nbCell+2];
	
	for(i=0;i<=nbCell+1;i++) /* we initialize the vectors FPlus, FrNumber and FMinus before computing. */
	{
		FPlus[i] = 0.0;
		FMinus[i] = 0.0;
		KPlus[i] = 0.0;
		KMinus[i] = 0.0;
		FrNumber[i] = getFroudeNumber(hu[i],h[i]);/* Before computing FPlus and FMinus, we need to compute the Froude number. */
		alpha[i] = min(1.0,max(-1.0,-FrNumber[i]));/* Then we compute the vector alpha for the flux FPlus*/
		hPower[i] = power3Over2(h[i]); /* we compute de vector hPower = h^(3/2). */
	}
	
	for(i=0;i<=2;i++) /* we initialize the vectors IPlus IMinus before computing. */
	{
		for(j=0;j<=nbCell+1;j++)
		{
			IPlus[i][j] = 0.0;
			IMinus[i][j] = 0.0;
		}
	}
	
	/* The last step before comuting integrals is to calculate KPlus and KMinus */
	/* First for the computation of integrals for the numerical flux FPlus...*/
	getKPlus(KPlus,h,Z,0,-1.0); /* here we compute kPlus[i] = D(Z(i+0.5))/h[i+1] where D(Z(i+0.5)) = Z(i+1)-Z(i) */
	getKMinus(KMinus,h,Z,0,1.0);/* here we compute kMinus[i] = D(Z(i+0.5))/h[i] where D(Z(i+0.5)) = Z(i+1)-Z(i) */
	/* ...and we compute the 3 integrals needed for computing de numerical flux FPlus. */
	computeIntegralsFPlus_h(alpha,FrNumber,IPlus,KPlus,KMinus);
	
	for(i=0;i<=nbCell+1;i++) /* Before computing integrals for the flux FMinus we need to recompute alpha*/
	{
		alpha[i] = max(-1.0,min(1.0,-FrNumber[i]));
	}

	/* and then for the computation of integrals for the numerical flux FMinus ...*/
	getKPlus(KPlus,h,Z,1,1.0);/* here we compute kPlus[i] = D(Z(i-0.5))/h[i] where D(Z(i-0.5)) = Z(i)-Z(i+1) */
	getKMinus(KMinus,h,Z,1,-1.0);/* here we compute kPlus[i] = D(Z(i-0.5))/h[i-1] where D(Z(i-0.5)) = Z(i)-Z(i+1) */
	/* and we compute the 3 integrals needed for computing de numerical flux FMinus. */
	computeIntegralsFMinus_h(alpha,FrNumber,IMinus,KPlus,KMinus);
	
	for(i=1;i<=nbCell;i++)
	{
		FPlus[i] = hPower[i]*IPlus[0][i] + hPower[i]*IPlus[1][i] + hPower[i+1]*IPlus[2][i];
		FMinus[i] = hPower[i]*IMinus[0][i] + hPower[i]*IMinus[1][i] + hPower[i-1]*IMinus[2][i];
	}

	/* In this case we don't mind about the flux at cell number 0 and n+1. That is why we put the value 0.0. */
	for(i=1;i<=nbCell;i++)
	{
		fh[i] = (2.0*sqrt(2.0*g)*(FPlus[i] - FMinus[i]))/PI;
	}
}

void flux_huKinetic(float h[nbCell+2],float hu[nbCell+2],float fhu[nbCell+2],float Z[nbCell+2])
{
	int i,j;
	float FPlus[nbCell+2]; /* it represents the flux at time tat the point i + 0.5 */
	float FMinus[nbCell+2]; /* it represents the flux at time tat the point i - 0.5 */
	float FrNumber[nbCell+2]; /* this vector will contain all Froude numbers (for each cell). */
	float IPlus[3][nbCell+2]; /* This vectors reresents the 3 integrals for the numercial flux FPlus. */
	float IMinus[3][nbCell+2]; /* This vectors reresents the 3 integrals for the numercial flux FMinus. */
	float alpha[nbCell+2]; /* alpha[i] = min(1,max(-1,-FrNumber[i])) */
	float hSquare[nbCell+2]; /* This vector contains h^2 */
	float KPlus[nbCell+2];
	float KMinus[nbCell+2];
	float K[nbCell+2]; /* To compute the ratio dZ/h[i] ou dZ/h[i+1] without the max in dZ.*/
	
	for(i=0;i<=nbCell+1;i++) /* we initialize the vectors FPlus, FrNumber and FMinus before computing. */
	{
		FPlus[i] = 0.0;
		FMinus[i] = 0.0;
		KPlus[i] = 0.0;
		KMinus[i] = 0.0;
		K[i] = 0.0;
		FrNumber[i] = getFroudeNumber(hu[i],h[i]);/* Before computing FPlus and FMinus, we need to compute the Froude number. */
		alpha[i] = min(1.0,max(-1.0,-FrNumber[i]));/* Then we compute the vector alpha */
		hSquare[i] = h[i]*h[i]; /* we compute de vector hPower = h^2. */
	}
	for(i=0;i<=2;i++) /* we initialize the vectors IPlus IMinus before computing. */
	{
		for(j=0;j<=nbCell+1;j++)
		{
			IPlus[i][j] = 0.0;
			IMinus[i][j] = 0.0;
		}
	}
	
	/* The last step before comuting integrals is to calculate KPlus and KMinus */
	/* First for the computation of integrals for the numerical flux FPlus...*/
	getKPlus(KPlus,h,Z,0,-1.0); /* here we compute kPlus[i] = D(Z(i+0.5))/h[i+1] where D(Z(i+0.5)) = Z(i+1)-Z(i) */
	getKMinus(KMinus,h,Z,0,1.0);/* here we compute kMinus[i] = D(Z(i+0.5))/h[i] where D(Z(i+0.5)) = Z(i+1)-Z(i) */
	getK(K,h,Z,0);
	/* ...and we compute the 3 integrals needed for computing de numerical flux FPlus. */
	computeIntegralsFPlus_hu(alpha,FrNumber,IPlus,KPlus,KMinus,K);
	
	for(i=0;i<=nbCell+1;i++) /* Before computing integrals for the flux FMinus we need to recompute alpha*/
	{
		alpha[i] = max(-1.0,min(1.0,-FrNumber[i]));
	}
	/* and then for the computation of integrals for the numerical flux FMinus ...*/
	getKPlus(KPlus,h,Z,1,1.0);/* here we compute kPlus[i] = D(Z(i-0.5))/h[i] where D(Z(i-0.5)) = Z(i)-Z(i+1) */
	getKMinus(KMinus,h,Z,1,-1.0);/* here we compute kPlus[i] = D(Z(i-0.5))/h[i-1] where D(Z(i-0.5)) = Z(i)-Z(i+1) */
	getK(K,h,Z,1);
	/* and we compute the 3 integrals needed for computing de numerical flux FMinus. */
	computeIntegralsFMinus_hu(alpha,FrNumber,IMinus,KPlus,KMinus,K);

	for(i=1;i<=nbCell;i++)
	{
		FPlus[i] = hSquare[i]*IPlus[0][i] + hSquare[i]*IPlus[1][i] + hSquare[i+1]*IPlus[2][i];
		FMinus[i] = hSquare[i]*IMinus[0][i] + hSquare[i]*IMinus[1][i] + hSquare[i-1]*IMinus[2][i];
	}
	
	/* In this case we don't mind about the flux at cell number 0 and n+1. That is why we put the value 0.0. */
	for(i=1;i<=nbCell;i++)
	{
		fhu[i] = 4.0*g*(FPlus[i] - FMinus[i])/PI;
	}
}

void computeIntegralsFPlus_h(float a[nbCell+2],float Fr[nbCell+2], float I[3][nbCell+2],float KPlus[nbCell+2],float KMinus[nbCell+2])
{
	int i;
	float beta = 0.0;
	
	for(i=1;i<=nbCell;i++)
	{
		I[0][i] = (1.0/3.0)*(1-a[i]*a[i])*sqrt(1-a[i]*a[i]) + 0.5*Fr[i]*acos(a[i]) - 0.5*Fr[i]*a[i]*sqrt(1-a[i]*a[i]);
		/* we compute beta for I[1][i] */
		beta = max(-1,min(sqrt(KMinus[i])-Fr[i],1));
		I[1][i] = (1.0/3.0)*(1-beta*beta)*sqrt(1-beta*beta) - (1.0/3.0)*(1-a[i]*a[i])*sqrt(1-a[i]*a[i]) + 0.5*Fr[i]*(acos(beta)-acos(a[i])) + 0.5*Fr[i]*(a[i]*sqrt(1-a[i]*a[i]) - beta*sqrt(1-beta*beta));
		/* we compute beta for I[1][i] */
		beta = max(-1,min(-sqrt(KPlus[i])-Fr[i+1],1));
		I[2][i] = -(1.0/3.0)*(1-beta*beta)*sqrt(1-beta*beta) + 0.5*Fr[i+1]*(PI - acos(beta)) + 0.5*Fr[i+1]*beta*sqrt(1-beta*beta);
	}	
}

void computeIntegralsFMinus_h(float a[nbCell+2],float Fr[nbCell+2], float I[3][nbCell+2],float KPlus[nbCell+2],float KMinus[nbCell+2])
{
	int i;
	float beta = 0.0;

	for(i=1;i<=nbCell;i++)
	{
		I[0][i] = -(1.0/3.0)*(1-a[i]*a[i])*sqrt(1-a[i]*a[i]) + 0.5*Fr[i]*(PI - acos(a[i])) + 0.5*Fr[i]*a[i]*sqrt(1-a[i]*a[i]);
		/* we compute beta for I[1][i] */
		beta = min(1,max(-sqrt(KPlus[i])-Fr[i],-1));
		I[1][i] = (1.0/3.0)*(1-beta*beta)*sqrt(1-beta*beta) - (1.0/3.0)*(1-a[i]*a[i])*sqrt(1-a[i]*a[i]) + 0.5*Fr[i]*(acos(beta)-acos(a[i])) + 0.5*Fr[i]*(a[i]*sqrt(1-a[i]*a[i]) - beta*sqrt(1-beta*beta));
		I[1][i] = -I[1][i];
		/* we compute beta for I[2][i] */
		beta = min(1,max(sqrt(KMinus[i])-Fr[i-1],-1));
		I[2][i] = (1.0/3.0)*(1-beta*beta)*sqrt(1-beta*beta) + 0.5*Fr[i-1]*acos(beta) - 0.5*Fr[i-1]*beta*sqrt(1-beta*beta);
	}
}

void computeIntegralsFPlus_hu(float a[nbCell+2],float Fr[nbCell+2], float I[3][nbCell+2],float KPlus[nbCell+2],float KMinus[nbCell+2],float K[nbCell+2])
{
	int i,j;
	float beta = 0.0;
	float omega = 0.0;
	for(i=1;i<=nbCell;i++)
	{
		I[0][i] = (1.0/3.0)*(1-a[i]*a[i])*sqrt(1-a[i]*a[i])*(2.0*Fr[i] + a[i]) + 0.5*acos(a[i])*(0.25 + Fr[i]*Fr[i]) + 0.5*a[i]*sqrt(1 - a[i]*a[i])*((1.0/6.0)*a[i]*a[i] - (5.0/12.0) - Fr[i]*Fr[i]);
		/* we compute beta for I[1][i] */
		/* we compute I[1][.] in 3 steps to be as clear as possible */
		beta = max(-1,min(sqrt(KMinus[i])-Fr[i],1));
		I[1][i] = (1.0/3.0)*(1-a[i]*a[i])*sqrt(1-a[i]*a[i])*(2.0*Fr[i] + a[i]) - (1.0/3.0)*(1-beta*beta)*sqrt(1-beta*beta)*(2.0*Fr[i] + beta);
		I[1][i] = I[1][i] + 0.5*(acos(a[i]) - acos(beta))*(0.25 + Fr[i]*Fr[i]);
		I[1][i] = I[1][i] + 0.5*beta*sqrt(1 - beta*beta)*(Fr[i]*Fr[i] - (1.0/6.0)*beta*beta + (5.0/12.0)) - 0.5*a[i]*sqrt(1 - a[i]*a[i])*(Fr[i]*Fr[i] - (1.0/6.0)*a[i]*a[i] + (5.0/12.0));
	}
	
	/* Finally we need to compute the approximate value of I[2][.] */
	for(i=1;i<=nbCell;i++)
	{
		/* we compute beta for I[1][i] */
		beta = max(-1,min(-sqrt(KPlus[i])-Fr[i+1],1));
		for(j=1;j<=N;j++)
		{
			omega = -1.0 + (j - 0.5)*(beta + 1.0)/N;
			I[2][i] = I[2][i] + f(omega,Fr[i+1],K[i]);
		}
		I[2][i] = ((1.0 + beta)/N)*I[2][i];
	}
}

void computeIntegralsFMinus_hu(float a[nbCell+2],float Fr[nbCell+2], float I[3][nbCell+2],float KPlus[nbCell+2],float KMinus[nbCell+2],float K[nbCell+2])
{
	int i,j;
	float beta = 0.0;
	float omega = 0.0;
	
	for(i=1;i<=nbCell;i++)
	{
		I[0][i] = -(1.0/3.0)*(1-a[i]*a[i])*sqrt(1-a[i]*a[i])*(2.0*Fr[i] + a[i]) - 0.5*acos(a[i])*(0.25 + Fr[i]*Fr[i]) - 0.5*a[i]*sqrt(1 - a[i]*a[i])*((1.0/6.0)*a[i]*a[i] - (5.0/12.0) - Fr[i]*Fr[i]);
		I[0][i] = I[0][i] + 0.5*PI*(0.25 + Fr[i]*Fr[i]) ;
		/* we compute beta for I[1][i] */
		/* we compute I[1][.] in 3 steps to be as clear as possible */
		beta = min(1,max(-sqrt(KPlus[i])-Fr[i],-1));
		I[1][i] = (1.0/3.0)*(1-beta*beta)*sqrt(1-beta*beta)*(2.0*Fr[i] + beta) - (1.0/3.0)*(1-a[i]*a[i])*sqrt(1-a[i]*a[i])*(2.0*Fr[i] + a[i]);
		I[1][i] = I[1][i] + 0.5*(acos(beta) - acos(a[i]))*(0.25 + Fr[i]*Fr[i]);
		I[1][i] = I[1][i] + 0.5*beta*sqrt(1 - beta*beta)*((1.0/6.0)*beta*beta - (5.0/12.0) - Fr[i]*Fr[i]) - 0.5*a[i]*sqrt(1 - a[i]*a[i])*((1.0/6.0)*a[i]*a[i] - (5.0/12.0) - Fr[i]*Fr[i]);
	}
	
	/* Finally we need to compute the approximate value of I[2][.] */
	for(i=1;i<=nbCell;i++)
	{
		/* we compute beta for I[1][i] */
		beta = min(1,max(sqrt(KMinus[i])-Fr[i-1],-1));
		for(j=1;j<=N;j++)
		{
			omega = beta + (j - 0.5)*(1.0 - beta)/N;
			I[2][i] = I[2][i] + f(omega,Fr[i-1],K[i]);
		}
		I[2][i] = ((beta - 1.0)/N)*I[2][i];
	}
}

float getFroudeNumber(float hu,float h)
{
	float res;

	if(h<=epsilon)
	{
		res = 0.0;
	}
	else
	{
		res = hu/(h*sqrt(2.0*g*h));
	}
	return res;
}

void getKPlus(float KPlus[nbCell+2],float h[nbCell+2],float Z[nbCell+2],int face,float coef)
{
	int i;
	
	if(face==0) /* if we are at the inteface i + 0.5 : */
	{
		for(i=1;i<=nbCell;i++)
		{
			if(h[i+1]<=epsilon)
			{
				KPlus[i] = 0.0;
			}
			else
			{
				KPlus[i] = max(0.0,coef*(Z[i+1]-Z[i]))/h[i+1];
			}
		}
	}
	else /* esle, if we are at the inteface i - 0.5 : */
	{
		for(i=1;i<=nbCell;i++)
		{
			if(h[i]<=epsilon)
			{
				KPlus[i] = 0.0;
			}
			else
			{
				KPlus[i] = max(0.0,coef*(Z[i-1]-Z[i]))/h[i];
			}
		}
	}
}

void getK(float K[nbCell+2],float h[nbCell+2],float Z[nbCell+2],int face)
{
	int i;
	
	if(face==0) /* if we are at the inteface i + 0.5 : */
	{
		for(i=1;i<=nbCell;i++)
		{
			if(h[i+1]<=epsilon)
			{
				K[i] = 0.0;
			}
			else
			{
				K[i] = (Z[i+1] - Z[i])/h[i+1];
			}
		}
	}
	else /* esle, if we are at the inteface i - 0.5 : */
	{
		for(i=1;i<=nbCell;i++)
		{
			if(h[i-1]<=epsilon)
			{
				K[i] = 0.0;
			}
			else
			{
				K[i] = (Z[i-1] - Z[i])/h[i-1];
			}
		}
	}
}

void getKMinus(float KMinus[nbCell+2],float h[nbCell+2],float Z[nbCell+2],int face,float coef)
{
	int i;
	
	if(face==0) /* if we are at the inteface i + 0.5 : */
	{
		for(i=1;i<=nbCell;i++)
		{
			if(h[i]<=epsilon)
			{
				KMinus[i] = 0.0;
			}
			else
			{
				KMinus[i] = max(0,coef*(Z[i+1]-Z[i]))/h[i];
			}
		}
	}
	else /* esle, if we are at the inteface i - 0.5 : */
	{
		for(i=1;i<=nbCell;i++)
		{
			if(h[i-1]<=epsilon)
			{
				KMinus[i] = 0.0;
			}
			else
			{
				KMinus[i] = max(0,coef*(Z[i-1]-Z[i]))/h[i-1];
			}
		}
	}
}

float f(float omega,float Fr, float K)
{
	float res;
	
	if((omega<=-1.0)||(omega>=1.0)) /* The function (1-omega²)^(1/2) is define on a compact support [-1;1] */
	{
		res = 0.0;
	}
	else
	{
		res = -(omega + Fr)*sqrt((omega + Fr)*(omega + Fr) + K)*sqrt(1.0-omega*omega);
	}
	return res;
}
