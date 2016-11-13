/*****************************************************************************************/
/* Various line shapes used in spectroscopy (see below)
 * 
 * References: 
 *  [1] P.L. Varghese and R.K. Hanson, Appl. Opt. Vol. 23, 
 *      Pag. 2376, (1984).
 *  [2] J. Humlicek, J Quant. Spectrosc. Radiat. Transfer Vol. 21, 
 *      Pag. 309 (1979)
 * 
 * general prototype:
 *
 * int lineshape(int n, double, *x, double *y, double x0, double a,
 *               double g, double s, double b)
 *
 * the return value is 1 in case of error 0 otherwise;
 * Inputs:
 * n             = number of points
 * x             = pointer to x[n]
 * y             = pointer to y[n]
 * x0            = center
 * a             = amplitude
 * g             = lorentz width
 * s             = gauss width
 * b             = collisional frequency
 *
 * Note that only positive values for s,g,b are meaningfull.  Negative
 * values are internally converted to absolute value so there is no
 * need to set limit on the parameters. An occasional fit resulting in
 * a negative parameter, say s=-1.2 should be interpreted as s=|-1.2|
 *
 * Normalization conditon is f(x0)=a.
 *
 * The functions xxxnorm are availble if an area normalization is required.
 * The area under the curve is given by a*xxxnorm(g,s,b).
 * y[i]+=f(x[i]) so y[i] should be cleared before calling the first time
 *
 * Of the 7 standard lineshapes discussed in [1] only 5 are
 * implemented here.  The other two require setting b as a complex
 * quantity (in order to include collisional shift). 
 * If you really need them should be straightforward to modify the code.
 *
 * Lineshapes:
 *
 * - Lorentzian -> lorentz()
 *
 *   computed as g^2/(dx^2+g^2), where dx=x-x0, so g=gamma and g=HWHM
 *
 * - Gaussian -> gauss()
 *
 *   computed as exp(-dx^2/s^2) so s=sqrt(2)*sigma or
 *   s=HWHM/sqrt(ln(2)) ~ HWHM*1.2
 * 
 * - Voigt -> voigt() or humlicek() 
 *
 *   computed as Re(W(Z))/Re(W(Z')) where W() is the complex
 *   probability function and Z=rx+i ry, Z'=i ry, rx=dx/s,
 *   ry=g/s. There are two algoritms to evaluate W(z). One good up
 *   to 6 or 7 digits (single precision) stolen from CERNLIBS, used
 *   by voigt() (see cwerf()) and another, approximate up to 4
 *   digits but faster, described in [2], used by humlicek().
 *
 * - Rautian -> rautian() or rautianh()   
 *   computed as Re(W(Z)/(1-sqrt(Pi)*rz*W(Z)))/n0
 *   with n0=Re(W(Z')/(1-sqrt(Pi)*rz*W(i(Z')),
 *        Z=rx+i(ry+rz),
 *        Z'=i(ry+rz),
 *        rz=b/s
 *   
 *  rautian() uses the slow, accurate algorithm for W(), rautianh() the 
 *  faster, approximate																	//*	
 *																						//*
 * - Galatry -> galatry()																//*
 *   see [1] for a description of the algorithm											//*	
*/																						//*
																						//*
#define M_PI          3.1415926535897932384626433832795									//*
#define M_2_SQRTPI    1.7724538509055160272981674833411									//*
																						//*	
//#include<math.h>
#include<cmath>                                                                         //*
																						//*
int lorentz(long n, double *x, double *y, double x0, double a,							//*
	    double g, double s, double b)													//*
{																						//*
  long i;																				//*		
  double g2=g*g;																		//*		
																						//*
  if(g2==0.)																			//*
    return 1;																			//*
  a*=g2;																				//*
  for(i=0;i<n;i++)																		//*
    y[i]+=a/((x[i]-x0)*(x[i]-x0)+g2);													//*
																						//*
  return 0;																				//*
}																						//*
																						//*
double lornorm(double g,double s, double b)												//*
{																						//*				
  return g*M_PI;																		//*	
}																						//*
																						//*
int gauss(long n, double *x, double *y, double x0, double a,							//*
	  double g, double s, double b)														//*
{																						//*	
  long i;																				//*
  double s2=s*s;																		//*	
																						//*
  if(s2==0.)																			//*
    return 1;																			//*
  for(i=0;i<n;i++)																		//*
    y[i]+=a*exp(-(x[i]-x0)*(x[i]-x0)/s2);												//*	
																						//*
  return 0;																				//*
}																						//*
																						//*		
double gaunorm(double g,double s, double b)												//*
{																						//*
  return s*sqrt(M_PI);																	//*												
}																						//*
																						//*
/*																						//*
   from CERN mathlib cwerf64.F															//*
																						//*
*/																						//*
																						//*
void cwerf(double x, double y, double *vr, double *vi)									//*
{																						//*
#define CVOIGT 1.12837916709551257						fffffff								//*
																						//*
  int n;																				//*
  double xa,ya;																			//*
  double zhr,zhi;																		//*
  double rr[38],ri[38];																	//*
  double tr,ti;																			//*
  double r,xl;																			//*
  double sr,si;																			//*
																						//*
  xa=fabs(x);																			//*
  ya=fabs(y);																			//*
  if (ya<7.4 && xa<8.3)																	//*
    {																					//*	
      zhr=ya+1.6;																		//*
      zhi=xa;																			//*
      rr[37]=0;																			//*
      ri[37]=0;																			//*
      for (n=36;n>0;n--)																//*
	{																					//*
	  tr=zhr+n*rr[n+1];																	//*
	  ti=zhi-n*ri[n+1];																	//*
	  xl=0.5/(tr*tr+ti*ti);																//*
	  rr[n]=tr*xl;																		//*
	  ri[n]=ti*xl;																		//*
	}																					//*
      xl=46768052394589056.;															//*
      sr=0;																				//*
      si=0;																				//*
      for (n=33;n>0;n--)																//*
	{																					//*
	  xl*=.3125;																		//*
	  tr=rr[n]*(sr+xl)-ri[n]*si;														//*
	  si=rr[n]*si+ri[n]*(sr+xl);														//*
	  sr=tr;																			//*
	}																					//*
      *vr=sr*CVOIGT;																	//*
      *vi=si*CVOIGT;																	//*
    }																					//*
  else																					//*
    {																					//*				
      zhr=ya;																			//*
      zhi=xa;																			//*
      rr[0]=0.;																			//*
      ri[0]=0.;																			//*
      for(n=9;n>0;n--)																	//*		
	{																					//*
	  tr=zhr+n*rr[0];																	//*
	  ti=zhi-n*ri[0];																	//*
	  xl=0.5/(tr*tr+ti*ti);																//*
	  rr[0]=tr*xl;																		//*
	  ri[0]=ti*xl;																		//*
	}																					//*
      *vr=rr[0]*CVOIGT;																	//*
      *vi=ri[0]*CVOIGT;																	//*
    }																					//*
  if (ya==0.)																			//*
    *vr=exp(-xa*xa);																	//*
  if(y<0.)																				//*	
    {																					//*
      tr=xa*xa-ya*ya;																	//*
      ti=2.*xa*ya;																		//*
      r=sqrt(tr*tr+ti*ti);																//*
      xl=exp(-r);																		//*
      sr=2.*xl*tr/r;																	//*
      si=-2.*xl*ti/r;																	//*
      *vr=2.*sr-*vr;																	//*	
      *vi=2.*si-*vi;																	//*
      if(x>0.)																			//*
	*vi=-*vi;																			//*		
    }																					//*
  else																					//*
    {																					//*
      if (x<0)																			//*
	*vi=-*vi;																			//*
    }																					//*	
  return ;																				//*
}																						//*
																						//*
int voigt(long n, double *x, double *y, double x0, double a,							//*
	  double g, double s, double b)														//*
{																						//*
  long i;																				//*
  double rx, ry;																		//*
  double vr,vi;																			//*	
																						//*
  if(g==0.)																				//*	
    return gauss(n,x,y,x0,a,g,s,b);														//*		
  if(s==0.)																				//*
    return lorentz(n,x,y,x0,a,g,s,b);													//*
																						//*	
  ry=fabs(g/s);																			//*
  cwerf(0.,ry,&vr,&vi);																	//*
  a/=vr;																				//*
																						//*
  for(i=0;i<n;i++)																		//*
    {																					//*
      rx=(x[i]-x0)/s;																	//*
      cwerf(rx,ry,&vr,&vi);																//*
      y[i]+=a*vr;																		//*
    }																					//*
																						//*
  return 0;																				//*												
}																						//*
																						//*
double voinorm(double g, double s, double b)											//*
{																						//*						
  double a, vr,vi;																		//*
																						//*
  if(g==0.)																				//*	
    return gaunorm(g,s,b);																//*
  if(s==0.)																				//*	
    return lornorm(g,s,b);																//*
																						//*	
  a=fabs(g/s);																			//*
  cwerf(0.,a,&vr,&vi);																	//*
  return sqrt(M_PI)*s/vr;																//*
}																						//*
																						//*	
//poor man's complex arithmetics														//*				
double approx1(double r, double i, double *ri)											//*	
{																						//*
  double rn, in, rd, id;																//*
  double m2;																			//*
																						//*	
  rn=.5641896*r;																		//*	
  in=.5641896*i;																		//*	
																						//*	
  rd=r*r-i*i+.5;																		//*	
  id=2*r*i;																				//*
																						//*
  m2=(rd*rd+id*id);																		//*
  *ri=(rn*id-in*rd)/m2;																	//*
																						//*
  return (rn*rd+in*id)/m2;																//*
}																						//*
																						//*
double approx2(double r, double i, double *ri)											//*
{																						//*
  double rn, in, rd, id, ru, iu, ra;													//*
  double m2;																			//*
																						//*
  ru=r*r-i*i;																			//*
  iu=2*r*i;																				//*
																						//*
  rn=.5641896*ru+1.410474;																//*
  in=.5641896*iu;																		//*	
  ra=rn*r-in*i;																			//*
  in=rn*i+in*r;																			//*
  rn=ra;																				//*
																						//*
  rd=ru+3;																				//*	
  id=iu;																				//*	
  ra=rd*ru-id*iu+0.75;																	//*
  id=rd*iu+ru*id;																		//*
																						//*	
  rd=ra;																				//*
  m2=(rd*rd+id*id);																		//*
  *ri=(rn*id-in*rd)/m2;																	//*
																						//*
  return (rn*rd+in*id)/m2;																//*
}																						//*
																						//*
double approx3(double r, double i, double *ri)											//*
{																						//*
  double rn, in, rd, id, ra;															//*
  double m2;																			//*
																						//*
  rn=.5642236*r+3.778987;																//*
  in=.5642236*i;																		//*
  ra=r*rn-i*in+11.96482;																//*
  in=r*in+i*rn;																			//*
  rn=ra;																				//*	
  ra=r*rn-i*in+20.20933;																//*
  in=r*in+i*rn;																			//*
  rn=ra;																				//*
  ra=r*rn-i*in+16.4955;																	//*
  in=r*in+i*rn;																			//*	
  rn=ra;																				//*
																						//*
  rd=r+6.699398;																		//*
  id=i;																					//*
  ra=r*rd-i*id+21.69274;																//*	
  id=r*id+i*rd;																			//*
  rd=ra;																				//*	
  ra=r*rd-i*id+39.27121;																//*	
  id=r*id+i*rd;																			//*
  rd=ra;																				//*	
  ra=r*rd-i*id+38.82363;																//*
  id=r*id+i*rd;																			//*
  rd=ra;																				//*
  ra=r*rd-i*id+16.4955;																	//*
  id=r*id+i*rd;																			//*
																						//*	
																						//*
  rd=ra;																				//*
  m2=(rd*rd+id*id);																		//*
  *ri=(rn*id-in*rd)/m2;																	//*
																						//*
  return (rn*rd+in*id)/(rd*rd+id*id);													//*
}																						//*
																						//*
double approx4(double r, double i, double *ri)											//*
{																						//*
  double rn, in, rd, id, ru, iu, ra;													//*
  double m2;																			//*
																						//*
  ru=-r*r+i*i;																			//*
  iu=-2*r*i;																			//*
																						//*
  rn=.56419*ru+1.320522;																//*
  in=.56419*iu;																			//*
  ra=rn*ru-in*iu+35.7668;																//*
  in=rn*iu+in*ru;																		//*
  rn=ra;																				//*
  ra=rn*ru-in*iu+219.031;																//*
  in=rn*iu+in*ru;																		//*
  rn=ra;																				//*
  ra=rn*ru-in*iu+1540.787;																//*
  in=rn*iu+in*ru;																		//*
  rn=ra;																				//*
  ra=rn*ru-in*iu+3321.99;																//*
  in=rn*iu+in*ru;																		//*
  rn=ra;																				//*	
  ra=rn*ru-in*iu+36183.31 ;																//*
  in=rn*iu+in*ru;																		//*
  rn=ra;																				//*
  ra=rn*r-in*i;																			//*		
  in=rn*i+in*r;																			//*
  rn=ra;																				//*
																						//*			
																						//*
  rd=ru+1.84144;																		//*
  id=iu;																				//*
  ra=rd*ru-id*iu+61.5704;																//*
  id=rd*iu+ru*id;																		//*
  rd=ra;																				//*	
  ra=rd*ru-id*iu+364.219;																//*
  id=rd*iu+ru*id;																		//*
  rd=ra;																				//*
  ra=rd*ru-id*iu+2186.18;																//*
  id=rd*iu+ru*id;																		//*
  rd=ra;																				//*
  ra=rd*ru-id*iu+9022.3;																//*
  id=rd*iu+ru*id;																		//*
  rd=ra;																				//*
  ra=rd*ru-id*iu+24322.8;																//*
  id=rd*iu+ru*id;																		//*
  rd=ra;																				//*
  ra=rd*ru-id*iu+32066.6;																//*	
  id=rd*iu+ru*id;																		//*
																						//*	
  rd=ra;																				//*		
  m2=(rd*rd+id*id);																		//*
  *ri=(rn*id-in*rd)/m2;																	//*
																						//*
  return (rn*ra+in*id)/m2;																//*
}																						//*
																						//*
double sqexp(double r, double i, double *im)											//*
{																						//*
  double r2,i2,e,s,c;																	//*	
																						//*
  r2=r*r-i*i;																			//*
  i2=2*r*i;																				//*
																						//*
  e=exp(r2);																			//*
#ifdef __GNU__																			//*
  sincos(i2,&s,&c);																		//*
#else																					//*
  s=sin(i2);																			//*	
  c=cos(i2);																			//*
#endif																					//*
  *im=e*s;																				//*
																						//*
  return e*c;																			//*
}																						//*
																						//*
void singhumlicek(double rx, double ry, double *vr, double *vi)							//*	
{																						//*
  double ax,sc,vt;																		//*
																						//*	
  ax=fabs(rx);																			//*
  sc=ax+ry;																				//*
  if(sc>15.)																			//*
    *vr=approx1(ry,-rx,vi);																//*
  else																					//*
    if((sc<=15.) && (sc>5.5))															//*	
      *vr=approx2(ry,-rx,vi);															//*
    else																				//*	
      if((sc<=5.5) && (ry>=(0.195*ax-0.176)))											//*
	*vr=approx3(ry,-rx,vi);																//*		
      else																				//*
	{																					//*	
	  *vr=(sqexp(ry,-rx,vi)-approx4(ry,-rx,&vt));										//*
	  *vi-=vt;																			//*
	}																					//*
}																						//*
																						//*
int humlicek(long n, double *x, double *y, double x0, double a,							//*
	     double g, double s, double b)													//*	
{																						//*
  long i;																				//*
  double rx, ry, vr, vi, sc, ax;														//*
  double dummy;																			//*
																						//*
  if(g==0.)																				//*
    return gauss(n,x,y,x0,a,g,s,b);														//*		
  if(s==0.)																				//*
    return lorentz(n,x,y,x0,a,g,s,b);													//*
																						//*
  ry=fabs(g/s);																			//*
  singhumlicek(0.,ry,&vr,&vi);															//*
  a/=vr;																				//*
																						//*
  if(ry>15.)																			//*	
    {																					//*
      for(i=0;i<n;i++)																	//*	
	{																					//*
	  rx=(x[i]-x0)/s;																	//*
	  y[i]+=a*approx1(ry,-rx,&dummy);													//*
	}																					//*
    }																					//*
  else if((ry<=15.) && (ry>5.5))														//*
    for(i=0;i<n;i++)																	//*
      {																					//*
	rx=(x[i]-x0)/s;																		//*
	sc=fabs(rx)+ry;																		//*
	if(sc>15.)																			//*
	  y[i]+=a*approx1(ry,-rx,&dummy);													//*													
	else																				//*
	  y[i]+=a*approx2(ry,-rx,&dummy);													//*
      }																					//*
  else if ((ry<=5.5) && (ry>0.75))														//*
    for(i=0;i<n;i++)																	//*
      {																					//*
	rx=(x[i]-x0)/s;																		//*
	sc=fabs(rx)+ry;																		//*
	if(sc>=15.)																			//*	
	  y[i]+=a*approx1(ry,-rx,&dummy);													//*	
	else																				//*
	  if(sc<5.5)																		//*
	    y[i]+=a*approx3(ry,-rx,&dummy);													//*	
	  else																				//*
	    y[i]+=a*approx2(ry,-rx,&dummy);													//*
      }																					//*	
  else																					//*
    {																					//*
      for(i=0;i<n;i++)																	//*
	{																					//*
	  rx=(x[i]-x0)/s;																	//*
	  ax=fabs(rx);																		//*
	  sc=ax+ry;																			//*
	  if(sc>15.)																		//*
	    y[i]+=a*approx1(ry,-rx,&dummy);													//*
	  else																				//*	
	    if((sc<=15.) && (sc>5.5))														//*
	      y[i]+=a*approx2(ry,-rx,&dummy);												//*	
	    else																			//*	
	      if((sc<=5.5) && (ry>=(0.195*ax-0.176)))										//*
		y[i]+=a*approx3(ry,-rx,&dummy);													//*
	      else																			//*
		y[i]+=a*(sqexp(ry,-rx,&dummy)-approx4(ry,-rx,&dummy));							//*																																																																																						
	}																					//*
    }																					//*
  return 0;																				//*
}																						//*
																						//*
																						//*
double raunorm(double g, double s, double b)											//*
{																						//*
  double ry, z, vr, vi, ar, ai, fd, fn;													//*	
																						//*
  if(s==0.)																				//*
    return lornorm(g,s,b);																//*
  if(g==0.)																				//*
    return gaunorm(g,s,b);																//*
  if(b==0.)																				//*
    return voinorm(g,s,b);																//*
																						//*
  ry=fabs((g+b)/s);																		//*
  z=sqrt(M_PI)*fabs(b/s);																//*
  cwerf(0,ry,&vr,&vi);																	//*	
  ar=(1.-z*vr);																			//*
  ai=vi*z;																				//*
  fd=ar*ar+ai*ai;																		//*
  fn=vr*ar-vi*ai;																		//*
  return sqrt(M_PI)*s*fd/fn;															//*
}																						//*
																						//*	
int rautian(int n, double *x, double *y, double x0, double a,							//*
	     double g, double s, double b)													//*
{																						//*
  long i;																				//*
  double rx, ry, z, ar, ai, fd, fn;														//*				
  double vr,vi;																			//*
																						//*
  if(b==0.)																				//*
    return voigt(n,x,y,x0,a,g,s,b);														//*
  if(s==0.)																				//*	
    return lorentz(n,x,y,x0,a,g,s,b);													//*
																						//*
  ry=fabs((g+b)/s);																		//*
  z=sqrt(M_PI)*fabs(b/s);																//*	
  cwerf(0.,ry,&vr,&vi);																	//*
  ar=(1.-z*vr);																			//*
  ai=vi*z;																				//*
  fd=ar*ar+ai*ai;																		//*
  fn=vr*ar-vi*ai;																		//*
  ar=fn/fd;																				//*
  a/=ar;																				//*
  for(i=0;i<n;i++)																		//*		
    {																					//*
      rx=(x[i]-x0)/s;																	//*
      cwerf(rx,ry,&vr,&vi);																//*
      ar=(1.-z*vr);																		//*	
      ai=vi*z;																			//*
      fd=ar*ar+ai*ai;																	//*
      fn=vr*ar-vi*ai;																	//*
      y[i]+=a*fn/fd;																	//*
    }																					//*
																						//*
  return 0;																				//*
}																						//*
																						//*	
int rautianh(int n, double *x, double *y, double x0, double a,							//*	
	     double g, double s, double b)													//*
{																						//*
  long i;																				//*
  double rx, ry, z, ar, ai, fd, fn;														//*
  double vr,vi;																			//*
																						//*
  if(b==0.)																				//*
    return humlicek(n,x,y,x0,a,g,s,b);													//*
  if(s==0.)																				//*
    return lorentz(n,x,y,x0,a,g,s,b);													//*
																						//*		
  ry=fabs((g+b)/s);																		//*
  z=sqrt(M_PI)*fabs(b/s);																//*
  singhumlicek(0.,ry,&vr,&vi);															//*
  ar=(1.-z*vr);																			//*
  ai=vi*z;																				//*
  fd=ar*ar+ai*ai;																		//*
  fn=vr*ar-vi*ai;																		//*
  ar=fn/fd;																				//*
  a/=ar;																				//*
  for(i=0;i<n;i++)																		//*
    {																					//*		
      rx=(x[i]-x0)/s;																	//*
      singhumlicek(rx,ry,&vr,&vi);														//*
      ar=(1.-z*vr);																		//*
      ai=vi*z;																			//*
      fd=ar*ar+ai*ai;																	//*
      fn=vr*ar-vi*ai;																	//*
      y[i]+=a*fn/fd;																	//*
    }																					//*
																						//*
  return 0;																				//*
}																						//*
																						//*
int galatry(long n, double *x, double *y, double x0, double a,							//*
	    double g, double s, double b)													//*
{																						//*	
  long i;																				//*
  int j, nd, rg, n2, n3;																//*
  double rx, ry, rz, ay, ax, ar, ai, br, bi, cr, ci, ct, nx;							//*
  double c[9],wr[9],wi[9];																//*
  double delta, thetar, thetai, thetai2;												//*
																						//*	
  if(b==0.)																				//*
    return voigt(n,x,y,x0,a,g,s,b);														//*
  if(s==0.)																				//*
    return lorentz(n,x,y,x0,a,g,s,b);													//*	
																						//*
  ry=fabs(g/s);																			//*
  rz=fabs(b/s);																			//*
																						//*
  //Region 1																			//*
  if((rz<0.04) && (ry<0.5))																//*
    {																					//*	
      c[0]=1.;																			//*
      c[1]=0.;																			//*
      c[2]=0.;																			//*	
      c[3]=rz/12;																		//*
      for(j=4;j<9;j++)																	//*	
	c[j]=-c[j-1]*rz/j;																	//*
      c[6]+=0.5*c[3]*c[3];																//*
      c[7]+=c[3]*c[4];																	//*
      c[8]+=c[3]*c[5]+0.5*c[4]*c[4];													//*
																						//*
      cwerf(0,ry,wr,wi);																//*
																						//*
      bi=2*ry;																			//*
      wr[1]=0.;																			//*
      wi[1]=M_2_SQRTPI-bi*wr[0];														//*
																						//*		
      nd=2;																				//*
      for(j=2;j<9;j++)																	//*
	{																					//*
	  wr[j]=bi*wi[j-1]-nd*wr[j-2];														//*	
	  wi[j]=-bi*wr[j-1]-nd*wi[j-2];														//*
	  nd+=2;																			//*		
	}																					//*
      ar=wr[0]-c[3]*wi[3]+c[4]*wr[4]+c[5]*wi[5]-										//*
	c[6]*wr[6]-c[7]*wi[7]+c[8]*wr[8];													//*	
																						//*
      //if all points are in region 1 set rg=1 else rg=0								//*
      br=2*s;																			//*
      for(i=0;i<n;i++)																	//*
	if(fabs(x[i]-x0)>br)																//*
	  break;																			//*
      rg=(i==n) ? 1 : 0;																//*
    }																					//*
  else																					//*
  //Region 2																			//*
    if(((rz>0.1) && (ry<4.*pow(rz,0.868))) || rz>5)										//*
      {																					//*
	n2=4+pow(rz,-1.05)*(1+3.*exp(-1.1*ry));												//*
	delta=0.5/rz/rz;																	//*
	thetar=delta+ry/rz;																	//*
																						//*
	ar=1/thetar;																		//*
	br=ar;																				//*
	for(j=1;j<=n2;j++)																	//*
	  {																					//*
	    br*=delta/(thetar+j);															//*
	    ar+=br;																			//*
	  }																					//*
	ar*=0.5*M_2_SQRTPI/rz;																//*
	rg=2;																				//*	
      }																					//*
  //Region 3																			//*	
    else																				//*
      {																					//*
	n3=2+37*exp(-0.6*ry);																//*		
	br=ry+n3*rz;																		//*
	nx=0.5*n3;																			//*
	ar=nx/br;																			//*
	for(j=1;j<n3;j++)																	//*
	  {																					//*								
	    br-=rz;																			//*
	    nx-=0.5;																		//*
	    ar=nx/(br+ar);																	//*
	  }																					//*
	ar=0.5*M_2_SQRTPI/(ry+ar);															//*		
	rg=3;																				//*
      }																					//*
  a/=ar;																				//*
																						//*
  switch(rg)																			//*
    {																					//*
    case 0:																				//*
      n3=2+37*exp(-0.6*ry);																//*	
      for(i=0;i<n;i++)																	//*
	{																					//*
	  rx=fabs((x[i]-x0)/s);																//*
	  if(rx>2.0)																		//*
	    {																				//*
	      br=ry+n3*rz;																	//*
	      bi=-rx;																		//*
	      nx=0.5*n3;																	//*
	      ax=nx/(br*br+bi*bi);															//*	
	      ar=br*ax;																		//*
	      ai=-bi*ax;																	//*
	      for(j=1;j<n3;j++)																//*
		{																				//*
		  br-=rz;																		//*
		  nx-=0.5;																		//*	
		  ax=nx/((br+ar)*(br+ar)+(bi+ai)*(bi+ai));										//*	
		  ar=(ar+br)*ax;																//*
		  ai=-(ai+bi)*ax;																//*
		}																				//*
	      ar=0.5*M_2_SQRTPI*(ry+ar)/((ry+ar)*(ry+ar)+(bi+ai)*(bi+ai));					//*		
	    }																				//*
	  else																				//*
	    {																				//*
	      cwerf(rx,ry,wr,wi);															//*
																						//*
	      bi=2*ry;																		//*
	      br=2*rx;																		//*
	      wr[1]=-br*wr[0]+bi*wi[0];														//*
	      wi[1]=M_2_SQRTPI-br*wi[0]-bi*wr[0];											//*		
	      nd=2;																			//*
	      for(j=2;j<9;j++)																//*		
		{																				//*
		  wr[j]=-br*wr[j-1]+bi*wi[j-1]-nd*wr[j-2];										//*
		  wi[j]=-br*wi[j-1]-bi*wr[j-1]-nd*wi[j-2];										//*
		  nd+=2;																		//*
		}																				//*
	      ar=wr[0]-c[3]*wi[3]+c[4]*wr[4]+c[5]*wi[5]-									//*
		c[6]*wr[6]-c[7]*wi[7]+c[8]*wr[8];												//*
	    }																				//*
	  y[i]+=a*ar;																		//*
	}																					//*
      break;																			//*
    case 1:																				//*
      bi=2*ry;																			//*
      for(i=0;i<n;i++)																	//*	
	{																					//*
	  rx=fabs((x[i]-x0)/s);																//*
	  cwerf(rx,ry,wr,wi);																//*
																						//*
	  br=2*rx;																			//*
	  wr[1]=-br*wr[0]+bi*wi[0];															//*
	  wi[1]=M_2_SQRTPI-br*wi[0]-bi*wr[0];												//*
																						//*
	  nd=2;																				//*
	  for(j=2;j<9;j++)																	//*
	    {																				//*
	      wr[j]=-br*wr[j-1]+bi*wi[j-1]-nd*wr[j-2];										//*
	      wi[j]=-br*wi[j-1]-bi*wr[j-1]-nd*wi[j-2];										//*	
	      nd+=2;																		//*
	    }																				//*
	  ar=wr[0]-c[3]*wi[3]+c[4]*wr[4]+c[5]*wi[5]-										//*
	    c[6]*wr[6]-c[7]*wi[7]+c[8]*wr[8];												//*	
	  y[i]+=a*ar;																		//*
	}																					//*
      break;																			//*
    case 2:																				//*
      for(i=0;i<n;i++)																	//*
	{																					//*	
	  rx=fabs((x[i]-x0)/s);																//*
	  thetai=-rx/rz;																	//*
	  thetai2=thetai*thetai;															//*
																						//*		
	  cr=thetar/(thetar*thetar+thetai2);												//*	
	  ci=-thetai/(thetar*thetar+thetai2);												//*	
	  ar=cr;																			//*
	  for(j=1;j<=n2;j++)																//*	
	    {																				//*
	      ay=delta/((thetar+j)*(thetar+j)+thetai2);										//*		
	      br=ay*(thetar+j);																//*
	      bi=-ay*thetai;																//*
	      ct=br*cr-bi*ci;																//*		
	      ci=bi*cr+br*ci;																//*
	      cr=ct;																		//*
	      ar+=cr;																		//*
	    }																				//*
	  ar*=0.5*M_2_SQRTPI/rz;															//*
	  y[i]+=a*ar;																		//*
	}																					//*
      break;																			//*
    case 3:																				//*
      for(i=0;i<n;i++)																	//*
	{																					//*
	  rx=fabs((x[i]-x0)/s);																//*
	  br=ry+n3*rz;																		//*
	  bi=-rx;																			//*
	  nx=0.5*n3;																		//*		
	  ax=nx/(br*br+bi*bi);																//*
	  ar=br*ax;																			//*
	  ai=-bi*ax;																		//*	
	  for(j=1;j<n3;j++)																	//*
	    {																				//*
	      br-=rz;																		//*
	      nx-=0.5;																		//*
	      ax=nx/((br+ar)*(br+ar)+(bi+ai)*(bi+ai));										//*
	      ar=(ar+br)*ax;																//*
	      ai=-(ai+bi)*ax;																//*
	    }																				//*
	  ar=0.5*M_2_SQRTPI*(ry+ar)/((ry+ar)*(ry+ar)+(bi+ai)*(bi+ai));						//*				
	  y[i]+=a*ar;																		//*
	}																					//*
    }																					//*
																						//*
  return 0;																				//*
}																						//*	
																						//*
double galnorm(double g, double s, double b)											//*
{																						//*
  int j, nd, n2, n3;																	//*
  double ry, rz, ar, br, bi, nx;														//*
  double c[9],wr[9],wi[9];																//*
  double delta, thetar;																	//*
																						//*
  if(b==0.)																				//*
    return voinorm(g,s,b);																//*	
																						//*
  ry=fabs(g/s);																			//*
  rz=fabs(b/s);																			//*
																						//*
  //Region 1																			//*	
  if((rz<0.04) && (ry<0.5))																//*	
    {																					//*
      c[0]=1.;																			//*
      c[1]=0.;																			//*
      c[2]=0.;																			//*
      c[3]=rz/12;																		//*
      for(j=4;j<9;j++)																	//*
	c[j]=-c[j-1]*rz/j;																	//*
      c[6]+=0.5*c[3]*c[3];																//*
      c[7]+=c[3]*c[4];																	//*	
      c[8]+=c[3]*c[5]+0.5*c[4]*c[4];													//*
																						//*
      cwerf(0,ry,wr,wi);																//*
																						//*	
      bi=2*ry;																			//*
      wr[1]=0.;																			//*
      wi[1]=M_2_SQRTPI-bi*wr[0];														//*
																						//*
      nd=2;																				//*
      for(j=2;j<9;j++)																	//*
	{																					//*
	  wr[j]=bi*wi[j-1]-nd*wr[j-2];														//*
	  wi[j]=-bi*wr[j-1]-nd*wi[j-2];														//*
	  nd+=2;																			//*	
	}																					//*
      ar=wr[0]-c[3]*wi[3]+c[4]*wr[4]+c[5]*wi[5]-										//*
	c[6]*wr[6]-c[7]*wi[7]+c[8]*wr[8];													//*	
    }																					//*	
  else																					//*
  //Region 2																			//*
    if(((rz>0.1) && (ry<4.*pow(rz,0.868))) || rz>5)										//*
      {																					//*
	n2=4+pow(rz,-1.05)*(1+3.*exp(-1.1*ry));												//*
	delta=0.5/rz/rz;																	//*
	thetar=delta+ry/rz;																	//*
																						//*
	ar=1/thetar;																		//*	
	br=ar;																				//*				
	for(j=1;j<=n2;j++)																	//*
	  {																					//*
	    br*=delta/(thetar+j);															//*
	    ar+=br;																			//*
	  }																					//*		
	ar*=0.5*M_2_SQRTPI/rz;																//*
      }																					//*		
  //Region 3																			//*
    else																				//*
      {																					//*
	n3=2+37*exp(-0.6*ry);																//*
	br=ry+n3*rz;																		//*
	nx=0.5*n3;																			//*			
	ar=nx/br;																			//*
	for(j=1;j<n3;j++)																	//*
	  {																					//*
	    br-=rz;																			//*
	    nx-=0.5;																		//*		
	    ar=nx/(br+ar);																	//*
	  }																					//*
	ar=0.5*M_2_SQRTPI/(ry+ar);															//*
      }																					//*	
  return sqrt(M_PI)*s/ar;																//*
}																						//*
																						//*									
/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/


