#include <cmath>
#include <cfloat>

namespace MathFonctions
{
    class HypergeometriqueConfluente
    {
		public: 
			static const double PI; 
			static const double M_log10e;
			static const double c_b8;
			static const double c_b53;
			static const double c_b65; 
			static const int c__7; 
			static const int c__1; 
			static const int c__9;
			static const int c__3; 
			static const int c__2;
			typedef struct {double r, i; } complex;

			//http://plasimo.phys.tue.nl/TBCI/online-docu/html/_t_o_m_s__707_8_c-source.html
			 /*      ALGORITHM 707, COLLECTED ALGORITHMS FROM ACM. */
			 /*      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
			 /*      VOL. 18, NO. 3, SEPTEMBER, 1992, PP. 345-349. */
			 /*     **************************************************************** */
			 /*     *                                                              * */
			 /*     *      SOLUTION TO THE CONFLUENT HYPERGEOMETRIC FUNCTION       * */
			 /*     *                                                              * */
			 /*     *                           by                                 * */
			 /*     *                                                              * */
			 /*     *                      MARK NARDIN,                            * */
			 /*     *                                                              * */
			 /*     *              W. F. PERGER and ATUL BHALLA                    * */
			 /*     *                                                              * */
			 /*     *                                                              * */
			 /*     *  Michigan Technological University, Copyright 1989           * */
			 /*     *                                                              * */
			 /*     *                                                              * */
			 /*     *  Description : A numerical evaluator for the confluent       * */
			 /*     *    hypergeometric function for complex arguments with large  * */
			 /*     *    magnitudes using a direct summation of the Kummer series. * */
			 /*     *    The method used allows an accuracy of up to thirteen      * */
			 /*     *    decimal places through the use of large real arrays       * */
			 /*     *    and a single final division.  LNCHF is a variable which   * */
			 /*     *    selects how the result should be represented.  A '0' will * */
			 /*     *    return the value in standard exponential form.  A '1'     * */
			 /*     *    will return the LOG of the result.  IP is an int      * */
			 /*     *    variable that specifies how many array positions are      * */
			 /*     *    desired (usually 10 is sufficient).  Setting IP=0 causes  * */
			 /*     *    the program to estimate the number of array positions.    * */
			 /*     *                                                              * */
			 /*     *    The confluent hypergeometric function is the solution to  * */
			 /*     *    the differential equation:                                * */
			 /*     *                                                              * */
			 /*     *             zf"(z) + (a-z)f'(z) - bf(z) = 0                  * */
			 /*     *                                                              * */
			 /*     *                                                              * */
			 /*     **************************************************************** */

			 static void conhyp_(complex *ret_val, const complex *a, const complex *b, const complex *z__, const int *lnchf, const int *ip);
 			 static double CHGM(double a, double b, double x);
			 static double taylor1f1(double a,double b, double z, double tol, double erreur);
	private:
			 static double gamma(double x);
			 static void chgf_(complex *ret_val, const complex *a, const complex *b, const complex *z__, const int *l, const int *lnchf);
			 static int aradd_(const double *a, const double *b, double *c__, const int *l, const double *rmax);
			 static int arsub_(const double *a, const double *b, double *c__, const int *l, const double *rmax);
			 static int armult_(const double *a, const double *b, double *c__, const int *l, const double *rmax);
			 static int arydiv_(const double *ar, const double *ai, const double *br, const double *bi, complex *c__, const int *l, const int *lnchf, const double *rmax, const int *bit); 
			 static int cmpadd_(const double *ar, const double *ai, const double *br, const double *bi, double *cr, double *ci, const int *l, const double *rmax);
			 static int cmpmul_(const double *ar, const double *ai, const double *br, const double *bi, double *cr, double *ci, const int *l, const double *rmax);
			 static int cmpsub_(double *ar, double *ai, double *br, double *bi, double *cr, double *ci, int *l, double *rmax);
			 static int bits_(void);
			 static double store_(double *x);
			 static int conv12_(complex *cn, double *cae), conv21_(double *cae, complex *cn);
			 static int emult_(double *n1, double *e1, double *n2, double *e2, double *nf, double *ef);
			 static int ecpmul_(double *a, double *b, double *c__);
			 static int ecpdiv_(double *a, double *b, double *c__);
			 static int ediv_(double *n1, double *e1, double *n2, double *e2, double *nf, double *ef);
			 static int eadd_(double *n1, double *e1, double *n2, double *e2, double *nf, double *ef);
			 static int esub_(double *n1, double *e1, double *n2, double *e2, double *nf, double *ef);

			 static inline double signe(const double x)
			 {
			   return ((x < 0.) ? -1.0 : 1.0) ;
			 }
			 
			 static inline double f__cabs(double real, double imag)
			 {
				double temp;

				if(real < 0)
					real = -real;
				if(imag < 0)
					imag = -imag;
				if(imag > real){
					temp = real;
					real = imag;
					imag = temp;
				}
				if((real+imag) == real)
					return(real);

				temp = imag/real;
				temp = real*sqrt(1.0 + temp*temp);  /*overflow!!*/
				return(temp);
			 }
			 
			 static inline void c_log(complex *r, complex *z)
			 {
				double zi, zr;
				r->i = atan2(zi = z->i, zr = z->r);
				r->r = log( f__cabs(zr, zi) );
			 }
			 
			 static inline void c_sin(complex *r, complex *z)
			 {
				double zi = z->i, zr = z->r;
				r->r = sin(zr) * cosh(zi);
				r->i = cos(zr) * sinh(zi);
			 } 
			 
			 static inline double pow_di(const double *ap, const int *bp)
			 {
				double pow, x;
				int n;
				unsigned long u;

				pow = 1;
				x = *ap;
				n = *bp;

				if(n != 0)
					{
					if(n < 0)
						{
						n = -n;
						x = 1/x;
						}
					for(u = n; ; )
						{
						if(u & 01)
							pow *= x;
						if(u >>= 1)
							x *= x;
						else
							break;
						}
					}
				return(pow);
			 }

			 static inline double pow_dd(const double *ap, const double *bp)
			 {
				return(pow(*ap, *bp) );
			 }

			 static inline double z__abs(const complex *z)
			 {
				return( f__cabs( z->r, z->i ) );
			 }

			 static inline double d_imag(const complex *z)
			 {
				return(z->i);
			 }

			 static inline void z_div(complex *c, complex *a, complex *b)
			 {
				double ratio, den;
				double abr, abi, cr;

				if( (abr = b->r) < 0.)
					abr = - abr;
				if( (abi = b->i) < 0.)
					abi = - abi;
				if( abr <= abi )
					{
					if(abi == 0)
					{
						c->i = c->r = DBL_MAX;
						return;
					}
					ratio = b->r / b->i ;
					den = b->i * (1 + ratio*ratio);
					cr = (a->r*ratio + a->i) / den;
					c->i = (a->i*ratio - a->r) / den;
					}

				else
					{
					ratio = b->i / b->r ;
					den = b->r * (1 + ratio*ratio);
					cr = (a->r + a->i*ratio) / den;
					c->i = (a->i - a->r*ratio) / den;
					}
				c->r = cr;
			 }

			 static inline double d_nint(double *x)
			 {
				return( (*x)>=0 ? floor(*x + .5) : -floor(.5 - *x) );
			 }

			 static inline double d_int(double *x)
			 {
				return( (*x>0) ? floor(*x) : -floor(- *x) );
			 }

			 static inline double d_sign(const double *a, const double *b)
			 {
				double x;
				x = (*a >= 0 ? *a : - *a);
				return( *b >= 0 ? x : -x);
			 }

			 static inline double d_lg10(const double *x)
			 {
				return( M_log10e * log(*x) );
			 }

			 static inline double dmax (const double a, const double b)
			 {
				return ((a) >= (b) ? (a) : (b));			 
			 }
    };
}

