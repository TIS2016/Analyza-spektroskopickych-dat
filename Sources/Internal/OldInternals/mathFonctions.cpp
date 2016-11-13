#include "MathFonctions.h"

using namespace std;

struct
{
    double y;
} stcom_;

namespace MathFonctions
{

	const double HypergeometriqueConfluente::PI              = 3.1415926535897932384626433832795;
	const double HypergeometriqueConfluente::M_log10e          = 0.43429448190325182765;
	const double HypergeometriqueConfluente::c_b8              = 2.0; 
	const double HypergeometriqueConfluente::c_b53             = 1.0; 
	const double HypergeometriqueConfluente::c_b65             = 10.0; 
	const int HypergeometriqueConfluente::c__7				  = 7; 
	const int HypergeometriqueConfluente::c__1                 = 1; 
	const int HypergeometriqueConfluente::c__9                 = 9; 
	const int HypergeometriqueConfluente::c__3                 = 3; 
	const int HypergeometriqueConfluente::c__2                 = 2;

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
 /*     *  Subprograms called: BITS, CHGF                              * */
 /*     *                                                              * */
 /*     **************************************************************** */

 void HypergeometriqueConfluente::conhyp_(complex *ret_val, const complex *a, const complex *b, const complex *z__, const int *lnchf, const int *ip)
{
    /* System generated locals */
    complex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions: prototypes.h */

    /* Local variables */
    static double term1, term2;
    static int i__;
    static double nterm, fx, ang, max__;

    if (z__abs(z__) != 0.) {
	ang = atan2(d_imag(z__), (double) z__->r);
    } else {
	ang = 1.;
    }
    if (abs(ang) < 1.570796) {
	ang = 1.;
    } else {
	ang = sin(abs(ang) - 1.570796325) + 1.;
    }
    max__ = 0.;
    nterm = 0.;
    fx = 0.;
    term1 = 0.;
L10:
    nterm += 1;
    z__4.r = a->r + nterm, z__4.i = a->i;
    z__3.r = z__4.r - 1, z__3.i = z__4.i;
    z__2.r = z__3.r * z__->r - z__3.i * z__->i, z__2.i = z__3.r * z__->i + 
	    z__3.i * z__->r;
    z__7.r = b->r + nterm, z__7.i = b->i;
    z__6.r = z__7.r - 1, z__6.i = z__7.i;
    z__5.r = nterm * z__6.r, z__5.i = nterm * z__6.i;
    z_div(&z__1, &z__2, &z__5);
    term2 = z__abs(&z__1);
    if (term2 == 0.) {
	goto L20;
    }
    if (term2 < 1.) {
	if (a->r + nterm - 1 > 1.) {
	    if (b->r + nterm - 1 > 1.) {
		if (term2 - term1 < 0.) {
		    goto L20;
		}
	    }
	}
    }
    fx += log(term2);
    if (fx > max__) {
	max__ = fx;
    }
    term1 = term2;
    goto L10;
L20:
    max__ = max__ * 2 / (bits_() * .69314718056);
    i__ = (int) (max__ * ang) + 7;
    if (i__ < 5) {
	i__ = 5;
    }
    if (*ip > i__) {
	i__ = *ip;
    }
    chgf_(&z__1, a, b, z__, &i__, lnchf);
     ret_val->r = z__1.r,  ret_val->i = z__1.i;
    return ;
}

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                   FUNCTION BITS                              * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Determines the number of significant figures  * */
 /*     *    of machine precision to arrive at the size of the array   * */
 /*     *    the numbers must must be stored in to get the accuracy    * */
 /*     *    of the solution.                                          * */
 /*     *                                                              * */
 /*     *  Subprogram called: STORE                                    * */
 /*     *                                                              * */
 /*     **************************************************************** */
 
 int HypergeometriqueConfluente::bits_(void)
{
    /* System generated locals */
    int ret_val;
    double d__1;

    /* Local variables */
    static int count;
    static double bit, bit2;

    bit = (double)1.;
    count = 0;
L10:
    ++count;
    d__1 = bit * (double)2.;
    bit2 = store_(&d__1);
    d__1 = bit2 + (double)1.;
    bit = store_(&d__1);
    if (bit - bit2 != (double)0.) {
	goto L10;
    }
    ret_val = count;
    return ret_val;
}

 double HypergeometriqueConfluente::store_(double *x)
 {
     /* System generated locals */
     double ret_val;


 /* *********************************************************** */
 /*   This function forces its argument X to be stored in a */
 /* memory location, thus providing a means of determining */
 /* floating point number characteristics (such as the machine */
 /* precision) when it is necessary to avoid computation in */
 /* high precision registers. */
 /* On input: */
 /*       X = Value to be stored. */
 /* X is not altered by this function. */
 /* On output: */
 /*       STORE = Value of X after it has been stored and */
 /*               possibly truncated or rounded to the double */
 /*               precision word length. */
 /* Modules required by STORE:  None */
 /* *********************************************************** */

     stcom_.y = *x;
     ret_val = stcom_.y;
     return ret_val;
 }

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                   FUNCTION CHGF                              * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Function that sums the Kummer series and      * */
 /*     *    returns the solution of the confluent hypergeometric      * */
 /*     *    function.                                                 * */
 /*     *                                                              * */
 /*     *  Subprograms called: ARMULT, ARYDIV, BITS, CMPADD, CMPMUL    * */
 /*     *                                                              * */
 /*     **************************************************************** */

 void HypergeometriqueConfluente::chgf_(complex *ret_val, const complex *a, const complex *b, const complex *z__, const int *l, const int *lnchf)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2;
    complex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */

    /* Local variables */
    static double rmax, numi[779], sumi[779], numr[779], sumr[779];
    static int i__;
    static complex final;
    static double ai, ci, ar, cr;
    static double xi, sigfig, xr, denomi[779], denomr[779], ai2;
    static double ci2;
    static double ar2, cr2, qi1[779], qi2[779], xi2, qr1[779], qr2[779], 
	    mx1, mx2, xr2;
    static int bit;
    static double cnt;

    bit = bits_();
    i__1 = bit / 2;
    rmax = pow_di(&c_b8, &i__1);
    i__1 = bit / 4;
    sigfig = pow_di(&c_b8, &i__1);
    ar2 = a->r * sigfig;
    ar = d_int(&ar2);
    d__1 = (ar2 - ar) * rmax;
    ar2 = d_nint(&d__1);
    ai2 = d_imag(a) * sigfig;
    ai = d_int(&ai2);
    d__1 = (ai2 - ai) * rmax;
    ai2 = d_nint(&d__1);
    cr2 = b->r * sigfig;
    cr = d_int(&cr2);
    d__1 = (cr2 - cr) * rmax;
    cr2 = d_nint(&d__1);
    ci2 = d_imag(b) * sigfig;
    ci = d_int(&ci2);
    d__1 = (ci2 - ci) * rmax;
    ci2 = d_nint(&d__1);
    xr2 = z__->r * sigfig;
    xr = d_int(&xr2);
    d__1 = (xr2 - xr) * rmax;
    xr2 = d_nint(&d__1);
    xi2 = d_imag(z__) * sigfig;
    xi = d_int(&xi2);
    d__1 = (xi2 - xi) * rmax;
    xi2 = d_nint(&d__1);
    sumr[0] = 1.;
    sumi[0] = 1.;
    numr[0] = 1.;
    numi[0] = 1.;
    denomr[0] = 1.;
    denomi[0] = 1.;
    i__1 = *l + 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	sumr[i__ + 1] = 0.;
	sumi[i__ + 1] = 0.;
	numr[i__ + 1] = 0.;
	numi[i__ + 1] = 0.;
	denomr[i__ + 1] = 0.;
	denomi[i__ + 1] = 0.;
/* L100: */
    }
    sumr[2] = 1.;
    numr[2] = 1.;
    denomr[2] = 1.;
    cnt = sigfig;
L110:
    if (sumr[2] < (double).5) {
	mx1 = sumi[*l + 2];
    } else if (sumi[2] < (double).5) {
	mx1 = sumr[*l + 2];
    } else {
/* Computing MAX */
	d__1 = sumr[*l + 2], d__2 = sumi[*l + 2];
	mx1 = dmax(d__1,d__2);
    }
    if (numr[2] < (double).5) {
	mx2 = numi[*l + 2];
    } else if (numi[2] < (double).5) {
	mx2 = numr[*l + 2];
    } else {
/* Computing MAX */
	d__1 = numr[*l + 2], d__2 = numi[*l + 2];
	mx2 = dmax(d__1,d__2);
    }
    if (mx1 - mx2 > (double)2.) {
	if (cr > 0.) {
	    z__3.r = ar, z__3.i = ai;
	    z__4.r = xr, z__4.i = xi;
	    z__2.r = z__3.r * z__4.r - z__3.i * z__4.i, z__2.i = z__3.r * 
		    z__4.i + z__3.i * z__4.r;
	    z__6.r = cr, z__6.i = ci;
	    z__5.r = cnt * z__6.r, z__5.i = cnt * z__6.i;
	    z_div(&z__1, &z__2, &z__5);
	    if (z__abs(&z__1) <= 1.) {
		goto L190;
	    }
	}
    }
    cmpmul_(sumr, sumi, &cr, &ci, qr1, qi1, l, &rmax);
    cmpmul_(sumr, sumi, &cr2, &ci2, qr2, qi2, l, &rmax);
    qr2[*l + 2] += -1;
    qi2[*l + 2] += -1;
    cmpadd_(qr1, qi1, qr2, qi2, sumr, sumi, l, &rmax);
    armult_(sumr, &cnt, sumr, l, &rmax);
    armult_(sumi, &cnt, sumi, l, &rmax);
    cmpmul_(denomr, denomi, &cr, &ci, qr1, qi1, l, &rmax);
    cmpmul_(denomr, denomi, &cr2, &ci2, qr2, qi2, l, &rmax);
    qr2[*l + 2] += -1;
    qi2[*l + 2] += -1;
    cmpadd_(qr1, qi1, qr2, qi2, denomr, denomi, l, &rmax);
    armult_(denomr, &cnt, denomr, l, &rmax);
    armult_(denomi, &cnt, denomi, l, &rmax);
    cmpmul_(numr, numi, &ar, &ai, qr1, qi1, l, &rmax);
    cmpmul_(numr, numi, &ar2, &ai2, qr2, qi2, l, &rmax);
    qr2[*l + 2] += -1;
    qi2[*l + 2] += -1;
    cmpadd_(qr1, qi1, qr2, qi2, numr, numi, l, &rmax);
    cmpmul_(numr, numi, &xr, &xi, qr1, qi1, l, &rmax);
    cmpmul_(numr, numi, &xr2, &xi2, qr2, qi2, l, &rmax);
    qr2[*l + 2] += -1;
    qi2[*l + 2] += -1;
    cmpadd_(qr1, qi1, qr2, qi2, numr, numi, l, &rmax);
    cmpadd_(sumr, sumi, numr, numi, sumr, sumi, l, &rmax);
    cnt += sigfig;
    ar += sigfig;
    cr += sigfig;
    goto L110;
L190:
    arydiv_(sumr, sumi, denomr, denomi, &final, l, lnchf, &rmax, &bit);
     ret_val->r = final.r,  ret_val->i = final.i;
    return ;
} 

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE ARADD                             * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Accepts two arrays of numbers and returns     * */
 /*     *    the sum of the array.  Each array is holding the value    * */
 /*     *    of one number in the series.  The parameter L is the      * */
 /*     *    size of the array representing the number and RMAX is     * */
 /*     *    the actual number of digits needed to give the numbers    * */
 /*     *    the desired accuracy.                                     * */
 /*     *                                                              * */
 /*     *  Subprograms called: none                                    * */
 /*     *                                                              * */
 /*     **************************************************************** */
 
 int HypergeometriqueConfluente::aradd_(const double *a, const double *b, 
        double *c__, const int *l, const double *rmax)
{
    /* System generated locals */
    int i__1;
    double d__1;

    /* Builtin functions */

    /* Local variables */
    static int i__, j, ediff;
    static double z__[779];

    /* Parameter adjustments */
    ++c__;
    ++b;
    ++a;

    /* Function Body */
    i__1 = *l + 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	z__[i__ + 1] = 0.;
/* L110: */
    }
    d__1 = a[*l + 1] - b[*l + 1];
    ediff = (int) d_nint(&d__1);
    if (abs(a[1]) < (double).5 || ediff <= -(*l)) {
	goto L111;
    }
    if (abs(b[1]) < (double).5 || ediff >= *l) {
	goto L113;
    }
    goto L115;
L111:
    i__1 = *l + 1;
    for (i__ = -1; i__ <= i__1; ++i__) {
	c__[i__] = b[i__];
/* L112: */
    }
    goto L311;
L113:
    i__1 = *l + 1;
    for (i__ = -1; i__ <= i__1; ++i__) {
	c__[i__] = a[i__];
/* L114: */
    }
    goto L311;
L115:
    z__[0] = a[-1];
    if ((d__1 = a[-1] - b[-1], abs(d__1)) < (double).5) {
	goto L200;
    }
    if (ediff > 0) {
	z__[*l + 2] = a[*l + 1];
	goto L233;
    }
    if (ediff < 0) {
	z__[*l + 2] = b[*l + 1];
	z__[0] = b[-1];
	goto L266;
    }
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (a[i__] > b[i__]) {
	    z__[*l + 2] = a[*l + 1];
	    goto L233;
	}
	if (a[i__] < b[i__]) {
	    z__[*l + 2] = b[*l + 1];
	    z__[0] = b[-1];
	    goto L266;
	}
/* L120: */
    }
    goto L300;
L200:
    if (ediff > 0) {
	goto L203;
    }
    if (ediff < 0) {
	goto L207;
    }
    z__[*l + 2] = a[*l + 1];
    for (i__ = *l; i__ >= 1; --i__) {
	z__[i__ + 1] = a[i__] + b[i__] + z__[i__ + 1];
	if (z__[i__ + 1] >= *rmax) {
	    z__[i__ + 1] -= *rmax;
	    z__[i__] = 1.;
	}
/* L201: */
    }
    if (z__[1] > (double).5) {
	for (i__ = *l; i__ >= 1; --i__) {
	    z__[i__ + 1] = z__[i__];
/* L202: */
	}
	z__[*l + 2] += 1.;
	z__[1] = 0.;
    }
    goto L300;
L203:
    z__[*l + 2] = a[*l + 1];
    i__1 = ediff + 1;
    for (i__ = *l; i__ >= i__1; --i__) {
	z__[i__ + 1] = a[i__] + b[i__ - ediff] + z__[i__ + 1];
	if (z__[i__ + 1] >= *rmax) {
	    z__[i__ + 1] -= *rmax;
	    z__[i__] = 1.;
	}
/* L204: */
    }
    for (i__ = ediff; i__ >= 1; --i__) {
	z__[i__ + 1] = a[i__] + z__[i__ + 1];
	if (z__[i__ + 1] >= *rmax) {
	    z__[i__ + 1] -= *rmax;
	    z__[i__] = 1.;
	}
/* L205: */
    }
    if (z__[1] > (double).5) {
	for (i__ = *l; i__ >= 1; --i__) {
	    z__[i__ + 1] = z__[i__];
/* L206: */
	}
	z__[*l + 2] += 1;
	z__[1] = 0.;
    }
    goto L300;
L207:
    z__[*l + 2] = b[*l + 1];
    i__1 = 1 - ediff;
    for (i__ = *l; i__ >= i__1; --i__) {
	z__[i__ + 1] = a[i__ + ediff] + b[i__] + z__[i__ + 1];
	if (z__[i__ + 1] >= *rmax) {
	    z__[i__ + 1] -= *rmax;
	    z__[i__] = 1.;
	}
/* L208: */
    }
    for (i__ = -ediff; i__ >= 1; --i__) {
	z__[i__ + 1] = b[i__] + z__[i__ + 1];
	if (z__[i__ + 1] >= *rmax) {
	    z__[i__ + 1] -= *rmax;
	    z__[i__] = 1.;
	}
/* L209: */
    }
    if (z__[1] > (double).5) {
	for (i__ = *l; i__ >= 1; --i__) {
	    z__[i__ + 1] = z__[i__];
/* L210: */
	}
	z__[*l + 2] += 1.;
	z__[1] = 0.;
    }
    goto L300;
L233:
    if (ediff > 0) {
	goto L243;
    }
    for (i__ = *l; i__ >= 1; --i__) {
	z__[i__ + 1] = a[i__] - b[i__] + z__[i__ + 1];
	if (z__[i__ + 1] < 0.) {
	    z__[i__ + 1] += *rmax;
	    z__[i__] = -1.;
	}
/* L234: */
    }
    goto L290;
L243:
    i__1 = ediff + 1;
    for (i__ = *l; i__ >= i__1; --i__) {
	z__[i__ + 1] = a[i__] - b[i__ - ediff] + z__[i__ + 1];
	if (z__[i__ + 1] < 0.) {
	    z__[i__ + 1] += *rmax;
	    z__[i__] = -1.;
	}
/* L244: */
    }
    for (i__ = ediff; i__ >= 1; --i__) {
	z__[i__ + 1] = a[i__] + z__[i__ + 1];
	if (z__[i__ + 1] < 0.) {
	    z__[i__ + 1] += *rmax;
	    z__[i__] = -1.;
	}
/* L245: */
    }
    goto L290;
L266:
    if (ediff < 0) {
	goto L276;
    }
    for (i__ = *l; i__ >= 1; --i__) {
	z__[i__ + 1] = b[i__] - a[i__] + z__[i__ + 1];
	if (z__[i__ + 1] < 0.) {
	    z__[i__ + 1] += *rmax;
	    z__[i__] = -1.;
	}
/* L267: */
    }
    goto L290;
L276:
    i__1 = 1 - ediff;
    for (i__ = *l; i__ >= i__1; --i__) {
	z__[i__ + 1] = b[i__] - a[i__ + ediff] + z__[i__ + 1];
	if (z__[i__ + 1] < 0.) {
	    z__[i__ + 1] += *rmax;
	    z__[i__] = -1.;
	}
/* L277: */
    }
    for (i__ = -ediff; i__ >= 1; --i__) {
	z__[i__ + 1] = b[i__] + z__[i__ + 1];
	if (z__[i__ + 1] < 0.) {
	    z__[i__ + 1] += *rmax;
	    z__[i__] = -1.;
	}
/* L278: */
    }
L290:
    if (z__[2] > (double).5) {
	goto L300;
    }
    i__ = 1;
L291:
    ++i__;
    if (z__[i__ + 1] < (double).5 && i__ < *l + 1) {
	goto L291;
    }
    if (i__ == *l + 1) {
	z__[0] = 1.;
	z__[*l + 2] = 0.;
	goto L300;
    }
/* L292: */
    i__1 = *l + 1 - i__;
    for (j = 1; j <= i__1; ++j) {
	z__[j + 1] = z__[j + i__];
/* L293: */
    }
    i__1 = *l;
    for (j = *l + 2 - i__; j <= i__1; ++j) {
	z__[j + 1] = 0.;
/* L294: */
    }
    z__[*l + 2] = z__[*l + 2] - i__ + 1;
L300:
    i__1 = *l + 1;
    for (i__ = -1; i__ <= i__1; ++i__) {
	c__[i__] = z__[i__ + 1];
/* L310: */
    }
L311:
    if (c__[1] < (double).5) {
	c__[-1] = 1.;
	c__[*l + 1] = 0.;
    }
    return 0;
}

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE ARSUB                             * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Accepts two arrays and subtracts each element * */
 /*     *    in the second array from the element in the first array   * */
 /*     *    and returns the solution.  The parameters L and RMAX are  * */
 /*     *    the size of the array and the number of digits needed for * */
 /*     *    the accuracy, respectively.                               * */
 /*     *                                                              * */
 /*     *  Subprograms called: ARADD                                   * */
 /*     *                                                              * */
 /*     **************************************************************** */
 
 int HypergeometriqueConfluente::arsub_(const double *a, const double *b, 
        double *c__, const int *l, const double *rmax)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;
    static double b2[779];

    /* Parameter adjustments */
    ++c__;
    ++b;
    ++a;

    /* Function Body */
    i__1 = *l + 1;
    for (i__ = -1; i__ <= i__1; ++i__) {
	b2[i__ + 1] = b[i__];
/* L100: */
    }
    b2[0] *= -1.;
    aradd_(&a[-1], b2, &c__[-1], l, rmax);
    return 0;
}

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE ARMULT                            * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Accepts two arrays and returns the product.   * */
 /*     *    L and RMAX are the size of the arrays and the number of   * */
 /*     *    digits needed to represent the numbers with the required  * */
 /*     *    accuracy.                                                 * */
 /*     *                                                              * */
 /*     *  Subprograms called: none                                    * */
 /*     *                                                              * */
 /*     **************************************************************** */

 int HypergeometriqueConfluente::armult_(const double *a, const double *b, double *c__, const int *l, const double *rmax)
{
    /* System generated locals */
    int i__1;
    double d__1;

    /* Builtin functions */

    /* Local variables */
    static double rmax2;
    static int i__;
    static double z__[779], carry, b2;

    /* Parameter adjustments */
    ++c__;
    ++a;

    /* Function Body */
    rmax2 = 1. / *rmax;
    z__[0] = d_sign(&c_b53, b) * a[-1];
    b2 = abs(*b);
    z__[*l + 2] = a[*l + 1];
    i__1 = *l;
    for (i__ = 0; i__ <= i__1; ++i__) {
	z__[i__ + 1] = 0.;
/* L100: */
    }
    if (b2 <= 1e-10 || a[1] <= 1e-10) {
	z__[0] = 1.;
	z__[*l + 2] = 0.;
	goto L198;
    }
    for (i__ = *l; i__ >= 1; --i__) {
	z__[i__ + 1] = a[i__] * b2 + z__[i__ + 1];
	if (z__[i__ + 1] >= *rmax) {
	    d__1 = z__[i__ + 1] / *rmax;
	    carry = d_int(&d__1);
	    z__[i__ + 1] -= carry * *rmax;
	    z__[i__] = carry;
	}
/* L110: */
    }
    if (z__[1] < (double).5) {
	goto L150;
    }
    for (i__ = *l; i__ >= 1; --i__) {
	z__[i__ + 1] = z__[i__];
/* L120: */
    }
    z__[*l + 2] += 1.;
    z__[1] = 0.;
L150:
L198:
    i__1 = *l + 1;
    for (i__ = -1; i__ <= i__1; ++i__) {
	c__[i__] = z__[i__ + 1];
/* L199: */
    }
    if (c__[1] < (double).5) {
	c__[-1] = 1.;
	c__[*l + 1] = 0.;
    }
    return 0;
}

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE CMPADD                            * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Takes two arrays representing one real and    * */
 /*     *    one imaginary part, and adds two arrays representing      * */
 /*     *    another complex number and returns two array holding the  * */
 /*     *    complex sum.                                              * */
 /*     *              (CR,CI) = (AR+BR, AI+BI)                        * */
 /*     *                                                              * */
 /*     *  Subprograms called: ARADD                                   * */
 /*     *                                                              * */
 /*     **************************************************************** */
 
 int HypergeometriqueConfluente::cmpadd_(const double *ar, const double *ai, const double *br, const double *bi, double *cr, double *ci, const int *l, const double *rmax)
{

    /* Parameter adjustments */
    ++ci;
    ++cr;
    ++bi;
    ++br;
    ++ai;
    ++ar;

    /* Function Body */
    aradd_(&ar[-1], &br[-1], &cr[-1], l, rmax);
    aradd_(&ai[-1], &bi[-1], &ci[-1], l, rmax);
    return 0;
}

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE CMPSUB                            * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Takes two arrays representing one real and    * */
 /*     *    one imaginary part, and subtracts two arrays representing * */
 /*     *    another complex number and returns two array holding the  * */
 /*     *    complex sum.                                              * */
 /*     *              (CR,CI) = (AR+BR, AI+BI)                        * */
 /*     *                                                              * */
 /*     *  Subprograms called: ARADD                                   * */
 /*     *                                                              * */
 /*     **************************************************************** */
 int HypergeometriqueConfluente::cmpsub_(double *ar, double *ai, double *br, double *bi, double *cr, double *ci, int *l, double *rmax)
{
    /* Parameter adjustments */
    ++ci;
    ++cr;
    ++bi;
    ++br;
    ++ai;
    ++ar;

    /* Function Body */
    arsub_(&ar[-1], &br[-1], &cr[-1], l, rmax);
    arsub_(&ai[-1], &bi[-1], &ci[-1], l, rmax);
    return 0;
} 

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE CMPMUL                            * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Takes two arrays representing one real and    * */
 /*     *    one imaginary part, and multiplies it with two arrays     * */
 /*     *    representing another complex number and returns the       * */
 /*     *    complex product.                                          * */
 /*     *                                                              * */
 /*     *  Subprograms called: ARMULT, ARSUB, ARADD                    * */
 /*     *                                                              * */
 /*     **************************************************************** */
 int HypergeometriqueConfluente::cmpmul_(const double *ar, const double *ai, const double *br, const double *bi, double *cr, double *ci, const int *l, const double *rmax)
{
    static double d1[779], d2[779];

    /* Parameter adjustments */
    ++ci;
    ++cr;
    ++ai;
    ++ar;

    /* Function Body */
    armult_(&ar[-1], br, d1, l, rmax);
    armult_(&ai[-1], bi, d2, l, rmax);
    arsub_(d1, d2, &cr[-1], l, rmax);
    armult_(&ar[-1], bi, d1, l, rmax);
    armult_(&ai[-1], br, d2, l, rmax);
    aradd_(d1, d2, &ci[-1], l, rmax);
    return 0;
}

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE ARYDIV                            * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Returns the double precision complex number   * */
 /*     *    resulting from the division of four arrays, representing  * */
 /*     *    two complex numbers.  The number returned will be in one  * */
 /*     *    two different forms.  Either standard scientific or as    * */
 /*     *    the log of the number.                                    * */
 /*     *                                                              * */
 /*     *  Subprograms called: CONV21, CONV12, EADD, ECPDIV, EMULT     * */
 /*     *                                                              * */
 /*     **************************************************************** */
 
 int HypergeometriqueConfluente::arydiv_(const double *ar, const double *ai, const double *br, const double *bi, complex *c__, const int *l, const int *lnchf, const double *rmax, const int *bit)
{
    /* System generated locals */
    double d__1;
    complex z__1;

    /* Builtin functions */

    /* Local variables */
    static int rexp;
    static double x;
    static double e1, e2, e3;
    static double n1, n2, n3, x1, x2, ae[4]	/* was [2][2] */, be[4]	/* 
	    was [2][2] */, ce[4]	/* was [2][2] */;
    static int ii10, ir10;
    static double ri10, phi, rr10, dum1, dum2;
    /* Parameter adjustments */
    ++bi;
    ++br;
    ++ai;
    ++ar;

    /* Function Body */
    rexp = *bit / 2;
    x = rexp * (ar[*l + 1] - 2);
    rr10 = x * d_lg10(&c_b8) / d_lg10(&c_b65);
    ir10 = (int) rr10;
    rr10 -= ir10;
    x = rexp * (ai[*l + 1] - 2);
    ri10 = x * d_lg10(&c_b8) / d_lg10(&c_b65);
    ii10 = (int) ri10;
    ri10 -= ii10;
    d__1 = ar[1] * *rmax * *rmax + ar[2] * *rmax + ar[3];
    dum1 = d_sign(&d__1, &ar[-1]);
    d__1 = ai[1] * *rmax * *rmax + ai[2] * *rmax + ai[3];
    dum2 = d_sign(&d__1, &ai[-1]);
    dum1 *= pow_dd(&c_b65, &rr10);
    dum2 *= pow_dd(&c_b65, &ri10);
    z__1.r = dum1, z__1.i = dum2;
    conv12_(&z__1, ae);
    ae[2] += ir10;
    ae[3] += ii10;
    x = rexp * (br[*l + 1] - 2);
    rr10 = x * d_lg10(&c_b8) / d_lg10(&c_b65);
    ir10 = (int) rr10;
    rr10 -= ir10;
    x = rexp * (bi[*l + 1] - 2);
    ri10 = x * d_lg10(&c_b8) / d_lg10(&c_b65);
    ii10 = (int) ri10;
    ri10 -= ii10;
    d__1 = br[1] * *rmax * *rmax + br[2] * *rmax + br[3];
    dum1 = d_sign(&d__1, &br[-1]);
    d__1 = bi[1] * *rmax * *rmax + bi[2] * *rmax + bi[3];
    dum2 = d_sign(&d__1, &bi[-1]);
    dum1 *= pow_dd(&c_b65, &rr10);
    dum2 *= pow_dd(&c_b65, &ri10);
    z__1.r = dum1, z__1.i = dum2;
    conv12_(&z__1, be);
    be[2] += ir10;
    be[3] += ii10;
    ecpdiv_(ae, be, ce);
    if (*lnchf == 0) {
	conv21_(ce, c__);
    } else {
	emult_(ce, &ce[2], ce, &ce[2], &n1, &e1);
	emult_(&ce[1], &ce[3], &ce[1], &ce[3], &n2, &e2);
	eadd_(&n1, &e1, &n2, &e2, &n3, &e3);
	n1 = ce[0];
	e1 = ce[2] - ce[3];
	x2 = ce[1];
	if (e1 > 74.) {
	    x1 = 1e75;
	} else if (e1 < -74.) {
	    x1 = 0.;
	} else {
	    x1 = n1 * pow_dd(&c_b65, &e1);
	}
	phi = atan2(x2, x1);
	d__1 = (log(n3) + e3 * log(10.)) * .5;
	z__1.r = d__1, z__1.i = phi;
	c__->r = z__1.r, c__->i = z__1.i;
    }
    return 0;
} 
 
 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE EMULT                             * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Takes one base and exponent and multiplies it * */
 /*     *    by another numbers base and exponent to give the product  * */
 /*     *    in the form of base and exponent.                         * */
 /*     *                                                              * */
 /*     *  Subprograms called: none                                    * */
 /*     *                                                              * */
 /*     **************************************************************** */
 int HypergeometriqueConfluente::emult_(double *n1, double *e1, double *n2, double *e2, double *nf, double *ef)
{
    *nf = *n1 * *n2;
    *ef = *e1 + *e2;
    if (abs(*nf) >= 10.) {
	*nf /= 10.;
	*ef += 1.;
    }
    return 0;
}

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE EDIV                              * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : returns the solution in the form of base and  * */
 /*     *    exponent of the division of two exponential numbers.      * */
 /*     *                                                              * */
 /*     *  Subprograms called: none                                    * */
 /*     *                                                              * */
 /*     **************************************************************** */
 int HypergeometriqueConfluente::ediv_(double *n1, double *e1, double *n2, double *e2, double *nf, double *ef)
 {
     *nf = *n1 / *n2;
     *ef = *e1 - *e2;
     if (abs(*nf) < 1. && *nf != 0.) {
         *nf *= 10.;
         *ef += -1.;
     }
     return 0;
 }

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE EADD                              * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Returns the sum of two numbers in the form    * */
 /*     *    of a base and an exponent.                                * */
 /*     *                                                              * */
 /*     *  Subprograms called: none                                    * */
 /*     *                                                              * */
 /*     **************************************************************** */
 
 int HypergeometriqueConfluente::eadd_(double *n1, double *e1, double *n2, double *e2, double *nf, double *ef)
{
    /* Builtin functions */

    /* Local variables */
    static double ediff;

    ediff = *e1 - *e2;
    if (ediff > 36.) {
	*nf = *n1;
	*ef = *e1;
    } else if (ediff < -36.) {
	*nf = *n2;
	*ef = *e2;
    } else {
	*nf = *n1 * pow_dd(&c_b65, &ediff) + *n2;
	*ef = *e2;
L400:
	if (abs(*nf) < 10.) {
	    goto L410;
	}
	*nf /= 10.;
	*ef += 1.;
	goto L400;
L410:
	if (abs(*nf) >= 1. || *nf == 0.) {
	    goto L420;
	}
	*nf *= 10.;
	*ef += -1.;
	goto L410;
    }
L420:
    return 0;
}

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE ESUB                              * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Returns the solution to the subtraction of    * */
 /*     *    two numbers in the form of base and exponent.             * */
 /*     *                                                              * */
 /*     *  Subprograms called: EADD                                    * */
 /*     *                                                              * */
 /*     **************************************************************** */
 
 int HypergeometriqueConfluente::esub_(double *n1, double *e1, double *n2, double *e2, double *nf, double *ef)
{
    /* System generated locals */
    double d__1;

    /* Local variables */

    d__1 = *n2 * -1.;
    eadd_(n1, e1, &d__1, e2, nf, ef);
    return 0;
}

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE CONV12                            * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Converts a number from complex notation to a  * */
 /*     *    form of a 2x2 real array.                                 * */
 /*     *                                                              * */
 /*     *  Subprograms called: none                                    * */
 /*     *                                                              * */
 /*     **************************************************************** */
 int HypergeometriqueConfluente::conv12_(complex *cn, double *cae)
{
    /* Builtin functions */

    /* Parameter adjustments */
    cae -= 3;

    /* Function Body */
    cae[3] = cn->r;
    cae[5] = 0.;
L300:
    if (abs(cae[3]) < 10.) {
	goto L310;
    }
    cae[3] /= 10.;
    cae[5] += 1.;
    goto L300;
L310:
    if (abs(cae[3]) >= 1. || cae[3] == 0.) {
	goto L320;
    }
    cae[3] *= 10.;
    cae[5] += -1.;
    goto L310;
L320:
    cae[4] = d_imag(cn);
    cae[6] = 0.;
L330:
    if (abs(cae[4]) < 10.) {
	goto L340;
    }
    cae[4] /= 10.;
    cae[6] += 1.;
    goto L330;
L340:
    if (abs(cae[4]) >= 1. || cae[4] == 0.) {
	goto L350;
    }
    cae[4] *= 10.;
    cae[6] += -1.;
    goto L340;
L350:
    return 0;
}

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE CONV21                            * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Converts a number represented in a 2x2 real   * */
 /*     *    array to the form of a complex number.                    * */
 /*     *                                                              * */
 /*     *  Subprograms called: none                                    * */
 /*     *                                                              * */
 /*     **************************************************************** */
 
 int HypergeometriqueConfluente::conv21_(double *cae, complex *cn)
{
    /* System generated locals */
    double d__1, d__2;
    complex z__1;

    /* Builtin functions */

    /* Parameter adjustments */
    cae -= 3;

    /* Function Body */
    if (cae[5] > 75. || cae[6] > 75.) {
	cn->r = 1e75, cn->i = 1e75;
    } else if (cae[6] < -75.) {
	d__1 = cae[3] * pow_dd(&c_b65, &cae[5]);
	z__1.r = d__1, z__1.i = 0.;
	cn->r = z__1.r, cn->i = z__1.i;
    } else {
	d__1 = cae[3] * pow_dd(&c_b65, &cae[5]);
	d__2 = cae[4] * pow_dd(&c_b65, &cae[6]);
	z__1.r = d__1, z__1.i = d__2;
	cn->r = z__1.r, cn->i = z__1.i;
    }
    return 0;
}

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE ECPMUL                            * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Multiplies two numbers which are each         * */
 /*     *    represented in the form of a two by two array and returns * */
 /*     *    the solution in the same form.                            * */
 /*     *                                                              * */
 /*     *  Subprograms called: EMULT, ESUB, EADD                       * */
 /*     *                                                              * */
 /*     **************************************************************** */

 int HypergeometriqueConfluente::ecpmul_(double *a, double *b, double *c__)
{
    static double c2[4]	/* was [2][2] */, e1, e2;
    static double n1, n2;
    /* Parameter adjustments */
    c__ -= 3;
    b -= 3;
    a -= 3;

    /* Function Body */
    emult_(&a[3], &a[5], &b[3], &b[5], &n1, &e1);
    emult_(&a[4], &a[6], &b[4], &b[6], &n2, &e2);
    esub_(&n1, &e1, &n2, &e2, c2, &c2[2]);
    emult_(&a[3], &a[5], &b[4], &b[6], &n1, &e1);
    emult_(&a[4], &a[6], &b[3], &b[5], &n2, &e2);
    eadd_(&n1, &e1, &n2, &e2, &c__[4], &c__[6]);
    c__[3] = c2[0];
    c__[5] = c2[2];
    return 0;
} 

 /*     **************************************************************** */
 /*     *                                                              * */
 /*     *                 SUBROUTINE ECPDIV                            * */
 /*     *                                                              * */
 /*     *                                                              * */
 /*     *  Description : Divides two numbers and returns the solution. * */
 /*     *    All numbers are represented by a 2x2 array.               * */
 /*     *                                                              * */
 /*     *  Subprograms called: EADD, ECPMUL, EDIV, EMULT               * */
 /*     *                                                              * */
 /*     **************************************************************** */
 
 int HypergeometriqueConfluente::ecpdiv_(double *a, double *b, double *c__)
{

	static double b2[4]	/* was [2][2] */, c2[4]	/* was [2][2] */, e1, 
	    e2, e3;
    static double n1, n2, n3;
    /* Parameter adjustments */
    c__ -= 3;
    b -= 3;
    a -= 3;

    /* Function Body */
    b2[0] = b[3];
    b2[2] = b[5];
    b2[1] = b[4] * -1.;
    b2[3] = b[6];
    ecpmul_(&a[3], b2, c2);
    emult_(&b[3], &b[5], &b[3], &b[5], &n1, &e1);
    emult_(&b[4], &b[6], &b[4], &b[6], &n2, &e2);
    eadd_(&n1, &e1, &n2, &e2, &n3, &e3);
    ediv_(c2, &c2[2], &n3, &e3, &c__[3], &c__[5]);
    ediv_(&c2[1], &c2[3], &n3, &e3, &c__[4], &c__[6]);
    return 0;
}

double HypergeometriqueConfluente::gamma(double x)
{
    int i,k,m;
    double ga,gr,r,z;

    static double g[] = {
        1.0,
        0.5772156649015329,
       -0.6558780715202538,
       -0.420026350340952e-1,
        0.1665386113822915,
       -0.421977345555443e-1,
       -0.9621971527877e-2,
        0.7218943246663e-2,
       -0.11651675918591e-2,
       -0.2152416741149e-3,
        0.1280502823882e-3,
       -0.201348547807e-4,
       -0.12504934821e-5,
        0.1133027232e-5,
       -0.2056338417e-6,
        0.6116095e-8,
        0.50020075e-8,
       -0.11812746e-8,
        0.1043427e-9,
        0.77823e-11,
       -0.36968e-11,
        0.51e-12,
       -0.206e-13,
       -0.54e-14,
        0.14e-14};

    if (x > 171.0) return 1e308;    // This value is an overflow flag.
    if (x == (int)x) {
        if (x > 0.0) {
            ga = 1.0;               // use factorial
            for (i=2;i<x;i++) {
               ga *= i;
            }
         }
         else
            ga = 1e308;
     }
     else {
        if (fabs(x) > 1.0) {
            z = fabs(x);
            m = (int)z;
            r = 1.0;
            for (k=1;k<=m;k++)
			{
                r *= (z-k);
            }
            z -= m;
        }
        else
            z = x;
        gr = g[24];
        for (k=23;k>=0;k--) {
            gr = gr*z+g[k];
        }
        ga = 1.0/(gr*z);
        if (fabs(x) > 1.0) {
            ga *= r;
            if (x < 0.0) {
                ga = -PI/(x*ga*sin(PI*x));
            }
        }
    }
    return ga;
}

double HypergeometriqueConfluente::CHGM(double a, double b, double x)
{


 double a0,a1;
 double hg,hg1,hg2;
 double r,r1,r2,rg;
 double sum1,sum2;
 double ta,tb,tba;
 double x0,xg;
 double y0(0.0),y1(0.0);
 int i,j,k,la,m,n,nl;

 a0 = a;
 a1 = a;
 x0 = x;
 hg = 0.0;

 if (b == 0.0 || b == -abs(floor(b)))
	 hg = DBL_MAX;
 else if (a==0.0 || x==0.0)
		hg = 1.0;
 else if (a==-1.0)
		hg = (1.0 + x/b)*exp(x);
 else if (a-b==1.0)
		hg = exp(x);
 else if (a==1.0 && b==2.0)
		hg = (exp(x) - 1.0)/x;
 else if (a == floor(a) && a < 0.0)
		{
			m = (int)floor(-a);
			r=1.0;
			for (k=1;k<m+1;++k)
			{
				r *= (a+(double)k-1.0)/(double)k/(b+(double)k-1.0)*x;
				hg+=r;
			}
		}
 if (hg!=0.)
	 return hg;

 if (x<0.)
 {
	a = b-a;
	a0 = a;
	x =abs(x);
 }
 if (a<2.)
	 nl=0;
 else
 {
	nl = 1;
	la = (int)floor(a);
	a-=(double)(la+1);
 }

 for (n=0;n<=nl;++n)
 {
	if (2.0<=a0)
		++a;
	if (x<= (30. + abs(b)) || a<0.)
	{
		hg = 1.;
		rg = 1.;
		j=1;
		while (j<=500)
		{
			rg *= (((a + (double)j - 1.0)/(((double)j)*(b + (double)j - 1)))*x);
			hg += rg;
			++j;
			if ( abs(rg/hg) < 1.0E-15 )
				break;
		}
	}
	else
	{
		ta = gamma(a);
		tb = gamma(b);
		xg = b - a;
		tba = gamma (xg);
		sum1 = 1.;
		sum2 = 1.;
		r1 = 1.;
		r2 = 1.;
		for (i=0;i<8;++i)
		{
	        r1 = - r1*(a + (double)i)*( a - b + ((double)i+1.0))/(x*((double)i+1.0));
		    r2 = - r2*(b - a + (double)i)*( a - ((double)i+1.0))/(x*((double)i+1.0));
			sum1 += r1;
			sum2 += r2;
		}
      hg1 = (tb/tba)*pow(x,-a)*cos(PI*a)*sum1;
      hg2 = (tb/ta)*exp(x)*pow(x,a-b)*sum2;
      hg = hg1 + hg2;
	}
	if (n==0)
		y0 = hg;
	if (n==1)
		y1 = hg;
 }
 hg = y0;
 if ( 2. <= a0 )
 {
	for (i = 0;i<la-1;++i)
	{
	  hg = ((2.*a - b + x)*y1 + (b - a)*y0)/a;
      y0 = y1;
      y1 = hg;
	  ++a;
	}
 }
 if ( x0 < 0.)
	 hg *= exp(x0);

 return hg;
}

double HypergeometriqueConfluente::taylor1f1(double a,double b, double z, double tol, double erreur)
{

	const int seuil = 1000;
	double a1_j = 1.0;
	double a1_jplus1 = 0.0;
	double b1 =1.0;
	double cpt = 1.0;
	double result = 0.0;

	for (int j=0;j<seuil;++j)
	{
	
		a1_jplus1 = (((a+cpt-1)/(b+cpt-1))*z/cpt)*a1_j;
		b1+=a1_jplus1;
		if (((abs(a1_j)/abs(b1))<tol) && ((abs(a1_jplus1)/abs(b1))<tol))
			break;
		a1_j = a1_jplus1;
		++cpt;
	}
	result = b1;
	if (cpt==(double)seuil)
	{
		result =-1.;
		erreur = -1;
	}
	return result;
}


}