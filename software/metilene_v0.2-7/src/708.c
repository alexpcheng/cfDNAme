/* 708.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

//#include "f2c.h"
#include "limits.h"
#include "basic-types.h"
#include "mathematics.h"
#include "float.h"
/* Table of constant values */

static long int c__1 = 1;
static long int c__0 = 0;
static double c_b188 = 1.f;



double r_sign(double *a, double *b)
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}


/*      ALGORITHM 708, COLLECTED ALGORITHMS FROM ACM. */
/*      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*      VOL. 18, NO. 3, SEPTEMBER, 1992, PP. 360-373z. */
/*     PROGRAM BTST (OUTPUT, TAPE6=OUTPUT) */
/* ----------------------------------------------------------------------- */

/*     SAMPLE PROGRAM USING BRATIO. GIVEN THE NONNEGATIVE VALUES */
/*     A, B, X, Y WHERE A AND B ARE NOT BOTH 0 AND X + Y = 1. THEN */

/*              CALL BRATIO (A, B, X, Y, W, W1, IERR) */

/*     COMPUTES THE VALUES */

/*                W = I (A,B)  AND W1 = 1 - I (A,B). */
/*                     X                     X */

/*     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS. */
/*     IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND */
/*     W AND W1 ARE COMPUTED. FOR MORE DETAILS SEE THE IN-LINE */
/*     DOCUMENTATION OF BRATIO. */

/*     THE LAST FUNCTION IN THIS PACKAGE, IPMPAR, MUST BE DEFINED */
/*     FOR THE PARTICULAR COMPUTER BEING USED. FOR DETAILS SEE THE */
/*     IN-LINE DOCUMENTATION OF IPMPAR. */

/*     NO DATA IS READ. THE OUTPUT FOR THE PROGRAM IS WRITTEN ON */
/*     UNIT 6. THE FIRST STATMENT OF THIS TEXT MAY BE USED TO */
/*     BEGIN THE PROGRAM FOR THE CDC 6000-7000 SERIES COMPUTERS. */
/* ----------------------------------------------------------------------- */
/* Main program */ 
/*
int main(void) {
   int l, ierr;
   double a = 5.3f;
   double   b = 10.1f; 
   double x = .01f, y, w, w1;
   
   for (l = 1; l <= 50; ++l) {
	y = .5f - x + .5f;
	bratio_(&a, &b, &x, &y, &w, &w1, &ierr);
	if (ierr != 0) {
      fprintf(stderr, "stop: %d!\n", ierr);
      exit(-1);
	}

    fprintf(stderr, "I(a,b,x,y) = I(%f,%f,%f,%f) = %.16f, %.16f\n", a, b, x, y, w, w1);
	x += .01f;
  }
}
*/
//{
//    /* Format strings */
//    static char fmt_1[] = "(\0021   X     Y\002,11x,\002W\002,14x,\002W1\002"
//	    "/)";
//    static char fmt_2[] = "(2f6.2,2e16.6)";

//    /* Builtin functions */
//    long int s_wsfe(cilist *), e_wsfe(void);
//    /* Subroutine */ int s_stop(char *, ftnlen);
//    long int do_fio(long int *, char *, ftnlen);

    /* Local variables */
//    static double a, b;
//    static long int l;
//    static double w, x, y, w1;
//    static long int ierr;
//    extern /* Subroutine */ int bratio_(double *, double *, double *, double *, double *
//	    , double *, long int *);

//    /* Fortran I/O blocks */
//    static cilist io___1 = { 0, 6, 0, fmt_1, 0 };
//    static cilist io___10 = { 0, 6, 0, fmt_2, 0 };


//    s_wsfe(&io___1);
//    e_wsfe();
/* L2: */

//    a = 5.3f;
//    b = 10.1f;
//    x = .01f;
//    for (l = 1; l <= 50; ++l) {
//	y = .5f - x + .5f;
//	bratio_(&a, &b, &x, &y, &w, &w1, &ierr);
//	if (ierr != 0) {
//	    s_stop("", (ftnlen)0);
//	}
//	s_wsfe(&io___10);
//	do_fio(&c__1, (char *)&x, (ftnlen)sizeof(double));
//	do_fio(&c__1, (char *)&y, (ftnlen)sizeof(double));
//	do_fio(&c__1, (char *)&w, (ftnlen)sizeof(double));
//	do_fio(&c__1, (char *)&w1, (ftnlen)sizeof(double));
//	e_wsfe();
//	x += .01f;
/* L10: */
 //   }
//    s_stop("", (ftnlen)0);
//    return 0;
//} /* MAIN__ */

/* Subroutine */ int bratio_(double *a, double *b, double *x, double *y, double *w, 
	double *w1, long int *ierr)
{
    /* System generated locals */
    double r__1, r__2;
    double d__1, d__2;

    /* Builtin functions */
    double pow(double, double);

    /* Local variables */
    static long int n;
    static double t, z__, a0, b0, x0, y0;
    static long int ind;
    extern double bup_(double *, double *, double *, double *, long int *, double *);
    static double eps;
    static long int ierr1;
    extern double bfrac_(double *, double *, double *, double *, double *, double *);
    extern /* Subroutine */ int bgrat_(double *, double *, double *, double *, double *,
	     double *, long int *);
    extern double apser_(double *, double *, double *, double *), bpser_(double *, 
	    double *, double *, double *), basym_(double *, double *, double *, double *), 
	    fpser_(double *, double *, double *, double *);
    static double lambda;

/* ----------------------------------------------------------------------- */

/*            EVALUATION OF THE INCOMPLETE BETA FUNCTION IX(A,B) */

/*                     -------------------- */

/*     IT IS ASSUMED THAT A AND B ARE NONNEGATIVE, AND THAT X .LE. 1 */
/*     AND Y = 1 - X.  BRATIO ASSIGNS W AND W1 THE VALUES */

/*                      W  = IX(A,B) */
/*                      W1 = 1 - IX(A,B) */

/*     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS. */
/*     IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND */
/*     W AND W1 ARE COMPUTED. OTHERWISE, IF AN ERROR IS DETECTED, */
/*     THEN W AND W1 ARE ASSIGNED THE VALUE 0 AND IERR IS SET TO */
/*     ONE OF THE FOLLOWING VALUES ... */

/*        IERR = 1  IF A OR B IS NEGATIVE */
/*        IERR = 2  IF A = B = 0 */
/*        IERR = 3  IF X .LT. 0 OR X .GT. 1 */
/*        IERR = 4  IF Y .LT. 0 OR Y .GT. 1 */
/*        IERR = 5  IF X + Y .NE. 1 */
/*        IERR = 6  IF X = A = 0 */
/*        IERR = 7  IF Y = B = 0 */

/* -------------------- */
/*     WRITTEN BY ALFRED H. MORRIS, JR. */
/*        NAVAL SURFACE WARFARE CENTER */
/*        DAHLGREN, VIRGINIA */
/*     REVISED ... NOV 1991 */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST */
/*            FLOATING POINT NUMBER FOR WHICH 1.0 + EPS .GT. 1.0 */

    eps = DBL_EPSILON;

/* ----------------------------------------------------------------------- */
    *w = 0.f;
    *w1 = 0.f;
    if (*a < 0.f || *b < 0.f) {
	goto L300;
    }
    if (*a == 0.f && *b == 0.f) {
	goto L310;
    }
    if (*x < 0.f || *x > 1.f) {
	goto L320;
    }
    if (*y < 0.f || *y > 1.f) {
	goto L330;
    }
    z__ = *x + *y - .5f - .5f;
    if (DABS(z__) > eps * 3.f) {
	goto L340;
    }

    *ierr = 0;
    if (*x == 0.f) {
	goto L200;
    }
    if (*y == 0.f) {
	goto L210;
    }
    if (*a == 0.f) {
	goto L211;
    }
    if (*b == 0.f) {
	goto L201;
    }

    eps = DMAX(eps,1e-15f);
    if (DMAX(*a,*b) < eps * .001f) {
	goto L230;
    }

    ind = 0;
    a0 = *a;
    b0 = *b;
    x0 = *x;
    y0 = *y;
    if (DMIN(a0,b0) > 1.f) {
	goto L30;
    }

/*             PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1 */

    if (*x <= .5f) {
	goto L10;
    }
    ind = 1;
    a0 = *b;
    b0 = *a;
    x0 = *y;
    y0 = *x;

L10:
/* Computing MIN */
    r__1 = eps, r__2 = eps * a0;
    if (b0 < DMIN(r__1,r__2)) {
	goto L80;
    }
/* Computing MIN */
    r__1 = eps, r__2 = eps * b0;
    if (a0 < DMIN(r__1,r__2) && b0 * x0 <= 1.f) {
	goto L90;
    }
    if (DMAX(a0,b0) > 1.f) {
	goto L20;
    }
    if (a0 >= DMIN(.2f,b0)) {
	goto L100;
    }
    d__1 = (double) x0;
    d__2 = (double) a0;
    if (pow(d__1, d__2) <= .9f) {
	goto L100;
    }
    if (x0 >= .3f) {
	goto L110;
    }
    n = 20;
    goto L130;

L20:
    if (b0 <= 1.f) {
	goto L100;
    }
    if (x0 >= .3f) {
	goto L110;
    }
    if (x0 >= .1f) {
	goto L21;
    }
    d__1 = (double) (x0 * b0);
    d__2 = (double) a0;
    if (pow(d__1, d__2) <= .7f) {
	goto L100;
    }
L21:
    if (b0 > 15.f) {
	goto L131;
    }
    n = 20;
    goto L130;

/*             PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1 */

L30:
    if (*a > *b) {
	goto L31;
    }
    lambda = *a - (*a + *b) * *x;
    goto L32;
L31:
    lambda = (*a + *b) * *y - *b;
L32:
    if (lambda >= 0.f) {
	goto L40;
    }
    ind = 1;
    a0 = *b;
    b0 = *a;
    x0 = *y;
    y0 = *x;
    lambda = DABS(lambda);

L40:
    if (b0 < 40.f && b0 * x0 <= .7f) {
	goto L100;
    }
    if (b0 < 40.f) {
	goto L140;
    }
    if (a0 > b0) {
	goto L50;
    }
    if (a0 <= 100.f) {
	goto L120;
    }
    if (lambda > a0 * .03f) {
	goto L120;
    }
    goto L180;
L50:
    if (b0 <= 100.f) {
	goto L120;
    }
    if (lambda > b0 * .03f) {
	goto L120;
    }
    goto L180;

/*            EVALUATION OF THE APPROPRIATE ALGORITHM */

L80:
    *w = fpser_(&a0, &b0, &x0, &eps);
    *w1 = .5f - *w + .5f;
    goto L220;

L90:
    *w1 = apser_(&a0, &b0, &x0, &eps);
    *w = .5f - *w1 + .5f;
    goto L220;

L100:
    *w = bpser_(&a0, &b0, &x0, &eps);
    *w1 = .5f - *w + .5f;
    goto L220;

L110:
    *w1 = bpser_(&b0, &a0, &y0, &eps);
    *w = .5f - *w1 + .5f;
    goto L220;

L120:
    r__1 = eps * 15.f;
    *w = bfrac_(&a0, &b0, &x0, &y0, &lambda, &r__1);
    *w1 = .5f - *w + .5f;
    goto L220;

L130:
    *w1 = bup_(&b0, &a0, &y0, &x0, &n, &eps);
    b0 += n;
L131:
    r__1 = eps * 15.f;
    bgrat_(&b0, &a0, &y0, &x0, w1, &r__1, &ierr1);
    *w = .5f - *w1 + .5f;
    goto L220;

L140:
    n = b0;
    b0 -= n;
    if (b0 != 0.f) {
	goto L141;
    }
    --n;
    b0 = 1.f;
L141:
    *w = bup_(&b0, &a0, &y0, &x0, &n, &eps);
    if (x0 > .7f) {
	goto L150;
    }
    *w += bpser_(&a0, &b0, &x0, &eps);
    *w1 = .5f - *w + .5f;
    goto L220;

L150:
    if (a0 > 15.f) {
	goto L151;
    }
    n = 20;
    *w += bup_(&a0, &b0, &x0, &y0, &n, &eps);
    a0 += n;
L151:
    r__1 = eps * 15.f;
    bgrat_(&a0, &b0, &x0, &y0, w, &r__1, &ierr1);
    *w1 = .5f - *w + .5f;
    goto L220;

L180:
    r__1 = eps * 100.f;
    *w = basym_(&a0, &b0, &lambda, &r__1);
    *w1 = .5f - *w + .5f;
    goto L220;

/*               TERMINATION OF THE PROCEDURE */

L200:
    if (*a == 0.f) {
	goto L350;
    }
L201:
    *w = 0.f;
    *w1 = 1.f;
    return 0;

L210:
    if (*b == 0.f) {
	goto L360;
    }
L211:
    *w = 1.f;
    *w1 = 0.f;
    return 0;

L220:
    if (ind == 0) {
	return 0;
    }
    t = *w;
    *w = *w1;
    *w1 = t;
    return 0;

/*           PROCEDURE FOR A AND B .LT. 1.E-3*EPS */

L230:
    *w = *b / (*a + *b);
    *w1 = *a / (*a + *b);
    return 0;

/*                       ERROR RETURN */

L300:
    *ierr = 1;
    return 0;
L310:
    *ierr = 2;
    return 0;
L320:
    *ierr = 3;
    return 0;
L330:
    *ierr = 4;
    return 0;
L340:
    *ierr = 5;
    return 0;
L350:
    *ierr = 6;
    return 0;
L360:
    *ierr = 7;
    return 0;
} /* bratio_ */

double fpser_(double *a, double *b, double *x, double *eps)
{
    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    double log(double), exp(double);

    /* Local variables */
    static double c__, s, t, an, tol;
    extern double exparg_(long int *);

/* ----------------------------------------------------------------------- */

/*                 EVALUATION OF I (A,B) */
/*                                X */

/*          FOR B .LT. MIN(EPS,EPS*A) AND X .LE. 0.5. */

/* ----------------------------------------------------------------------- */

/*                  SET  FPSER = X**A */

    ret_val = 1.f;
    if (*a <= *eps * .001f) {
	goto L10;
    }
    ret_val = 0.f;
    t = *a * log(*x);
    if (t < exparg_(&c__1)) {
	return ret_val;
    }
    ret_val = exp(t);

/*                NOTE THAT 1/B(A,B) = B */

L10:
    ret_val = *b / *a * ret_val;
    tol = *eps / *a;
    an = *a + 1.f;
    t = *x;
    s = t / an;
L20:
    an += 1.f;
    t = *x * t;
    c__ = t / an;
    s += c__;
    if (DABS(c__) > tol) {
	goto L20;
    }

    ret_val *= *a * s + 1.f;
    return ret_val;
} /* fpser_ */

double apser_(double *a, double *b, double *x, double *eps)
{
    /* Initialized data */

    static double g = .577215664901533f;

    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    double log(double);

    /* Local variables */
    static double c__, j, s, t, aj, bx;
    extern double psi_(double *);
    static double tol;

/* ----------------------------------------------------------------------- */
/*     APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR */
/*     A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. USED WHEN */
/*     A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED. */
/* ----------------------------------------------------------------------- */
/* -------------------- */
/* -------------------- */
    bx = *b * *x;
    t = *x - bx;
    if (*b * *eps > .02f) {
	goto L10;
    }
    c__ = log(*x) + psi_(b) + g + t;
    goto L20;
L10:
    c__ = log(bx) + g + t;

L20:
    tol = *eps * 5.f * DABS(c__);
    j = 1.f;
    s = 0.f;
L30:
    j += 1.f;
    t *= *x - bx / j;
    aj = t / j;
    s += aj;
    if (DABS(aj) > tol) {
	goto L30;
    }

    ret_val = -(*a) * (c__ + s);
    return ret_val;
} /* apser_ */

double bpser_(double *a, double *b, double *x, double *eps)
{
    /* System generated locals */
    long int i__1;
    double ret_val;
    double d__1, d__2;

    /* Builtin functions */
    double log(double), exp(double), pow(double, double);

    /* Local variables */
    static double c__;
    static long int i__, m;
    static double n, t, u, w, z__, a0, b0, apb, tol, sum;
    extern double gam1_(double *), gamln1_(double *), betaln_(double *, double *),
	     algdiv_(double *, double *);

/* ----------------------------------------------------------------------- */
/*     POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B .LE. 1 */
/*     OR B*X .LE. 0.7.  EPS IS THE TOLERANCE USED. */
/* ----------------------------------------------------------------------- */

    ret_val = 0.f;
    if (*x == 0.f) {
	return ret_val;
    }
/* ----------------------------------------------------------------------- */
/*            COMPUTE THE FACTOR X**A/(A*BETA(A,B)) */
/* ----------------------------------------------------------------------- */
    a0 = DMIN(*a,*b);
    if (a0 < 1.f) {
	goto L10;
    }
    z__ = *a * log(*x) - betaln_(a, b);
    ret_val = exp(z__) / *a;
    goto L70;
L10:
    b0 = DMAX(*a,*b);
    if (b0 >= 8.f) {
	goto L60;
    }
    if (b0 > 1.f) {
	goto L40;
    }

/*            PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1 */

    d__1 = (double) (*x);
    d__2 = (double) (*a);
    ret_val = pow(d__1, d__2);
    if (ret_val == 0.f) {
	return ret_val;
    }

    apb = *a + *b;
    if (apb > 1.f) {
	goto L20;
    }
    z__ = gam1_(&apb) + 1.f;
    goto L30;
L20:
    u = (double) (*a) + (double) (*b) - 1.;
    z__ = (gam1_(&u) + 1.f) / apb;

L30:
    c__ = (gam1_(a) + 1.f) * (gam1_(b) + 1.f) / z__;
    ret_val = ret_val * c__ * (*b / apb);
    goto L70;

/*         PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8 */

L40:
    u = gamln1_(&a0);
    m = b0 - 1.f;
    if (m < 1) {
	goto L50;
    }
    c__ = 1.f;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b0 += -1.f;
/* L41: */
	c__ *= b0 / (a0 + b0);
    }
    u = log(c__) + u;

L50:
    z__ = *a * log(*x) - u;
    b0 += -1.f;
    apb = a0 + b0;
    if (apb > 1.f) {
	goto L51;
    }
    t = gam1_(&apb) + 1.f;
    goto L52;
L51:
    u = (double) a0 + (double) b0 - 1.;
    t = (gam1_(&u) + 1.f) / apb;
L52:
    ret_val = exp(z__) * (a0 / *a) * (gam1_(&b0) + 1.f) / t;
    goto L70;

/*            PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8 */

L60:
    u = gamln1_(&a0) + algdiv_(&a0, &b0);
    z__ = *a * log(*x) - u;
    ret_val = a0 / *a * exp(z__);
L70:
    if (ret_val == 0.f || *a <= *eps * .1f) {
	return ret_val;
    }
/* ----------------------------------------------------------------------- */
/*                     COMPUTE THE SERIES */
/* ----------------------------------------------------------------------- */
    sum = 0.f;
    n = 0.f;
    c__ = 1.f;
    tol = *eps / *a;
L100:
    n += 1.f;
    c__ = c__ * (.5f - *b / n + .5f) * *x;
    w = c__ / (*a + n);
    sum += w;
    if (DABS(w) > tol) {
	goto L100;
    }
    ret_val *= *a * sum + 1.f;
    return ret_val;
} /* bpser_ */

double bup_(double *a, double *b, double *x, double *y, long int *n, double *eps)
{
    /* System generated locals */
    long int i__1;
    double ret_val, r__1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double d__;
    static long int i__, k;
    static double l, r__, t, w;
    static long int mu;
    static double ap1;
    static long int nm1, kp1;
    static double apb;
    extern double brcmp1_(long int *, double *, double *, double *, double *), 
	    exparg_(long int *);

/* ----------------------------------------------------------------------- */
/*     EVALUATION OF IX(A,B) - IX(A+N,B) WHERE N IS A POSITIVE INTEGER. */
/*     EPS IS THE TOLERANCE USED. */
/* ----------------------------------------------------------------------- */

/*          OBTAIN THE SCALING FACTOR EXP(-MU) AND */
/*             EXP(MU)*(X**A*Y**B/BETA(A,B))/A */

    apb = *a + *b;
    ap1 = *a + 1.f;
    mu = 0;
    d__ = 1.f;
    if (*n == 1 || *a < 1.f) {
	goto L10;
    }
    if (apb < ap1 * 1.1f) {
	goto L10;
    }
    mu = (r__1 = exparg_(&c__1), DABS(r__1));
    k = exparg_(&c__0);
    if (k < mu) {
	mu = k;
    }
    t = (double) mu;
    d__ = exp(-t);

L10:
    ret_val = brcmp1_(&mu, a, b, x, y) / *a;
    if (*n == 1 || ret_val == 0.f) {
	return ret_val;
    }
    nm1 = *n - 1;
    w = d__;

/*          LET K BE THE INDEX OF THE MAXIMUM TERM */

    k = 0;
    if (*b <= 1.f) {
	goto L40;
    }
    if (*y > 1e-4f) {
	goto L20;
    }
    k = nm1;
    goto L30;
L20:
    r__ = (*b - 1.f) * *x / *y - *a;
    if (r__ < 1.f) {
	goto L40;
    }
    k = nm1;
    t = (double) nm1;
    if (r__ < t) {
	k = r__;
    }

/*          ADD THE INCREASING TERMS OF THE SERIES */

L30:
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = (double) (i__ - 1);
	d__ = (apb + l) / (ap1 + l) * *x * d__;
	w += d__;
/* L31: */
    }
    if (k == nm1) {
	goto L50;
    }

/*          ADD THE REMAINING TERMS OF THE SERIES */

L40:
    kp1 = k + 1;
    i__1 = nm1;
    for (i__ = kp1; i__ <= i__1; ++i__) {
	l = (double) (i__ - 1);
	d__ = (apb + l) / (ap1 + l) * *x * d__;
	w += d__;
	if (d__ <= *eps * w) {
	    goto L50;
	}
/* L41: */
    }

/*               TERMINATE THE PROCEDURE */

L50:
    ret_val *= w;
    return ret_val;
} /* bup_ */

double bfrac_(double *a, double *b, double *x, double *y, double *lambda, double *eps)
{
    /* System generated locals */
    double ret_val, r__1;

    /* Local variables */
    static double c__, e, n, p, r__, s, t, w, c0, c1, r0, an, bn, yp1, anp1, 
	    bnp1, beta, alpha;
    extern double brcomp_(double *, double *, double *, double *);

/* ----------------------------------------------------------------------- */
/*     CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B .GT. 1. */
/*     IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B. */
/* ----------------------------------------------------------------------- */
/* -------------------- */
    ret_val = brcomp_(a, b, x, y);
    if (ret_val == 0.f) {
	return ret_val;
    }

    c__ = *lambda + 1.f;
    c0 = *b / *a;
    c1 = 1.f / *a + 1.f;
    yp1 = *y + 1.f;

    n = 0.f;
    p = 1.f;
    s = *a + 1.f;
    an = 0.f;
    bn = 1.f;
    anp1 = 1.f;
    bnp1 = c__ / c1;
    r__ = c1 / c__;

/*        CONTINUED FRACTION CALCULATION */

L10:
    n += 1.f;
    t = n / *a;
    w = n * (*b - n) * *x;
    e = *a / s;
    alpha = p * (p + c0) * e * e * (w * *x);
    e = (t + 1.f) / (c1 + t + t);
    beta = n + w / s + e * (c__ + n * yp1);
    p = t + 1.f;
    s += 2.f;

/*        UPDATE AN, BN, ANP1, AND BNP1 */

    t = alpha * an + beta * anp1;
    an = anp1;
    anp1 = t;
    t = alpha * bn + beta * bnp1;
    bn = bnp1;
    bnp1 = t;

    r0 = r__;
    r__ = anp1 / bnp1;
    if ((r__1 = r__ - r0, DABS(r__1)) <= *eps * r__) {
	goto L20;
    }

/*        RESCALE AN, BN, ANP1, AND BNP1 */

    an /= bnp1;
    bn /= bnp1;
    anp1 = r__;
    bnp1 = 1.f;
    goto L10;

/*                 TERMINATION */

L20:
    ret_val *= r__;
    return ret_val;
} /* bfrac_ */

double brcomp_(double *a, double *b, double *x, double *y)
{
    /* Initialized data */

    static double const__ = .398942280401433f;

    /* System generated locals */
    long int i__1;
    double ret_val, r__1;

    /* Builtin functions */
    double log(double), exp(double), sqrt(double);

    /* Local variables */
    static double c__, e, h__;
    static long int i__, n;
    static double t, u, v, z__, a0, b0, x0, y0, apb, lnx, lny;
    extern double gam1_(double *), rlog1_(double *), bcorr_(double *, double *), 
	    gamln1_(double *);
    static double lambda;
    extern double betaln_(double *, double *), algdiv_(double *, double *), 
	    alnrel_(double *);

/* ----------------------------------------------------------------------- */
/*               EVALUATION OF X**A*Y**B/BETA(A,B) */
/* ----------------------------------------------------------------------- */
/* ----------------- */
/*     CONST = 1/SQRT(2*PI) */
/* ----------------- */

    ret_val = 0.f;
    if (*x == 0.f || *y == 0.f) {
	return ret_val;
    }
    a0 = DMIN(*a,*b);
    if (a0 >= 8.f) {
	goto L100;
    }

    if (*x > .375f) {
	goto L10;
    }
    lnx = log(*x);
    r__1 = -(*x);
    lny = alnrel_(&r__1);
    goto L20;
L10:
    if (*y > .375f) {
	goto L11;
    }
    r__1 = -(*y);
    lnx = alnrel_(&r__1);
    lny = log(*y);
    goto L20;
L11:
    lnx = log(*x);
    lny = log(*y);

L20:
    z__ = *a * lnx + *b * lny;
    if (a0 < 1.f) {
	goto L30;
    }
    z__ -= betaln_(a, b);
    ret_val = exp(z__);
    return ret_val;
/* ----------------------------------------------------------------------- */
/*              PROCEDURE FOR A .LT. 1 OR B .LT. 1 */
/* ----------------------------------------------------------------------- */
L30:
    b0 = DMAX(*a,*b);
    if (b0 >= 8.f) {
	goto L80;
    }
    if (b0 > 1.f) {
	goto L60;
    }

/*                   ALGORITHM FOR B0 .LE. 1 */

    ret_val = exp(z__);
    if (ret_val == 0.f) {
	return ret_val;
    }

    apb = *a + *b;
    if (apb > 1.f) {
	goto L40;
    }
    z__ = gam1_(&apb) + 1.f;
    goto L50;
L40:
    u = (double) (*a) + (double) (*b) - 1.;
    z__ = (gam1_(&u) + 1.f) / apb;

L50:
    c__ = (gam1_(a) + 1.f) * (gam1_(b) + 1.f) / z__;
    ret_val = ret_val * (a0 * c__) / (a0 / b0 + 1.f);
    return ret_val;

/*                ALGORITHM FOR 1 .LT. B0 .LT. 8 */

L60:
    u = gamln1_(&a0);
    n = b0 - 1.f;
    if (n < 1) {
	goto L70;
    }
    c__ = 1.f;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b0 += -1.f;
	c__ *= b0 / (a0 + b0);
/* L61: */
    }
    u = log(c__) + u;

L70:
    z__ -= u;
    b0 += -1.f;
    apb = a0 + b0;
    if (apb > 1.f) {
	goto L71;
    }
    t = gam1_(&apb) + 1.f;
    goto L72;
L71:
    u = (double) a0 + (double) b0 - 1.;
    t = (gam1_(&u) + 1.f) / apb;
L72:
    ret_val = a0 * exp(z__) * (gam1_(&b0) + 1.f) / t;
    return ret_val;

/*                   ALGORITHM FOR B0 .GE. 8 */

L80:
    u = gamln1_(&a0) + algdiv_(&a0, &b0);
    ret_val = a0 * exp(z__ - u);
    return ret_val;
/* ----------------------------------------------------------------------- */
/*              PROCEDURE FOR A .GE. 8 AND B .GE. 8 */
/* ----------------------------------------------------------------------- */
L100:
    if (*a > *b) {
	goto L101;
    }
    h__ = *a / *b;
    x0 = h__ / (h__ + 1.f);
    y0 = 1.f / (h__ + 1.f);
    lambda = *a - (*a + *b) * *x;
    goto L110;
L101:
    h__ = *b / *a;
    x0 = 1.f / (h__ + 1.f);
    y0 = h__ / (h__ + 1.f);
    lambda = (*a + *b) * *y - *b;

L110:
    e = -lambda / *a;
    if (DABS(e) > .6f) {
	goto L111;
    }
    u = rlog1_(&e);
    goto L120;
L111:
    u = e - log(*x / x0);

L120:
    e = lambda / *b;
    if (DABS(e) > .6f) {
	goto L121;
    }
    v = rlog1_(&e);
    goto L130;
L121:
    v = e - log(*y / y0);

L130:
    z__ = exp(-(*a * u + *b * v));
    ret_val = const__ * sqrt(*b * x0) * z__ * exp(-bcorr_(a, b));
    return ret_val;
} /* brcomp_ */

double brcmp1_(long int *mu, double *a, double *b, double *x, double *y)
{
    /* Initialized data */

    static double const__ = .398942280401433f;

    /* System generated locals */
    long int i__1;
    double ret_val, r__1;

    /* Builtin functions */
    double log(double), sqrt(double), exp(double);

    /* Local variables */
    static double c__, e, h__;
    static long int i__, n;
    static double t, u, v, z__, a0, b0, x0, y0, apb, lnx, lny;
    extern double gam1_(double *), esum_(long int *, double *), rlog1_(double *),
	     bcorr_(double *, double *), gamln1_(double *);
    static double lambda;
    extern double betaln_(double *, double *), algdiv_(double *, double *), 
	    alnrel_(double *);

/* ----------------------------------------------------------------------- */
/*          EVALUATION OF  EXP(MU) * (X**A*Y**B/BETA(A,B)) */
/* ----------------------------------------------------------------------- */
/* ----------------- */
/*     CONST = 1/SQRT(2*PI) */
/* ----------------- */

    a0 = DMIN(*a,*b);
    if (a0 >= 8.f) {
	goto L100;
    }

    if (*x > .375f) {
	goto L10;
    }
    lnx = log(*x);
    r__1 = -(*x);
    lny = alnrel_(&r__1);
    goto L20;
L10:
    if (*y > .375f) {
	goto L11;
    }
    r__1 = -(*y);
    lnx = alnrel_(&r__1);
    lny = log(*y);
    goto L20;
L11:
    lnx = log(*x);
    lny = log(*y);

L20:
    z__ = *a * lnx + *b * lny;
    if (a0 < 1.f) {
	goto L30;
    }
    z__ -= betaln_(a, b);
    ret_val = esum_(mu, &z__);
    return ret_val;
/* ----------------------------------------------------------------------- */
/*              PROCEDURE FOR A .LT. 1 OR B .LT. 1 */
/* ----------------------------------------------------------------------- */
L30:
    b0 = DMAX(*a,*b);
    if (b0 >= 8.f) {
	goto L80;
    }
    if (b0 > 1.f) {
	goto L60;
    }

/*                   ALGORITHM FOR B0 .LE. 1 */

    ret_val = esum_(mu, &z__);
    if (ret_val == 0.f) {
	return ret_val;
    }

    apb = *a + *b;
    if (apb > 1.f) {
	goto L40;
    }
    z__ = gam1_(&apb) + 1.f;
    goto L50;
L40:
    u = (double) (*a) + (double) (*b) - 1.;
    z__ = (gam1_(&u) + 1.f) / apb;

L50:
    c__ = (gam1_(a) + 1.f) * (gam1_(b) + 1.f) / z__;
    ret_val = ret_val * (a0 * c__) / (a0 / b0 + 1.f);
    return ret_val;

/*                ALGORITHM FOR 1 .LT. B0 .LT. 8 */

L60:
    u = gamln1_(&a0);
    n = b0 - 1.f;
    if (n < 1) {
	goto L70;
    }
    c__ = 1.f;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b0 += -1.f;
	c__ *= b0 / (a0 + b0);
/* L61: */
    }
    u = log(c__) + u;

L70:
    z__ -= u;
    b0 += -1.f;
    apb = a0 + b0;
    if (apb > 1.f) {
	goto L71;
    }
    t = gam1_(&apb) + 1.f;
    goto L72;
L71:
    u = (double) a0 + (double) b0 - 1.;
    t = (gam1_(&u) + 1.f) / apb;
L72:
    ret_val = a0 * esum_(mu, &z__) * (gam1_(&b0) + 1.f) / t;
    return ret_val;

/*                   ALGORITHM FOR B0 .GE. 8 */

L80:
    u = gamln1_(&a0) + algdiv_(&a0, &b0);
    r__1 = z__ - u;
    ret_val = a0 * esum_(mu, &r__1);
    return ret_val;
/* ----------------------------------------------------------------------- */
/*              PROCEDURE FOR A .GE. 8 AND B .GE. 8 */
/* ----------------------------------------------------------------------- */
L100:
    if (*a > *b) {
	goto L101;
    }
    h__ = *a / *b;
    x0 = h__ / (h__ + 1.f);
    y0 = 1.f / (h__ + 1.f);
    lambda = *a - (*a + *b) * *x;
    goto L110;
L101:
    h__ = *b / *a;
    x0 = 1.f / (h__ + 1.f);
    y0 = h__ / (h__ + 1.f);
    lambda = (*a + *b) * *y - *b;

L110:
    e = -lambda / *a;
    if (DABS(e) > .6f) {
	goto L111;
    }
    u = rlog1_(&e);
    goto L120;
L111:
    u = e - log(*x / x0);

L120:
    e = lambda / *b;
    if (DABS(e) > .6f) {
	goto L121;
    }
    v = rlog1_(&e);
    goto L130;
L121:
    v = e - log(*y / y0);

L130:
    r__1 = -(*a * u + *b * v);
    z__ = esum_(mu, &r__1);
    ret_val = const__ * sqrt(*b * x0) * z__ * exp(-bcorr_(a, b));
    return ret_val;
} /* brcmp1_ */

/* Subroutine */ int bgrat_(double *a, double *b, double *x, double *y, double *w, double 
	*eps, long int *ierr)
{
    /* System generated locals */
    long int i__1;
    double r__1;

    /* Builtin functions */
    double log(double), exp(double);

    /* Local variables */
    static double c__[30], d__[30];
    static long int i__;
    static double j, l;
    static long int n;
    static double p, q, r__, s, t, u, v, z__, n2, t2, dj, cn, nu, bm1;
    static long int nm1;
    static double lnx, sum;
    extern double gam1_(double *);
    static double bp2n, coef;
    extern /* Subroutine */ int grat1_(double *, double *, double *, double *, double *,
	     double *);
    extern double algdiv_(double *, double *), alnrel_(double *);

/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR IX(A,B) WHEN A IS LARGER THAN B. */
/*     THE RESULT OF THE EXPANSION IS ADDED TO W. IT IS ASSUMED */
/*     THAT A .GE. 15 AND B .LE. 1.  EPS IS THE TOLERANCE USED. */
/*     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS. */
/* ----------------------------------------------------------------------- */

    bm1 = *b - .5f - .5f;
    nu = *a + bm1 * .5f;
    if (*y > .375f) {
	goto L10;
    }
    r__1 = -(*y);
    lnx = alnrel_(&r__1);
    goto L11;
L10:
    lnx = log(*x);
L11:
    z__ = -nu * lnx;
    if (*b * z__ == 0.f) {
	goto L100;
    }

/*                 COMPUTATION OF THE EXPANSION */
/*                 SET R = EXP(-Z)*Z**B/GAMMA(B) */

    r__ = *b * (gam1_(b) + 1.f) * exp(*b * log(z__));
    r__ = r__ * exp(*a * lnx) * exp(bm1 * .5f * lnx);
    u = algdiv_(b, a) + *b * log(nu);
    u = r__ * exp(-u);
    if (u == 0.f) {
	goto L100;
    }
    grat1_(b, &z__, &r__, &p, &q, eps);

/* Computing 2nd power */
    r__1 = 1.f / nu;
    v = r__1 * r__1 * .25f;
    t2 = lnx * .25f * lnx;
    l = *w / u;
    j = q / r__;
    sum = j;
    t = 1.f;
    cn = 1.f;
    n2 = 0.f;
    for (n = 1; n <= 30; ++n) {
	bp2n = *b + n2;
	j = (bp2n * (bp2n + 1.f) * j + (z__ + bp2n + 1.f) * t) * v;
	n2 += 2.f;
	t *= t2;
	cn /= n2 * (n2 + 1.f);
	c__[n - 1] = cn;
	s = 0.f;
	if (n == 1) {
	    goto L21;
	}
	nm1 = n - 1;
	coef = *b - n;
	i__1 = nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s += coef * c__[i__ - 1] * d__[n - i__ - 1];
/* L20: */
	    coef += *b;
	}
L21:
	d__[n - 1] = bm1 * cn + s / n;
	dj = d__[n - 1] * j;
	sum += dj;
	if (sum <= 0.f) {
	    goto L100;
	}
	if (DABS(dj) <= *eps * (sum + l)) {
	    goto L30;
	}
/* L22: */
    }

/*                    ADD THE RESULTS TO W */

L30:
    *ierr = 0;
    *w += u * sum;
    return 0;

/*               THE EXPANSION CANNOT BE COMPUTED */

L100:
    *ierr = 1;
    return 0;
} /* bgrat_ */

/* Subroutine */ int grat1_(double *a, double *x, double *r__, double *p, double *q, 
	double *eps)
{
    /* System generated locals */
    double r__1;

    /* Builtin functions */
    double log(double), exp(double), sqrt(double);

    /* Local variables */
    static double c__, g, h__, j, l, t, w, z__, an, am0, an0, a2n, b2n, cma;
    extern double erf_(double *);
    static double tol, sum;
    extern double gam1_(double *);
    static double a2nm1, b2nm1;
    extern double rexp_(double *), erfc1_(long int *, double *);

/* ----------------------------------------------------------------------- */
/*        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS */
/*                      P(A,X) AND Q(A,X) */

/*     IT IS ASSUMED THAT A .LE. 1.  EPS IS THE TOLERANCE TO BE USED. */
/*     THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A). */
/* ----------------------------------------------------------------------- */
    if (*a * *x == 0.f) {
	goto L130;
    }
    if (*a == .5f) {
	goto L120;
    }
    if (*x < 1.1f) {
	goto L10;
    }
    goto L50;

/*             TAYLOR SERIES FOR P(A,X)/X**A */

L10:
    an = 3.f;
    c__ = *x;
    sum = *x / (*a + 3.f);
    tol = *eps * .1f / (*a + 1.f);
L11:
    an += 1.f;
    c__ = -c__ * (*x / an);
    t = c__ / (*a + an);
    sum += t;
    if (DABS(t) > tol) {
	goto L11;
    }
    j = *a * *x * ((sum / 6.f - .5f / (*a + 2.f)) * *x + 1.f / (*a + 1.f));

    z__ = *a * log(*x);
    h__ = gam1_(a);
    g = h__ + 1.f;
    if (*x < .25f) {
	goto L20;
    }
    if (*a < *x / 2.59f) {
	goto L40;
    }
    goto L30;
L20:
    if (z__ > -.13394f) {
	goto L40;
    }

L30:
    w = exp(z__);
    *p = w * g * (.5f - j + .5f);
    *q = .5f - *p + .5f;
    return 0;

L40:
    l = rexp_(&z__);
    w = l + .5f + .5f;
    *q = (w * j - l) * g - h__;
    if (*q < 0.f) {
	goto L110;
    }
    *p = .5f - *q + .5f;
    return 0;

/*              CONTINUED FRACTION EXPANSION */

L50:
    a2nm1 = 1.f;
    a2n = 1.f;
    b2nm1 = *x;
    b2n = *x + (1.f - *a);
    c__ = 1.f;
L51:
    a2nm1 = *x * a2n + c__ * a2nm1;
    b2nm1 = *x * b2n + c__ * b2nm1;
    am0 = a2nm1 / b2nm1;
    c__ += 1.f;
    cma = c__ - *a;
    a2n = a2nm1 + cma * a2n;
    b2n = b2nm1 + cma * b2n;
    an0 = a2n / b2n;
    if ((r__1 = an0 - am0, DABS(r__1)) >= *eps * an0) {
	goto L51;
    }
    *q = *r__ * an0;
    *p = .5f - *q + .5f;
    return 0;

/*                SPECIAL CASES */

L100:
    *p = 0.f;
    *q = 1.f;
    return 0;

L110:
    *p = 1.f;
    *q = 0.f;
    return 0;

L120:
    if (*x >= .25f) {
	goto L121;
    }
    r__1 = sqrt(*x);
    *p = erf_(&r__1);
    *q = .5f - *p + .5f;
    return 0;
L121:
    r__1 = sqrt(*x);
    *q = erfc1_(&c__0, &r__1);
    *p = .5f - *q + .5f;
    return 0;

L130:
    if (*x <= *a) {
	goto L100;
    }
    goto L110;
} /* grat1_ */

double basym_(double *a, double *b, double *lambda, double *eps)
{
    /* Initialized data */

    static long int num = 20;
    static double e0 = 1.12837916709551f;
    static double e1 = .353553390593274f;

    /* System generated locals */
    long int i__1, i__2, i__3, i__4;
    double ret_val, r__1, r__2;

    /* Builtin functions */
    double sqrt(double), exp(double);

    /* Local variables */
    static double c__[21], d__[21], f, h__;
    static long int i__, j, m, n;
    static double r__, s, t, u, w, z__, a0[21], b0[21], j0, j1, h2, r0, r1, t0, 
	    t1, w0, z0, z2, hn, zn;
    static long int im1, mm1, np1, imj, mmj;
    static double sum, znm1, bsum, dsum;
    extern double erfc1_(long int *, double *), rlog1_(double *), bcorr_(double *
	    , double *);

/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B. */
/*     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED. */
/*     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT */
/*     A AND B ARE GREATER THAN OR EQUAL TO 15. */
/* ----------------------------------------------------------------------- */
/* ------------------------ */
/*     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP */
/*            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN. */
/*            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1. */

/* ------------------------ */
/*     E0 = 2/SQRT(PI) */
/*     E1 = 2**(-3/2) */
/* ------------------------ */
/* ------------------------ */
    ret_val = 0.f;
    if (*a >= *b) {
	goto L10;
    }
    h__ = *a / *b;
    r0 = 1.f / (h__ + 1.f);
    r1 = (*b - *a) / *b;
    w0 = 1.f / sqrt(*a * (h__ + 1.f));
    goto L20;
L10:
    h__ = *b / *a;
    r0 = 1.f / (h__ + 1.f);
    r1 = (*b - *a) / *a;
    w0 = 1.f / sqrt(*b * (h__ + 1.f));

L20:
    r__1 = -(*lambda) / *a;
    r__2 = *lambda / *b;
    f = *a * rlog1_(&r__1) + *b * rlog1_(&r__2);
    t = exp(-f);
    if (t == 0.f) {
	return ret_val;
    }
    z0 = sqrt(f);
    z__ = z0 / e1 * .5f;
    z2 = f + f;

    a0[0] = r1 * .66666666666666663f;
    c__[0] = a0[0] * -.5f;
    d__[0] = -c__[0];
    j0 = .5f / e0 * erfc1_(&c__1, &z0);
    j1 = e1;
    sum = j0 + d__[0] * w0 * j1;

    s = 1.f;
    h2 = h__ * h__;
    hn = 1.f;
    w = w0;
    znm1 = z__;
    zn = z2;
    i__1 = num;
    for (n = 2; n <= i__1; n += 2) {
	hn = h2 * hn;
	a0[n - 1] = r0 * 2.f * (h__ * hn + 1.f) / (n + 2.f);
	np1 = n + 1;
	s += hn;
	a0[np1 - 1] = r1 * 2.f * s / (n + 3.f);

	i__2 = np1;
	for (i__ = n; i__ <= i__2; ++i__) {
	    r__ = (i__ + 1.f) * -.5f;
	    b0[0] = r__ * a0[0];
	    i__3 = i__;
	    for (m = 2; m <= i__3; ++m) {
		bsum = 0.f;
		mm1 = m - 1;
		i__4 = mm1;
		for (j = 1; j <= i__4; ++j) {
		    mmj = m - j;
/* L30: */
		    bsum += (j * r__ - mmj) * a0[j - 1] * b0[mmj - 1];
		}
/* L31: */
		b0[m - 1] = r__ * a0[m - 1] + bsum / m;
	    }
	    c__[i__ - 1] = b0[i__ - 1] / (i__ + 1.f);

	    dsum = 0.f;
	    im1 = i__ - 1;
	    i__3 = im1;
	    for (j = 1; j <= i__3; ++j) {
		imj = i__ - j;
/* L40: */
		dsum += d__[imj - 1] * c__[j - 1];
	    }
/* L41: */
	    d__[i__ - 1] = -(dsum + c__[i__ - 1]);
	}

	j0 = e1 * znm1 + (n - 1.f) * j0;
	j1 = e1 * zn + n * j1;
	znm1 = z2 * znm1;
	zn = z2 * zn;
	w = w0 * w;
	t0 = d__[n - 1] * w * j0;
	w = w0 * w;
	t1 = d__[np1 - 1] * w * j1;
	sum += t0 + t1;
	if (DABS(t0) + DABS(t1) <= *eps * sum) {
	    goto L60;
	}
/* L50: */
    }

L60:
    u = exp(-bcorr_(a, b));
    ret_val = e0 * t * u * sum;
    return ret_val;
} /* basym_ */


double exparg_(long int *l)
{
    /* System generated locals */
    double lnb = .69314718055995f;
    int m;
    if(*l==0) {
      m = DBL_MAX_EXP;
      return m * lnb * .99999;
    }
    m = DBL_MIN_EXP -1;
    return m * lnb * .99999;
} /* exparg_ */

double esum_(long int *mu, double *x)
{
    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double w;

/* ----------------------------------------------------------------------- */
/*                    EVALUATION OF EXP(MU + X) */
/* ----------------------------------------------------------------------- */
    if (*x > 0.f) {
	goto L10;
    }

    if (*mu < 0) {
	goto L20;
    }
    w = *mu + *x;
    if (w > 0.f) {
	goto L20;
    }
    ret_val = exp(w);
    return ret_val;

L10:
    if (*mu > 0) {
	goto L20;
    }
    w = *mu + *x;
    if (w < 0.f) {
	goto L20;
    }
    ret_val = exp(w);
    return ret_val;

L20:
    w = (double) (*mu);
    ret_val = exp(w) * exp(*x);
    return ret_val;
} /* esum_ */

double rexp_(double *x)
{
    /* Initialized data */

    static double p1 = 9.14041914819518e-10f;
    static double p2 = .0238082361044469f;
    static double q1 = -.499999999085958f;
    static double q2 = .107141568980644f;
    static double q3 = -.0119041179760821f;
    static double q4 = 5.95130811860248e-4f;

    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double w;

/* ----------------------------------------------------------------------- */
/*            EVALUATION OF THE FUNCTION EXP(X) - 1 */
/* ----------------------------------------------------------------------- */
/* ----------------------- */
    if (DABS(*x) > .15f) {
	goto L10;
    }
    ret_val = *x * (((p2 * *x + p1) * *x + 1.f) / ((((q4 * *x + q3) * *x + q2)
	     * *x + q1) * *x + 1.f));
    return ret_val;

L10:
    w = exp(*x);
    if (*x > 0.f) {
	goto L20;
    }
    ret_val = w - .5f - .5f;
    return ret_val;
L20:
    ret_val = w * (.5f - 1.f / w + .5f);
    return ret_val;
} /* rexp_ */

double alnrel_(double *a)
{
    /* Initialized data */

    static double p1 = -1.29418923021993f;
    static double p2 = .405303492862024f;
    static double p3 = -.0178874546012214f;
    static double q1 = -1.62752256355323f;
    static double q2 = .747811014037616f;
    static double q3 = -.0845104217945565f;

    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    double log(double);

    /* Local variables */
    static double t, w, x, t2;

/* ----------------------------------------------------------------------- */
/*            EVALUATION OF THE FUNCTION LN(1 + A) */
/* ----------------------------------------------------------------------- */
/* -------------------------- */
    if (DABS(*a) > .375f) {
	goto L10;
    }
    t = *a / (*a + 2.f);
    t2 = t * t;
    w = (((p3 * t2 + p2) * t2 + p1) * t2 + 1.f) / (((q3 * t2 + q2) * t2 + q1) 
	    * t2 + 1.f);
    ret_val = t * 2.f * w;
    return ret_val;

L10:
    x = (double) (*a) + 1.;
    ret_val = log(x);
    return ret_val;
} /* alnrel_ */

double rlog1_(double *x)
{
    /* Initialized data */

    static double a = .0566749439387324f;
    static double b = .0456512608815524f;
    static double p0 = .333333333333333f;
    static double p1 = -.224696413112536f;
    static double p2 = .00620886815375787f;
    static double q1 = -1.27408923933623f;
    static double q2 = .354508718369557f;

    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    double log(double);

    /* Local variables */
    static double h__, r__, t, w, w1;

/* ----------------------------------------------------------------------- */
/*             EVALUATION OF THE FUNCTION X - LN(1 + X) */
/* ----------------------------------------------------------------------- */
/* ------------------------ */
/* ------------------------ */
    if (*x < -.39f || *x > .57f) {
	goto L100;
    }
    if (*x < -.18f) {
	goto L10;
    }
    if (*x > .18f) {
	goto L20;
    }

/*              ARGUMENT REDUCTION */

    h__ = *x;
    w1 = 0.f;
    goto L30;

L10:
    h__ = (double) (*x) + .3;
    h__ /= .7f;
    w1 = a - h__ * .3f;
    goto L30;

L20:
    h__ = (double) (*x) * .75 - .25;
    w1 = b + h__ / 3.f;

/*               SERIES EXPANSION */

L30:
    r__ = h__ / (h__ + 2.f);
    t = r__ * r__;
    w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.f);
    ret_val = t * 2.f * (1.f / (1.f - r__) - r__ * w) + w1;
    return ret_val;


L100:
    w = *x + .5f + .5f;
    ret_val = *x - log(w);
    return ret_val;
} /* rlog1_ */

double erf_(double *x)
{
    /* Initialized data */

    static double c__ = .564189583547756f;
    static double a[5] = { 7.7105849500132e-5f,-.00133733772997339f,
	    .0323076579225834f,.0479137145607681f,.128379167095513f };
    static double b[3] = { .00301048631703895f,.0538971687740286f,
	    .375795757275549f };
    static double p[8] = { -1.36864857382717e-7f,.564195517478974f,
	    7.21175825088309f,43.1622272220567f,152.98928504694f,
	    339.320816734344f,451.918953711873f,300.459261020162f };
    static double q[8] = { 1.f,12.7827273196294f,77.0001529352295f,
	    277.585444743988f,638.980264465631f,931.35409485061f,
	    790.950925327898f,300.459260956983f };
    static double r__[5] = { 2.10144126479064f,26.2370141675169f,
	    21.3688200555087f,4.6580782871847f,.282094791773523f };
    static double s[4] = { 94.153775055546f,187.11481179959f,99.0191814623914f,
	    18.0124575948747f };

    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    double exp(double), r_sign(double *, double *);

    /* Local variables */
    static double t, x2, ax, bot, top;

/* ----------------------------------------------------------------------- */
/*             EVALUATION OF THE REAL ERROR FUNCTION */
/* ----------------------------------------------------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
    ax = DABS(*x);
    if (ax > .5f) {
	goto L10;
    }
    t = *x * *x;
    top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.f;
    bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.f;
    ret_val = *x * (top / bot);
    return ret_val;

L10:
    if (ax > 4.f) {
	goto L20;
    }
    top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax 
	    + p[5]) * ax + p[6]) * ax + p[7];
    bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax 
	    + q[5]) * ax + q[6]) * ax + q[7];
    ret_val = .5f - exp(-(*x) * *x) * top / bot + .5f;
    if (*x < 0.f) {
	ret_val = -ret_val;
    }
    return ret_val;

L20:
    if (ax >= 5.8f) {
	goto L30;
    }
    x2 = *x * *x;
    t = 1.f / x2;
    top = (((r__[0] * t + r__[1]) * t + r__[2]) * t + r__[3]) * t + r__[4];
    bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.f;
    ret_val = (c__ - top / (x2 * bot)) / ax;
    ret_val = .5f - exp(-x2) * ret_val + .5f;
    if (*x < 0.f) {
	ret_val = -ret_val;
    }
    return ret_val;

L30:
    ret_val = r_sign(&c_b188, x);
    return ret_val;
} /* erf_ */

double erfc1_(long int *ind, double *x)
{
    /* Initialized data */

    static double c__ = .564189583547756f;
    static double a[5] = { 7.7105849500132e-5f,-.00133733772997339f,
	    .0323076579225834f,.0479137145607681f,.128379167095513f };
    static double b[3] = { .00301048631703895f,.0538971687740286f,
	    .375795757275549f };
    static double p[8] = { -1.36864857382717e-7f,.564195517478974f,
	    7.21175825088309f,43.1622272220567f,152.98928504694f,
	    339.320816734344f,451.918953711873f,300.459261020162f };
    static double q[8] = { 1.f,12.7827273196294f,77.0001529352295f,
	    277.585444743988f,638.980264465631f,931.35409485061f,
	    790.950925327898f,300.459260956983f };
    static double r__[5] = { 2.10144126479064f,26.2370141675169f,
	    21.3688200555087f,4.6580782871847f,.282094791773523f };
    static double s[4] = { 94.153775055546f,187.11481179959f,99.0191814623914f,
	    18.0124575948747f };

    /* System generated locals */
    double ret_val, r__1;

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    static double e, t;
    static double w;
    static double ax, bot, top;
    extern double exparg_(long int *);

/* ----------------------------------------------------------------------- */
/*         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION */

/*          ERFC1(IND,X) = ERFC(X)            IF IND = 0 */
/*          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE */
/* ----------------------------------------------------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */

/*                     ABS(X) .LE. 0.5 */

    ax = DABS(*x);
    if (ax > .5f) {
	goto L10;
    }
    t = *x * *x;
    top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.f;
    bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.f;
    ret_val = .5f - *x * (top / bot) + .5f;
    if (*ind != 0) {
	ret_val = exp(t) * ret_val;
    }
    return ret_val;

/*                  0.5 .LT. ABS(X) .LE. 4 */

L10:
    if (ax > 4.f) {
	goto L20;
    }
    top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax 
	    + p[5]) * ax + p[6]) * ax + p[7];
    bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax 
	    + q[5]) * ax + q[6]) * ax + q[7];
    ret_val = top / bot;
    goto L40;

/*                      ABS(X) .GT. 4 */

L20:
    if (*x <= -5.6f) {
	goto L50;
    }
    if (*ind != 0) {
	goto L30;
    }
    if (*x > 100.f) {
	goto L60;
    }
    if (*x * *x > -exparg_(&c__1)) {
	goto L60;
    }

L30:
/* Computing 2nd power */
    r__1 = 1.f / *x;
    t = r__1 * r__1;
    top = (((r__[0] * t + r__[1]) * t + r__[2]) * t + r__[3]) * t + r__[4];
    bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.f;
    ret_val = (c__ - t * top / bot) / ax;

/*                      FINAL ASSEMBLY */

L40:
    if (*ind == 0) {
	goto L41;
    }
    if (*x < 0.f) {
	ret_val = exp(*x * *x) * 2.f - ret_val;
    }
    return ret_val;
L41:
    w = (double) (*x) * (double) (*x);
    t = w;
    e = w - (double) t;
    ret_val = (.5f - e + .5f) * exp(-t) * ret_val;
    if (*x < 0.f) {
	ret_val = 2.f - ret_val;
    }
    return ret_val;

/*             LIMIT VALUE FOR LARGE NEGATIVE X */

L50:
    ret_val = 2.f;
    if (*ind != 0) {
	ret_val = exp(*x * *x) * 2.f;
    }
    return ret_val;

/*             LIMIT VALUE FOR LARGE POSITIVE X */
/*                       WHEN IND = 0 */

L60:
    ret_val = 0.f;
    return ret_val;
} /* erfc1_ */

double gam1_(double *a)
{
    /* Initialized data */

    static double p[7] = { .577215664901533f,-.409078193005776f,
	    -.230975380857675f,.0597275330452234f,.0076696818164949f,
	    -.00514889771323592f,5.89597428611429e-4f };
    static double q[5] = { 1.f,.427569613095214f,.158451672430138f,
	    .0261132021441447f,.00423244297896961f };
    static double r__[9] = { -.422784335098468f,-.771330383816272f,
	    -.244757765222226f,.118378989872749f,9.30357293360349e-4f,
	    -.0118290993445146f,.00223047661158249f,2.66505979058923e-4f,
	    -1.32674909766242e-4f };
    static double s1 = .273076135303957f;
    static double s2 = .0559398236957378f;

    /* System generated locals */
    double ret_val;

    /* Local variables */
    static double d__, t, w, bot, top;

/*     ------------------------------------------------------------------ */
/*     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 .LE. A .LE. 1.5 */
/*     ------------------------------------------------------------------ */
/*     ------------------- */
/*     ------------------- */
/*     ------------------- */
/*     ------------------- */
/*     ------------------- */
    t = *a;
    d__ = *a - .5f;
    if (d__ > 0.f) {
	t = d__ - .5f;
    }
    if (t < 0.f) {
	goto L30;
    } else if (t == 0) {
	goto L10;
    } else {
	goto L20;
    }

L10:
    ret_val = 0.f;
    return ret_val;

L20:
    top = (((((p[6] * t + p[5]) * t + p[4]) * t + p[3]) * t + p[2]) * t + p[1]
	    ) * t + p[0];
    bot = (((q[4] * t + q[3]) * t + q[2]) * t + q[1]) * t + 1.f;
    w = top / bot;
    if (d__ > 0.f) {
	goto L21;
    }
    ret_val = *a * w;
    return ret_val;
L21:
    ret_val = t / *a * (w - .5f - .5f);
    return ret_val;

L30:
    top = (((((((r__[8] * t + r__[7]) * t + r__[6]) * t + r__[5]) * t + r__[4]
	    ) * t + r__[3]) * t + r__[2]) * t + r__[1]) * t + r__[0];
    bot = (s2 * t + s1) * t + 1.f;
    w = top / bot;
    if (d__ > 0.f) {
	goto L31;
    }
    ret_val = *a * (w + .5f + .5f);
    return ret_val;
L31:
    ret_val = t * w / *a;
    return ret_val;
} /* gam1_ */

double gamln1_(double *a)
{
    /* Initialized data */

    static double p0 = .577215664901533f;
    static double p1 = .844203922187225f;
    static double p2 = -.168860593646662f;
    static double p3 = -.780427615533591f;
    static double p4 = -.402055799310489f;
    static double p5 = -.0673562214325671f;
    static double p6 = -.00271935708322958f;
    static double q1 = 2.88743195473681f;
    static double q2 = 3.12755088914843f;
    static double q3 = 1.56875193295039f;
    static double q4 = .361951990101499f;
    static double q5 = .0325038868253937f;
    static double q6 = 6.67465618796164e-4f;
    static double r0 = .422784335098467f;
    static double r1 = .848044614534529f;
    static double r2 = .565221050691933f;
    static double r3 = .156513060486551f;
    static double r4 = .017050248402265f;
    static double r5 = 4.97958207639485e-4f;
    static double s1 = 1.24313399877507f;
    static double s2 = .548042109832463f;
    static double s3 = .10155218743983f;
    static double s4 = .00713309612391f;
    static double s5 = 1.16165475989616e-4f;

    /* System generated locals */
    double ret_val;

    /* Local variables */
    static double w, x;

/* ----------------------------------------------------------------------- */
/*     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25 */
/* ----------------------------------------------------------------------- */
/* ---------------------- */
/* ---------------------- */
    if (*a >= .6f) {
	goto L10;
    }
    w = ((((((p6 * *a + p5) * *a + p4) * *a + p3) * *a + p2) * *a + p1) * *a 
	    + p0) / ((((((q6 * *a + q5) * *a + q4) * *a + q3) * *a + q2) * *a 
	    + q1) * *a + 1.f);
    ret_val = -(*a) * w;
    return ret_val;

L10:
    x = *a - .5f - .5f;
    w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) / (((((s5 * 
	    x + s4) * x + s3) * x + s2) * x + s1) * x + 1.f);
    ret_val = x * w;
    return ret_val;
} /* gamln1_ */

double psi_(double *xx)
{
    /* Initialized data */

    static double piov4 = .785398163397448f;
    static double dx0 = 1.461632144968362341262659542325721325;
    static double p1[7] = { .0089538502298197f,4.77762828042627f,
	    142.441585084029f,1186.45200713425f,3633.51846806499f,
	    4138.10161269013f,1305.60269827897f };
    static double q1[6] = { 44.8452573429826f,520.752771467162f,
	    2210.0079924783f,3641.27349079381f,1908.310765963f,
	    6.91091682714533e-6f };
    static double p2[4] = { -2.12940445131011f,-7.01677227766759f,
	    -4.48616543918019f,-.648157123766197f };
    static double q2[4] = { 32.2703493791143f,89.2920700481861f,
	    54.6117738103215f,7.77788548522962f };

    /* System generated locals */
    double ret_val, r__1, r__2;

    /* Builtin functions */
    double cos(double), sin(double), log(double);

    /* Local variables */
    static long int i__, m, n;
    static double w, x, z__;
    static long int nq;
    static double den, aug, sgn, xmx0, xmax1, upper;
    static double xsmall;

/* --------------------------------------------------------------------- */

/*                 EVALUATION OF THE DIGAMMA FUNCTION */

/*                           ----------- */

/*     PSI(XX) IS ASSIGNED THE VALUE 0 WHEN THE DIGAMMA FUNCTION CANNOT */
/*     BE COMPUTED. */

/*     THE MAIN COMPUTATION INVOLVES EVALUATION OF RATIONAL CHEBYSHEV */
/*     APPROXIMATIONS PUBLISHED IN MATH. COMP. 27, 123-127(1973) BY */
/*     CODY, STRECOK AND THACHER. */

/* --------------------------------------------------------------------- */
/*     PSI WAS WRITTEN AT ARGONNE NATIONAL LABORATORY FOR THE FUNPACK */
/*     PACKAGE OF SPECIAL FUNCTION SUBROUTINES. PSI WAS MODIFIED BY */
/*     A.H. MORRIS (NSWC). */
/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

/*     PIOV4 = PI/4 */
/*     DX0 = ZERO OF PSI TO EXTENDED PRECISION */

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF */
/*     PSI(X) / (X - X0),  0.5 .LE. X .LE. 3.0 */

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF */
/*     PSI(X) - LN(X) + 1 / (2*X),  X .GT. 3.0 */

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

/*     MACHINE DEPENDENT CONSTANTS ... */

/*        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT */
/*                 WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED */
/*                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE */
/*                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH */
/*                 PSI MAY BE REPRESENTED AS ALOG(X). */

/*        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X) */
/*                 MAY BE REPRESENTED BY 1/X. */

/* --------------------------------------------------------------------- */
    xmax1 = (double) INT_MAX;
/* Computing MIN */
    r__1 = xmax1, r__2 = 0.5 / (0.5 * DBL_EPSILON);
    xmax1 = DMIN(r__1,r__2);
    xsmall = 1e-9f;
/* --------------------------------------------------------------------- */
    x = *xx;
    aug = 0.f;
    if (x >= .5f) {
	goto L200;
    }
/* --------------------------------------------------------------------- */
/*     X .LT. 0.5,  USE REFLECTION FORMULA */
/*     PSI(1-X) = PSI(X) + PI * COTAN(PI*X) */
/* --------------------------------------------------------------------- */
    if (DABS(x) > xsmall) {
	goto L100;
    }
    if (x == 0.f) {
	goto L400;
    }
/* --------------------------------------------------------------------- */
/*     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE */
/*     FOR  PI*COTAN(PI*X) */
/* --------------------------------------------------------------------- */
    aug = -1.f / x;
    goto L150;
/* --------------------------------------------------------------------- */
/*     REDUCTION OF ARGUMENT FOR COTAN */
/* --------------------------------------------------------------------- */
L100:
    w = -x;
    sgn = piov4;
    if (w > 0.f) {
	goto L120;
    }
    w = -w;
    sgn = -sgn;
/* --------------------------------------------------------------------- */
/*     MAKE AN ERROR EXIT IF X .LE. -XMAX1 */
/* --------------------------------------------------------------------- */
L120:
    if (w >= xmax1) {
	goto L400;
    }
    nq = (long int) w;
    w -= (double) nq;
    nq = (long int) (w * 4.f);
    w = (w - (double) nq * .25f) * 4.f;
/* --------------------------------------------------------------------- */
/*     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X. */
/*     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST */
/*     QUADRANT AND DETERMINE SIGN */
/* --------------------------------------------------------------------- */
    n = nq / 2;
    if (n + n != nq) {
	w = 1.f - w;
    }
    z__ = piov4 * w;
    m = n / 2;
    if (m + m != n) {
	sgn = -sgn;
    }
/* --------------------------------------------------------------------- */
/*     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X) */
/* --------------------------------------------------------------------- */
    n = (nq + 1) / 2;
    m = n / 2;
    m += m;
    if (m != n) {
	goto L140;
    }
/* --------------------------------------------------------------------- */
/*     CHECK FOR SINGULARITY */
/* --------------------------------------------------------------------- */
    if (z__ == 0.f) {
	goto L400;
    }
/* --------------------------------------------------------------------- */
/*     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND */
/*     SIN/COS AS A SUBSTITUTE FOR TAN */
/* --------------------------------------------------------------------- */
    aug = sgn * (cos(z__) / sin(z__) * 4.f);
    goto L150;
L140:
    aug = sgn * (sin(z__) / cos(z__) * 4.f);
L150:
    x = 1.f - x;
L200:
    if (x > 3.f) {
	goto L300;
    }
/* --------------------------------------------------------------------- */
/*     0.5 .LE. X .LE. 3.0 */
/* --------------------------------------------------------------------- */
    den = x;
    upper = p1[0] * x;

    for (i__ = 1; i__ <= 5; ++i__) {
	den = (den + q1[i__ - 1]) * x;
	upper = (upper + p1[i__]) * x;
/* L210: */
    }

    den = (upper + p1[6]) / (den + q1[5]);
    xmx0 = (double) x - dx0;
    ret_val = den * xmx0 + aug;
    return ret_val;
/* --------------------------------------------------------------------- */
/*     IF X .GE. XMAX1, PSI = LN(X) */
/* --------------------------------------------------------------------- */
L300:
    if (x >= xmax1) {
	goto L350;
    }
/* --------------------------------------------------------------------- */
/*     3.0 .LT. X .LT. XMAX1 */
/* --------------------------------------------------------------------- */
    w = 1.f / (x * x);
    den = w;
    upper = p2[0] * w;

    for (i__ = 1; i__ <= 3; ++i__) {
	den = (den + q2[i__ - 1]) * w;
	upper = (upper + p2[i__]) * w;
/* L310: */
    }

    aug = upper / (den + q2[3]) - .5f / x + aug;
L350:
    ret_val = aug + log(x);
    return ret_val;
/* --------------------------------------------------------------------- */
/*     ERROR RETURN */
/* --------------------------------------------------------------------- */
L400:
    ret_val = 0.f;
    return ret_val;
} /* psi_ */

double betaln_(double *a0, double *b0)
{
    /* Initialized data */

    static double e = .918938533204673f;

    /* System generated locals */
    long int i__1;
    double ret_val, r__1;

    /* Builtin functions */
    double log(double);

    /* Local variables */
    static double a, b, c__, h__;
    static long int i__, n;
    static double u, v, w, z__;
    extern double gamln_(double *), bcorr_(double *, double *), algdiv_(double *, 
	    double *), alnrel_(double *), gsumln_(double *, double *);

/* ----------------------------------------------------------------------- */
/*     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION */
/* ----------------------------------------------------------------------- */
/*     E = 0.5*LN(2*PI) */
/* -------------------------- */
/* -------------------------- */
    a = DMIN(*a0,*b0);
    b = DMAX(*a0,*b0);
    if (a >= 8.f) {
	goto L60;
    }
    if (a >= 1.f) {
	goto L20;
    }
/* ----------------------------------------------------------------------- */
/*                   PROCEDURE WHEN A .LT. 1 */
/* ----------------------------------------------------------------------- */
    if (b >= 8.f) {
	goto L10;
    }
    r__1 = a + b;
    ret_val = gamln_(&a) + (gamln_(&b) - gamln_(&r__1));
    return ret_val;
L10:
    ret_val = gamln_(&a) + algdiv_(&a, &b);
    return ret_val;
/* ----------------------------------------------------------------------- */
/*                PROCEDURE WHEN 1 .LE. A .LT. 8 */
/* ----------------------------------------------------------------------- */
L20:
    if (a > 2.f) {
	goto L30;
    }
    if (b > 2.f) {
	goto L21;
    }
    ret_val = gamln_(&a) + gamln_(&b) - gsumln_(&a, &b);
    return ret_val;
L21:
    w = 0.f;
    if (b < 8.f) {
	goto L40;
    }
    ret_val = gamln_(&a) + algdiv_(&a, &b);
    return ret_val;

/*                REDUCTION OF A WHEN B .LE. 1000 */

L30:
    if (b > 1e3f) {
	goto L50;
    }
    n = a - 1.f;
    w = 1.f;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a += -1.f;
	h__ = a / b;
	w *= h__ / (h__ + 1.f);
/* L31: */
    }
    w = log(w);
    if (b < 8.f) {
	goto L40;
    }
    ret_val = w + gamln_(&a) + algdiv_(&a, &b);
    return ret_val;

/*                 REDUCTION OF B WHEN B .LT. 8 */

L40:
    n = b - 1.f;
    z__ = 1.f;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b += -1.f;
	z__ *= b / (a + b);
/* L41: */
    }
    ret_val = w + log(z__) + (gamln_(&a) + (gamln_(&b) - gsumln_(&a, &b)));
    return ret_val;

/*                REDUCTION OF A WHEN B .GT. 1000 */

L50:
    n = a - 1.f;
    w = 1.f;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a += -1.f;
	w *= a / (a / b + 1.f);
/* L51: */
    }
    ret_val = log(w) - n * log(b) + (gamln_(&a) + algdiv_(&a, &b));
    return ret_val;
/* ----------------------------------------------------------------------- */
/*                   PROCEDURE WHEN A .GE. 8 */
/* ----------------------------------------------------------------------- */
L60:
    w = bcorr_(&a, &b);
    h__ = a / b;
    c__ = h__ / (h__ + 1.f);
    u = -(a - .5f) * log(c__);
    v = b * alnrel_(&h__);
    if (u <= v) {
	goto L61;
    }
    ret_val = log(b) * -.5f + e + w - v - u;
    return ret_val;
L61:
    ret_val = log(b) * -.5f + e + w - u - v;
    return ret_val;
} /* betaln_ */

double gsumln_(double *a, double *b)
{
    /* System generated locals */
    double ret_val, r__1;

    /* Builtin functions */
    double log(double);

    /* Local variables */
    static double x;
    extern double gamln1_(double *), alnrel_(double *);

/* ----------------------------------------------------------------------- */
/*          EVALUATION OF THE FUNCTION LN(GAMMA(A + B)) */
/*          FOR 1 .LE. A .LE. 2  AND  1 .LE. B .LE. 2 */
/* ----------------------------------------------------------------------- */
    x = (double) (*a) + (double) (*b) - 2.;
    if (x > .25f) {
	goto L10;
    }
    r__1 = x + 1.f;
    ret_val = gamln1_(&r__1);
    return ret_val;
L10:
    if (x > 1.25f) {
	goto L20;
    }
    ret_val = gamln1_(&x) + alnrel_(&x);
    return ret_val;
L20:
    r__1 = x - 1.f;
    ret_val = gamln1_(&r__1) + log(x * (x + 1.f));
    return ret_val;
} /* gsumln_ */

double bcorr_(double *a0, double *b0)
{
    /* Initialized data */

    static double c0 = .0833333333333333f;
    static double c1 = -.00277777777760991f;
    static double c2 = 7.9365066682539e-4f;
    static double c3 = -5.9520293135187e-4f;
    static double c4 = 8.37308034031215e-4f;
    static double c5 = -.00165322962780713f;

    /* System generated locals */
    double ret_val, r__1;

    /* Local variables */
    static double a, b, c__, h__, t, w, x, s3, s5, x2, s7, s9, s11;

/* ----------------------------------------------------------------------- */

/*     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE */
/*     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A). */
/*     IT IS ASSUMED THAT A0 .GE. 8 AND B0 .GE. 8. */

/* ----------------------------------------------------------------------- */
/* ------------------------ */
    a = DMIN(*a0,*b0);
    b = DMAX(*a0,*b0);

    h__ = a / b;
    c__ = h__ / (h__ + 1.f);
    x = 1.f / (h__ + 1.f);
    x2 = x * x;

/*                SET SN = (1 - X**N)/(1 - X) */

    s3 = x + x2 + 1.f;
    s5 = x + x2 * s3 + 1.f;
    s7 = x + x2 * s5 + 1.f;
    s9 = x + x2 * s7 + 1.f;
    s11 = x + x2 * s9 + 1.f;

/*                SET W = DEL(B) - DEL(A + B) */

/* Computing 2nd power */
    r__1 = 1.f / b;
    t = r__1 * r__1;
    w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 * 
	    s3) * t + c0;
    w *= c__ / b;

/*                   COMPUTE  DEL(A) + W */

/* Computing 2nd power */
    r__1 = 1.f / a;
    t = r__1 * r__1;
    ret_val = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a + 
	    w;
    return ret_val;
} /* bcorr_ */

double algdiv_(double *a, double *b)
{
    /* Initialized data */

    static double c0 = .0833333333333333f;
    static double c1 = -.00277777777760991f;
    static double c2 = 7.9365066682539e-4f;
    static double c3 = -5.9520293135187e-4f;
    static double c4 = 8.37308034031215e-4f;
    static double c5 = -.00165322962780713f;

    /* System generated locals */
    double ret_val, r__1;

    /* Builtin functions */
    double log(double);

    /* Local variables */
    static double c__, d__, h__, t, u, v, w, x, s3, s5, x2, s7, s9, s11;
    extern double alnrel_(double *);

/* ----------------------------------------------------------------------- */

/*     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B .GE. 8 */

/*                         -------- */

/*     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY */
/*     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X). */

/* ----------------------------------------------------------------------- */
/* ------------------------ */
    if (*a <= *b) {
	goto L10;
    }
    h__ = *b / *a;
    c__ = 1.f / (h__ + 1.f);
    x = h__ / (h__ + 1.f);
    d__ = *a + (*b - .5f);
    goto L20;
L10:
    h__ = *a / *b;
    c__ = h__ / (h__ + 1.f);
    x = 1.f / (h__ + 1.f);
    d__ = *b + (*a - .5f);

/*                SET SN = (1 - X**N)/(1 - X) */

L20:
    x2 = x * x;
    s3 = x + x2 + 1.f;
    s5 = x + x2 * s3 + 1.f;
    s7 = x + x2 * s5 + 1.f;
    s9 = x + x2 * s7 + 1.f;
    s11 = x + x2 * s9 + 1.f;

/*                SET W = DEL(B) - DEL(A + B) */

/* Computing 2nd power */
    r__1 = 1.f / *b;
    t = r__1 * r__1;
    w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 * 
	    s3) * t + c0;
    w *= c__ / *b;

/*                    COMBINE THE RESULTS */

    r__1 = *a / *b;
    u = d__ * alnrel_(&r__1);
    v = *a * (log(*b) - 1.f);
    if (u <= v) {
	goto L30;
    }
    ret_val = w - v - u;
    return ret_val;
L30:
    ret_val = w - u - v;
    return ret_val;
} /* algdiv_ */

double gamln_(double *a)
{
    /* Initialized data */

    static double d__ = .418938533204673f;
    static double c0 = .0833333333333333f;
    static double c1 = -.00277777777760991f;
    static double c2 = 7.9365066682539e-4f;
    static double c3 = -5.9520293135187e-4f;
    static double c4 = 8.37308034031215e-4f;
    static double c5 = -.00165322962780713f;

    /* System generated locals */
    long int i__1;
    double ret_val, r__1;

    /* Builtin functions */
    double log(double);

    /* Local variables */
    static long int i__, n;
    static double t, w;
    extern double gamln1_(double *);

/* ----------------------------------------------------------------------- */
/*            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A */
/* ----------------------------------------------------------------------- */
/*     WRITTEN BY ALFRED H. MORRIS */
/*          NAVAL SURFACE WARFARE CENTER */
/*          DAHLGREN, VIRGINIA */
/* -------------------------- */
/*     D = 0.5*(LN(2*PI) - 1) */
/* -------------------------- */
/* -------------------------- */
/* ----------------------------------------------------------------------- */
    if (*a > .8f) {
	goto L10;
    }
    ret_val = gamln1_(a) - log(*a);
    return ret_val;
L10:
    if (*a > 2.25f) {
	goto L20;
    }
    t = *a - .5f - .5f;
    ret_val = gamln1_(&t);
    return ret_val;

L20:
    if (*a >= 10.f) {
	goto L30;
    }
    n = *a - 1.25f;
    t = *a;
    w = 1.f;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t += -1.f;
/* L21: */
	w = t * w;
    }
    r__1 = t - 1.f;
    ret_val = gamln1_(&r__1) + log(w);
    return ret_val;

L30:
/* Computing 2nd power */
    r__1 = 1.f / *a;
    t = r__1 * r__1;
    w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / *a;
    ret_val = d__ + w + (*a - .5f) * (log(*a) - 1.f);
    return ret_val;
} /* gamln_ */

