/* test.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "toms.h"
#include "g2c.h"

/* Common Block Declarations */

struct {
    integer nfcall;
} refnum_;

#define refnum_1 refnum_

/* Table of constant values */

static integer c__1 = 1;

/*    **********                                                        000000
10*/
/*                                                                      000000
20*/
/*    THIS IS A SAMPLE PROGRAM FOR THE EASY-TO-USE VERSION BRENT1       000000
30*/
/*    OF SUBROUTINE BRENTM. THIS PROGRAM SOLVES THE DISCRETE            000000
40*/
/*    BOUNDARY VALUE PROBLEM DEFINED BY THE SYSTEM OF NONLINEAR         000000
50*/
/*    EQUATIONS                                                         000000
60*/
/*                                                                      000000
70*/
/*    2*X(I) - X(I-1) - X(I+1)                                          000000
80*/
/*                                                                      000000
90*/
/*           + 0.5*(H**2)*(X(I) + I*H + 1)**3 = 0 , I = 1,...,N         000001
00*/
/*                                                                      000001
10*/
/*    WHERE H = 1/(N+1), AND X(0) = X(N+1) = 0.                         000001
20*/
/*                                                                      000001
30*/
/*    **********                                                        000001
40*/
/* Main program */
int main()
{
    /* Initialized data */

    static integer nwrite = 6;

    /* Format strings */
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */

    /* Local variables */
    doublereal fvec[10];
    integer info, nfev;
    doublereal temp, h__;
    integer i__, n;
    doublereal x[10];
    doublereal fnorm1, fnorm2, wa[130];
    integer lwa;
    extern /* Subroutine */ int bvp_(integer *, doublereal *, doublereal *, 
      integer *);
    doublereal tol;


/*                                                                      00
000210*/
/*    LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.                    00
000220*/
/*                                                                      00
000230*/
/*                                                                      00
000250*/
    lwa = 130;
    tol = 1e-10;
    n = 10;

/*     STARTING VALUES. */

    i__1 = n + 1;
    h__ = 1. / (doublereal) i__1;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
  temp = (doublereal) i__ * h__;
  x[i__ - 1] = temp * (temp - 1.);
/* L10: */
    }

/*     INITIAL MAX-NORM OF THE RESIDUALS. */

    fnorm1 = 0.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
  bvp_(&n, x, fvec, &i__);
/* Computing MAX */
  d__2 = fnorm1, d__3 = (d__1 = fvec[i__ - 1], abs(d__1));
  fnorm1 = max(d__2,d__3);
/* L20: */
    }

    refnum_1.nfcall = 0;
    brent1_(makeFunctor((VOF*)0, bvp_), &n, x, fvec, &tol, &info, wa, &lwa);
    nfev = refnum_1.nfcall / n;

/*     FINAL MAX-NORM OF THE RESIDUALS. */

    fnorm2 = 0.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
  bvp_(&n, x, fvec, &i__);
/* Computing MAX */
  d__2 = fnorm2, d__3 = (d__1 = fvec[i__ - 1], abs(d__1));
  fnorm2 = max(d__2,d__3);
/* L30: */
    }
    cout << "DIMENSION " << n << endl;
    cout << "INITIAL MAX-NORM OF THE RESIDUALS " << fnorm1 << endl;
    cout << "FINAL MAX-NORM OF THE RESIDUALS " << fnorm2 << endl;
    cout << "NUMBER OF FUNCTION EVALUATIONS " << nfev << endl;
    cout << "EXIT PARAMETER " << info << endl;
    cout << "FINAL APPROXIMATE SOLUTION" << endl;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        cout << x[i__ - 1] << endl;
    }

/*     LAST CARD OF SAMPLE PROGRAM. */

    return 0;
} /* MAIN__ */

/* Subroutine */ int bvp_(integer *n, doublereal *x, doublereal *fvec, 
  integer *iflag)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    doublereal temp, temp1, temp2, h__;

/*     ********** */

/*     SUBROUTINE BVP DEFINES THE BOUNDARY VALUE PROBLEM. */

/*     ********** */
    /* Parameter adjustments */
    --fvec;
    --x;

    /* Function Body */
    i__1 = *n + 1;
    h__ = 1. / (doublereal) i__1;
/* Computing 3rd power */
    d__1 = x[*iflag] + (doublereal) (*iflag) * h__ + 1., d__2 = d__1;
    temp = d__2 * (d__1 * d__1) * .5;
    temp1 = 0.;
    if (*iflag != 1) {
  temp1 = x[*iflag - 1];
    }
    temp2 = 0.;
    if (*iflag != *n) {
  temp2 = x[*iflag + 1];
    }
/* Computing 2nd power */
    d__1 = h__;
    fvec[*iflag] = x[*iflag] * 2. - temp1 - temp2 + temp * (d__1 * d__1);
    ++refnum_1.nfcall;
    return 0;

/*     LAST CARD OF SUBROUTINE BVP. */

} /* bvp_ */

