//algorithm 554

//this function is from netlib/toms, see readme.txt for license
//I convertied it to C++ by means of f2c and made some minor chages
//Evgenii Rudnyi rudnyi@comp.chem.msu.su

/* brent1.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "toms.h"
#include "f2c.h"
#include "f2clib.h"

int brentm_(const VOF &fcn, integer *n, double *x, double *fvec, double
            *ftol, double *xtol, integer *maxfev, integer *mopt, integer *info,
            integer *nfev, double *q, integer *ldq, double *sigma, double
            *wa1, double *wa2);


/* Subroutine */ int brent1_(const VOF &fcn, integer *n, doublereal *x, doublereal *
  fvec, doublereal *tol, integer *info, doublereal *wa, integer *lwa)
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    doublereal emax;
    integer nfev;
    doublereal ftol, temp;
    integer mopt;
    doublereal xtol;
    integer i__, maxfev;

/*     ********** */

/*     SUBROUTINE BRENT1 */

/*     THE PURPOSE OF THIS SUBROUTINE IS TO FIND A ZERO OF */
/*     A SYSTEM OF N NONLINEAR EQUATIONS IN N VARIABLES BY A */
/*     METHOD DUE TO R. BRENT. THIS IS DONE BY USING THE */
/*     MORE GENERAL NONLINEAR EQUATION SOLVER BRENTM. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE BRENT1(FCN,N,X,FVEC,TOL,INFO,WA,LWA) */

/*     WHERE */

/*       FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH */
/*         CALCULATES COMPONENTS OF THE FUNCTION. FCN SHOULD BE */
/*         DECLARED IN AN EXTERNAL STATEMENT IN THE USER CALLING */
/*         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS. */

/*         SUBROUTINE FCN(N,X,FVEC,IFLAG) */
/*         INTEGER N,IFLAG */
/*         DOUBLE PRECISION X(N),FVEC(N) */
/*         ---------- */
/*         CALCULATE THE IFLAG-TH COMPONENT OF THE FUNCTION */
/*         AND RETURN THIS VALUE IN FVEC(IFLAG). */
/*         ---------- */
/*         RETURN */
/*         END */

/*         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS */
/*         THE USER WANTS TO TERMINATE EXECUTION OF BRENT1. */
/*         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER. */

/*       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF EQUATIONS AND VARIABLES. */

/*       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN */
/*         AN INITIAL ESTIMATE OF THE SOLUTION VECTOR. ON OUTPUT X */
/*         CONTAINS THE FINAL ESTIMATE OF THE SOLUTION VECTOR. */

/*       FVEC IS AN ARRAY OF LENGTH N. ON OUTPUT IT CONTAINS */
/*         THE FINAL RESIDUALS. */

/*       TOL IS A NONNEGATIVE INPUT VARIABLE. THE ALGORITHM CONVERGES */
/*         IF EITHER ALL THE RESIDUALS ARE AT MOST TOL IN MAGNITUDE, */
/*         OR IF THE ALGORITHM ESTIMATES THAT THE RELATIVE ERROR */
/*         BETWEEN X AND THE SOLUTION IS AT MOST TOL. */

/*       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS. IF */
/*         THE USER HAS TERMINATED EXECUTION, INFO WILL BE SET TO */
/*         THE (NEGATIVE) VALUE OF IFLAG. SEE DESCRIPTION OF FCN. */
/*         OTHERWISE */

/*         INFO = 0   IMPROPER INPUT PARAMETERS. */

/*         INFO = 1   ALL RESIDUALS ARE AT MOST TOL IN MAGNITUDE. */

/*         INFO = 2   ALGORITHM ESTIMATES THAT THE RELATIVE ERROR */
/*                    BETWEEN X AND THE SOLUTION IS AT MOST TOL. */

/*         INFO = 3   CONDITIONS FOR INFO = 1 AND INFO = 2 BOTH HOLD. */

/*         INFO = 4   NUMBER OF FUNCTION EVALUATIONS HAS REACHED OR */
/*                    EXCEEDED 50*(N+3). */

/*         INFO = 5   APPROXIMATE JACOBIAN MATRIX IS SINGULAR. */

/*         INFO = 6   ITERATION IS NOT MAKING GOOD PROGRESS. */

/*         INFO = 7   ITERATION IS DIVERGING. */

/*         INFO = 8   ITERATION IS CONVERGING, BUT TOL IS TOO */
/*                    SMALL, OR THE CONVERGENCE IS VERY SLOW */
/*                    DUE TO A JACOBIAN SINGULAR NEAR THE OUTPUT */
/*                    X OR DUE TO BADLY SCALED VARIABLES. */

/*       WA IS A WORK ARRAY OF LENGTH LWA. */

/*       LWA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN */
/*         N*(N+3). */

/*     SUBPROGRAMS REQUIRED */

/*       USER-SUPPLIED ...... FCN, BRENTM */

/*       FORTRAN-SUPPLIED ... DLOG */

/*     ********** */
    /* Parameter adjustments */
    --fvec;
    --x;
    --wa;

    /* Function Body */
    *info = 0;

/*     CHECK THE INPUT PARAMETERS FOR ERRORS. */

    if (*n <= 0 || *tol < zero || *lwa < *n * (*n + 3)) {
  goto L30;
    }

/*     DETERMINE AN OPTIMAL VALUE FOR MOPT. */

    emax = zero;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
  i__2 = i__ + 1;
  i__3 = *n + (i__ << 1) + 1;
  temp = log((doublereal) i__2) / (doublereal) i__3;
  if (temp < emax) {
      goto L20;
  }
  mopt = i__;
  emax = temp;
/* L10: */
    }
L20:

/*     CALL BRENTM. */

    maxfev = (*n + 3) * 50;
    ftol = *tol;
    xtol = *tol;
    brentm_(fcn, n, &x[1], &fvec[1], &ftol, &xtol, &maxfev, &mopt, info,
       &nfev, &wa[*n * 3 + 1], n, &wa[1], &wa[*n + 1], &wa[(*n << 1) + 
      1]);
L30:
    return 0;

/*     LAST CARD OF SUBROUTINE BRENT1. */

} /* brent1_ */

/* Subroutine */ int brentm_(const VOF &fcn, integer *n, doublereal *x, doublereal *
  fvec, doublereal *ftol, doublereal *xtol, integer *maxfev, integer *
  mopt, integer *info, integer *nfev, doublereal *q, integer *ldq, 
  doublereal *sigma, doublereal *wa1, doublereal *wa2)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal p05 = .05;
    static doublereal scale = 10.;

    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    logical conv;
    doublereal temp;
    integer nier6, nier7, nier8;
    doublereal h__;
    integer i__, j, k, m, iflag;
    doublereal delta, difit;
    integer nsing;
    doublereal fnorm, xnorm, difit1, fnorm1;
    integer nfcall;
    doublereal epsmch, sknorm, eta, eps, fky, fkz;

/*     ********** */

/*     SUBROUTINE BRENTM */

/*     THE PURPOSE OF THIS SUBROUTINE IS TO FIND A ZERO TO */
/*     A SYSTEM OF N NONLINEAR EQUATIONS IN N VARIABLES BY A */
/*     METHOD DUE TO R. BRENT. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE BRENTM(FCN,N,X,FVEC,FTOL,XTOL,MAXFEV,MOPT, */
/*                         INFO,NFEV,Q,LDQ,SIGMA,WA1,WA2) */

/*     WHERE */

/*       FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH */
/*         CALCULATES COMPONENTS OF THE FUNCTION. FCN SHOULD BE */
/*         DECLARED IN AN EXTERNAL STATEMENT IN THE USER CALLING */
/*         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS. */

/*         SUBROUTINE FCN(N,X,FVEC,IFLAG) */
/*         INTEGER N,IFLAG */
/*         DOUBLE PRECISION X(N),FVEC(N) */
/*         ---------- */
/*         CALCULATE THE IFLAG-TH COMPONENT OF THE FUNCTION */
/*         AND RETURN THIS VALUE IN FVEC(IFLAG). */
/*         ---------- */
/*         RETURN */
/*         END */

/*         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS */
/*         THE USER WANTS TO TERMINATE EXECUTION OF BRENTM. */
/*         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER. */

/*       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF */
/*         EQUATIONS AND VARIABLES. */

/*       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN */
/*         AN ESTIMATE TO THE SOLUTION OF THE SYSTEM OF EQUATIONS. */
/*         ON OUTPUT X CONTAINS THE FINAL ESTIMATE TO THE SOLUTION */
/*         OF THE SYSTEM OF EQUATIONS. */

/*       FVEC IS AN ARRAY OF LENGTH N. ON OUTPUT IT CONTAINS */
/*         THE FINAL RESIDUALS. */

/*       FTOL IS A NONNEGATIVE INPUT VARIABLE. CONVERGENCE */
/*         OCCURS IF ALL RESIDUALS ARE AT MOST FTOL IN MAGNITUDE. */

/*       XTOL IS A NONNEGATIVE INPUT VARIABLE. CONVERGENCE */
/*         OCCURS IF THE RELATIVE ERROR BETWEEN TWO SUCCESSIVE */
/*         ITERATES IS AT MOST XTOL. */

/*       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION */
/*         OCCURS IF THE NUMBER OF FUNCTION EVALUATIONS IS AT */
/*         LEAST MAXFEV BY THE END OF AN ITERATION. IN BRENTM, */
/*         A FUNCTION EVALUATION CORRESPONDS TO N CALLS TO FCN. */

/*       MOPT IS A POSITIVE INTEGER INPUT VARIABLE. MOPT SPECIFIES */
/*         THE NUMBER OF TIMES THAT THE APPROXIMATE JACOBIAN IS */
/*         USED DURING EACH ITERATION WHICH EMPLOYS ITERATIVE */
/*         REFINEMENT. IF MOPT IS 1, NO ITERATIVE REFINEMENT WILL */
/*         BE DONE. MAXIMUM EFFICIENCY IS USUALLY OBTAINED IF */
/*         MOPT MAXIMIZES LOG(K+1)/(N+2*K+1) FOR K = 1,...,N. */

/*       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS. IF */
/*         THE USER HAS TERMINATED EXECUTION, INFO WILL BE SET TO */
/*         THE (NEGATIVE) VALUE OF IFLAG. SEE DESCRIPTION OF FCN. */
/*         OTHERWISE */

/*         INFO = 0   IMPROPER INPUT PARAMETERS. */

/*         INFO = 1   ALL RESIDUALS ARE AT MOST FTOL IN MAGNITUDE. */

/*         INFO = 2   RELATIVE ERROR BETWEEN TWO SUCCESSIVE ITERATES */
/*                    IS AT MOST XTOL. */

/*         INFO = 3   CONDITIONS FOR INFO = 1 AND INFO = 2 BOTH HOLD. */

/*         INFO = 4   NUMBER OF FUNCTION EVALUATIONS HAS REACHED OR */
/*                    EXCEEDED MAXFEV. */

/*         INFO = 5   APPROXIMATE JACOBIAN MATRIX IS SINGULAR. */

/*         INFO = 6   ITERATION IS NOT MAKING GOOD PROGRESS. */

/*         INFO = 7   ITERATION IS DIVERGING. */

/*         INFO = 8   ITERATION IS CONVERGING, BUT XTOL IS TOO */
/*                    SMALL, OR THE CONVERGENCE IS VERY SLOW */
/*                    DUE TO A JACOBIAN SINGULAR NEAR THE OUTPUT */
/*                    X OR DUE TO BADLY SCALED VARIABLES. */

/*       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF */
/*         FUNCTION EVALUATIONS USED IN PRODUCING X. IN BRENTM, */
/*         A FUNCTION EVALUATION CORRESPONDS TO N CALLS TO FCN. */

/*       Q IS AN N BY N ARRAY. IF JAC DENOTES THE APPROXIMATE */
/*         JACOBIAN, THEN ON OUTPUT Q IS (A MULTIPLE OF) AN */
/*         ORTHOGONAL MATRIX SUCH THAT JAC*Q IS A LOWER TRIANGULAR */
/*         MATRIX. ONLY THE DIAGONAL ELEMENTS OF JAC*Q NEED */
/*         TO BE STORED, AND THESE CAN BE FOUND IN SIGMA. */

/*       LDQ IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N */
/*         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY Q. */

/*       SIGMA IS A LINEAR ARRAY OF LENGTH N. ON OUTPUT SIGMA */
/*         CONTAINS THE DIAGONAL ELEMENTS OF THE MATRIX JAC*Q. */
/*         SEE DESCRIPTION OF Q. */

/*       WA1 AND WA2 ARE LINEAR WORK ARRAYS OF LENGTH N. */

/*     SUBPROGRAMS REQUIRED */

/*       USER-SUPPLIED ...... FCN */

/*       FORTRAN-SUPPLIED ... DABS,DMAX1,DSQRT,DSIGN */

/*     ********** */
    /* Parameter adjustments */
    --wa2;
    --wa1;
    --sigma;
    --fvec;
    --x;
    q_dim1 = *ldq;
    q_offset = q_dim1 + 1;
    q -= q_offset;

    /* Function Body */

/*     WARNING. */

/*     THIS IS AN IBM CODE. TO RUN THIS CODE ON OTHER MACHINES IT */
/*     IS NECESSARY TO CHANGE THE VALUE OF THE MACHINE PRECISION */
/*     EPSMCH. THE MACHINE PRECISION IS THE SMALLEST FLOATING */
/*     POINT NUMBER EPSMCH SUCH THAT */

/*           1 + EPSMCH .GT. 1 */

/*     IN WORKING PRECISION. IF IN DOUBT ABOUT THE VALUE OF */
/*     EPSMCH, THEN THE FOLLOWING PROGRAM SEGMENT DETERMINES */
/*     EPSMCH ON MOST MACHINES. */

/*     EPSMCH = 0.5D0 */
/*   1 CONTINUE */
/*     IF (1.D0+EPSMCH .EQ. 1.D0) GO TO 2 */
/*     EPSMCH = 0.5D0*EPSMCH */
/*     GO TO 1 */
/*   2 CONTINUE */
/*     EPSMCH = 2.D0*EPSMCH */

/*     THE IBM DOUBLE PRECISION EPSMCH. */

    epsmch = DBL_EPSILON*100.;
    *info = 0;
    iflag = 0;
    *nfev = 0;
    nfcall = 0;

/*     CHECK THE INPUT PARAMETERS FOR ERRORS. */

    if (*n <= 0 || *ftol < zero || *xtol < zero || *maxfev <= 0 || *mopt <= 0 
      || *ldq < *n) {
  goto L220;
    }

/*     INITIALIZE SOME OF THE VARIABLES. */

    nier6 = -1;
    nier7 = -1;
    nier8 = 0;
    fnorm = zero;
    difit = zero;
    xnorm = zero;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
  d__2 = xnorm, d__3 = (d__1 = x[i__], abs(d__1));
  xnorm = max(d__2,d__3);
/* L10: */
    }
    eps = sqrt(epsmch);
    delta = scale * xnorm;
    if (xnorm == zero) {
  delta = scale;
    }

/*     ENTER THE PRINCIPAL ITERATION. */

L20:

/*     TO PRINT THE ITERATES, PLACE WRITE STATEMENTS */
/*     FOR THE VECTOR X HERE. */

    nsing = *n;
    fnorm1 = fnorm;
    difit1 = difit;
    fnorm = zero;

/*     COMPUTE THE STEP H FOR THE DIVIDED DIFFERENCE WHICH */
/*     APPROXIMATES THE K-TH ROW OF THE JACOBIAN MATRIX. */

    h__ = eps * xnorm;
    if (h__ == zero) {
  h__ = eps;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
  i__2 = *n;
  for (i__ = 1; i__ <= i__2; ++i__) {
      q[i__ + j * q_dim1] = zero;
/* L30: */
  }
  q[j + j * q_dim1] = h__;
  wa1[j] = x[j];
/* L40: */
    }

/*     ENTER A SUBITERATION. */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
  iflag = k;
        fcn(n, &wa1[1], &fvec[1], &iflag);
  fky = fvec[k];
  ++nfcall;
  *nfev = nfcall / *n;
  if (iflag < 0) {
      goto L230;
  }
/* Computing MAX */
  d__1 = fnorm, d__2 = abs(fky);
  fnorm = max(d__1,d__2);

/*        COMPUTE THE K-TH ROW OF THE JACOBIAN MATRIX. */

  i__2 = *n;
  for (j = k; j <= i__2; ++j) {
      i__3 = *n;
      for (i__ = 1; i__ <= i__3; ++i__) {
    wa2[i__] = wa1[i__] + q[i__ + j * q_dim1];
/* L50: */
      }
            fcn(n, &wa2[1], &fvec[1], &iflag);
      fkz = fvec[k];
      ++nfcall;
      *nfev = nfcall / *n;
      if (iflag < 0) {
    goto L230;
      }
      sigma[j] = fkz - fky;
/* L60: */
  }
  fvec[k] = fky;

/*        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE K-TH RO
W */
/*        OF THE JACOBIAN MATRIX TO A MULTIPLE OF THE K-TH UNIT VECTOR
. */

  eta = zero;
  i__2 = *n;
  for (i__ = k; i__ <= i__2; ++i__) {
/* Computing MAX */
      d__2 = eta, d__3 = (d__1 = sigma[i__], abs(d__1));
      eta = max(d__2,d__3);
/* L70: */
  }
  if (eta == zero) {
      goto L150;
  }
  --nsing;
  sknorm = zero;
  i__2 = *n;
  for (i__ = k; i__ <= i__2; ++i__) {
      sigma[i__] /= eta;
/* Computing 2nd power */
      d__1 = sigma[i__];
      sknorm += d__1 * d__1;
/* L80: */
  }
  sknorm = sqrt(sknorm);
  if (sigma[k] < zero) {
      sknorm = -sknorm;
  }
  sigma[k] += sknorm;

/*        APPLY THE TRANSFORMATION AND COMPUTE THE MATRIX Q. */

  i__2 = *n;
  for (i__ = 1; i__ <= i__2; ++i__) {
      wa2[i__] = zero;
/* L90: */
  }
  i__2 = *n;
  for (j = k; j <= i__2; ++j) {
      temp = sigma[j];
      i__3 = *n;
      for (i__ = 1; i__ <= i__3; ++i__) {
    wa2[i__] += temp * q[i__ + j * q_dim1];
/* L100: */
      }
/* L110: */
  }
  i__2 = *n;
  for (j = k; j <= i__2; ++j) {
      temp = sigma[j] / (sknorm * sigma[k]);
      i__3 = *n;
      for (i__ = 1; i__ <= i__3; ++i__) {
    q[i__ + j * q_dim1] -= temp * wa2[i__];
/* L120: */
      }
/* L130: */
  }

/*        COMPUTE THE SUBITERATE. */

  sigma[k] = sknorm * eta;
  temp = fky / sigma[k];
  if (h__ * abs(temp) > delta) {
      d__1 = delta / h__;
      temp = d_sign(&d__1, &temp);
  }
  i__2 = *n;
  for (i__ = 1; i__ <= i__2; ++i__) {
      wa1[i__] += temp * q[i__ + k * q_dim1];
/* L140: */
  }
L150:
  ;
    }

/*     COMPUTE THE NORMS OF THE ITERATE AND CORRECTION VECTOR. */

    xnorm = zero;
    difit = zero;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
  d__2 = xnorm, d__3 = (d__1 = wa1[i__], abs(d__1));
  xnorm = max(d__2,d__3);
/* Computing MAX */
  d__2 = difit, d__3 = (d__1 = x[i__] - wa1[i__], abs(d__1));
  difit = max(d__2,d__3);
  x[i__] = wa1[i__];
/* L160: */
    }

/*     UPDATE THE BOUND ON THE CORRECTION VECTOR. */

/* Computing MAX */
    d__1 = delta, d__2 = scale * xnorm;
    delta = max(d__1,d__2);

/*     DETERMINE THE PROGRESS OF THE ITERATION. */

    conv = fnorm < fnorm1 && difit < difit1 && nsing == 0;
    ++nier6;
    ++nier7;
    ++nier8;
    if (conv) {
  nier6 = 0;
    }
    if (fnorm < fnorm1 || difit < difit1) {
  nier7 = 0;
    }
    if (difit > eps * xnorm) {
  nier8 = 0;
    }

/*     TESTS FOR CONVERGENCE. */

    if (fnorm <= *ftol) {
  *info = 1;
    }
    if (difit <= *xtol * xnorm && conv) {
  *info = 2;
    }
    if (fnorm <= *ftol && *info == 2) {
  *info = 3;
    }
    if (*info != 0) {
  goto L230;
    }

/*     TESTS FOR TERMINATION. */

    if (*nfev >= *maxfev) {
  *info = 4;
    }
    if (nsing == *n) {
  *info = 5;
    }
    if (nier6 == 5) {
  *info = 6;
    }
    if (nier7 == 3) {
  *info = 7;
    }
    if (nier8 == 4) {
  *info = 8;
    }
    if (*info != 0) {
  goto L230;
    }

/*     ITERATIVE REFINEMENT IS USED IF THE ITERATION IS CONVERGING. */

    if (! conv || difit > p05 * xnorm || *mopt == 1) {
  goto L220;
    }

/*     START ITERATIVE REFINEMENT. */

    i__1 = *mopt;
    for (m = 2; m <= i__1; ++m) {
  fnorm1 = fnorm;
  fnorm = zero;
  i__2 = *n;
  for (k = 1; k <= i__2; ++k) {
      iflag = k;
            fcn(n, &wa1[1], &fvec[1], &iflag);
      fky = fvec[k];
      ++nfcall;
      *nfev = nfcall / *n;
      if (iflag < 0) {
    goto L230;
      }
/* Computing MAX */
      d__1 = fnorm, d__2 = abs(fky);
      fnorm = max(d__1,d__2);

/*           ITERATIVE REFINEMENT IS TERMINATED IF IT DOES NOT */
/*           GIVE A REDUCTION OF THE RESIDUALS. */

      if (fnorm < fnorm1) {
    goto L170;
      }
      fnorm = fnorm1;
      goto L220;
L170:
      temp = fky / sigma[k];
      i__3 = *n;
      for (i__ = 1; i__ <= i__3; ++i__) {
    wa1[i__] += temp * q[i__ + k * q_dim1];
/* L180: */
      }
/* L190: */
  }

/*        COMPUTE THE NORMS OF THE ITERATE AND CORRECTION VECTOR. */

  xnorm = zero;
  difit = zero;
  i__2 = *n;
  for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
      d__2 = xnorm, d__3 = (d__1 = wa1[i__], abs(d__1));
      xnorm = max(d__2,d__3);
/* Computing MAX */
      d__2 = difit, d__3 = (d__1 = x[i__] - wa1[i__], abs(d__1));
      difit = max(d__2,d__3);
      x[i__] = wa1[i__];
/* L200: */
  }

/*        STOPPING CRITERIA FOR ITERATIVE REFINEMENT. */

  if (fnorm <= *ftol) {
      *info = 1;
  }
  if (difit <= *xtol * xnorm) {
      *info = 2;
  }
  if (fnorm <= *ftol && *info == 2) {
      *info = 3;
  }
  if (*nfev >= *maxfev && *info == 0) {
      *info = 4;
  }
  if (*info != 0) {
      goto L230;
  }
/* L210: */
    }
L220:

/*     END OF THE ITERATIVE REFINEMENT. */

    goto L20;

/*     TERMINATION, EITHER NORMAL OR USER IMPOSED. */

L230:
    if (iflag < 0) {
  *info = iflag;
    }
    return 0;

/*     LAST CARD OF SUBROUTINE BRENTM. */

} /* brentm_ */

