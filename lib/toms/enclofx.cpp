//algorithm 748

//this function is from netlib/toms, see readme.txt for license
//I convertied it to C++ by means of f2c and made some minor chages
//Evgenii Rudnyi rudnyi@comp.chem.msu.su

/* enclofx.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "toms.h"
#include "f2c.h"

int brackt_(const D2D &f, double *a, double *b, double *c__, double *fa,
            double *fb, double *tol, integer *neps, double *eps, double
            *d__, double *fd);

int isign_(double *x);

int tole_(double *b, double *tol, integer *neps, double *eps);

int newqua_(double *a, double *b, double *d__, double *fa, double *fb,
            double *fd, double *c__, integer *k);

int pzero_(double *a, double *b, double *d__, double *e, double *fa,
           double *fb, double *fd, double *fe, double *c__);

int rmp_(double *rel);

/* Table of constant values */

static integer c__2 = 2;
static integer c__3 = 3;

/* *** enclofx.f */
/* Subroutine */ int rroot_(const D2D &f, integer *neps, doublereal *eps,
  doublereal *a, doublereal *b, doublereal *root)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    doublereal prof, c__, d__, e, u, a0, b0;
    integer itnum;
    doublereal fa, fb, fd, fe, fu;
    doublereal tol;

/* FINDS EITHER AN EXACT SOLUTION OR AN APPROXIMATE SOLUTION OF THE */
/* EQUATION F(X)=0 IN THE INTERVAL [A,B]. AT THE BEGINING OF EACH */
/* ITERATION, THE CURRENT ENCLOSING INTERVAL IS RECORDED AS [A0,B0]. */
/* THE FIRST ITERATION IS SIMPLY A SECANT STEP. STARTING WITH THE */
/* SECOND ITERATION, THREE STEPS ARE TAKEN IN EACH ITERATION. FIRST */
/* TWO STEPS ARE EITHER QUADRATIC INTERPOLATION OR CUBIC INVERSE */
/* INTERPOLATION. THE THIRD STEP IS A DOUBLE-SIZE SECANT STEP. IF THE */
/* DIAMETER OF THE ENCLOSING INTERVAL OBTAINED AFTER THOSE THREE STEPS */
/* IS LARGER THAN 0.5*(B0-A0), THEN AN ADDITIONAL BISECTION STEP WILL */
/* BE TAKEN. */
/*  NPROB -- INTEGER. INDICATING THE PROBLEM TO BE SOLVED; */
/*  NEPS  -- INTEGER. USED TO DETERMINE THE TERMINATION CRITERION; */
/*  EPS   -- DOUBLE PRECISION. USED IN THE TERMINATION CRITERION; */
/*  A,B   -- DOUBLE PRECISION. INPUT AS THE INITIAL INTERVAL AND */
/*           OUTPUT AS THE ENCLOSING INTERVAL AT THE TERMINATION; */
/*  ROOT  -- DOUBLE PRECISION. OUTPUT SOLUTION OF THE EQUATION. */

/* INITIALIZATION. SET THE NUMBER OF ITERATION AS 0. CALL SUBROUTINE */
/* "FUNC" TO OBTAIN THE INITIAL FUNCTION VALUES F(A) AND F(B). SET */
/* DUMB VALUES FOR THE VARIABLES "E" AND "FE". */

    itnum = 0;
//EBR
//    func_(nprob, a, &fa);
//    func_(nprob, b, &fb);
    fa = f(*a);
    fb = f(*b);
//EBR
  if (isign_(&fa)*isign_(&fb) > 0)
    throw err_rroot();

    e = 1e5;
    fe = 1e5;

/* ITERATION STARTS. THE ENCLOSING INTERVAL BEFORE EXECUTING THE */
/* ITERATION IS RECORDED AS [A0, B0]. */

L10:
    a0 = *a;
    b0 = *b;

/* UPDATES THE NUMBER OF ITERATION. */

    ++itnum;

/* CALCULATES THE TERMINATION CRITERION. STOPS THE PROCEDURE IF THE */
/* CRITERION IS SATISFIED. */

    if (abs(fb) <= abs(fa)) {
  tole_(b, &tol, neps, eps);
    } else {
  tole_(a, &tol, neps, eps);
    }
    if (*b - *a <= tol) {
  goto L400;
    }

/* FOR THE FIRST ITERATION, SECANT STEP IS TAKEN. */

    if (itnum == 1) {
  c__ = *a - fa / (fb - fa) * (*b - *a);

/* CALL SUBROUTINE "BRACKT" TO GET A SHRINKED ENCLOSING INTERVAL AS */
/* WELL AS TO UPDATE THE TERMINATION CRITERION. STOP THE PROCEDURE */
/* IF THE CRITERION IS SATISFIED OR THE EXACT SOLUTION IS OBTAINED. */

        brackt_(f, a, b, &c__, &fa, &fb, &tol, neps, eps, &d__, &fd);
  if (fa == 0. || *b - *a <= tol) {
      goto L400;
  }
  goto L10;
    }

/* STARTING WITH THE SECOND ITERATION, IN THE FIRST TWO STEPS, EITHER */
/* QUADRATIC INTERPOLATION IS USED BY CALLING THE SUBROUTINE "NEWQUA" */
/* OR THE CUBIC INVERSE INTERPOLATION IS USED BY CALLING THE SUBROUTINE */
/* "PZERO". IN THE FOLLOWING, IF "PROF" IS NOT EQUAL TO 0, THEN THE */
/* FOUR FUNCTION VALUES "FA", "FB", "FD", AND "FE" ARE DISTINCT, AND */
/* HENCE "PZERO" WILL BE CALLED. */

//EBR
//    prof = (fa - fb) * (fa - fd)
//         * (fa - fe) * (fb - fd)
//         * (fb - fe) * (fd - fe);
    prof = min((fa - fb), 1.) * min((fa - fd), 1.)
         * min((fa - fe), 1.) * min((fb - fd), 1.)
         * min((fb - fe), 1.) * min((fd - fe), 1.);


    if (itnum == 2 || prof == 0.) {
  newqua_(a, b, &d__, &fa, &fb, &fd, &c__, &c__2);
    } else {
  pzero_(a, b, &d__, &e, &fa, &fb, &fd, &fe, &c__);
  if ((c__ - *a) * (c__ - *b) >= 0.) {
      newqua_(a, b, &d__, &fa, &fb, &fd, &c__, &c__2);
  }
    }
    e = d__;
    fe = fd;

/* CALL SUBROUTINE "BRACKT" TO GET A SHRINKED ENCLOSING INTERVAL AS */
/* WELL AS TO UPDATE THE TERMINATION CRITERION. STOP THE PROCEDURE */
/* IF THE CRITERION IS SATISFIED OR THE EXACT SOLUTION IS OBTAINED. */

    brackt_(f, a, b, &c__, &fa, &fb, &tol, neps, eps, &d__, &fd);
    if (fa == 0. || *b - *a <= tol) {
  goto L400;
    }
//EBR
//    prof = (fa - fb) * (fa - fd)
//         * (fa - fe) * (fb - fd)
//         * (fb - fe) * (fd - fe);
    prof = min((fa - fb), 1.) * min((fa - fd), 1.)
         * min((fa - fe), 1.) * min((fb - fd), 1.)
         * min((fb - fe), 1.) * min((fd - fe), 1.);
    if (prof == 0.) {
  newqua_(a, b, &d__, &fa, &fb, &fd, &c__, &c__3);
    } else {
  pzero_(a, b, &d__, &e, &fa, &fb, &fd, &fe, &c__);
  if ((c__ - *a) * (c__ - *b) >= 0.) {
      newqua_(a, b, &d__, &fa, &fb, &fd, &c__, &c__3);
  }
    }

/* CALL SUBROUTINE "BRACKT" TO GET A SHRINKED ENCLOSING INTERVAL AS */
/* WELL AS TO UPDATE THE TERMINATION CRITERION. STOP THE PROCEDURE */
/* IF THE CRITERION IS SATISFIED OR THE EXACT SOLUTION IS OBTAINED. */

    brackt_(f, a, b, &c__, &fa, &fb, &tol, neps, eps, &d__, &fd);
    if (fa == 0. || *b - *a <= tol) {
  goto L400;
    }
    e = d__;
    fe = fd;

/* TAKES THE DOUBLE-SIZE SECANT STEP. */

    if (abs(fa) < abs(fb)) {
  u = *a;
  fu = fa;
    } else {
  u = *b;
  fu = fb;
    }
    c__ = u - fu / (fb - fa) * 2. * (*b - *a);
    if ((d__1 = c__ - u, abs(d__1)) > (*b - *a) * .5) {
  c__ = *a + (*b - *a) * .5;
    }

/* CALL SUBROUTINE "BRACKT" TO GET A SHRINKED ENCLOSING INTERVAL AS */
/* WELL AS TO UPDATE THE TERMINATION CRITERION. STOP THE PROCEDURE */
/* IF THE CRITERION IS SATISFIED OR THE EXACT SOLUTION IS OBTAINED. */

    brackt_(f, a, b, &c__, &fa, &fb, &tol, neps, eps, &d__, &fd);
    if (fa == 0. || *b - *a <= tol) {
  goto L400;
    }

/* DETERMINES WHETHER AN ADDITIONAL BISECTION STEP IS NEEDED. AND TAKES */
/* IT IF NECESSARY. */

    if (*b - *a < (b0 - a0) * .5) {
  goto L10;
    }
    e = d__;
    fe = fd;

/* CALL SUBROUTINE "BRACKT" TO GET A SHRINKED ENCLOSING INTERVAL AS */
/* WELL AS TO UPDATE THE TERMINATION CRITERION. STOP THE PROCEDURE */
/* IF THE CRITERION IS SATISFIED OR THE EXACT SOLUTION IS OBTAINED. */

    d__1 = *a + (*b - *a) * .5;
    brackt_(f, a, b, &d__1, &fa, &fb, &tol, neps, eps, &d__, &fd);
    if (fa == 0. || *b - *a <= tol) {
  goto L400;
    }
    goto L10;

/* TERMINATES THE PROCEDURE AND RETURN THE "ROOT". */

L400:
    *root = *a;
    return 0;
} /* rroot_ */

/* Subroutine */ int brackt_(const D2D &f, doublereal *a, doublereal *b,
  doublereal *c__, doublereal *fa, doublereal *fb, doublereal *tol, 
  integer *neps, doublereal *eps, doublereal *d__, doublereal *fd)
{
    doublereal fc;

/* GIVEN CURRENT ENCLOSING INTERVAL [A,B] AND A NUMBER C IN (A,B), IF */
/* F(C)=0 THEN SETS THE OUTPUT A=C. OTHERWISE DETERMINES THE NEW */
/* ENCLOSING INTERVAL: [A,B]=[A,C] OR [A,B]=[C,B]. ALSO UPDATES THE */
/* TERMINATION CRITERION CORRESPONDING TO THE NEW ENCLOSING INTERVAL. */
/*  NPROB   -- INTEGER. INDICATING THE PROBLEM TO BE SOLVED; */
/*  A,B     -- DOUBLE PRECISION. [A,B] IS INPUT AS THE CURRENT */
/*             ENCLOSING INTERVAL AND OUTPUT AS THE SHRINKED NEW */
/*             ENCLOSING INTERVAL; */
/*  C       -- DOUBLE PRECISION. USED TO DETERMINE THE NEW ENCLOSING */
/*             INTERVAL; */
/*  D       -- DOUBLE PRECISION. OUTPUT: IF THE NEW ENCLOSING INTERVAL */
/*             IS [A,C] THEN D=B, OTHERWISE D=A; */
/*  FA,FB,FD-- DOUBLE PRECISION. FA=F(A), FB=F(B), AND FD=F(D); */
/*  TOL     -- DOUBLE PRECISION. INPUT AS THE CURRENT TERMINATION */
/*             CRITERION AND OUTPUT AS THE UPDATED TERMINATION */
/*             CRITERION ACCORDING TO THE NEW ENCLOSING INTERVAL; */
/*  NEPS    -- INTEGER. USED TO DETERMINE THE TERMINATION CRITERION; */
/*  EPS     -- DOUBLE PRECISION. USED IN THE TERMINATION CRITERION. */

/* ADJUST C IF (B-A) IS VERY SMALL OR IF C IS VERY CLOSE TO A OR B. */

    *tol *= .7;
    if (*b - *a <= *tol * 2.) {
  *c__ = *a + (*b - *a) * .5;
    } else if (*c__ <= *a + *tol) {
  *c__ = *a + *tol;
    } else {
  if (*c__ >= *b - *tol) {
      *c__ = *b - *tol;
  }
    }

/* CALL SUBROUTINE "FUNC" TO OBTAIN F(C) */

//EBR
//    func_(nprob, c__, &fc);
    fc = f(*c__);

/* IF F(C)=0, THEN SET A=C AND RETURN. THIS WILL TERMINATE THE */
/* PROCEDURE IN SUBROUTINE "RROOT" AND GIVE THE EXACT SOLUTION OF */
/* THE EQUATION F(X)=0. */

    if (fc == 0.) {
  *a = *c__;
  *fa = 0.;
  *d__ = 0.;
  *fd = 0.;
  return 0;
    }

/* IF F(C) IS NOT ZERO, THEN DETERMINE THE NEW ENCLOSING INTERVAL. */

    if (isign_(fa) * isign_(&fc) < 0) {
  *d__ = *b;
  *fd = *fb;
  *b = *c__;
  *fb = fc;
    } else {
  *d__ = *a;
  *fd = *fa;
  *a = *c__;
  *fa = fc;
    }

/* UPDATE THE TERMINATION CRITERION ACCORDING TO THE NEW ENCLOSING */
/* INTERVAL. */

    if (abs(*fb) <= abs(*fa)) {
  tole_(b, tol, neps, eps);
    } else {
  tole_(a, tol, neps, eps);
    }

/* END OF THE SUBROUTINE. */

    return 0;
} /* brackt_ */

int isign_(doublereal *x)
{
    /* System generated locals */
    integer ret_val;

/* INDICATES THE SIGN OF THE VARIABLE "X". */
/*  X     -- DOUBLE PRECISION. */
/*  ISIGN -- INTEGER. */
    if (*x > 0.) {
  ret_val = 1;
    } else if (*x == 0.) {
  ret_val = 0;
    } else {
  ret_val = -1;
    }
    return ret_val;
} /* isign_ */

/* Subroutine */ int tole_(doublereal *b, doublereal *tol, integer *neps, 
  doublereal *eps)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;

/* DETERMINES THE TERMINATION CRITERION. */
/*  B    -- DOUBLE PRECISION. */
/*  NEPS -- INTEGER. */
/*  EPS  -- DOUBLE PRECISION. */
/*  TOL  -- DOUBLE PRECISION. OUTPUT AS THE TERMINATION CRITERION. */
/*           TOL =2*(2*EPS*|B| + 10D-{NEPS}),  IF NEPS IS NOT 1000; */
/*    AND    TOL =2*(2*EPS*|B|),               IF NEPS = 1000. */
    if (*neps == 1000) {
  *tol = 0.;
    } else {
  *tol = 1.;
  i__1 = *neps;
  for (i__ = 1; i__ <= i__1; ++i__) {
      *tol /= 10.;
/* L10: */
  }
    }
    *tol += abs(*b) * 2. * *eps;
    *tol *= 2.;
    return 0;
} /* tole_ */

/* Subroutine */ int newqua_(doublereal *a, doublereal *b, doublereal *d__, 
  doublereal *fa, doublereal *fb, doublereal *fd, doublereal *c__, 
  integer *k)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    doublereal a0, a1, a2, pc;
    integer ierror;
    doublereal pdc;

/* USES K NEWTON STEPS TO APPROXIMATE THE ZERO IN (A,B) OF THE */
/* QUADRATIC POLYNOMIAL INTERPOLATING F(X) AT A, B, AND D. SAFEGUARD */
/* IS USED TO AVOID OVERFLOW. */
/*  A,B,D,FA,FB,FD -- DOUBLE PRECISION. D LIES OUTSIDE THE INTERVAL */
/*                    [A,B]. FA=F(A), FB=F(B), AND FD=F(D). F(A)F(B)<0. */
/*  C              -- DOUBLE PRECISION. OUTPUT AS THE APPROXIMATE ZERO */
/*                    IN (A,B) OF THE QUADRATIC POLYNOMIAL. */
/*  K              -- INTEGER. INPUT INDICATING THE NUMBER OF NEWTON */
/*                    STEPS TO TAKE. */

/* INITIALIZATION. FIND THE COEFFICIENTS OF THE QUADRATIC POLYNOMIAL. */

    ierror = 0;
    a0 = *fa;
    a1 = (*fb - *fa) / (*b - *a);
    a2 = ((*fd - *fb) / (*d__ - *b) - a1) / (*d__ - *a);

/* SAFEGUARD TO AVOID OVERFLOW. */

L10:
    if (a2 == 0. || ierror == 1) {
  *c__ = *a - a0 / a1;
  return 0;
    }

/* DETERMINE THE STARTING POINT OF NEWTON STEPS. */

    if (isign_(&a2) * isign_(fa) > 0) {
  *c__ = *a;
    } else {
  *c__ = *b;
    }

/* START THE SAFEGUARDED NEWTON STEPS. */

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
  if (ierror == 0) {
      pc = a0 + (a1 + a2 * (*c__ - *b)) * (*c__ - *a);
      pdc = a1 + a2 * (*c__ * 2. - (*a + *b));
      if (pdc == 0.) {
    ierror = 1;
      } else {
    *c__ -= pc / pdc;
      }
  }
/* L20: */
    }
    if (ierror == 1) {
  goto L10;
    }
    return 0;
} /* newqua_ */

/* Subroutine */ int pzero_(doublereal *a, doublereal *b, doublereal *d__, 
  doublereal *e, doublereal *fa, doublereal *fb, doublereal *fd, 
  doublereal *fe, doublereal *c__)
{
    long double d21, d31, d32, q11, q21, q31, q22, q32, q33;

/* USES CUBIC INVERSE INTERPOLATION OF F(X) AT A, B, D, AND E TO */
/* GET AN APPROXIMATE ROOT OF F(X). THIS PROCEDURE IS A SLIGHT */
/* MODIFICATION OF AITKEN-NEVILLE ALGORITHM FOR INTERPOLATION */
/* DESCRIBED BY STOER AND BULIRSCH IN "INTRO. TO NUMERICAL ANALYSIS" */
/* SPRINGER-VERLAG. NEW YORK (1980). */
/*  A,B,D,E,FA,FB,FD,FE -- DOUBLE PRECISION. D AND E LIE OUTSIDE */
/*                         THE INTERVAL [A,B]. FA=F(A), FB=F(B), */
/*                         FD=F(D), AND FE=F(E). */
/*  C                   -- DOUBLE PRECISION. OUTPUT OF THE SUBROUTINE. */

    q11 = (*d__ - *e) * *fd / (*fe - *fd);
    q21 = (*b - *d__) * *fb / (*fd - *fb);
    q31 = (*a - *b) * *fa / (*fb - *fa);
    d21 = (*b - *d__) * *fd / (*fd - *fb);
    d31 = (*a - *b) * *fb / (*fb - *fa);
    q22 = (d21 - q11) * *fb / (*fe - *fb);
    q32 = (d31 - q21) * *fa / (*fd - *fa);
    d32 = (d31 - q21) * *fd / (*fd - *fa);
    q33 = (d32 - q22) * *fa / (*fe - *fa);

/* CALCULATE THE OUTPUT C. */

    *c__ = q31 + q32 + q33;
    *c__ = *a + *c__;
    return 0;
} /* pzero_ */

/* Subroutine */ int rmp_(doublereal *rel)
{
    doublereal beta, a, b;

/* CALCULATES THE RELATIVE MACHINE PRECISION (RMP). */
/*  REL -- DOUBLE PRECISION. OUTPUT OF RMP. */

    beta = 2.;
    a = 1.;
L10:
    b = a + 1.;
    if (b > 1.) {
  a /= beta;
  goto L10;
    }
    *rel = a * beta;
    return 0;
} /* rmp_ */

