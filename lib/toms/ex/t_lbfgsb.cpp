/* driver1.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include <toms.h>
#include "g2c.h"
#include "f2clib.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;

/*                             DRIVER 1 */
/*     -------------------------------------------------------------- */
/*                SIMPLE DRIVER FOR L-BFGS-B (version 2.3) */
/*     -------------------------------------------------------------- */

/*        L-BFGS-B is a code for solving large nonlinear optimization */
/*        problems with simple bounds on the variables. */

/*        The code can also be used for unconstrained problems and is */
/*        as efficient for these problems as the earlier limited memory */
/*        code L-BFGS. */

/*        This is the simplest driver in the package. It uses all the */
/*        default settings of the code. */


/*     References: */

/*        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited */
/*        memory algorithm for bound constrained optimization'', */
/*        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. */

/*        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN */
/*        Subroutines for Large Scale Bound Constrained Optimization'' */
/*        Tech. Report, NAM-11, EECS Department, Northwestern University, */
/*        1994. */

/*        (Postscript files of these papers are available via anonymous */
/*        ftp to ece.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.) */

/*                              *  *  * */

/*        NEOS, November 1994. (Latest revision April 1997.) */
/*        Optimization Technology Center. */
/*        Argonne National Laboratory and Northwestern University. */
/*        Written by */
/*                           Ciyou Zhu */
/*        in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */

/*     NOTE: The user should adapt the subroutine 'timer' if 'etime' is */
/*           not available on the system.  An example for system */
/*           AIX Version 3.2 is available at the end of this driver. */

/*     ************** */

extern "C" /* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_16[] = "(/5x,\002Solving sample problem.\002/5x,\002 (f "
      "= 0.0 at the optimal solution.)\002/)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    char task[60];
    doublereal f, g[1024];
    integer i__;
    doublereal l[1024];
    integer m, n;
    doublereal u[1024], x[1024], factr;
    char csave[60];
    doublereal dsave[29];
    integer isave[44];
    logical lsave[4];
    doublereal pgtol, t1, t2, wa[42227];
    integer iprint, nbd[1024], iwa[3072];

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 6, 0, fmt_16, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };


/*     This simple driver demonstrates how to call the L-BFGS-B code to */
/*     solve a sample problem (the extended Rosenbrock function */
/*     subject to bounds on the variables). The dimension n of this */
/*     problem is variable. */
/*     nmax  is the dimension of the largest problem to be solved. */
/*     mmax  is the maximum number of limited memory corrections. */
/*     lenwa is the corresponding real workspace required. */
/*     Declare the variables needed by the code. */
/*     A description of all these variables is given at the end of */
/*     the driver. */
/*     Declare a few additional variables for this sample problem. */
/*     We wish to have output at every iteration. */
    iprint = 1;
/*     We specify the tolerances in the stopping criteria. */
    factr = 1e7;
    pgtol = 1e-5;
/*     We specify the dimension n of the sample problem and the number */
/*     m of limited memory corrections stored.  (n and m should not */
/*     exceed the limits nmax and mmax respectively.) */
    n = 25;
    m = 5;
/*     We now provide nbd which defines the bounds on the variables: */
/*     l   specifies the lower bounds, */
/*     u   specifies the upper bounds. */
/*     First set bounds on the odd-numbered variables. */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; i__ += 2) {
  nbd[i__ - 1] = 2;
  l[i__ - 1] = 1.;
  u[i__ - 1] = 100.;
/* L10: */
    }
/*     Next set bounds on the even-numbered variables. */
    i__1 = n;
    for (i__ = 2; i__ <= i__1; i__ += 2) {
  nbd[i__ - 1] = 2;
  l[i__ - 1] = -100.;
  u[i__ - 1] = 100.;
/* L12: */
    }
/*     We now define the starting point. */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
  x[i__ - 1] = 3.;
/* L14: */
    }
/*     We now write the heading of the output. */
    s_wsfe(&io___11);
    e_wsfe();
/*     We start the iteration by initializing task. */
    s_copy(task, "START", 60L, 5L);
/*     ------- The beginning of the loop ---------- */
L111:
/*     This is the call to the L-BFGS-B code. */
//    setulb_(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &
//      iprint, csave, lsave, isave, dsave, 60L, 60L);
    setulb_(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &
      iprint, csave, lsave, isave, dsave);
    if (s_cmp(task, "FG", 2L, 2L) == 0) {
/*        The minimization routine has returned to request the */
/*        function f and gradient g values at the current x. */
/*        Compute function value f for the sample problem. */
/* Computing 2nd power */
        d__1 = x[0] - 1.;
  f = d__1 * d__1 * .25;
  i__1 = n;
  for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing 2nd power */
      d__2 = x[i__ - 2];
/* Computing 2nd power */
      d__1 = x[i__ - 1] - d__2 * d__2;
      f += d__1 * d__1;
/* L20: */
  }
  f *= 4.;
/*        Compute gradient g for the sample problem. */
/* Computing 2nd power */
  d__1 = x[0];
  t1 = x[1] - d__1 * d__1;
  g[0] = (x[0] - 1.) * 2. - x[0] * 16. * t1;
  i__1 = n - 1;
  for (i__ = 2; i__ <= i__1; ++i__) {
      t2 = t1;
/* Computing 2nd power */
      d__1 = x[i__ - 1];
      t1 = x[i__] - d__1 * d__1;
      g[i__ - 1] = t2 * 8. - x[i__ - 1] * 16. * t1;
/* L22: */
  }
  g[n - 1] = t1 * 8.;
/*        Go back to the minimization routine. */
  goto L111;
    } else if (s_cmp(task, "NEW_X", 5L, 5L) == 0) {
/*        The minimization routine has returned with a new iterate, */
/*        and we have opted to continue the iteration. */
  goto L111;
    } else {
/*        We terminate execution when task is neither FG nor NEW_X. */
/*        We print the information contained in the string task */
/*        if the default output is not used and the execution is */
/*        not stopped intentionally by the user. */
  if (iprint <= -1 && s_cmp(task, "STOP", 4L, 4L) != 0) {
      s_wsle(&io___23);
      do_lio(&c__9, &c__1, task, 60L);
      e_wsle();
  }
    }
/*     ---------- The end of the loop ------------- */
    s_stop("", 0L);
    return 0;
} /* MAIN__ */

/* Main program alias */ int main () { MAIN__ (); return 0; }
