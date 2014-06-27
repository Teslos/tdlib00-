//algorithm 666

//this function is from netlib/toms, see readme.txt for license
//I convertied it to C++ by means of f2c and made some minor chages
//Evgenii Rudnyi rudnyi@comp.chem.msu.su

/* intsub.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "toms.h"
#include "f2c.h"
#include "f2clib.h"

int intsub_(const DL2D &fnc, long *n, double *xo, double *h__, double
            *delta, double *epsilo, long *icon, long *inf1, double *as,
            double *vas, long *inf2, double *wa, long *lwa);

int chapol_(const DL2D &fnc, long *n, long *nv, long *ndg, double *xo,
            double *h__, double *epsmch, double *delta, double *epsilo,
            long *icon, double *b, double *c__, double *po, double *spo,
            double *cp, double *a, long *inf1, double *x, double *fx,
            double *wm, double *wa1, double *wa2);

int genbis_(const DL2D &fnc, long *n, long *nv, long *ns, long *ndg,
            double *b, double *c__, double *cp, double *epsmch, double
            *epsilo, double *cpr, double *dm, double *dg, double *bp,
            double *fvb, long *inf2, double *a, double *arp, double *wm,
            double *wa1, double *wa2, double *wa3);

long lsinnc_(long *n, double *sx, double *epsilo);

long iscomc_(long *n, double *sx, double *sy, long *incy);

double sgn_(double *sv);

int saverc_(long *n, double *sx, long *incx, double *sy);

int sscalc_(long *n, double *v1, double *sx, long *lx, long *incx,
            double *v2, double *sy, long *ly, long *incy, double *sz,
            long *lz, long *incz);

int smnlgc_(long *n, double *sx, long *incx, double *smn, double *slg);


int chabis(const DL2D &fnc, long n, double *x, double *h,
           double *vas, double delta)
{
  if (!delta)
//    delta = 0.0625;
    delta = 1e-6;
  double epsilo = 1e-6;
  long icon = 1;
  long inf1, inf2;
  long lwa = 2*n + (6*n+1) * long(pow(2., int(n)));
  double *wa = new double[lwa + 1]; //plus 1 just in case
  double *as = new double[n];
  intsub_(fnc, &n, x, h, &delta, &epsilo, &icon, &inf1, as, vas, &inf2,
          wa, &lwa);
  delete [] wa;
  for (int i = 0; i < n; ++i)
    x[i] = as[i];
  delete [] as;
  if (inf2 == 0 && inf1 != 4)
    throw err_chabis("chabis: no convergence", inf2 + inf1*10 + 100);
  return inf2 + inf1*10 + 100;
}


/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b41 = 1.;
static doublereal c_b96 = 2.;

/* Subroutine */ int intsub_(const DL2D &fnc, integer *n, doublereal *xo, doublereal
  *h__, doublereal *delta, doublereal *epsilo, integer *icon, integer *
  inf1, doublereal *as, doublereal *vas, integer *inf2, doublereal *wa, 
  integer *lwa)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer j, nn, ns, nv;
    doublereal epsmch;
    integer ld1, ld2, ld3, ld4, ld5, ld6, ld7, ld8, ld9, ld10, ld11, ndg;
//    doublereal tol;

/* ********************************************************************** 
*/
/*                                                                     * 
*/
/*     SUBROUTINE INTSUB                                               * 
*/
/*                                                                     * 
*/
/*       This is an interface subroutine between the  main program and * 
*/
/*       the subroutines  CHAPOL  and  GENBIS,  through which a single * 
*/
/*       workspace array WA is passed. INTSUB evokes CHAPOL and GENBIS * 
*/
/*       in such a way that the addresses of  the corresponding arrays * 
*/
/*       within the subroutines are computed in WA. In this way, their * 
*/
/*       arrays can appear with variable dimensions.                   * 
*/
/*                                                                     * 
*/
/*                                                                     * 
*/
/*     The subroutine statement is :                                   * 
*/
/*                                                                     * 
*/
/*     SUBROUTINE INTSUB( FNC, N, XO, H, DELTA, EPSILO, ICON, INF1,    * 
*/
/*    +                   AS, VAS, INF2, WA, LWA )                     * 
*/
/*                                                                     * 
*/
/*     where                                                           * 
*/
/*                                                                     * 
*/
/*       FNC is the name of the user-supplied function which evaluates * 
*/
/*         components of the given function. FNC should be declared in * 
*/
/*         an external  statement in  the user - calling  program, and * 
*/
/*         should be written as follows :                              * 
*/
/*                                                                     * 
*/
/*         REAL FUNCTION FNC( X, IFLAG )                               * 
*/
/*         INTEGER IFLAG                                               * 
*/
/*         REAL X(N)                                                   * 
*/
/*         ------------------------------------------------------      * 
*/
/*         Calculate the IFLAG-th component of the function at X.      * 
*/
/*         ------------------------------------------------------      * 
*/
/*         RETURN                                                      * 
*/
/*         END                                                         * 
*/
/*                                                                     * 
*/
/*       N is a positive integer input variable which defines the num- * 
*/
/*         ber of equations and variables.                             * 
*/
/*                                                                     * 
*/
/*       XO is an input array of length  N  which defines  the initial * 
*/
/*         guess of the solution.                                      * 
*/
/*                                                                     * 
*/
/*       H is an input array of length N, with positive entries, which * 
*/
/*         determines the stepsizes in each coordinate direction.      * 
*/
/*                                                                     * 
*/
/*       DELTA is a positive input variable which determines the accu- * 
*/
/*         racy of  the computation of the roots, of the components of * 
*/
/*         the  given function, which are located on  the edges of the * 
*/
/*         Initial N-Polyhedron and are used for the construction of a * 
*/
/*         Characteristic  N-Polyhedron.  If  DELTA  is less  than the * 
*/
/*         machine precision  EPSMCH,  DELTA takes the value  0.0625 . * 
*/
/*         The value of EPSMCH is computed within INTSUB.              * 
*/
/*                                                                     * 
*/
/*       EPSILO  is a nonnegative  input  variable. Termination occurs * 
*/
/*         when the algorithm estimates that the Infinity Norm, of the * 
*/
/*         function  values  at an  approximate  solution, is at  most * 
*/
/*         EPSILO. If EPSILO is less than the machine precision EPSMCH,* 
*/
/*         EPSILO becomes equal to  EPSMCH.                            * 
*/
/*                                                                     * 
*/
/*       ICON  is an integer input variable. If  ICON  is equal to  1, * 
*/
/*         then GENBIS is evoked even if a Characteristic N-Polyhedron * 
*/
/*         has not been constructed.                                   * 
*/
/*                                                                     * 
*/
/*       INF1 is an integer output variable set as follows :           * 
*/
/*                                                                     * 
*/
/*         INF1 = 0   Improper input parameters.                       * 
*/
/*                    CHAPOL and GENBIS have not been evoked.          * 
*/
/*                                                                     * 
*/
/*         INF1 = 1   A  Characteristic  N-Polyhedron  has  been  con- * 
*/
/*                    structed.                                        * 
*/
/*                                                                     * 
*/
/*         INF1 = 2   A Characteristic N-Polyhedron has not  been con- * 
*/
/*                    structed.                                        * 
*/
/*                                                                     * 
*/
/*         INF1 = 3   More than two vertices of the constructed  Char- * 
*/
/*                    acteristic  N-Polyhedron are located on the same * 
*/
/*                    edge of the Initial N-Polyhedron.                * 
*/
/*                                                                     * 
*/
/*         INF1 = 4   An approximate solution according to the  preci- * 
*/
/*                    sion  EPSILO has been found during the construc- * 
*/
/*                    tion of the Characteristic N-Polyhedron.         * 
*/
/*                    GENBIS has not been evoked.                      * 
*/
/*                                                                     * 
*/
/*       AS is an output array of length N  which determines the final * 
*/
/*         approximate solution.                                       * 
*/
/*                                                                     * 
*/
/*       VAS is an output array of length N  which specifies the func- * 
*/
/*          tion values at  AS.                                        * 
*/
/*                                                                     * 
*/
/*       INF2 is an integer output variable set as follows :           * 
*/
/*                                                                     * 
*/
/*         INF2 = 0   GENBIS has not been evoked.                      * 
*/
/*                                                                     * 
*/
/*         INF2 = 1   The solution has been found  within the required * 
*/
/*                    accuracy of  EPSILO.                             * 
*/
/*                                                                     * 
*/
/*         INF2 = 2   Iterations have reached their upper limit.       * 
*/
/*                    The answer may not be accurate.                  * 
*/
/*                                                                     * 
*/
/*         INF2 = 3   The  length  of  the longest  diagonal  has been * 
*/
/*                    found to be less than  2.0 * N * EPSILO.         * 
*/
/*                                                                     * 
*/
/*       WA is a workspace array of length LWA.                        * 
*/
/*                                                                     * 
*/
/*       LWA is a positive integer input variable  not less  than  the * 
*/
/*         value of ( 2*N + (6*N+1) * 2**N ).                          * 
*/
/*                                                                     * 
*/
/*                                                                     * 
*/
/*     Subprograms required :                                          * 
*/
/*                                                                     * 
*/
/*       USER-Supplied ..... FNC                                       * 
*/
/*                                                                     * 
*/
/*       CHABIS-Supplied ... CHAPOL, GENBIS                            * 
*/
/*                                                                     * 
*/
/*                                                                     * 
*/
/*                                        CHABIS. Version of June 1988 * 
*/
/*                                        Michael N. Vrahatis          * 
*/
/*                                                                     * 
*/
/* ********************************************************************** 
*/

/*     Check the input parameters for errors.                          * 
*/

    /* Parameter adjustments */
    --vas;
    --as;
    --h__;
    --xo;
    --wa;

    /* Function Body */
    *inf1 = 0;
    *inf2 = 0;
    nv = pow_ii(&c__2, n);
    if (*n <= 1 || *lwa < (*n << 1) + (*n * 6 + 1) * nv) {
  goto L30;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
  if (h__[j] <= 0.) {
      goto L30;
  }
/* L10: */
    }

/*     Compute the machine precision.                                  * 
*/

    epsmch = DBL_EPSILON;

//L20:
//    epsmch /= 2.;
//    tol = epsmch + 1.;
//    if (tol > 1.) {
//        goto L20;
//    }
//    epsmch *= 2.;

/*     Determine if the input values of EPSILO and DELTA are less than * 
*/
/*     the machine precision.                                          * 
*/

    if (*epsilo < epsmch) {
  *epsilo = epsmch;
    }
    if (*delta < epsmch) {
  *delta = .0625;
    }

/*     Call CHAPOL.                                                    * 
*/

    ndg = nv / 2;
    nn = *n * nv;
    ld1 = nn + 1;
    ld2 = ld1 + nn;
    ld3 = ld2 + nn;
    ld4 = ld3 + nn;
    ld5 = ld4 + nn;
    ld6 = ld5 + nv;
    ld7 = ld6 + nn;
    ld8 = ld7 + *n;
    chapol_(fnc, n, &nv, &ndg, &xo[1], &h__[1], &epsmch, delta, epsilo,
      icon, &wa[1], &wa[ld1], &wa[ld3], &wa[ld4], &wa[ld2], &wa[ld5], 
      inf1, &as[1], &vas[1], &wa[ld6], &wa[ld7], &wa[ld8]);
    if (*inf1 == 4) {
  goto L30;
    }
    if (*inf1 == 2 && *icon != 1) {
  goto L30;
    }

/*     Call GENBIS.                                                    * 
*/

    ns = *n * ndg;
    ld9 = ld6 + ns;
    ld10 = ld9 + ndg;
    ld11 = ld10 + *n;
    genbis_(fnc, n, &nv, &ns, &ndg, &wa[1], &wa[ld1], &wa[ld2], &epsmch,
             epsilo, &wa[ld3], &wa[ld6], &wa[ld9], &as[1], &vas[1], inf2, &wa[
      ld5], &wa[ld10], &wa[ld4], &wa[ld11], &xo[1], &h__[1]);
L30:
    return 0;

/*     Last statement of the interface subroutine INTSUB.              * 
*/

} /* intsub_ */

/* ---------------------------------------------------------------------* */

/* Subroutine */ int chapol_(const DL2D &fnc, integer *n, integer *nv, integer *ndg,
  doublereal *xo, doublereal *h__, doublereal *epsmch, doublereal *
  delta, doublereal *epsilo, integer *icon, doublereal *b, doublereal *
  c__, doublereal *po, doublereal *spo, doublereal *cp, doublereal *a, 
  integer *inf1, doublereal *x, doublereal *fx, doublereal *wm, 
  doublereal *wa1, doublereal *wa2)
{
    /* System generated locals */
    integer b_dim1, b_offset, c_dim1, c_offset, po_dim1, po_offset, spo_dim1, 
      spo_offset, cp_dim1, cp_offset, wm_dim1, wm_offset, i__1, i__2, 
      i__3, i__4, i__5;
    doublereal d__1;

    /* Local variables */
    doublereal rlen, ritr;
    logical test;
    integer i__, j, k, l;
    doublereal s, t, dstar;
    integer i1, j1, j2, j3, j4, j5, i2, i3, i4;
    doublereal fe, el, er, rn, ro, rt;
    doublereal sc1, sc2;
    integer itr, nnv;
    doublereal sro, srn;

/* ********************************************************************** 
*/
/*                                                                     * 
*/
/*                                                                     * 
*/
/*       EPSILO  is a nonnegative  input  variable. Termination occurs * 
*/
/*         when the algorithm estimates that the Infinity Norm, of the * 
*/
/*         function  values  at an  approximate  solution, is at  most * 
*/
/*         EPSILO. If EPSILO is less than the machine precision EPSMCH,* 
*/
/*         EPSILO  becomes equal to  EPSMCH.                           * 
*/
/*                                                                     * 
*/
/*       ICON is an integer input variable. If  ICON is equal to 1 and * 
*/
/*         the  Characteristic  N-Polyhedron construction fails,  then * 
*/
/*         CHAPOL  will complete  the unconstructed  Characteristic N- * 
*/
/*         Polyhedron vertices, using the  Initial N-Polyhedron verti- * 
*/
/*         ces. Note that this action is taken when the Characteristic * 
*/
/*         N-Polyhedron  is incomplete and  the subroutine  GENBIS  is * 
*/
/*         utilized.                                                   * 
*/
/*                                                                     * 
*/
/*       B is an output  NV by N  matrix  which defines  the  N-binary * 
*/
/*         matrix.                                                     * 
*/
/*                                                                     * 
*/
/*       C is an output  NV by N  matrix which defines the  N-complete * 
*/
/*         matrix.                                                     * 
*/
/*                                                                     * 
*/
/*       PO is an output NV by N matrix,with entries in each row which * 
*/
/*         are  the corresponding  components  of the vertices  of the * 
*/
/*         Initial N-Polyhedron.                                       * 
*/
/*                                                                     * 
*/
/*       SPO is an output  NV  by  N  matrix, with entries in each row * 
*/
/*         which are the  corresponding  components  of the vectors of * 
*/
/*         signs of the function relative to the vertices of  PO.      * 
*/
/*                                                                     * 
*/
/*       CP is an output NV by N matrix,with entries in each row which * 
*/
/*         are  the corresponding  components  of the vertices of  the * 
*/
/*         Characteristic N-Polyhedron.                                * 
*/
/*                                                                     * 
*/
/*       A is an output array of length NV and determines which verti- * 
*/
/*         ces of the  Characteristic N-Polyhedron have not been  con- * 
*/
/*         structed. When all its entries are equal to zero, the Char- * 
*/
/*         acteristic N-Polyhedron has been constructed.               * 
*/
/*                                                                     * 
*/
/*       INF1 is an integer output variable set as follows :           * 
*/
/*                                                                     * 
*/
/*         INF1 = 0   Improper input parameters.                       * 
*/
/*                    CHAPOL and GENBIS have not been evoked.          * 
*/
/*                                                                     * 
*/
/*         INF1 = 1   A  Characteristic  N-Polyhedron  has  been  con- * 
*/
/*                    structed.                                        * 
*/
/*                                                                     * 
*/
/*         INF1 = 2   A Characteristic N-Polyhedron has not  been con- * 
*/
/*                    structed.                                        * 
*/
/*                                                                     * 
*/
/*         INF1 = 3   More than two vertices of the constructed  Char- * 
*/
/*                    acteristic  N-Polyhedron are located on the same * 
*/
/*                    edge of the Initial N-Polyhedron.                * 
*/
/*                                                                     * 
*/
/*         INF1 = 4   An approximate solution, according to the preci- * 
*/
/*                    sion EPSILO, has been found during the construc- * 
*/
/*                    tion of the Characteristic N-Polyhedron.         * 
*/
/*                    GENBIS has not been evoked.                      * 
*/
/*                                                                     * 
*/
/*       X is an output array of length  N  which  determines  a point * 
*/
/*         which will be a vertex of the  Characteristic N-Polyhedron. * 
*/
/*         In the case of INF1 = 4,  X contains the coordinates of the * 
*/
/*         approximate solution.                                       * 
*/
/*                                                                     * 
*/
/*       FX is an output array of length  N which determines the func- * 
*/
/*         tion values at  X .                                         * 
*/
/*                                                                     * 
*/
/*       WM  is an  NV by N  work matrix.                              * 
*/
/*                                                                     * 
*/
/*       WA1, WA2 are work arrays of length N.                         * 
*/
/*                                                                     * 
*/
/*                                                                     * 
*/
/*     Subprograms required :                                          * 
*/
/*                                                                     * 
*/
/*       USER-Supplied ...... FNC                                      * 
*/
/*                                                                     * 
*/
/*       BLAS-Supplied ...... SCOPY                                    * 
*/
/*                                                                     * 
*/
/*       CHABIS-Supplied .... SGN, LSINNC, ISCOMC, SAVERC, SSCALC,     * 
*/
/*                            SMNLGC                                   * 
*/
/*                                                                     * 
*/
/*       FORTRAN-Supplied ... FLOAT, LOG, INT, ABS, SIGN               * 
*/
/*                                                                     * 
*/
/*                                                                     * 
*/
/*                                        CHABIS. Version of June 1988 * 
*/
/*                                        Michael N. Vrahatis          * 
*/
/*                                                                     * 
*/
/* ********************************************************************** 
*/

/*     Set  DSTAR, a small positive number, which determines a  toler- * 
*/
/*     ance  of the roots  of  the components  of the given  function, * 
*/
/*     which are  located on  the edges of the  Initial  N-Polyhedron. * 
*/
/*     DSTAR must be greater than, or equal to, DELTA.                 * 
*/

    /* Parameter adjustments */
    --wa2;
    --wa1;
    --fx;
    --x;
    --h__;
    --xo;
    wm_dim1 = *nv;
    wm_offset = wm_dim1 + 1;
    wm -= wm_offset;
    --a;
    cp_dim1 = *nv;
    cp_offset = cp_dim1 + 1;
    cp -= cp_offset;
    spo_dim1 = *nv;
    spo_offset = spo_dim1 + 1;
    spo -= spo_offset;
    po_dim1 = *nv;
    po_offset = po_dim1 + 1;
    po -= po_offset;
    c_dim1 = *nv;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;
    b_dim1 = *nv;
    b_offset = b_dim1 + 1;
    b -= b_offset;

    /* Function Body */
    dstar = *delta + *epsmch * 2.;

/*     Construct the  N-Complete  Matrix and the matrix  PO. Also, set * 
*/
/*     all the entries of CP equal to the corresponding entries of PO. * 
*/

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
  i__2 = *n - j;
  k = pow_ii(&c__2, &i__2);
  l = k << 1;
  i__2 = *nv;
  for (i__ = 1; i__ <= i__2; ++i__) {
      i1 = i__ - 1;
      t = (doublereal) (i1 / k - (i1 / l << 1));
      b[i__ + j * b_dim1] = t;
      c__[i__ + j * c_dim1] = t * 2. - 1.;
      po[i__ + j * po_dim1] = xo[j] + t * h__[j];
      cp[i__ + j * cp_dim1] = po[i__ + j * po_dim1];
/* L10: */
  }
/* L20: */
    }

/*     Construct the indexing array, A.                                * 
*/

    i__1 = *nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
  a[i__] = (doublereal) i__;
/* L30: */
    }

/*     Construct  the matrix of signs  of the given function, relative * 
*/
/*     to vertices of the Initial N-Polyhedron.                        * 
*/

    i__1 = *nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
  dcopy_(n, &po[i__ + po_dim1], nv, &x[1], &c__1);
  i__2 = *n;
  for (j = 1; j <= i__2; ++j) {
            fx[j] = fnc(&x[1], j);
      spo[i__ + j * spo_dim1] = sgn_(&fx[j]);
/* L40: */
  }
  test = lsinnc_(n, &fx[1], epsilo);
  if (test) {
      goto L270;
  }

/*        Start the construction of a Characteristic  N-Polyhedron.   
 * */

  dcopy_(n, &spo[i__ + spo_dim1], nv, &wa1[1], &c__1);
  k = iscomc_(n, &wa1[1], &c__[c_offset], nv);
  if (k > 0) {
      if (a[k] > 0.) {
    dcopy_(n, &x[1], &c__1, &cp[k + cp_dim1], nv);
    a[k] = 0.;
      }
  }
/* L50: */
    }

/*     Determine if a  Characteristic N-Polyhedron has been completed. * 
*/

    test = lsinnc_(nv, &a[1], epsmch);
    if (test) {
  goto L250;
    }

/*     Compute  the roots of the components of the function, which are * 
*/
/*     located on the edges of the  Initial N-Polyhedron, and use them * 
*/
/*     to construct a  Characteristic N-Polyhedron.                    * 
*/

//cout << "roots" << endl;
    i__1 = *n;
    for (j1 = 1; j1 <= i__1; ++j1) {
  j2 = pow_ii(&c__2, &j1) - 2;
  i__2 = *n - j1;
  j3 = pow_ii(&c__2, &i__2);
  i__2 = j2;
  for (j4 = 0; j4 <= i__2; j4 += 2) {
      i__3 = j3;
      for (j5 = 1; j5 <= i__3; ++j5) {
    i1 = j4 * j3 + j5;
    i2 = i1 + j3;
    dcopy_(n, &po[i1 + po_dim1], nv, &x[1], &c__1);
    el = po[i1 + j1 * po_dim1];
    er = po[i2 + j1 * po_dim1];
    rlen = er - el;

/*              Compute  the number  of iterations which are  
required * */
/*              to obtain a root to the predetermined accuracy
, DELTA. * */

    ritr = log(rlen / *delta) / log(2.);
    k = (integer) ritr;
    itr = k + 1;
    if (ritr == (doublereal) k) {
        --itr;
    }

    i__4 = *n;
    for (j = 1; j <= i__4; ++j) {
        if (spo[i1 + j * spo_dim1] * spo[i2 + j * spo_dim1] > 0.) 
          {
      goto L100;
        }
        ro = el;
        sro = spo[i1 + j * spo_dim1];
        s = rlen / 2.;
        rn = ro + s;
        rt = er;

/*                 Begin the iterations.                  
             * */

        i__5 = itr;
        for (i3 = 1; i3 <= i__5; ++i3) {
      x[j1] = rn;
//cout << "2" << endl;
                        fe = fnc(&x[1], j);
      srn = sgn_(&fe);
      ro = rn;
      s /= 2.;

/*                    Compute the new approximate root
.                * */

      rn = ro + sro * srn * s;

/*                    Specify another stopping test  u
sing the differ- * */
/*                    ence between two successive iter
ations.          * */

      if ((d__1 = rn - ro, abs(d__1)) <= *delta) {
          goto L70;
      }
/* L60: */
        }
L70:
        rt = rn;
        if (er - rt <= *delta) {
      goto L100;
        }

/*                 Construct a  Characteristic  N-Polyhedr
on using the * */
/*                 computed roots.                        
             * */

        for (k = -1; k <= 1; k += 2) {
      x[j1] = rt - (doublereal) k * dstar;
//cout << rt << " " << k << " " << dstar << endl;
                        i__5 = *n;
      for (l = 1; l <= i__5; ++l) {
//cout << "3" << endl;
                            fx[l] = fnc(&x[1], l);
          wa1[l] = sgn_(&fx[l]);
/* L80: */
      }
      test = lsinnc_(n, &fx[1], epsilo);
      if (test) {
          goto L270;
      }
      i__ = iscomc_(n, &wa1[1], &c__[c_offset], nv);
      if (i__ != 0) {
          dcopy_(n, &x[1], &c__1, &cp[i__ + cp_dim1], nv);
          a[i__] = 0.;

/*                       Determine if  a  Characte
ristic  N-Polyhedron * */
/*                       has been completed.      
                     * */

          test = lsinnc_(nv, &a[1], epsmch);
          if (test) {
        goto L150;
          }
      }
/* L90: */
        }
L100:
        ;
    }
/* L110: */
      }
/* L120: */
  }
/* L130: */
    }

/*     Complete the  Characteristic  N-Polyhedron vertices, using  the * 
*/
/*     Initial  N-Polyhedron vertices.  Note that this action is taken * 
*/
/*     when  the  Characteristic  N-Polyhedron is incomplete,  and the * 
*/
/*     the subroutine  GENBIS  is utilized.                            * 
*/

    if (*icon != 1) {
  goto L260;
    }
    saverc_(n, &po[po_offset], nv, &x[1]);
    i__1 = *nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
  l = 1 - i__ + *nv;
  if (a[i__] > 0.) {
      if (a[l] == 0.) {
    sc1 = 2.;
    sc2 = -1.;
    sscalc_(n, &sc1, &x[1], &c__1, &c__1, &sc2, &cp[cp_offset], &
      l, nv, &cp[cp_offset], &i__, nv);
      }
  }
/* L140: */
    }
    goto L260;

/*     Reconstruct the  Characteristic N-Polyhedron in such a way that * 
*/
/*     its vertices are the  vertices of  a scaled translation  of the * 
*/
/*     unit N-cube.                                                    * 
*/

L150:
    *inf1 = 1;
    nnv = *n * *nv;
    dcopy_(&nnv, &cp[cp_dim1 + 1], &c__1, &wm[wm_dim1 + 1], &c__1);
    smnlgc_(n, &cp[cp_offset], nv, &wa1[1], &wa2[1]);
    i__1 = *nv;
    for (i1 = 1; i1 <= i__1; ++i1) {
  i__2 = *n;
  for (j = 1; j <= i__2; ++j) {
      x[j] = wa1[j] + b[i1 + j * b_dim1] * wa2[j];
/* L160: */
  }
  k = iscomc_(n, &x[1], &wm[wm_offset], nv);
  if (k != 0) {
      goto L190;
  }
  k = iscomc_(n, &x[1], &po[po_offset], nv);
  if (k != 0) {
      goto L190;
  }
  i__2 = *n;
  for (j = 1; j <= i__2; ++j) {
            fx[j] = fnc(&x[1], j);
/* L170: */
  }
  test = lsinnc_(n, &fx[1], epsilo);
  if (test) {
      goto L270;
  }

/*        Change vertices.                                            
 * */

  i__2 = *n;
  for (j = 1; j <= i__2; ++j) {
      fx[j] = d_sign(&c_b41, &fx[j]);
/* L180: */
  }
  i__ = iscomc_(n, &fx[1], &c__[c_offset], nv);
  if (i__ == 0) {
      goto L190;
  }
  dcopy_(n, &x[1], &c__1, &cp[i__ + cp_dim1], nv);
L190:
  ;
    }

/*     Reconstruct the Characteristic  N-Polyhedron in such a way that * 
*/
/*     not more  than two of its vertices are located on the same edge * 
*/
/*     of the  Initial  N-Polyhedron.                                  * 
*/

    i__1 = *ndg;
    for (i1 = 1; i1 <= i__1; ++i1) {
  i2 = 1 - i1 + *nv;
  i__2 = *n;
  for (k = 1; k <= i__2; ++k) {
      if ((d__1 = cp[i1 + k * cp_dim1] - cp[i2 + k * cp_dim1], abs(d__1)
        ) < *epsmch) {
    goto L210;
      }
/* L200: */
  }
  goto L240;
L210:
  *inf1 = 3;
  for (i3 = 1; i3 <= 2; ++i3) {
      i4 = i1;
      if (i3 == 2) {
    i4 = i2;
      }
      dcopy_(n, &cp[i4 + cp_dim1], nv, &x[1], &c__1);
      if ((d__1 = x[k] - xo[k], abs(d__1)) < *epsmch) {
    x[k] = xo[k] + h__[k];
      } else {
    x[k] -= h__[k];
      }
      i__2 = *n;
      for (j = 1; j <= i__2; ++j) {
                fx[j] = fnc(&x[1], j);
    wa1[j] = d_sign(&c_b41, &fx[j]);
/* L220: */
      }
      test = lsinnc_(n, &fx[1], epsilo);
      if (test) {
    goto L270;
      }

/*           Change vertices.                                     
     * */

      i__ = iscomc_(n, &wa1[1], &c__[c_offset], nv);
      if (i__ == 0) {
    goto L240;
      }
      if (i4 == i__) {
    *inf1 = 1;
    dcopy_(n, &x[1], &c__1, &cp[i__ + cp_dim1], nv);
    goto L240;
      }
/* L230: */
  }
L240:
  ;
    }

    return 0;
L250:
    *inf1 = 1;
    return 0;
L260:
    *inf1 = 2;
    return 0;
L270:
    *inf1 = 4;
    return 0;

/*     Last statement of the subroutine CHAPOL.                        * 
*/

} /* chapol_ */

/* ---------------------------------------------------------------------* */

/* Subroutine */ int genbis_(const DL2D &fnc, integer *n, integer *nv, integer *ns,
  integer *ndg, doublereal *b, doublereal *c__, doublereal *cp, 
  doublereal *epsmch, doublereal *epsilo, doublereal *cpr, doublereal *
        dm, doublereal *dg, doublereal *bp, doublereal *fvb, integer *inf2,
  doublereal *a, doublereal *arp, doublereal *wm, doublereal *wa1, 
  doublereal *wa2, doublereal *wa3)
{
    /* System generated locals */
    integer b_dim1, b_offset, c_dim1, c_offset, cp_dim1, cp_offset, cpr_dim1, 
      cpr_offset, wm_dim1, wm_offset, i__1, i__2, i__3, i__4, i__5, 
      i__6;

    /* Local variables */
    doublereal zeta;
    integer nirp;
    doublereal ritr;
    logical test;
    integer i__, j, k, m, i1, j1, i3, j2, j3, j4, j5, i2;
    doublereal sc1, sc2, gdg, gdm;
    integer itr, nnv, irp1, irp2;

/* ********************************************************************** 
*/
/*                                                                     * 
*/
/*       EPSMCH  is a positive input variable which determines the ma- * 
*/
/*         chine precision.                                            * 
*/
/*                                                                     * 
*/
/*       EPSILO  is a nonnegative  input  variable. Termination occurs * 
*/
/*         when the algorithm estimates that the Infinity Norm, of the * 
*/
/*         function  values  at an  approximate  solution, is at  most * 
*/
/*         EPSILO. If EPSILO is less than the machine precision EPSMCH,* 
*/
/*         EPSILO  becomes equal to  EPSMCH.                           * 
*/
/*                                                                     * 
*/
/*       CPR is an  output  NV by N  matrix, with entries  in each row * 
*/
/*         which are the corresponding elements of  CP, and determines * 
*/
/*         the Characteristic  N-Polyhedron used for the iterative re- * 
*/
/*         finement.                                                   * 
*/
/*                                                                     * 
*/
/*       DM is an output array of length NS which contains the lengths * 
*/
/*         of all the proper 1-simplexes of the Characteristic N-Poly- * 
*/
/*         hedron.                                                     * 
*/
/*                                                                     * 
*/
/*       DG is an output array of length NDG which contains the lengths* 
*/
/*         of all the diagonals of the Characteristic N-Polyhedron.    * 
*/
/*                                                                     * 
*/
/*       BP is an output array of length N which determines an approxi-* 
*/
/*         mate solution.                                              * 
*/
/*                                                                     * 
*/
/*       FVB is an output array of length N which determines the func- * 
*/
/*         tion values at  BP.                                         * 
*/
/*                                                                     * 
*/
/*       INF2 is an integer output variable set as follows :           * 
*/
/*                                                                     * 
*/
/*         INF2 = 0   GENBIS has not been evoked.                      * 
*/
/*                                                                     * 
*/
/*         INF2 = 1   The solution has been found  within the required * 
*/
/*                    accuracy of  EPSILO.                             * 
*/
/*                                                                     * 
*/
/*         INF2 = 2   Iterations have reached their upper limit.       *
*/
/*                    The answer may not be accurate.                  * 
*/
/*                                                                     * 
*/
/*         INF2 = 3   The  length  of  the longest  diagonal  has been * 
*/
/*                    found to be less than   2.0 * N * EPSILO.        * 
*/
/*                                                                     * 
*/
/*       A is an  output  array of length  NV  which  determines which * 
*/
/*         vertices of the Characteristic N-Polyhedron are substituted * 
*/
/*         during the bisection of its proper 1-simplexes.             * 
*/
/*                                                                     * 
*/
/*       ARP is an output array of length  N  which is used during the * 
*/
/*         relaxation procedure.                                       * 
*/
/*                                                                     * 
*/
/*       WM  is an  NV by N  work matrix.                              * 
*/
/*                                                                     * 
*/
/*       WA1, WA2, WA3  are work arrays of length N.                   * 
*/
/*                                                                     * 
*/
/*                                                                     * 
*/
/*     Subprograms required :                                          * 
*/
/*                                                                     * 
*/
/*       USER-Supplied ...... FNC                                      * 
*/
/*                                                                     * 
*/
/*       BLAS-Supplied ...... SCOPY, SNRM2, ISAMAX                     * 
*/
/*                                                                     * 
*/
/*       CHABIS-Supplied .... SSCALC, LSINNC, ISCOMC, SMNLGC           * 
*/
/*                                                                     * 
*/
/*       FORTRAN-Supplied ... FLOAT, LOG, INT, SIGN                    * 
*/
/*                                                                     * 
*/
/*                                                                     * 
*/
/*                                        CHABIS. Version of June 1988 * 
*/
/*                                        Michael N. Vrahatis          * 
*/
/*                                                                     * 
*/
/* ********************************************************************** 
*/

/*     Set  ZETA, a small positive real variable, which determines the *
*/
/*     accuracy of a stopping criterion during  the bisection process. * 
*/
/*     The algorithm terminates  if it is determined  that the longest * 
*/
/*     diagonal  of a  refined  Characteristic  N-Polyhehedron  has  a * 
*/
/*     length less than  ZETA.                                         * 
*/

    /* Parameter adjustments */
    --wa3;
    --wa2;
    --wa1;
    --arp;
    --fvb;
    --bp;
    wm_dim1 = *nv;
    wm_offset = wm_dim1 + 1;
    wm -= wm_offset;
    --a;
    cpr_dim1 = *nv;
    cpr_offset = cpr_dim1 + 1;
    cpr -= cpr_offset;
    cp_dim1 = *nv;
    cp_offset = cp_dim1 + 1;
    cp -= cp_offset;
    c_dim1 = *nv;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;
    b_dim1 = *nv;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    --dm;
    --dg;

    /* Function Body */
    zeta = (doublereal) (*n) * 2. * *epsilo;

/*     Set the maximum number of iterations of the relaxation process. * 
*/

//    nirp = 2;
// EBR
    nirp = 1;

/*     Set all the entries of the matrix  CPR equal to the correspond- * 
*/
/*     ing entries of  CP.                                             * 
*/

    nnv = *n * *nv;
    dcopy_(&nnv, &cp[cp_dim1 + 1], &c__1, &cpr[cpr_dim1 + 1], &c__1);

/*     Compute the lengths of all  the proper 1-simplexes of the Char- * 
*/
/*     acteristic N-Polyhedron and store them in the array, DM.        * 
*/

    i3 = 0;
    i__1 = *n;
    for (j1 = 1; j1 <= i__1; ++j1) {
  j2 = pow_ii(&c__2, &j1) - 2;
  i__2 = *n - j1;
  j3 = pow_ii(&c__2, &i__2);
  i__2 = j2;
  for (j4 = 0; j4 <= i__2; j4 += 2) {
      i__3 = j3;
      for (j5 = 1; j5 <= i__3; ++j5) {
    ++i3;
    i1 = j4 * j3 + j5;
    i2 = i1 + j3;
    sc1 = 1.;
    sc2 = -1.;
    sscalc_(n, &sc1, &cpr[cpr_offset], &i1, nv, &sc2, &cpr[
      cpr_offset], &i2, nv, &wa1[1], &c__1, &c__1);
    dm[i3] = dnrm2_(n, &wa1[1], &c__1);
/* L10: */
      }
/* L20: */
  }
/* L30: */
    }

/*     Compute the length of the longest proper 1-simplexes.           * 
*/

    i__ = idamax_(ns, &dm[1], &c__1);
    gdm = dm[i__];

/*     Compute the maximum number of iterations.                       * 
*/

    ritr = log(gdm * 2. / ((doublereal) (*n) * *epsilo)) / log(2.);
    k = (integer) ritr;
    itr = k + 1;
    if (ritr == (doublereal) k) {
  --itr;
    }

//EBR
//  cout << "TTT" << itr << endl;
   itr = 3;
/*     Begin the iterations.                                           * 
*/

    irp1 = 0;
    i__1 = itr;
    for (i3 = 1; i3 <= i__1; ++i3) {

/*        Bisect the diagonals.                                       
 * */

  if (irp1 != 0) {
      goto L80;
  }
  i__2 = *ndg;
  for (i1 = 1; i1 <= i__2; ++i1) {
      i2 = 1 - i1 + *nv;
//EBR
  size_t internal_counter = 0;
L40:
      sc1 = .5;
      sscalc_(n, &sc1, &cpr[cpr_offset], &i1, nv, &sc1, &cpr[cpr_offset]
        , &i2, nv, &bp[1], &c__1, &c__1);
      i__3 = *n;
//cout << "1 " << bp[1] << " " << bp[2] << " " << bp[3] << endl;
      for (j = 1; j <= i__3; ++j) {
                fvb[j] = fnc(&bp[1], j);
    wa1[j] = d_sign(&c_b41, &fvb[j]);
/* L50: */
      }

/*           Stopping test.                                       
     * */
/*           The solution has been found within the required accur
acy. * */

      test = lsinnc_(n, &fvb[1], epsilo);
      if (test) {
    goto L270;
      }

/*           Change vertices.                                     
     * */

      i__ = iscomc_(n, &wa1[1], &c__[c_offset], nv);
      if (i__ == 0) {
    goto L60;
      }
      dcopy_(n, &bp[1], &c__1, &cpr[i__ + cpr_dim1], nv);
//            if (i__ == i1 || i__ == i2) {
//EBR
            if ((i__ == i1 || i__ == i2) && internal_counter < 3) {
//EBR
  ++internal_counter;
                goto L40;
      }
L60:
      ;
  }

/*        Compute the length of the longest diagonal, GDG.            
 * */

  i__2 = *ndg;
  for (i1 = 1; i1 <= i__2; ++i1) {
      i2 = 1 - i1 + *nv;
      sc1 = 1.;
      sc2 = -1.;
      sscalc_(n, &sc1, &cpr[cpr_offset], &i1, nv, &sc2, &cpr[cpr_offset]
        , &i2, nv, &wa1[1], &c__1, &c__1);
      dg[i1] = dnrm2_(n, &wa1[1], &c__1);
/* L70: */
  }
  i__ = idamax_(ndg, &dg[1], &c__1);
  gdg = dg[i__];

/*        Determine if  GDG  is less than  ZETA.                      
 * */

  if (gdg < zeta) {
      *inf2 = 3;
      goto L250;
  }

L80:
  irp1 = 0;

/*        Construct the indexing array, A.                            
 * */

  i__2 = *nv;
  for (i__ = 1; i__ <= i__2; ++i__) {
      a[i__] = (doublereal) i__;
/* L90: */
  }

/*        Bisect the proper 1-simplexes.                              
 * */

  i__2 = *n;
  for (j1 = 1; j1 <= i__2; ++j1) {
      j2 = pow_ii(&c__2, &j1) - 2;
      i__3 = *n - j1;
      j3 = pow_ii(&c__2, &i__3);
      i__3 = j2;
      for (j4 = 0; j4 <= i__3; j4 += 2) {
    i__4 = j3;
    for (j5 = 1; j5 <= i__4; ++j5) {
        i1 = j4 * j3 + j5;
        i2 = i1 + j3;
        irp2 = 0;
        sc1 = .5;
        sscalc_(n, &sc1, &cpr[cpr_offset], &i1, nv, &sc1, &cpr[
          cpr_offset], &i2, nv, &bp[1], &c__1, &c__1);
        goto L110;

/*                 Execute the relaxation procedure.      
             * */

L100:
        ++irp2;
        irp1 = 1;
        sc1 = 2.;
        sc2 = -1.;
                    sscalc_(n, &sc1, &cpr[cpr_offset], &i__, nv, &sc2, &arp[1]
          , &c__1, &c__1, &bp[1], &c__1, &c__1);
L110:
        i__5 = *n;
//cout << "2 " << bp[1] << " " << bp[2] << " " << bp[3] << endl;
                    for (j = 1; j <= i__5; ++j) {
      wa3[j] = 0.;
                        fvb[j] = fnc(&bp[1], j);
      if (fvb[j] == 0.) {
          wa3[j] = (doublereal) j;
      }
      wa1[j] = d_sign(&c_b41, &fvb[j]);
/* L120: */
        }

/*                 Stopping test.                         
             * */
/*                 The  solution has been  found within  t
he  required * */
/*                 accuracy.                              
             * */

        test = lsinnc_(n, &fvb[1], epsilo);
        if (test) {
      goto L270;
        }

/*                 Change vertices.                       
             * */

        i__ = iscomc_(n, &wa1[1], &c__[c_offset], nv);
        if (i__ == 0) {
      goto L150;
        }
        i__5 = *n;
        for (k = 1; k <= i__5; ++k) {
      j = (integer) wa3[k];
      if (wa3[k] == 0.) {
          goto L130;
      }
      i__6 = *n - j;
      m = i__ - (integer) (c__[i__ + j * c_dim1] * pow_di(&
        c_b96, &i__6));
      sc1 = 1.;
      sc2 = -1.;
      sscalc_(n, &sc1, &cpr[cpr_offset], &i__, nv, &sc2, &
        cpr[cpr_offset], &m, nv, &wa2[1], &c__1, &
        c__1);
      test = lsinnc_(n, &wa2[1], epsmch);
      if (test) {
          goto L140;
      }
L130:
      ;
        }
        dcopy_(n, &cpr[i__ + cpr_dim1], nv, &arp[1], &c__1);
        dcopy_(n, &bp[1], &c__1, &cpr[i__ + cpr_dim1], nv);
        a[i__] = 0.;
L140:
                    if (i__ != i1 && i__ != i2 && irp2 < nirp) {
      goto L100;
        }
L150:
        ;
    }
/* L160: */
      }
/* L170: */
  }

/*        Compute the length of the longest diagonal, GDG.            
 * */

  i__2 = *ndg;
  for (i1 = 1; i1 <= i__2; ++i1) {
      i2 = 1 - i1 + *nv;
      sc1 = 1.;
      sc2 = -1.;
      sscalc_(n, &sc1, &cpr[cpr_offset], &i1, nv, &sc2, &cpr[cpr_offset]
        , &i2, nv, &wa1[1], &c__1, &c__1);
      dg[i1] = dnrm2_(n, &wa1[1], &c__1);
/* L180: */
  }
  i__ = idamax_(ndg, &dg[1], &c__1);
  gdg = dg[i__];

/*        Determine if  GDG  is less than  ZETA.                      
 * */

  if (gdg < zeta) {
      *inf2 = 3;
      goto L250;
  }

  k = 0;
  i__2 = *nv;
  for (i__ = 1; i__ <= i__2; ++i__) {
      if (a[i__] != 0.) {
    ++k;
      }
/* L190: */
  }
  if (k == 0 || irp1 == 0) {
      irp1 = 0;
      goto L230;
  }

/*        Reconstruct the refined Characteristic  N-Polyhedron in such
 * */
/*        a way that its vertices are  the vertices of a scaled trans-
 * */
/*        lation of the unit N-cube.                                  
 * */

  dcopy_(&nnv, &cpr[cpr_dim1 + 1], &c__1, &wm[wm_dim1 + 1], &c__1);
  smnlgc_(n, &cpr[cpr_offset], nv, &wa1[1], &wa2[1]);
  i__2 = *nv;
  for (i1 = 1; i1 <= i__2; ++i1) {
      i__3 = *n;
      for (j = 1; j <= i__3; ++j) {
    bp[j] = wa1[j] + b[i1 + j * b_dim1] * wa2[j];
/* L200: */
      }
      i__ = iscomc_(n, &bp[1], &wm[wm_offset], nv);
      if (i__ != 0) {
    goto L220;
      }
      i__3 = *n;
//cout << "3 " << bp[1] << " " << bp[2] << " " << bp[3] << endl;
            for (j = 1; j <= i__3; ++j) {
                fvb[j] = fnc(&bp[1], j);
    wa3[j] = d_sign(&c_b41, &fvb[j]);
/* L210: */
      }

/*           Stopping test.                                       
     * */
/*           The solution has been found within the required accur
acy. * */

      test = lsinnc_(n, &fvb[1], epsilo);
      if (test) {
    goto L270;
      }

/*           Change vertices.                                     
     * */

      i__ = iscomc_(n, &wa3[1], &c__[c_offset], nv);
      if (i__ == 0) {
    goto L220;
      }
      dcopy_(n, &bp[1], &c__1, &cpr[i__ + cpr_dim1], nv);
L220:
      ;
  }
L230:
  ;
    }

    *inf2 = 2;

/*     Find the longest diagonal of the refined Characteristic N-Poly- * 
*/
/*     hedron and take its midpoint as the final approximate solution. * 
*/

    i__1 = *ndg;
    for (i1 = 1; i1 <= i__1; ++i1) {
  i2 = 1 - i1 + *nv;
  sc1 = 1.;
  sc2 = -1.;
  sscalc_(n, &sc1, &cpr[cpr_offset], &i1, nv, &sc2, &cpr[cpr_offset], &
    i2, nv, &wa1[1], &c__1, &c__1);
  dg[i1] = dnrm2_(n, &wa1[1], &c__1);
/* L240: */
    }
    i__ = idamax_(ndg, &dg[1], &c__1);
L250:
    i1 = 1 - i__ + *nv;
    sc1 = .5;
    sscalc_(n, &sc1, &cpr[cpr_offset], &i__, nv, &sc1, &cpr[cpr_offset], &i1, 
      nv, &bp[1], &c__1, &c__1);
    i__1 = *n;
//cout << "4 " << bp[1] << " " << bp[2] << " " << bp[3] << endl;
    for (j = 1; j <= i__1; ++j) {
        fvb[j] = fnc(&bp[1], j);
/* L260: */
    }

/*     Stopping test.                                                  * 
*/
/*     The solution has been found within the required accuracy.       * 
*/

    test = lsinnc_(n, &fvb[1], epsilo);
    if (test) {
  goto L270;
    }

    return 0;
L270:
    *inf2 = 1;
    return 0;

/*     Last statement of the subroutine GENBIS.                        * 
*/

} /* genbis_ */

/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */
/*                                                                     * */
/*             B L A S  -  L I K E     U T I L I T I E S               * */
/*                                                                     * */
/*                   CHABIS.  Version of June 1988                     * */
/*                                                                     * */
/*                        Michael N. Vrahatis                          * */
/*                                                                     * */
/* ---------------------------------------------------------------------* */
/* ---------------------------------------------------------------------* */

logical lsinnc_(integer *n, doublereal *sx, doublereal *epsilo)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    logical ret_val;

    /* Local variables */
    integer i__;


/*     Determine if  the Infinity Norm of a single precision  N-vector * 
*/
/*     stored in SX( ) is less than, or equal to, the single precision * 
*/
/*     positive constant EPSILO; if so, set  LSINNC  equal to " TRUE ";* 
*/
/*     otherwise, set  LSINNC equal to " FALSE ".                      * 
*/


    /* Parameter adjustments */
    --sx;

    /* Function Body */
    ret_val = FALSE_;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
  if ((d__1 = sx[i__], abs(d__1)) > *epsilo) {
      return ret_val;
  }
/* L10: */
    }
    ret_val = TRUE_;

    return ret_val;
} /* lsinnc_ */

/* ---------------------------------------------------------------------* */

integer iscomc_(integer *n, doublereal *sx, doublereal *sy, integer *incy)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;

    /* Local variables */
    integer i__, j, k;


/*     Find the value of I, where I can be I = 1,2,...,INCY, such that * 
*/
/*     the entries of a single precision N-vector, stored in  SX( J ), * 
*/
/*     J = 1,2,...,N, coincide with the corresponding single precision * 
*/
/*     entries  SY( I+(J-1)*INCY )  for all J = 1,2,...,N. If there is * 
*/
/*     no such value, then return the value of zero.                   * 
*/


    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    i__1 = *incy;
    for (i__ = 1; i__ <= i__1; ++i__) {
  k = i__;
  i__2 = *n;
  for (j = 1; j <= i__2; ++j) {
      if (sx[j] != sy[k]) {
    goto L20;
      }
      k += *incy;
/* L10: */
  }
  ret_val = i__;
  return ret_val;
L20:
  ;
    }
    ret_val = 0;

    return ret_val;
} /* iscomc_ */

/* ---------------------------------------------------------------------* */

doublereal sgn_(doublereal *sv)
{
    /* System generated locals */
    doublereal ret_val;


/*     Transfer the sign of the single precision variable SV in such a * 
*/
/*     way that if  SV is equal to zero, the sign of  SV will be equal * 
*/
/*     to zero.                                                        * 
*/


    ret_val = d_sign(&c_b41, sv);
    if (*sv == 0.) {
  ret_val = 0.;
    }

    return ret_val;
} /* sgn_ */

/* ---------------------------------------------------------------------* */

/* Subroutine */ int saverc_(integer *n, doublereal *sx, integer *incx, 
  doublereal *sy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, j, k, m, k1, k2;
    doublereal sum;


/*     Find the average value for each of the  N vector segments (each * 
*/
/*     of which contains  INCX  elements)  of the single precision  N- * 
*/
/*     vector stored in SX( ). Store these values in the single preci- * 
*/
/*     sion  N-vector  SY( ).                                          * 
*/


    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    m = *incx % 4;
    k = -(*incx);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
  k += *incx;
  sum = 0.;
  if (m == 0) {
      goto L20;
  }
  k1 = k + 1;
  k2 = k + m;
  i__2 = k2;
  for (i__ = k1; i__ <= i__2; ++i__) {
      sum += sx[i__];
/* L10: */
  }
  if (*incx < 4) {
      goto L40;
  }
L20:
  k1 = k + m + 1;
  k2 = k + *incx;
  i__2 = k2;
  for (i__ = k1; i__ <= i__2; i__ += 4) {
      sum = sum + sx[i__] + sx[i__ + 1] + sx[i__ + 2] + sx[i__ + 3];
/* L30: */
  }
L40:
  sy[j] = sum / (doublereal) (*incx);
/* L50: */
    }

    return 0;
} /* saverc_ */

/* ---------------------------------------------------------------------* */

/* Subroutine */ int sscalc_(integer *n, doublereal *v1, doublereal *sx, 
  integer *lx, integer *incx, doublereal *v2, doublereal *sy, integer *
  ly, integer *incy, doublereal *sz, integer *lz, integer *incz)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, k, ix, iy, iz;


/*     Multiply the single precision variable  V1 by the single preci- * 
*/
/*     sion N-vector stored in  SX( ), with leading dimension  LX  and * 
*/
/*     storage increment  INCX ; to the product, add the single preci- * 
*/
/*     sion variable  V2 multiplied by  the single precision  N-vector * 
*/
/*     stored in SY( ),with leading dimension LY and storage increment * 
*/
/*     INCY. Store the result in the single precision N-vector  SZ( ), * 
*/
/*     with leading dimension  LZ  and storage increment  INCZ.        * 
*/


    /* Parameter adjustments */
    --sz;
    --sy;
    --sx;

    /* Function Body */
    k = *n % 4;
    ix = *lx;
    iy = *ly;
    iz = *lz;
    if (k == 0) {
  goto L20;
    }
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
  sz[iz] = *v1 * sx[ix] + *v2 * sy[iy];
  ix += *incx;
  iy += *incy;
  iz += *incz;
/* L10: */
    }
    if (*n < 4) {
  return 0;
    }
L20:
    i__1 = *n - 1;
    for (i__ = k; i__ <= i__1; i__ += 4) {
  sz[iz] = *v1 * sx[ix] + *v2 * sy[iy];
  ix += *incx;
  iy += *incy;
  iz += *incz;
  sz[iz] = *v1 * sx[ix] + *v2 * sy[iy];
  ix += *incx;
  iy += *incy;
  iz += *incz;
  sz[iz] = *v1 * sx[ix] + *v2 * sy[iy];
  ix += *incx;
  iy += *incy;
  iz += *incz;
  sz[iz] = *v1 * sx[ix] + *v2 * sy[iy];
  ix += *incx;
  iy += *incy;
  iz += *incz;
/* L30: */
    }

    return 0;
} /* sscalc_ */

/* ---------------------------------------------------------------------* */

/* Subroutine */ int smnlgc_(integer *n, doublereal *sx, integer *incx, 
  doublereal *smn, doublereal *slg)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    doublereal smin, smax;
    integer i__, j, k;
    doublereal x;

/*                                                                     * 
*/
/*     Find the maximum value, and the difference between  the maximum * 
*/
/*     and minimum values, within each of the  N vector segments (each * 
*/
/*     of which contains  INCX  elements) of  the single precision  N- * 
*/
/*     vector stored in  SX( ). Store  these values in the single pre- * 
*/
/*     cision  N-vectors  SMN( ) and  SLG( ),  respectively.           * 
*/


    /* Parameter adjustments */
    --slg;
    --smn;
    --sx;

    /* Function Body */
    k = -(*incx);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
  k += *incx;
  smin = sx[k + 1];
  smax = smin;
  i__2 = k + *incx;
  for (i__ = k + 2; i__ <= i__2; ++i__) {
      x = sx[i__];
      smin = min(smin,x);
      smax = max(smax,x);
/* L10: */
  }
  smn[j] = smin;
  slg[j] = smax - smin;
/* L20: */
    }

    return 0;
} /* smnlgc_ */

