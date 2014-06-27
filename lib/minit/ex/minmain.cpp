/*
 * This is a front end to test the linear programming package
 * minit. This algorithm was written by Das and Salazar in algol.
 *
 * Solve the LP      T
 *    max C X  subj. to
 *              AX <= B
 *               X >= 0
 */

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <minit.h>

main(int argc, char **argv)
{
  int  ier,m,n,m1,m2; /* Dimensions of various matrices */
  Cmatrix A;
  vec_double b, c, x, w;
  double z;
  float tmp;
  char format[64], *usage;
  int suppress;
  int int_part, frac_part;
  double precision;
  register int i,j,arg;
  FILE *in;

  n = m1 = m2 = 0;
  in = stdin;
  suppress = 1; /* Suppress printing interactive messages. */
  precision = 0.0;

  usage = "\
Usage:\tminit [-pprec] [-i[filename]]\n\
or:\tminit [-pprec] < filename\n\
For detailed help: minit -help\n";

  for (arg = 1; arg < argc; arg++)
  {
    if (argv[arg][0] != '-' ||
      (argv[arg][1] != 'p' && argv[arg][1] != 'i'))
    {
        if (!strcmp(argv[arg],"-help") || !strcmp(argv[arg],"help"))
      (void) system("man minit");
        else
      (void) fputs(usage, stderr);
        exit(1);
    }

    if (argv[arg][1] == 'p')
      (void) sscanf(&(argv[arg][2]),"%f",&precision);
    else  /* if argv[arg][1] == 'i') */
    {
      if (in != stdin)  /* in already assigned. */
      {
        (void) fputs("Only one input file permitted.\n",
          stderr);
        exit(1);
      }
      if((in = fopen(&(argv[arg][2]),"r")) == NULL)
      {
        (void) fprintf(stderr,
          "No input file %s; interactive mode\n",
          &(argv[arg][2]));

        in = stdin;
        suppress = 0;
      }
    }
  }

  if (!suppress)
    (void) fputs("Enter value of n (No. of unknowns): ", stderr);

  (void) fscanf(in,"%d",&n);

  if (!suppress)
    (void) fputs("Enter no. of inequality constraints: ", stderr);

  (void) fscanf(in,"%d",&m1);

  if (!suppress)
    (void) fputs("Enter no. of equality constraints: ", stderr);

  (void) fscanf(in,"%d",&m2);

  m = m1+m2;

  /* Allocate space for the arrays */
  A.resize(m, n);
  
  b.resize(m);
  c.resize(n);
  x.resize(n);
  w.resize(m);
  
  if (!suppress)
    (void) fprintf(stderr,"Enter (%d by %d) elements of A:\n",
      m,n);
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
    {
      (void) fscanf(in,"%f",&tmp);
      A[i][j] = tmp;
    }

  if (!suppress)
    (void) fprintf(stderr,
      "Enter %d elements of b, entered as a linear array:\n ",
        m);
  for(i=0;i<m;i++)
  {
    (void) fscanf(in,"%f",&tmp);
    b[i] = tmp;
  }

  if (!suppress)
    (void) fprintf(stderr,
      "Enter %d elements of c, entered as a linear array:\n",
        n);
  for(j=0;j<n;j++)
  {
    (void) fscanf(in,"%f",&tmp);
    c[j] = tmp;
  }

  if (suppress)
  {
    (void) fputs("INPUT DATA:\n",stderr);
    (void) fprintf(stderr,"No. of variables (n) = %d\n", n);
    (void) fprintf(stderr,
      "No. of inequality constraints (m1) = %d\n", m1);
    (void) fprintf(stderr,
      "No. of equality constraints (m2) = %d\n", m2);
  }

  (void) fputs("------------------------------",stderr);
  (void) fputs(" EXECUTING MINIT ",stderr);
  (void) fputs("------------------------------",stderr);
  (void) fputc('\n',stderr);

  minit doit;
  doit.zx3lp(A,b,c,n,m1,m2,z,x,w,ier);
  cout << "ier is " << ier << endl; 
  switch(ier)
  {
    case 133:
    break;

    case 1:
    (void) fprintf(stderr,"No solution.\n");
    exit(1);

    case 131:
    (void) fprintf(stderr,"Primal has unbounded solution.\n");
    exit(3);
  }

  /* Make sure precision is acceptible. */
  int_part = precision;
  frac_part = precision - int_part;
  if (int_part <= frac_part || int_part > 25)
    precision = 10.3; /* The default precision. */

  (void) sprintf(format,"\nValue of z: %%%gg\n",precision);
  (void) printf(format,z);
  (void) puts("Primal solution:");

  (void) sprintf(format,"%%%gg ",precision);
  for(j=0;j<n;j++)
    (void) printf(format,x[j]);

  (void) puts("\nDual solution:");
  for(i=0;i<m;i++)
    (void) printf(format,w[i]);

  (void) putchar('\n');
  exit(0);
}
