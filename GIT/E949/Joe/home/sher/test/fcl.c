#include<math.h>
#define PI          3.14159265358979323846
#define LIKELIHOODRATIO
/* #define PBOCK */
/* #define CPRIME */
#define ESTBIG 1.0E10
#define NCOMBMAX 1.0E15

/* An explicit probability sum for computing confidence levels.
   Uses ecl.f for estimators and a Poisson probability subroutine.
   Note the underscores at the end of externally visible symbols.
   These are there because +ppu is needed to compile fortran code
   on the HP when linking with CERNLIB.  Take them off if you
   are building on another platform or not using the +ppu option
   of f77.
   Warning -- this sum can be either very fast or very slow.
   Also -- it may not be protected against all strange inputs.
   (don't try negative expected s or b, for instance.)
*/

/*flag = 0 -- do only CLB -- flag=1 do cls, clsb too */

/* Function prototypes and constants for fcl.c */

void fcl_(int *n, double *s, double *b, int *d,
         double *cls, double *clsb, double *clb, int *flag);

#define MAXCHANS 400
#define MAXOUT 100

void fcl_estlist(int n, double *v, double *s, double *b, double dest,
                 double estlist[MAXCHANS][MAXOUT],
                 double plist[MAXCHANS][MAXOUT],
                 int nplist[MAXCHANS],int *sizelist,
		 double *ncombs);

double fcl_compute(int n, double dest, 
                   double estlist[MAXCHANS][MAXOUT],
                   double plist[MAXCHANS][MAXOUT],
                   int nplist[MAXCHANS]);

void e_estimf_(double *s, double *b, int *d, double *estim);

void e_destnf_(int *n, double *s, double *b, int *d, double *dest);

void e_poissf_(double *u, int *n, double *prob);

double e_factf_(int *n);

double e_f2f_(int *n);

/*------subroutines---------*/

void fcl_(int *ni, double *si, double *bi, int *di,
         double *cls, double *clsb, double *clb, int *flag)
{
  double estlist[MAXCHANS][MAXOUT];
  double plist[MAXCHANS][MAXOUT];
  double v[MAXCHANS];
  int nplist[MAXCHANS];
  double dest;
  int i,sizelist;
  int k, *d;
  double *s, *b;
  double ncombs;
  double *vp;
  int n;

  n = *ni;
  s = si;
  b = bi;
  d = di;

  if (n>MAXCHANS) 
    { printf("fcl: Too many channels %d (max: %d)\n",n,MAXCHANS);
      printf("fcl: exiting. \n");
      exit(0);}

  /* chop down the problem size a little bit to keep the number
     of combinations within reason.  Assume the best channels
     are saved for last */

  for (i=0;i<n;i++)
    { v[i] = s[i] + b[i]; }
  vp = v;
  vp--;
  s--; b--; d--; n++;
  do 
    {
       n--;
       e_destnf_(&n,++s,++b,++d,&dest);
       vp++;
       fcl_estlist(n,v,s,b,dest,estlist,plist,nplist,&sizelist,&ncombs);
    }
  while(ncombs > NCOMBMAX);

  if (*flag != 0) 
    {
      for (i=0;i<n;i++)
        { v[i] = s[i] + b[i]; }
      fcl_estlist(n,v,s,b,dest,estlist,plist,nplist,&sizelist,&ncombs);
      if (sizelist == 0) 
        { *clsb = 1.0; }
      else
        { *clsb = fcl_compute(sizelist,dest,estlist,plist,nplist); }
    }

  for (i=0;i<n;i++)
    { v[i] = b[i]; }
  fcl_estlist(n,v,s,b,dest,estlist,plist,nplist,&sizelist,&ncombs);
  if (sizelist == 0) 
    { *clb = 1.0; }
  else
    { *clb = fcl_compute(sizelist,dest,estlist,plist,nplist); }

  if (*flag != 0) 
    {
      if (*clb >= 0.0)
        { *cls = *clsb / *clb; }
      else
        { *cls = 2.0; }
    }
}

/*-------------------------------------------------------------------*/
/* Make a list of possible outcomes of each channel such that the
   estimator is <= the data estimator.  Give up on a channel if it
   takes more than MAXOUT to do it -- reason -- if you have a candidate
   in a channel with very low b (and some s) in the data, then to get
   the same estimator out of a channel with very high b and not much s
   in the data you need lots of events.  The probability of having 
   high estimators in that channel is nil, so the cumulative prob.
   of having less than that is very near 1, and you can neglect the
   high-background channel in this case.
*/

void fcl_estlist(int n, double *v, double *s, double *b, double dest,
                 double estlist[MAXCHANS][MAXOUT],
                 double plist[MAXCHANS][MAXOUT],
                 int nplist[MAXCHANS], int *sizelist,
                 double *ncombs)
{
  int i,j,k;
  double psum;
  *ncombs = 1;
  *sizelist = 0;
  for (i=0;i<n;i++)
    {
      psum = 0.0;
      for (j=0;j<MAXOUT;j++)
	{
	  e_estimf_(&(s[i]),&(b[i]),&j,&(estlist[*sizelist][j]));
	  e_poissf_(&(v[i]),&j,&(plist[*sizelist][j]));
	  psum += plist[*sizelist][j];
	  nplist[*sizelist] = j+1;
          if (estlist[*sizelist][j] > dest) break;
	}
      if (nplist[*sizelist] < MAXOUT) 
	{ (*sizelist)++; *ncombs *= j; }
      if (psum > 1.0 && *sizelist > 0)
	{ k = *sizelist - 1;
          for (j=0;j<nplist[k];j++)
	    { plist[k][j] /= psum; }
        }
    }
}

/*-------------------------------------------------------------------*/

double fcl_compute(int n, double dest, 
                   double estlist[MAXCHANS][MAXOUT],
                   double plist[MAXCHANS][MAXOUT],
                   int nplist[MAXCHANS])
{
  int i,j;
  double psum = 0.0;

  n--;
  if (n==0)
    {
      for (i=0;i<nplist[n];i++)
	{ if (dest<estlist[n][i]) break;
          psum += plist[n][i];
        }
    }
  else
    {
      for (i=0;i<nplist[n];i++)
        {
	  if (dest<estlist[n][i]) break;
	  psum += plist[n][i] *
	    fcl_compute(n,dest-estlist[n][i],
                        estlist,plist,nplist);
        }
    }
  return psum;
}


/* test statistic computer -- no syst. errors */


void e_estimf_(double *s, double *b, int *d, double *estim)
{

#ifdef PBOCK
if (*s == 0.0) 
  { *estim = 0.0; }
else
  { *estim = ((double) d)/(CPBOCK + *b/(*s)); }
#endif

#ifdef CPRIME
if (*s == 0.0) 
  { *estim = 0.0; }
else
  { *estim = ((double) *d)/(*s/2.0 + *b/(*s)); }
#endif

#ifdef LIKELIHOODRATIO

if (*b <= 0.0 && *d > 0)
  { *estim = ESTBIG; }
else if (*b == 0.0 && *s == 0.0)
  { *estim = 0.0; }
else if (*b == 0.0 && *d == 0) 
  { *estim = 0.0; }
else
  { *estim = ((double) *d)*log(1.0 + *s/(*b)); }

#endif
}

void e_destnf_(int *n, double *s, double *b, int *d, double *dest)
{
  double estim;
  int i;
  double *sl,*bl;
  int *dl;
  *dest = 0.0;
  sl = s;
  bl = b;
  dl = d;
  for(i=0;i<*n;i++)
   {
     e_estimf_(sl++,bl++,dl++,&estim);
     *dest += estim;
   }
}

void e_poissf_(double *u, int *n, double *prob)
{
/* make sure 0^0=1 in this case.  Zero expected events gives 100%
   probability to zero observed events. */
if (*n<100) 
  { 
    if (*u<1.0e-6 && *n == 0) 
      { *prob = 1.0; }
    else 
      { *prob = pow(*u,(double) *n)*exp(-*u)/e_factf_(n); }
  }
else
  { *prob = exp(-((double) *n - *u)*((double) *n - *u)/(2.0*(*u)))/
                sqrt(2.0*PI*(*u)); }
}

double e_factf_(int *n)
  {
    int i,j;
    static double fl[100];
    static int first; /* assume initialized to zero */
    if (first == 0)
      { first = 1;
        for (i=0;i<100;i++)
          { fl[i] = e_f2f_(&i); }
      }
    if (*n >= 100)
      { return(e_f2f_(n)); }
    else
      { return(fl[*n]); }
  }

double e_f2f_(int *n)
  {
    double prod;
    int i;
    prod = 1.0;
    for (i=1;i<=*n;i++)
      { prod *= (double) i; }
    return(prod);
  }
