/*  lake.c  */
#define VERSION_INFORMATION  "svnid: $ Id: $"

/*
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
*/

/*
COMPILE (for UNIX or Linux or Cygwin):
   make 	     build and compile the program
   make test	     run a test case as an example
   more xxx.log      view the test case


CITATION
   Use this reference to cite this code in a publication.

   PhD Thesis, 1990, Northwestern University, 
   Evanston, IL, USA, Pete R. Jemian.
   http://jemian.org/pjthesis.pdf

USE
   These instructions are very limited.
   Additional help might be available in the Theory chapter
   of my PhD thesis:
   http://jemian.org/pjthesis.pdf

   Compile the program as above and run the "lake" executable.
   No "install" methods are provided.  The Makefile
   will build "lake" executable as a shared-object.  (Not static).

   Input data will be provided in an ASCII TEXT file
   as three columns (Q  I  dI) separated by white space.
   Units must be compatible.  (I and dI must have same units)
   Q: scattering vector
   I measured SAS intensity
   dI estimated uncertainties of I (usually standard deviation)
   Note that dI MUST be provided and MUST not be zero.
   
NOTES
    revised to compile on:
       Linux gcc
       Solaris Sun Workshop C
       Windows 95/98/NT Metrowerks CodeWarrior
       Macintosh Metrowerks CodeWarrior
       Macintosh Symantec C
    4 February 1994
       direct translation of LAKE.FOR
       which would not run correctly under System 7
       (Some compiler error that would not pass array DC(*)
       into subroutine SMEAR.)
    2002-Sept-18:
    	users of Windows version report errors associated with file names
	Probably a bug was edited in.


HISTORY
    lake.c was derived from the FORTRAN program:
    Lake.FOR  25 May 1991
    ref: J.A. Lake; ACTA CRYST 23 (1967) 191-194.
    by: Pete R. Jemian, Late-Nite(tm) Software

    Also see: O. Glatter; ACTA CRYST 7 (1974) 147-153
      W.E. Blass & G.W.Halsey (1981).  "Deconvolution of
        Absorption Spectra."  New York City: Academic Press
      P.A. Jansson.  (1984) "Deconvolution with Applications
        in Spectroscopy."  New York City: Academic Press.
      G.W.Halsey & W.E. Blass.  "Deconvolution Examples"
        in "Deconvolution with Applications in Spectroscopy."
        Ed. P.A. Jansson.  (see above)


PURPOSE
      This program applies the iterative desmearing technique of Lake
        to small-angle scattering data.  The way that the program works
        is that the user selects a file of data (x,y,dy) to be desmeared.
        If a file was not chosen, the program will end.  Otherwise the
        user is then asked to specify the slit-length (in the units of the
        x-axis); the X at which to begin fitting the last data points to a
        power-law of X, the output file name, and the number of iterations
        to be run.  Then the data file is opened, the data is read, and the
        data file is closed.  The program begins iterating and shows an
        indicator of progress on the screen in text format.
      It is a mistake to run this program on data that has been desmeared
        at least once (by this program) as you will see.  The problem is
        that the program expects that the input data has been smeared, NOT
        partially desmeared.  Lake's technique should be made to iterate
        with the original, smeared data and subsequent trial solutions
        of desmeared data.
      The integration technique used by this program to smear the data
        is the trapezoid-rule where the step-size is chosen by the
        spacing of the data points themselves.  A linear
        interpolation of the data is performed.  To avoid truncation
        effects, a power-law extrapolation of the intensity
        is made for all values beyond the range of available
        data.  This region is also integrated by the trapezoid
        rule.  The integration covers the region from l = 0
        up to l = lo. (see routine SMEAR).
        This technique allows the slit-length weighting function
        to be changed without regard to the limits of integration
        coded into this program.
*/



/*****************************************************************
 ***                  global definitions                       ***
 *****************************************************************/
#define LakeUnit    (1)
#define LakeFast    (2)
#define LakeChi2    (3)
#define InfItr      (10000)

/*****************************************************************
 ***                  global variables                         ***
 *****************************************************************/
    int     NumPts;     /* number of data points read in */
    double  *q;         /* q : scattering vector, horizontal axis */
    double  *E, *dE;    /* E, dE : experimental intensity and estimated error */
    double  *S;         /* S : re-smeared intensity */
    double  *C, *dC;    /* C, dC : corrected intensity and estimated error */
    double  *resid;     /* resid : normalized residuals, = (E-S)/dE */
    int     NumItr = InfItr;    /* maximum number of iterations to run */
    char    InFile[256], OutFil[256];
    double  sLengt = 1.0;  /* slit length, as defined by Lake */
    double  sFinal = 1.0;  /* to start evaluating the constants for extrapolation */
    int     mForm = 1;     /* model final data as a constant */
    int     LakeForm = 2;  /* shows the fastest convergence most times */
    double  fSlope;     /* linear coefficient of data fit */
    double  fIntercept; /* constant coefficient of data fit */

/*****************************************************************
 ***                  subroutines & functions                  ***
 *****************************************************************/
double Plengt (double x);
double FindIc (double x, double y, int NumPts, double *q, double *C);

/* from toolbox.c */
void AskString (char *question, char *answer);
void AskDouble (char *question, double *answer);
void AskInt (char *question, int *answer);
char AskYesOrNo (char *question, char standard);


/*****************************************************************
 ***                  #includes                                ***
 *****************************************************************/
#ifdef THINK_C
#include <console.h>
extern long _fcreator;
#endif
#ifdef __MWERKS__
#include <sioux.h>
#if TARGET_OS_WIN32
#include <WinSIOUX.h>
#endif
#endif
#include <math.h>
#include "recipes.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


void GetInf (
    char *InFile,
    char *OutFil,
    double *sLengt,
    double *sFinal,
    int *NumItr,
    int MaxItr,
    int *mForm,
    int *LakeForm)
        /*  Get information about the desmearing parameters.
         *  This is designed to be independent of wavelength
         *    or radiation-type (i.e. neutrons, X rays, etc.)
         */
{
    InFile[0] = (char) 0;
    AskString ("What is the input data file name? <Quit>", InFile);
    if (!InFile[0]) return;
    do
        AskString ("What is the output data file?", OutFil);
    while (!strcmp (InFile, OutFil));
    AskDouble ("What is the slit length (x-axis units)?", sLengt);
    printf (
        "Extrapolation forms to avoid truncation-error.\n"
        "   1 = flat background, I(q) = B\n"
        "   2 = linear, I(q) = b + q * m\n"
        "   3 = power law, I(q) = b * q^m\n"
        "   4 = Porod law, I(q) = Cp + Bkg / q^4\n"
    );
    AskInt ("Which form?", mForm);
    AskDouble ("What X to begin evaluating extrapolation (x-axis units)?", 
                sFinal);
    AskInt ("How many iteration(s)? (10000 = infinite)", NumItr);
    printf (
        " Weighting methods for iterative corrections:\n"
        " Correction = weight * (MeasuredI - SmearedI)\n"
        "   #1) weight = 1.0\n"
        "   #2) weight = CorrectedI / SmearedI\n"
        "   #3) weight = 2*SQRT(ChiSqr(0) / ChiSqr(i))\n"
    );
    AskInt ("Which method?", LakeForm);
}

void FixErr (int n, double *x, double *y, double *dy, double *z, double *dz)
/*
 *  Estimate the error on Z based on data point scatter and
 *  previous error values and smooth that estimate.
 */
{
    int i, j;
    double slope, intercept, zNew, w1, w2;

    /* Error proportional to input error */
    for (i = 1; i <= n; i++) 
        dz[i] = z[i] * dy[i] / y[i];

    /*
     *  Error based on scatter of desmeared data points.
     *    Determine this by fitting a line to the points
     *    i-1, i, i+1 and take the difference.  Add this to dz.
     */
    SumClr ();
    SumAdd (x[1], z[1]);
    SumAdd (x[2], z[2]);
    SumAdd (x[3], z[3]);
    SumLR (&slope, &intercept);
    dz[1] += fabs (intercept + slope*x[1] - z[1]);
    dz[2] += fabs (intercept + slope*x[2] - z[2]);
    for (i = 3; i < n; i++) {
      SumClr ();
      SumAdd (x[i-1], z[i-1]);
      SumAdd (x[i],   z[i]);
      SumAdd (x[i+1], z[i+1]);
      SumLR (&slope, &intercept);
      zNew = intercept + slope * x[i];
      dz[i] += fabs (zNew - z[i]);
    }
    dz[n] += fabs (intercept + slope*x[n] - z[n]);

    /*
     *  Smooth the error by a 3-point moving average filter.
     *    Do this 5 times.  Don't smooth the end points.
     *    Weight the data points by distance^2 (as a penalty)
     *    using the function weight(u,v)=(1 - |1 - u/v|)**2
     *    By its definition, weight(x0,x0) == 1.0.  I speed
     *    computation using this definition.  Why I can't use
     *    this definition of weight as a statement function
     *    with some compilers is beyond me!
     *  Smoothing is necessary to increase the error estimate
     *    for some grossly under-estimated errors.
     */
    for (j = 1; j <= 5; j++) {
        for (i = 2; i < n; i++) {
            w1 = SQR(1 - fabs (1 - (x[i-1]/x[i])));
            w2 = SQR(1 - fabs (1 - (x[i+1]/x[i])));
            dz[i] = (w1 * dz[i-1] + dz[i] + w2 * dz[i+1])
                    / (w1 + 1.0 + w2);
        }
    }
}

void Prep (double *x, double *y, double *dy, int NumPts)
    /*
     *  Calculate the constants for an extrapolation fit
     *  from all the data that satisfy x(i) >= sFinal.
     */
{
    double h4;
    int i;
    SumClr ();
    for (i = 1; i < NumPts; i++) {
        if (x[i] >= sFinal) {
            switch (mForm) {
                case 1:
                    SwtAdd (x[i], y[i], dy[i]);     /* weighted */
                    break;
                case 2:
                    SumAdd (x[i], y[i]);            /* un-weighted */
                    break;
                case 3:
                    SumAdd (log(x[i]), log(y[i]));  /* un-weighted */
                    break;
                case 4:
                    h4 = SQR(SQR(x[i]));
                    SwtAdd (h4, y[i]*h4, dy[i]*h4); /* weighted */
                    break;
            }
        }
    }
    switch (mForm) {
        case 1:
            MeanXY (&fSlope, &fIntercept);
            fSlope = 0.0;   /* flat background */
            break;
        case 2:
        case 3:
        case 4:
            SumLR (&fSlope, &fIntercept);
            break;
    }
}

void Smear (double *S, int NumPts, double *q, double *C, double *dC)
    /*
     *  Smear the data of C(q) into S using the slit-length
     *    weighting function "Plengt" and a power-law extrapolation
     *    of the data to avoid truncation errors.  Assume that
     *    Plengt goes to zero for l > lo (the slit length).
     *  Also assume that the slit length  function is symmetrical
     *    about l = zero.
     *  This routine is written so that if "Plengt" is changed
     *    (for example) to a Gaussian, that no further modification
     *    is necessary to the integration procedure.  That is,
     *    this routine will integrate the data out to "lo".
     */
{
    double *x, *w, hLo, hNow, sum, ratio;
    int i, j, k;

    x = dvector ((long) 1, (long) NumPts);
    w = dvector ((long) 1, (long) NumPts);

    Prep (q, C, dC, NumPts-2);  /* get coefficients */
    switch (mForm) {
        case 1:
#ifdef __MWERKS__
            printf ("%25s fit: I = %g\n", 
                    "constant background", fIntercept);
#else
            printf ("%25s fit: I = %lg\n", 
                    "constant background", fIntercept);
#endif
            break;
        case 2:
#ifdef __MWERKS__
            printf ("%25s fit: I = (%g) + q*(%g)\n", 
                    "linear", fIntercept, fSlope);
#else
            printf ("%25s fit: I = (%lg) + q*(%lg)\n", 
                    "linear", fIntercept, fSlope);
#endif
            break;
        case 3:
#ifdef __MWERKS__
            printf ("%25s fit: I = (%g) * q^(%g)\n", 
                    "Power law", fIntercept, fSlope);
#else
            printf ("%25s fit: I = (%lg) * q^(%lg)\n", 
                    "Power law", fIntercept, fSlope);
#endif
            break;
        case 4:
#ifdef __MWERKS__
            printf ("%25s fit: I = (%g) + (%g)/q^4\n", 
                    "Power law", fIntercept, fSlope);
#else
            printf ("%25s fit: I = (%lg) + (%lg)/q^4\n", 
                    "Power law", fIntercept, fSlope);
#endif
            break;
    }

    hLo = q[1];
    ratio = sLengt / (q[NumPts] - hLo);
    for (i = 1; i <= NumPts; i++) {
      x[i] = ratio * (q[i] - hLo);  /* values for "l" */
      w[i] = Plengt (x[i]);         /* probability at "l" */
    }

    w[1] *= x[2] - x[1];
    for (i = 2; i < NumPts; i++)
        w[i] *= x[i+1] - x[i-1];        /* step sizes */
    w[NumPts] = x[NumPts] - x[NumPts-1];

    for (i = 1; i <= NumPts; i++) {     /* evaluate each integral ... */
      Spinner (i);
      hNow = q[i];                      /* ... using trapezoid rule */
      sum = w[1] * FindIc (hNow, x[1], NumPts, q, C);
      for (k = 2; k < NumPts; k++)
        sum += w[k] * FindIc (hNow, x[k], NumPts, q, C);
      S[i] = sum + w[NumPts] * FindIc(hNow, x[NumPts], NumPts, q, C);
    }
    free_dvector (x, 1, NumPts);
    free_dvector (w, 1, NumPts);
}


double Plengt (double x)
    /*
     *  Here is the definition of the slit-length weighting function.
     *    It is defined for a rectangular slit of length 2*sLengt
     *    and probability 1/(2*sLengt).  It is zero elsewhere.
     *  It is not necessary to change the limit of the integration
     *    if the functional form here is changed.  You may, however,
     *    need to ask the user for more parameters.  Pass these
     *    around to the various routines through the use of the
     *    /PrepCm/ COMMON block.
     */
{
    return (fabs(x) > sLengt) ? 0.0 : 0.5 / sLengt;
}

#define GetIt(x,x1,y1,x2,y2)  (y1 + (y2-y1) * (x-x1) / (x2-x1))

double FindIc (double x, double y, int NumPts, double *q, double *C)
    /*
     *  Determine the "corrected" intensity at u = SQRT (x*x + y*y)
     *  Note that only positive values of "u" will be searched!
     */
{
    double u, value;
    int iTest, iLo, iHi;

    u = sqrt (x*x + y*y);               /* circularly symmetric */
    /*
     * dhunt(q, (unsigned long) NumPts, u, &iTest);
     * iTest++;
     */
    BSearch (u, q, NumPts, &iTest);     /* find index */
    iLo = iTest - 1;
    iHi = iLo + 1;
    if (iTest < 1) {
        printf ("\n\n Bad value of U or array Q in routine FindIc\n");
      exit (0);
    }
    if (iTest <= NumPts) {
      if (u == q[iLo])
        value = C[iLo];     /* exactly! */
      else                  /* linear interpolation */
        value = GetIt(u, q[iLo],C[iLo], q[iHi],C[iHi]);
    } else {                /* functional extrapolation */
      switch (mForm) {
        case 1: 
            value = fIntercept;
            break;
        case 2: 
            value = fIntercept + fSlope * u;
            break;
        case 3: 
            value = exp(fIntercept + fSlope * log(u));
            break;
        case 4: 
            value = fIntercept + fSlope / SQR(SQR(u));
            break;
      }
    }
    return value;
}


void DesmearData ()
{
    double ChiSqr, ChiSq0, weighting;
    int i, j, iteration;
    char reply[256], trimReply[256];

    printf ("Number of points read: %d\n", NumPts);
    printf ("Output file: %s\n", OutFil);
#ifdef __MWERKS__
    printf ("Slit length: %g (x-axis units)\n", sLengt);
    printf ("Final form approx. will begin at: %g (x-axis units)\n", sFinal);
#else
    printf ("Slit length: %lg (x-axis units)\n", sLengt);
    printf ("Final form approx. will begin at: %lg (x-axis units)\n", sFinal);
#endif
    printf ("Final form approximation: ");
    switch (mForm) {
        case 1: printf ("flat background, I(q) = B\n"); break;
        case 2: printf ("linear, I(h) = b + q * m\n"); break;
        case 3: printf ("power law, I = b * q^m\n"); break;
        case 4: printf ("Porod, I = Cp + B / q^4\n"); break;
    }
    printf ("Number of iterations: ");
    if (NumItr < InfItr)
        printf ("%d\n", NumItr);
    else
        printf ("infinite\n");
    printf ("Iterative weighting: ");
    switch (LakeForm) {
        case LakeUnit: printf ("unity\n"); break;
        case LakeFast: printf ("fast\n"); break;
        case LakeChi2: printf ("ChiSqr\n"); break;
    }

    /*
     *  To start Lake's method, assume that the 0-th approximation
     *    of the corrected intensity is the measured intensity.
     */
    for (i = 1; i <= NumPts; i++) {
        C[i] = E[i];
        dC[i] = dE[i];
    }
    printf ("\n Smearing to get first approximation...\n");
    fflush (stdout);
    Smear (S, NumPts, q, C, dC);

    ChiSqr = 0.0;                           /* find the ChiSqr */
    for (j = 1; j <= NumPts; j++)
      ChiSqr += SQR((S[j] - E[j])/dE[j]);
    ChiSq0 = ChiSqr;                        /* remember the first one */

    iteration = 0;
    do {
#ifdef __MWERKS__
#ifdef __INTEL__
        WinSIOUXclrscr();
#endif
#endif
        iteration++;
        if (NumItr < InfItr)
            printf ("\n Iteration #%d of %d iterations\n", 
                        iteration, NumItr);
        else
            printf ("\n Iteration #%d\n", iteration);
        printf ("Applying the iterative correction ...\n");
        fflush (stdout);
        switch (LakeForm) {
            case LakeUnit: weighting = 1.0; break;
            case LakeChi2: weighting = 2*sqrt(ChiSq0/ChiSqr); break;
        }
        for (j = 1; j <= NumPts; j++) {
            if (LakeForm == LakeFast)
                weighting = C[j] / S[j];
            C[j] += weighting * (E[j] - S[j]);
        }
        printf ("Examining scatter to calculate the errors ...\n");
        fflush (stdout);
        FixErr (NumPts, q, E, dE, C, dC);
        printf ("Smearing again ...\n");
        fflush (stdout);
        Smear (S, NumPts, q, C, dC);
        ChiSqr = 0.0;
        for (j = 1; j <= NumPts; j++) {
            resid[j] = (S[j] - E[j]) / dE[j];
            ChiSqr += SQR(resid[j]);
        }
        printf ("\n Residuals plot for iteration #%d\n", iteration);
        ResPlot (NumPts-1, resid);
#ifdef __MWERKS__
        printf ("ChiSqr = %g for #%d points\n", ChiSqr, NumPts);
#else
        printf ("ChiSqr = %lg for #%d points\n", ChiSqr, NumPts);
#endif
        fflush (stdout);
        if (NumItr >= InfItr) {
            if (AskYesOrNo ("Save this data?", 'n') == 'y') {
                printf ("What file name? (%s) ==> ", OutFil);
                fgets (reply, 256, stdin);
                /*
                 * need a string trim function here
                   On Windows under CodeWarrior, the user
                   presses <return> and the code receives
                   the EOL symbols.  Got to fix this.
                 */
                /* printf ("(%s,%d) <%s>\n", __FUNCTION__, __LINE__, reply); */
                strtrim(reply, trimReply, NULL);
                strcpy(reply, trimReply);
                /* printf ("(%s,%d) <%s>\n", __FUNCTION__, __LINE__, reply); */
                if (!reply[0] || (reply[0] == '\n'))
                    SavDat (OutFil, q, C, dC, NumPts);
                else
                    SavDat (reply,  q, C, dC, NumPts);
            }
            if (AskYesOrNo ("Continue iterating?", 'y') == 'n')
                iteration = InfItr;
            fflush (stdout);
        }
    } while (iteration < NumItr);
    if (NumItr < InfItr)
        SavDat (OutFil, q, C, dC, NumPts);
    printf ("Plot of log(desmeared intensity) vs. q ...\n");
    fflush (stdout);
    for (i = 1; i <= NumPts; i++)
        C[i] = log (fabs (C[i]));
    Plot (NumPts, q, C);
    fflush (stdout);
    printf ("Same, but now log-log ...\n");
    for (i = 1; i <= NumPts; i++)
        q[i] = log (fabs (q[i]));
    Plot (NumPts, q, C);
    fflush (stdout);
}

main ()
{
#ifdef THINK_C
    cecho2file("Lake.log", 0, stdout);
    _fcreator = 'QKPT';       /* output file owner is Kaleidagraph */
#endif
#ifdef __MWERKS__
    SIOUXSettings.autocloseonquit = 0;      /* don't automatically exit */
    SIOUXSettings.asktosaveonclose = 0;     /* don't ask to save window */
    SIOUXSettings.showstatusline = 0;       /* don't show the status line */
#endif
    printf ("\n");
    printf (" %s, by Pete R. Jemian\n", VERSION_INFORMATION);
    printf (" This executable was built %s, %s\n\n", __DATE__, __TIME__);
    printf (" SAS data desmearing using the iterative technique of JA Lake.\n");
    printf (" P.R. Jemian; Ph.D. Thesis (1990) Northwestern University, Evanston, IL, USA.\n\n");
    printf (" J.A. Lake; ACTA CRYST 23 (1967) 191-194.\n\n");
    /*
     * console_options.pause_atexit = 0;
     */

    do {
        GetInf (InFile, OutFil, &sLengt, &sFinal, 
            &NumItr, InfItr, &mForm, &LakeForm);
        if (!InFile[0]) break;          /* no file name given, so exit */
        if (NumItr == 0) NumItr = InfItr;
        fflush (stdout);
        printf ("Input file: %s\n", InFile);
        GetDat (InFile, &q, &E, &dE, &NumPts);
        C     = dvector ((long) 1, (long) NumPts);
        dC    = dvector ((long) 1, (long) NumPts);
        S     = dvector ((long) 1, (long) NumPts);
        resid = dvector ((long) 1, (long) NumPts);
        if (NumPts == 0) break;
        if (sFinal > q[NumPts-1]) break;

        fflush (stdout);
        DesmearData ();
        fflush (stdout);

        free_dvector (q,     1, NumPts);
        free_dvector (E,     1, NumPts);
        free_dvector (dE,    1, NumPts);
        free_dvector (C,     1, NumPts);
        free_dvector (dC,    1, NumPts);
        free_dvector (S,     1, NumPts);
        free_dvector (resid, 1, NumPts);
    } while (InFile[0]);
    printf ("\n");
    return (0);
}
