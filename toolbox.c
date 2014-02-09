/** @file toolbox.c
 *  @brief part of my general mathematical "toolbox"
 *
 *  Pete R. Jemian, 15 May 1989.
 *  Some of these routines are taken
 *    (with reference) from book(s) but most, I have
 *    developed on my own.  They are modular in construction
 *    so that they may be improved, as needed.
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "recipes.h"

/**
 * request string input from the command line
 */
void AskString (char *question, char *answer)
{
    char reply[256];
    printf ("%s <%s> ==> ", question, answer);
    fgets (reply, 256, stdin);
    if (reply[0])
        sscanf (reply, "%s\n", answer);
}

/**
 * request a double input from the command line
 */
void AskDouble (char *question, double *answer)
{
    char reply[256];
#ifdef __MWERKS__
    printf ("%s <%g> ==> ", question, *answer);
#else
    printf ("%s <%lg> ==> ", question, *answer);
#endif
    fgets (reply, 256, stdin);
    if (reply[0])
        sscanf (reply, "%le", answer);
}

/**
 * request an integer input from the command line
 */
void AskInt (char *question, int *answer)
{
    char reply[256];
    printf ("%s <%d> ==> ", question, *answer);
    fgets (reply, 256, stdin);
    if (reply[0])
        sscanf (reply, "%d", answer);
}

/**
 * request yes or no from the command line
 */
char AskYesOrNo (char *question, char standard)
{
    char reply[256];
    do {
        printf ("%s (Y=yes, N=no) <%c> ==> ", question, standard);
        fgets (reply, 256, stdin);
        switch (reply[0]) {
            case (char) 0:
                reply[0] = standard;
                break;
            case 'y':
            case 'Y':
                reply[0] = 'y';
                break;
            case 'n':
            case 'N':
                reply[0] = 'n';
                break;
            default:
                reply[0] = 0;
        }
    } while (!reply);
    return reply[0];
}


/**
 * spins a stick to indicate program is still working
 * call many times
 */
void Spinner (int i)
{
    switch (i % 4) {
        case 1: printf ("%c", '-'); break;
        case 2: printf ("%c", '/'); break;
        case 3: printf ("%c", '|'); break;
        case 0: printf ("%c", '\\'); break;
    }
    printf ("%c", 0x08);   /* backspace */
    fflush(stdout);
}

/**
 * Does this line of text contain valid data?
 */
static int isDataLine (char *line)
{
  int isData = 0;
  if (strlen (line))       /* only consider non-blank lines */
    if (line[0] != '#')    /* but don't count comment lines */
      isData = 1;
  return (isData);
}


/**
 * read three-column data from a text file
 */
void GetDat (
     char *InFile,
     double **x,
     double **y, 
     double **dy, 
     int *n)
{
    FILE *path;
    int i;
    char aLine[256];
    double xx, yy, esd;

    *n = 0;
    path = fopen (InFile, "r");
    if (!path) return;
    i = 0;
    while (!feof (path)) {      /* count # lines in the file */
        fgets (aLine, 256, path);
        if (isDataLine (aLine)) i++;
    }
    *n = i - 1;
    *x = dvector ((long) 1, (long) *n);
    *y = dvector ((long) 1, (long) *n);
    *dy = dvector ((long) 1, (long) *n);
    freopen (InFile, "r", path);
    i = 0;
    while (!feof (path) && (i < *n)) {      /* go through each line in file */
        fgets (aLine, 256, path);
        if (isDataLine (aLine)) {
            i++;
            sscanf (aLine, "%le%le%le", &xx, &yy, &esd);
            (*x)[i] = xx;
            (*y)[i] = yy;
            (*dy)[i] = esd;
        }
    }
    fclose (path);
}

#ifdef __MWERKS__
#define SAVE_FORMAT "%g\t%g\t%g\n"
#else
#define SAVE_FORMAT "%lg\t%lg\t%lg\n"
#endif

/**
 * save three-columns of data to a text file
 */
void SavDat (char *OutFil, double *x, double *y, double *dy, int n)
{
    FILE *path;
    int i;

    /* printf ("(%s,%d) <%s>\n", __FUNCTION__, __LINE__, OutFil); */
    path = fopen (OutFil, "w");
    if (!path) return;
    printf ("Saving data in file: %s\n", OutFil);
    for (i = 1; i <= n; i++)
        fprintf (path, SAVE_FORMAT, x[i], y[i], dy[i]);
    fclose (path);
}


/**
 * swap a and b
 */
void Iswap (int *a, int *b)
{
    int c = *a;
    *a = *b;
    *b = c;
}

#define MaxRow (20)
#define MaxCol (75)
#define Blank = ' ';
#define Symbol = 'O';
#define hBordr = '-';
#define vBordr = '|';

/**
 * low-resolution ASCII text plot
 *
 * Make a scatter plot on the default display device (UNIT=*).
 * MaxRow and MaxCol correspond to the display dimensions.
 */
void Plot (int n, double *x, double *y)
{
    char screen[MaxRow+2][MaxCol+2+1];
    int i, r, c, nRow, nCol;
    double  xMin, xMax, yMin, yMax, ColDel, RowDel;

    if (n < 2) return;          /* can't draw anything important */
    /* wipe the "screen" clean & paint a border */
    for (r = 0; r < MaxRow+2; r++) {
        for (c = 0; c < MaxCol+2; c++)
            screen[r][c] = ' ';
        screen[r][MaxCol+2] = (char) 0;
    }
    for (r = 1; r <= MaxRow; r++) {
        screen[r][0]        = '|';
        screen[r][MaxCol+1] = '|';
    }
    for (c = 1; c <= MaxCol; c++) {
        screen[0][c]        = '-';
        screen[MaxRow+1][c] = '-';
    }
    /* get the data limits */
    xMin = x[1];
    xMax = x[1];
    yMin = y[1];
    yMax = y[1];
    for (i = 2; i <= n; i++) {
        xMin = DMIN (xMin, x[i]);
        xMax = DMAX (xMax, x[i]);
        yMin = DMIN (yMin, y[i]);
        yMax = DMAX (yMax, y[i]);
    }
    ColDel = (MaxCol - 1) / (xMax - xMin);
    RowDel = (MaxRow - 1) / (yMax - yMin);

    /*
     * plot the data points
     * data scaling functions are offset by +1 for plot frame
     */
    for (i = 1; i <= n; i++) {
        c = ColDel * (x[i] - xMin) + 1;
        r = RowDel * (y[i] - yMin) + 1;
        screen[r][c] = 'O';
    }

    /* convey the "screen" to the default output */
#ifdef __MWERKS__
    printf ("row:    min=%g   step=%g   max=%g\n", xMin, ColDel, xMax);
    printf ("column: min=%g   step=%g   max=%g\n", yMin, RowDel, yMax);
#else
    printf ("row:    min=%lg   step=%lg   max=%lg\n", xMin, ColDel, xMax);
    printf ("column: min=%lg   step=%lg   max=%lg\n", yMin, RowDel, yMax);
#endif
    for (r = MaxRow + 1; r >= 0; r--)
        printf ("%s\n", screen[r]);
}



/**
 *  Draw a plot of the standardized residuals on the screen.
 *  Mark the rows of + and - one standard deviation.
 */
void ResPlot (int n, double *x)
{
    char screen[MaxRow+2][MaxCol+2+1];
    int i, r, c, nRow, nCol, mPlus, mMinus, nPack;
    double  xMin, xMax, ColDel, RowDel;

    if (n < 2) return;          /* can't draw anything important */
    /* Find out how many points to pack per column and how many columns */
    nPack = ((double) n / MaxCol - 1.0/n) + 1;
    nCol  = ((double) n - 1./n)/nPack + 1;
    /* wipe the "screen" clean & paint a border */
    for (r = 0; r < MaxRow+2; r++) {
        for (c = 0; c < nCol+2; c++)
            screen[r][c] = ' ';
        screen[r][nCol+2] = (char) 0;
    }
    for (r = 1; r <= MaxRow; r++) {
        screen[r][0]      = '|';
        screen[r][nCol+1] = '|';
    }
    for (c = 1; c <= nCol; c++) {
        screen[0][c]        = '-';
        screen[MaxRow+1][c] = '-';
    }

    /* get the data limits */
    xMin = x[1];
    xMax = x[1];
    for (i = 2; i <= n; i++) {
        xMin = DMIN (xMin, x[i]);
        xMax = DMAX (xMax, x[i]);
    }
    RowDel = (MaxRow - 1) / (xMax - xMin);

    /* show the standard deviation bars */
    mMinus = RowDel*(-1.0 - xMin) + 1;
    mPlus  = RowDel*( 1.0 - xMin) + 1;
    for (c = 1; c < nCol+1; c++) {
        screen[mMinus][c] = '=';
        screen[mPlus][c]  = '=';
    }

    /*
     * draw the data (overdrawing the residuals bars if necessary)
     * data scaling functions (offset by +1 for the plot frame)
     */
    for (i = 1; i <= n; i++) {
        c = ((double) i - 1./n)/nPack + 1;
        r = RowDel * (x[i] - xMin) + 1;
        screen[r][c] = 'O';
    }

    /* convey the "screen" to the default output */
#ifdef __MWERKS__
    printf ("row:    min=%g   step=%g   max=%g\n", xMin, RowDel, xMax);
    printf ("column: min=%d    max=%d\n",  1,            n);
#else
    printf ("row:    min=%lg   step=%lg   max=%lg\n", xMin, RowDel, xMax);
    printf ("column: min=%d    max=%d\n",  1,            n);
#endif
    for (r = MaxRow + 1; r >= 0; r--)
        printf ("%s\n", screen[r]);
}

int iLo, iHi;   /* used to bracket the current index in FindIc */

/**
 *  Search the array "x" for (iLo) <= z < x(iHi)
 *  On exit, iLo and iHi will exactly bracket the datum
 *    and iTest will be the same as iLo.
 *  If z is below [above] the range, iTest = -1 [NumPts+1]. 
 */
void BSearch (double z, double *x, int NumPts, int *iTest)
{
    *iTest = -1;                /* assume that z < x[1] and test */
    if (z < x[1]) return;
    *iTest = NumPts + 1;        /* assume z > x[n] and test */
    if (z > x[NumPts]) return;
    if (iLo < 1 || iHi > NumPts || iLo >= iHi) {
      iLo = 1;
      iHi = NumPts;
    }
    while (z < x[iLo])
        iLo /= 2;
    while (z > x[iHi])      /* expand up? */
      iHi = (iHi + 1 + NumPts) / 2;
    *iTest = iHi;
    while (iHi - iLo > 1) {
        *iTest = (iLo + iHi) / 2;
        if (z >= x[*iTest])
            iLo = *iTest;
        else
            iHi = *iTest;
    }
}


/**
 * remove stuff from a string
 */
void strtrim (char *inStr, char *outStr, char* list)
{
  char *in = inStr;
  char *out = outStr;
  /* strcpy(outStr, inStr); */
  
  /*
   * For the time being, 
   * cut out any white space from the string
   */
  *outStr = 0;
  while (*in) {
    if (!isspace(*in)) {
      *out = *in;
      out++;
      *out = 0;   /* terminate the output string */
    }
    in++;
  }
}


