


/*.FE{P 10}{Tabellierung von Polynomsplines}
           {Tabellierung von Polynomsplines}*/

/* -------------------- DEKLARATIONEN splintab.h -------------------- */

int sptab       /* Tabellierung eines kubischen Polynomsplines .......*/
         (
          int  n,         /* Anzahl der Splinestuecke ................*/
          REAL xanf,      /* linker Intervallrand ....................*/
          REAL xend,      /* rechter Intervallrand ...................*/
          REAL deltx,     /* Schrittweite ............................*/
          int  anzahl,    /* maximale Tabellengroesse ................*/
          REAL x[],       /* Stuetzstellen ...........................*/
          REAL a[],       /* Splinekoeffizienten von (x-x[i])^0 ......*/
          REAL b[],       /* Splinekoeffizienten von (x-x[i])^1 ......*/
          REAL c[],       /* Splinekoeffizienten von (x-x[i])^2 ......*/
          REAL d[],       /* Splinekoeffizienten von (x-x[i])^3 ......*/
          REAL xtab[],    /* x-Koordinaten der Splinetabelle .........*/
          REAL ytab[],    /* y-Koordinaten der Splinetabelle .........*/
          int  *lentab    /* tatsaechliche Tabellengroesse ...........*/
         );               /* Fehlercode ..............................*/

int partab    /* Tabellierung eines parametr. kub. Polynomsplines ....*/
          (
           int  n,         /* Anzahl der Splinestuecke ...............*/
           REAL tanf,      /* linker Intervallrand ...................*/
           REAL tend,      /* rechter Intervallrand ..................*/
           REAL delt,      /* Schrittweite ...........................*/
           int  anzahl,    /* maximale Tabellengroesse ...............*/
           REAL t[],       /* Stuetzstellen ..........................*/
           REAL ax[],      /* x-Splinekoeffizienten von (t-t[i])^0 ...*/
           REAL bx[],      /* x-Splinekoeffizienten von (t-t[i])^1 ...*/
           REAL cx[],      /* x-Splinekoeffizienten von (t-t[i])^2 ...*/
           REAL dx[],      /* x-Splinekoeffizienten von (t-t[i])^3 ...*/
           REAL ay[],      /* y-Splinekoeffizienten von (t-t[i])^0 ...*/
           REAL by[],      /* y-Splinekoeffizienten von (t-t[i])^1 ...*/
           REAL cy[],      /* y-Splinekoeffizienten von (t-t[i])^2 ...*/
           REAL dy[],      /* y-Splinekoeffizienten von (t-t[i])^3 ...*/
           REAL xtab[],    /* x-Koordinaten der Splinetabelle ........*/
           REAL ytab[],    /* y-Koordinaten der Splinetabelle ........*/
           int  *lentab    /* tatsaechliche Tabellengroesse ..........*/
          );               /* Fehlercode .............................*/

int hmtab         /* Tabellierung eines Hermite-Polynomsplines .......*/
         (
          int  n,         /* Anzahl der Splinestuecke ................*/
          REAL xanf,      /* linker Intervallrand ....................*/
          REAL xend,      /* rechter Intervallrand ...................*/
          REAL deltx,     /* Schrittweite ............................*/
          int  anzahl,    /* maximale Tabellengroesse ................*/
          REAL x[],       /* Stuetzstellen ...........................*/
          REAL a[],       /* Splinekoeffizienten von (x-x[i])^0 ......*/
          REAL b[],       /* Splinekoeffizienten von (x-x[i])^1 ......*/
          REAL c[],       /* Splinekoeffizienten von (x-x[i])^2 ......*/
          REAL d[],       /* Splinekoeffizienten von (x-x[i])^3 ......*/
          REAL e[],       /* Splinekoeffizienten von (x-x[i])^4 ......*/
          REAL f[],       /* Splinekoeffizienten von (x-x[i])^5 ......*/
          REAL xtab[],    /* x-Koordinaten der Splinetabelle .........*/
          REAL ytab[],    /* y-Koordinaten der Splinetabelle .........*/
          int  *lentab    /* tatsaechliche Tabellengroesse ...........*/
         );               /* Fehlercode ..............................*/

int pmtab   /* Tabellierung eines parametr. Hermite-Polynomsplines ...*/
         (
          int  n,         /* Anzahl der Splinestuecke ................*/
          REAL tanf,      /* linker Intervallrand ....................*/
          REAL tend,      /* rechter Intervallrand ...................*/
          REAL delt,      /* Schrittweite ............................*/
          int  anzahl,    /* maximale Tabellengroesse ................*/
          REAL t[],       /* Stuetzstellen ...........................*/
          REAL ax[],      /* x-Splinekoeffizienten von (t-t[i])^0 ....*/
          REAL bx[],      /* x-Splinekoeffizienten von (t-t[i])^1 ....*/
          REAL cx[],      /* x-Splinekoeffizienten von (t-t[i])^2 ....*/
          REAL dx[],      /* x-Splinekoeffizienten von (t-t[i])^3 ....*/
          REAL ex[],      /* x-Splinekoeffizienten von (t-t[i])^4 ....*/
          REAL fx[],      /* x-Splinekoeffizienten von (t-t[i])^5 ....*/
          REAL ay[],      /* y-Splinekoeffizienten von (t-t[i])^0 ....*/
          REAL by[],      /* y-Splinekoeffizienten von (t-t[i])^1 ....*/
          REAL cy[],      /* y-Splinekoeffizienten von (t-t[i])^2 ....*/
          REAL dy[],      /* y-Splinekoeffizienten von (t-t[i])^3 ....*/
          REAL ey[],      /* y-Splinekoeffizienten von (t-t[i])^4 ....*/
          REAL fy[],      /* y-Splinekoeffizienten von (t-t[i])^5 ....*/
          REAL xtab[],    /* x-Koordinaten der Splinetabelle .........*/
          REAL ytab[],    /* y-Koordinaten der Splinetabelle .........*/
          int  *lentab    /* tatsaechliche Tabellengroesse ...........*/
         );               /* Fehlercode ..............................*/

int strtab   /* Tabellierung eines transf.-param. kub. Polynomspl. ...*/
          (
           int  n,         /* Anzahl der Splinestuecke ...............*/
           REAL panf,      /* Anfangswinkel ..........................*/
           REAL pend,      /* Endwinkel ..............................*/
           REAL phin[],    /* Winkelkoordinaten der Stuetzpunkte .....*/
           REAL a[],       /* Splinekoeff. von (phi-phin[i])^0 .......*/
           REAL b[],       /* Splinekoeff. von (phi-phin[i])^1 .......*/
           REAL c[],       /* Splinekoeff. von (phi-phin[i])^2 .......*/
           REAL d[],       /* Splinekoeff. von (phi-phin[i])^3 .......*/
           REAL phid,      /* Drehwinkel des Koordinatensystems ......*/
           REAL px,        /* x-Koordinate, ..........................*/
           REAL py,        /* y-Koordinate des Verschiebungsvektors ..*/
           REAL x[],       /* Stuetzstellen ..........................*/
           REAL y[],       /* Stuetzwerte ............................*/
           int  nl,        /* maximale Tabellengroesse ...............*/
           int  *nt,       /* tatsaechliche Tabellengroesse ..........*/
           REAL xtab[],    /* x-Koordinaten der Splinetabelle ........*/
           REAL ytab[]     /* y-Koordinaten der Splinetabelle ........*/
          );               /* Fehlercode .............................*/

/* ------------------------- ENDE splintab.h ------------------------ */
