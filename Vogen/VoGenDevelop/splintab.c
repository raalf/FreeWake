/* ------------------------ MODUL splintab.c ------------------------ */

/***********************************************************************
*                                                                      *
* Funktionen zur Tabellierung von Polynomsplinefunktionen              *
* -------------------------------------------------------              *
*                                                                      *
* Programmiersprache: ANSI C                                           *
* Compiler:           Turbo C 2.0                                      *
* Rechner:            IBM PS/2 70 mit 80387                            *
* Bemerkung:          Umsetzung einer aequivalenten TP-Unit und eines  *
*                     aequivalenten QuickBASIC-Moduls                  *
* Autor:              Elmar Pohl (QuickBASIC)                          *
* Bearbeiter:         Juergen Dietel, Rechenzentrum der RWTH Aachen    *
* Datum:              MO 8. 7. 1991                                    *
*                                                                      *
***********************************************************************/

#include "basis.h"      /* wegen intervall, COS, SIN, REAL, ZERO, TWO */
#include "splintab.h"   /* wegen sptab, partab, hmtab, pmtab, strtab  */
#include <math.h>


/* ------------------------------------------------------------------ */

static void sptabh      /* Hilfsfunktion fuer sptab() ................*/
/*.IX{sptabh}*/
                  (
                   REAL xa,
                   REAL xe,
                   REAL dx,
                   REAL xi,
                   REAL a,
                   REAL b,
                   REAL c,
                   REAL d,
                   REAL xt[],
                   REAL yt[],
                   int  anzahl,
                   int  *lt
                  )

/***********************************************************************
* Diese Hilfsfunktion fuer die Funktion sptab() tabelliert eine        *
* kubische Splinefunktion im Inneren eines Stuetzstellenintervalls von *
* xa bis xe mit der Schrittweite dx.                                   *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* xa:      linker Rand des zu tabellierenden Intervalls                *
* xe:      rechter Rand des zu tabellierenden Intervalls               *
* dx:      Schrittweite                                                *
* xi:      linker Rand des Stuetzstellenintervalls                     *
* a,b,c,d: Koeffizienten des Splines im Stuetzstellenintervall         *
* xt:      [0..anzahl]-Feld mit alten x-Werten der Tabelle             *
* yt:      [0..anzahl]-Feld mit alten y-Werten der Tabelle             *
* anzahl:  obere Indexgrenze von xtab und ytab                         *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* xt: [0..anzahl]-Feld mit alten und neuen x-Werten der Tabelle        *
* yt: [0..anzahl]-Feld mit alten und neuen y-Werten der Tabelle        *
* lt: Index des letzten berechneten Tabellenpunktes                    *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* REAL
***********************************************************************/

{
  REAL x0,
       x1;

  xt += *lt, yt += *lt;
  for (x0 = xa; x0 < xe; x0 += dx)
  {
    if (*lt >= anzahl)
      return;
    (*lt)++;
    x1 = x0 - xi;
    *++xt = x0;
    *++yt = ((d * x1 + c) * x1 + b) * x1 + a;
  }
}


/* ------------------------------------------------------------------ */

int sptab       /* Tabellierung eines kubischen Polynomsplines .......*/
/*.IX{sptab}*/
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
         )                /* Fehlercode ..............................*/

/***********************************************************************
* eine Tabelle mit Funktionswerten einer kubischen Polynomspline-      *
* funktion erstellen.                                                  *
* Durch Rundungsfehler wegen einer unguenstig gewaehlten Schrittweite  *
* kann es vorkommen, dass an manchen Stuetzstellen und bei xend        *
* doppelt tabelliert wird.                                             *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* n:       Index der letzten Stuetzstelle in x                         *
* xanf:\   legen das x-Intervall fest, in dem die Funktion             *
* xend:/   tabelliert werden soll.                                     *
* deltx:   Tabellenschrittweite. Die Tabellenwerte werden              *
*          fuer  x = xanf(deltx)xend  berechnet.                       *
* x:       [0..n]-Feld mit den Stuetzstellen                           *
* a,b,c,d: [0..n-1]-Felder mit den Splinekoeffizienten.                *
*          Von a wird eventuell auch das Element a[n] benoetigt.       *
* anzahl:  obere Indexgrenze von xtab und ytab                         *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* xtab:   [0..anzahl]-Feld mit x-Werten der Tabelle. Die Dimension     *
*         anzahl sollte so gross gewaehlt werden, dass fuer alle       *
*         Tabellenwerte Platz ist.                                     *
* ytab:   [0..anzahl]-Feld mit y-Werten der Tabelle.                   *
*         ytab[i] ist der Funktionswert zu xtab[i], i=0(1)lentab       *
* lentab: Index des letzten berechneten Tabellenwertes                 *
*                                                                      *
* Funktionswert:                                                       *
* ==============                                                       *
* 0: kein Fehler                                                       *
* 1: xanf > xend                                                       *
* 2: deltx <= 0                                                        *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* sptabh, REAL, intervall, ZERO, TWO                                   *
***********************************************************************/

{
  int anf,    /* Index des zu xanf gehoerigen Splinepolynoms (0..n-1) */
      end,    /* Index des zu xend gehoerigen Splinepolynoms (0..n-1) */
      anf2,   /* Nummer der naechsten Stuetzstelle <= xanf (-1..n)    */
      end2,   /* Nummer der naechsten Stuetzstelle <= xend (-1..n)    */
      i;      /* Laufvariable                                         */

  if (xanf > xend)
    return 1;
  if (deltx <= ZERO)
    return 2;
  anf2 = anf = intervall(n, xanf, x);
  end2 = end = intervall(n, xend, x);
  if (xanf < x[0])        /* xanf links vom Interpolationsintervall?  */
    anf2--;
  if (xend > x[n])        /* xend rechts vom Interpolationsintervall? */
    end2++;
  *lentab = -1;
  if (anf2 < end2)
  {
    sptabh(xanf, x[anf2 + 1], deltx, x[anf], a[anf], b[anf],
           c[anf], d[anf], xtab, ytab, anzahl, lentab);
    for (i = anf2 + 1; i < end2; i++)
      sptabh(x[i], x[i + 1], deltx, x[i], a[i], b[i], c[i], d[i],
             xtab, ytab, anzahl, lentab);
    xanf = x[end2];                /* bei xanf die Tabelle fortsetzen */
    if (end2 == n)                         /* dafuer sorgen, dass der */
      if (*lentab < anzahl)                /* Funktionswert bei x[n]  */
        xtab[++(*lentab)] = x[n],          /* exakt eingetragen wird, */
        ytab[*lentab]     = a[n],          /* falls er im Tabellie-   */
        xanf += deltx;                     /* rungsintervall liegt    */
    sptabh(xanf, xend, deltx, x[end], a[end], b[end],
           c[end], d[end], xtab, ytab, anzahl, lentab);
  }
  else
    sptabh(xanf, xend, deltx, x[anf], a[anf], b[anf], c[anf],
           d[anf], xtab, ytab, anzahl, lentab);

  /* ----- den rechten Rand xend gesondert behandeln, da er oben ---- */
  /* ----- normalerweise unter den Tisch faellt                  ---- */
  sptabh(xend, xend + deltx / TWO, deltx, x[end], a[end], b[end],
         c[end], d[end], xtab, ytab, anzahl, lentab);

  return 0;
}



/* ------------------------------------------------------------------ */

static void partabh    /* Hilfsfunktion fuer partab() ................*/
/*.IX{partabh}*/
                   (
                    REAL ta,
                    REAL te,
                    REAL dt,
                    REAL ti,
                    REAL ax,
                    REAL bx,
                    REAL cx,
                    REAL dx,
                    REAL ay,
                    REAL by,
                    REAL cy,
                    REAL dy,
                    REAL xt[],
                    REAL yt[],
                    int  anzahl,
                    int  *lt
                   )

/***********************************************************************
* Diese Hilfsfunktion fuer die Funktion partab() tabelliert eine       *
* parametrische kubische Splinefunktion im Inneren eines Parameter-    *
* intervalls von ta bis te mit der Schrittweite dt.                    *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* ta:           linker Rand des zu tabellierenden Intervalls           *
* te:           rechter Rand des zu tabellierenden Intervalls          *
* dt:           Schrittweite                                           *
* ti:           linker Rand des Stuetzstellenintervalls                *
* ax,bx,cx,dx:\ Koeffizienten des Splines im Stuetzstellenintervall    *
* ay,by,cy,dy:/                                                        *
* xt:           [0..anzahl]-Feld mit alten x-Werten der Tabelle        *
* yt:           [0..anzahl]-Feld mit alten y-Werten der Tabelle        *
* anzahl:       obere Indexgrenze von xtab und ytab                    *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* xt: [0..anzahl]-Feld mit alten und neuen x-Werten der Tabelle        *
* yt: [0..anzahl]-Feld mit alten und neuen y-Werten der Tabelle        *
* lt: Index des letzten berechneten Tabellenpunktes                    *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* REAL
***********************************************************************/

{
  REAL t0,
       t1;

  xt += *lt, yt += *lt;
  for (t0 = ta; t0 < te; t0 += dt)
  {
    if (*lt >= anzahl)
      return;
    (*lt)++;
    t1 = t0 - ti;
    *++xt = ((dx * t1 + cx) * t1 + bx) * t1 + ax;
    *++yt = ((dy * t1 + cy) * t1 + by) * t1 + ay;
  }
}


/* ------------------------------------------------------------------ */

int partab    /* Tabellierung eines parametr. kub. Polynomsplines ....*/
/*.IX{partab}*/
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
          )                /* Fehlercode .............................*/

/***********************************************************************
* eine Tabelle mit Funktionswerten einer parametrischen kubischen      *
* Polynomsplinefunktion erstellen.                                     *
* Durch Rundungsfehler wegen einer unguenstig gewaehlten Schrittweite  *
* kann es vorkommen, dass an manchen Stuetzstellen und bei tend        *
* doppelt tabelliert wird.                                             *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* n:            Index des letzten Parameterwertes in t                 *
* tanf:\        legen das Parameterintervall fest, in dem die          *
* tend:/        Funktion tabelliert werden soll.                       *
* delt:         Tabellenschrittweite. Die Tabellenwerte werden         *
*               fuer t = tanf(delt)tend  berechnet.                    *
* t:            [0..n]-Feld mit den Parameterwerten des Splines        *
* ax,bx,cx,dx:\                                                        *
* ax,by,cy,dy:/ [0..n-1]-Felder mit den Splinekoeffizienten.           *
*               Von ax und ay werden eventuell auch die Elemente       *
*               ax[n] und ay[n] benoetigt.                             *
* anzahl:       obere Indexgrenze von XTAB und YTAB                    *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* xtab:   [0..anzahl]-Feld mit x-Werten der Tabelle. Die Dimension     *
*         anzahl muss so gross gewaehlt werden, dass fuer alle         *
*         Tabellenwerte Platz ist.                                     *
* ytab:   [0..anzahl]-Feld mit y-Werten der Tabelle                    *
* lentab: Index des letzten berechneten Tabellenwertes                 *
*                                                                      *
* Funktionswert:                                                       *
* ==============                                                       *
* 0: kein Fehler                                                       *
* 1: tanf > tend                                                       *
* 2: delt <= 0                                                         *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* partabh, REAL, intervall, ZERO, TWO                                  *
***********************************************************************/

{
  int anf,    /* Index des zu tanf gehoerigen Splinepolynoms (0..n-1) */
      end,    /* Index des zu tend gehoerigen Splinepolynoms (0..n-1) */
      anf2,   /* Nummer der naechsten Stuetzstelle <= tanf (-1..n)    */
      end2,   /* Nummer der naechsten Stuetzstelle <= tend (-1..n)    */
      i;      /* Laufvariable                                         */

  if (tanf > tend)
    return 1;
  if (delt <= ZERO)
    return 2;
  anf2 = anf = intervall(n, tanf, t);
  end2 = end = intervall(n, tend, t);
  if (tanf < t[0])        /* tanf links vom Interpolationsintervall?  */
    anf2--;
  if (tend > t[n])        /* tend rechts vom Interpolationsintervall? */
    end2++;
  *lentab = -1;
  if (anf2 < end2)
  {
    partabh(tanf, t[anf2 + 1], delt, t[anf], ax[anf], bx[anf],
            cx[anf], dx[anf], ay[anf], by[anf], cy[anf],
            dy[anf], xtab, ytab, anzahl, lentab);
    for (i = anf2 + 1; i < end2; i++)
       partabh(t[i], t[i + 1], delt, t[i], ax[i], bx[i], cx[i],
               dx[i], ay[i], by[i], cy[i], dy[i], xtab, ytab,
               anzahl, lentab);
    tanf = t[end2];                /* bei tanf die Tabelle fortsetzen */
    if (end2 == n)                         /* dafuer sorgen, dass der */
      if (*lentab < anzahl)                /* Funktionswert bei t[n]  */
        xtab[++(*lentab)] = ax[n],         /* exakt eingetragen wird, */
        ytab[*lentab]     = ay[n],         /* falls er im Tabellie-   */
        tanf += delt;                      /* rungsintervall liegt    */
    partabh(tanf, tend, delt, t[end], ax[end], bx[end],
            cx[end], dx[end], ay[end], by[end], cy[end],
            dy[end], xtab, ytab, anzahl, lentab);
  }
  else
    partabh(tanf, tend, delt, t[anf], ax[anf], bx[anf],
            cx[anf], dx[anf], ay[anf], by[anf], cy[anf],
            dy[anf], xtab, ytab, anzahl, lentab);

  /* ----- den rechten Rand xend gesondert behandeln, da er oben ---- */
  /* ----- normalerweise unter den Tisch faellt                  ---- */
  partabh(tend, tend + delt / TWO, delt, t[end], ax[end], bx[end],
          cx[end], dx[end], ay[end], by[end], cy[end], dy[end],
          xtab, ytab, anzahl, lentab);

  return 0;
}



/* ------------------------------------------------------------------ */

static void hmtabh      /* Hilfsfunktion fuer hmtab() ................*/
/*.IX{hmtabh}*/
                  (
                   REAL xa,
                   REAL xe,
                   REAL dx,
                   REAL xi,
                   REAL a,
                   REAL b,
                   REAL c,
                   REAL d,
                   REAL e,
                   REAL f,
                   REAL xt[],
                   REAL yt[],
                   int  anzahl,
                   int  *lt
                  )

/***********************************************************************
* Diese Hilfsfunktion fuer die Funktion hmtab() tabelliert eine        *
* Hermite-Splinefunktion im Inneren eines Stuetzstellenintervalls      *
* von xa bis xe mit der Schrittweite dx.                               *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* xa:          linker Rand des zu tabellierenden Intervalls            *
* xe:          rechter Rand des zu tabellierenden Intervalls           *
* dx:          Schrittweite                                            *
* xi:          linker Rand des Stuetzstellenintervalls                 *
* a,b,c,d,e,f: Koeffizienten des Splines im Stuetzstellenintervall     *
* xt:          [0..anzahl]-Feld mit alten x-Werten der Tabelle         *
* yt:          [0..anzahl]-Feld mit alten y-Werten der Tabelle         *
* anzahl:      obere Indexgrenze von xtab und ytab                     *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* xt: [0..anzahl]-Feld mit alten und neuen x-Werten der Tabelle        *
* yt: [0..anzahl]-Feld mit alten und neuen y-Werten der Tabelle        *
* lt: Index des letzten berechneten Tabellenpunktes                    *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* REAL
***********************************************************************/

{
  REAL x0,
       x1;

  xt += *lt, yt += *lt;
  for (x0 = xa; x0 < xe; x0 += dx)
  {
    if (*lt >= anzahl)
      return;
    (*lt)++;
    x1 = x0 - xi;
    *++xt = x0;
    *++yt = ((((f * x1 + e) * x1 + d) * x1 + c) * x1 + b) * x1 + a;
  }
}


/* ------------------------------------------------------------------ */

int hmtab         /* Tabellierung eines Hermite-Polynomsplines .......*/
/*.IX{hmtab}*/
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
         )                /* Fehlercode ..............................*/

/***********************************************************************
* eine Tabelle mit Funktionswerten einer Hermite-Polynomsplinefunktion *
* erstellen. Durch Rundungsfehler wegen einer unguenstig gewaehlten    *
* Schrittweite kann es vorkommen, dass an manchen Stuetzstellen und    *
* bei xend doppelt tabelliert wird.                                    *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* n:           Index der letzten Stuetzstelle in x                     *
* xanf:\       legen das x-Intervall fest, in dem die Funktion         *
* xend:/       tabelliert werden soll.                                 *
* deltx:       Tabellenschrittweite. Die Tabellenwerte werden          *
*              fuer x = xanf(deltx)xend  berechnet.                    *
* x:           [0..n]-Feld mit den Stuetzstellen                       *
* a,b,c,d,e,f: [0..n-1]-Felder mit den Splinekoeffizienten.            *
*              Von a wird eventuell auch das Element a[n] benoetigt.   *
* anzahl:      obere Indexgrenze von xtab und ytab                     *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* xtab:   [0..anzahl]-Feld mit x-Werten der Tabelle. Die Dimension     *
*         anzahl sollte so gross gewaehlt werden, dass fuer alle       *
*         Tabellenwerte Platz ist.                                     *
* ytab:   [0..anzahl]-Feld mit y-Werten der Tabelle.                   *
*         ytab[i] ist der Funktionswert zu xtab[i], i=0(1)lentab       *
* lentab: Index des letzten berechneten Tabellenwertes                 *
*                                                                      *
* Funktionswert:                                                       *
* ==============                                                       *
* 0: kein Fehler                                                       *
* 1: xanf > xend                                                       *
* 2: deltx <= 0                                                        *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* hmtabh, REAL, intervall, ZERO, TWO                                   *
***********************************************************************/

{
  int anf,    /* Index des zu xanf gehoerigen Splinepolynoms (0..n-1) */
      end,    /* Index des zu xend gehoerigen Splinepolynoms (0..n-1) */
      anf2,   /* Nummer der naechsten Stuetzstelle <= xanf (-1..n)    */
      end2,   /* Nummer der naechsten Stuetzstelle <= xend (-1..n)    */
      i;      /* Laufvariable                                         */

  if (xanf > xend)
    return 1;
  if (deltx <= ZERO)
    return 2;
  anf2 = anf = intervall(n, xanf, x);
  end2 = end = intervall(n, xend, x);
  if (xanf < x[0])        /* xanf links vom Interpolationsintervall?  */
    anf2--;
  if (xend > x[n])        /* xend rechts vom Interpolationsintervall? */
    end2++;
  *lentab = -1;
  if (anf2 < end2)
  {
    hmtabh(xanf, x[anf2 + 1], deltx, x[anf], a[anf], b[anf], c[anf],
           d[anf], e[anf], f[anf], xtab, ytab, anzahl, lentab);
    for (i = anf2 + 1; i < end2; i++)
      hmtabh(x[i], x[i + 1], deltx, x[i], a[i], b[i], c[i], d[i], e[i],
             f[i], xtab, ytab, anzahl, lentab);
    xanf = x[end2];                /* bei xanf die Tabelle fortsetzen */
    if (end2 == n)                         /* dafuer sorgen, dass der */
      if (*lentab < anzahl)                /* Funktionswert bei x[n]  */
        xtab[++(*lentab)] = x[n],          /* exakt eingetragen wird, */
        ytab[*lentab]     = a[n],          /* falls er im Tabellie-   */
        xanf += deltx;                     /* rungsintervall liegt    */
    hmtabh(xanf, xend, deltx, x[end], a[end], b[end], c[end], d[end],
           e[end], f[end], xtab, ytab, anzahl, lentab);
  }
  else
    hmtabh(xanf, xend, deltx, x[anf], a[anf], b[anf], c[anf], d[anf],
           e[anf], f[anf], xtab, ytab, anzahl, lentab);

  /* ----- den rechten Rand xend gesondert behandeln, da er oben ---- */
  /* ----- normalerweise unter den Tisch faellt                  ---- */
  hmtabh(xend, xend + deltx / TWO, deltx, x[end], a[end], b[end],
         c[end], d[end], e[end], f[end], xtab, ytab, anzahl, lentab);

  return 0;
}



/* ------------------------------------------------------------------ */

static void pmtabh      /* Hilfsfunktion fuer pmtab() ................*/
/*.IX{pmtabh}*/
                  (
                   REAL ta,
                   REAL te,
                   REAL dt,
                   REAL ti,
                   REAL ax,
                   REAL bx,
                   REAL cx,
                   REAL dx,
                   REAL ex,
                   REAL fx,
                   REAL ay,
                   REAL by,
                   REAL cy,
                   REAL dy,
                   REAL ey,
                   REAL fy,
                   REAL xt[],
                   REAL yt[],
                   int  anzahl,
                   int  *lt
                  )

/***********************************************************************
* Diese Hilfsfunktion fuer die Funktion partab() tabelliert eine       *
* parametrische Hermite-Splinefunktion im Inneren eines Parameter-     *
* intervalls von ta bis te mit der Schrittweite dt.                    *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* ta:                 linker Rand des zu tabellierenden Intervalls     *
* te:                 rechter Rand des zu tabellierenden Intervalls    *
* dt:                 Schrittweite                                     *
* ti:                 linker Rand des Stuetzstellenintervalls          *
* ax,bx,cx,dx,ex,fx:\ Koeffizienten des Splines im                     *
* ay,by,cy,dy,ey,fy:/ Stuetzstellenintervall                           *
* xt:                 [0..anzahl]-Feld mit alten x-Werten der Tabelle  *
* yt:                 [0..anzahl]-Feld mit alten y-Werten der Tabelle  *
* anzahl:             obere Indexgrenze von xtab und ytab              *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* xt: [0..anzahl]-Feld mit alten und neuen x-Werten der Tabelle        *
* yt: [0..anzahl]-Feld mit alten und neuen y-Werten der Tabelle        *
* lt: Index des letzten berechneten Tabellenpunktes                    *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* REAL
***********************************************************************/
{
  REAL t0,
       t;

  xt += *lt, yt += *lt;
  for (t0 = ta; t0 < te; t0 += dt)
  {
    if (*lt >= anzahl)
      return;
    (*lt)++;
    t = t0 - ti;
    *++xt = ((((fx * t + ex) * t + dx) * t + cx) * t + bx) * t + ax;
    *++yt = ((((fy * t + ey) * t + dy) * t + cy) * t + by) * t + ay;
  }
}


/* ------------------------------------------------------------------ */

int pmtab   /* Tabellierung eines parametr. Hermite-Polynomsplines ...*/
/*.IX{pmtab}*/
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
         )                /* Fehlercode ..............................*/

/***********************************************************************
* eine Tabelle mit Funktionswerten einer parametrischen Hermite-       *
* Polynomsplinefunktion erstellen.                                     *
* Durch Rundungsfehler wegen einer unguenstig gewaehlten Schrittweite  *
* kann es vorkommen, dass an manchen Stuetzstellen und bei tend        *
* doppelt tabelliert wird.                                             *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* n:                  Index des letzten Parameterwertes in t           *
* tanf:\              legen das Parameterintervall fest, in dem die    *
* tend:/              Funktion tabelliert werden soll.                 *
* delt:               Tabellenschrittweite. Die Tabellenwerte werden   *
*                     fuer t = tanf(delt)tend  berechnet.              *
* t:                  [0..n]-Feld mit den Parameterwerten des Splines  *
* ax,bx,cx,dx,ex,fx:\                                                  *
* ax,by,cy,dy,ey,fy:/ [0..n-1]-Felder mit den Splinekoeffizienten.     *
*                     Von ax und ay werden eventuell auch die Elemente *
*                     ax[n] und ay[n] benoetigt.                       *
* anzahl:             obere Indexgrenze von XTAB und YTAB              *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* xtab:   [0..anzahl]-Feld mit x-Werten der Tabelle. Die Dimension     *
*         anzahl muss so gross gewaehlt werden, dass fuer alle         *
*         Tabellenwerte Platz ist.                                     *
* ytab:   [0..anzahl]-Feld mit y-Werten der Tabelle                    *
* lentab: Index des letzten berechneten Tabellenwertes                 *
*                                                                      *
* Funktionswert:                                                       *
* ==============                                                       *
* 0: kein Fehler                                                       *
* 1: tanf > tend                                                       *
* 2: delt <= 0                                                         *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* pmtabh, REAL, intervall, ZERO, TWO                                   *
***********************************************************************/

{
  int anf,    /* Index des zu tanf gehoerigen Splinepolynoms (0..n-1) */
      end,    /* Index des zu tend gehoerigen Splinepolynoms (0..n-1) */
      anf2,   /* Nummer der naechsten Stuetzstelle <= tanf (-1..n)    */
      end2,   /* Nummer der naechsten Stuetzstelle <= tend (-1..n)    */
      i;      /* Laufvariable                                         */

  if (tanf > tend)
    return 1;
  if (delt <= ZERO)
    return 2;
  anf2 = anf = intervall(n, tanf, t);
  end2 = end = intervall(n, tend, t);
  if (tanf < t[0])        /* tanf links vom Interpolationsintervall?  */
    anf2--;
  if (tend > t[n])        /* tend rechts vom Interpolationsintervall? */
    end2++;
  *lentab = -1;
  if (anf2 < end2)
  {
    pmtabh(tanf, t[anf2 + 1], delt, t[anf], ax[anf], bx[anf],
           cx[anf], dx[anf], ex[anf], fx[anf], ay[anf], by[anf],
           cy[anf], dy[anf], ey[anf], fy[anf], xtab, ytab, anzahl,
           lentab);
    for (i = anf2 + 1; i < end2; i++)
      pmtabh(t[i], t[i + 1], delt, t[i], ax[i], bx[i], cx[i], dx[i],
             ex[i], fx[i], ay[i], by[i], cy[i], dy[i], ey[i], fy[i],
             xtab, ytab, anzahl, lentab);
    tanf = t[end2];                /* bei tanf die Tabelle fortsetzen */
    if (end2 == n)                         /* dafuer sorgen, dass der */
      if (*lentab < anzahl)                /* Funktionswert bei t[n]  */
        xtab[++(*lentab)] = ax[n],         /* exakt eingetragen wird, */
        ytab[*lentab]     = ay[n],         /* falls er im Tabellie-   */
        tanf += delt;                      /* rungsintervall liegt    */
    pmtabh(tanf, tend, delt, t[end], ax[end], bx[end], cx[end],
           dx[end], ex[end], fx[end], ay[end], by[end], cy[end],
           dy[end], ey[end], fy[end], xtab, ytab, anzahl, lentab);
  }
  else
    pmtabh(tanf, tend, delt, t[anf], ax[anf], bx[anf], cx[anf],
           dx[anf], ex[anf], fx[anf], ay[anf], by[anf], cy[anf],
           dy[anf], ey[anf], fy[anf], xtab, ytab, anzahl, lentab);

  /* ----- den rechten Rand xend gesondert behandeln, da er oben ---- */
  /* ----- normalerweise unter den Tisch faellt                  ---- */
  pmtabh(tend, tend + delt / TWO, delt, t[end], ax[end], bx[end],
         cx[end], dx[end], ex[end], fx[end], ay[end], by[end],
         cy[end], dy[end], ey[end], fy[end], xtab, ytab, anzahl,
         lentab);

  return 0;
}



/* ------------------------------------------------------------------ */

static void strtabh    /* Hilfsfunktion fuer strtab() ................*/
/*.IX{strtabh}*/
                   (
                    REAL pa,
                    REAL pe,
                    REAL dp,
                    REAL pi,
                    REAL a,
                    REAL b,
                    REAL c,
                    REAL d,
                    REAL phid,
                    REAL px,
                    REAL py,
                    REAL xt[],
                    REAL yt[],
                    int  anzahl,
                    int  *lt
                   )

/***********************************************************************
* Diese Hilfsfunktion fuer die Funktion strtab() tabelliert eine       *
* transformiert-parametrische kubische Splinefunktion im Inneren eines *
* Stuetzstellenintervalls von pa bis pe mit der Schrittweite dp.       *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* pa:      linker Rand des zu tabellierenden Intervalls                *
* pe:      rechter Rand des zu tabellierenden Intervalls               *
* dp:      Schrittweite                                                *
* pi:      linker Rand des Stuetzstellenintervalls                     *
* a,b,c,d: Koeffizienten des Splines im Stuetzstellenintervall         *
* phid:    Drehwinkel des Koordinatensystems                           *
* py,py:   Verschiebungsvektor                                         *
* xt:      [0..anzahl]-Feld mit alten x-Werten der Tabelle             *
* yt:      [0..anzahl]-Feld mit alten y-Werten der Tabelle             *
* anzahl:  obere Indexgrenze von xtab und ytab                         *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* xt: [0..anzahl]-Feld mit alten und neuen x-Werten der Tabelle        *
* yt: [0..anzahl]-Feld mit alten und neuen y-Werten der Tabelle        *
* lt: Index des letzten berechneten Tabellenpunktes                    *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* REAL, SIN, COS
***********************************************************************/

{
  REAL p0,
       p1,
       s,
       rho;

  xt += *lt, yt += *lt;
  for (p0 = pa; p0 < pe; p0 += dp)
  {
    if (*lt >= anzahl)
      return;
    (*lt)++;
    p1    = p0 - pi;
    s     = ((d * p1 + c) * p1 + b) * p1 + a;
    rho   = p0 + phid;
    *++xt = s * COS(rho) + px;
    *++yt = s * SIN(rho) + py;
  }
}


/* ------------------------------------------------------------------ */

int strtab   /* Tabellierung eines transf.-param. kub. Polynomspl. ...*/
/*.IX{strtab}*/
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
          )                /* Fehlercode .............................*/

/***********************************************************************
* eine transformiert-parametrische kubische Polynomsplinefunktion in   *
* der Darstellung                                                      *
*   s(phi) = a[i] + b[i](phi-phin[i]) + c[i](phi-phin[i])^2 +          *
*                                     + d[i](phi-phin[i])^3            *
* fuer phi aus  [phin[i], phin[i+1]], i=0(1)n-1,  tabellieren.         *
* Es wird ein Tabelle mit den Wertepaaren                              *
*        xtab = xtab(phi) = s(phi) * cos(phi + phid) + px,             *
*        ytab = ytab(phi) = s(phi) * sin(phi + phid) + py,             *
* phi aus  [panf, pend]  erstellt, wobei gilt:                         *
*   - Ist panf < phin[0], wird fuer alle Werte xtab(phi) mit           *
*     phi < phin[0] das Randpolynom p[0] ausgewertet.                  *
*   - Ist pend > phin[n], wird fuer alle Werte xtab(phi) mit           *
*     phi > phin[n] das Randpolynom p[n-1] ausgewertet.                *
*   - Die Intervallgrenzen panf und pend und alle dazwischenliegenden  *
*     Stuetzstellen phin[i] werden auf jeden Fall tabelliert.          *
*   - In jedem Teilintervall [phin[i],phin[i+1]] wird die Tabelle mit  *
*     aequidistanter Schrittweite h erzeugt, wobei h jeweils von der   *
*     Intervallaenge und der gewaehlten Tabellenlaenge nl abhaengt.    *
*   - Der Eingabeparameter nl gibt die ungefaehre Tabellenlaenge vor;  *
*     die tatsaechliche Tabellenlaenge ist nt + 1. (nt ist der letzte  *
*     benutzte Tabellenindex.) Fuer nt gilt:                           *
*                     0 < nt < nl + n + 1.                             *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* n:       Index des letzten Stuetzpunktes                             *
* panf:\   legen das Parameterintervall fest, in dem die Funktion      *
* pend:/   tabelliert werden soll.                                     *
* phin:    [0..n]-Feld mit Parameterwerten der Splineknoten            *
* a,b,c,d: [0..n-1]-Felder mit den Splinekoeffizienten                 *
* phid:    Drehwinkel des Koordinatensystems                           *
* py,py:   Verschiebungsvektor                                         *
* nl:      obere Indexgrenze von xtab und ytab                         *
* x,y:     [0..n]-Vektoren mit den urspruenglichen Stuetzpunkten (wird *
*          benoetigt, um an den Stuetzstellen phin[i] exakte Tabellen- *
*          werte zu bekommen)                                          *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* nt:   Index des letzten berechneten Tabellenpunktes                  *
* xtab: [0..nl]-Feld mit x-Werten der Tabelle                          *
* ytab: [0..nl]-Feld mit y-Werten der Tabelle                          *
*                                                                      *
* Funktionswert:                                                       *
* ==============                                                       *
* 0: kein Fehler                                                       *
* 1: panf >= pend                                                      *
* 2: n < 1                                                             *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* strtabh, REAL, intervall, TWO                                        *
***********************************************************************/

{
  int  anf,   /* Index des zu panf gehoerigen Splinepolynoms (0..n-1) */
       end,   /* Index des zu pend gehoerigen Splinepolynoms (0..n-1) */
       anf2,  /* Nummer der naechsten Stuetzstelle <= panf (-1..n)    */
       end2,  /* Nummer der naechsten Stuetzstelle <= pend (-1..n)    */
       i;     /* Laufvariable                                         */
  REAL h;     /* Schrittweite                                         */

  if (pend <= panf)
    return 1;
  if (n < 1)
    return 2;

  anf2 = anf = intervall(n, panf, phin);
  end2 = end = intervall(n, pend, phin);
  if (panf < phin[0])     /* panf links vom Interpolationsintervall?  */
    anf2--;
  if (pend > phin[n])     /* pend rechts vom Interpolationsintervall? */
    end2++;
                                      /* die Schrittweite gerade so   */
  h = (pend - panf) / (nl - n - 1);   /* klein waehlen, dass auch die */
  *nt = -1;                           /* Stuetzpunkte in der Tabelle  */
                                      /* Platz finden                 */
  if (anf2 < end2)        /* Stuetzstellen im Tabellierungsintervall? */
  {
    if (panf == phin[anf])
      if (*nt < nl)                                /* den Stuetzpunkt */
        xtab[++(*nt)] = x[anf],                    /* (x[anf],y[anf]  */
        ytab[*nt]     = y[anf],                    /* noetigenfalls   */
        panf += h;                                 /* exakt eintragen */
    strtabh(panf, phin[anf2 + 1], h,         /* die Tabellenwerte von */
            phin[anf], a[anf],               /* panf bis phin[anf2+1] */
            b[anf], c[anf], d[anf],          /* berechnen             */
            phid, px, py, xtab, ytab,
            nl, nt);

    for (i = anf2 + 1; i < end2; i++)        /* die Tabellenwerte von */
    {                                        /* phin[anf2+1] bis      */
      panf = phin[i];                        /* phin[end2] berechnen  */
      if (*nt < nl)
        xtab[++(*nt)] = x[i],                    /* den Stuetzpunkt   */
        ytab[*nt]     = y[i],                    /* (x[i],y[i]) exakt */
        panf += h;                               /* eintragen         */
      strtabh(panf, phin[i + 1], h,
              phin[i], a[i], b[i],
              c[i], d[i], phid, px,
              py, xtab, ytab, nl, nt);
    }

    panf = phin[end2];
    if (*nt < nl)
      xtab[++(*nt)] = x[end2],                   /* den Stuetzpunkt   */
      ytab[*nt]     = y[end2],                   /* (x[end2],y[end2]) */
      panf += h;                                 /* exakt eintragen   */
    strtabh(panf, pend, h, phin[end],        /* die Tabellenwerte von */
            a[end], b[end], c[end],          /* phin[end2] bis pend   */
            d[end], phid, px, py,            /* berechnen             */
            xtab, ytab, nl, nt);
  }
  else              /* keine Stuetzstellen im Tabellierungsintervall? */
    strtabh(panf, pend, h, phin[anf], a[anf], b[anf], c[anf],
            d[anf], phid, px, py, xtab, ytab, nl, nt);

  /* ----- den rechten Rand xend gesondert behandeln, da er oben ---- */
  /* ----- normalerweise unter den Tisch faellt                  ---- */
  strtabh(pend, pend + h / TWO, h, phin[end], a[end], b[end],
          c[end], d[end], phid, px, py, xtab, ytab, nl, nt);

  return 0;
}

/* ------------------------- ENDE splintab.c ------------------------ */
