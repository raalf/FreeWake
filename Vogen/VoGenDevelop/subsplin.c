/* ------------------------ MODUL subsplin.c ------------------------ */

/***********************************************************************
*                                                                      *
* Funktionen zur Berechnung von Akima- und Renner-Subsplines           *
* ----------------------------------------------------------           *
*                                                                      *
* Programmiersprache: ANSI C                                           *
* Compiler:           Borland C++ 2.0                                  *
* Rechner:            IBM PS/2 70 mit 80387                            *
* Bemerkung:          Umsetzung einer aequivalenten TP-Unit, aber      *
*                     teilweise Verwendung neuer Algorithmen           *
* Autor:              Elmar Pohl (QuickBASIC)                          *
* Bearbeiter:         Juergen Dietel, Rechenzentrum der RWTH Aachen    *
* Datum:              DO 10. 9. 1992                                   *
*                                                                      *
***********************************************************************/

#include "basis.h"     /* wegen sqr, FABS, SQRT, MACH_EPS, ONE, ZERO, */
                       /*       TWO, THREE, REAL, SIX                 */
#include "vmblock.h"   /* wegen vmalloc, vmcomplete, vmfree, vminit,  */
                       /*       VEKTOR, VVEKTOR                       */
#include "subsplin.h"  /* wegen akima, renner                         */

#include <math.h>

/*.FE{P 13.1}{Akima--Subsplines}{Akima--Subsplines}*/

/* ------------------------------------------------------------------ */

static int eckrund
/*.IX{eckrund}*/
                  (
				   int  *n,        /* Nummer des letzten Stuetzpunkts */
                   int  nmax,      /* obere Indexgrenze von x und y   */
                   REAL x[],       /* Stuetzstellen                   */
                   REAL y[],       /* Stuetzwerte                     */
                   REAL beta       /* Mass der Rundung                */
                  )                /* Fehlercode                      */

/***********************************************************************
* Diese Hilfsfunktion fuer akima() rundet Ecken bei Akima-Subsplines   *
* ab, indem sie jeden Eckpunkt durch zwei neue, leicht verschobene     *
* Punkte ersetzt.                                                      *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* n     Nummer des letzten Stuetzpunktes. Die Numerierung beginnt      *
*       bei 0. Das Minimum ist 4.                                      *
* x     [0..n]-Vektor mit den Stuetzstellen.                           *
*       Die x[i] muessen streng monoton steigen.                       *
* y     [0..n]-Vektor mit den zu interpolierenden Werten.              *
*       Die Felder x und x muessen gross genug dimensioniert sein,     *
*       denn durch die Eckrundung koennen maximal [(n+1)/2] Punkte     *
*       hinzukommen.                                                   *
* beta  bestimmt die Staerke der Rundung. Jeder Eckpunkt wird um eine  *
*       von beta abhaengige Distanz in Richtung seines linken und      *
*       seines rechten Nachbarpunkts verschoben und dann durch die     *
*       beiden verschobenen Punkte ersetzt.                            *
*       Falls beta nicht im Intervall ]0,1[ liegt, wird keine          *
*       Eckenrundung vorgenommen.                                      *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* n     neuer, eventuell groesserer Wert fuer n                        *
* x,y   Felder wie bei der Eingabe, jedoch eventuell mit neu einge-    *
*       fuegten Punkten.                                               *
*                                                                      *
* Funktionswert:                                                       *
* ==============                                                       *
* 0: kein Fehler (oder beta nicht im Intervall ]0,1[)                  *
* 2: Die x[i] sind nicht streng monoton steigend.                      *
* 5: kein Platz mehr fuer Eckenrundungen                               *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* REAL, MACH_EPS, min, FABS, ZERO                                      *
***********************************************************************/

{
  REAL qimin2,          /* Sekantensteigung in P[i-2]                 */
       qimin1,          /* Sekantensteigung in P[i-1]                 */
       qi,              /* Sekantensteigung in P[i]                   */
       qiplus1,         /* Sekantensteigung in P[i+1]                 */
       L, R, B,         /* Hilfsvariablen zur Berechnung der Gewichte */
       lambda,          /* Gewicht fuer die Verschiebung nach hinten  */
       my;              /* Gewicht fuer die Verschiebung nach vorne   */
  int  i, j;            /* Laufvariablen                              */


  if (beta <= ZERO || beta >= ONE)  /* keine Eckenrundung gewuenscht? */
    return 0;

  if (x[1] == x[0] || x[2] == x[1] ||     /* keine strenge Monotonie? */
      x[3] == x[2])
    return 2;


  qimin2 = (y[1] - y[0]) / (x[1] - x[0]);
  qimin1 = (y[2] - y[1]) / (x[2] - x[1]);
  qi     = (y[3] - y[2]) / (x[3] - x[2]);


  for (i = 2; i <= *n - 2; i++)
  {
    if (x[i + 2] == x[i + 1])             /* keine strenge Monotonie? */
      return 2;

    qiplus1 = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);

    if (FABS(qiplus1 - qi) +                            /* Ist hier   */
        FABS(qimin1 - qimin2) <  MACH_EPS &&            /* eine Ecke? */
        FABS(qi - qimin1)     >= MACH_EPS)
    {

      if (*n == nmax - 1)        /* kein Platz mehr fuer neue Punkte? */
        return 5;

      for (j = *n; j >= i; j--)     /* die Punkte P[i] ... P[n] einen */
        x[j + 1] = x[j],            /* Platz nach oben ruecken lassen */
        y[j + 1] = y[j];

      /* ---------- P[i] und P[i+1] ein Stueck verschieben ---------- */

      L = 2 * (x[i]     - x[i - 1]);
      R = 2 * (x[i + 2] - x[i + 1]);
      B = beta * min(L, R);
      lambda   = B / L;
      my       = B / R;
      x[i]     = x[i]     - lambda * (x[i]     - x[i - 1]);
      y[i]     = y[i]     - lambda * (y[i]     - y[i - 1]);
      x[i + 1] = x[i + 1] + my     * (x[i + 2] - x[i + 1]);
      y[i + 1] = y[i + 1] + my     * (y[i + 2] - y[i + 1]);


      (*n)++;                    /* Wir haben jetzt einen Punkt mehr. */

      qimin2 = (y[i]     - y[i - 1]) / (x[i]     - x[i - 1]);
      qimin1 = (y[i + 1] - y[i])     / (x[i + 1] - x[i]);
      qi     = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
    }

    else                                          /* keine Ecke hier? */
      qimin2 = qimin1,
      qimin1 = qi,
      qi     = qiplus1;
  }


  return 0;
}



/* ------------------------------------------------------------------ */

int akima
/*.IX{akima}*/
         (
		  int  *n,                 /* Nummer des letzten Stuetzpunkts */
          int  nmax,               /* obere Indexgrenze fuer x, y,... */
          REAL x[],                /* Stuetzstellen                   */
          REAL y[],                /* Stuetzwerte                     */
          int  perio,              /* periodische Interpolation?      */
          REAL beta,               /* Mass fuer die Eckenrundung      */
          REAL b[],                /* Splinekoeffizienten             */
          REAL c[],
          REAL d[]
         )                         /* Fehlercode                      */

/***********************************************************************
* die Koeffizienten eines interpolierenden Akima-Subsplines berechnen. *
*                                                                      *
* Die Akima-Subsplinefunktion kann aehnlich wie eine kubische Spline-  *
* funktion mit der Funktion spwert() aus dem Modul spliwert ausgewer-  *
* tet und mit der Funktion sptab() aus splintab tabelliert werden.     *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* n      Nummer der letzten Stuetzstelle. Die Numerierung beginnt      *
*        bei 0. Das Minimum ist 4.                                     *
* nmax   obere Indexgrenze der Vektoren x, y, b, c, d                  *
*        (zur Pruefung, ob noch Eckenrundungen moeglich sind)          *
* x      [0..n]-Vektor mit den Stuetzstellen.                          *
*        Die x[i] muessen streng monoton steigen.                      *
* y      [0..n]-Vektor mit den zu interpolierenden Werten.             *
*        Fuer periodische Interpolation (perio gesetzt) muss ausserdem *
*        gelten:  y[0] = y[n].                                         *
* perio  Flagge fuer periodische oder nichtperiodische Interpolation.  *
*        perio nicht gesetzt (FALSE): nichtperiodische Interpolation   *
*        perio gesetzt       (TRUE):  periodische Interpolation        *
* beta   Masszahl fuer die Rundung von Ecken.                          *
*        Falls Ecken auftreten, werden sie abgerundet, indem jeder     *
*        Eckpunkt durch zwei neue Punkte ersetzt wird. Diese neuen     *
*        Punkte entstehen dadurch, dass der Eckpunkt um eine von beta  *
*        abhaengige Distanz laengs seines links bzw. rechts benachbar- *
*        ten Sehnenvektors verschoben wird.                            *
*        Die Felder x, y, a..d muessen gross genug dimensioniert wer-  *
*        den, um alle neu einzufuegenden Punkte aufnehmen zu koennen,  *
*        denn es koennen maximal [(n+1)/2] Punkte hinzukommen.         *
*        Falls beta nicht im Intervall ]0,1[ liegt, wird keine         *
*        Eckenrundung vorgenommen.                                     *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* b \   [0..n-1]-Felder mit den Splinekoeffizienten nach dem Ansatz    *
* c  >     s(x)  =  a[i] + b[i] * (x - x[i]) + c[i] * (x - x[i]) ^ 2   *
* d /                    + d[i] * (x - x[i]) ^ 3                       *
*       a entspricht y und hat daher noch ein zusaetzliches            *
*       Element a[n].                                                  *
*                                                                      *
* und bei beta aus]0,1[:                                               *
* n     neuer Wert fuer n mit eingefuegten Zwischenpunkten             *
* x,y   wie bei der Eingabe, aber mit eventuell eingefuegten           *
*       Zwischenpunkten                                                *
*                                                                      *
* Funktionswert:                                                       *
* ==============                                                       *
* 0: kein Fehler                                                       *
* 1: n < 4                                                             *
* 2: Die x[i] sind nicht streng monoton steigend.                      *
* 3: bei gesetztem perio:  y[0] != y[n]                                *
* 4: Speichermangel                                                    *
* 5: kein Platz mehr fuer Eckenrundungen                               *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* eckrund, REAL, MACH_EPS, sqr, FABS, vminit, vmalloc, vmcomplete,     *
* vmfree, VEKTOR, ZERO, TWO, THREE                                     *
***********************************************************************/

{
  void *vmblock;   /* Liste der dynamisch vereinbarten Vektoren und   */
                   /* Matrizen                                        */
  REAL *tl,        /* [0..n]-Vektor mit den linksseitigen Steigungen  */
       *tr,        /* [0..n]-Vektor mit den rechtsseitigen Steigungen */
       *m,         /* [-2..n+1]-Vektor mit den Sehnensteigungen       */
       hi,         /* Laenge des betrachteten Stuetzstellenintervalls */
       nenner,     /* Nenner bei der Berechnung von alpha             */
       alpha;      /* Faktor bei der Berechnung der rechts- und       */
                   /* linksseitigen Steigungen                        */
  int  i,          /* Laufvariable                                    */
       fehler;     /* Fehlercode von eckrund()                        */


  if (*n < 4)
    return 1;


  fehler = eckrund(n, nmax, x, y, beta);    /* Ecken eventuell runden */
  if (fehler != 0)
    return fehler;


  /* ----------- Speicher fuer die Hilfsfelder anfordern: ----------- */
  /* ----------- 1 [0..n+3]-Vektor und 2 [0..n]-Vektoren  ----------- */

  vmblock = vminit();                 /* Speicherblock initialisieren */
  tl = (REAL *)vmalloc(vmblock, VEKTOR, *n + 1, 0);
  tr = (REAL *)vmalloc(vmblock, VEKTOR, *n + 1, 0);
  m  = (REAL *)vmalloc(vmblock, VEKTOR, *n + 4, 0);
  if (! vmcomplete(vmblock))   /* Ging eine der Speicheranforderungen */
    return 4;                  /* fuer den Block schief?              */

  m += 2;                /* Nun kann m auch negativ indiziert werden. */


  for (i = 0; i < *n; i++)   /* die Steigungen m[0]..m[n-1] berechnen */
  {                          /* und dabei die Stuetzstellen auf       */
    hi = x[i + 1] - x[i];    /* strenge Monotonie pruefen             */
    if (hi <= ZERO)          /* nicht erfuellt?                       */
    {
      vmfree(vmblock);       /* Speicherplatz freigeben               */
      return 2;              /* Fehler melden                         */
    }
    m[i] = (y[i + 1] - y[i]) / hi;
  }


  if (perio)                  /* periodische Akima-Interpolation?     */
  {
    if (y[*n] != y[0])        /* ungeeignete Stuetzpunkte?            */
    {
      vmfree(vmblock);        /* Speicherplatz freigeben              */
      return 3;               /* Fehler melden                        */
    }
    m[-2]     = m[*n - 2];    /* die Randsteigungen durch periodische */
    m[-1]     = m[*n - 1];    /* Fortsetzung gewinnen                 */
    m[*n]     = m[0];
    m[*n + 1] = m[1];
  }

  else                                /* normale Akima-Interpolation? */
                                      /* die Randsteigungen berechnen */
    m[-2]     = THREE * m[0]      - TWO * m[1],
    m[-1]     = TWO   * m[0]      -       m[1],
    m[*n]     = TWO   * m[*n - 1] -       m[*n - 2],
    m[*n + 1] = THREE * m[*n - 1] - TWO * m[*n - 2];


  for (i = 0; i <= *n; i++)    /* die links- und rechtsseitigen Stei- */
  {                            /* gungen tl[i] und tr[i], i=0(1)n,    */
                               /* berechnen                           */
    nenner = FABS(m[i + 1] - m[i]) + FABS(m[i - 1] - m[i - 2]);
    if (nenner >= MACH_EPS)
      alpha = FABS(m[i - 1] - m[i - 2]) / nenner,
      tr[i] = tl[i] = m[i - 1] + alpha * (m[i] - m[i - 1]);
    else
      tl[i] = m[i - 1],
      tr[i] = m[i];
  }


  for (i = 0; i < *n; i++)    /* die Subsplinekoeffizienten berechnen */
    hi   = x[i + 1] - x[i],
    b[i] = tr[i],
    c[i] = (THREE * m[i] - TWO * tr[i] - tl[i + 1]) / hi,
    d[i] = (tr[i] + tl[i + 1] - TWO * m[i]) / sqr(hi);


  vmfree(vmblock);                         /* Speicherplatz freigeben */
  return 0;                                /* Erfolg melden           */
}



/*.FE{P 13.2}{Renner--Subsplines}{Renner--Subsplines}*/

/* ------------------------------------------------------------------ */

/**********************************************************************
* Typ fuer Punkte in der Ebene in kartesischen Koordinaten und fuer   *
* Richtungsvektoren                                                   *
**********************************************************************/

typedef struct { REAL x, y; } punkttyp;
/*.IX{punkttyp}*/



/* ------------------------------------------------------------------ */

static REAL abskreu2
/*.IX{abskreu2}*/
                    (
					 REAL x[],      /* Stuetzstellen \ Stuetzpunkte P */
                     REAL y[],      /* Stuetzwerte   /                */
                     int  i         /* Index des betrachteten Punktes */
                    )               /* Flaecheninhalt bzw. Fehlercode */

/***********************************************************************
* Diese Hilfsfunktion fuer eckrundp() berechnet die Flaeche des von    *
* den beiden normierten Sehnenvektoren                                 *
*                 (x[i+1]-x[i],   y[i+1]-y[i])0    und                 *
*                 (x[i+2]-x[i+1], y[i+2]-y[i+1])0                      *
* aufgespannten Parallelogramms und gibt sie als Funktionswert         *
* zurueck. Falls einer der Sehnenvektoren der Nullvektor ist, wird     *
* dies durch den Funktionswert -1 angezeigt.                           *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* REAL, ZERO, SQRT, sqr, FABS, ONE, punkttyp                           *
***********************************************************************/

{
  punkttyp si,               /* Sehnenvektor s[i]   = P[i+1] - P[i]   */
           sip1;             /* Sehnenvektor s[i+1] = P[i+2] - P[i+1] */


  si.x   = x[i + 1] - x[i];
  si.y   = y[i + 1] - y[i];
  sip1.x = x[i + 2] - x[i + 1];
  sip1.y = y[i + 2] - y[i + 1];

  if ((si.x   == ZERO && si.y   == ZERO) ||  /* zwei zusammenfallende */
      (sip1.x == ZERO && sip1.y == ZERO))    /* Stuetzpunkte?         */
    return -ONE;


  return FABS(si.x * sip1.y - sip1.x * si.y) /
         (SQRT(sqr(si.x)   + sqr(si.y)) *
          SQRT(sqr(sip1.x) + sqr(sip1.y)));
}



/* ------------------------------------------------------------------ */

static int eckrundp
/*.IX{eckrundp}*/
                   (
				    int  *n,       /* Nummer des letzten Stuetzpunkts */
                    int  nmax,     /* obere Indexgrenze von x und y   */
                    REAL x[],      /* Abszissen,                      */
                    REAL y[],      /* Ordinaten der Stuetzpunkte P    */
                    REAL beta      /* Mass der Rundung                */
                   )               /* Fehlercode                      */

/***********************************************************************
* Diese Hilfsfunktion fuer renner() rundet Ecken bei Renner-Subsplines *
* ab, indem sie jeden Eckpunkt durch zwei neue, leicht verschobene     *
* Punkte ersetzt.                                                      *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* n     Nummer des letzten Stuetzpunktes. Die Numerierung beginnt      *
*       bei 0. Das Minimum ist 4.                                      *
* nmax  obere Indexgrenze der Vektoren x, y                            *
*       (zur Pruefung, ob noch Eckenrundungen moeglich sind)           *
* x     [0..n]-Vektor mit den x-Koordinaten der Stuetzpunkte           *
* y     [0..n]-Vektor mit den y-Koordinaten der Stuetzpunkte           *
*       Aufeinanderfolgende Punkte duerfen nicht zusammenfallen.       *
*       Die Felder x und y muessen gross genug dimensioniert sein,     *
*       denn durch die Eckrundung koennen maximal [(n+1)/2] Punkte     *
*       hinzukommen.                                                   *
* beta  bestimmt die Staerke der Rundung. Jeder Eckpunkt wird um eine  *
*       von beta abhaengige Distanz in Richtung seines linken und      *
*       seines rechten Nachbarpunkts verschoben und dann durch die     *
*       beiden verschobenen Punkte ersetzt.                            *
*       Falls beta nicht im Intervall ]0,1[ liegt, wird keine          *
*       Eckenrundung vorgenommen.                                      *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* n     neuer, eventuell groesserer Wert fuer n                        *
* x,y   Felder wie bei der Eingabe, jedoch eventuell mit neu einge-    *
*       fuegten Punkten.                                               *
*                                                                      *
* Funktionswert:                                                       *
* ==============                                                       *
* 0: kein Fehler (oder beta nicht im Intervall ]0,1[)                  *
* 2: Zwei aufeinanderfolgende Punkte fallen zusammen.                  *
* 4: kein Platz mehr fuer Eckenrundungen                               *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* abskreu2, REAL, min, MACH_EPS, SQRT, sqr, ZERO, TWO                  *
***********************************************************************/

{
  REAL simin2kr,        /* Flaecheninhalt F(s0[i-2],s0[i-1])          */
       simin1kr,        /* Flaecheninhalt F(s0[i-1],s0[i])            */
       sikr,            /* Flaecheninhalt F(s0[i],  s0[i+1])          */
       L, R, B,         /* Hilfsvariablen zur Berechnung der Gewichte */
       lambda,          /* Gewicht fuer die Verschiebung nach hinten  */
       my;              /* Gewicht fuer die Verschiebung nach vorne   */
  int  i, j;            /* Laufvariablen                              */


  if (beta <= ZERO || beta >= ONE)  /* keine Eckenrundung gewuenscht? */
    return 0;


  if ((simin2kr = abskreu2(x, y, 0)) < ZERO ||
      (simin1kr = abskreu2(x, y, 1)) < ZERO)
    return 2;

  for (i = 2; i <= *n - 2; i++)         /* nach Ecken Ausschau halten */
  {
    if ((sikr = abskreu2(x, y, i)) < ZERO)
      return 2;

    if (simin2kr + sikr <  MACH_EPS &&         /* Liegt beim Punkt    */
        simin1kr        >= MACH_EPS)           /* P[i] eine Ecke vor? */
    {
      if (*n == nmax - 1)        /* kein Platz mehr fuer neue Punkte? */
        return 4;

      for (j = *n; j >= i; j--)     /* die Punkte P[i] ... P[n] einen */
        x[j + 1] = x[j],            /* Platz nach oben ruecken lassen */
        y[j + 1] = y[j];

      /* ---------- P[i] und P[i+1] ein Stueck verschieben ---------- */

      L        = TWO * SQRT(sqr(x[i]     - x[i - 1]) +
                            sqr(y[i]     - y[i - 1]));
      R        = TWO * SQRT(sqr(x[i + 2] - x[i + 1]) +
                            sqr(y[i + 2] - y[i + 1]));
      B        = beta * min(L, R);
      lambda   = B / L;
      my       = B / R;
      x[i]     = x[i]     - lambda * (x[i]     - x[i - 1]);
      y[i]     = y[i]     - lambda * (y[i]     - y[i - 1]);
      x[i + 1] = x[i + 1] + my     * (x[i + 2] - x[i + 1]);
      y[i + 1] = y[i + 1] + my     * (y[i + 2] - y[i + 1]);

      (*n)++;          /* die frohe Botschaft verkuenden, dass ein    */
                       /* neuer Punkt das Licht der Welt erblickt hat */
      simin2kr = abskreu2(x, y, i - 3);
      simin1kr = abskreu2(x, y, i - 2);
    }

    else                               /* keine Ecke beim Punkt P[i]? */
      simin2kr = simin1kr,
      simin1kr = sikr;

  }


  return 0;
}



/* ------------------------------------------------------------------ */

static REAL abskreuz
/*.IX{abskreuz}*/
                    (
					 punkttyp s[],    /* Sehnenvektoren               */
                     int      i       /* Startindex der beiden zu be- */
                                      /* trachtenden Sehnenvektoren   */
                    )                 /* Flaecheninhalt               */

/***********************************************************************
* Diese Hilfsfunktion fuer renner() berechnet die Flaeche des von den  *
* beiden Sehnenvektoren s[i] und s[i+1] aufgespannten Parallelogramms. *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* REAL, punkttyp, FABS                                                 *
***********************************************************************/

{
  return FABS(s[i].x * s[i + 1].y - s[i].y * s[i + 1].x);
}



/* ------------------------------------------------------------------ */

int renner
/*.IX{renner}*/
          (
		   int  *n,            /* Nummer des letzten Stuetzpunkts     */
           int  nmax,          /* obere Indexgrenze fuer x, y,...     */
           REAL x[],           /* x-Komponenten und                   */
           REAL y[],           /* y-Komponenten der Stuetzpunkte      */
           REAL beta,          /* Mass fuer die Eckenrundung          */
           REAL T[],           /* Parameterwerte zu den Stuetzpunkten */
           REAL bx[],          /* x-Komponenten der                   */
           REAL cx[],          /* Splinekoeffizienten                 */
           REAL dx[],
           REAL by[],          /* y-Komponenten der                   */
           REAL cy[],          /* Splinekoeffizienten                 */
           REAL dy[]
          )                    /* Fehlercode                          */

/***********************************************************************
* die Koeffizienten eines interpolierenden Renner-Subsplines berechnen *
*                                                                      *
* Die Interpolationsfunktion kann aehnlich wie eine parametrische      *
* kubische Splinefunktion mit der Funktion pspwert() aus dem Modul     *
* spliwert ausgewertet und mit der Funktion partab() aus splintab      *
* tabelliert werden. Allerdings setzen diese Funktionen einen monoton  *
* steigenden Kurvenparameter voraus, waehrend der Kurvenparameter bei  *
* Renner-Subsplines an jedem Knoten bei 0 beginnt. Man summiert also   *
* vorher die T[i] zu monoton steigenden Parameterwerten Tsum[i] auf,   *
* z. B. durch folgendes Programmstueck, wobei Tsum vom selben Typ      *
* wie T ist:                                                           *
*                                                                      *
*             for (Tsum[0] = 0, i = 1; i <= n; i++)                    *
*               Tsum[i] = Tsum[i-1] + T[i-1];                          *
*                                                                      *
* Dann setzt man in pspwert() und partab() Tsum anstelle von T ein.    *
* Wegen der besonderen Eigenschaften der Renner-Parameter ist uebri-   *
* gens Tsum[i] annaehernd gleich der Bogenlaenge vom Anfangsknoten bis *
* zum Knoten mit der Nummer i.                                         *
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* n     Nummer des letzten Knotens. Die Numerierung beginnt bei 0.     *
*       Das Minimum ist 4.                                             *
* nmax  obere Indexgrenze der Vektoren x, y, T, bx, cx, dx, by, cy, dy *
*       (zur Pruefung, ob noch Eckenrundungen moeglich sind)           *
* x \   [0..n]-Vektoren mit den zu interpolierenden Punkten            *
* y /           P[i] = (x[i],y[i]), i=0(1)n                            *
*       Je zwei direkt aufeinanderfolgende Punkte muessen voneinander  *
*       verschieden sein.                                              *
* beta  Masszahl fuer die Rundung von Ecken.                           *
*       Falls Ecken auftreten, werden sie abgerundet, indem jeder Eck- *
*       punkt durch zwei neue Punkte ersetzt wird. Diese neuen Punkte  *
*       entstehen dadurch, dass der Eckpunkt um eine von beta abhaen-  *
*       gige Distanz laengs seines links bzw. rechts benachbarten      *
*       Sehnenvektors verschoben wird.                                 *
*       Die Felder x, y, t, ax..dx und ay..dy muessen gross genug di-  *
*       mensioniert werden, um alle neu einzufuegenden Punkte aufneh-  *
*       men zu koennen, denn es koennen maximal [(n+1)/2] Punkte hin-  *
*       zukommen.                                                      *
*       Falls beta nicht im Intervall ]0,1[ liegt, wird keine          *
*       Eckenrundung vorgenommen.                                      *
*                                                                      *
* Ausgabeparameter:                                                    *
* =================                                                    *
* T      [0..n]-Vektor mit den Laengen der Parameterintervalle zu den  *
*        einzelnen Kurvensegmenten von P[i] bis P[i+1], i=0(1)n        *
* bx \   [0..n-1]-Felder mit den Koeffizienten der x-Komponente        *
* cx  >  der Subspline-Funktion                                        *
* dx /                                                                 *
* by \   [0..n-1]-Felder mit den Koeffizienten der y-Komponente        *
* cy  >  der Subspline-Funktion                                        *
* dy /                                                                 *
*        Die Subspline-Funktion ist fuer i=0(1)n-1 gegeben durch       *
*                                                                      *
*                (ax[i])       (bx[i])         (cx[i])         (dx[i]) *
*        S(t) =  (     ) + t * (     ) + t^2 * (     ) + t^3 * (     ) *
*                (ay[i])       (by[i])         (cy[i])         (dy[i]) *
*                                                                      *
*        fuer 0 <= t <= T[i].                                          *
*        ax entspricht x, ay entspricht y, daher haben beide Vektoren  *
*        noch ein zusaetzliches Element ax[n] bzw ay[n].               *
*                                                                      *
* und bei beta aus ]0,1[:                                              *
* n      neuer Wert fuer n mit eingefuegten Zwischenpunkten            *
* x,y    wie bei der Eingabe, aber mit eventuell eingefuegten          *
*        Zwischenpunkten                                               *
*                                                                      *
* Funktionswert:                                                       *
* ==============                                                       *
* 0: kein Fehler                                                       *
* 1: n < 4                                                             *
* 2: Zwei aufeinanderfolgende Punkte sind gleich.                      *
* 3: Speichermangel                                                    *
* 4: kein Platz mehr fuer Eckenrundungen                               *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* punkttyp, eckrundp, abskreuz, REAL, MACH_EPS, sqr, FABS, vminit,     *
* vmalloc, vmcomplete, vmfree, VVEKTOR, SQRT, TWO, THREE, SIX          *
***********************************************************************/

{
  void     *vmblock;      /* Liste der dynamisch vereinbarten         */
                          /* Vektoren und Matrizen                    */
  int      m,             /* n nach der Eckenrundung                  */
           fehler,        /* Fehlercode von eckrundp()                */
           i;             /* Laufvariable                             */
  punkttyp *tL,           /* [0..n]-Vektor mit den linksseitigen      */
                          /* Tangenteneinheitsvektoren, zugleich Zei- */
                          /* ger auf den gesamten dynamisch belegten  */
                          /* Speicher                                 */
           *tR,           /* [0..n]-Vektor mit den rechtsseitigen     */
                          /* Tangenteneinheitsvektoren                */
           *s,            /* [-2..n+1]-Vektor mit den Sehnenvektoren  */
           *s0,           /* [-2..n+1]-Vektor mit den Sehneneinheits- */
                          /* vektoren                                 */
           hilf;          /* noch nicht normierter Tangentenvektor    */
                          /* tL[i] = tR[i], spaeter Summe von         */
                          /* tR[i] und tL[i+1]                        */
  REAL     norm,          /* Euklidische Norm von Vektoren            */
           alpha,         /* Gewichtungsfaktor bei der Berechnung der */
                          /* links- und rechtsseitigen Tangentenein-  */
                          /* heitsvektoren aus den benachbarten       */
                          /* Sehnenvektoren                           */
           nenner,        /* Nenner bei der Berechnung von alpha      */
           A, B, C,       /* Hilfsvariablen zur Berechnung der        */
                          /* Laengen der Parameterintervalle T[i]     */
           TiInv,         /* Kehrwert von T[i]                        */
           TiInvSq;       /* Quadrat des Kehrwerts von T[i]           */


  if (*n < 4)                               /* zuwenige Stuetzpunkte? */
    return 1;


  fehler = eckrundp(n, nmax, x, y, beta);   /* Ecken eventuell runden */
  if (fehler != 0)
    return fehler;

  m = *n;


  /* ---------- Speicher fuer die Hilfsfelder anfordern:   ---------- */
  /* ---------- 2 [0..m]-Vektoren und 2 [-2..m+1]-Vektoren ---------- */

  vmblock = vminit();                 /* Speicherblock initialisieren */
  tL = (punkttyp *)vmalloc(vmblock, VVEKTOR, m + 1, sizeof(punkttyp));
  tR = (punkttyp *)vmalloc(vmblock, VVEKTOR, m + 1, sizeof(punkttyp));
  s  = (punkttyp *)vmalloc(vmblock, VVEKTOR, m + 4, sizeof(punkttyp));
  s0 = (punkttyp *)vmalloc(vmblock, VVEKTOR, m + 4, sizeof(punkttyp));
  if (! vmcomplete(vmblock))   /* Ging eine der Speicheranforderungen */
  {                            /* fuer den Block schief?              */
    vmfree(vmblock);
    return 3;
  }
                /* dafuer sorgen, dass s und s0 auch mit negativen    */
  s  += 2;      /* Indizes verwendet werden koennen, wie es der       */
  s0 += 2;      /* Algorithmus erfordert, so dass man sich eine nur   */
                /* Verwirrung stiftende Indexverschiebung sparen kann */


  for (i = 0; i < m; i++)               /* die inneren Sehnenvektoren */
  {                                     /* s[i], i=0(1)m-1, berechnen */
    s[i].x = x[i + 1] - x[i];
    s[i].y = y[i + 1] - y[i];
    if (FABS(s[i].x) < MACH_EPS &&        /* Fallen zwei aufeinander- */
        FABS(s[i].y) < MACH_EPS)          /* folgende Stuetzpunkte    */
    {                                     /* zusammen?  =>  Fehler!   */
      vmfree(vmblock);
      return 2;
    }
  }


  /* ------- weitere Sehnenvektoren s[-2], s[-1], s[m], s[m+1] ------ */
  /* ------- bereitstellen                                     ------ */

  if (FABS(x[0] - x[m]) < MACH_EPS &&      /* geschlossene Kurve?     */
      FABS(y[0] - y[m]) < MACH_EPS)
    s[-2]    = s[m - 2],                   /* die Sehnenvektoren      */
    s[-1]    = s[m - 1],                   /* am Rand durch zyklische */
    s[m]     = s[0],                       /* Wiederholung bilden     */
    s[m + 1] = s[1];

  else                                               /* offene Kurve? */
    s[-2].x    = THREE * s[0].x     - TWO * s[1].x,      /* die Seh-  */
    s[-2].y    = THREE * s[0].y     - TWO * s[1].y,      /* nenvek-   */
    s[-1].x    = TWO   * s[0].x     -       s[1].x,      /* toren am  */
    s[-1].y    = TWO   * s[0].y     -       s[1].y,      /* Rand aus  */
    s[m].x     = TWO   * s[m - 1].x -       s[m - 2].x,  /* den be-   */
    s[m].y     = TWO   * s[m - 1].y -       s[m - 2].y,  /* nachbar-  */
    s[m + 1].x = THREE * s[m - 1].x - TWO * s[m - 2].x,  /* ten Vek-  */
    s[m + 1].y = THREE * s[m - 1].y - TWO * s[m - 2].y;  /* toren     */
                                                         /* herleiten */


  for (i = -2; i <= m + 1; i++)                /* die Sehneneinheits- */
  {                                            /* vektoren berechnen  */
    norm = SQRT(sqr(s[i].x) + sqr(s[i].y));
    if (norm < MACH_EPS)                      /* verschwindende Norm? */
      s0[i].x = s0[i].y = ZERO;               /* Nullvektor eintragen */
    else
      s0[i].x = s[i].x / norm,                /* normierten Vektor    */
      s0[i].y = s[i].y / norm;                /* notieren             */
  }


  for (i = 0; i <= m; i++)      /* die links- und rechtsseitigen Tan- */
  {                             /* genteneinheitsvektoren berechnen   */
    nenner = abskreuz(s0, i - 2) + abskreuz(s0, i);
    if (nenner < MACH_EPS)                 /* verschwindender Nenner? */
      tL[i] = s0[i - 1],                   /* Bei (x[i],y[i]) liegt   */
      tR[i] = s0[i];                       /* eine Ecke vor.          */
    else
      alpha  =  abskreuz(s0, i - 2) / nenner,
      hilf.x =  s[i - 1].x + alpha * (s[i].x - s[i - 1].x),
      hilf.y =  s[i - 1].y + alpha * (s[i].y - s[i - 1].y),
      norm   =  SQRT(sqr(hilf.x) + sqr(hilf.y)),
      hilf.x /= norm,
      hilf.y /= norm,
      tR[i]  =  tL[i] = hilf;
  }


  for (i = 0; i < m; i++)    /* die Laengen der Parameterintervalle   */
  {                          /* und die Splinekoeffizienten berechnen */
    hilf.x = tR[i].x + tL[i + 1].x;
    hilf.y = tR[i].y + tL[i + 1].y;

    A    = (REAL)16.0 - sqr(hilf.x) - sqr(hilf.y),
    B    = SIX * (s[i].x * hilf.x + s[i].y * hilf.y),
    C    = (REAL)36.0 * (sqr(s[i].x) + sqr(s[i].y)),
    T[i] = (-B + SQRT(sqr(B) + A * C)) / A;

    bx[i]   = tR[i].x;        /* die Subsplinekoeffizienten berechnen */
    by[i]   = tR[i].y;
    TiInv   = 1 / T[i];
    TiInvSq = sqr(TiInv);
    cx[i]   = (THREE * TiInv * s[i].x - TWO * tR[i].x - tL[i + 1].x) *
              TiInv;
    cy[i]   = (THREE * TiInv * s[i].y - TWO * tR[i].y - tL[i + 1].y) *
              TiInv;
    dx[i]   = (hilf.x - TWO * TiInv * s[i].x) * TiInvSq;
    dy[i]   = (hilf.y - TWO * TiInv * s[i].y) * TiInvSq;
  }


  vmfree(vmblock);
  return 0;
}

/* ------------------------- ENDE subsplin.c ------------------------ */
