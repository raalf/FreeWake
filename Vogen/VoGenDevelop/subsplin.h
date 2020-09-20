/*.KA{P 13}{Akima-- und Renner--Subsplines}
           {Akima-- und Renner--Subsplines}*/
/* -------------------- DEKLARATIONEN subsplin.h -------------------- */

int akima(int  *n,                 /* Nummer des letzten Stuetzpunkts */
          int  nmax,               /* obere Indexgrenze fuer x, y,... */
          REAL x[],                /* Stuetzstellen                   */
          REAL y[],                /* Stuetzwerte                     */
          int  perio,              /* periodische Interpolation?      */
          REAL beta,               /* Mass fuer die Eckenrundung      */
          REAL b[],                /* Splinekoeffizienten             */
          REAL c[],
          REAL d[]
         );                        /* Fehlercode                      */

int renner(int  *n,            /* Nummer des letzten Stuetzpunkts     */
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
          );                   /* Fehlercode                          */

/* ------------------------- ENDE subsplin.h ------------------------ */
