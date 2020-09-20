/* ------------------------- MODUL vmblock.c ------------------------ */

/***********************************************************************
*                                                                      *
* Verwaltung eines Satzes von dynamischen Vektoren und Matrizen        *
* -------------------------------------------------------------        *
*                                                                      *
* Idee:   In vielen Unterprogrammen der Numerikbibliothek werden immer *
*         wieder dynamisch vereinbarte Vektoren und Matrizen           *
*         benoetigt. Dabei tritt jedoch manchmal das Problem auf, dass *
*         nur fuer einen Teil der Vektoren und Matrizen Speicher       *
*         vorhanden ist, so dass der schon belegte Speicher            *
*         zurueckgegeben und auf den Speichermangel geeignet reagiert  *
*         werden muss. Dies kostet viel Muehe und stellt eine haeufige *
*         Fehlerquelle dar, wenn man es jedesmal neu formulieren muss. *
*         Zur Vereinfachung dieser Arbeit wurde daher dieses C-Modul   *
*         geschrieben. Es verwaltet alle zusammengehoerigen            *
*         Speicheranforderungen fuer Vektoren und Matrizen in einer    *
*         einfach verketteten Liste. Dazu werden dem Benutzer folgende *
*         vier Funktionen zur Verfuegung gestellt:                     *
*                                                                      *
*         - vminit(),    das einen typlosen Listenanfangszeiger        *
*                        liefert, mit dem alle weiteren Funktionen     *
*                        arbeiten,                                     *
*                                                                      *
*         - vmalloc()    zur Anforderung von Speicher fuer einen neuen *
*                        Vektor oder eine neue Matrix,                 *
*                                                                      *
*         - vmcomplete() zur nachtraeglichen Pruefung, ob alle         *
*                        bisherigen in der Liste vorgenommenen         *
*                        Speicheranforderungen zum Erfolg fuehrten,    *
*                        und                                           *
*                                                                      *
*         - vmfree(),    das den von einer Vektor-Matrix-Liste         *
*                        beanspruchten Speicher wieder vollstaendig    *
*                        freigibt.                                     *
*                                                                      *
*         Ausserdem werden noch die fuenf Makros                       *
*                                                                      *
*         - VEKTOR  (fuer REAL-Vektoren),                              *
*         - VVEKTOR (fuer variable Vektoren),                          *
*         - MATRIX  (fuer REAL-Matrizen),                              *
*         - IMATRIX (fuer int-Matrizen) und                            *
*         - PMATRIX (fuer Punktmatrizen im R3) exportiert,             *
*                                                                      *
*         mit denen der Benutzer beim Aufruf von vmalloc() den Typ der *
*         anzufordernden Datenstruktur waehlen kann.                   *
*                                                                      *
*         Achtung: 1. Der von einer Vektor-Matrix-Liste                *
*                     beanspruchte Speicher darf nur durch vmfree()    *
*                     freigegeben werden!                              *
*                  2. vmfree() gibt immer nur den gesamten schon       *
*                     angeforderten Speicher frei, der zu einer Liste  *
*                     gehoert, laesst sich also nicht auf einzelne     *
*                     Vektoren oder Matrizen der Liste anwenden!       *
*                                                                      *
* Aufruf: Der Benutzer vereinbart einen typlosen Zeiger, der           *
*         zuallererst durch einen Aufruf von vminit() initialisiert    *
*         werden muss und von da an den einzigen gueltigen Zugang zur  *
*         Speicherliste darstellt. Ueber diesen Zeiger koennen nun mit *
*         Hilfe von vmalloc() Vektoren und Matrizen dynamisch angelegt *
*         werden. Wurden alle Speicheranforderungen getaetigt, sollte  *
*         man mit vmcomplete() pruefen, ob sie auch gelungen sind, und *
*         dann entsprechend reagieren. Wenn die zur Liste gehoerenden  *
*         Vektoren und Matrizen nicht mehr benoetigt werden, empfiehlt *
*         es sich, denn davon beanspruchten Speicher durch Aufruf von  *
*         vmfree() der Allgemeinheit wieder zur Verfuegung zu stellen. *
*         Beispiel:                                                    *
*             ...                                                      *
*             void *vmblock;                                           *
*             REAL *vektor1;    /+ REAL-Vektor mit n1 Elementen     +/ *
*             REAL *vektor2;    /+ REAL-Vektor mit n2 Elementen     +/ *
*             int  *vektor3;    /+ int-Vektor mit n3 Elementen      +/ *
*             REAL **matrix1;   /+ Matrix mit m1 Zeilen, n1 Spalten +/ *
*             int  **matrix2;   /+ Matrix mit m2 Zeilen, n2 Spalten +/ *
*             REAL ***pmatrix;  /+ Matrix mit m2*n2 Punkten im R3   +/ *
*             ...                                                      *
*             vmblock = vminit();                                      *
*             vektor1 = (REAL *)vmalloc(vmblock, VEKTOR,  n1, 0);      *
*             vektor2 = (REAL *)vmalloc(vmblock, VEKTOR,  n2, 0);      *
*             vektor2 = (int *) vmalloc(vmblock, VVEKTOR, n2,          *
*                                       sizeof(int));                  *
*             ...                                                      *
*             matrix1 = (REAL **) vmalloc(vmblock, MATRIX,  m1, n1);   *
*             matrix2 = (int  **) vmalloc(vmblock, IMATRIX, m2, n2);   *
*             pmatrix = (REAL ***)vmalloc(vmblock, PMATRIX, m2, n2);   *
*             ...                                                      *
*             if (! vmcomplete(vmblock))  /+ teilweise misslungen? +/  *
*             {                                                        *
*               vmfree(vmblock);          /+ Block ganz freigeben  +/  *
*               return 99;                /+ Fehler melden         +/  *
*             }                                                        *
*             ...                                                      *
*             vmfree(vmblock);                                         *
*             ...                                                      *
*                                                                      *
* Programmiersprache: ANSI C                                           *
* Compiler:           Borland C++ 2.0                                  *
* Rechner:            IBM PS/2 70 mit 80387                            *
* Autor:              Juergen Dietel, Rechenzentrum der RWTH Aachen    *
* Datum:              DO 10. 9. 1992                                   *
*                                                                      *
***********************************************************************/

#include "basis.h"     /* wegen size_t, NULL, malloc, free, calloc,   */
                       /*       boolean, FALSE, TRUE, REAL            */
#include "vmblock.h"   /* wegen vmalloc, vmcomplete, vmfree, vminit,  */
                       /*       VEKTOR, VVEKTOR, MATRIX, IMATRIX,     */
                       /*       PMATRIX                               */



/*--------------------------------------------------------------------*/

typedef struct VML           /* Element einer Vektor-Matrix-Liste     */
{
  void       *vmzeiger;      /* Zeiger auf den Vektor bzw. die Matrix */
  int        typ;            /* Typ des Zeigers: Vektor oder Matrix   */
                             /* (moegliche Werte: VEKTOR, VVEKTOR,    */
                             /*                   MATRIX, IMATRIX,    */
                             /*                   PMATRIX)            */
  size_t     groesse;        /* im Ankerelement die Anzahl der fuer   */
                             /* Vektoren oder Matrizen genutzten      */
                             /* Listenelemente, sonst ungenutzt       */
                             /* ausser bei Matrizen, wo groesse fuer  */
                             /* die Zeilenanzahl "missbraucht" wird   */
  size_t     spalten;        /* Spaltenanzahl bei Punktmatrizen       */
  struct VML *naechst;       /* Zeiger auf das naechste Listenelement */
} vmltyp;

#define VMALLOC  (vmltyp *)malloc(sizeof(vmltyp))    /* Speicher fuer */
                                                     /* ein neues     */
                                                     /* Listenelement */
                                                     /* anfordern     */

#define LISTE    ((vmltyp *)vmblock)              /* zur Abkuerzung   */
                                                  /* der Schreibweise */
#define MAGIC    410      /* soll ein gueltiges Ankerelement anzeigen */



/*--------------------------------------------------------------------*/

void *vminit(void    /* eine leere Vektor-Matrix-Liste erzeugen ......*/
/*.IX{vminit}*/
            )                   /* Adresse der Liste .................*/

/***********************************************************************
* eine leere Vektor-Matrix-Liste erzeugen. Diese besteht aus einem     *
* Ankerelement, das nur dazu benoetigt wird, die Anzahl der            *
* angeforderten Listenelemente und einen magischen Wert fuer           *
* Plausibilitaetskontrollen aufzunehmen. So kann man spaeter in        *
* vmcomplete() pruefen, ob der Funktionswert NULL von vmalloc() daher  *
* kommt, dass die Liste nicht mehr verlaengert werden konnte oder dass *
* fuer den Vektor bzw. die Matrix kein Speicher mehr vorhanden war.    *
* Als Funktionswert wird die Adresse des Ankers zurueckgegeben, im     *
* Fehlerfall natuerlich NULL. Um in den nachfolgenden Aufrufen von     *
* vmalloc(), vmcomplete() und vmfree() pruefen zu koennen, ob der      *
* uebergebene typlose Zeiger wirklich auf eine Vektor-Matrix-Liste     *
* zeigt, wird die Komponente `typ' des Ankerelements dazu missbraucht, *
* einen magischen Wert aufzunehmen, der ein gueltiges Ankerelement     *
* anzeigen soll.                                                       *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* vmltyp, VMALLOC, MAGIC, NULL, malloc                                 *
***********************************************************************/

{
  vmltyp *liste;             /* Zeiger auf das Ankerelement der Liste */


  if ((liste = VMALLOC) == NULL) /* Speicher fuer den Anker anfordern */
    return NULL;                 /* misslungen? => Fehler melden      */
  liste->vmzeiger = NULL;        /* damit vmfree() sich nicht vertut  */
  liste->typ      = MAGIC;       /* einen gueltigen Anker anzeigen    */
  liste->groesse  = 0;           /* noch keine Listenelemente da      */
  liste->naechst  = NULL;        /* noch kein Nachfolger              */


  return (void *)liste;
}



/*--------------------------------------------------------------------*/

static REAL **malloc_realmatrix
                               (
                                size_t m,
                                size_t n
                               )

/***********************************************************************
* Speicherplatz fuer eine rechteckige [0..m-1,0..n-1]-Matrix mit       *
* Elementen vom Typ REAL anfordern und ihre Anfangsadresse als         *
* Funktionswert zurueckgeben, falls die Anforderung zum Erfolg         *
* fuehrte, sonst NULL. Dabei wird fuer jede Zeile der Matrix ein       *
* eigener Zeiger angelegt.                                             *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* REAL, size_t, NULL, calloc, free                                     *
***********************************************************************/

{
  REAL   **matrix;                   /* Zeiger auf die Zeilenvektoren */
  size_t i;                          /* laufender Zeilenindex         */


  matrix = (REAL **)                        /* fuer jede der m Zeilen */
           calloc(m, sizeof(*matrix));      /* einen Zeiger           */

  if (matrix == NULL)           /* nicht genug Speicher?              */
    return NULL;                /* Speichermangel melden              */

  for (i = 0; i < m; i++)       /* fuer jeden Zeilenzeiger eine Zeile */
  {                             /* mit je n Elementen anfordern       */
    matrix[i] = (REAL *)
                calloc(n, sizeof(**matrix));
    if (matrix[i] == NULL)      /* nicht genug Speicher?              */
    {
      while (i != 0)            /* schon reservierte Zeilen freigeben */
        free(matrix[--i]);
      free(matrix);             /* Vektor der Zeilenzeiger freigeben  */
      return NULL;              /* Speichermangel melden              */
    }
  }


  return matrix;
}



/*--------------------------------------------------------------------*/

static int **malloc_intmatrix
                             (
                              size_t m,
                              size_t n
                             )

/***********************************************************************
* Speicherplatz fuer eine rechteckige [0..m-1,0..n-1]-Matrix mit       *
* Elementen vom Typ int anfordern und ihre Anfangsadresse als          *
* Funktionswert zurueckgeben, falls die Anforderung zum Erfolg         *
* fuehrte, sonst NULL. Dabei wird fuer jede Zeile der Matrix ein       *
* eigener Zeiger angelegt.                                             *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* size_t, NULL, calloc, free                                           *
***********************************************************************/

{
  int    **matrix;                   /* Zeiger auf die Zeilenvektoren */
  size_t i;                          /* laufender Zeilenindex         */


  matrix = (int **)                         /* fuer jede der m Zeilen */
           calloc(m, sizeof(*matrix));      /* einen Zeiger           */

  if (matrix == NULL)           /* nicht genug Speicher?              */
    return NULL;                /* Speichermangel melden              */

  for (i = 0; i < m; i++)       /* fuer jeden Zeilenzeiger eine Zeile */
  {                             /* mit je n Elementen anfordern       */
    matrix[i] = (int *)
                calloc(n, sizeof(**matrix));
    if (matrix[i] == NULL)      /* nicht genug Speicher?              */
    {
      while (i != 0)            /* schon reservierte Zeilen freigeben */
        free(matrix[--i]);
      free(matrix);             /* Vektor der Zeilenzeiger freigeben  */
      return NULL;              /* Speichermangel melden              */
    }
  }


  return matrix;
}



/*--------------------------------------------------------------------*/

static void free_matrix
                       (
                        void   **matrix,
                        size_t m
                       )

/***********************************************************************
* eine wie in malloc_realmatrix() oder malloc_intmatrix() erzeugte     *
* Matrix mit m Zeilen freigeben                                        *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* size_t, NULL, free                                                   *
***********************************************************************/

{
  if (matrix != NULL)
  {
    while (m != 0)
      free(matrix[--m]);
    free(matrix);
  }
}



/*--------------------------------------------------------------------*/

static void free_pmatrix
                        (
                         void   ***matrix,
                         size_t m,
                         size_t n
                        )

/***********************************************************************
* eine wie in malloc_punktmatrix() erzeugte Matrix mit m Zeilen und    *
* n Spalten freigeben                                                  *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* size_t, NULL, free, free_matrix                                      *
***********************************************************************/

{
  if (matrix != NULL)
  {
    while (m != 0)
      free_matrix(matrix[--m], n);
    free(matrix);
  }
}



/*--------------------------------------------------------------------*/

static REAL ***malloc_punktmatrix
                                 (
                                  size_t m,
                                  size_t n
                                 )

/***********************************************************************
* Speicherplatz fuer eine [0..m-1,0..n-1,0..2]-Matrix mit Elementen    *
* vom Typ REAL anfordern und ihre Anfangsadresse als Funktionswert     *
* zurueckgeben, falls die Anforderung zum Erfolg fuehrte, sonst NULL.  *
* Dabei wird fuer jede Zeile der Matrix ein eigener Zeiger angelegt.   *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* size_t, NULL, calloc, free_matrix, free                              *
***********************************************************************/

{
  REAL   ***matrix;                  /* Zeiger auf die Zeilenvektoren */
  size_t i;                          /* laufender Zeilenindex         */


  matrix = (REAL ***)                       /* fuer jede der m Zeilen */
           calloc(m, sizeof(*matrix));      /* einen Zeiger           */

  if (matrix == NULL)         /* nicht genug Speicher?                */
    return NULL;              /* Speichermangel melden                */

  for (i = 0; i < m; i++)     /* fuer jeden Zeilenzeiger eine         */
  {                           /* (n,3)-Matrix anfordern               */
    matrix[i] = malloc_realmatrix(n, 3);
    if (matrix[i] == NULL)    /* nicht genug Speicher?                */
    {
      while (i != 0)          /* schon reservierte Matrizen freigeben */
        free_matrix((void **)(matrix[--i]), n);
      free(matrix);           /* Vektor der Zeilenzeiger freigeben    */
      return NULL;            /* Speichermangel melden                */
    }
  }


  return matrix;
}



/*--------------------------------------------------------------------*/

void *vmalloc        /* dynamischen Vektor bzw. Matrix erzeugen ......*/
/*.IX{vmalloc}*/
             (
              void   *vmblock,  /* Adresse einer Vektor-Matrix-Liste  */
              int    typ,       /* Art des Vektors/der Matrix ........*/
              size_t zeilen,    /* Elementanzahl/Zeilenanzahl ........*/
              size_t spalten    /* Spaltenanzahl/Elementgroesse ......*/
             )                  /* Adresse des geschaffenen Objekts ..*/

/***********************************************************************
* ein durch `typ' bestimmtes Element (Vektor oder Matrix), dessen      *
* Groesse durch `zeilen' und `spalten' festgelegt wird, erzeugen und   *
* vorne in die bei `vmblock' beginnende einfach verkettete Liste       *
* einfuegen. Die Adresse des neuen Vektors bzw. der neuen Matrix wird  *
* als Funktionswert zurueckgegeben. Bei einem REAL-Vektor (Typ VEKTOR) *
* enthaelt der Parameter `zeilen' die Laenge, `spalten' wird nicht     *
* benutzt. Beim Typ VVEKTOR (variabler Vektor) muss in `spalten' die   *
* Groesse eines einzelnen Vektorelements stehen. Bei einer Matrix (Typ *
* MATRIX, IMATRIX oder PMATRIX) enthaelt `zeilen' die Zeilen- und      *
* `spalten' die Spaltenanzahl der Matrix.                              *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* vmltyp, VMALLOC, LISTE, MAGIC, malloc_realmatrix, malloc_intmatrix,  *
* malloc_punktmatrix, REAL, VEKTOR, VVEKTOR, MATRIX, IMATRIX, PMATRIX, *
* NULL, size_t, malloc, calloc                                         *
***********************************************************************/

{
  vmltyp *element;               /* Zeiger auf das neue Listenelement */


  if (LISTE == NULL)                 /* ungueltige Liste?             */
    return NULL;                     /* Fehler melden                 */

  if (LISTE->typ != MAGIC)           /* ungueltiges Ankerelement?     */
    return NULL;                     /* Fehler melden                 */


  LISTE->groesse++;                  /* neue Anforderung registrieren */


  if ((element = VMALLOC) == NULL)   /* neues Listenelement anfordern */
    return NULL;                     /* misslungen? => Fehler melden  */

  switch (typ)         /* Speicher fuer die gewuenschte Datenstruktur */
  {                    /* anfordern (Vektor oder Matrix) und ihre     */
                       /* Adresse in das neue Listenelement eintragen */

    case VEKTOR:          /* ---------- REAL-Vektor?       ---------- */
      element->vmzeiger = calloc(zeilen, sizeof(REAL));
      break;

    case VVEKTOR:         /* ---------- beliebiger Vektor? ---------- */
      element->vmzeiger = calloc(zeilen, spalten);
      break;

    case MATRIX:          /* ---------- REAL-Matrix?       ---------- */
      element->vmzeiger = (void *)malloc_realmatrix(zeilen, spalten);
      element->groesse  = zeilen;      /* fuer vmfree() unter groesse */
      break;                           /* die Zeilenanzahl eintragen  */

    case IMATRIX:         /* ---------- int-Matrix?        ---------- */
      element->vmzeiger = (void *)malloc_intmatrix(zeilen, spalten);
      element->groesse  = zeilen;      /* fuer vmfree() unter groesse */
      break;                           /* die Zeilenanzahl eintragen  */

    case PMATRIX:         /* ---------- Punktmatrix?       ---------- */
      element->vmzeiger = (void *)malloc_punktmatrix(zeilen, spalten);
      element->groesse  = zeilen;      /* fuer vmfree() unter groesse */
      element->spalten  = spalten;     /* und spalten die Zeilen- und */
      break;                           /* Spaltenanzahl eintragen     */
    default:              /* ---- ungueltiger Datenstrukturtyp? ----  */
      element->vmzeiger = NULL;        /* Nullzeiger eintragen        */
  }

  element->typ = typ;                  /* Datenstrukturtyp im         */
                                       /* Listenelement notieren      */
  element->naechst = LISTE->naechst;   /* neues Element einfuegen vor */
                                       /* dem ersten Element und ...  */

  LISTE->naechst = element;            /* ... hinter dem Ankerelement */


  return element->vmzeiger;            /* neue Vektor/Matrix-Adresse  */
}                                      /* zurueckgeben                */



/*--------------------------------------------------------------------*/

boolean vmcomplete   /* Vektor-Matrix-Liste auf Speichermangel testen */
/*.IX{vmcomplete}*/
                  (
                   void *vmblock  /* Adresse der Liste ...............*/
                  )               /* kein Speichermangel? ............*/

/***********************************************************************
* die bei `vmblock' beginnende Liste durchlaufen und dabei pruefen, ob *
* allen darin notierten Vektoren und Matrizen Speicher zugeordnet      *
* werden konnte und ob fuer alle angeforderten Vektoren bzw. Matrizen  *
* ueberhaupt Listenelemente angelegt werden konnten. Der Funktionswert *
* zeigt an, ob diese Pruefung erfolgreich war (TRUE) oder nicht        *
* (FALSE).                                                             *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* vmltyp, LISTE, MAGIC, boolean, FALSE, TRUE, NULL, size_t             *
***********************************************************************/

{
  vmltyp *element;       /* Zeiger zum Durchlaufen der Liste          */
  size_t elementzahl;    /* Zahl der Listenelemente, die tatsaechlich */
                         /* angelegt werden konnten (im Gegensatz zu  */
                         /* denjenigen, die angefordert wurden und    */
                         /* deren Anzahl im Ankerelement der Liste    */
                         /* steht, und zwar in `groesse')             */


  if (LISTE == NULL)                /* ungueltige Liste?              */
    return FALSE;                   /* Fehler melden                  */

  if (LISTE->typ != MAGIC)          /* ungueltiges Ankerelement?      */
    return FALSE;                   /* Fehler melden                  */


  for (elementzahl = 0,             /* die Liste durchlaufen und      */
       element = LISTE->naechst;    /* dabei mitzaehlen, wieviele     */
       element != NULL;             /* Elemente sie wirklich enthaelt */
       element = element->naechst,
       elementzahl++)
    if (element->vmzeiger == NULL) /* Element ohne Speicherzuordnung? */
      return FALSE;                /* Misserfolg melden               */

  if (elementzahl != LISTE->groesse)  /* Konnten nicht genuegend      */
                                      /* Elemente erzeugt werden?     */
    return FALSE;                     /* Misserfolg melden            */


  return TRUE;                      /* erfolgreiche Pruefung melden   */
}



/*--------------------------------------------------------------------*/

void vmfree          /* Speicher einer Vektor-Matrix-Liste freigeben  */
/*.IX{vmfree}*/
           (
            void *vmblock       /* Adresse der Liste .................*/
           )

/***********************************************************************
* saemtlichen dynamischen Speicher der bei `vmblock' beginnenden Liste *
* freigeben                                                            *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* vmltyp, LISTE, MAGIC, free_matrix, free_pmatrix, VEKTOR, VVEKTOR,    *
* MATRIX, IMATRIX, PMATRIX, NULL, free                                 *
***********************************************************************/

{
  vmltyp *hilf;             /* Zwischenspeicher fuer einen Zeigerwert */


  if (LISTE == NULL)                     /* ungueltige Liste?         */
    return;                              /* nichts tun                */

  if (LISTE->typ != MAGIC)               /* ungueltiges Ankerelement? */
    return;                              /* nichts tun                */


  for ( ; LISTE != NULL; vmblock = (void *)hilf)
  {

    switch (LISTE->typ)
    {
      case VEKTOR:
      case VVEKTOR: if (LISTE->vmzeiger != NULL)
                      free(LISTE->vmzeiger);
                    break;
      case MATRIX:
      case IMATRIX: free_matrix((void **)LISTE->vmzeiger,
                                LISTE->groesse);
                    break;
      case PMATRIX: free_pmatrix((void ***)LISTE->vmzeiger,
                                 LISTE->groesse, LISTE->spalten);
    }

    hilf = LISTE->naechst;               /* Nachfolgerzeiger retten   */
    free(LISTE);                         /* Listenelement freigeben   */
  }
}

/* ------------------------- ENDE vmblock.c ------------------------- */
