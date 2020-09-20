


/*.FE{P 0.2}{Dynamische Vektoren und Matrizen}
            {Dynamische Vektoren und Matrizen}*/

/* --------------------- DEKLARATIONEN vmblock.h -------------------- */

#define VEKTOR   0                 /* fuer einen REAL-Vektor          */
/*.IX{VEKTOR}*/
#define VVEKTOR  1                 /* fuer einen Vektor mit Elementen */
/*.IX{VVEKTOR}*/
                                   /* von angegebener Groesse         */
#define MATRIX   2                 /* fuer eine REAL-Matrix           */
/*.IX{MATRIX}*/
#define IMATRIX  3                 /* fuer eine int-Matrix            */
/*.IX{IMATRIX}*/
#define PMATRIX  4                 /* fuer eine Punktmatrix im R3     */
/*.IX{PMATRIX}*/


void *vminit(void    /* eine leere Vektor-Matrix-Liste erzeugen ......*/
            );                  /* Adresse der Liste .................*/


void *vmalloc        /* dynamischen Vektor bzw. Matrix erzeugen ......*/
             (
              void   *vmblock,  /* Adresse einer Vektor-Matrix-Liste  */
              int    typ,       /* Art des Vektors/der Matrix ........*/
              size_t zeilen,    /* Elementanzahl/Zeilenanzahl ........*/
              size_t spalten    /* Spaltenanzahl/Elementgroesse ......*/
             );                 /* Adresse des geschaffenen Objekts ..*/


boolean vmcomplete   /* Vektor-Matrix-Liste auf Speichermangel testen */
                  (
                   void *vmblock  /* Adresse der Liste ...............*/
                  );              /* kein Speichermangel? ............*/


void vmfree          /* Speicher einer Vektor-Matrix-Liste freigeben  */
           (
            void *vmblock       /* Adresse der Liste .................*/
           );

/* ------------------------- ENDE vmblock.h ------------------------- */
