/*******************errmes.c**********************                  
*                                                *      
*      Vorentwurfs-Generator "VoGen"             *          
*      Version 2.00                              *
*      Date of last modification 21.04.2020      *
*                                                *          
*      Writen by: Jan Himisch                    *      
*             jan.himisch@dlr.de                 *      
*                                                *      
*                                                *
*************************************************/

#include <stdlib.h>
#include <stdio.h>
void errmessage(int Numb)
{
  switch (Numb)
  {
  case 0: 
      {
      printf("\n Irgendwas stimmt hier nicht!! \n");
      exit (EXIT_FAILURE);
      }
  case 1: 
      {
      exit (EXIT_FAILURE);
      }

  case 11: 
      {
      printf("\n Falscher Aufruf!\n ./VoGen parameterfile\n");
      exit (EXIT_FAILURE);
      }
  case 12: 
      {
      printf(" \n Fehler beim �ffnen einer Datei \n ");
      exit (EXIT_FAILURE);
      }
   case 13: 
      {
      printf(" \n Im Parameter-File muss erst der Parameter \" number of targetCL \" gesetzt werden \n");
      printf(" \n bevor ein Ziel CL gesetzt werden darf\n");
      exit (EXIT_FAILURE);
      }
  case 14: 
      {
      printf(" \n Im Parameter-File muss erst der Parameter \" number of target alpha \" gesetzt werden \n");
      printf(" \n bevor ein Ziel Alpha gesetzt werden darf\n");
      exit (EXIT_FAILURE);
      }
  case 15: 
      {
      printf(" \n Im Parameter-File muss erst der Parameter \" Anzahl Fl�gel \" gesetzt werden \n");
      printf(" \n bevor ein Fl�gel Parameter gesetzt werden darf\n");
      exit (EXIT_FAILURE);
      }
  case 16: 
      {
      printf(" \n Im Parameter-File muss erst der Parameter \" Anzahl Kreuzungen \" gesetzt werden \n");
      printf(" \n bevor ein Kreuzungs Paramter gesetzt werden darf\n");
      exit (EXIT_FAILURE);
      }
  case 17:
      {
      printf(" \n Der Parameter  \"Parameters for Optimization\" muss geoeffnet und \n geschlossen werden");
      exit (EXIT_FAILURE);
      }
  case 18:
      {
      printf(" \n Gr�ssere Werte als 9 werden f�r den Parameter  \"Ordnung des Alfa Einflusses (1-4)\" nicht unterst�tzt \n");
      exit (EXIT_FAILURE);
      }
  case 19:
      {
      printf("\n Fehler beim anlegen der Lifting Line Input Datei. (Genuegend SpeciherPlatz?, Schreibrechte, ...)\n");
      exit (EXIT_FAILURE);
      }
  case 20:
      {
      printf("\n Das setzen einer Kreuzung und der Symmetriebedingung im selben Schnitt ist nicht m�glich\n");
      exit (EXIT_FAILURE);
      }
  case 21:
      {
      printf("\n Die Mindestzahl der Definitionswerte f�r die Camber-Line muss in jedem Schnitt mindestens 2 betragen\n");
      exit (EXIT_FAILURE);
      }
  case 22:
      {
      printf("\n Wenn die XML-Aussgabe aktiviert ist muss f�r jeden Schnitt der Profiname angegeben sein! \n");
      exit (EXIT_FAILURE);
      }
  case 23:
      {
      printf("\n Wenn die XML-Aussgabe aktiviert ist muss auch die Dichte, Temperatur und die Geschwindigkeit gesetzt sein! \n");
      exit (EXIT_FAILURE);
      }
  case 24:
      {
      printf("\n Fehler bei der LiftinLine Ausf�hrung. Bitte Lifting-Line-Monitor-File �berpr�fen.\n");
      exit (EXIT_FAILURE);
      }
  case 25:
      {
      printf("\n Fehler bei der Polint Ausf�hrung. Bitte Polint-Output �berpr�fen.\n");
      exit (EXIT_FAILURE);
      }
  case 26:
      {
      printf("\n Fehler beim Einlesen der Parameterdatei.\n");
      exit (EXIT_FAILURE);
      }
  case 27:
      {
      printf("\n Fehler beim Oefnen der Parameterdatei.\n");
      exit (EXIT_FAILURE);
      }
  case 28:
      {
      printf("\n Fehler beim Einlesen der Parameterdatei fuer die Spline Parametriesierung\n");
      exit (EXIT_FAILURE);
      }
  case 29:
      {
      printf("\n Der Eta-Wert bei der Spline-Parametrisierung muss immer zwischen 0 und 1 liegen!!\n");
      exit (EXIT_FAILURE);
      }
  case 30:
      {
      printf("\n Fehler beim erstellen des Akima Splines!\n");
      exit (EXIT_FAILURE);
      }
  case 31:
      {
      printf("\n Fehler beim einlesen des 14-Files. Unter Umstaenden falsche Lifting-Line Version\n Mindestens Release2.3b\n");
      exit (EXIT_FAILURE);
      }
  case 32:
      {
      printf("\n Fehler beim oeffnen des tecplot Loesungs-Files\n");
      exit (EXIT_FAILURE);
      }
  case 33:
      {
      printf("\n  Fehler beim oeffnen des *.14 Loesungs-Files\n");
      exit (EXIT_FAILURE);
      }
  case 34:
      {
      printf("\n  Geschlossener Ringfluegel ist statisch Ueberbestimmt. Keine BM-Bestimmung moeglich!\n");
      exit (EXIT_FAILURE);
      }
  case 35:
      {
      printf("\n Wenn das Einlesen der relativen Dicke aktiviert ist, muss diese auch f�r jeden Schnitt angegeben sein!\n");
      exit (EXIT_FAILURE);
      }
  case 36:
      {
      printf("\n Irgendwas stimmt bei der Bestimmung von CL CQ nicht!! \n");
      exit (EXIT_FAILURE);
      }
  case 37:
      {
      printf(" \n Alle Marker m�ssen in der Form { X_Y_Z} direkt hinter der Option \" Zonen \" angegeben sein.\n");
      exit (EXIT_FAILURE);
      }
  case 38:
      {
      printf(" \n Konnte BSURF-Result-File nicht �ffnen.\n");
      exit (EXIT_FAILURE);
      }
  case 39:
      {
      printf("\n Fehler beim einlesen des BSURF-Result-Files. \n");
      exit (EXIT_FAILURE);
      }
  case 40:
      {
      printf("\n Beim Austauschen der LAstverteilung in der Funktion read_Bsurf passiert etwas was ich nicht verstehe!\n");
      exit (EXIT_FAILURE);
      }
  case 41:
      {
      printf("\n If structure coupled simulations are performed currently only one CL or one alpha can be set!\n");
      exit (EXIT_FAILURE);
      }
  case 42:
      {
      printf("\n Fehler beim erstellen der Datei fxyz4beam.txt f�r die Structurkopplung!\n");
      exit (EXIT_FAILURE);
      }
  case 43:
      {
      printf("\n Der ProjektionsFl�gel muss definiert sein!\n");
      exit (EXIT_FAILURE);
      }
  case 44:
      {
      printf("\n Die jigh2flightdatei passt nicht zur aktuellen Rechnung!\n");
      exit (EXIT_FAILURE);
      }
  case 45:
      {
      printf("\n Konnte jigh2flightdatei nicht oeffnen.\n");
      exit (EXIT_FAILURE);
      }
  case 46:
      {
      printf("\n Bei der Structurkopplung ist das Sichern von mehr als 100 iterationen nicht vorgesehen.\n");
      exit (EXIT_FAILURE);
      }
  case 47:
      {
      printf("\n Kann bsurf ergebniss nicht Einlesen. Egebnisfile �berpr�fen, oder Read Bsurf Load Dist. ausschalte.\n");
      exit (EXIT_FAILURE);
      }
  case 48:
      {
      printf("\n Wenn die Geschwindigkeit f�r die Polintkopplung errechnet werden soll,\nmuss die Dichte die bezugs Fl�sche und das Gewicht angegeben sein!\n");
      exit (EXIT_FAILURE);
      }
  case 49:
      {
      printf("\n Wenn die Geschwindigkeit f�r die Polintkopplung errechnet werden soll,\nd�rfen keine Ziel-Alfa angegeben sein!\n");
      exit (EXIT_FAILURE);
      }
  case 50: 
      {
      printf(" \n Im Parameter-File muss erst der Parameter \" Anzahl Klappen \" gesetzt werden \n");
      printf(" \n bevor ein Klappen Paramter gesetzt werden darf\n");
      exit (EXIT_FAILURE);
      }
  case 51: 
      {
      printf(" \n Fehler beim Setzen der Klappenbedingungen \n");
      exit (EXIT_FAILURE);
      }
  case 52: 
      {
      printf(" \n Die Klappen d�rfen sich nicht �berschneiden\n");
      exit (EXIT_FAILURE);
      }
  case 53: 
      {
      printf(" \n Klappen-Schnitt 2 muss gr��er sein als Klappen-Schnitt 1\n");
      exit (EXIT_FAILURE);
      }
  case 54: 
      {
      printf(" \n Naja die Klappe sollte schon nicht �ber den Fl�gel hinaus ragen!\n");
      exit (EXIT_FAILURE);
      }
  case 55: 
      {
      printf(" \n Innerhalb einer Polare darf dieselbe Kombination von Mach-Zahl, Re-Zahl & Klappenwinkel nur einmal vorkommen\n");
      exit (EXIT_FAILURE);
      }
  case 56: 
      {
      printf(" \n Die Inputdateien  (Name of parametric input File :) fuer\n die Optimierung muessen f�r unterschiedliche Optimierungssaetze unterschiedlich sein!\n");
      exit (EXIT_FAILURE);
      }
  case 57: 
      {
      printf(" \n Die Anzahl der Parametersets im parafile entsprechen nicht der Anzahl der spezifizierten Cas und Alfas!\n");
      exit (EXIT_FAILURE);
      }
  case 58: 
      {
      printf(" \n Es darf nur ein Cl oder Alfa f�r den Vergleich mit der Tau-Lastverteilung gegeben sein!\n");
      exit (EXIT_FAILURE);
      }
  case 59: 
      {
      printf(" \n Beim Spiegeln einer Kreuzung wir das Kreuzen eines nicht gespigelten Fl�gels (Create Sym.(off=0/on=1/newWing=2):0) \n mit einem durch die Option 2 gespiegelten Fl�gel (Create Sym.(off=0/on=1/newWing=2):2)\n nicht unterst�tzt!!!\n");
      exit (EXIT_FAILURE);
      }
  case 60: 
      {
      printf(" \n Seit der Version 1.12 darf in einer Parametergruppe fuer die Optimierung nur Parameter eines Fluegels definiert werden\n Fuer parameter unterschiedlicher Fluegel muessen mehrere Parametersets definiert werden\n");
      exit (EXIT_FAILURE);
      }
  case 61:
      {
      printf("\n Fehler beim Oefnen der CW0Wertdatei.\n");
      exit (EXIT_FAILURE);
      }
  case 62:
      {
      printf("\n Fehler beim Einlesen der CW0Wertdatei.\n");
      exit (EXIT_FAILURE);
      }
  case 63:
      {
      printf("\n Zur Zeit wird nur die LiftingLine Version 2.3beta(V2p3b) oder die entsprechende Kreisflug-Modifikation(V2p3bCirc) unterstuezt\n");
      exit (EXIT_FAILURE);
      }
  case 64:
      {
      printf("\n Der Kreisflugmodus wird nur mit der (V2p3bCirc) unterstuezt\n");
      exit (EXIT_FAILURE);
      }
  case 65:
      {
      printf("\n Currently only the x-variable can be defined as dev varaible in the mixer mode\n");
      exit (EXIT_FAILURE);
      }
  case 66:
      {
      printf("\n Currently only the flaps indicatet by f are aloud as targets of mixer-functions\n");
      exit (EXIT_FAILURE);
      }
  case 67:
      {
      printf("\n Currently only Lifting Line or Free Wake can be started not both at same time!\n");
      exit (EXIT_FAILURE);
      }
  case 68:
      {
      printf("\n Near Lifting Line nor Free Wake Start >0 --> nothing to do!\n");
      exit (EXIT_FAILURE);
      }
  case 69:
      {
      printf("\n To run one-by-one for more than one CL + Alpha is only feasible if a new input file can be generated (don't use old input)\n");
      exit (EXIT_FAILURE);
      }
  case 70:
      {
      printf("\n Error while creating FreeWake Input File (DiskSpace?, disk-rights?, ...)\n");
      exit (EXIT_FAILURE);
      }
  case 71:
      {
      printf("\n Number of airfoils is not allowed to be larger than 15 if using FreeWake\n");
      exit (EXIT_FAILURE);
      }
  }  
}
