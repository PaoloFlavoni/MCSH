 /*PROGRAMMA ADATTATO ALLO STUDIO DI UN PENDOLO FORZATO*/
 /* Programma che studia un paio di segnale complessi dipendenti dal
   tempo e relativi al moto di un sistema dinamico autonomo (si
   intende che i segnali sono discretizzati ad istanti di tempo
   equidistanziati, quindi sono descritti da due successioni finite di
   valori nel piano complesso).
   Questo programma determina le principali componenti periodiche dei
   segnali utilizzando il metodo dell'analisi in frequenza.  I valori
   complessi delle successioni finite relativi ai segnali vengono
   letti da un paio di file esterni che devono essere scritti da un
   altro programma eseguito precedentemente a quello presente in
   questo file.
   Un esempio di programma che scrive i segnali (relativi al modello
   di Henon-Heiles) e' contenuto nel file scrivisegnale_HH.c */
#include <stdio.h>
#include <math.h>
/* Inclusione del file che contiene le definizioni generali */
#include "lib_ainf.h"
/* Inclusione del file che contiene tutte le function necessarie per
   implementare l'analisi in frequenza */
#include "lib_ainf.c"
/* ### */
int main()
{
  int j;                         /* Indice che funge da contatore all'interno
				    dei cicli for */
  int sign_count;                /* Indice che funge da contatore all'interno
				    del grande ciclo for che distingue tra i
				    segnali da analizzare */
  unsigned long int i;           /* Indice che funge da contatore all'interno
				    dei cicli for */
  double omega[NFREQMAX];        /* omega[0], ..., omega[nfreqstudio-1] sono
				    le velocita' angolari corrispondenti ai
				    massimi relativi dell'ampiezza
				    dell'integrale fondamentale dell'analisi
				    in frequenza */
  FILE *in;                      /* Puntatore al file di input */
  char nomefile[30];             /* Stringa dove viene scritto il nome
				    del/dei file di input */
  /* Inizio delle istruzioni eseguibili */
  pig = 4. * atan(1.);
  dupig = 2. * pig;
  /* Apertura del file di dati concernenti alcune informazioni generali
     riguardanti l'analisi in frequenza che vogliamo effettuare */
  in = fopen("freqmax_PF.dat","r");
  fscanf(in,"%ld", &npt);
  fscanf(in,"%le", &tstep);
  fclose(in);  /* Chiusura del file dati */
  /* Definizione di tgrande.
     Si osservi che se il sistema in considerazione e' autonomo
     (=senza dipendenza esplicita dal tempo nella Hamiltoniana
     oppure nella mappa simplettica), allora possiamo studiare
     l'integrale fondamentale dell'analisi in frequenza
     nell'intervallo tra - tgrande e + tgrande anche quando l'orbita
     e' conosciuta tra t0 - tgrande e t0 + tgrande con t0 diverso da
     0. */
  tgrande = 0.5 * npt * tstep; 
    /* Apertura del file di dati contenente il
       segnale corrispondenti all'orbita che vogliamo studiare con
       l'analisi in frequenza */
    sprintf(nomefile, "pm_pendolo_eulero.out");
    in = fopen(nomefile,"r");
    printf("Analisi del file %s\n", nomefile);
    /* Lettura delle npt + 1 coppie di dati che costituiscono il
        segnale */
    leggi_segnale(in, npt, 0);
    fclose(in);  /* Chiusura del file dati */
    /* Calcolo della velocita' angolare corrispondente al massimo
       assoluto dell'ampiezza del sign_count-esimo integrale
       fondamentale dell'analisi in frequenza.  */
    omega[0] = max_segnale(TRUE, ERRSOLUZ);
    /* Operazioni di stampa dei risultati riguardo al punto di massimo
       del sign_count-esimo integrale fondamentale dell'analisi in
       frequenza. */
    printf("Indice          omega         \n");
    printf("  %3d %24.17le \n\n", 0, omega[0]);
  
}