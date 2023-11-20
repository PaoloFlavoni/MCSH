/* nyquist.c : file sorgente di un programma che che studia un paio di
               segnali complessi dipendenti dal tempo e relativi al
               moto di un sistema dinamico autonomo (si intende che i
               segnali sono discretizzati ad istanti di tempo
               equidistanziati, quindi sono descritti da due
               successioni finite di valori nel piano complesso).
	       Lo scopo di questo programma e' quello di scrivere su
	       due file esterni i dati necessari a tracciare i grafici
	       (per entrambi i segnali) dei valori assoluti
	       dell'integrale fondamentale dell'analisi in frequenza
	       in funzione della velocita' angolare. La variabile
	       indipendente, cioe' la velocita' angolare, viene fatta
	       variare per valori discretizzati e equidistanziati
	       compresi tra meno la frequenza di Nyquist e la frequenza
	       di Nyquist stessa, con segno positivo.
	       Modificando opportunamente i parametri che definiscono
	       i segnali da analizzare (creati da un programma che
	       deve essere eseguito preliminarmente rispetto a quello
	       presente in questo file), questo programma puo'
	       evidenziare gli effetti di aliasing. */
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
  int npt_graph;                 /* Numero dei punti che vogliamo mettere in
				    grafico */
  int i_graph;                   /* Indice che funge da contatore sul
				    numero dei punti che vogliamo
				    mettere in grafico */
  double nu;                     /* Valore della velocita' angolare
				    rispetto alla quale si calcola
				    l'integrale fondamentale
				    dell'analisi in frequenza */
  double valint[2];              /* Valore (nel piano complesso)
				    dell'integrale fondamentale
				    dell'analisi in frequenza */
  FILE *in;                      /* Puntatore al file di input */
  FILE *out;                     /* Puntatore al file di output */
  char nomefile[30];             /* Stringa dove viene scritto il nome
				    del/dei file di input */
  /* Inizio delle istruzioni eseguibili */
  pig = 4. * atan(1.);
  dupig = 2. * pig;
  /* Apertura del file di dati concernenti alcune informazioni generali
     riguardanti l'analisi in frequenza che vogliamo effettuare */
  in = fopen("nyquist.dat","r");
  fscanf(in,"%ld", &npt);
  fscanf(in,"%le", &tstep);
  fscanf(in,"%d", &npt_graph);
  fclose(in);  /* Chiusura del file dati */
  /* Definizione di tgrande.
     Si osservi che se il sistema in considerazione e' autonomo
     (=senza dipendenza esplicita dal tempo nella Hamiltoniana
     oppure nella mappa simplettica), allora possiamo studiare
     l'integrale fondamentale dell'analisi in frequenza
     nell'intervallo tra - tgrande e + tgrande anche quando l'orbita
     e' conosciuta tra t0 - tgrande e t0 + tgrande con t0 diverso da
     0. */
  tgrande = .5 * npt * tstep;

  /* Inizio del ciclo sul numero di segnali */
 
    /* Apertura del file di dati contenente il
       segnale corrispondenti all'orbita che vogliamo studiare con
       l'analisi in frequenza */
    sprintf(nomefile, "orbite_1.dat");
    in = fopen(nomefile,"r");
    sprintf(nomefile, "nyquist_modello_8192.out");
    out = fopen(nomefile,"w");
    printf("       Produzione del file %s\n", nomefile);
    /* Lettura delle npt + 1 coppie di dati che costituiscono il
       sign_count-esimo segnale */
    leggi_segnale(in, npt, 0);
    fclose(in);  /* Chiusura del file dati */
    for(i_graph=0; i_graph<=npt_graph; i_graph++) {
      /* Calcolo della velocita' angolare nu */
      nu = - pig / tstep + dupig / (tstep * npt_graph) * i_graph;
      /* Calcolo del valore di A(nu), cosi' come definito
	 dall'integrale fondamentale dell'analisi in frequenza; tale
	 calcolo viene effettuato con un'opportuna chiamata alla
	 funzione prodscal, che e' inclusa nel file lib_ainf.c */
      prodscal(nu, valint);
      /* Stampa su file esterno dei valori nu e di |A(nu)| */
      fprintf(out, "  %le  %le\n",
	      nu, sqrt(valint[0] * valint[0] + valint[1] * valint[1]) );
    }      
    fclose(out);  /* Chiusura del file di output */
  
}
