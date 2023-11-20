/* Programma che studia un segnale complesso dipendente dal
   tempo e relativo al moto di un sistema dinamico autonomo (si
   intende che il segnale è discretizzato ad istanti di tempo
   equidistanziati, quindi è descritto da una successione finita di
   valori nel piano complesso).
   Questo programma determina le principali componenti periodiche del
   segnale utilizzando il metodo dell'analisi in frequenza. I valori
   complessi della successione finita relativa al segnale vengono
   letti da un file esterno che devo essere scritto da un
   altro programma eseguito precedentemente a quello presente in
   questo file.
   IL SEGUENTE PROGRAMMA È STATO MODIFICATO PER IL CASO SPECIFICO DEL 
   PENDOLO FORZATO. */
#include <stdio.h>
#include <math.h>
/* Inclusione del file che contiene le definizioni generali */
#include "lib_ainf.h"
/* Inclusione del file che contiene tutte le function necessarie per
   implementare l'analisi in frequenza */
#include "lib_ainf.c"
double segnale_orig[2*NPTMAX+2];
/* ### */
int main()
{
  int j;                         /* Indice che funge da contatore all'interno
				    dei cicli for */
  unsigned long int i;           /* Indice che funge da contatore all'interno
				    dei cicli for */
  int nfreqstudio;               /* Numero delle componenti del segnale
				    che vogliamo determinare per mezzo
				    dell'analisi in frequenza */
  double omega[NFREQMAX];        /* omega[0], ..., omega[nfreqstudio-1] sono
				    le velocita' angolari corrispondenti ai
				    massimi relativi dell'ampiezza
				    dell'integrale fondamentale dell'analisi
				    in frequenza */
  double coord_base_non_ort[NFREQMAX][2];
                                 /* Componenti reali e immaginarie delle
				    coordinate del segnale rispetto ai
				    vettori della base non ortonormale */
  int nfreq_fond;                /* Numero delle frequenze fondamentali */
  int modkmax;                   /* Massimo modulo delle combinazioni
				    lineari a coefficienti interi che
				    devono approssimare le velocita'
				    angolari delle varie componenti
				    del segnale a partire da poche
				    frequenze fondamentali */
  double omega_fond[NFREQMAX];   /* omega_fond[0], ...,
				    omega_fond[nfreq_fond-1] sono le
				    velocita' angolari fondamentali da
				    cui si dovrebbero ottenere, per
				    combinazione lineare a
				    coefficienti interi, tutte le
				    altre che corrispondono ai massimi
				    relativi dell'ampiezza
				    dell'integrale fondamentale
				    dell'analisi in frequenza */
  int kbest[NFREQMAX];           /* Coefficienti interi della combinazione
				    lineare delle frequenze fondamentali che
				    costituisce la migliore approssimazione
				    della velocita' angolare che si sta
				    considerando */
  double err;                    /* Differenza tra la velocita' angolare che
				    si sta considerando e la combinazione
				    lineare delle frequenze fondamentali */
  double errmax;                 /* Massima differenza tra il segnale e la
				    sua ricostruzione effettuata secondo
				    il metodo dell'analisi in frequenza */
  FILE *in, *otp;                      /* Puntatore al file di input */
  char nomefile[30];             /* Stringa dove viene scritto il nome
				    del/dei file di input */
  /* Inizio delle istruzioni eseguibili */
  pig = 4. * atan(1.);
  dupig = 2. * pig;
  /* Apertura del file di dati concernenti alcune informazioni generali
     riguardanti l'analisi in frequenza che vogliamo effettuare */
  in = fopen("scomposizione.dat","r");
  fscanf(in,"%ld %d", &npt, &nfreqstudio);
  fscanf(in,"%le", &tstep);
  fscanf(in,"%d %d", &nfreq_fond, &modkmax);
  for(j=0; j<nfreq_fond; j++) 
  {
    fscanf(in,"%le", &omega_fond[j]);
  }
  fclose(in);  /* Chiusura del file dati */
  /* Definizione di tgrande.*/
  tgrande = .5 * npt * tstep;
  
    /* Apertura del file di dati contenente il
       segnale corrispondenti all'orbita che vogliamo studiare con
       l'analisi in frequenza */
    sprintf(nomefile, "pm_pendolo_eulero.out");
    in = fopen(nomefile,"r");
    printf("Analisi del file %s\n", nomefile);
    leggi_segnale(in, npt, 0);
    fclose(in);  /* Chiusura del file dati */
    for(j = 0; j<npt; j++)
    {
      segnale_orig[2*j] = segnale[2*j];
      segnale_orig[2*j + 1] = segnale[2*j + 1]; 
    }
    /* Calcolo di nfreqstudio velocita' angolari corrispondenti ai
       massimi relativi dell'ampiezza dell'integrale fondamentale
       dell'analisi in frequenza. Inoltre, vengono calcolate le
       coordinate del segnale espresse rispetto alla base NON
       ortonormale costituita dai vettori exp(i omega[0] t), ...,
       exp(i omega[nfreqstudio-1] t). */
    scomponi(TRUE, nfreqstudio, omega, coord_base_non_ort);
    /* Operazioni di stampa dei risultati della scomposizione del
       segnale */
    printf("Indice          omega         ");
    for(j=0; j<nfreq_fond-2; j++)  printf("  ");
    if(nfreq_fond > 0)  printf("comb. lin.");
    for(j=0; j<nfreq_fond-2; j++)  printf("  ");
    if(nfreq_fond > 0)  printf("    err    ");
    printf("   modulo      fase\n");
    for(j=0; j<nfreqstudio; j++) {
      int l;
      /* Calcolo della combinazione lineare a coefficienti interi
	 che a partire dalle frequenze fondamentali fornisce la
	 migliore approssimazione di omega[j] */
      mira_freq(nfreq_fond, modkmax, omega[j], omega_fond, kbest, &err);
      /* Stampa dei dati principali riguardanti la scomposizione
	 del segnale secondo il metodo dell'analisi in frequenza */
      printf("  %3d %24.17le ", j, omega[j]);
      for(l=0; l<nfreq_fond; l++)  printf(" %3d", kbest[l]);
      if(nfreq_fond > 0)  printf("  %10.4le ", err);
      printf("%10.4le %10.3le\n",
	     norma_in_C(&coord_base_non_ort[j][0]),
	     atan2(coord_base_non_ort[j][1], coord_base_non_ort[j][0]));
    }
    /* Riapertura del file di dati contenente il segnale che stiamo
       studiando con l'analisi in frequenza */
    in = fopen(nomefile,"r");
    /* Rilettura del segnale iniziale */
    leggi_segnale(in, npt, 0);
    fclose(in);  /* Chiusura del file dati */
    /* Sottrazione al segnale iniziale di tutte le sue componenti
       che sono state individuate con l'analisi in frequenza */
    prova_ainf(omega, nfreqstudio, coord_base_non_ort);
    /* Calcolo della massima differenza tra il segnale e la sua ricostruzione
       effettuata secondo il metodo dell'analisi in frequenza */
    errmax = 0.;
    otp = fopen("grafico_c.dat","w");
    for(i=0; i<npt; i++) 
    {
      double valtmp;
      valtmp = segnale[2*i] * segnale[2*i] + segnale[2*i+1] * segnale[2*i+1];
      if(errmax < valtmp) 
      {
	errmax = valtmp;
      }
      if(TRUE)
      {
         fprintf(otp,"%le  %le\n", segnale_orig[2*i] - segnale[2*i], segnale_orig[2*i + 1] - segnale[2*i + 1]);      
         //fprintf(otp,"%le  %le\n\n", segnale_orig[2*i +1] , segnale_orig[2*i]);
      }
    }
    fclose(otp);
    errmax = sqrt(errmax);
    /* Stampa della massima differenza */
    printf("\n\n  Ampiezza massima segnale dopo sottrazione finale:\n");
    printf("%24.17le\n\n\n", errmax);
}