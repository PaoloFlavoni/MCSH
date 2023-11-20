/* Libreria che include le routines necessarie a svolgere l'analisi in
   frequenza secondo le prescrizioni fornite da Laskar.
   Per un'introduzione all'analisi in frequenza si rimanda a:
     Laskar, J.: "Frequency map analysis of an Hamiltonian system", in
     "Workshop on Non--Linear Dynamics in Particle Accelarators,
     Arcidosso Sept. 1994", AIP Conf. Proc., Vol. 344, pag. 130-159
     (1995);
     Laskar, J.: "Introduction to frequency map analysis", in
     Proceedings of the NATO ASI school: "Hamiltonian Systems with
     Three or More Degrees of Freedom", S'Agaro (Spain), June 19-30,
     1995, C. Sim\`o (managing ed.), pag. 134-150, Kluwer (1999).

   Per poter utilizzare correttamente le routines seguenti, vi sono
   poche restrizioni fondamentali. Quella essenziale e' che l'oggetto
   in studio sia un segnale discretizzato ad istanti di tempo
   equidistanziati, quindi esso sara' descritto da una successione
   finita di valori nel piano complesso. Inoltre, quando si considera
   un segnale relativo al moto di un sistema dinamico, le routines
   funzionano in modo appropriato quando il sistema stesso e'
   autonomo; altrimenti sono necessarie alcune modifiche (non
   sostanziali, ma probabimente noiose, forse e' piu' semplice rendere
   autonomo il sistema) nei calcoli.

   Qui di seguito, e' descritta la struttura delle chiamate delle
   routines del programma di esempio nel file scomposizione.c (che
   ricostruisce le componenti periodiche principali di un segnale). Il
   suddetto programma utilizza tutte le routines qui incluse.

      main
        |
	|-- leggi_segnale
	|
	|-- scomponi
	|      |
	|      |-- max_segnale
	|      |        |
	|      |        |-- pesow
	|      |        |
	|      |        |-- dfour1
	|      |        |     |-- scambia 
	|      |        |
	|      |        |-- max_fft
	|      |        |
	|      |        |-- delimita_max
	|      |        |     |-- deriv_prodscal
	|      |        |     |          |-- pesow
	|      |        |     | 
	|      |        |     |-- scambia 
	|      |        |
	|      |        |-- newton
	|      |              |-- deriv_prodscal
	|      |                         |-- pesow
	|      |
	|      |-- ortonormalizza
	|      |         |-- prodscal_base_nonort
	|      |
	|      |-- proiezione_base_ort
	|      |            |-- prodscal
	|      |                   |-- pesow
	|      |
	|      |-- sottrai_al_segnale
	|      |
	|      |-- calcola_coord_base_non_ort
	|
	|-- mira_freq
	|
	|-- norma_in_C
        |
	|-- leggi_segnale  (seconda chiamata)
	|
	|-- prova_ainf

   Nel caso in cui siamo interessati al solo studio della frequenza
   principale (come puo' accadere quando vogliamo riprodurre la mappa
   dell'analisi in frequenza - si veda ad esempio il programma nel
   file mappa_freq.c - oppure quando consideriamo la variazione nel
   tempo della frequenza principale stessa - si veda ad esempio il
   programma in studia_var_freq.c -), allora la struttura delle
   chiamate sara' del tipo:

      main
        |
	|-- leggi_segnale
	|
	|-- max_segnale
	         |
	         |-- pesow
	         |
	         |-- dfour1
                 |     |-- scambia 
	         |
	         |-- max_fft
	         |
	         |-- delimita_max
                 |     |-- deriv_prodscal
                 |     |          |-- pesow
                 |     | 
                 |     |-- scambia 
                 |
                 |-- newton
                       |-- deriv_prodscal
                                  |-- pesow

    Quest'ultima struttura di chiamate e' quindi analoga a quella
    precedente a parte il fatto che si usano solo (e ripetutamente) le
    routines leggi_segnale, max_segnale, e ovviamente tutte le
    routines "piu' interne" di max_segnale.

    Le routines sono disposte in ordine alfabetico nel seguito del
    file.  */
#include <stdio.h>
#include <math.h>
/* Inclusione del file che contiene le definizioni generali */
#include "lib_ainf.h"
/* ### */
void calcola_coord_base_non_ort(mat_ort_nonort, n, coord_base_ort,
				coord_base_non_ort)
     double mat_ort_nonort[NFREQMAX][NFREQMAX];
     int n;
     double coord_base_ort[NFREQMAX][2], coord_base_non_ort[NFREQMAX][2];
     /* Routine che calcola le coordinate del segnale espresse
	rispetto alla base non ortonormale. Si ricordi che la base non
	ortonormale e' costituita dai vettori seguenti
	   exp(i vetfreq[0] t), ..., exp(i vetfreq[n-1] t),
	dove i e' l'unita' immaginaria, t e' il tempo e vetfreq e' il vettore
	di frequenze che viene passato alla routine ortonormalizza (e ciascuna
	delle sue componenti corrisponde a un massimo relativo dell'integrale
	fondamentale dell'analisi in frequenza).

	Variabili in entrata:
	  mat_ort_nonort[0][0], ...,  mat_ort_nonort[0][n-1]
                   .              .             .
                   .              .             .
	  mat_ort_nonort[n-1][0], ...,  mat_ort_nonort[n-1][n-1]
	                         => Matrice di passaggio dalle coordinate
				    relative alla base ortonormale a quelle
				    che si riferiscono alla base non ortonormale
	  n                      => Numero delle frequenze determinate finora
	                            tra quelle che corrispondono a dei massimi
				    relativi dell'integrale fondamentale
				    dell'analisi in frequenza
	  coord_base_ort[0][0] ,   coord_base_ort[0][1]
                  .                        .
                  .                        .
	  coord_base_ort[n-1][0] , coord_base_ort[n-1][1]
	                         => Componenti reali e immaginarie delle
				    proiezioni del segnale sui primi n
				    vettori della base ortonormale

	Variabili in uscita:
	  coord_base_non_ort[0][0] ,   coord_base_non_ort[0][1]
                    .                            .
                    .                            .
	  coord_base_non_ort[n-1][0] , coord_base_non_ort[n-1][1]
	                         => Componenti reali e immaginarie delle
				    coordinate del segnale rispetto agli
				    n vettori della base NON ortonormale
     */
{
  int i, j, l;
  /* printf("Entry calcola_coord_base_non_ort\n"); */
  for(i=0; i<n; i++) {
    for(l=0; l<2; l++) {
      coord_base_non_ort[i][l] = 0.;
    }
    for(j=i; j<n; j++) {
      for(l=0; l<2; l++) {
	coord_base_non_ort[i][l] += mat_ort_nonort[j][i] * coord_base_ort[j][l] ;
      }
    }
  }
  /* printf("Exit calcola_coord_base_non_ort\n"); */
}
/* ### */
void delimita_max(funz, flag_bisez, passo_x,
		  xa, funz_xa, d_funz_xa, d2_funz_xa, xb)
     double (*funz)();
     short unsigned int flag_bisez;
     double passo_x, *xa, *funz_xa, *d_funz_xa, *d2_funz_xa, *xb;
     /* Routine che determina due valori xa e xb (posti a distanza
	passo_x l'uno dall'altro) tale che la funzione funz (che
	compare tra gli argomenti della routine) soddisfa le due
	seguenti condizioni: (1) funz(xa) >= funz(xb) e (2) la
	derivata di funz in xa ha lo stesso segno della differenza
	xb - xa; quindi, la funzione funz ammette un massimo relativo
	nell'intervallo compreso tra xa e xb.

	Variabili in entrata:
	  funz                   => L'indirizzo corrispondente a funz rimanda
	                            a una routine che calcola il valore
				    della funzione cui siamo interessati
				    e i valori delle sue derivate prima e
				    seconda
	  flag_bisez             => Variabile flag che se e' uguale a FALSE (=0)
	                            allora si suppone che sia la prima volta
				    che questa routine viene chiamata, quindi
				    ci sono alcune operazioni in piu' che vanno
				    eseguite all'inizio e alla fine; se invece
				    e' uguale a TRUE (=1) si suppone che questa
				    non sia la prima chiamata della routine e
				    che essa venga utilizzata all'interno di
				    una procedura di ricerca del massimo
				    relativo con metodo di bisezione.
	  passo_x                => Valore dell'ampiezza dell'intervallo
	                            entro il quale vogliamo assicurarci che
				    esiste un massimo relativo
	  *xa                    => Approssimazione iniziale del valore
	                            del primo estremo xa

	Inoltre, se flag_bisez e' uguale a TRUE sono variabili in entrata
	anche le seguenti:
	  *funz_xa               => Valore assunto dalla funzione funz in
	                            corrispondenza all'estremo xa
	  *d_funz_xa             => Valore della derivata della funzione funz
	                            in corrispondenza all'estremo xa
	  *xb                    => Valore del secondo estremo xb.
	                            In verita', xb viene convenientemente
				    ridefinito all'inizio della routine senza
				    tener conto del valore che esso aveva in
				    entrata della routine stessa.
	                            Cio' nonostante, all'interno di una
				    procedura per bisezione al momento della
				    chiamata si suppone che:
				    (1) |xa - xb| = 2 * passo_x;
				    (2) funz(xa) >= funz(xb);
				    (3) d_funz_xa * (xb - xa) > 0;
				    ovvero, sappiamo che il massimo relativo
				    e' compreso tra xa e xb, ma in entrata
				    della routine l'ampiezza di questo intervallo
				    e' il doppio di quel che vorremmo.

	Variabili in uscita:
	  *xa                    => Valore del primo estremo xa dell'intervallo
	                            in cui sappiamo che esiste un massimo
				    relativo della funzione funz
	  *d_funz_xa             => Valore della derivata della funzione funz
	                            in corrispondenza all'estremo xa
	  *d2_funz_xa            => Valore della derivata seconda della funzione
	                            funz in corrispondenza all'estremo xa
	  *xb                    => Valore del secondo estremo xb dell'intervallo
	                            in cui sappiamo che esiste un massimo
				    relativo della funzione funz. Infatti,
				    in uscita xa, xb e d_funz_xa sono tali che:
				    (1) |xa - xb| = passo_x;
				    (2) funz(xa) >= funz(xb);
				    (3) d_funz_xa * (xb - xa) > 0;

	     *** Attenzione ***     Si suppone che funz sia una funzione
	                            per cui e' estremamente improbabile
				    che capiti di valutarla in un punto
				    in cui una delle derivate prima o
				    seconda si annulla. Infatti, e' facile
				    pensare un esempio in cui questa routine
				    sbaglia ad individuare un intervallo con
				    un massimo relativo, proprio nel caso in
				    cui le derivate esaminate dalla routine
				    non sono sempre diverse da zero.
				    Sarebbe possibile inserire dei tests in
				    modo da cautelarsi opportunamente da questa
				    eventualita', ma al fine di non complicare
				    e di non rallentare il programma, i suddetti
				    tests non sono stati inseriti.
     */
{
  short int flag_direzione;
  double funz_xb, d_funz_xb, d2_funz_xb, xb_old;
  /* printf("Entry delimita_max.\n"); */
  if(flag_bisez == FALSE) {
    /* Se questa e' la prima chiamata della routine, dobbiamo calcolare
       i valori (per ora ignoti) di funz(xa), d_funz_xa e d2_funz_xa */
    *funz_xa = (*funz)(*xa, d_funz_xa, d2_funz_xa);
  }
  else {
    xb_old = *xb;
  }
  if(*d_funz_xa >= 0.)  flag_direzione = 1;
  else  flag_direzione = -1;
  while(TRUE) {
    *xb = *xa + flag_direzione * passo_x;
    funz_xb = (*funz)(*xb, &d_funz_xb, &d2_funz_xb);
    if(*funz_xa >= funz_xb)  break;
    else {
      scambia(xa, xb);
      *d_funz_xa = d_funz_xb;
      *d2_funz_xa = d2_funz_xb;
      if(flag_direzione * (*d_funz_xa) < 0.)  break;
      if(flag_bisez == TRUE) {
	/* Abbiamo appena supposto che questa routine e' stata chiamata
           all'interno di una procedura di ricerca del massimo
           relativo secondo il metodo di bisezione. Di conseguenza, a
           questo punto siamo certi che il massimo relativo e'
           compreso tra il valore corrente di xa (cioe' il punto medio
           tra i valori che avevano gli estremi xa e xb al momento
           della chiamata della routine) e il valore iniziale
           dell'estremo xb. */
	*xb = xb_old;
	break;
      }
    }
  }
  /* printf("Exit delimita_max.\n"); */
}
/* ### */
double deriv_prodscal(nu, deriv1, deriv2)
     double nu, *deriv1, *deriv2;
     /* Routine che calcola il modulo al quadrato dell'integrale
	fondamentale dell'analisi in frequenza e le sue derivate prima
	e seconda in corrispondenza alla velocita' angolare nu.

	Variabili in entrata:
	  nu                     => Valore della velocita' angolare nu rispetto
	                            al quale calcoliamo le derivate dell'integrale
				    fondamentale dell'analisi in frequenza

	Variabili (di tipo globale) in entrata:
	  tgrande                => Tempo pari a meta' dell'intervallo di
	                            integrazione, cioe' l'integrale fondamentale
				    dell'analisi in frequenza viene calcolato
				    tra - tgrande e + tgrande
	  npt                    => Numero dei punti su cui calcoliamo
	                            approssimativamente l'integrale fondamentale
				    dell'analisi in frequenza
	  tstep                  => Tempo di scansione dei punti che ci consentono
	                            di calcolare approssimativamente l'integrale
				    fondamentale dell'analisi in frequenza;
				    tstep altro non e' che 2 * tgrande / npt
	  segnale[0], ..., segnale[2*npt+1]
	                         => Ognuna delle coppie di elementi consecutivi
				    rappresenta le componenti reali e
				    immaginarie del segnale che vogliamo
				    studiare con l'analisi in frequenza

	Variabili in uscita:
	                            Il valore restituito dalla routine
				    deriv_prodscal e' il quadrato del
				    modulo dell'integrale fondamentale
				    dell'analisi in frequenza.
	  *deriv1, *deriv2       => Valori delle derivate prima e seconda
	                            (eseguite rispetto alla velocita' angolare
				    nu) del modulo al quadrato dell'integrale
				    fondamentale dell'analisi in frequenza.

	     *** Attenzione ***     Il calcolo numerico dell'integrale
	                            fondamentale dell'analisi in frequenza
				    dipende dalla particolare funzione peso
				    che e' utilizzata nella routine pesow.
				    In particolare, SI ASSUME che la funzione
				    peso E' NULLA agli estremi dell'intervallo
				    integrazione; se cio' non e' verificato
				    occorre che la presente routine sia modificata
				    in modo tale che sia considerato il contributo
				    degli estremi al calcolo numerico
				    dell'integrale che viene qui effettuato
				    utilizzando la formula dei trapezi.  */
{
  int i, j;                      /* Indici usati come contatori nei cicli for */
  double t;                      /* Tempo utilizzato come variabile di
				    integrazione nel calcolo dell'integrale
				    fondamentale dell'analisi in frequenza */
  double reint, imint, redint1, imdint1, redint2, imdint2;
                                 /* Variabili che contengono rispettivamente
				    i valori reali e immaginari dell'integrale
				    fondamentale dell'analisi in frequenza e
				    delle sue derivate prime e seconde */
  double costmp, sintmp, pesotmp;
                                 /* Variabili che contengono i valori temporanei
				    calcolati rispettivamente a proposito delle
				    funzioni coseno, seno e peso (cioe' il filtro
				    dell'integrale fondamentale dell'analisi in
				    frequenza) */
  double retmpterm, imtmpterm;   /* Variabili che contengono alcuni valori
				    temporanei che servono a calcolare gli
				    integrali */
  /* printf("Entry deriv_prodscal\n"); */
  reint = imint = redint1 = imdint1 = redint2 = imdint2 = 0.;
  for(i=1; i<npt; i++) {
    /* Calcolo preliminare di alcune variabili (t, costmp, sintmp, pesotmp,
       retmpterm, imtmpterm) che vengono utili per la valutazione numerica
       degli integrali seguenti */
    t = i * tstep - tgrande ;
    costmp = cos(nu * t);
    sintmp = sin(nu * t);
    pesotmp = pesow(t);
    retmpterm = (segnale[2*i] * costmp + segnale[2*i+1] * sintmp) * pesotmp;
    imtmpterm = (- segnale[2*i] * sintmp + segnale[2*i+1] * costmp) * pesotmp;
    /* Aggiornamento della PARTE REALE dell'integrale
       fondamentale dell'analisi in frequenza */
    reint += retmpterm;
    /* Aggiornamento della PARTE IMMAGINARIA dell'integrale
       fondamentale dell'analisi in frequenza */
    imint += imtmpterm;
    /* Aggiornamento della PARTE REALE DELLA DERIVATA PRIMA
       dell'integrale fondamentale dell'analisi in frequenza */
    imtmpterm *= t;
    redint1 += imtmpterm;
    /* Aggiornamento della PARTE IMMAGINARIA DELLA DERIVATA PRIMA
       dell'integrale fondamentale dell'analisi in frequenza */
    retmpterm *= t;
    imdint1 -= retmpterm;
    /* Aggiornamento della PARTE REALE DELLA DERIVATA SECONDA
       dell'integrale fondamentale dell'analisi in frequenza */
    retmpterm *= t;
    redint2 -= retmpterm;
    /* Aggiornamento della PARTE IMMAGINARIA DELLA DERIVATA SECONDA
       dell'integrale fondamentale dell'analisi in frequenza */
    imtmpterm *= t;
    imdint2 -= imtmpterm;
  }
  /* Calcolo definitivo delle parti reali ed immaginarie dell'integrale
     fondamentale dell'analisi in frequenza e delle sue derivate prima
     e seconda */
  reint *= tstep / (2. * tgrande);  imint *= tstep / (2. * tgrande);
  redint1 *= tstep / (2. * tgrande);  imdint1 *= tstep / (2. * tgrande);
  redint2 *= tstep / (2. * tgrande);  imdint2 *= tstep / (2. * tgrande);
  /* Calcolo finale delle derivate prima e seconda del modulo al quadrato
     dell'integrale fondamentale dell'analisi in frequenza */
  *deriv1 = 2. * (reint * redint1 + imint * imdint1);
  *deriv2 = 2. * (reint * redint2 + redint1 * redint1 +
		 imint * imdint2 + imdint1 * imdint1);
  /* printf("Exit deriv_prodscal\n"); */
  return(reint * reint + imint * imint);
}
/* ### */
void dfour1(ncoppie, segno)
     unsigned long int ncoppie;
     int segno;
     /* Routine che esegue la Fast Fourier Transform secondo l'algoritmo
	di Danielson e Lanczos. Le istruzioni che seguono sono quasi
	esattamente la copia di quelle che si possono trovare nel paragrafo
	12.2 del Numerical Recipes, dove sono riportate molte altre
	informazioni concernenti la trasformata di Fourier.

	Variabili in entrata:
	  ncoppie                => Numero di coppie di dati cui vogliamo
	                            applicare la trasformata di Fourier
	     *** Attenzione ***     Questa routine funziona correttamente
	                            solo se ncoppie e' una potenza di 2,
				    ma non viene effettuato nessun controllo
				    sul valore di ncoppie
	     *** Consiglio ***      Se i dati significativi a disposizione
	                            sono tali che ncoppie non e' una potenza
				    di 2 si puo' mettere comunque la routine
				    in grado di funzionare (in modo rozzo)
				    riempiendo di zeri il vettore dati_fft finche'
				    ncoppie non diventa uguale a una potenza di 2.
	  segno                  => Se e' posto uguale a 1 viene eseguita la
	                            trasformata di Fourier diretta, se e'
				    invece uguale a -1 viene effettuata
				    quella inversa
	     *** Attenzione ***     Per trasformata di Fourier diretta qui
	                            s'intende quella che sul Numerical Recipes
				    e' considerata come inversa e viceversa

	Variabili (di tipo globale) in entrata:
	  dati_fft[0], ..., dati_fft[2*ncoppie-1]
	                         => Ognuna delle coppie di elementi consecutivi
				    rappresenta le componenti reali e
				    immaginarie del segnale (che si intende
				    campionato a intervalli di tempo regolari)
				    a cui vogliamo applicare la trasformata
				    di Fourier discreta

	Variabili (di tipo globale) in uscita:
	  dati_fft[0], ..., dati_fft[2*ncoppie-1]
	                         => Trasformata di Fourier discreta dei valori
				    che erano nel vettore dati_fft stesso
				    al momento della chiamata della routine

				    Supponiamo che il segnale in entrata sia
				    campionato ad intervalli di tempo D, cioe':
				      dati_fft[0],dati_fft[1] si riferiscono a t=0
				      dati_fft[2],dati_fft[3] si riferiscono a t=D
				         .       .         .             .    .
					 .       .         .             .    .
					 .       .         .             .    .
				      dati_fft[2*ncoppie-2],dati_fft[2*ncoppie-1]
				      si riferiscono a t = (ncoppie - 1) * D
				    Quando segno = 1, posto C = 2 * pigreco / D,
				    per convenzione in uscita si intende che
				     dati_fft[0],dati_fft[1] sono relativi alla
				     vel. ang. 0
				     dati_fft[2],dati_fft[3] sono relativi alla
				     vel. ang. C
				          .         .         .       .     .
					  .         .         .       .     .
				     dati_fft[ncoppie],dati_fft[ncoppie+1] sono
				     relativi
				     alla vel. ang. ncoppie / 2 * C
				     dati_fft[ncoppie+2],dati_fft[ncoppie+3] sono
				     relativi
				     alla vel. ang. - (ncoppie / 2 - 1) * C
				          .         .         .       .     .
					  .         .         .       .     .
				     dati_fft[2*ncoppie-2],dati_fft[2*ncoppie-1]
				     sono relativi alla vel. ang. - C.

				    Se invece, segno = - 1 allora i dati in
				    entrata e in uscita sono organizzati nel modo
				    opposto (cioe' a dati in entrata che sono
				    disposti in corrispondenza a multipli della
				    frequenza angolare C nell'ultimo modo indicato,
				    corrisponderanno in uscita dati che sono
				    relativi a multipli del tempo D, tra 0 e
				    (ncoppie - 1) * D).

	     *** Attenzione ***     Affinche' la routine funzioni correttamente
	                            deve essere definita la variabile globale
				    dupig che contiene il valore di 2 * pigreco
  */
{
  unsigned long int n, mmax, m, j, istep, i;
  double tmpr, tmpi, wtmp, wr, wpr, wpi, wi, theta;
  /* printf("Entry dfour1.\n"); */
  n = 2 * ncoppie;
  j = 0;
  for(i=0; i<n; i+=2) {
    if(j > i) {
      scambia(&dati_fft[j], &dati_fft[i]);
      scambia(&dati_fft[j+1], &dati_fft[i+1]);
    }
    m = n / 2;
    while(m >= 2 && j >= m) {
      j -= m;
      m /= 2;
    }
    j += m;
  }
  mmax = 2;
  while(n > mmax) {
    istep = 2 * mmax;
    theta = - segno * (dupig / mmax);
    wtmp = sin(0.5 * theta);
    wpr = -2. * wtmp * wtmp;
    wpi = sin(theta);
    wr = 1.;
    wi = 0.;
    for(m=0; m<mmax; m+=2) {
      for(i=m; i<n; i+=istep) {
	j = i + mmax;
	tmpr = wr * dati_fft[j] - wi * dati_fft[j+1];
	tmpi = wr * dati_fft[j+1] + wi * dati_fft[j];
	dati_fft[j] = dati_fft[i] - tmpr;
	dati_fft[j+1] = dati_fft[i+1] - tmpi;
	dati_fft[i] += tmpr;
	dati_fft[i+1] += tmpi;
      }
      wtmp = wr;
      wr = wr * wpr - wi * wpi + wr;
      wi = wi * wpr + wtmp * wpi + wi;
    }
    mmax = istep;
  }
  /* printf("Exit dfour1.\n"); */
  return;
}
/* ### */
int leggi_segnale(in, numfin, numletture)
     FILE *in;
     unsigned long int numfin;
     int numletture;
     /* Routine che legge (in modo binario) da file il segnale che vogliamo
	studiare con l'analisi in frequenza.

	Variabili in entrata:
	  in                     => Puntatore al file di input
	  numfin                 => Numero di dati (costituiti ciascuno da
	                            una coppia di valori, il primo dei quali
				    rappresenta la componente reale e il secondo
				    quella immaginaria) del segnale da leggere
	  numletture             => Numero di chiamate effettuate a questa routine
	                            al fine di leggere il segnale (o una parte di
				    esso)

	Variabile in uscita:
	                            Il valore restituito dalla routine
				    leggi_segnale e' un intero che, se e' uguale
				    a 0, significa che sono state lette
				    correttamente numfin + 1 coppie di dati che
				    costituiscono il nostro segnale da studiare;
				    se invece l'intero restituito dalla routine e'
				    maggiore di 0, allora non tutta l'operazione
				    di lettura ha funzionato correttamente
				    
	Variabili (di tipo globale) in uscita:
	  segnale[0], ..., segnale[2*numfin+1]
	                         => Ognuna delle coppie di elementi consecutivi
				    rappresenta le componenti reali e
				    immaginarie del segnale che e' stato letto
				    dal file e che vogliamo studiare con l'analisi
				    in frequenza
     */
{
  unsigned long int i;
  /* printf("Entry leggi_segnale.\n"); */
  if(numletture == 0) {
    fscanf(in,"%lf %lf", segnale, segnale+1);
  }
  else {
    /* Se non stiamo considerando la prima chiamata della presente routine,
       allora la prima coppia di dati e' uguale all'ultima coppia di dati
       letta durante la chiamata precedente */
    segnale[0] = segnale[2*numfin];
    segnale[1] = segnale[2*numfin+1];
  }
  // printf("%24.16le  %24.16le\n", segnale[0], segnale[1]);
  for(i=1; i<=numfin; i++) {
    fscanf(in,"%lf %lf", segnale+2*i, segnale+2*i+1);
    // printf("i = %d %24.16le  %24.16le\n", i, segnale[2*i], segnale[2*i+1]);
  }
  /* printf("Exit leggi_segnale.\n"); */
  return(0);
}
/* ### */
double max_fft(ncoppie)
     unsigned long int ncoppie;
     /* Routine che determina la velocita' angolare corrispondente alla
	coppia (costituita da un valore reale e uno immaginario) di massima
	ampiezza tra gli elementi del vettore dati_ftt.
	Affinche' la determinazione della velocita' angolare abbia senso,
	si deve intendere che dati_fft contiene la Fast Fourier Transform
	di un segnale che dipende dal tempo ed e' stato discretizzato
	ad intervalli di tempo tstep.
	Per avere maggiori informazioni su come viene memorizzata la
	Fast Fourier Transform nel vettore dati_fft, si faccia riferimento
	ai commenti relativi alla routine dfour1.

	Variabili in entrata:
	  ncoppie                => Numero di coppie di dati che costituiscono
	                            la trasformata di Fourier

	Variabili (di tipo globale) in entrata:
	  dati_fft[0], ..., dati_fft[2*ncoppie-1]
	                         => Ognuna delle coppie di elementi consecutivi
				    rappresenta le componenti reali e
				    immaginarie della trasformata di Fourier
				    discreta di un segnale, che si intende
				    campionato a intervalli di tempo regolari
	  tstep                  => Intervallo di tempo di scansione del
	                            segnale prima che venga effettuata la
				    trasformata di Fourier

	Variabile in uscita:
	                            Il valore restituito dalla routine max_fft
	                            e' la velocita' angolare che corrisponde alla
	                            coppia di dati di massima ampiezza (dove
				    per ampiezza intendiamo il modulo del numero
				    complesso che ha per parte reale il primo
				    valore della coppia di dati e per parte
				    immaginaria il secondo).
  */
{
  unsigned long int i, imax;
  double amp, maxamp, omegamax;
  /* printf("Entry max_fft.\n"); */
  /* Determinazione dell'ampiezza massima della trasformata */
  maxamp = 0.;
  for(i=0; i<ncoppie; i++) {
    amp = dati_fft[2*i] * dati_fft[2*i] + dati_fft[2*i+1] * dati_fft[2*i+1];
    if(maxamp < amp) {
      imax = i;
      maxamp = amp;
    }    
  }
  /* Calcolo della velocita' angolare che corrisponde all'ampiezza massima
     della trasformata di Fourier */
  // printf("imax = %d\n", imax);
  // printf("ncoppie = %d\n", ncoppie);
  if(imax < ncoppie / 2)  omegamax = (imax * dupig) / (ncoppie * tstep);
  else  omegamax = - ((ncoppie - imax) * dupig) / (ncoppie * tstep);
  /* printf("Exit max_fft.\n"); */
  return(omegamax);
}
/* ### */
double max_segnale(flag_prudente, err)
     unsigned short int flag_prudente;
     double err;
     /* Routine che individua una velocita' angolare che corrisponde a
	un massimo relativo dell'ampiezza dell'integrale fondamentale
	dell'analisi in frequenza.

	     *** Attenzione ***     Per come e' stata scritta questa routine,
	                            essa cerca di determinare la velocita'
				    angolare che corrisponde al MASSIMO ASSOLUTO;
				    il fatto che ci riesca o meno dipende da
				    quanto bene sono differenziate le quote
				    corrispondenti ai massimi dell'integrale
				    fondamentale dell'analisi in frequenza e
				    dipende anche da quanto e' buona
				    l'approssimazione fornita dalla Fast Fourier
				    transform.
				    Alla peggio, questa routine dovrebbe
				    determinare una velocita' angolare
				    corrispondente a uno dei principali massimi
				    relativi.

	Variabili in entrata:
	  flag_prudente          => Variabile flag che se e' uguale a TRUE (=1)
	                            imposta la ricerca dei massimi relativi in
				    modalita' prudente (che e' meno veloce, ma
				    piu' indicata quando ci sono segnali in cui
				    la componente caotica e'rilevante); se invece
				    e' uguale a FALSE (=0) imposta la ricerca dei
				    massimi relativi in modalita' veloce
	  err                    => Errore ammesso nella determinazione della
	                            velocita' angolare che corrisponde a
				    un massimo relativo dell'ampiezza
				    dell'integrale fondamentale dell'analisi in
				    frequenza

	Variabili (di tipo globale) in entrata:
	  tgrande                => Tempo pari a meta' dell'intervallo di
	                            integrazione, cioe' l'integrale fondamentale
				    dell'analisi in frequenza viene calcolato
				    tra - tgrande e + tgrande
	  npt                    => Numero dei punti su cui calcoliamo
	                            approssimativamente l'integrale fondamentale
				    dell'analisi in frequenza
	  tstep                  => Tempo di scansione dei punti che ci consentono
	                            di calcolare approssimativamente l'integrale
				    fondamentale dell'analisi in frequenza;
				    tstep altro non e' che 2 * tgrande / npt
	  segnale[0], ..., segnale[2*npt+1]
	                         => Ognuna delle coppie di elementi consecutivi
				    rappresenta le componenti reali e
				    immaginarie del segnale iniziale che
				    vogliamo studiare con l'analisi in
				    frequenza

	Variabile in uscita:
	                            Il valore restituito dalla routine max_segnale
				    e' la velocita' angolare che corrisponde
				    al massimo dell'ampiezza dell'integrale
				    fondamentale dell'analisi in frequenza */
{
  unsigned long int i;           /* Indice che funge da contatore all'interno
				    dei cicli for */ 
  unsigned long int ncoppie;     /* Numero di coppie di dati cui vogliamo
	                            applicare la trasformata di Fourier.
				    ncoppie viene posta uguale alla massima
				    potenza intera di 2 tra i numeri minori
				    o uguali di npt */
  double omega;                  /* Valore della velocita' angolare che corrisponde
				    a un massimo relativo dell'ampiezza
				    dell'integrale fondamentale dell'analisi in
				    frequenza */
  double omegaappr;              /* Valore approssimato della velocita' angolare
				    omega; omegaappr serve a determinare proprio
				    omega */
  double passo_x;                /* Valore dell'ampiezza dell'intervallo (di
	                            frequenze) entro il quale assicuriamo che
				    esiste un massimo relativo dell'integrale
				    fondamentale dell'analisi in frequenza */
  double xa;                     /* Valore del primo estremo dell'intervallo
	                            in cui esiste un massimo relativo
				    dell'integrale fondamentale dell'analisi
				    in frequenza */
  double funz_xa;                /* Valore assunto dal quadrato del modulo
				    dell'integrale fondamentale dell'analisi
				    in frequenza quando la velocita' angolare
				    e' uguale a xa */
  double d_funz_xa;              /* Valore assunto dalla derivata del quadrato
				    del modulo dell'integrale fondamentale
				    dell'analisi in frequenza quando la
				    velocita' angolare e' uguale a xa */
  double d2_funz_xa;             /* Valore assunto dalla derivata seconda
				    del quadrato del modulo dell'integrale
				    fondamentale dell'analisi in frequenza
				    quando la velocita' angolare e' uguale a xa */
  double xb;                     /* Valore del secondo estremo dell'intervallo
	                            in cui esiste un massimo relativo
				    dell'integrale fondamentale dell'analisi
				    in frequenza */
  double dx;                     /* Quando si cosidera un singolo passo del
				    metodo iterativo di Newton, dx e' il
				    valore della correzione della velocita'
				    angolare che corrisponde a un massimo
				    relativo dell'ampiezza dell'integrale
				    fondamentale dell'analisi in frequenza */
  /* printf("Entry max_segnale.\n"); */
  /* Definizione di ncoppie (ncoppie e' la massima potenza intera di 2
     tra quelle minori o uguali di npt). */
  ncoppie = 1;
  while(ncoppie <= npt / 2) {
    ncoppie *= 2;
  }
  /* Preparazione del vettore di dati_fft cui applicheremo la trasformata di
     Fourier.
     *** Osservazione *** : Gli esperimenti numerici, hanno evidenziato
     che la ricerca del punto di massimo assoluto e' piu' efficiente quando
     il segnale viene moltiplicato per la funzione peso prima dell'applicazione
     della Fast Fourier transform.
  */
  for(i=0; i<ncoppie; i++) {
    double t, pesotmp;
    t = i * tstep - tgrande ;
    pesotmp = pesow(t);
    dati_fft[2*i] = segnale[2*i] * pesotmp;
    dati_fft[2*i+1] = segnale[2*i+1] * pesotmp;
  }
  /* Calcolo della trasformata di Fourier */
  dfour1(ncoppie, 1);
  /* In prima approssimazione assumiamo la velocita' angolare che corrisponde
     alla massima ampiezza della trasformata di Fourier */
  omegaappr = max_fft(ncoppie);
  // printf("omegaappr = %24.16le\n", omegaappr);
  if(flag_prudente == TRUE) {
    /* Ricerca del massimo relativo in modalita' prudente:
       dapprima (1) viene determinato un intervallo in cui sicuramente esiste
       il massimo, poi (2) l'intervallo "delimitante" il massimo viene
       progressivamente diviso per due con un metodo di bisezione, finche'
       il metodo di Newton fornisce delle correzioni che sono sicuramente
       all'interno dell'intervallo "delimitante", allora (3) viene innescata
       la ricerca del massimo con il metodo di Newton */
    xa = omegaappr;
    /* Un valore pari a 1/5 della distanza di scansione dei punti della
       trasformata di Fourier (cioe' pig / (tgrande * 2) ) e' una buona
       ampiezza dell'intervallo nel quale vogliamo assicurare l'esistenza
       di un massimo.
       Per comprendere il motivo di questa scelta si osservi l'algoritmo
       utilizzato in delimita_max. */
    passo_x = pig / (tgrande * 10.);
    /* Determiniamo un primo intervallo in cui sicuramente esiste il
       massimo */
    delimita_max(deriv_prodscal, FALSE, passo_x,
		 &xa, &funz_xa, &d_funz_xa, &d2_funz_xa, &xb);
    /* Cominciamo la ricerca del massimo relativo con il metodo di
       bisezione */
    while(TRUE) {
      passo_x *= .5;
      dx = d_funz_xa / d2_funz_xa;
      if(fabs(dx) < passo_x && d2_funz_xa < 0.) {
	/* Il metodo di Newton e' ora competitivo con quello di bisezione */
	omegaappr = xa - dx;
	/* Calcolo definitivo di omega secondo il metodo di Newton */
	omega = newton(deriv_prodscal, omegaappr, err);
	break;
      }
      /* Determiniamo ora un altro intervallo in cui esiste sicuramente il
	 massimo relativo e la lunghezza del nuovo intervallo deve essere
	 la meta' di quello precedente */
      delimita_max(deriv_prodscal, TRUE, passo_x,
		   &xa, &funz_xa, &d_funz_xa, &d2_funz_xa, &xb);
      if(passo_x < err) {
	/* Il metodo di bisezione e' riuscito a delimitare il massimo
	   in un intervallo la cui ampiezza e' inferiore all'errore
	   ammesso */
	omega = .5 * (xa + xb);
	break;
      }
    }
  }
  else {
    /* Ricerca del massimo relativo in modalita' veloce */
    /* Temporaneamente consideriamo un intervallo di integrazione
       rispetto al tempo di lunghezza ridotta; quindi ridefiniamo npt
       e tgrande */
    npt /= 4;
    tgrande /= 4.;
    /* In seconda approssimazione assumiamo la velocita' angolare che
       corrisponde alla massima ampiezza dell'integrale fondamentale
       dell'analisi in frequenza per l'intervallo di integrazione di
       lunghezza ridotta */
    omegaappr = newton(deriv_prodscal, omegaappr, err);
    /* Ritorniamo all'intervallo di integrazione di lunghezza originaria,
       quindi ridefiniamo npt e tgrande */
    npt *= 4;
    tgrande *= 4.;
    /* Calcolo definitivo di omega con il metodo di Newton */
    omega = newton(deriv_prodscal, omegaappr, err);
  }
  /* printf("Exit max_segnale.\n"); */
  return(omega);
}
/* ### */
void mira_freq(nfreq_fond, modkmax, omega, omega_fond, kbest, errmin)
     int nfreq_fond, modkmax;
     double omega, omega_fond[];
     int kbest[];
     double *errmin;
     /* Routine che calcola la combinazione lineare a coefficienti
	interi che a partire dalle frequenze fondamentali fornisce
	la migliore approssimazione di omega.

	Variabili in entrata:
	  nfreq_fond             => Numero delle frequenze fondamentali
	  modkmax                => Massimo modulo delle combinazioni lineari
				    a coefficienti interi che devono approssimare
				    omega a partire dalle frequenze fondamentali
	  omega                  => Valore della velocita' angolare di cui
	                            vogliamo determinare la migliore
				    approssimazione con una combinazione
				    lineare a coefficienti interi delle
				    frequenze fondamentali
	  omega_fond[0], ..., omega_fond[nfreq_fond-1]
	                         => Valori delle frequenze fondamentali

	Variabili in uscita:
	  kbest[0], ..., kbest[nfreq_fond-1]
	                         => Coefficienti interi della combinazione
				    lineare delle frequenze fondamentali che
				    costituisce la migliore approssimazione
				    di omega
	  *errmin                => Valore della differenza tra omega e la
	                            combinazione lineare delle frequenze
				    fondamentali che costituisce la migliore
				    approssimazione della stessa omega
     */
{
  int i, modk, indbox, box[NFREQMAX], k[NFREQMAX];
  double err;
  /* printf("Entry mira_freq.\n"); */
  if(nfreq_fond == 0) {
    return;
  }
  /* Definizioni iniziali di kbest e di *errmin */
  for(i=0; i<nfreq_fond; i++) {
    kbest[i] = 0;
  }
  *errmin = fabs(omega);
  /* Ciclo su modk: si effettuano i confronti tra omega e le combinazioni
     lineari a coefficienti interi k[0], ..., k[nfreq_fond-1] delle
     frequenze fondamentali. k[0], ..., k[nfreq_fond-1] sono tali che
     la |k[0]| + ... |k[nfreq_fond-1]| = modk . */
  for(modk=1; modk<=modkmax; modk++) {
    /* Definizione iniziale del vettore di coefficienti interi k e di
       alcune quantita' utili per i successivi aggiornamenti di k */
    for(i=0; i<nfreq_fond; i++) {
      k[i] = 0;
    }
    k[0] = modk;
    indbox = 0;
    box[0] = 0;
    while(indbox >= 0) {
      /* Calcolo della differenza tra omega e la combinazione lineare a
	 coefficienti interi k[0], ..., k[nfreq_fond-1] delle frequenze
	 fondamentali */
      err = omega;
      for(i=0; i<nfreq_fond; i++) {
	err -= k[i] * omega_fond[i];
      }
      err = fabs(err);
      if(err < *errmin) {
	/* Aggiornamento della migliore approssimazione */
	*errmin=err;
	for(i=0; i<nfreq_fond; i++) {
	  kbest[i]=k[i];
	}
      }
      /* Aggiornamento del vettore a coefficienti interi k */
      k[box[indbox]] *= -1;
      if(k[box[indbox]] > 0) {
	if(box[indbox] == nfreq_fond-1)  {
	  indbox--;
	  k[nfreq_fond-1]=0;
	  if(indbox < 0)  break;
	  k[box[indbox]] *= -1;
	  if(k[box[indbox]] > 0) {
	    k[box[indbox]]--;
	  }
	}
	else {
	  k[box[indbox]]--;
	}
	indbox++;
	box[indbox] = box[indbox-1] + 1;
	k[box[indbox]] = modk;
	for(i=0; i<box[indbox]; i++) {
	  if(k[i] >= 0) {
	    k[box[indbox]] -= k[i];
	  }
	  else {
	    k[box[indbox]] += k[i];
	  }
	}
	if(k[box[indbox-1]] == 0) {
	  box[indbox-1] = box[indbox];
	  indbox--;
	}
      }
      /* Fine aggiornamento del vettore a coefficienti interi k */
    }
  }
  /* printf("Exit mira_freq.\n"); */
}
/* ### */
double newton(funz, xiniz, err)
     double (*funz)();
     double xiniz, err;
     /* Routine che risolve un'equazione del tipo f(x)=0 con il metodo di
	Newton.

	Variabili in entrata:
	  funz                   => L'indirizzo corrispondente a funz rimanda
	                            a una routine che calcola il valore
				    della funzione f(x) cui siamo interessati
				    e della sua derivata in corrispondenza
				    al valore di x che gli viene passato come
				    argomento della funzione stessa
	  xiniz                  => Approssimazione iniziale della soluzione
	  err                    => Errore massimo consentito nella determinazione
	                            della soluzione stessa

	Variabile in uscita:
	                            Il valore restituito dalla routine newton
				    corrisponde all'approssimazione numerica
				    della soluzione dell'equazione f(x)=0
     */
{
  int j;
  double x, f, dx, df;
  /* printf("Entry newton.\n"); */
  x = xiniz;
  for(j=0; j<MAXITERNEWTON; j++) {
    (*funz)(x, &f, &df);
    dx = f / df;
    x -= dx;
    if(fabs(dx) < err) {
      /* printf("Exit newton.\n"); */
      return(x);
    }
  }
  printf("  Il metodo di Newton non e' convergente!\n");
  /* printf("Exit newton.\n"); */
  return(x);
}
/* ### */
double norma_in_C(z)
     double z[2];
     /* Routine che calcola la norma del numero complesso avente z[0]
	come parte reale e z[1] come parte immaginaria */
{
  double tmp;
  /* printf("Entry norma_in_C.\n"); */
  tmp = sqrt(z[0] * z[0] + z[1] * z[1]);
  /* printf("Exit norma_in_C.\n"); */
  return(tmp);
}
/* ### */
void ortonormalizza(mat_nonort_ort, mat_ort_nonort, n, vetfreq)
     double mat_nonort_ort[NFREQMAX][NFREQMAX], mat_ort_nonort[NFREQMAX][NFREQMAX];
     int n;
     double vetfreq[NFREQMAX];
     /* Routine che calcola una riga (la n-esima) della matrice che consente
	il passaggio dalle coordinate relative alla base non ortonormale
	costituita dai vettori
	   exp(i vetfreq[0] t), ..., exp(i vetfreq[n-1] t),
	(dove i e' l'unita' immaginaria e t il tempo) alle coordinate che
	si riferiscono alla base ortonormale che si ottiene dalla precedente
	seguendo l'algoritmo di Gram-Schmidt.
	Inoltre, viene contemporaneamente calcolata anche la n-esima riga
	della matrice inversa.

	Variabili in entrata:
	  mat_nonort_ort[0][0], ...,  mat_nonort_ort[0][n-2]
                   .              .             .
                   .              .             .
	  mat_nonort_ort[n-2][0], ...,  mat_nonort_ort[n-2][n-2]
	                         => Matrice di passaggio dalle coordinate
				    relative alla base non ortonormale
				    a quelle che si riferiscono alla base
				    ortonormale
	  mat_ort_nonort[0][0], ...,  mat_ort_nonort[0][n-2]
                   .              .             .
                   .              .             .
	  mat_ort_nonort[n-2][0], ...,  mat_ort_nonort[n-2][n-2]
	                         => Matrice di passaggio dalle coordinate
				    relative alla base ortonormale a quelle
				    che si riferiscono alla base non ortonormale
	  n                      => Indice della riga della matrice che dobbiamo
	                            calcolare
	  vetfreq[0], ...,  vetfreq[n-1]
	                         => Valori delle frequenze che corrispondono
				    a dei massimi relativi dell'integrale
				    fondamentale dell'analisi in frequenza

	Variabili in uscita:
	  mat_nonort_ort[0][0], ...,  mat_nonort_ort[0][n-1]
                   .              .             .
                   .              .             .
	  mat_nonort_ort[n-1][0], ...,  mat_nonort_ort[n-1][n-1]
	                         => Rispetto ai valori in entrata, viene
				    aggiunta la n-esima riga. Si osservi che
				    l'algoritmo di Gram-Schmidt prevede che
				    tutti i termini della n-esima colonna
				    sopra la diagonale sono posti uguali a 0,
				    quindi non occorre calcolarli.
	  mat_ort_nonort[0][0], ...,  mat_ort_nonort[0][n-1]
                   .              .             .
                   .              .             .
	  mat_ort_nonort[n-1][0], ...,  mat_ort_nonort[n-1][n-1]
	                         => Rispetto ai valori in entrata, viene
				    aggiunta la n-esima riga. Si osservi che
				    l'algoritmo di Gram-Schmidt prevede che
				    tutti i termini della n-esima colonna
				    sopra la diagonale sono posti uguali a 0,
				    quindi non occorre calcolarli.
     */
{
  int i, j;
  double norma;
  /* printf("Entry ortonormalizza\n"); */
  /* Allo scopo di semplificare la notazione nei commenti di questa
     routine, d'ora in poi sia u_0, ..., u_n-1 la base non ortonormale
     dei vettori e sia e_0, ..., e_n-1 quella ortonormale */
  /* Ricordiamo che e_n-1 e' definito da e_n-1 = v_n-1 / || v_n-1 || , dove
       v_n-1 = u_n-1 - (u_n-1 , e_n-2) e_n-2 - ... - (u_n-1 , e_0) e_0 */
  /* Cominciamo a scrivere la formula che e' immediata conseguenza della
     precedente, ovvero
      u_n-1 - || v_n-1 || e_n-1 = (u_n-1 , e_n-2) e_n-2 + ... + (u_n-1 , e_0) e_0 ,
     nella n-esima riga di mat_nonort_ort */
  for(i=0; i<n-1; i++) {
    /* Calcolo di (u_n-1 , e_i) .
       Si ricordi che dobbiamo ricondurci a dei prodotti scalari con vettori
       della base non ortonormale, cioe' del tipo (u_n-1 , u_j), perche' sono
       gli unici che si possono calcolare banalmente */
    mat_nonort_ort[n-1][i] = 0.;
    for(j=0; j<=i; j++) {
      mat_nonort_ort[n-1][i] +=
	mat_ort_nonort[i][j] * prodscal_base_nonort(vetfreq[n-1], vetfreq[j]);
    }
  }
  /* Calcolo di || v_n-1 || .
     Si ricordi che || u_n-1 || = 1 e siccome 
      u_n-1 = || v_n-1 || e_n-1 + (u_n-1 , e_n-2) e_n-2 + ... + (u_n-1 , e_0) e_0 ,
     possiamo sfruttare il teorema di Pitagora in n dimensioni */
  norma = 1.;
  for(i=0; i<n-1; i++) {
    norma -= mat_nonort_ort[n-1][i] * mat_nonort_ort[n-1][i];
  }
  norma = sqrt(norma);
  /* Completiamo la scrittura della formula
     u_n-1 = || v_n-1 || e_n-1 + (u_n-1 , e_n-2) e_n-2 + ... + (u_n-1 , e_0) e_0
     nella n-esima riga di mat_nonort_ort */
  mat_nonort_ort[n-1][n-1] = norma;
  /* Scriviamo ora
       e_n-1 = a_n-1,n-1 u_n-1 + a_n-1,n-2 u_n-2 + ... + a_n-1,0 u_0 ,
     dove i coefficienti a_n-1,i altro non sono che i valori finali di
     mat_ort_nonort[n-1][i]. A questo scopo, dobbiamo sostituire nella
     formula
     e_n-1 = [u_n-1 - (u_n-1 , e_n-2) e_n-2 - ... - (u_n-1 , e_0) e_0]/|| v_n-1 ||
     le espressioni lineari di e_i in funzione della base non ortonormale, cioe'
       e_i = a_i,i u_i + a_i,i-1 u_i-1 + ... + a_i,0 u_0 ,
     dove i=0, ..., n-2. */
  mat_ort_nonort[n-1][n-1] = 1. / norma;
  for(i=0; i<n-1; i++) {
    mat_ort_nonort[n-1][i] = 0.;
  }
  for(j=0; j<n-1; j++) {
    for(i=j; i<n-1; i++) {
      mat_ort_nonort[n-1][j] += mat_nonort_ort[n-1][i] * mat_ort_nonort[i][j];
    }
  }
  for(j=0; j<n-1; j++) {
    mat_ort_nonort[n-1][j] /= (-norma);
  }
  /* printf("Exit ortonormalizza\n"); */
}
/* ### */
double pesow(t)
     double t;
     /* Routine che calcola la funzione peso che compare nell'integrale
	fondamentale dell'analisi in frequenza.
	Nella versione attuale il peso viene calcolato utilizzando il
	"filtro di Hanning".
	     *** Consiglio ***      Se vogliamo utilizzare il metodo dell'analisi
	                            in frequenza senza utilizzare la funzione
				    peso, basta scrivere tmp = 1. in luogo
				    dell'istruzione che esegue il calcolo della
				    variabile tmp

	     *** Attenzione ***     Se vogliamo variare la presente routine
	                            pesow, allora dobbiamo modificare
				    coerentemente anche il calcolo nella
				    routine prodscal_base_nonort che dipende
	                            dal peso che si e' assunto nella definizione
				    dell'integrale fondamentale dell'analisi
				    in frequenza.
	Variabile in entrata:
	  t                      => Valore del tempo che compare come variabile
	                            di integrazione nel calcolo dell'integrale
				    fondamentale dell'analisi in frequenza.

	Variabile (di tipo globale) in entrata:
	  tgrande                => Tempo pari a meta' dell'intervallo di
	                            integrazione
     */
{
  double tmp;
  tmp = 1. + cos(pig * t / tgrande);
  return(tmp);
}
/* ### */
void prodscal(omega, valint)
     double omega, valint[2];
     /* Routine che calcola solamente le componenti reale e
	immaginaria dell'integrale fondamentale dell'analisi in
	frequenza in corrispondenza alla velocita' angolare omega.
	Questa routine altro non e' che una versione accorciata di
	deriv_prodscal.

	Variabili in entrata:
	  omega                  => Valore della velocita' angolare omega per
	                            cui calcoliamo la parte reale e immaginaria
				    dell'integrale fondamentale dell'analisi
				    in frequenza

	Variabili (di tipo globale) in entrata:
	  tgrande                => Tempo pari a meta' dell'intervallo di
	                            integrazione, cioe' l'integrale fondamentale
				    dell'analisi in frequenza viene calcolato
				    tra - tgrande e + tgrande
	  npt                    => Numero dei punti su cui calcoliamo
	                            approssimativamente l'integrale fondamentale
				    dell'analisi in frequenza
	  tstep                  => Tempo di scansione dei punti che ci consentono
	                            di calcolare approssimativamente l'integrale
				    fondamentale dell'analisi in frequenza;
				    tstep altro non e' che 2 * tgrande / npt
	  segnale[0], ..., segnale[2*npt+1]
	                         => Ognuna delle coppie di elementi consecutivi
				    rappresenta le componenti reali e
				    immaginarie del segnale che vogliamo
				    studiare con l'analisi in frequenza

	Variabili in uscita:
	  valint[0], valint[1]   => Valori reale e immaginario dell'integrale
				    fondamentale dell'analisi in frequenza

	     *** Attenzione ***     Il calcolo numerico dell'integrale
	                            fondamentale dell'analisi in frequenza
				    dipende dalla particolare funzione peso
				    che e' utilizzata nella routine pesow.
				    In particolare, SI ASSUME che la funzione
				    peso E' NULLA agli estremi dell'intervallo
				    integrazione; se cio' non e' verificato
				    occorre che la presente routine sia modificata
				    in modo tale che sia considerato il contributo
				    degli estremi al calcolo numerico
				    dell'integrale che viene qui effettuato
				    utilizzando la formula dei trapezi.
  */
{
  int i, j;                      /* Indici usati come contatori nei cicli for */
  double t;                      /* Tempo utilizzato come variabile di
				    integrazione nel calcolo dell'integrale
				    fondamentale dell'analisi in frequenza */
  double costmp, sintmp, pesotmp;
                                 /* Variabili che contengono i valori temporanei
				    calcolati rispettivamente a proposito delle
				    funzioni coseno, seno e peso (cioe' il filtro
				    dell'integrale fondamentale dell'analisi in
				    frequenza) */
  double retmpterm, imtmpterm;   /* Variabili che contengono alcuni valori
				    temporanei che servono a calcolare gli
				    integrali */
  /* printf("Entry prodscal\n"); */
  valint[0] = valint[1] = 0.;
  for(i=1; i<npt; i++) {
    /* Calcolo preliminare di alcune variabili (t, costmp, sintmp, pesotmp,
       retmpterm, imtmpterm) che vengono utili per la valutazione numerica
       degli integrali seguenti */
    t = i * tstep - tgrande ;
    costmp = cos(omega * t);
    sintmp = sin(omega * t);
    pesotmp = pesow(t);
    retmpterm = (segnale[2*i] * costmp + segnale[2*i+1] * sintmp) * pesotmp;
    imtmpterm = (- segnale[2*i] * sintmp + segnale[2*i+1] * costmp) * pesotmp;
    /* Aggiornamento della PARTE REALE dell'integrale
       fondamentale dell'analisi in frequenza */
    valint[0] += retmpterm;
    /* Aggiornamento della PARTE IMMAGINARIA dell'integrale
       fondamentale dell'analisi in frequenza */
    valint[1] += imtmpterm;
  }
  /* Calcolo definitivo delle parti reali ed immaginarie dell'integrale
     fondamentale dell'analisi in frequenza e delle sue derivate prima
     e seconda */
  valint[0] *= tstep / (2. * tgrande);  valint[1] *= tstep / (2. * tgrande);
  /* printf("Exit prodscal\n"); */
}
/* ### */
double prodscal_base_nonort(omega1, omega2)
     double omega1, omega2;
     /* Routine che calcola il prodotto interno (relativo all'integrale
	fondamentale dell'analisi in frequenza) tra le funzioni
	        exp(i omega1 t)   e   exp(i omega2 t)

	Variabili in entrata:
	  omega1, omega2         => Velocita' angolari delle funzioni di
	                            cui vogliamo fare il prodotto scalare

	Variabili (di tipo globale) in entrata:
	  tgrande                => Tempo pari a meta' dell'intervallo di
	                            integrazione, cioe' l'integrale fondamentale
				    dell'analisi in frequenza viene calcolato
				    tra - tgrande e + tgrande

	Variabile in uscita:
	                            Il valore restituito dalla routine
				    prodscal_base_nonort e' proprio il prodotto
				    interno tra le funzioni exp(i omega1 t) e
				    exp(i omega2 t)

	     *** Attenzione ***     Il calcolo del prodotto interno dipende
	                            dal peso che si e' assunto nella definizione
				    dell'integrale fondamentale dell'analisi
				    in frequenza; quindi la presente routine
				    prodscal_base_nonort deve essere variata
				    quando viene cambiata pesow.
				    Nella versione attuale di questa routine,
				    il calcolo si riferisce al "filtro di
				    Hanning".
     */
{
  double tmp;
  /* printf("Entry prodscal_base_nonort\n"); */
  if(omega1 == omega2) {
    tmp = 1.;
  }
  else {
    tmp = (omega1 - omega2) * tgrande;
    tmp = sin(tmp) / tmp + .5 * sin(tmp + pig) / (tmp + pig)
      + .5 * sin(tmp - pig) / (tmp - pig);
  }
  /* printf("Exit prodscal_base_nonort\n"); */
  return(tmp);
}
/* ### */
void proiezione_base_ort(omega, mat_nonort_ort_n, coord_base_ort_n)
     double omega;
     double mat_nonort_ort_n;
     double coord_base_ort_n[2];
     /* Routine che calcola la proiezione del segnale sull'n-esimo vettore
	della base ortonormale. Si ricordi che la base ortonormale e'
	quella che si ottiene secondo l'algoritmo di Gram-Schmidt dai
	vettori seguenti
	   exp(i vetfreq[0] t), ..., exp(i vetfreq[n-1] t),
	dove i e' l'unita' immaginaria, t e' il tempo e vetfreq e' il vettore
	di frequenze che viene passato alla routine ortonormalizza (e ciascuna
	delle sue componenti corrisponde a un massimo relativo dell'integrale
	fondamentale dell'analisi in frequenza).

	Variabili in entrata:
	  omega                  => n-esima frequenza trovata che corrisponde
	                            a un massimo relativo dell'integrale
				    fondamentale dell'analisi in frequenza
	  mat_nonort_ort_n       => n-esimo elemento sulla diagonale della
	                            matrice di passaggio dalle coordinate
				    relative alla base non ortonormale
				    a quelle che si riferiscono alla base
				    ortonormale

	Variabili in uscita:
	  coord_base_ort_n[0] , coord_base_ort_n[1]
	                         => Componenti reale e immaginaria della
				    proiezione del segnale sull'n-esimo
				    vettore della base ortonormale
     */
{
  int j;
  double valint[2];
  /* printf("Entry proiezione_base_ort\n"); */
  /* Allo scopo di semplificare la notazione nei commenti di questa
     routine, d'ora in poi sia u_0, ..., u_n-1 la base non ortonormale
     dei vettori e sia e_0, ..., e_n-1 quella ortonormale. Inoltre,
     sia f=f(t) il segnale di cui vogliamo calcolare la proiezione
     sull'n-esimo vettore della base ortonormale, cioe' (f , e_n-1).
        *** Attenzione ***   Si intende che quando viene chiamata questa
	                     routine al segnale sono gia' state sottratte
			     le componenti appartenenti allo spazio lineare
			     generato da e_0, ..., e_n-2, cioe'
			     (f , e_0) = ... = (f , e_n-2) = 0 */
  /* Calcolo del prodotto interno (f , u_n-1) */
  prodscal(omega, valint);
  /* Calcoliamo ora il prodotto interno (f , e_n-1).
     A questo scopo, utilizzando la formula 
       u_n-1 = b_n-1,n-1 e_n-1 + b_n-1,n-2 e_n-2 + ... + b_n-1,0 e_0 ,
     possiamo scrivere il prodotto interno (f , u_n-1) come segue:
       (f , u_n-1) = b_n-1,n-1 (f , e_n-1) + ... + b_n-1,0  (f , e_0) ,
     ma siccome abbiamo assunto l'ipotesi che il segnale f e' tale che
       (f , e_0) = ... = (f , e_n-2) = 0 ,
     allora ne segue che
       (f , e_n-1) = (f , u_n-1) / b_n-1,n-1 .
     Ricordando che (A) i valori delle parti reale e immaginaria di (f , e_n-1)
     sono memorizzati rispettivamente nelle variabili coord_base_ort_n[0]
     e coord_base_ort_n[1] e che (B) il valore di b_n-1,n-1 e' memorizzato
     in mat_nonort_ort_n, dovrebbe essere evidente che le seguenti
     istruzioni permettono di calcolare (f , e_n-1). */
  for(j=0; j<2; j++) {
    coord_base_ort_n[j] = valint[j] / mat_nonort_ort_n;
  }
  /* printf("Exit proiezione_base_ort\n"); */
}
/* ### */
void prova_ainf(vetfreq, n, coord_base_non_ort)
     double vetfreq[NFREQMAX];
     int n;
     double coord_base_non_ort[NFREQMAX][2];
     /* Routine che sottrae al segnale iniziale la sua stessa proiezione
	sullo spazio formato dai primi n vettori della base ortonormale.
	Il calcolo (per ragioni di efficienza e perche' e' un test piu'
	significativo) e' effettuato utilizzando le componenti del segnale
	rispetto ai vettori della base NON ortonormale.
	Si ricordi che la base non ortonormale e' costituita dai vettori
	seguenti
	   exp(i vetfreq[0] t), ..., exp(i vetfreq[n-1] t),
	dove i e' l'unita' immaginaria e t e' il tempo.

	Variabili in entrata:
	  vetfreq[0], ...,  vetfreq[n-1]
	                         => Valori delle frequenze che corrispondono
				    a dei massimi relativi dell'integrale
				    fondamentale dell'analisi in frequenza
	  n	          	 => Numero delle componenti del segnale
				    che sono state individuate seguendo
				    il metodo dell'analisi in frequenza
	  coord_base_non_ort[0][0] ,   coord_base_non_ort[0][1]
                    .                            .
                    .                            .
	  coord_base_non_ort[n-1][0] , coord_base_non_ort[n-1][1]
	                         => Componenti reali e immaginarie delle
				    coordinate del segnale rispetto agli
				    n vettori della base NON ortonormale

	Variabili (di tipo globale) in entrata:
	  tgrande                => Tempo pari a meta' dell'intervallo di
	                            integrazione, cioe' l'integrale fondamentale
				    dell'analisi in frequenza viene calcolato
				    tra - tgrande e + tgrande
	  npt                    => Numero dei punti su cui calcoliamo
	                            approssimativamente l'integrale fondamentale
				    dell'analisi in frequenza
	  tstep                  => Tempo di scansione dei punti che ci consentono
	                            di calcolare approssimativamente l'integrale
				    fondamentale dell'analisi in frequenza;
				    tstep altro non e' che 2 * tgrande / npt
	  segnale[0], ..., segnale[2*npt+1]
	                         => Ognuna delle coppie di elementi consecutivi
				    rappresenta le componenti reali e
				    immaginarie del segnale iniziale che
				    vogliamo studiare con l'analisi in
				    frequenza
				    
	Variabili (di tipo globale) in uscita:
	  segnale[0], ..., segnale[2*npt+1]
	                         => Ognuna delle coppie di elementi consecutivi
				    rappresenta le componenti reali e
				    immaginarie del segnale che abbiamo ricevuto
				    in entrata e a cui sono state sottratte tutte
				    le sue n componenti che abbiamo individuato
				    con l'analisi in frequenza
     */
{
  unsigned long int i;
  int j, l;
  double t, costmp, sintmp, term_sottr[2];
  /* printf("Entry prova_ainf.\n"); */
  for(j=0; j<n; j++) {
    for(i=0; i<=npt; i++) {
      /* Calcolo del tempo t, di cos(vetfreq[j] t) e sin(vetfreq[j] t) */
      t = i * tstep - tgrande ;
      costmp = cos(vetfreq[j] * t);
      sintmp = sin(vetfreq[j] * t);
      term_sottr[0] = coord_base_non_ort[j][0] * costmp -
	coord_base_non_ort[j][1] * sintmp;
      term_sottr[1] = coord_base_non_ort[j][0] * sintmp +
	coord_base_non_ort[j][1] * costmp;
      for(l=0; l<2; l++) {
	segnale[2*i+l] -= term_sottr[l];
      }
    }
  }
  /* printf("Exit prova_ainf.\n"); */
}
/* ### */
void scambia(x,y)
     double *x,*y;
{
  double tmp;
  /* printf("Entry scambia.\n"); */
  tmp = *x;
  *x = *y;
  *y = tmp;
  /* printf("Exit scambia.\n"); */
  return;
}
/* ### */
void scomponi(flag_prudente, n, vetfreq, coord_base_non_ort)
     unsigned short int flag_prudente;
     int n;
     double vetfreq[NFREQMAX], coord_base_non_ort[NFREQMAX][2];
     /* Routine che individua n velocita' angolari (vetfreq[0], ...,
	vetfreq[n-1]) corrispondenti a dei massimi relativi
	dell'ampiezza dell'integrale fondamentale dell'analisi in
	frequenza.  Inoltre, la presente routine calcola le coordinate
	del segnale espresse rispetto alla base NON ortonormale
	costituita dai seguenti vettori
	   exp(i vetfreq[0] t), ..., exp(i vetfreq[n-1] t),
	dove i e' l'unita' immaginaria e t e' il tempo.

	     *** Attenzione ***     In uscita da questa routine il segnale
	                            che viene restituito e' quello di entrata
				    meno la sua stessa proiezione (calcolata in
				    modo iterativo e con l'ausilio della base
				    ortonormale) sullo spazio lineare generato
				    da exp(i vetfreq[0] t), ...,
				    exp(i vetfreq[n-1] t.

	Variabili in entrata:
	  flag_prudente          => Variabile flag che se e' uguale a TRUE (=1)
	                            imposta la ricerca dei massimi relativi in
				    modalita' prudente (che e' meno veloce, ma
				    piu' indicata quando ci sono segnali in cui
				    la componente caotica e'rilevante); se invece
				    e' uguale a FALSE (=0) imposta la ricerca dei
				    massimi relativi in modalita' veloce
          n                      => Numero delle componenti del segnale
				    che vogliamo determinare per mezzo
				    dell'analisi in frequenza

	Variabili (di tipo globale) in entrata:
	  segnale[0], ..., segnale[2*npt+1]
	                         => Ognuna delle coppie di elementi consecutivi
				    rappresenta le componenti reali e
				    immaginarie del segnale iniziale che
				    vogliamo studiare con l'analisi in
				    frequenza

	Variabili in uscita:
	  vetfreq[0], ..., vetfreq[n-1]
				 => Velocita' angolari corrispondenti a
				    dei massimi relativi dell'ampiezza
				    dell'integrale fondamentale
				    dell'analisi in frequenza
	  coord_base_non_ort[0][0] ,   coord_base_non_ort[0][1]
                    .                            .
                    .                            .
	  coord_base_non_ort[n-1][0] , coord_base_non_ort[n-1][1]
	                         => Componenti reali e immaginarie delle
				    coordinate del segnale rispetto agli
				    n vettori della base NON ortonormale
				    
	Variabili (di tipo globale) in uscita:
	  segnale[0], ..., segnale[2*npt+1]
	                         => Ognuna delle coppie di elementi consecutivi
				    rappresenta le componenti reali e
				    immaginarie del segnale che abbiamo ricevuto
				    in entrata e a cui e' stata sottratta
				    la proiezione sullo spazio lineare formato
				    dagli n vettori della base

     */
{
  int j;                         /* Indice che funge da contatore all'interno
				    dei cicli for */ 
  double mat_nonort_ort[NFREQMAX][NFREQMAX];
                                 /* Matrice di passaggio dalle coordinate
				    relative alla base non ortonormale
				    a quelle che si riferiscono alla base
				    ortonormale */
  double mat_ort_nonort[NFREQMAX][NFREQMAX];
                                 /* Matrice di passaggio dalle coordinate
				    relative alla base ortonormale a quelle
				    che si riferiscono alla base non ortonormale */
  double coord_base_ort[NFREQMAX][2];
                                 /* Componenti reali e immaginarie
				    delle proiezioni del segnale sui
				    vettori della base ortonormale */
  /* printf("Entry scomponi.\n"); */
  for(j=0; j<n; j++) {
    /* Determinazione di un massimo relativo dell'ampiezza dell'integrale
       fondamentale dell'analisi in frequenza */
    vetfreq[j]=max_segnale(flag_prudente, ERRSOLUZ);
    // printf("max_segnale in omega = %24.16le\n", vetfreq[j]);
    /* Aggiornamento delle matrici che consentono il passaggio dalle coordinate
       della base non ortonormale a quella ortonormale */
    ortonormalizza(mat_nonort_ort, mat_ort_nonort, j+1, vetfreq);
    /* Calcolo della proiezione del segnale sul j-esimo vettore
       della base ortonormale */
    proiezione_base_ort(vetfreq[j], mat_nonort_ort[j][j], coord_base_ort[j]);
    /* Sottrazione al segnale della sua stessa proiezione sul j-esimo
       vettore della base ortonormale */
    sottrai_al_segnale(mat_ort_nonort, vetfreq, j+1, coord_base_ort[j]);
  }
  /* Calcolo delle coordinate del segnale rispetto alla base non
     ortonormale */
  calcola_coord_base_non_ort(mat_ort_nonort, n,
			     coord_base_ort, coord_base_non_ort);
  /* printf("Exit scomponi.\n"); */
}
/* ### */
void sottrai_al_segnale(mat_ort_nonort, vetfreq, n, coord_vet_n)
     double mat_ort_nonort[NFREQMAX][NFREQMAX], vetfreq[NFREQMAX];
     int n;
     double coord_vet_n[2];
     /* Routine che sottrae al segnale la sua stessa proiezione
	sull'n-esimo vettore della base ortonormale. Si ricordi che la
	base ortonormale e' quella che si ottiene secondo l'algoritmo
	di Gram-Schmidt dai vettori seguenti
	   exp(i vetfreq[0] t), ..., exp(i vetfreq[n-1] t),
	dove i e' l'unita' immaginaria, t e' il tempo e vetfreq e' il
	vettore di frequenze tale che ciascuna delle sue componenti
	corrisponde a un massimo relativo dell'integrale fondamentale
	dell'analisi in frequenza.

	Variabili in entrata:
	  mat_ort_nonort[0][0], ...,  mat_ort_nonort[0][n-1]
                   .              .             .
                   .              .             .
	  mat_ort_nonort[n-1][0], ...,  mat_ort_nonort[n-1][n-1]
	                         => Matrice di passaggio dalle coordinate
				    relative alla base ortonormale a quelle
				    che si riferiscono alla base non ortonormale
	  vetfreq[0], ...,  vetfreq[n-1]
	                         => Valori delle frequenze che corrispondono
				    a dei massimi relativi dell'integrale
				    fondamentale dell'analisi in frequenza
	  n                      => Numero delle frequenze determinate finora
	                            tra quelle che corrispondono a dei massimi
				    relativi dell'integrale fondamentale
				    dell'analisi in frequenza
	  coord_vet_n[0] , coord_vet_n[1]
	                         => Componenti reale e immaginaria della
				    proiezione del segnale sull'n-esimo
				    vettore della base ortonormale
				    
	Variabili (di tipo globale) in entrata:
	  tgrande                => Tempo pari a meta' dell'intervallo di
	                            integrazione, cioe' l'integrale fondamentale
				    dell'analisi in frequenza viene calcolato
				    tra - tgrande e + tgrande
	  npt                    => Numero dei punti su cui calcoliamo
	                            approssimativamente l'integrale fondamentale
				    dell'analisi in frequenza
	  tstep                  => Tempo di scansione dei punti che ci consentono
	                            di calcolare approssimativamente l'integrale
				    fondamentale dell'analisi in frequenza;
				    tstep altro non e' che 2 * tgrande / npt
	  segnale[0], ..., segnale[2*npt+1]
	                         => Ognuna delle coppie di elementi consecutivi
				    rappresenta le componenti reali e
				    immaginarie del segnale che e' stato
				    letto inizialmente e a cui e' stata sottratta
				    la sua stessa proiezione sullo spazio
				    generato dai primi n-1 vettori della
				    base ortonormale
				    
	Variabili (di tipo globale) in uscita:
	  segnale[0], ..., segnale[2*npt+1]
	                         => Ognuna delle coppie di elementi consecutivi
				    rappresenta le componenti reali e
				    immaginarie del segnale che abbiamo ricevuto
				    in entrata e a cui e' stata sottratta la sua
				    stessa proiezione sull'n-esimo vettore
				    della base ortonormale */
{
  unsigned long int i;
  int j, l;
  double coefcompl[2], t, costmp, sintmp, term_sottr[2];
  /* printf("Entry sottrai_al_segnale\n"); */
  /* Allo scopo di semplificare la notazione nei commenti di questa
     routine, d'ora in poi sia u_0, ..., u_n-1 la base non ortonormale
     dei vettori e sia e_0, ..., e_n-1 quella ortonormale. Inoltre,
     sia f=f(t) il segnale che vogliamo scomporre. */
  /* Dobbiamo calcolare f - (f , e_n-1) e_n-1.
     A questo scopo, siccome e' piu' facile effettuare i calcoli con i
     vettori della base non ortonormale, utilizziamo l'equazione
       e_n-1 = a_n-1,n-1 u_n-1 + a_n-1,n-2 u_n-2 + ... + a_n-1,0 u_0 ,
     dove i coefficienti a_n-1,j altro non sono che i valori di
     mat_ort_nonort[n-1][j]. */
  for(j=0; j<n; j++) {
    /* Calcoliamo a_n-1,j (f , e_n-1) */
    for(l=0; l<2; l++) {
      coefcompl[l] =  mat_ort_nonort[n-1][j] * coord_vet_n[l];
    }
    for(i=0; i<=npt; i++) {
      /* Calcolo del tempo t, di cos(vetfreq[j] t) e sin(vetfreq[j] t) */
      t = i * tstep - tgrande ;
      costmp = cos(vetfreq[j] * t);
      sintmp = sin(vetfreq[j] * t);
      /* Calcolo di a_n-1,j (f , e_n-1) u_j(t).
	 Si ricordi che Re(u_j(t)) = cos(vetfreq[j] t)
	 e Im(u_j(t)) = sin(vetfreq[j] t) */
      term_sottr[0] = coefcompl[0] * costmp - coefcompl[1] * sintmp;
      term_sottr[1] = coefcompl[0] * sintmp + coefcompl[1] * costmp;
      /* Sottrazione a f(t) del termine a_n-1,j (f , e_n-1) u_j(t) */
      for(l=0; l<2; l++) {
	segnale[2*i+l] -= term_sottr[l];
      }
    }
  }
  /* printf("Exit sottrai_al_segnale\n"); */
}
