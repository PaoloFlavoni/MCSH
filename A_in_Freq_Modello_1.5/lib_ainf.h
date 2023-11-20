/* Definizione di alcune costanti utili */
#define TRUE 1
#define FALSE 0
#define MAXITERNEWTON 15 /* Massimo numero di iterazioni all'internodel metodo di Newton */
#define ERRSOLUZ 1.e-15  /* Errore ammesso nella determinazione dei punti dimassimo relativo quando viene chiamata la routinemax_segnale */
/*	     *** Attenzione ***     Dalla definizione delle prossime due costanti
	                            dipende la memoria necessaria al programma.
	     *** Consiglio ***      Non e' conveniente porre un valore di NFREQMAX
	                            molto elevato (e' difficile pensare che
				    possano venire utili applicazioni con
				    NFREQMAX > 100).
				    Pertanto, si puo' con buona approssimazione
				    stimare la memoria necessaria al programma
				    con il seguente calcolo:
				          4 x NPTMAX x sizeof(TIPO) ,
				    dove TIPO e' il tipo di variabile con cui
				    trattiamo il segnale ed effettuiamo i calcoli
				    (cioe' puo' essere float o double oppure
				    long double, etc.) e sizeof e' il numero di
				    bytes relativo alla taglia di TIPO */
#define NPTMAX 33554432  /* Limite superiore sulla variabile npt (vedi sotto) */
#define NFREQMAX 100  /* Limite superiore sulla variabile nfreqstudio (vedi sotto) */
/* Dichiarazione delle variabili globali */
double pig, dupig;               /* Contengono rispettivamente i valori
				    di pigreco e 2 * pigreco */
double tgrande;                  /* Tempo pari a meta' dell'intervallo di
	                            integrazione, cioe' l'integrale fondamentale
				    dell'analisi in frequenza viene calcolato
				    tra - tgrande e + tgrande */
unsigned long int npt;           /* Numero dei punti dell'orbita discretizzata
				    per cui calcoliamo l'integrale
				    fondamentale dell'analisi in frequenza */
double tstep;                    /* Tempo di scansione dei punti che ci consentono
	                            di calcolare approssimativamente l'integrale
				    fondamentale dell'analisi in frequenza;
				    tstep e' uguale a 2 * tgrande / npt */
double segnale[2*NPTMAX+2];      /* segnale[0], ..., segnale[2*npt+1]:
	                            Ognuna delle coppie di elementi consecutivi
				    rappresenta le componenti reali e
				    immaginarie del segnale che vogliamo
				    studiare con l'analisi in frequenza */
double dati_fft[2*NPTMAX];       /* dati_fft[0], ..., dati_fft[2*ncoppie-1] sono
				    ncoppie di numeri complessi (se l'indice
				    e' pari allora si intende che e' la componente
				    reale, altrimenti e' quella immaginaria)
				    cui vogliamo applicare la trasformata di
				    Fourier discreta. Dopo l'esecuzione della
				    routine dfour1, il vettore dati_fft ospita
				    la trasformata di Fourier degli stessi
				    valori di dati_fft precedenti alla chiamata */
/* Dichiarazione delle routines della libreria */
void calcola_coord_base_non_ort(double mat_ort_nonort[][NFREQMAX], int n,
				double coord_base_ort[][2],
				double coord_base_non_ort[][2]);
void delimita_max(double (*funz)(), short unsigned int flag_bisez,
		  double passo_x, double *xa, double *funz_xa,
		  double *d_funz_xa, double *d2_funz_xa, double *xb);
double deriv_prodscal(double nu, double *deriv1, double *deriv2);
void dfour1(unsigned long int ncoppie, int segno);
int leggi_segnale(FILE *in, unsigned long int numfin, int numletture);
double max_fft(unsigned long int ncoppie);
double max_segnale(unsigned short int flag_veloce, double err);
void mira_freq(int nfreq_fond, int modkmax, double omega,
	       double omega_fond[], int kbest[], double *errmin);
double newton(double (*funz)(), double xiniz, double err);
double norma_in_C(double z[]);
void ortonormalizza(double mat_nonort_ort[][NFREQMAX],
		    double mat_ort_nonort[][NFREQMAX],
		    int n, double vetfreq[]);
double pesow(double t);
void prodscal(double omega, double valint[]);
double prodscal_base_nonort(double omega1, double omega2);
void proiezione_base_ort(double omega, double mat_nonort_ort_n,
			 double coord_base_ort_n[]);
void prova_ainf(double vetfreq[], int n, double coord_base_non_ort[][2]);
void scambia(double *x, double *y);
void scomponi(unsigned short int flag_veloce, int n, double vetfreq[],
	      double coord_base_non_ort[][2]);
void sottrai_al_segnale(double mat_ort_nonort[][NFREQMAX], double vetfreq[],
			int n, double coord_vet_n[]);