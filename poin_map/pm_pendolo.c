/*pf.c: Programma che consente di calcolare la mappa di POincarè di un pendolo forzato
	    per un orbita di librazione*/


#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define pi 4*atan(1)
/*Paradigma di programmazione Top-Down*/

/*Definizione delle funzioni*/
void kin_flow(double x[2], double y[2],double omega[2], double dt, double epsilon);
void pot_flow(double x[2], double y[2], double dt,  double omega[2], double epsilon);
void saba1(double x[2], double y[2],double omega[2],double c_dt[2], double d_dt[2], double epsilon);

/* ### */

int main()
{
    unsigned long int i, j, N, num_sezioni, num_intervalli, i_graph, k;
    double  epsilon, Tf, dt;
    double en_t, en_0, errmax_en, delta_en;
    double x[2], y[2], xiniz[2], yiniz[2], omega[2];
    double c[2], d[2], c_dt[2], d_dt[2], dt_graph;
    char riga[200];
    FILE *ipf, *opf;
    printf("Questo è pigreco %lf\n\n", pi);
    ipf = fopen("pm_pendolo.inp", "r");
    printf("Leggo da file esterno i parametri e le condizioni iniziali\n\n");
    /*Inizio delle operazioni di lettura da file di input*/
    /*Come al solito leggo prima la riga di commenti e poi i valori*/
    fgets(riga, 200, ipf);
    fgets(riga, 200, ipf);
    sscanf(riga, "%lf %lf", &omega[0], &omega[1]);
    printf("Valore di omega_1 %f, Valore di omega_2 %f\n\n", omega[0], omega[1]);
    fgets(riga, 200, ipf);
    fgets(riga, 200, ipf);
    sscanf(riga, "%lf", &epsilon);
    printf("Valore piccolo parametro %f\n\n", epsilon);
    fgets(riga, 200, ipf);
    fgets(riga, 200, ipf);
    sscanf(riga, "%ld", &N);
    printf("Valore di N %ld\n\n", N);
    fgets(riga, 200, ipf);
    fgets(riga, 200, ipf);
    sscanf(riga, "%ld", &num_sezioni);
    printf("Valore del numero di Sezioni %ld\n\n", num_sezioni);
    /*Diamo le condizoni iniziali sul primo e sul secondo momento*/
    xiniz[0] = 0. ;
    yiniz[0] = 0. ;
    /* Pongo i valori degli elementi dei vettori c e d uguali alle
     costanti relative al metodo leap-frog, che e' anche il SABA1,
     utilizzando la notazione in accordo con l'articolo di Laskar e
     Robutel del 2001 su Celestial Mechanics and Dynamical
     Astronomy) */
    c[0] = 0.5;
    d[0] = 1.;
    /*Pongo num internalli pari al seguente prodotto*/
    num_intervalli = N*num_sezioni;
    printf("Valore Numero di Intervalli %ld\n\n", num_intervalli);
    /*Pongo pari dt pari alla seguente quantità*/
    dt = 2*pi/(N*omega[1]);
    printf("Valore tempuscolo %f\n\n",dt);
    /* Calcolo i valori degli elementi dei vettori c_dt e d_dt in modo
    che corrispondano ai valori dei passi di integrazione che
    vengono effettuati (lungo il flusso delle due parti in cui
    viene divisa la Hamiltoniana) durante l'applicazione del metodo
    leap-frog */
    for(i=0; i<1; i++) 
    {
        c_dt[i] = c[i] * dt;
        d_dt[i] = d[i] * dt;
    }
    /*Effettuo l'integrazione delle equazioni del moto*/
    opf = fopen("pm_pendolo.out", "w");
    printf("Effettuo l'integrazione numerica delle equazioni del moto\n\n");
    /*Effettuo una lettura a vuoto del file di input*/
    fgets(riga, 200, ipf);
    /*Uso il while in maniera da poter avere diverse condizionin inziali da vedere*/
    while (fgets(riga, 200, ipf) != NULL)
    {
         /*Apro file di output*/
        sscanf(riga, "%le %le", &xiniz[1], &yiniz[1]);
        printf("Condizioni iniziali %le %le\n", xiniz[1], yiniz[1]);
        for(i=0;i<2;i++)
        {
            x[i] = xiniz[i];
            y[i] = yiniz[i];
        }
        for(j=1; j<= num_intervalli; j++)
        {
            saba1(x,y,omega,c_dt,d_dt,epsilon);
            if (j%N == 0)
            {
                /*Stampo su file esterno i valori della mappa*/
                fprintf(opf, "%le %le \n", x[1], y[1]);
                //k += 1;
            }
        }
    /*
    printf("QUESTO È IL VALORE DI K %ld", k);
    effettivamente l'algoritmo fa quello che ci si aspetta, ovvero raccoglie un numero di punti pari a num_sezioni
    */
    
    }
     fclose(ipf);
     fclose(opf);
}

/*Funzione che integra lungo il flusso della parte cinetica*/
void kin_flow(double x[2], double y[2],double omega[2], double dt, double epsilon)
{
    y[0] += - epsilon * x[1] * omega[1] * cos(omega[1] * x[0]) * dt;
    y[1] += (-omega[0] * omega[0] * sin(x[1]) + epsilon * sin(omega[1] * x[0])) * dt;
   
} 
/*Funzione che integra lungo il flusso del potenziale*/
void pot_flow(double x[2], double y[2], double dt,  double omega[2], double epsilon)
{
   x[0] += dt;
   x[1] += y[1]*dt ; 
    if (x[1] > pi)
        x[1] -=  2 * pi;
    if (x[1] < - pi)
        x[1] += 2 *  pi;
}
/*Integratore Leap-Frog*/
void saba1(double x[2], double y[2],double omega[2],double c_dt[2], double d_dt[2], double epsilon)
{
    kin_flow( x, y, omega, c_dt[0], epsilon);//Flusso lungo A per un tempo dt/2
    pot_flow(x, y, d_dt[0], omega, epsilon);//Flusso lungo B per un tempo dt
    kin_flow( x, y, omega, c_dt[0], epsilon);//Flusso lungo A per un tempo dt/2
}
