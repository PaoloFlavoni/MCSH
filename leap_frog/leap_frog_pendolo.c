/* leap_frog_pendolo.c: Programma che consente di osservare
                        l'andamento dell'errore sull'energia
                        per il metodo SABA1. Test effettuato 
                        controllando la conservazione dell'energia
                        per il pendolo forzato.*/



#include<stdio.h>
#include<math.h>
#include<stdlib.h>

/*Paradigma di programmazione Top-Down*/

/*Definizione delle funzioni*/
double energy(double x[2], double y[2],double omega[2], double epsilon);
void kin_flow(double x[2], double y[2],double omega[2], double dt);
void pot_flow(double x[2], double y[2], double dt,  double omega[2], double epsilon);
void saba1(double x[2], double y[2],double omega[2],double c_dt[2], double d_dt[2], double epsilon);

/* ### */

int main()
{
    int i, N, num_intervalli, i_graph;
    double  epsilon, Tf, dt;
    double en_t, en_0, errmax_en, delta_en;
    double x[2], y[2], xiniz[2], yiniz[2], omega[2];
    double c[2], d[2], c_dt[2], d_dt[2], dt_graph;
    char riga[200];
    FILE *ipf, *opf;
    ipf = fopen("leap_frog_pendolo.inp", "r");
    printf("Leggo da file esterno i parametri e le condizioni iniziali\n\n");
    /*Inizio delle operazioni di lettura da file di input*/
    /*Come al solito leggo prima la riga di commenti e poi i valori*/
    fgets(riga, 200, ipf);
    fgets(riga, 200, ipf);
    sscanf(riga, "%lf %lf", &omega[0], &omega[1]);
    fgets(riga, 200, ipf);
    fgets(riga, 200, ipf);
    sscanf(riga, "%lf", &epsilon);
    fgets(riga, 200, ipf);
    fgets(riga, 200, ipf);
    sscanf(riga, "%lf %lf", &xiniz[0], &yiniz[0]);
    fgets(riga, 200, ipf);
    fgets(riga, 200, ipf);
    sscanf(riga, "%lf", &Tf);
    fgets(riga, 200, ipf);
    fgets(riga, 200, ipf);
    sscanf(riga, "%lf", &dt);
    fgets(riga, 200, ipf);
    fgets(riga, 200, ipf);
    sscanf(riga, "%lf", &dt_graph);
    fclose(ipf);
    /*Diamo le condizoni iniziali sul primo e sul secondo momento*/
    xiniz[1] = 0. ;
    yiniz[1] = 1. ;
    /*Calcolo il valore inizale dell'energia*/
    en_0 = energy(xiniz, yiniz, omega, epsilon);
    /* Pongo i valori degli elementi dei vettori c e d uguali alle
     costanti relative al metodo leap-frog, che e' anche il SABA1,
     utilizzando la notazione in accordo con l'articolo di Laskar e
     Robutel del 2001 su Celestial Mechanics and Dynamical
     Astronomy) */
    c[0] = 0.5;
    d[0] = 1.;
    /*Poniamo i valori delle coordinate canoniche uguali ai 
    corrispondenti valori iniziali letti dal file di input*/
    for(i=0; i<2; i++)
    {
        x[i] = xiniz[i];
        y[i] = yiniz[i];
    }
    /*POngo num internalli pari al seguente rapporto*/
    num_intervalli = Tf / dt;
    /*Pongo i_graph uguale al numero di passi cge devono intercorrere tra il tracciamento
    di un punto per il grafico e il successivo*/
    i_graph = dt_graph / dt;
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
    /*Apro file di output*/
    opf = fopen("leap_frog_pendolo.out", "w");
    /*Stampo su file esterno il valore iniziale dell'errore, che è zero ovviamente*/
    fprintf(opf, " %le %le\n", 0.0, 0.0);
    /*Coerentemente pongo a zero il valore dell'errore massimo sull'energia*/
    errmax_en = 0.;
    /*Effettuo l'integrazione delle equazioni del moto*/
    printf("Effettuo l'integrazione numerica delle equazioni del moto\n\n");
    for(i=1; i<=num_intervalli; i++)
    {
        saba1(x , y, omega, c_dt, d_dt, epsilon);
            if(i % i_graph == 0)
            {
                en_t = energy(x, y, omega, epsilon);
                delta_en = en_t - en_0;
                fprintf(opf, "%le  %le\n", i * dt, delta_en );
                if(errmax_en < fabs(delta_en))
                {
                    errmax_en = fabs(delta_en);
                }
            }
    }
    fclose(opf);
    /* Stampa finale a proposito della conservazione dell'energia */
    printf(" En. iniz.    max_t | E(t) - E(0) |   max_t | (E(t)-E(0))/E(0) |\n");
    printf("%le       %le              %le\n\n",
	 en_0, errmax_en, fabs(errmax_en/en_0) );
    /*in un certo
    intervallo dei valori dei passi di integrazione si ha che dimezzando il
    passo di integrazione l'erroresi riduce di un fattore 4, perché il
    metodo leap-frog  è quadratico*/
    //spoiler : pare funzioni
}
    

double energy(double x[2], double y[2],double omega[2], double epsilon)
{
    double ris;
    ris = 0. ;
    ris = omega[1] * y[0] + y[1] * y[1] /2 - cos(x[1]) - epsilon * x[1] * sin(x[0]);
    return ris;
}
    

/*Funzione che integra lungo il flusso della parte cinetica*/
void kin_flow(double x[2], double y[2],double omega[2], double dt)
{
   x[0] += dt;
   x[1] += omega[0]*y[1]*dt;
}
/*Funzione che integra lungo il flusso del potenziale*/
void pot_flow(double x[2], double y[2], double dt,  double omega[2], double epsilon)
{
    y[0] +=   epsilon * x[1] * omega[1]* (cos( omega[1] * x[0])) * dt;
    y[1] +=  (epsilon * (sin( omega[0] * x[0])) - omega[0]*sin(x[1]) ) * dt;
}
/*Integratore Leap-Frog*/
void saba1(double x[2], double y[2],double omega[2],double c_dt[2], double d_dt[2], double epsilon)
{
    kin_flow( x, y, omega, c_dt[0]);//Flusso lungo A per un tempo dt/2
    pot_flow(x, y, d_dt[0], omega, epsilon);//Flusso lungo B per un tempo dt
    kin_flow( x, y, omega, c_dt[0]);//Flusso lungo A per un tempo dt/2
}
