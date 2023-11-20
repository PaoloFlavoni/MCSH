#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define pi 4*atan(1)




int main()
{
    unsigned long int i, j, N, num_sezioni, num_intervalli, i_graph, k;
    double  epsilon, Tf, dt, t;
    double x_p, x_d, y_p, y_d, omega[2];
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
    num_intervalli = N*num_sezioni;
    printf("Valore Numero di Intervalli %ld\n\n", num_intervalli);
    /*Pongo pari dt pari alla seguente quantità*/
    dt = 2*pi/(N*omega[1]);
    opf = fopen("pm_pendolo_eulero.out", "w");
    printf("Effettuo l'integrazione numerica delle equazioni del moto\n\n");
    /*Effettuo una lettura a vuoto del file di input*/
    fgets(riga, 200, ipf);
    /*Uso il while in maniera da poter avere diverse condizionin inziali da vedere*/
    t = 0 ;
    while (fgets(riga, 200, ipf) != NULL)
    {
         /*Apro file di output*/
        sscanf(riga, "%le %le", &x_d, &y_d);
        printf("Condizioni iniziali %le %le\n", x_d, y_d);
    for(i = 1; i <= num_intervalli; i++)
    {   
        
        x_d = x_d + y_d*dt ;
        if (x_d > pi)
        x_d -=  2 * pi;
        if (x_d < - pi)
        x_d += 2 *  pi;
        y_d = y_d - omega[0] * omega[0] * sin(x_d) * dt + epsilon * cos(omega[1] * t) * dt ;
        t += dt ;
        if (i%N == 0)
            {
                //printf("valore di t %f\n\n", t);
                /*Stampo su file esterno i valori della mappa*/
                fprintf(opf, "%le %le \n", x_d, y_d);
                //k += 1;
            }
     
    }
    }
    fclose(ipf);
    fclose(opf);
}

