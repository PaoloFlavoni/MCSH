#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define pi 4*atan(1)




int main()
{
    unsigned long int i, j, num_intervalli, k, i_graph;
    double  epsilon, Tf, dt, t, dt_graph;
    double x, y, omega[2];
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
    sscanf(riga, "%lf", &Tf);
    printf("Valore del tempo d'integrazione Tf %f\n\n", Tf);
    fgets(riga, 200, ipf);
    fgets(riga, 200, ipf);
    sscanf(riga, "%lf", &dt);
    printf("Valore del tempuscolo d'integrazione dt %f\n\n", dt);   
    num_intervalli = Tf/dt;
    printf("Valore Numero di Intervalli %ld\n\n", num_intervalli);
    fgets(riga, 200, ipf);
    fgets(riga, 200, ipf);
    sscanf(riga, "%lf", &dt_graph);
    printf("Valore di scansionamento del segnale %f\n\n", dt_graph);
    printf("Effettuo l'integrazione numerica delle equazioni del moto\n\n");
    /*Effettuo una lettura a vuoto del file di input*/
    fgets(riga, 200, ipf);
    /*Uso il while in maniera da poter avere diverse condizionin inziali da vedere*/
    t = 0 ;
    /* Pongo i_graph uguale al numero di passi che devono in tercorrere tra
    il tracciamento di un punto per il segnale e il successivo*/
    i_graph = dt_graph / dt;
    while (fgets(riga, 200, ipf) != NULL)
    {
         /*Apro file di output*/
        opf = fopen("pm_pendolo_eulero.out", "w");
        sscanf(riga, "%le %le", &x, &y);
        printf("Condizioni iniziali %le %le\n", x, y);
    for(i = 1; i <= num_intervalli; i++)
    {   
        
        x = x + y * dt ;
        if (x > pi)
        x -=  2 * pi;
        if (x < - pi)
        x += 2 *  pi;
        y = y - omega[0] * omega[0] * sin(x) * dt + epsilon * cos(omega[1] * t) * dt ;
        t += dt ;
        if(i % i_graph == 0)
        {
                    fprintf(opf, "%le %le \n", y, x);
        }
        if(i%(num_intervalli/100) == 0)
        {
            printf("INTERVALLI RIMANENTI  %li  \n\n", num_intervalli - i);
        }
    
    }
    }
    fclose(ipf);
    fclose(opf);
}