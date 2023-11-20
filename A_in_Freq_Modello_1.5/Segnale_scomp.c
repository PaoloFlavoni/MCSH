/*Programma che controlla se la divisione in parte reale e complessa del segnale 
pm_pendolo_eulero.out, fatta tramite il programma scomposizione.c, usando 
Mathematica è corretta*/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

double Re(double t);
double Im(double t);

int main()
{
   unsigned long int i, num_intervalli, i_get;
   double pi, dt, x_2dot, x_2, t, k, omega, T, dt_get;
   char riga[200];
   FILE *inf, *otf;
   pi = 4*atan(1);
   inf = fopen("input.dat", "r");
   fgets(riga, 200, inf);
   fgets(riga, 200, inf);
   sscanf(riga, "%lf", &omega);
   printf("Valore di omega %lf\n\n", omega);
   fgets(riga, 200, inf);
   fgets(riga, 200, inf);
   sscanf(riga, "%lf", &dt);
   printf("Valore del tempuscolo d'integrazione %lf \n\n", dt);
   fgets(riga, 200, inf);
   fgets(riga, 200, inf);
   sscanf(riga, "%lf", &T);
   printf("Valore del Tempo d'integrazione %lf \n\n", T);
   fgets(riga, 200, inf);
   fgets(riga, 200, inf);
   sscanf(riga, "%lf", &dt_get);
   printf("Valore Tempo di scansionamento del segnale %lf \n\n", dt_get);
   num_intervalli = T/dt;
   printf("Valore Numero di Intervalli %ld \n\n", num_intervalli);   
   printf("Effettuo l'integrazione numerica delle equazioni del moto\n\n");
   otf = fopen("orbite_1.dat", "w");
   fgets(riga, 200, inf);
   fgets(riga, 200, inf);
   sscanf(riga, "%lf", &k);
   printf("Valore della Costante Elastica %f\n\n", k);
   t = 0;
   fgets(riga, 200, inf);
   /* Pongo i_graph uguale al numero di passi che devono intercorrere
     tra il tracciamento di un punto per il grafico e il successivo */
   i_get = dt_get / dt;
   while(fgets(riga, 200, inf) != NULL)
   {
      sscanf(riga, "%le %le", &x_2, &x_2dot);
      printf("Condizioni iniziali %le %le\n", x_2, x_2dot);
      fprintf(otf, " %le %le \n", x_2, x_2dot);
      for(i=1; i<= num_intervalli; i++)
      {
         
         x_2 = x_2 + dt * x_2dot;
         x_2dot = x_2dot - dt * ((x_2 - pi/2)- (x_2 - pi/2)*(x_2 - pi/2)) + k*dt*Im(t) ;
         //fprintf(otf, "%le %le  \n",x_2, x_2dot);
         t += dt;
         if(i % i_get == 0)
         {
            fprintf(otf, " %le %le \n", x_2, x_2dot);
         }
         
      }
      fprintf(otf, "\n\n");
   }
   
   fclose(inf);
   fclose(otf);
}


 double Re(double t)
 {
    double y;

    y = 1.70337 + 0.0057001 * cos(1.841 - 0.64191560007675122 * t) + 
 0.0035484 * cos(1.456 - 0.53492966673062602 * t) + 
 0.085243 * cos(0.4499 - 0.32095780003837561 * t) + 
 0.12178 * cos(1.837 - 0.21397186669225041 * t) + 
 0.16103 * cos(2.454 - 0.21397186669225041 * t) + 
 1.0003 * cos(2.379 - 0.10698593334612520 * t) + 
 0.86358 * cos(3.034 + 0.10698593334612520 * t) + 
 0.042899 * cos(2.199 + 0.32095780003837561 * t) + 
 0.032066 * cos(0.9495 + 0.42794373338450081 * t) + 
 0.019363 * cos(2.643 + 0.42794373338450081 * t) + 
 0.012128 * cos(2.838 + 0.53492966673062602 * t) + 
 0.0034397 * cos(0.2395 + 0.64191560007675122 * t);
    return(y);
 }

 double Im(double t)
 {
    double x;

    x = 1.45053 - 0.0057001 * sin(1.841 - 0.64191560007675122 * t) + 
 0.0035484 * sin(1.456 - 0.53492966673062602 * t) - 
 0.085243 * sin(0.4499 - 0.32095780003837561 * t) + 
 0.12178 * sin(1.837 - 0.21397186669225041 * t) - 
 0.16103 * sin(2.454 - 0.21397186669225041 * t) + 
 1.0003 * sin(2.379 - 0.10698593334612520 * t) + 
 0.86358 * sin(3.034 + 0.10698593334612520 * t) - 
 0.042899 * sin(2.199 + 0.32095780003837561 * t) + 
 0.032066 * sin(0.9495 + 0.42794373338450081 * t) - 
 0.019363 * sin(2.643 + 0.42794373338450081 * t) + 
 0.012128 * sin(2.838 + 0.53492966673062602 * t) - 
 0.0034397 * sin(0.2395 + 0.64191560007675122 * t);
        return(x);
 }