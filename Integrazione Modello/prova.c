#include<stdio.h>
#include<math.h>
#include<stdlib.h>


int main()
{
    unsigned long int i,N, num_intervalli, num_sezioni;
    double dt, pi, omega, t;
    char riga[200];
    FILE *otp;
    pi = 4*atan(1);
    omega = 1.41421356237309504880;
    t = 0;
    N = 10;
    num_sezioni = 10;
    num_intervalli = N*num_sezioni;
    dt = 2*pi/(omega*N);
    otp = fopen("test.txt", "w");
    for(i = 1; i<= num_intervalli; i++)
    {
        t += dt;
        if(i%N == 0)
        {
            fprintf(otp, "%ld %le\n", i, t);
        }
    }
    fclose(otp);
}