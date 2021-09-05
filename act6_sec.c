#include<stdio.h>
#include "omp.h"
#include <math.h>

#define PI 3.141592

void runge1(FILE *fptr);
void runge2(FILE *fptr);
void runge3(FILE *fptr);
void runge4(FILE *fptr);

int main(){
	const double startTime = omp_get_wtime();
	
	(void)runge1(fopen("RungeSEC_f1.txt", "w"));
	(void)runge2(fopen("RungeSEC_f2.txt", "w"));
	(void)runge3(fopen("RungeSEC_f3.txt", "w"));
	(void)runge4(fopen("RungeSEC_f4.txt", "w"));
		
	const double endTime = omp_get_wtime();
    printf("tomo = (%lf) segundos\n",(endTime - startTime));
}

void runge1(FILE *fptr){
	long N = 50000;
	//fprintf(fptr, "Numero de pasos:%d Atendido por thread:%d\n", N,omp_get_thread_num());
 	fprintf(fptr, "Datos que encuentra el metodo de Euler(variable ind.\t variable dep.)\n");
 	   	double h,t,w;
       	double k1= 0.0, k2= 0.0, k3=0.0, k4=0.0;
       	double w0=PI/4, t0=0,a=0.0,b=PI;
       	double ab=0.0;
       	int i;
       	w=w0;
       	fprintf(fptr, "%f\t %f\n", a, w);
       	for(i=0;i<N;i++){
        	h=(b-a)/N;
           	t=a+(h*i);
           	ab=t*t;
           	k1=h*(t*exp(3*t)-2*w);
           	k2=h*((t+0.5*h)*exp(3*(t+0.5*h))-2*(w+0.5*k1));
           	k3=h*((t+0.5*h)*exp(3*(t+0.5*h))-2*(w+0.5*k2));
           	k4=h*((t+h)*exp(3*(t+h))-2*(w+k3));
           	w=w+((1.0/6.0)*(k1+(2*k2)+(2*k3)+k4));
           	fprintf(fptr, "%f\t %f \t numero de thread:%d\n", t+h, w,omp_get_thread_num());
        }
   fclose(fptr);
}

void runge2(FILE *fptr){
	long N = 50000;
	//fprintf(fptr, "Numero de pasos:%d Atendido por thread:%d\n", N,omp_get_thread_num());
 	fprintf(fptr, "Datos que encuentra el metodo de Euler(variable ind.\t variable dep.)\n");
 	   	double h,t,w;
       	double k1= 0.0, k2= 0.0, k3=0.0, k4=0.0;
       	double w0=PI/4, t0=0,a=0.0,b=PI;
       	double ab=0.0;
       	int i;
       	w=w0;
       	fprintf(fptr, "%f\t %f\n", a, w);
       	for(i=0;i<N;i++){
        	h=(b-a)/N;
           	t=a+(h*i);
           	ab=t*t;
           	k1=h*(1.0+pow(t-w,2));
           	k2=h*(1.0+pow((t+0.5*h)-(w+0.5*k1),2));
           	k3=h*(1.0+pow((t+0.5*h)-(w+0.5*k2),2));
           	k4=h*(1.0+pow((t+h)-(w+k3),2));
           	w=w+((1.0/6.0)*(k1+(2*k2)+(2*k3)+k4));
           	fprintf(fptr, "%f\t %f \t numero de thread:%d\n", t+h, w,omp_get_thread_num());
        }
   fclose(fptr);
}

void runge3(FILE *fptr){
	long N = 50000;
	//fprintf(fptr, "Numero de pasos:%d Atendido por thread:%d\n", N,omp_get_thread_num());
 	fprintf(fptr, "Datos que encuentra el metodo de Euler(variable ind.\t variable dep.)\n");
 	   	double h,t,w;
       	double k1= 0.0, k2= 0.0, k3=0.0, k4=0.0;
       	double w0=PI/4, t0=0,a=0.0,b=PI;
       	double ab=0.0;
       	int i;
       	w=w0;
       	fprintf(fptr, "%f\t %f\n", a, w);
       	for(i=0;i<N;i++){
        	h=(b-a)/N;
           	t=a+(h*i);
           	ab=t*t;
           	k1=h*(1.0+w/t);
           	k2=h*(1.0+(w+0.5*k1)/(t+0.5*h));
           	k3=h*(1.0+(w+0.5*k2)/(t+0.5*h));
           	k4=h*(1.0+(w+k3)/(t+h));
            w=w+((1.0/6.0)*(k1+(2*k2)+(2*k3)+k4));
           	fprintf(fptr, "%f\t %f \t numero de thread:%d\n", t+h, w,omp_get_thread_num());
        }
   fclose(fptr);
}

void runge4(FILE *fptr){
	long N = 50000;
	//fprintf(fptr, "Numero de pasos:%d Atendido por thread:%d\n", N,omp_get_thread_num());
 	fprintf(fptr, "Datos que encuentra el metodo de Euler(variable ind.\t variable dep.)\n");
 	   	double h,t,w;
       	double k1= 0.0, k2= 0.0, k3=0.0, k4=0.0;
       	double w0=PI/4, t0=0,a=0.0,b=PI;
       	double ab=0.0;
       	int i;
       	w=w0;
       	fprintf(fptr, "%f\t %f\n", a, w);
       	for(i=0;i<N;i++){
        	h=(b-a)/N;
           	t=a+(h*i);
           	ab=t*t;
           	k1=h*(cos(2*t*w)+sin(3*t*w));
           	k2=h*(cos(2*(t+0.5*h)*(w+0.5*k1))+sin(3*(t+0.5*h)*(w+0.5*k1)));
           	k3=h*(cos(2*(t+0.5*h)*(w+0.5*k2))+sin(3*(t+0.5*h)*(w+0.5*k2)));
           	k4=h*(cos(2*(t+h)*(w+k3))+sin(3*(t+h)*(w+k3)));
           	w=w+((1.0/6.0)*(k1+(2*k2)+(2*k3)+k4));
           	fprintf(fptr, "%f\t %f \t numero de thread:%d\n", t+h, w,omp_get_thread_num());
        }
   fclose(fptr);
}