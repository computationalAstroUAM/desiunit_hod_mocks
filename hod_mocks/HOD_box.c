/* This program reads the OuterRim halo catalogue directly in it's GenericIO format.
 */
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <cassert>
#include "/mnt/lustre/eboss/OuterRim/genericio/GenericIO.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

int main() {
   /* my first program in C */
   printf("Hello, World! \n");
   
   return 0;
}

/*
#define rho_crit (27.755e10)


#define outbase ("output/galaxies_V1.0_M1_%.2f_DLM%.2f_Lth%.2f.dat")

#define input ("/mnt/lustre/eboss/OuterRim/OuterRim_sim/ascii/vol200Mpc/OuterRim_STEP266_z0.865/OuterRim_STEP266_fofproperties_200Mpc.txt")
float zsnap = 0.865;

double LBOX = 200.;
double OMEGA_M = 0.25;

int Nz;
double rad_to_deg(double rad){
	return rad*180.0/M_PI;
}

double deg_to_rad(double dec){
	return dec*M_PI/180.0;
}
struct cosm_params { double Omega_m; double Omega_L; };


//void read_photoz(char *);
double rand_gauss (void);
int countlines(char *);
void NFW_to_pos(double,double *, double *,double *,double );
void vir_to_vel(double ,double,double *,double *,double *);
//void box_to_lightcone(double ,double ,double ,double ,double ,double ,double *,double *,double *);
double rand_photoz(double zspec);
void invert_z_r();
double spline_inter(double *x, double *y, int N, double *y2, double xi);
int HOD_alpha1(double M, double M1);

int CX,CY,CZ;

double *mask;

//void initialise_HMF();
//double HMF1_to_HMF2(double );
//double watson(double );

int main(int argc, char **argv){
  	fprintf(stderr,"Run call:\n\t");
	int i ;
   	for (i=0;i<argc;i++){
			fprintf(stderr,"%s  ",argv[i]);
   	}
	fprintf(stderr,"\n");
	if (argc!=4){
		fprintf(stderr,"usage: %s M1 DeltaLM Lth\n",argv[0]);
		return -1;
	} 


	int seed;
//	double *mask;
	
	int j,Nhalos,Nsat;
	FILE *f_in,*f_out;
	char output[1024],line[1024];
	double M,x,y,z,vx,vy,vz,L;
	double Dx,Dy,Dz,Dvx,Dvy,Dvz;
	double xgal,ygal,zgal;
	double M1,DeltaLM,logLth,logM1,Lth;


	logM1 = atof(argv[1]);
	DeltaLM = atof(argv[2]);
	logLth = atof(argv[3]);

	sprintf(output,outbase,logM1,DeltaLM,logLth);

	if (logM1>0.0)
		M1 = pow(10.0,(double) logM1 );
	else M1=0.;

	DeltaLM = DeltaLM*log(10.0);

	if (logLth>0.0)
		Lth = pow(10.0,(double) logLth );
	else Lth=0.;

	seed = 42;
	srand(seed);

	

	Nhalos = countlines(input);

	fprintf(stderr,"found N=%d halos! in input file %s\n", Nhalos,input);
	
        if ((f_in = fopen(input,"r"))==NULL){
                fprintf(stderr,"Could not open input file %s\n", input);
                exit(0);
        }
	fprintf(stderr,"Writing output in %s...\n",output);

        if ((f_out = fopen(output,"w"))==NULL){
                fprintf(stderr,"Could not open input file %s\n", output);
                exit(0);
        }

	fgets(line,1000,f_in);
	if (line[0]=='#')
		fgets(line,1000,f_in);

	fprintf(stderr,"Starting HOD with M1=%e\n",M1);
	int a;
	double b,c,d,e;
	for (i=0;i<Nhalos;i++){	
                sscanf(line,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&a, &b, &M, &x,&y,&z, &c,&d,&e, &vx,&vy,&vz);
		Nsat = HOD_alpha1(M,M1);

		L = exp(log(M) +  DeltaLM*rand_gauss());
		if (L>Lth)
			fprintf(f_out,"%f %f %f %f %f %f %e %e %d\n",x,y,z,vx,vy,vz,M,L,Nsat);
		
		for (j=0;j<Nsat;j++){
			L = exp(log(M) +  DeltaLM*rand_gauss());
			if (L>Lth){
				NFW_to_pos(M,&Dx,&Dy,&Dz,zsnap);
				vir_to_vel(M,zsnap,&Dvx,&Dvy,&Dvz);
				xgal=x+Dx;
				ygal=y+Dy;
				zgal=z+Dz;
				if (xgal>=LBOX)
					xgal-=LBOX;
				if (xgal<0)
					xgal+=LBOX;
				if (ygal>=LBOX)
					ygal-=LBOX;
				if (ygal<0)
					ygal+=LBOX;
				if (zgal>=LBOX)
					zgal-=LBOX;
				if (zgal<0)
					zgal+=LBOX;
				fprintf(f_out,"%f %f %f %f %f %f %e %e %d\n",xgal,ygal,zgal,vx+Dvx,vy+Dvy,vz+Dvz,M,L,-1);
			}
		}
          	fgets(line,1000,f_in);
	}
	
	fclose(f_in);
	fclose(f_out);

	fprintf(stderr,"\nDone and writen to %s\n\n",output);

	return 0;
}

unsigned long factorial(unsigned long f)
{
    if ( f == 0 ) 
        return 1;
    return(f * factorial(f - 1));
}

int poisson(double l){
	int k =0;
	double prob=0.0;
	double r = (double )rand()/RAND_MAX;
	do {
		prob+=pow(l,k)*exp(-l)/factorial(k);
		if (r<prob)
			return k;
		k++;
	
	} while(prob<0.9999);
	return k;
}

int HOD_alpha1(double M, double M1){
	if (M1==0.0)
		return 0;
	double xsat = M/M1;
	return poisson(xsat);
}

int countlines(char *fname) {
        char stemp[2005];
        int i=0;
        FILE *f;
        if ((f = fopen(fname, "r"))==NULL) {
                perror("fopen:");
                fprintf(stderr,"Error while reading the file %s\n",fname);
                return -1;
        }


        while(fgets(stemp,2000,f)!=NULL) {
		if (stemp[0]!='#')
                	i++;
        }
        fclose(f);
        return i;
}

double rand_gauss (void) {
  double v1,v2,s;

  do {
    v1 = 2.0 * ((double) rand()/RAND_MAX) - 1;
    v2 = 2.0 * ((double) rand()/RAND_MAX) - 1;

    s = v1*v1 + v2*v2;
  } while ( s >= 1.0 );

  if (s == 0.0)
    return 0.0;
  else
    return (v1*sqrt(-2.0 * log(s) / s));
}


double R_from_mass(double Mass,double OVD) {
        return  (double) pow((3./(4.*rho_crit*OVD*M_PI)*Mass),(1./3.));
}


double E2(double z){
	return (OMEGA_M * (1.0+z)*(1.0+z)*(1.0+z) + (1.0 - OMEGA_M));
}

double Omega(double z){
	return ( OMEGA_M * (1.0+z)* (1.0+z)* (1.0+z) ) / (OMEGA_M * (1.0+z)*(1.0+z)*(1.0+z) + (1.0 - OMEGA_M) );
}


double Delta_vir(double z){
	double d;
	d = 1.0-Omega(z);
	return 18*M_PI*M_PI + 82*d -39*d*d;
}

double I_NFW(double x){
	return ((1.0/(1.0+x)+log(1.0+x)-1.0));
}

double concentration_MICE(double M,double z){
	return (9.0/(1.0+z)*pow(M/2.0e13,-0.13));
}

double concentration_klypin(double M, double z){
	double C0, gamma, M0;
	if (z<0.25){
		C0=  9.5;	
		gamma=  0.09;
		M0=  3.0e5*1.0e12;
	}
	else if (z<0.75){
		C0=  6.75;	
		gamma=  0.088;
		M0=  5000*1.0e12;
	
	}
	else if(z<1.22){
		C0=  5.0;	
		gamma= 0.086;
		M0=  450*1.0e12;

	}	
	else{
		C0=  4.05;	
		gamma=  0.085;
		M0=  90.0*1.0e12;
	}

	return C0 * pow(M/1.0e12,-gamma)*(1.0+pow(M/M0,0.4));
}

void NFW_to_pos(double M,double *Dx,double *Dy,double *Dz,double z){
	double x_max = 50.0;
	double y_rand = ((double) rand()/RAND_MAX * I_NFW(x_max));
	double tol = 0.01*I_NFW(x_max);
	double low=0.0;
	double high = x_max;
	double mid=0.5*(low+high);	
	double y_try = I_NFW(mid); 
	double R, phi, costh;
	do{
		if (y_try>y_rand)
			high = mid;
		else 
			low = mid;
		mid = 	0.5*(low+high);
		y_try = I_NFW(mid);
  	}while(fabs(y_try - y_rand)<tol);
	
	R = mid*R_from_mass(M,Delta_vir(M))/concentration_klypin(M,z);
	phi = (double )rand()/RAND_MAX * 2*M_PI;
	costh = (double)rand()/RAND_MAX*2-1.0;

	*Dx = R * sqrt(1.0 - costh*costh) * cos(phi);
	*Dy = R * sqrt(1.0 - costh*costh) * sin(phi);
	*Dz = R * costh;
}
void vir_to_vel(double M,double z,double *Dvx,double *Dvy,double *Dvz){
	//double sigma = sqrt(0.2*Gnewton*M/M_to_R(M));
	double sigma;
	sigma = 476*0.9*pow( Delta_vir(z) * E2(z),1.0/6.0)*pow(M/1.0e15,1.0/3.0);
	
	*Dvx = sigma*rand_gauss();
	*Dvy = sigma*rand_gauss();
	*Dvz = sigma*rand_gauss();
}

double sigma_z(double z){
	        if (z<0.975)
                return (0.03 - 0.012*(z-0.6) + 1.23*(z-0.6)*(z-0.6)*(z-0.6)*(z-0.6));
        else
                return (0.0494 + (0.0494-0.0412)/0.5*(z-0.975));

}

double rand_photoz_polinomial(double z){
	return (z+(1.0+z)*rand_gauss()*sigma_z(z));
}


double corner(int c, double x){
        if (c==0)
                return x;
        else
                return LBOX-x;
}

*/
