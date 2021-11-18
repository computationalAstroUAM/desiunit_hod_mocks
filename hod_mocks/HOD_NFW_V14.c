#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
//#include <omp.h>
//#include <time.h>
//#include <gsl/gsl_integration.h>
//#include "chealpix.h"

#define rho_crit (27.755e10)

//V0.4 include vt
//V0.6 skewed Gaussian
//V0.7 Binomial
//V1.0 K constant. Got the I_NFW wrong!
//V1.3, re-fixed NFW (same I:NFW as V<1.0, xmax=1.0 vs xmax=50)
//V1.4, NFW truncated at rvir/K

#define MORE 1

#ifdef MORE
	//#define outbase ("/mnt/lustre/savila/HOD_NFW/output_V1/galaxies_1000Mpc_V%smore_NFW_mu%.3f_Ac%.4f_As%.5f_vfact%.2f_beta%.3f_K%.2f_vt%.0fpm%.0f_mock%d%d%d.dat")
        #define outbase ("../../DESI_outputs/output_V1/pruebas/galaxies_1000Mpc_V%smore_NFW_mu%.3f_Ac%.4f_As%.5f_vfact%.2f_beta%.3f_K%.2f_vt%.0fpm%.0f_BVG_product.dat") //BVG
#else
  	#define outbase ("../../DESI_outputs/output_V1/pruebas/galaxies_1000Mpc_V%s_NFW_mu%.3f_Ac%.4f_As%.5f_vfact%.2f_beta%.3f_K%.2f_vt%.0fpm%.0f_BVG_product.dat") //BVG
#endif
//#define inbase ("/mnt/lustre/eboss/OuterRim/OuterRim_sim/ascii/OuterRim_STEP266_z0.865/subvols27/OuterRim_STEP266_fofproperties_%d%d%d.txt")
#define inbase ("../../DESI_outputs/out_100p_X_Y_Z_VX_VY_VZ_logM.txt")  //BVG
//#define partfile ("/mnt/lustre/eboss/OuterRim/OuterRim_sim/ascii/vol500Mpc/OuterRim_STEP266_z0.865/OuterRim_STEP266_particles_500Mpc.txt")


float zsnap = 0.865;
//float zsnap = 0.9873; //BVG

double LBOX = 1000.;
//double OMEGA_M = 0.272;
double OMEGA_M = 0.3089; //UNITSIM

//double Mmin = 11.0;
//double M1 = 12.0;
//double alpha = 1.0;
//double siglogM = 0.25;
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
void NFW_to_pos(double,double *, double *,double *,double , double K);
void vir_to_vel(double ,double,double *,double *,double *);
//void box_to_lightcone(double ,double ,double ,double ,double ,double ,double *,double *,double *);
double rand_photoz(double zspec);
void invert_z_r();
double spline_inter(double *x, double *y, int N, double *y2, double xi);
int HOD_erf(double logM, double mu, double sig, double As);
int HOD_gauss(double logM, double mu, double sig, double As);
int HOD_gaussPL(double logM, double mu, double sig, double Ac,double gamma);
int HOD_powerlaw(double M, double M0, double M1, double alpha, double As);

int CX,CY,CZ;

double *mask;

//void initialise_HMF();
//double HMF1_to_HMF2(double );
//double watson(double );

double BETA;
          //BVG

int main(int argc, char **argv){
  	fprintf(stderr,"Run call:\n\t");
	int i ;
   	for (i=0;i<argc;i++){
			fprintf(stderr,"%s  ",argv[i]);
   	}
	fprintf(stderr,"\n");
	fprintf(stderr,"Number of params: %d \n",argc);
	if (argc!=9){      //BVG
		fprintf(stderr,"usage: %s mu Ac As vfact beta K  vt sigma(vt) \n",argv[0]); //BVG
		return -1;
	}


	int seed;
//	double *mask;
	
	int j,Nhalos,Nsat,Ncent;
	FILE *f_in,*f_out;
	char output[1024],line[1024], input[1024];
	double M,x,y,z,vx,vy,vz;
	double Dx,Dy,Dz,Dvx,Dvy,Dvz;
	double xgal,ygal,zgal;
	double M0,M1,alpha,As;
	double mu,sig,Ac;

	mu = atof(argv[1]);
	Ac = atof(argv[2]);
	As = atof(argv[3]);
	double vfact = atof(argv[4]);   //alpha_v
	BETA  = atof(argv[5]); //Poisson = 0 ,  nearest int = -2
	double K = atof(argv[6]);
	double vt = atof(argv[7]);  //=0  =500
	double vtdisp = atof(argv[8]);  //=0  =200
	//BVG   Binomial, only if BETA = 0, BETA2 > 0. If noy using binomial, then set BETA2 = 0
	double ux,uy,uz,vtrand,Dr;

	//int imock = atoi(argv[9]);  
	//int jmock = atoi(argv[10]);
	//int kmock = atoi(argv[11]);
	#ifdef HOD1
	//Old Version:
	//char *version="1.41";
        //M0 = 0.;
        char *version="1.41b";
	M0 = pow(10.,mu);
        M1 = pow(10.,mu + 1.3);
        alpha  = 1.0;
	sig = 0.15;
	#else 
		#ifdef HOD2 
		char *version="1.42";
		M0 = pow(10.,mu - 0.1);
	        M1 = pow(10.,mu + 0.3);
	        alpha  = 0.8;
	        sig = 0.12;
		#else
			char *version="1.4";
        		M0 = pow(10.,mu - 0.05);
        		M1 = pow(10.,mu + 0.35);
        		alpha  = 0.9;
        		sig = 0.08;
        		double gamma= -1.4;
		#endif
	#endif
	//sprintf(output,outbase,version,mu,Ac,As,vfact,BETA,K,vt,vtdisp,imock,jmock,kmock);
	//sprintf(input,inbase,imock,jmock,kmock);
	sprintf(output,outbase,version,mu,Ac,As,vfact,BETA,K,vt,vtdisp);    //BVG
	sprintf(input,inbase);

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
	double logM;
	long unsigned id;
	for (i=0;i<Nhalos;i++){	
                sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lu", &x,&y,&z, &vx,&vy,&vz, &logM, &id);
		M = pow(10.,logM);
		Nsat = HOD_powerlaw(M, M0, M1, alpha, As);
		#ifdef HOD1
		Ncent= HOD_erf(logM,mu,sig,Ac);
		#else
			#ifdef HOD2
			Ncent = HOD_gauss(logM,mu,sig,Ac);
			#else
			Ncent = HOD_gaussPL(logM,mu,sig,Ac,gamma);
			#endif
		#endif
		if (Ncent == 1)
			#ifdef MORE
			fprintf(f_out,"%f %f %f %f %f %f %e %d %f %f %f %f %f %f %lu\n",x,y,z,vx,vy,vz,M,Nsat,0.,0.,0.,0.,0.,0.,id);
			#else
			fprintf(f_out,"%f %f %f %f %f %f %e %d\n",x,y,z,vx,vy,vz,M,Nsat);
			#endif
		
	
		for (j=0;j<Nsat;j++){
				NFW_to_pos(M,&Dx,&Dy,&Dz,zsnap,K);
				vir_to_vel(M,zsnap,&Dvx,&Dvy,&Dvz);
				vtrand = rand_gauss()*vtdisp+vt;
				Dr=sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
				ux = -Dx/Dr;
				uy = -Dy/Dr;
				uz = -Dz/Dr;
				Dvx = (vfact*Dvx + ux*vtrand);  //Eliminado el factor (1+zsnap) en esta linea y las dos siguientes
				Dvy = (vfact*Dvy + uy*vtrand);
				Dvz = (vfact*Dvz + uz*vtrand);

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
				#ifdef MORE
				fprintf(f_out,"%f %f %f %f %f %f %e %d %f %f %f %f %f %f %lu\n",xgal,ygal,zgal,vx+Dvx,vy+Dvy,vz+Dvz,M,Nsat, Dvx,Dvy,Dvz,Dx,Dy,Dz,id); //BVG
				#else
				fprintf(f_out,"%f %f %f %f %f %f %e %d\n",xgal,ygal,zgal,vx+Dvx,vy+Dvy,vz+Dvz,M, Nsat);  //BVG
				#endif
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
	double r = (double) rand()/(RAND_MAX+1.0);
	do {
		prob+=pow(l,k)*exp(-l)/factorial(k);
		if (r<prob)
			return k;
		k++;
	
	} while(prob<0.999999999999999999999);
	return k;
}

/*
int poisson(double l){
	int n =(int) l;
	double rest = l-n;
	if( ( (double) rand()/(RAND_MAX+1.))<rest)
		return n+1;
	else
		return n;
}
*/






int next_integer(float x){
	int low = (int) floor(x);
	float rand01 = ((float) rand()/(RAND_MAX+1.));

	if (rand01>(x-low))
		return low;
	else
		return low+1;

}



/*int neg_binomial(double x, double beta){

	double sigma = sqrt(x)*(1.0+beta);
	double p = x/(sigma*sigma);
	double r = x*x/(sigma*sigma-x);

	double P=0;
	int N=-1;
	double rand01 =  ((float) rand()/(RAND_MAX+1.));
	do{
		N++;
		P+= tgamma(N+r)/( tgamma(r)*tgamma(N+1) ) * pow(p,r)*pow(1-p,N) ;	

	} while(P<rand01);

	return N;
}*/

/*BVG code for negative binomial and binomial. x is the mean of the distribution, and beta the parameter that increases the variance of the distribution*/


/*overflow happens in gamma(x) when x>171.7*/

double product(double a,double b){ /*LONG DOUBLE BAD RESULTS*/
        int c = (round)(a-b);
        double s = 1.0;
        int j;
        if (c%1 == 0)
                for(j=0;j<c+1;j++){
                        s *= (j+b);
                }
	else
                printf("Error: no integer number c");
	return s;


}

/*
double P_NB3(int C,double D, double E){
        double a;
        a = product(C+D-1,D)/( tgamma(C+1) ) * pow(E,D)*pow(1-E,C) ;
        return a;

}
*/

int neg_binomial(double x,double beta){   /*the definition of r can be extended to reals*/

        double r = x/beta;
	double p = r/(r+x);
/*	double xmax = 171.7;*/
        double P=0;
        int N=-1;
        double rand01 =  ((float) rand()/(RAND_MAX+1.));
        do{
                N++;
		P+= product(N+r-1,r)/( tgamma(N+1) ) * pow(p,r)*pow(1-p,N) ;
/*                P+= tgamma(N+r)/( tgamma(r)*tgamma(N+1) ) * pow(p,r)*pow(1-p,N) ; */
        } while(P<rand01);

        return N;
}


int n(double y, double z){  //BVG
        int outn = 0;
	if (1/z >= trunc(y+1.0))
                outn = 1/z;
        else
                outn = trunc(y+1.0);
        return outn;
}

double p(double y, double z){   //BVG
        if (1/z >= trunc(y+1.0))
                return y*z;
        else
                return y/(trunc(y+1.0));

}

int binomial(double x, double beta){                                //BVG

        double a = -beta;
        double P=0;
        int N=-1;
        double rand01 =  ((float) rand()/(RAND_MAX+1.));
        do{
                N++;
                P+= tgamma(n(x,a)+1)/( tgamma(n(x,a)+1-N)*tgamma(N+1) ) * pow(p(x,a),N)*pow(1-p(x,a),n(x,a)-N) ;

        } while(P<rand01);

        return N;
}



int HOD_powerlaw(double M, double M0, double M1, double alpha, double As){
	if (M1<=0.0)
		return 0;
	if (M<M0)
		return 0;
	double xsat = (M-M0)/M1;
	if (BETA<-1.0)
		return next_integer((As*pow(xsat, alpha)));
	if (BETA==0.)                                                 //BVG
		return poisson(As*pow(xsat, alpha));                  //BVG
	if (BETA<0. && BETA>=-1.0)                                     //BVG
	        return binomial(As*pow(xsat, alpha),BETA);            //BVG
	if (BETA>0.)                                                  //BVG
		return neg_binomial(As*pow(xsat, alpha),BETA);        //BVG
}
int HOD_erf(double logM, double mu, double sig, double As){
        double r = As*0.5*(1.0+erf((logM-mu)/sig));
        if (  (((double) (rand()*(RAND_MAX+1.0)+rand())/((RAND_MAX+1.0)*(RAND_MAX+1.0)-1)))<r )
                return 1;
        return 0;
}

int HOD_gauss(double logM, double mu, double sig, double As){

        double r = As/(sig*sqrt(2.0*M_PI))*exp(-(logM-mu)*(logM-mu)/(2*sig*sig));
	if (  (((double) (rand()*(RAND_MAX+1.0)+rand())/((RAND_MAX+1.0)*(RAND_MAX+1.0)-1)))<r )
		return 1;
	return 0;
}
int HOD_gaussPL(double logM, double mu, double sig, double Ac,double gamma){
        double r;
        if (logM<mu)
                r = Ac/(sig*sqrt(2.0*M_PI))*exp(-(logM-mu)*(logM-mu)/(2*sig*sig));
        else
                r = Ac/(sig*sqrt(2.0*M_PI))*pow(10.0,gamma*(logM-mu));

	if (r>1.0){
		fprintf(stderr,"\n\n ## ##  ##  ##  ## ## \n");
		fprintf(stderr,"\n\n******************************\n");
		fprintf(stderr,"\n*  WARNING: <Ncen>=%.2f >1     *\n",r);
		fprintf(stderr,"\n******************************\n\n");
	}
        if (  (((double) (rand()*(RAND_MAX+1.0)+rand())/((RAND_MAX+1.0)*(RAND_MAX+1.0)-1)))<r )
                return 1;

        return 0;
}

/*
void read_Mth(char *fname){
        int i;
        FILE *f;
        char line[1024];
        float ftemp,z;
        if ((f = fopen(fname, "r"))==NULL) {
                perror("fopen:");
                fprintf(stderr,"Error while reading the file %s\n",fname);
                exit(0);
        }
        for (i=0;i<Nth;i++){
                fgets(line,1000,f);
                sscanf(line,"%f %f",&z,&ftemp);
                Mth_array[i]=ftemp;
        }
        fclose(f);
}


void read_weights(char *fname){
	int i;
	FILE *f;
	char line[1024];
	double ftemp,z;
        if ((f = fopen(fname, "r"))==NULL) {
                perror("fopen:");
                fprintf(stderr,"Error while reading the file %s\n",fname);
                exit(0);
        }
	for (i=0;i<Nw;i++){
		fgets(line,1000,f);
                sscanf(line,"%lf %lf",&z,&ftemp);
		weight_array[i]=ftemp;
	}
	fclose(f);
}
*/

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
/*
*/


//Version V1.0 is gone had 1/(ax^3+bx^2+cx) ish shape

double I_NFW(double x){
	return ((1.0/(1.0+x)+log(1.0+x)-1.0));
}


double concentration_MICE(double M,double z){
	return (9.0/(1.0+z)*pow(M/2.0e13,-0.13));
}
// WMAP7 cosmology, virial overdensity
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

void NFW_to_pos(double M,double *Dx,double *Dy,double *Dz,double z, double k){
	double x_max = concentration_klypin(M,z);
	// this was in olds versions (V<1.0), but NFW diverge! double x_max = 50.0;
	double y_rand = ((double) rand()/RAND_MAX * I_NFW(x_max));
	double tol = 0.001*I_NFW(x_max);
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
	
	R = mid*R_from_mass(M,Delta_vir(M))/(concentration_klypin(M,z)*k);
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






