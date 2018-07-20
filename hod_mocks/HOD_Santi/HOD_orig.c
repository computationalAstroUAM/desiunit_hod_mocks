#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <gsl/gsl_integration.h>
#include "chealpix.h"

#define rho_crit (27.755e10)

//#define Nbins 8
//double zcuts[Nbins+1]={0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
//double M1z[Nbins];
//={13.4, 13.5, 13.9, 13.9, 13.7, 0, 0, 0};

#define Nbins 8
double zcuts[Nbins+3]={0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2};
double sz[Nbins+2]={0.031,0.029,0.029,0.029,0.029,0.029,0.029,0.030,0.04,0.05};

double M1z_log[Nbins] ={13.4, 13.6, 14.2, 14.5, 14.0,-6.0,-6.0,-8.0};
double M1z[Nbins];

#define NPIX (12582912*4*4)
long nsidem=4096;

double Lbox=3072.0;

//#define mask_file ("input/Y1LSSmask_v1_il22_4096ring.dat") // V1.4.0
#define mask_file ("input/Y1LSSmask_v1_il22seeil4.04096ring.dat")

#define weight_file ("input/nbarY1redDESY1_1435redmean_z_bpz_short.dat")

#define Mth_file ("input/Mth_mean24.txt")
//#define outbase ("output/HALOGENlamps_V1.4.0_mock%04d.dat")
#define outbase ("output/HALOGENlamps_V1.4.1_mock%04d.dat")

#define inbase ("/bigdata/savila/BIG_CATS/Mocks/BigMD_Mocks/HALOGEN/output/HALOGEN_V1.3.0_lightcone_mock%03d_corner%d%d%d_zspec")

double rand_gauss (void);

#define Nw 40
double zwmin = 0.6;
double zwmax = 1.0;
double weight_array[Nw];

#define Nth 40
double zthmin = 0.6;
double zthmax = 1.0;
double Mth_array[Nw];

double LBOX = 3072.;
double OMEGA_M = 0.25;
double *z_array, *r_array, *y2_global;
double Binz = 0.1;
double z_max =2.5;
int Nz;
double rad_to_deg(double rad){
	return rad*180.0/M_PI;
}

double deg_to_rad(double dec){
	return dec*M_PI/180.0;
}
struct cosm_params { double Omega_m; double Omega_L; };
int which_bin(double z){
        if (z<zcuts[1])
                return 0;
        else if (z<zcuts[2])
                return 1;
        else if (z<zcuts[3])
                return 2;
        else if (z<zcuts[4])
                return 3;
        else if (z<zcuts[5])
                return 4;
        else if (z<zcuts[6])
                return 5;
        else if (z<zcuts[7])
                return 6;
        else if (z<zcuts[8])
                return 7;
        else if (z<zcuts[9])
                return 8;
        else
                return 9;
}


int which_M1_bin(double z){
        if (z<zcuts[1])
                return 0;
        else if (z<zcuts[2])
                return 1;
        else if (z<zcuts[3])
                return 2;
        else if (z<zcuts[4])
                return 3;
        else if (z<zcuts[5])
                return 4;
        else if (z<zcuts[6])
                return 5;
        else if (z<zcuts[7])
                return 6;
        else
                return 7;
}




//void read_photoz(char *);
int countlines(char *);
void NFW_to_pos(double,double *, double *,double *,double );
void vir_to_vel(double ,double,double *,double *,double *);
void box_to_lightcone(double ,double ,double ,double ,double ,double ,double *,double *,double *);
double rand_photoz(double zspec);
void invert_z_r();
double spline_inter(double *x, double *y, int N, double *y2, double xi);
int HOD_alpha1(double M, double M1);
void read_weights(char *);
void read_Mth(char *);
double * read_mask(char *filename);
double which_weight(double *ra, double *dec, double z, double M);

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
	if (argc!=5){
		fprintf(stderr,"usage: %s halos M1 ztrue galaxies\n",argv[0]);
		return -1;
	} 


	int seed;
//	double *mask;
	
	int j,Nhalos,Nsat;
	FILE *f_in,*f_out;
	//char maskfile[1024],input[1024],output[1024],line[1024];
	char input[1024],output[1024],line[1024];
	double M,x,y,z,vx,vy,vz;
	double Dx,Dy,Dz,Dvx,Dvy,Dvz;
	double xgal,ygal,zgal;
	double M1,ztrue;


	sprintf(input,argv[1]);
	M1 = atof(argv[2]);
	ztrue = atof(argv[3]);
	sprintf(output,argv[4]);	

	if (M1>0)
		M1 = pow(10.0,(double) M1 );

	
	seed = 42;
	srand(seed);

//	invert_z_r();
	

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

	
	#ifdef DEBUG 
	fprintf(stderr,"\nDEBUG mode ON **\n");
	Nhalos=100;
	#endif

	for (i=0;i<Nhalos;i++){
	
                sscanf(line,"%lf %lf %lf %lf %lf %lf %le",&x,&y,&z,&vx,&vy,&vz,&M);
		#ifdef DEBUG 
                fprintf(stderr,"\n- ihalo %d: %lf %lf %lf %le %lf %lf %lf %lf %lf %lf %lf %lf\n\t",i,ra,dec,zspec,M,zwrong,ztrue,x,y,z,vx,vy,vz);
		#endif
	#ifdef DEBUG 
		Dz= (double)(zthmax-zthmin)/Nth;
		int bin =(int) ((float)(photoz-zthmin)/Dz);
		fprintf(stderr,"z=%f, zbin=%d, Mth=%e ",photoz,bin,Mth_array[bin]);
	#endif
		if (M1>0){
			Nsat = HOD_alpha1(M,M1);
		} else if (M1==0.0){
			Nsat=0;
		}
		else{
			Nsat = 0;
			M = exp(log(M) + M1*rand_gauss());
		}
		//w = which_weight(&ra,&dec,photoz,M);

		fprintf(f_out,"%f %f %f %f %f %f %e %d\n",x,y,z,vx,vy,vz,M,Nsat);
		#ifdef DEBUG 
		fprintf(stderr,"\t%f %f %f %f %e %f %d\n",ra,dec,photoz,w,M,zspec,Nsat);
		#endif
		
		for (j=0;j<Nsat;j++){
			NFW_to_pos(M,&Dx,&Dy,&Dz,ztrue);
			vir_to_vel(M,ztrue,&Dvx,&Dvy,&Dvz);
			xgal=x+Dx;
			ygal=y+Dy;
			zgal=z+Dz;
			if (xgal>=Lbox)
				xgal-=Lbox;
			if (xgal<0)
				xgal+=Lbox;
			if (ygal>=Lbox)
				ygal-=Lbox;
			if (ygal<0)
				ygal+=Lbox;
			if (zgal>=Lbox)
				zgal-=Lbox;
			if (zgal<0)
				zgal+=Lbox;
			fprintf(f_out,"%f %f %f %f %f %f %e %d\n",xgal,ygal,zgal,vx+Dvx,vy+Dvy,vz+Dvz,M,-1);
		}
          	fgets(line,1000,f_in);
	}
	
	fclose(f_in);
	fclose(f_out);

	fprintf(stderr,"\nDone and writen to %s\n\n",output);

	return 0;

}

#ifdef NOMASK
double which_weight(double *ra, double *dec, double z, double M){
	return 1.0;
}
#else


double which_weight(double *ra, double *dec, double z, double M){
	long ipr;
	int bin;
	double Dz;
	//z-limits
	if (z<zwmin)
		return 0.0;	
	if (z>=zwmax)
		return 0.0;	
	#ifdef DEBUG 
	fprintf(stderr,"passed limits;  ");
	#endif

	//Mth
  #ifndef NOMTH
	Dz= (double)(zthmax-zthmin)/Nth;
	bin =(int) ((double) (z-zthmin)/Dz);
	#ifdef DEBUG 
	fprintf(stderr,"zbin=%d, Mth=%e ",bin,Mth_array[bin]);
	#endif


	if (M<Mth_array[bin])
		return 0.0;	
	#ifdef DEBUG 
	fprintf(stderr,"passed th; ");
	#endif
  #endif

	//mask
        if (*ra>180.0)
                (*ra) = (*ra) -360.0;
        (*dec) = -1*(*dec);

         if (sin(deg_to_rad(*dec))>-0.09 && (*ra)>30){
               (*ra)-=80.0;
               (*dec) = rad_to_deg( (double) asin( sin(deg_to_rad( (*dec) ) )+0.05));
         }
         else if ( sin( deg_to_rad( (*dec) ) )>-0.61 &&  sin( deg_to_rad( (*dec) ) )<-0.36  && (*ra)>15 ){
                (*ra) -= 80.0;
                (*dec) = rad_to_deg( (double) asin( sin( deg_to_rad(*dec)) -0.27 ) );
         }
         else
                (*ra) +=10.0;

	#ifdef DEBUG 
	fprintf(stderr,"patched: [%lf,%lf];  ",*ra,*dec);
	#endif	
	ang2pix_ring(nsidem,(M_PI*0.5-deg_to_rad(*dec)),deg_to_rad(*ra),&ipr);
	if (mask[ipr]==0.0){
		#ifdef DEBUG 
		fprintf(stderr,"out of mask; ");
		#endif
		return 0.0;
        } else if (((double) rand()/RAND_MAX)>mask[ipr]){
		#ifdef DEBUG 
		fprintf(stderr,"failed pix %f",mask[ipr]);
		#endif
		return 0.0;
	}
	#ifdef DEBUG 
	fprintf(stderr,"passed mask pix=%f; ",mask[ipr]);
	#endif
	//weight
	Dz= (double)(zwmax-zwmin)/Nw;
	bin =(int) ((double)(z-zwmin))/Dz;
	return weight_array[bin];	
}
#endif

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
	if (M1==1.0)
		return 0;
	double xsat = M/M1;
	return poisson(xsat);
}


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


/*
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
                //if (stemp[0]!='#')
                        i++;
        }
        fclose(f);
        return i;
}
*/

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
double rand_photoz(double z){
	return (z+(1.0+z)*rand_gauss()*sz[which_bin(z)]);
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

void box_to_lightcone(double x, double y,double z,double vx,double vy,double vz,double *ra,double *dec,double *zspec){
                double L=LBOX;
                double MaxX=LBOX;
                double MaxY=LBOX;
                double MaxZ=LBOX;
                double MinX=0;
                double MinY=0;
                double MinZ=0;

		double r,r2,ztrue,H,a,prod_v_r,Dx,Dy,Dz,D;

                if (x>MaxX)
                        x-=L;
                if (x<MinX)
                        x+=L;
                if (y>MaxY)
                        y-=L;
                if (y<MinY)
                        y+=L;
                if (z>MaxZ)
                       z-=L;
                if (z<MinZ)
                       z+=L;



		x = corner(CX,x);
		y = corner(CY,y);
		z = corner(CZ,z);			

                r2 = x*x + y*y + z*z;
                ztrue  =  spline_inter(r_array,z_array,Nz,y2_global,sqrt(r2));
                H = 100.0*sqrt(  OMEGA_M*(1+z)*(1+z)*(1+z) + (1-OMEGA_M) ); //[h km/s /Mpc]
                a = 1/(1+ztrue);
                prod_v_r = x*vx + y*vy + z*vz;

                Dx = 1.0/(a*H) * prod_v_r * x /r2;
                Dy = 1.0/(a*H) * prod_v_r * y /r2;
                Dz = 1.0/(a*H) * prod_v_r * z /r2;
                D=sqrt(Dx*Dx + Dy*Dy + Dz*Dz);

                x += Dx;
                if (x>MaxX)
                        x-=L;
                if (x<MinX)
                        x+=L;

                y += Dy;
                if (y>MaxY)
                        y-=L;
                if (y<MinY)
                        y+=L;

                z += Dz;
                if (z>MaxZ)
                       z-=L;
                if (z<MinZ)
                       z+=L;

                r = sqrt( (double) x*x + y*y + z*z );
                *zspec =  spline_inter(r_array,z_array,Nz,y2_global,r);
                *ra = rad_to_deg((double) atan2(y,x));
                *dec = rad_to_deg((double)asin(z/r));
}
double z_to_r_integrand(double x, void *p){
        struct cosm_params * params =  (struct cosm_params *) p;
        double Omega_m = (params->Omega_m);
        double Omega_L = (params->Omega_L);
        return  2997.92458/sqrt(Omega_m*(1+x)*(1+x)*(1+x)+Omega_L);
}



double z_to_r(double z){
        double r, epsrel=0.01,epsabs=0.01,err;
        size_t neval;

        gsl_function F;
        struct cosm_params params = {OMEGA_M, 1.0-OMEGA_M};
        F.function = &z_to_r_integrand;
        F.params = &params;

        gsl_integration_qng((const gsl_function *) &F, 0, z, epsabs, epsrel, &r, &err, &neval);
        //fprintf(stderr,"integral = %f +- %f (%d)\n",r,err,(int)neval);
        return r;
}


double spline_inter(double *x, double *y, int N, double *y2, double xi){

        double h,b,a;
        int i_low=0,i_high=N-1,i_mid;

        while(i_high-i_low>1){
                i_mid=(i_high+i_low)/2;
                if (x[i_mid]>xi)
                        i_high = i_mid;
                else
                        i_low = i_mid;
        }

        h = x[i_high] - x[i_low];
        a = (x[i_high]-xi)/h;
        b = (xi-x[i_low])/h;
        return (a*y[i_low]+b*y[i_high]+((a*a*a-a)*y2[i_low]+(b*b*b-b)*y2[i_high])*h*h/6.0);
}

void get_cubic_splines(double *x, double *y, int N, double *y2){
        int i;
        double *u,s,p;

        u = (double *) calloc(N,sizeof(double));

        y2[0]=0.0; //Natural spline
        u[0]=0.0;
        //fprintf(stderr,"%f\n",y2[0]);


        for (i=1;i<N-1;i++){
                s = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
                p = s * y2[i-1] + 2.0;
                y2[i] = (s - 1.0) / p;
                u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
                u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-s*u[i-1])/p;
                //fprintf(stderr,"%f %f %f\n",s,p,y2[i]);
        }

        y2[N-1]=0.0; //Natural spline

        for (i=N-2;i>=0;i--){
                y2[i] = y2[i]*y2[i+1]+u[i];
                //fprintf(stderr,"%f\n",y2[i]);
        }
        free(u);
}


void invert_z_r(){
        int i;
        Nz = (((int) (z_max/Binz) ) + 1);
        z_array = (double *) malloc(sizeof(double)*Nz);
        r_array = (double *) malloc(sizeof(double)*Nz);
        y2_global  = (double *) malloc(sizeof(double)*Nz);
        for (i=0;i<Nz;i++) {
                z_array[i] = 0.0 + Binz*i;
                r_array[i] = z_to_r(z_array[i]);
        }
        //Compute   z(r)   splines
        get_cubic_splines(r_array,z_array,Nz,y2_global);
}

double * read_mask(char *filename){
  double *mask;
  long N=0,id;
  double idf;
  FILE *f;
  char line[1000];
  double weight;
  double Area=0;
  mask = (double *) calloc(NPIX,sizeof(double));
  //N=countlines(filename);  

  f=fopen(filename,"r");
  //for(i=0;i<N;i++){
  while(fgets(line,1000,f)!=NULL){
                //sscanf(line,"%ld %lf",&id,&weight); V1.4.0
                sscanf(line,"%lf %lf",&idf,&weight); //V1.4.1
		id = (long) idf;
                if(id>NPIX){
                        fprintf(stderr,"ERROR: found pix=%ld\n>NPIX=%d at line %ld",id,NPIX,N);
                }

                mask[id]=weight;
                N++;
		Area+=weight;
  }
  fclose(f);
  fprintf(stderr,"Ngood_pix=%ld Area(int)=%f  Area(weight)=%f\n",N,(double)N/M_PI*4*180*180,Area/M_PI*4*180*180);
  return mask;
}
