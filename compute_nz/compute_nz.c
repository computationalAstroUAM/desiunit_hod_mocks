#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>


float Omm;

#define Nz 70
double zmin = 0.5;
double zmax = 1.2;

typedef struct {
        float ra;
        float dec;
        float zphoto;
        float zspec;
        float M;
} particle;

int read_catalog(char *catalog, particle **P);
double z_to_r(double z);

int main(int argc, char **argv){
        int i;
        fprintf(stderr,"Run call:\n\t");
        for (i=0;i<argc;i++){
                        fprintf(stderr,"%s  ",argv[i]);
        }
        fprintf(stderr,"\n");
        if (argc!=5){
                fprintf(stderr,"usage: %s mock out fsky Omm\n",argv[0]);
                return -1;
        }
        char catalog[1024],output[1024];
        int Nhalos,ibin;
        int *histogram;
        double z,Dz,Area,V,zup,zdown,fsky;
        FILE *f;
        particle *Parts;

        sprintf(catalog,"%s",argv[1]);
        sprintf(output,"%s",argv[2]);
        fsky = atof(argv[3]);
        Omm = atof(argv[4]);


        Nhalos = read_catalog(argv[1],&Parts);
        fprintf(stderr,"Read %d halos in %s\n",Nhalos,argv[1]);

        histogram = (int *) calloc(Nz,sizeof(int));
        Dz=(double)(zmax-zmin)/(Nz);

        for(i=0;i<Nhalos;i++){
                ibin =  (int) floor((float) (Parts[i].zphoto-zmin)/Dz);
                if (ibin>=0 && ibin<Nz)
                        histogram[ibin]++;
                //if(i<30)
                  //      fprintf(stderr,"i=%d, z=%f, ibin=%d, hist=%d\n",i,Parts[i].zphoto,ibin,histogram[ibin]);
        }

        if ((f=fopen(output,"w") )== NULL){
                fprintf(stderr,"Couldnt open output file %s\n",output);
                exit(0);
        }

        fprintf(f,"#z N dN/dz dN/dA/dz dN/dV #used fsky=%f  Omega_matter=%f\n",fsky,Omm);
        for(ibin=0;ibin<Nz;ibin++){
                z = zmin + (ibin+0.5)*Dz;
                zup = zmin + (ibin+1.0)*Dz;
                zdown = zmin + (ibin)*Dz;
                V= fsky*(pow(z_to_r(zup),3.0) - pow(z_to_r(zdown),3.0))*4.0*M_PI/3.0;
                //fprintf(stderr,"z=%f, V=%e, N=%d, n=%f\n",z,V,histogram[ibin],histogram[ibin]/(360.0*360.0/M_PI*fsky*Dz));
                fprintf(f,"%f %d %f %e %e\n",z,histogram[ibin],histogram[ibin]/Dz,histogram[ibin]/(360.0*360.0/M_PI*fsky*Dz),histogram[ibin]/V);
        }
        fclose(f);
        return 0;
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



struct cosm_params { double Omega_m; double Omega_L; };
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
        struct cosm_params params = {Omm, 1.0-Omm};
        F.function = &z_to_r_integrand;
        F.params = &params;

        gsl_integration_qng((const gsl_function *) &F, 0, z, epsabs, epsrel, &r, &err, &neval);
        return r;
}


int read_catalog(char *catalog, particle **P){
        char line[1005];
        int i=0,Nhalos;
        FILE *file;
        float a,b,c,d,e;
        fprintf(stderr,"read %s\n",catalog);

        Nhalos=countlines(catalog);
        *P = (particle *) malloc(sizeof(particle)*(Nhalos));

        if ((file = fopen(catalog,"r"))==NULL){
                fprintf(stderr,"Could not open input file %s\n", catalog);
                exit(0);
        }

        for(i=0;i<(Nhalos);i++){
                fgets(line,1000,file);
                if (line[0]!='#'){
                        sscanf(line,"%f %f %f %f %f",&a,&b,&c,&d,&e);
                        (*P)[i].ra =a;
                        (*P)[i].dec =b;
                        (*P)[i].zphoto =c;
                        (*P)[i].M=d;
                        (*P)[i].zspec =e;
                }

        }
        fclose(file);
        return Nhalos;
}

