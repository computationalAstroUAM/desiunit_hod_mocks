#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>



int countlines(char *);
void write_catalog(char *,int,double*,double*,double*);
void read_catalog(char *,int ,double**,double**, double**,double**,int);

int main(int argc, char **argv){

        if (argc!=6){
                fprintf(stderr,"usage: %s catalog Lbox z Omega_m output\n",argv[0]);
                return -1;
        }

        char catalog[200], cat_out[200];
        double Lbox,redshift,Om,h;
        int Vcomp,Nhalos;

        strcpy(catalog,argv[1]);
        Lbox=atof(argv[2]);
        redshift=atof(argv[3]);
        Om=atof(argv[4]);
	strcpy(cat_out,argv[5]);

        h=1.0;
        Vcomp=1;

        fprintf(stderr,"Welcome to Redshift convert\n");
        fprintf(stderr,"Parameters: L=%f, z=%f, Om=%f",Lbox,redshift,Om);

        Vcomp=1;

        fprintf(stderr,"Input file: %s\n",catalog);
        fprintf(stderr,"Output file: %s\n",cat_out);

        Nhalos = countlines(catalog);

        fprintf(stderr,"Number of halos: %d\n",Nhalos);

        double *x,*y,*z,*v;
        x = (double *) calloc(Nhalos,sizeof(double));
        if (x==NULL){
                fprintf(stderr,"Not enough memory for x\n");
                exit(0);
        }
        y = (double *) calloc(Nhalos,sizeof(double));
        if (y==NULL){
                fprintf(stderr,"Not enough memory for y\n");
                exit(0);
        }
        z = (double *) calloc(Nhalos,sizeof(double));
        if (z==NULL){
                fprintf(stderr,"Not enough memory for z\n");
                exit(0);
        }
        v = (double *) calloc(Nhalos,sizeof(double));
        if (v==NULL){
                fprintf(stderr,"Not enough memory for z\n");
                exit(0);
        }

        fprintf(stderr,"Allocated\n");

        read_catalog(catalog,Nhalos,&x,&y,&z,&v,Vcomp);

        fprintf(stderr,"Catalog read!\n");
        fprintf(stderr,"x[0]=%f, Y[0]=%f, Z[0]=%f v[0]=%f\n",x[0],y[0],z[0],v[0]);
        fprintf(stderr,"x[1]=%f, Y[1]=%f, Z[1]=%f v[1]=%f\n",x[1],y[1],z[1],v[1]);

        fprintf(stderr,"         ...\n");
        fprintf(stderr,"x[last]=%f, Y[last]=%f, Z[last]=%f v[last]=%f\n",x[Nhalos-1],y[Nhalos-1],z[Nhalos-1],v[Nhalos-1]);

        double H=100*sqrt( Om*pow(1+redshift,3.0) + (1-Om) );
        double fact=1.0/H*(1+redshift);
        int i;
        fprintf(stderr,"  fact=%f\n",fact);

        for(i=0;i<Nhalos;i++){
                if (Vcomp==1){
                        x[i]=x[i]+fact*v[i];
                        if (x[i]>Lbox)
                                x[i]=x[i]-Lbox;
                        if (x[i]<0)
                                x[i]=x[i]+Lbox;
                }
                if (Vcomp==2){
                        y[i]=y[i]+fact*v[i];
                        if (y[i]>Lbox)
                                y[i]=y[i]-Lbox;
                        if (y[i]<0)
                                y[i]=y[i]+Lbox;
                }
                if (Vcomp==3){
                        z[i]=z[i]+fact*v[i];
                        if (z[i]>Lbox)
                                z[i]=z[i]-Lbox;
                        if (z[i]<0)
                                z[i]=z[i]+Lbox;
                }
                if(i<10){
                        fprintf(stderr,"i=%d, Delta(r)=%f\n",i,fact*v[i]);
                }
                if (i%1000000==0)
                        fprintf(stderr,"%d million done!\n",i/1000000);
        }
        fprintf(stderr,"Transformation done!\n");

        write_catalog(cat_out,Nhalos,x,y,z);
        fprintf(stderr,"Catalog written, Bye,Bye!\n");
        return 0;
}

void write_catalog(char *filename,int Nhalos,double *x,double *y,double *z){
        FILE *f;
        int i;
        if ((f=fopen(filename,"w") )== NULL){
                fprintf(stderr,"Couldnt open output file %s\n",filename);
                exit(0);
        }
        for(i=0;i<Nhalos;i++){
                fprintf(f,"%f %f %f\n",x[i],y[i],z[i]);
                if (i%1000000==0)
                        fprintf(stderr,"%d million done!\n",i/1000000);
        }
        fclose(f);

}
void read_catalog(char *catalog,int Nhalos,double **x,double **y, double **z,double **v,int Vcomp){
        char line[1005];
        int i=0;
        FILE *file;
        double a,b,c,d,e,f;


        fprintf(stderr,"Vcomp=%d\n",Vcomp);
        if ((file = fopen(catalog,"r"))==NULL){
                fprintf(stderr,"Could not open input file %s\n", catalog);
                exit(0);
        }

        for(i=0;i<Nhalos;i++){
                fgets(line,1000,file);
                sscanf(line,"%lf %lf %lf %lf %lf %lf",&a,&b,&c,&d,&e,&f);
                (*x)[i]=a;
                (*y)[i]=b;
                (*z)[i]=c;
                if (Vcomp==1)
                        (*v)[i]=d;
                if (Vcomp==2)
                        (*v)[i]=e;
                if (Vcomp==3)
                        (*v)[i]=f;
                if (i%1000000==0){
                        fprintf(stderr,"%d million done!\n",i/1000000);
                        fprintf(stderr,"x,y,x,v = %f,%f,%f,%f\n",a,b,c,d);
                        fprintf(stderr,"x,y,x,v = %f,%f,%f,%f\n",(*x)[i],(*y)[i],(*z)[i],(*v)[i]);
                }

        }
        fclose(file);
}


int countlines(char *fname) {
        char stemp[200];
        int i=0;
        FILE *f;

        if ((f = fopen(fname, "r"))==NULL) {
                perror("fopen:");
                fprintf(stderr,"Error while reading the file %s\n",fname);
                return -1;
        }


        while(fgets(stemp,200,f)!=NULL) {
                i++;
        }
        fclose(f);
        return i;
}

