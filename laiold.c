#define NNN 100

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"
#include "nrutil.h"

#ifdef IBM_PC
#include <alloc.h>
#endif

/**************************/
/*                        */
/*  Measurment Routines   */
/*                        */
/**************************/


int qmax;//str factor array size
int *distrib, *distrib_brpt, *distrib_ci_n, *distrib_ci_p;
int *distrib_all, *distrib_all_ch;
int *distrib_pe, *distrib_pe_ch, *distrib_pe_d_all, *distrib_pe_d_ch;
int *distrib_cmdend_linecm, *distrib_cmall_cmdend;
int *distrib_dend, *distrib_dend_r0cm, *distrib_dend_ch;
double *hist_eta_loc;

double *rg, *g1, *g2, *g3, *g4, *g5, *cor, *energy, *energy_cl;
double *x_end, *x_mid, *x_cm, *r_ee, *r2_ee, *l2_;
double *y_end, *y_mid, *y_cm, a2m, dnpm, dnlpm, count;
double *z_end, *z_mid, *z_cm, am2, dlpm, dllpm, dallmono;
double *str, *strav1d, *stravalld, *strav1p, *stravallp,  *stravall, *strav2;
double *order_parameters, *order_parameters_track;
int *chdb;
//float diag_store[4];

#define square(x) ((x)*(x))

int get_full_energy()//хз, правильно ли работает...
{
    int i;
    double energ=0;
    for(i=0; i<Nmono; i++)
    {
        if(type[i]==1) energ += en_of_groups(i, x[i], y[i], z[i]);
        energ += get_energy(x[i], y[i], z[i]);
    }

    return (int)(((int)energ)%100000000000)/10000000;
}

void tred2(float *a, int n, float *d, float *e)
{
    int l, k, j, i;
    float scale, hh, h, g, f;

    for (i=n; i>=2; i--)
    {
        l=i-1;
        h=scale=0.0;
        if (l>1)
        {
            for (k=1; k<=l; k++) scale += fabs(*(a+k*4+i));
            if (scale==0.0) e[i]=*(a+l*4+i);
            else
            {
                for (k=1; k<=l; k++)
                {
                    *(a+k*4+i) /= scale;
                    h += *(a+k*4+i)**(a+k*4+i);
                }
                f=*(a+4*l+i);
                g=(f >=0.0 ? -sqrt(h) : sqrt(h));
                e[i]=scale*g;
                h -=f*g;
                *(a+4*l+i)=f-g;
                f=0.0;
                for (j=1; j<=l; j++)
                {
                    *(a+4*i+j)=*(a+4*j+i)/h;
                    g=0.0;
                    for (k=1; k<=j; k++)	g += *(a+4*k+j)**(a+4*k+i);
                    for (k=j+1; k<=l; k++)	g += *(a+4*j+k)**(a+4*k+i);
                    e[j]=g/h;
                    f += e[j]**(a+4*j+i);
                }
                hh=f/(h+h);
                for (j=1; j<=l; j++)
                {
                    f=*(a+4*j+i);
                    e[j]=g=e[j]-hh*f;
                    for (k=1; k<=j; k++) *(a+4*k+j) -= (f*e[k]+g**(a+4*k+i));
                }
            }
        }
        else	e[i]=*(a+4*l+i);
        d[i]=h;
    }
    d[1]=0.0;
    e[1]=0.0;
    for (i=1; i<=n; i++)
    {
        l=i-1;
        if (d[i])
        {
            for (j=1; j<=l; j++)
            {
                g=0.0;
                for (k=1; k<=l; k++) g +=*(a+4*k+i)**(a+4*j+k);
                for (k=1; k<=l; k++) *(a+4*j+k) -= g**(a+4*i+k);
            }
        }
        d[i]=*(a+4*i+i);
        *(a+4*i+i)=1.0;
        for (j=1; j<=l; j++) *(a+4*i+j)=*(a+4*j+i)=0.0;
    }
}

void tqli(float *d, float *e, int n, float *z)
{
    float pythag(float a, float b);
    int m, l, iter, i, k;
    float s, r, p, g, f, dd, c, b;

    for (i=2; i<=n; i++) e[i-1]=e[i];
    e[n]=0.0;
    for (l=1; l<=n; l++)
    {
        iter=0;
        do
        {
            for (m=l; m<=n-1; m++)
            {
                dd=fabs(d[m])+fabs(d[m+1]);
                if ((float)(fabs(e[m])+dd) == dd) break;
            }
            if (m !=l)
            {
                if (iter++ ==100) nrerror("Too many iterations in tqli");
                g=(d[l+1]-d[l])/(2.0*e[l]);
                r=pythag(g,1.0);
                g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
                s=c=1.0;
                p=0.0;
                for (i=m-1; i>=l; i--)
                {
                    f=s*e[i];
                    b=c*e[i];
                    e[i+1]=(r=pythag(f,g));
                    if (r == 0.0)
                    {
                        d[i+1] -= p;
                        e[m]=0.0;
                        break;
                    }
                    s=f/r;
                    c=g/r;
                    g=d[i+1]-p;
                    r=(d[i]-g)*s+2.0*c*b;
                    d[i+1]=g+(p=s*r);
                    g=c*r-b;
                    for (k=1; k<=n; k++)
                    {
                        f=*(z+4*k+i+1);
                        *(z+4*k+i+1)=s**(z+4*k+i)+c*f;
                        *(z+4*k+i)=c**(z+4*k+i)-s*f;
                    }
                }
                if (r==0.0 && i>=l) continue;
                d[l] -= p;
                e[l]=g;
                e[m]=0.0;
            }
        }
        while (m!=l);
    }
}


float pythag(float a, float b)
{
    float absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

long  put_chdb( char *fname )
{
    int i,j;
    FILE *fp;
    fp=fopen(fname,"w");
    if(fp==(NULL)) return(-1);
    for (i=0; i<maxtime; i++)
    {
        for (j=0; j<ngen; j++)
        {
            fprintf(fp,"%d ",chdb[i*ngen+j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    return(0);
}


long  put_IT( char *fname )
{
    int i;
    double dd;
    FILE *fp;
    dd=dlpm*dnpm*di*di/(double)(nrun);
    fp=fopen(fname,"w");
    if(fp==(NULL)) return(-1);
    for (i=0; i<maxtime; i++)
        fprintf(fp,"%lf %lf %lf %lf\n",diag_store[4*i]*dd, diag_store[4*i+1]*dd, diag_store[4*i+2]*dd, diag_store[4*i+3]*dd);
    fclose(fp);
    return(0);
}


long  put_distrib( char *fname )
{
    double dens;
    int i, g[10];
    FILE *fp;
    fp=fopen(fname,"w");
    if(fp==(NULL)) return(-1);

//	fprintf(fp,"time\tcounterN\tcounterN\tcounterP\tcounterP\tdistrMU\tdensDe\tcmDe-CMU\tdenLnLn\tdenChLnLn\tdenLnDe\tdenChLnDe\tcmDe-cmLn\tcmDe-cmSys\tdensAll\tdensChAll\tBrPt\tBrPt0\tBrPt1\tBrPt2\tBrPt3\tBrPt4\tBrPt5\tBrPt6\tBrPt7\tBrPt8\tBrPt9\n");
    fprintf(fp,"x\tcounterN\tcounterN\tcounterP\tcounterP\tdistrMU\tdensDe\tcmDe-CMU\tdenChDe\tdenLnLn\tdenChLnLn\tdenLnDe\tdenChLnDe\tcmDe-cmLn\tcmDe-cmSys\tdensAll\tdensChAll\tBrPt\tBrPt0\tBrPt1\tBrPt2\tBrPt3\tBrPt4\tBrPt5\tBrPt6\tBrPt7\tBrPt8\tBrPt9\n");

    for (i=0; i<10; i++)
    {
        if (i>=ngen) g[i]=ngen-1;
        else g[i]=i;
    }
    for (i=0; i<4*sysi; i++)
    {
        //4.189=4*Pi/3
        fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                (double)i/4,
                (double)distrib_ci_n[(int)i]/(maxtime*(Nall-Nmono)),//distrib counterions
                (double)distrib_ci_n[(int)i]*pow(0.7,3)/(maxtime*(Nall-Nmono)*(pow((double)(i+1)/4,3)-pow((double)i/4,3))),//density counterions
                (double)distrib_ci_p[(int)i]/(maxtime*(Nall-Nmono)),//distrib counterions
                (double)distrib_ci_p[(int)i]*pow(0.7,3)/(maxtime*(Nall-Nmono)*(pow((double)(i+1)/4,3)-pow((double)i/4,3))),//density counterions
                (double)distrib_dend[(int)i]/(maxtime*lpm*npm),//distrib monomer units
                (double)distrib_dend[(int)i]*pow(0.7,3)/(maxtime*lpm*npm*(pow((double)(i+1)/4,3)-pow((double)i/4,3))),//density of dendrimer

                (double)distrib_dend_r0cm[(int)i],//dist center mass of dendrimer and central monomer

                (double)distrib_dend_ch[(int)i]*pow(0.7,3)/(maxtime*Nch*(pow((double)(i+1)/4,3)-pow((double)i/4,3))),//density of charged units of dendrimer
                (double)distrib_pe[(int)i]*pow(0.7,3)/(maxtime*llpm*(pow((double)(i+1)/4,3)-pow((double)i/4,3))),//density of polyelectrolyte (TO CENTER OF MASS OF polyelectrolyte)
                (double)distrib_pe_ch[(int)i]*pow(0.7,3)/(maxtime*Nch*(pow((double)(i+1)/4,3)-pow((double)i/4,3))),//density of charged units of polyelectrolyte (TO CENTER OF MASS OF polyelectrolyte)
                //10
                (double)distrib_pe_d_all[(int)i]*pow(0.7,3)/(maxtime*llpm*(pow((double)(i+1)/4,3)-pow((double)i/4,3))),//density of polyelectrolyte (TO CENTER OF MASS OF DENDRIMER)
                (double)distrib_pe_d_ch[(int)i]*pow(0.7,3)/(maxtime*Nch*(pow((double)(i+1)/4,3)-pow((double)i/4,3))),//density of charged units of polyelectrolyte (TO CENTER OF MASS OF DENDRIMER)
                (double)distrib_cmdend_linecm[(int)i],//hist cm of dendrimer and cm oflinear
                (double)distrib_cmall_cmdend[(int)i],//hist cm of dendrimer and cm of system
                (double)distrib_all[(int)i]*pow(0.7,3)/(maxtime*(npm*lpm+nlpm*llpm)*(pow((double)(i+1)/4,3)-pow((double)i/4,3))),//density of all monomer units
                (double)distrib_all_ch[(int)i]*pow(0.7,3)/(maxtime*(Nch+Nch)*(pow((double)(i+1)/4,3)-pow((double)i/4,3))),//density of all charged monomer units
                (double)distrib_brpt[(int)i]*2/(maxtime*Nmono),//distrib branch points
                (double)distrib[(int)(g[0]*4*sysi+i)]/(maxtime*pow(2,g[0]+1)),//distrib br pt of generation 0
                (double)distrib[(int)(g[1]*4*sysi+i)]/(maxtime*pow(2,g[1]+1)),//etc.
                (double)distrib[(int)(g[2]*4*sysi+i)]/(maxtime*pow(2,g[2]+1)),
                //20
                (double)distrib[(int)(g[3]*4*sysi+i)]/(maxtime*pow(2,g[3]+1)),
                (double)distrib[(int)(g[4]*4*sysi+i)]/(maxtime*pow(2,g[4]+1)),
                (double)distrib[(int)(g[5]*4*sysi+i)]/(maxtime*pow(2,g[5]+1)),
                (double)distrib[(int)(g[6]*4*sysi+i)]/(maxtime*pow(2,g[6]+1)),
                (double)distrib[(int)(g[7]*4*sysi+i)]/(maxtime*pow(2,g[7]+1)),
                (double)distrib[(int)(g[8]*4*sysi+i)]/(maxtime*pow(2,g[8]+1)),
                (double)distrib[(int)(g[9]*4*sysi+i)]/(maxtime*pow(2,g[9]+1))
               );
    }
    fclose(fp);
    return(0);
}


void str_factor( long itime )
{
//	double *str, *strav1d,*stravalld,*strav1p,*stravallp, *stravall, *strav2;
    const double pi=3.1415926539, twopi=6.2831853078;
    int i,j, iq, jq, mon;
    double theta, phi, qx[NNN], qy[NNN], qz[NNN];
    double *qvec;
    double cossqr, sinsqr, addcos, addsin, arg;
//	int NNN=20;

    qvec = (double *) calloc(qmax,sizeof(double));
    for (iq=0; iq<qmax; iq++)
    {
        qvec[iq]=(double)iq*4/qmax;
    }

    for (jq=0; jq<NNN; jq++)
    {
        do  //dlya ravnomernosti _napravleniy_
        {
            theta=(pi*rand())/RAND_MAX; //random theta from 0 to pi
        }
        while ( (double)rand() / RAND_MAX >= sin(theta) );
        phi=(twopi*rand())/RAND_MAX;	//random phi from 0 to twopi = 6.28...
        qx[jq]=sin(theta)*cos(phi);
        qy[jq]=sin(theta)*sin(phi);
        qz[jq]=cos(theta);
    }


    for (i=0; i<npm; i++)// dlya 1 dendrimera
    {

        for (iq=0; iq<qmax; iq++)
        {
            cossqr=0;
            sinsqr=0;
            for (jq=0; jq<NNN; jq++)
            {
                addcos=0;
                addsin=0;
                for (j=0; j<lpm; j++)
                {
                    mon=i*lpm+j;
                    arg=qvec[iq]*(x[mon]*qx[jq]+y[mon]*qy[jq]+z[mon]*qz[jq])*di;
                    addcos+=cos(arg);
                    addsin+=sin(arg);
                }
                cossqr+=(addcos*addcos);
                sinsqr+=(addsin*addsin);
            }
            str[iq]=(cossqr+sinsqr)/(lpm*NNN);
            strav1d[iq]+=str[iq];
//			strav2[iq]+=str[iq]*str[iq];
        }
    }

    for (iq=0; iq<qmax; iq++)//dlya vseh dendrimerov
    {
        cossqr=0;
        sinsqr=0;
        for (jq=0; jq<NNN; jq++)
        {
            addcos=0;
            addsin=0;
            for (j=0; j<lpm*npm; j++)
            {
                arg=qvec[iq]*(x[j]*qx[jq]+y[j]*qy[jq]+z[j]*qz[jq])*di;
                addcos+=cos(arg);
                addsin+=sin(arg);
            }
            cossqr+=(addcos*addcos);
            sinsqr+=(addsin*addsin);
        }
        str[iq]=(cossqr+sinsqr)/(lpm*npm*NNN);
        stravalld[iq]+=str[iq];
    }

    for (i=0; i<nlpm; i++)// dlya 1 polimera
    {

        for (iq=0; iq<qmax; iq++)
        {
            cossqr=0;
            sinsqr=0;
            for (jq=0; jq<NNN; jq++)
            {
                addcos=0;
                addsin=0;
                for (j=0; j<llpm; j++)
                {
                    mon=npm*lpm+i*llpm+j;
                    arg=qvec[iq]*(x[mon]*qx[jq]+y[mon]*qy[jq]+z[mon]*qz[jq])*di;
                    addcos+=cos(arg);
                    addsin+=sin(arg);
                }
                cossqr+=(addcos*addcos);
                sinsqr+=(addsin*addsin);
            }
            str[iq]=(cossqr+sinsqr)/(llpm*NNN);
            strav1p[iq]+=str[iq];
        }
    }


    for (iq=0; iq<qmax; iq++)// dlya vseh polimerov
    {
        cossqr=0;
        sinsqr=0;
        for (jq=0; jq<NNN; jq++)
        {
            addcos=0;
            addsin=0;
            for (j=0; j<llpm*nlpm; j++)
            {
                mon=npm*lpm+j;
                arg=qvec[iq]*(x[mon]*qx[jq]+y[mon]*qy[jq]+z[mon]*qz[jq])*di;
                addcos+=cos(arg);
                addsin+=sin(arg);
            }
            cossqr+=(addcos*addcos);
            sinsqr+=(addsin*addsin);
        }
        str[iq]=(cossqr+sinsqr)/(llpm*nlpm*NNN);
        stravallp[iq]+=str[iq];
    }

    for (iq=0; iq<qmax; iq++)// dlya vseh 4astic
    {
        cossqr=0;
        sinsqr=0;
        for (jq=0; jq<NNN; jq++)
        {
            addcos=0;
            addsin=0;
            for (j=0; j<llpm*nlpm+npm*lpm; j++)
            {
                arg=qvec[iq]*(x[j]*qx[jq]+y[j]*qy[jq]+z[j]*qz[jq])*di;
                addcos+=cos(arg);
                addsin+=sin(arg);
            }
            cossqr+=(addcos*addcos);
            sinsqr+=(addsin*addsin);
        }
        str[iq]=(cossqr+sinsqr)/((lpm*npm+llpm*nlpm)*NNN);
        stravall[iq]+=str[iq];
    }


}


//output of the structure factor
int put_strfact(char *fname)
{

    int  i;
    FILE *fp;
    fp=fopen(fname,"w");
    if(fp==(NULL)) return(-1);
    fprintf( fp,"q\t1De\tallDe\t1Pe\tAllPe\tAll\n");
    for( i=0; i<qmax; i++)
    {
        fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",(double)i*4/qmax,strav1d[i]*Cmsr/mcsmax/npm,stravalld[i]*Cmsr/mcsmax,
                strav1p[i]*Cmsr/mcsmax/nlpm,stravallp[i]*Cmsr/mcsmax,stravall[i]*Cmsr/mcsmax);
    }
    fclose(fp);
    return(0);

}

void count_order_parameters(int step)
{
    int i, j, cur_cont;
    int lev_cont_num[15];
    double c;

    for(i=0; i<15; i++)
    {
        lev_cont_num[i]=0;
    }

    for(i=0; i<clasterNum; i++)//пробегаем по всем контактам, суммирум сколько связей какого уровня есть в системе
    {
        for(j=0; j<clasterNum; j++)
        {
            cur_cont=claster_cont[j*clasterNum+i];
            if(cur_cont) lev_cont_num[cur_cont]++;
        }
    }
    //число заполненных кластеров данного уровня
    for(i=1; i<15; i++)
    {
        c=pow(2,i-1)*clasterNum;//число цзлов в 1 кластере * длина кластера * ширини кластера * число кластеров данного уровня
        order_parameters[i]+=lev_cont_num[i]/c;
        order_parameters_track[i+15*step]=lev_cont_num[i]/c;
    }

}

int put_order_parameters(char *fname)
{

    int  i, j;
    FILE *fp;
    fp=fopen(fname,"w");
    if(fp==(NULL)) return(-1);
    for( i=1; i<15; i++)
    {
        fprintf(fp,"%lf\n",(double)order_parameters[i]*Cmsr/mcsmax);
    }
    fclose(fp);

    fp=fopen("ord_track.dat","w");
    if(fp==(NULL)) return(-1);
    for( j=1; j<mcsmax/Cmsr; j++)
    {
        for( i=1; i<15; i++)
        {
            fprintf(fp,"%lf\t",(double)order_parameters_track[i+15*j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);


    return(0);

}

void local_orientaton()
{
    double cosa, sum;
    long x1, x2, y1, y2, z1, z2;
    long NNeib, i, j;


//	distrib_eta_loc cosa = (x1*x2+y1*y2+z1*z2)/(sqrt(x1*x1+y1*y1+z1*z1)*sqrt(x2*x2+y2*y2+z2*z2))


    for(NNeib=1; NNeib<llpm-1; NNeib++) //vibiraem N`
    {
        for(sum=0, i=lpm*npm; i<lpm*npm+llpm-NNeib-1; i++)//probegaem po vsem zven`yam
        {
            for(j=1; j<=NNeib; j++)//sommiruem kosinusi
            {
                x1=xabs[i]-xabs[i+1];
                x2=xabs[i+j]-xabs[i+j+1];

                y1=yabs[i]-yabs[i+1];
                y2=yabs[i+j]-yabs[i+j+1];

                z1=zabs[i]-zabs[i+1];
                z2=zabs[i+j]-zabs[i+j+1];

                cosa = (x1*x2*di+y1*y2+z1*z2)/(sqrt(x1*x1+y1*y1+z1*z1)*sqrt(x2*x2+y2*y2+z2*z2));
//				if(NNeib==1) printf("%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",x1,x2,y1,y2,z1,z2,cosa);

                sum+=3*cosa*cosa-1;
            }

        }

        hist_eta_loc[NNeib]=sum/(2*NNeib*(llpm-NNeib));
    }

}

void put_local_orientaton(char *fname)
{
    FILE *fp;
    int i;

    fp=fopen(fname,"w");
    if(fp==(NULL)) return;
    for(i=0; i<llpm; i++)
    {
        fprintf(fp,"%d\t%lf\n",i,hist_eta_loc[i]);
    }


    fclose(fp);

}

long  put_measurment( char *fname )
{
    long  i, it, t;
    double rx, w, fac, Cv, w2;
    FILE *fp;

    fp=fopen(fname,"w");
    if(fp==(NULL)) return(-1);
    fprintf( fp,"acceptance: %lf \n", accp/(nrun*mcsmax*Nmono));
    /*
    for( i=0; i<maxtime1; i++)
    {
       t = (i+1)*Cmsr;
       //put measurement
       fprintf(fp,"%5.0lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf\n",
       t,
       energy_cl[i]/(nrun*96),
       energy[i]/(nrun*Nall),
       rg[i]*di*di*dlpm*dnpm/nrun,
       rgLn[i]/nrun,
       r2Ln[i]/nrun,
       cmcdLn[i],
       cmcd[i]
       );
    }
       /////////////////
    */
    for(i=1; i<=maxtime1; i++)
    {
        t = (i)*Cmsr;
        fprintf(fp,"%ld\t%lf\t%lf\n",t, rgLn[i]/nrun, r2Ln[i]/nrun);
    }
    fprintf(fp,"%ld\t%lf\t%lf\n", (maxtime1+1)*Cmsr, rgLn[0]/nrun, r2Ln[0]/nrun);

    fclose(fp);
    return(0);
}

void step_measure( long itime )
{
    long i, j, ic, k , it, jt, dt, i0, j0, mid1, mid2, a1, a2, a3;
    double cmsysx, cmsysy, cmsysz, cmx, cmy, cmz, cmlx, cmly, cmlz, rg0, rgLn0, rd, r0cm2, rcmcm2, e, ecl, xee, yee, zee;
    long dr1, rgx, rgy, rgz;
    double a12, dummy;
    float IT[4][4], diag[4], subdiag[4];
//initialization

    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            IT[i][j]=0;
        }
        diag[i]=0;
        subdiag[i]=0;
    }

//energy measurement
    it = itime&maxtime1;
    if (itime!=maxtime)
    {
        energy[it]=0;
        energy_cl[it]=0;
    }

    if(itime-maxtime)
    {


        if((nlpm==1) && (npm==1))
        {
            cmsysx=cmsysy=cmsysz=0.0;
            cmx=cmy=cmz=0.0;
            for(j=0; j<lpm; j++)
            {
                cmx += xabs[j];
                cmy += yabs[j];
                cmz += zabs[j];
                cmsysx += xabs[j];
                cmsysy += yabs[j];
                cmsysz += zabs[j];
            }
            cmx *= dlpm;
            cmy *= dlpm;
            cmz *= dlpm;

            cmlx=cmly=cmlz=0.0;
            for(j=lpm; j<lpm+llpm; j++)
            {
                cmlx += xabs[j];
                cmly += yabs[j];
                cmlz += zabs[j];
                cmsysx += xabs[j];
                cmsysy += yabs[j];
                cmsysz += zabs[j];
            }
            cmlx *= dllpm;
            cmly *= dllpm;
            cmlz *= dllpm;
            cmsysx *= dallmono;
            cmsysy *= dallmono;
            cmsysz *= dallmono;

            rd=0.0;
            rd += square(cmsysx-cmx);
            rd += square(cmsysy-cmy);
            rd += square(cmsysz-cmz);
            rd=sqrt(rd)*di;

            distrib_cmall_cmdend[(int)(rd*4)]++;

            rd=0.0;
            rd += square(cmlx-cmx);
            rd += square(cmly-cmy);
            rd += square(cmlz-cmz);
            rd=sqrt(rd)*di;
            distrib_cmdend_linecm[(int)(rd*4)]++;


        }
        else
        {
            cmsysx=cmsysy=cmsysz=0.0;//cm of system
            for(j=0; j<lpm*npm+llpm*nlpm; j++)
            {
                cmsysx += xabs[j];
                cmsysy += yabs[j];
                cmsysz += zabs[j];
            }
            cmsysx *= dallmono;
            cmsysy *= dallmono;
            cmsysz *= dallmono;
        }

        for(i=0; i<lpm*npm+llpm*nlpm; i++)//all monomer units
        {
            rd=0.0;
            rd += square(xabs[i]-cmsysx);
            rd += square(yabs[i]-cmsysy);
            rd += square(zabs[i]-cmsysz);
            rd=sqrt(rd)*di;
            distrib_all[(int)(rd*4)]++;
            if(charge[i]) distrib_all_ch[(int)(rd*4)]++;
        }

        for(ic=0; ic<nlpm; ic++)//polyelectrolyte
        {
            i=npm*lpm+ic*llpm;

            cmlx=cmly=cmlz=0.0;
            for(j=0; j<llpm; j++)
            {
                cmlx += xabs[j+i];
                cmly += yabs[j+i];
                cmlz += zabs[j+i];
            }
            cmlx *= dllpm;
            cmly *= dllpm;
            cmlz *= dllpm;

            for (j=0; j<llpm; j++)
            {
                rd=0.0;
                rd += square(xabs[i+j]-cmlx);
                rd += square(yabs[i+j]-cmly);
                rd += square(zabs[i+j]-cmlz);
                rd=sqrt(rd)*di;
                distrib_pe[(int)(rd*4)]++;
                if(charge[i+j]) distrib_pe_ch[(int)(rd*4)]++;
            }
        }


    }

    for(ic=0; ic<npm; ic++)//dendrimer
    {
        i = ic * lpm;
        k = (it * npm) + ic;
        //center of mass of dendrimer
        cmx=cmy=cmz=0.0;
        for(j=0; j<lpm; j++)
        {
            cmx += xabs[i+j];
            cmy += yabs[i+j];
            cmz += zabs[i+j];
        }
        cmx *= dlpm;
        cmy *= dlpm;
        cmz *= dlpm;
        x_cm[k] = cmx;
        y_cm[k] = cmy;
        z_cm[k] = cmz;
        cmcd[it]=((cmx-xabs[1])*(cmx-xabs[1])+
                  (cmy-yabs[1])*(cmy-yabs[1])+
                  (cmz-zabs[1])*(cmz-zabs[1]))*di*di;

        if(itime-maxtime)
        {

            //distance between center of mass and center monomer
            rd=0.0;
            rd += square(xabs[i]-cmx);
            rd += square(yabs[i]-cmy);
            rd += square(zabs[i]-cmz);
            rd=sqrt(rd)*di;
            distrib_dend_r0cm[(int)(rd*4)]++;


            //gyration radius: rg0
            //radial distribution functions: rd
            rg0 = 0.0;
            for(j=0; j<lpm; j++)
            {
                rd=0.0;
                rd += square(xabs[i+j]-cmx);
                rd += square(yabs[i+j]-cmy);
                rd += square(zabs[i+j]-cmz);
                rg0 += rd;
                rd=sqrt(rd)*di;
                if ((gen[i+j]>=0)&&(type[i+j]!=2))
                {
                    //branch points of specified generation
                    distrib[(int)(gen[i+j]*sysi*4+rd*4)]++;
                    //all branch points
                    distrib_brpt[(int)(rd*4)]++;
                }
                if(charge[i+j]) distrib_dend_ch[(int)(rd*4)]++;
                //all monomer units
                distrib_dend[(int)(rd*4)]++;
            }
            rg[it] += rg0;
            for(j=Nmono; j<Nall; j++)
            {
                rd=0.0;
                rd += square(xabs[i+j]-cmx);
                rd += square(yabs[i+j]-cmy);
                rd += square(zabs[i+j]-cmz);
                rd=sqrt(rd)*di;
                //counterions
                if (charge[j]>0) distrib_ci_p[(int)(rd*4)]++;
                if (charge[j]<0) distrib_ci_n[(int)(rd*4)]++;
                /*
                if (charge[j]>0) distrib_ci_p[(int)(rd*4)]++;
                if (charge[j]<0) distrib_ci_n[(int)(rd*4)]++;

                */
            }
            //polyelectrolyte
            for (j=lpm; j<Nmono; j++)
            {
                rd=0.0;
                rd += square(xabs[i+j]-cmx);
                rd += square(yabs[i+j]-cmy);
                rd += square(zabs[i+j]-cmz);
                rd=sqrt(rd)*di;
                distrib_pe_d_all[(int)(rd*4)]++;
                if(charge[i+j]) distrib_pe_d_ch[(int)(rd*4)]++;

            }

        }
    }

    //for linear relative to it`s cm

    for(rgLn[it]=0, r2Ln[it]=0, ic=0; ic<nlpm; ic++)
    {
        i = ic * llpm + npm;
        cmlx=cmly=cmlz=0;
        for(j=0; j<llpm; j++)
        {
            cmlx += xabs[i+j];
            cmly += yabs[i+j];
            cmlz += zabs[i+j];
        }
        cmlx*=dllpm;
        cmly*=dllpm;
        cmlz*=dllpm;
        cmcdLn[it] = ((cmx-xabs[1])*(cmx-xabs[1])+(cmy-yabs[1])*(cmy-yabs[1])+(cmz-zabs[1])*(cmz-zabs[1]))*di*di;

        //gyration radius: rgLn
        rgLn0=0;
        for(j=0; j<llpm; j++)
        {
            rgLn0+= (xabs[i+j]-cmx)*(xabs[i+j]-cmx);
            rgLn0+= (yabs[i+j]-cmy)*(yabs[i+j]-cmy);
            rgLn0+= (zabs[i+j]-cmz)*(zabs[i+j]-cmz);
        }
        rgLn[it]+=rgLn0*di*di*dllpm*dnlpm;

        //r2 ends
        r2Ln[it]+=(xabs[i]-xabs[i+llpm-1])*(xabs[i]-xabs[i+llpm-1]);
        r2Ln[it]+=(yabs[i]-yabs[i+llpm-1])*(yabs[i]-yabs[i+llpm-1]);
        r2Ln[it]+=(zabs[i]-zabs[i+llpm-1])*(zabs[i]-zabs[i+llpm-1]);
        r2Ln[it]=sqrt(r2Ln[it])*di*dnlpm;

    }



    //energy
    if(itime-maxtime)
    {
        for(i=0,e=0.0,ecl=0.0; i<Nall; i++)
        {
            e += nrg[i];
            ecl += clnrg[i];
        }
        energy[it] += (e/2.0);
        energy_cl[it] += (ecl/2.0);
    }
    //structure factor
//	str_factor(itime);
    //local orientation
    local_orientaton();
    // net order parameters
    count_order_parameters(itime);
    //inertia tensor
    /*
    for(ic=0; ic<npm; ic++)
    {
    	i = ic * lpm;
    	k = (it * npm) + ic;
    	cmx=cmy=cmz=0.0;
    	for(j=0; j<lpm; j++)
    	{
    		cmx += xabs[i+j]; cmy += yabs[i+j]; cmz += zabs[i+j];
    	}
    	cmx *= dlpm; cmy *= dlpm; cmz *= dlpm;//here
    	x_cm[k] = cmx;
    	y_cm[k] = cmy;
    	z_cm[k] = cmz;
    	if(itime-maxtime)
    	{
    		rg0 = 0.0;
    		for(j=lpm-64; j<lpm; j++)
    		{
    			rgx=xabs[i+j]-cmx;
    			rgy=yabs[i+j]-cmy;
    			rgz=zabs[i+j]-cmz;
    			IT[1][1]+=abs(rgy*rgy+rgz*rgz);
    			IT[1][2]+=abs(rgx*rgy);
    			IT[1][3]+=abs(rgx*rgz);
    			IT[2][2]+=abs(rgx*rgx+rgz*rgz);
    			IT[2][3]+=abs(rgy*rgz);
    			IT[3][3]+=abs(rgx*rgx+rgy*rgy);

    		}
    		IT[2][1]=IT[1][2];
    		IT[3][1]=IT[1][3];
    		IT[3][2]=IT[2][3];
    		tred2(IT, 3, diag, subdiag);
    		tqli(diag, subdiag, 3, IT);
    		for (j=0; j<4; j++)
    			diag_store[it*4+j]+=diag[j];
    	}
    }
    */





}



void charge_distrib( long itime )
{
    long i,tmp;
    for (i=0; i<Nmono; i++)
    {
        if (gen[i]>=0)
            if (charge[i]==1)
                chdb[itime*ngen+gen[i]]++;
    }
}


void init_measure(void)
{
    long itime, i, j;

    qmax=10*sysi;
    dnpm = 1.0/(double)npm;
    dlpm = 1.0/(double)lpm;
    dnlpm = 1.0/(double)nlpm;
    dllpm = 1.0/(double)llpm;
    dallmono = 1.0/((double)(npm*lpm+nlpm*llpm));

    distrib = (int*) calloc( ngen*sysi*4, sizeof(int));
    distrib_brpt = (int*) calloc( sysi*4, sizeof(int));
    distrib_ci_n = (int*) calloc( sysi*4, sizeof(int));
    distrib_ci_p = (int*) calloc( sysi*4, sizeof(int));


    distrib_all = (int*) calloc( sysi*4, sizeof(int));//relative to cm of system
    distrib_all_ch = (int*) calloc( sysi*4, sizeof(int));//relative to cm of system


    distrib_dend = (int*) calloc( sysi*4, sizeof(int));
    distrib_dend_r0cm = (int*) calloc( sysi*4, sizeof(int));
    distrib_dend_ch = (int*) calloc( sysi*4, sizeof(int));

    distrib_cmdend_linecm = (int*) calloc( sysi*4, sizeof(int));
    distrib_cmall_cmdend = (int*) calloc( sysi*4, sizeof(int));

    distrib_pe = (int*) calloc( sysi*4, sizeof(int));
    distrib_pe_ch = (int*) calloc( sysi*4, sizeof(int));
    distrib_pe_d_all = (int*) calloc( sysi*4, sizeof(int));//relative to cm of dendrimer
    distrib_pe_d_ch = (int*) calloc( sysi*4, sizeof(int));//relative to cm of dendrimer

    hist_eta_loc = (int*) calloc(llpm, sizeof(double));

    itime = maxtime;
    rg = (double *) calloc( itime, sizeof(double) );
    r2_ee = (double *) calloc( itime, sizeof(double) );
    energy = (double *) calloc(itime+2,sizeof(double) );
    energy_cl = (double *) calloc(itime+2,sizeof(double) );
    itime = maxtime*npm;
    x_cm = (double *) calloc( itime, sizeof(double) );
    y_cm = (double *) calloc( itime, sizeof(double) );
    z_cm = (double *) calloc( itime, sizeof(double) );
    x_end = (double *) calloc( itime, sizeof(double) );
    y_end = (double *) calloc( itime, sizeof(double) );
    z_end = (double *) calloc( itime, sizeof(double) );
    x_mid = (double *) calloc( itime, sizeof(double) );
    y_mid = (double *) calloc( itime, sizeof(double) );
    z_mid = (double *) calloc( itime, sizeof(double) );
    r_ee  = (double *) calloc( itime, sizeof(double) );
    cmcd  = (double *) calloc( itime, sizeof(double) );
    cmcdLn  = (double *) calloc( itime, sizeof(double) );
    rgLn = (double *) calloc( itime, sizeof(double) );
    r2Ln = (double *) calloc( itime, sizeof(double) );
    str = (double *) calloc(qmax, sizeof(double));
    strav1d = (double *) calloc(qmax, sizeof(double));
    stravalld = (double *) calloc(qmax, sizeof(double));
    strav1p = (double *) calloc(qmax, sizeof(double));
    stravallp = (double *) calloc(qmax, sizeof(double));
    stravall = (double *) calloc(qmax, sizeof(double));
    strav2 = (double *) calloc(qmax, sizeof(double));
    diag_store  = (double *) calloc( itime*4, sizeof(double) );
    chdb  = (int *) calloc( itime*ngen, sizeof(double) );

    if(distrib==NULL) goto error;
    if(distrib_brpt==NULL) goto error;
    if(distrib_all==NULL) goto error;
    if(distrib_ci_n==NULL) goto error;
    if(distrib_ci_p==NULL) goto error;

    if(distrib_all==NULL) goto error;
    if(distrib_all_ch==NULL) goto error;

    if(distrib_pe==NULL) goto error;
    if(distrib_pe_ch==NULL) goto error;
    if(distrib_pe_d_all==NULL) goto error;
    if(distrib_pe_d_ch==NULL) goto error;

    if(distrib_cmdend_linecm==NULL) goto error;
    if(distrib_cmall_cmdend==NULL) goto error;

    if(hist_eta_loc==NULL) goto error;

    if(distrib_dend==NULL) goto error;
    if(distrib_dend_r0cm==NULL) goto error;
    if(distrib_dend_ch==NULL) goto error;


    if(rg==NULL) goto error;
    if(r2_ee==NULL) goto error;
    if(energy==NULL) goto error;
    if(energy_cl==NULL) goto error;
    if(x_cm==NULL) goto error;
    if(y_cm==NULL) goto error;
    if(z_cm==NULL) goto error;
    if(x_end==NULL) goto error;
    if(y_end==NULL) goto error;
    if(z_end==NULL) goto error;
    if(x_mid==NULL) goto error;
    if(y_mid==NULL) goto error;
    if(z_mid==NULL) goto error;
    if(r_ee==NULL) goto error;
    if(cmcd==NULL) goto error;
    if(diag_store==NULL) goto error;
    if(chdb==NULL) goto error;
    goto quit;

error:
    printf("*** ERROR during allocation of the Stat arrays\n");
    exit(-1);

quit:
    for (i=0; i<sysi*ngen*4; i++)
    {
        distrib[i]=0;
    }
    for (i=0; i<sysi*4; i++)
    {
        distrib_brpt[i]=0;
        distrib_all[i]=0;
        distrib_ci_n[i]=0;
        distrib_ci_p[i]=0;
        distrib_pe[i]=0;
    }
    for (i=0; i<maxtime*4; i++)	diag_store[i]=0;
    for (i=0; i<maxtime*ngen; i++)	chdb[i]=0;
    return;
}

void put_data( long itime )
{
    char tmp[128];
    long  i, j;
    double rx, ry, rz, rxabs, ryabs, rzabs;
    FILE *fp;
    sprintf(tmp,"beads%03ld.txt",itime);
    fp=fopen(tmp,"w");
    if(fp==(NULL)) return ;
    for(i=0; i<Nall; i++)
    {
        rx=di*(x[i]+0.5);
        ry=di*(y[i]+0.5);
        rz=di*(z[i]+0.5);
        rxabs=di*(xabs[i]+0.5);
        ryabs=di*(yabs[i]+0.5);
        rzabs=di*(zabs[i]+0.5);
        if (charge[i]== 0)fprintf(fp,"%g %g %g %1d\n",rxabs,ryabs,rzabs,0);
        if (charge[i]==-1)fprintf(fp,"%g %g %g %1d\n",rxabs,ryabs,rzabs,1);
        if (charge[i]==+1)fprintf(fp,"%g %g %g %1d\n",rxabs,ryabs,rzabs,3);
    }
    fclose(fp);
    return;
}
