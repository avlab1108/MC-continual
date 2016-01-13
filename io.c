
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

#ifdef IBM_PC
#include <alloc.h>
#endif

/*************************
 *                       *
 * Input/Output routines *
 *                       *
 *************************/

long get_params(char *fname)
{
    FILE *fp;
    char buf[80];
    int debug_flag=1;

    fp=fopen(fname,"r");
    if(fp==(NULL)) return(-1);

    fgets(buf,80,fp);
    sscanf(buf,"%ld ",&npm);
    if(debug_flag) printf("1|%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%d", &ngen);
    if(debug_flag) printf("1|%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%ld ",&lpm);
    if(debug_flag) printf("64|%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%ld ",&nlpm);
    if(debug_flag) printf("0|%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%ld ",&llpm);
    if(debug_flag) printf("0|%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%ld ",&Nmono);
    if(debug_flag) printf("64|%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%ld ",&Nch);
    if(debug_flag) printf("0|%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%ld ",&Nobst);
    if(debug_flag) printf("0|%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%ld ",&mcsmax);
    if(debug_flag) printf("26000|%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%ld ",&nrun);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&sysi);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&temperature);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&range);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&hvar);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&spring);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&u0);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&uhydro);
    if(debug_flag) printf("%s|%f\n", buf,uhydro);
    fgets(buf,80,fp);
    sscanf(buf,"%d ",&thydro);
    if(debug_flag) printf("%s|%d\n", buf,thydro);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&ucentr);
    if(debug_flag) printf("%s|%lf\n", buf,ucentr);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&us0);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&us1);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&rmin_group);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&group_cont);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%d ",&clasterNum);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&bjerrum);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&bmin);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%lf ",&bmax);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%ld ",&dim);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%ld ",&iseed);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%s ", fin1);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%s ", fin2);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%s ", fin3);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%s ", fout1);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%s ", fout2);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%s ", fout3);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%s ", fout4);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%s ", fout5);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%s ", fout6);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%s ", fout7);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%s ", fout8);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%s ", fout9);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%ld", &irange);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%ld", &Cmsr);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%d", &react);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%d", &hydro_react);
    if(debug_flag) printf("%s\n", buf);
    fgets(buf,80,fp);
    sscanf(buf,"%ld", &jump);
    if(debug_flag) printf("%s\n", buf);

//printf("\n%ld\n%ld\n%ld\n%ld\n%ld\n%ld\n%ld\n%ld\n%ld\n%ld\n",npm,ngen,lpm,nlpm,llpm,Nmono,Nch,Nobst,mcsmax,nrun);
//printf("%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n" ,sysi, temperature,range,hvar,spring,u0,us0,us1,bjerrum,bmin,bmax);
//printf("%ld\n%ld\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n" ,dim,iseed, fin1,fin2,fin3,fout1,fout2,fout3,fout4,fout5,fout6,fout7);

    ngen++;

    fclose(fp);

    return(0);
}


long  get_monomers(char *fname)
{
    long  i, j=0;
    double rx, ry, rz;
    int tp, gn, ch, chvar;
    FILE *fp;
    char buf[80];
    int flag;


    fp=fopen(fname,"r");
    if(fp==(NULL)) return(-1);

    for(i=0; i< Nmono; i++)
    {
        fgets(buf,80,fp);
        sscanf(buf,"%lf %lf %lf %d %d %d %d", &rx, &ry, &rz, &tp, &ch, &chvar, &gn);
        xold[i]=xabs[i]=rx/di;
        yold[i]=yabs[i]=ry/di;
        zold[i]=zabs[i]=rz/di;
        x[i] =(long)(xabs[i]);// & Lmax1;
        y[i] =(long)(yabs[i]);// & Lmax1;
        z[i] =(long)(zabs[i]);// & Lmax1;
        type[i]=tp;
        chargevar[i]=chvar;
        if (chvar==0)    charge[i]=ch;
        else
        {
            charge[i]=0;
            chvarlist[j]=i;
            j++;
        }
        gen[i]=gn;
    }
    Nvar=j;
    printf("Nvar: %ld\n",Nvar);
    fclose(fp);

    if(Nvar)
    {
        for(i=0; i<Nch; i++)
        {
            flag=0;
            while (flag==0)
            {
                j=rand()%Nvar;
                if (charge[chvarlist[j]]==0)
                {
                    charge[chvarlist[j]]=1;
                    flag=1;
                }
            }
        }
    }

    return(0);
}


int get_bonds(char *fname)
{
    long  i, j;
    int b[6];
    FILE *fp;
    char buf[80];
    fp=fopen(fname,"r");
    if(fp==(NULL)) return(-1);

    for(i=0; i< Nmono; i++)
    {
        fgets(buf,80,fp);
        for (j=0; j<react; j++)
        {
            b[j]=-1;
            bnd[react*i+j]=b[j];
        }
        if (sscanf(buf,"%ld %d %d %d %d %d %d", &j, &b[0], &b[1], &b[2], &b[3], &b[4], &b[5])!=EOF)
        {
            for (j=0; j<react; j++)
            {
                bnd[react*i+j]=b[j];
            }
        }
    }
    fclose(fp);
    return(0);
}

int get_hydro(char *fname) //водородные связи
{
    long  i, j;
    int monN, tmp;
    FILE *fp;
    char buf[20];
    fp=fopen(fname,"r");
    if(fp==(NULL)) return(-1);
    for(i=0; i<Nmono; i++)
    {
        for(j=0; j<hydro_react; j++)
        {
            hydro_bnd[i*hydro_react+j]=-10; // -10 - связи нет и быть не может, -1 - возможна связь, положительное число - звено, с которым установлена связь
        }
    }
    j=0;
    while (fscanf(fp, "%d\t", &tmp) !=EOF)
    {
        if(tmp<Nmono) //номер не больше,чем число мономеров
        {
            if(j==0) //номер звена, для которого потом записываются все связи
            {
                if(tmp>=0) monN=tmp;
                else printf("Negative monomer number = %d\n", tmp);
            }
            else hydro_bnd[monN*hydro_react+j-1]=tmp; //

            if(j==(hydro_react)) j=0;
            else j++;
        }
        else
        {
            printf("Big monomer number = %d\n", tmp);
        }
    }
    fclose(fp);
    return(0);
}


long  get_obstacles(char *fname)
{
    long  i, j;
    double rx, ry, rz;
    int tp, ch, gn;
    FILE *fp;
    char buf[80];

    fp=fopen(fname,"r");
    if(fp==(NULL)) return(-1);

    for(i=0; i< Nobst; i++)
    {
        j = Nmono + i;
        fgets(buf,80,fp);
        sscanf(buf,"%lf %lf %lf %d %d %d", &rx, &ry, &rz, &tp, &ch, &gn);
        xold[j]=xabs[j]=rx/di;
        yold[j]=yabs[j]=ry/di;
        zold[j]=zabs[j]=rz/di;
        x[j] =(long)(xabs[j]) & Lmax1;
        y[j] =(long)(yabs[j]) & Lmax1;
        z[j] =(long)(zabs[j]) & Lmax1;
        xold[j]=xabs[j]=x[j];
        yold[j]=yabs[j]=y[j];
        zold[j]=zabs[j]=z[j];
        type[j]=tp;
        charge[j]=ch;
        gen[j]=gn;
    }
    fclose(fp);
    return(0);
}

long  get_claster_cont(char *fname)
{
    long  i, j, tmp;
    FILE *fp;
    char buf[3];

    fp=fopen(fname,"r");
    if(fp==(NULL)) return(-1);

    for(i=0; i< clasterNum; i++)
    {
        for(j=0; j< clasterNum; j++)
        {
            fscanf(fp,"%ld", &tmp);
            claster_cont[i*clasterNum+j]=claster_cont[j*clasterNum+i]=tmp;
            //printf("%d ",claster_cont[i*clasterNum+j]);

//        claster_cont[i*clasterNum+j]=claster_cont[j*clasterNum+i]=0;
        }
//    printf("\n");
    }
    fclose(fp);

    return(0);
}

long  put_monomers(char *fname)
{
    long  i, j;
    double rx, ry, rz, rxabs, ryabs, rzabs;
    int tp, ch, chvar, gn;
    FILE *fp;

    fp=fopen(fname,"w");
    if(fp==(NULL)) return(-1);

    for(i=0; i<Nall; i++)
    {

        rxabs=di*(xabs[i]+0.5);
        ryabs=di*(yabs[i]+0.5);
        rzabs=di*(zabs[i]+0.5);
        rx=di*(x[i]+0.5);
        ry=di*(y[i]+0.5);
        rz=di*(z[i]+0.5);
        tp=type[i];
        ch=charge[i];
        chvar=chargevar[i];
        gn=gen[i];
        fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\n", rx, ry, rz, tp, ch, chvar, gn );
    }
    fclose(fp);
    return(0);
}


void  put_hydro_bnd(char *fname)
{
    int  i, j, notnull;
    FILE *fp;
    fp=fopen(fname,"w");
    if(fp==(NULL))
    {
        printf("Enable to write to file %s",fname);
        return;
    }

    for(i=0; i<Nmono; i++)
    {
        notnull=0;
        for(j=0; j<hydro_react; j++)
        {
            if(hydro_bnd[i*hydro_react+j]!=-10) notnull=1;
        }

        if(notnull)
        {
            fprintf(fp,"%d",i);
            for(j=0; j<hydro_react; j++)
            {
                fprintf(fp,"\t%d", hydro_bnd[i*hydro_react+j]);
            }
            fprintf(fp,"\n");
        }
    }
    fclose(fp);
    return ;
}

long  put_claster_cont(char *fname)
{
    long  i, j;
    FILE *fp;
    fp=fopen(fname,"w");
    if(fp==(NULL)) return(-1);

    for(i=0; i<clasterNum; i++)
    {
        for(j=0; j<clasterNum; j++)
        {
            fprintf(fp,"%d ", claster_cont[i*clasterNum+j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    return(0);
}

long  put_pair_cont(char *fname)
{
    long  i, j, m1, m2, lev;
//double tmp;
    FILE *fp;
    fp=fopen(fname,"w");
    if(fp==(NULL)) return(-1);

    for(i=0; i<clasterNum; i++)
    {
        for(j=0; j<clasterNum; j++)
        {
            m1=monomer_clasters[i];
            m2=monomer_clasters[j];
            lev=claster_cont[i*clasterNum+j];
            if(dist_real(x[m1]-x[m2],y[m1]-y[m2],z[m1]-z[m2])<claster_distance[lev])
            {
                fprintf(fp,"1 ");
            }
            else fprintf(fp,"0 ");

        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    return(0);
}

void create_wrl(char *fname)
{
    long  i, j;
    double rx, ry, rz, rxabs, ryabs, rzabs;
    FILE *fp;

    fp=fopen(fname,"w");
    if(fp==(NULL)) return;

    fprintf(fp,"#VRML V2.0 utf8\n");
    /*
    fprintf(fp,"DEF topview Viewpoint {\n");
    fprintf(fp,"position 32 22 50\n");
    fprintf(fp,"orientation 0 1 0 0.0\n");
    fprintf(fp,"description \"top\"\n");
    fprintf(fp,"}\n");
    */

    fprintf(fp,"DirectionalLight { direction -1 -1 -1 }\n");
    fprintf(fp,"DirectionalLight { direction 1 1 1 }\n\n");

    fprintf(fp,"DEF Ball0 Shape {\n");
    fprintf(fp,"geometry Sphere { radius 0.5 }\n");
    fprintf(fp,"appearance Appearance {\n");
    fprintf(fp,"material Material {\n");
    fprintf(fp,"diffuseColor  .00 .000  .000\n");
    fprintf(fp,"shininess 1.0\n");
    fprintf(fp,"transparency  .4\n}}}\n\n");

    fprintf(fp,"DEF Ball1 Shape {\n");
    fprintf(fp,"geometry Sphere { radius 0.2 }\n");
    fprintf(fp,"appearance Appearance {\n");
    fprintf(fp,"material Material {\n");
    fprintf(fp,"diffuseColor  1.000 0.00  .000\n");
    fprintf(fp,"shininess 1.0\n");
    fprintf(fp,"transparency  .0\n}}}\n\n");

    fprintf(fp,"DEF Ball2 Shape {\n");
    fprintf(fp,"geometry Sphere { radius 0.2 }\n");
    fprintf(fp,"appearance Appearance {\n");
    fprintf(fp,"material Material {\n");
    fprintf(fp,"diffuseColor  .00 1.00  .000\n");
    fprintf(fp,"shininess 1.0\n");
    fprintf(fp,"transparency  .0\n}}}\n\n");

    fprintf(fp,"DEF Ball3 Shape {\n");
    fprintf(fp,"geometry Sphere { radius 0.5 }\n");
    fprintf(fp,"appearance Appearance {\n");
    fprintf(fp,"material Material {\n");
    fprintf(fp,"diffuseColor  0.000 0.00  1.000\n");
    fprintf(fp,"shininess 1.0\n");
    fprintf(fp,"transparency  .0\n}}}\n\n");

    fprintf(fp,"DEF Ball4 Shape {\n");
    fprintf(fp,"geometry Sphere { radius 0.5 }\n");
    fprintf(fp,"appearance Appearance {\n");
    fprintf(fp,"material Material {\n");
    fprintf(fp,"diffuseColor  1.000 1.000  .000\n");
    fprintf(fp,"shininess 1.0\n");
    fprintf(fp,"transparency  .0\n}}}\n\n");

    fprintf(fp,"DEF Ball5 Shape {\n");
    fprintf(fp,"geometry Sphere { radius 0.5 }\n");
    fprintf(fp,"appearance Appearance {\n");
    fprintf(fp,"material Material {\n");
    fprintf(fp,"diffuseColor  1.000 0.0  1.000\n");
    fprintf(fp,"shininess 1.0\n");
    fprintf(fp,"transparency  .0\n}}}\n\n");

    fprintf(fp,"DEF Ball6 Shape {\n");
    fprintf(fp,"geometry Sphere { radius 0.5 }\n");
    fprintf(fp,"appearance Appearance {\n");
    fprintf(fp,"material Material {\n");
    fprintf(fp,"diffuseColor  .000 1  1.000\n");
    fprintf(fp,"shininess 1.0\n");
    fprintf(fp,"transparency  .0\n}}}\n\n");

    fprintf(fp,"DEF Ball7 Shape {\n");
    fprintf(fp,"geometry Sphere { radius 0.5 }\n");
    fprintf(fp,"appearance Appearance {\n");
    fprintf(fp,"material Material {\n");
    fprintf(fp,"diffuseColor  0.8  0.80  0.800\n");
    fprintf(fp,"shininess 1.0\n");
    fprintf(fp,"transparency  .0\n}}}\n\n");

    fprintf(fp,"DEF Ball8 Shape {\n");
    fprintf(fp,"geometry Sphere { radius 0.5 }\n");
    fprintf(fp,"appearance Appearance {\n");
    fprintf(fp,"material Material {\n");
    fprintf(fp,"diffuseColor  0.000 0.00  0.0\n");
    fprintf(fp,"shininess 1.0\n");
    fprintf(fp,"transparency  .0\n}}}\n\n");

    fprintf(fp,"DEF Ball9 Shape {\n");
    fprintf(fp,"geometry Sphere { radius 0.5 }\n");
    fprintf(fp,"appearance Appearance {\n");
    fprintf(fp,"material Material {\n");
    fprintf(fp,"diffuseColor  .000 0.00  .000\n");
    fprintf(fp,"shininess 1.0\n");
    fprintf(fp,"transparency  .0\n}}}\n\n");

    fprintf(fp,"DEF Ball10 Shape {\n");
    fprintf(fp,"geometry Sphere { radius 0.5 }\n");
    fprintf(fp,"appearance Appearance {\n");
    fprintf(fp,"material Material {\n");
    fprintf(fp,"diffuseColor  0 0 0\n");
    fprintf(fp,"shininess 1.0\n");
    fprintf(fp,"transparency  .0\n}}}\n\n");

    for(i=0; i<Nall; i++)
    {
        rx=di*(x[i]+0.5);
        ry=di*(y[i]+0.5);
        rz=di*(z[i]+0.5);
        rxabs=di*(xabs[i]+0.5);
        ryabs=di*(yabs[i]+0.5);
        rzabs=di*(zabs[i]+0.5);

        fprintf(fp,"Transform {\n");
        fprintf(fp,"translation\t%lf\t%lf\t%lf\n",rx, ry, rz);
        fprintf(fp,"children [\n");
        //деление на 4 части
/*        if(i<Nall/4)
        {
          fprintf(fp,"USE Ball1\n");
        }
        else if((i<Nall/2))
        {
          fprintf(fp,"USE Ball1\n");
        }
        else if((i<3*Nall/4))
        {
          fprintf(fp,"USE Ball2\n");
        }
        else
        {
          fprintf(fp,"USE Ball2\n");
        }
*/

        if(type[i]==-10)
        {
            fprintf(fp,"USE Ball2\n");
        }
        else
        {
            fprintf(fp,"USE Ball1\n");
        }


        fprintf(fp,"]\n}\n\n");
    }

    fprintf(fp, "Shape {\n");
    fprintf(fp, "appearance Appearance {\n");
    fprintf(fp, "material Material {\n");
    fprintf(fp, "emissiveColor 0 0 0\n}\n}");
    fprintf(fp, "geometry IndexedLineSet {\n");
    fprintf(fp, "coord Coordinate {\n");
    fprintf(fp, "point [\n");
    for (i=0; i<Nmono; i++)
    {
        rx=di*(x[i]+0.5);
        ry=di*(y[i]+0.5);
        rz=di*(z[i]+0.5);
        rxabs=di*(xabs[i]+0.5);
        ryabs=di*(yabs[i]+0.5);
        rzabs=di*(zabs[i]+0.5);
        if (i!=Nall) fprintf(fp,"\t%lf\t%lf\t%lf,\n",rx, ry, rz);
        else fprintf(fp,"\t%lf\t%lf\t%lf\n",rx, ry, rz);
    }
    fprintf(fp,"]\n}\n");
    fprintf(fp,"coordIndex [\n");
    for (i=0; i<Nmono; i++)
    {
        for (j=0; j<react; j++)
        {
            if (bnd[i*react+j]!=-1)
            {
                fprintf(fp,"\t%ld, %d, -1,\n",i,bnd[i*react+j]);
            }
        }
    }
    fprintf(fp,"]\n}\n}\n");

    fclose(fp);
}


/* Set the initial constants and allocate the arrays */

long init_mem()
{
    long i, j, w;

    w = (long)sysi;
    if( bitc(w) != 1 )
    {
        printf("*** Error incorrect ratio sysi/range : %ld\n", w);
        exit(-1);
    }
    i = bitc(w-1);
    irangeL = 16-i;
    irange=(long)1<<irangeL;
    irange1 = irange - 1;
    dirange2 = pow((double)irange,2.0);
    c_irange1 = ~irange1;

    Ncell = w * w * w;
    Lcell = w;
    Lmax = w * irange;
    Lmax1 = Lmax - 1;
    LmaxL = bitc(Lmax1);
    Lmax2 = Lmax >> 1;

    di = sysi / (double)Lmax;
    cellL = LmaxL - irangeL;
    /* cellL=bitc(w-1); */


    if(bitc(mcsmax)!=1)
    {
        for(w=0; mcsmax; w++, mcsmax>>=1);
        mcsmax = (long)1<<(w-1);
        printf("*** Warning: incorrect mcsmax : changed to %ld\n", mcsmax);
    }

    if(bitc(Cmsr)!=1)
    {
        for(w=0; Cmsr; w++, Cmsr>>=1);
        Cmsr = (long)1<<(w-1);
        printf("*** Warning: incorrect Cmsr : changed to %ld\n", Cmsr);
    }

    if(bitc(jump)!=1)
    {
        for(w=0; jump; w++, jump>>=1);
        jump = (long)1<<(w-1);
        printf("*** Warning: incorrect jump step : changed to %ld\n", jump);
    }
    Lmsr = bitc(Cmsr-1);
    maxtime = mcsmax/Cmsr;
    maxtime1 = maxtime - 1;

    Nall  = Nmono + Nobst;
    mask_z = ((long)1<<cellL) - 1;
    mask_y = mask_z<<cellL;
    mask_x = mask_y<<cellL;
    one_z = (long)1;
    one_y = (long)1<<cellL;
    one_x = (long)1<<(2*cellL);

    bnum = (long)1<<B2L;
    bnum1 = bnum - 1;
    NB=(irangeL*2) - B2L;
    irangesq = (bnum<<NB)-1; /* In fact : (irange^2)-1 */

    printf("NB:%ld :%ld %ld\n",NB, irange, irangesq);
    ihvar = hvar / di;
    if(bitc(ihvar)!=1)
    {
        for(w=0; ihvar; w++, ihvar>>=1);
        ihvar = (long)1<<(w-1);
        printf("*** Warning: incorrect ihvar : changed to %lf\n", ihvar*di);
    }
    ihvar1 = ihvar - 1;
    c_bnum1 = bnum1 ^ (long)(-1);
    EqL = ((bmin + bmax)/2)/di;
    brange = ((bmax - bmin)/2.0)/di;

    x = (long *) calloc( (size_t) Nall, sizeof(long) );
    y = (long *) calloc((size_t)  Nall, sizeof(long) );
    z = (long *) calloc((size_t)  Nall, sizeof(long) );
    xabs = (double *) calloc((size_t)  Nall, sizeof(double) );
    yabs = (double *) calloc((size_t)  Nall, sizeof(double) );
    zabs = (double *) calloc((size_t)  Nall, sizeof(double) );
    xold = (double *) calloc((size_t)  Nall, sizeof(double) );
    yold = (double *) calloc((size_t)  Nall, sizeof(double) );
    zold = (double *) calloc((size_t)  Nall, sizeof(double) );
    type = (int *) calloc((size_t)  Nall, sizeof(int) );
    gen = (int *) calloc((size_t)  Nall, sizeof(int) );
    charge = (int *) calloc((size_t)  Nall, sizeof(int) );
    chargevar = (int *) calloc((size_t)  Nall, sizeof(int) );
    monomer_clasters = (int *) calloc((size_t)  clasterNum, sizeof(int) );
    mon_in_claster = (int *) calloc((size_t)  10, sizeof(int) );
    claster_distance = (int *) calloc((size_t)  log(clasterNum)/log(2), sizeof(int) );
    claster_energy = (int *) calloc((size_t)  log(clasterNum)/log(2), sizeof(int) );
    claster_cont = (int *) calloc((size_t)  clasterNum*clasterNum, sizeof(int) );

    order_parameters = (double *) calloc((size_t)  15, sizeof(double) );
    order_parameters_track = (double *) calloc((size_t)  15*mcsmax/Cmsr, sizeof(double) );
    chvarlist = (int *) calloc((size_t)  Nall, sizeof(int) );
    enp  = (double *) calloc((size_t)  bnum+100, sizeof(double) );
    enh  = (double *) calloc((size_t)  bnum+100, sizeof(double) );
    enc  = (double *) calloc((size_t)  bnum+100, sizeof(double) );
    enc1  = (double *) calloc((size_t)  bnum+100, sizeof(double) );
    enc2  = (double *) calloc((size_t)  bnum+100, sizeof(double) );
    ran0 = (long *) calloc((size_t)  Nall+16, sizeof(long) );
    cpp0 = (long *) calloc((size_t)  (3*Nall)+16, sizeof(long) );
    npp0 = (long *) calloc((size_t)  Nall+16, sizeof(long) );
    jmpa0 = (long *) calloc((size_t)  Nvar+16, sizeof(long) );
    jmpb0 = (long *) calloc((size_t)  Nvar+16, sizeof(long) );
    nrg  = (double *) calloc((size_t)  Nall, sizeof(double) );
    nrg_new = (double *) calloc((size_t)  Nall, sizeof(double) );
    clnrg  = (double *) calloc((size_t)  Nall, sizeof(double) );
    clnrg_new = (double *) calloc((size_t)  2*Nall, sizeof(double) );
    who = (long *) calloc((size_t)  Nall, sizeof(long) );
    clwho = (long *) calloc((size_t)  2*Nall, sizeof(long) );
    mycell = (long *) calloc((size_t)  Nall, sizeof(long) );
    cellp = (long *) calloc((size_t)  Lcell, sizeof(long) );
    cellm = (long *) calloc((size_t)  Lcell, sizeof(long) );


    list = (struct part *) calloc((size_t)  Nall, sizeof(struct part) );
    cell = (struct part *) calloc((size_t)  Ncell, sizeof(struct part) );
    bnd  = (long *) calloc((size_t)  react*2*Nall, sizeof(int) );
    hydro_bnd  = (int *) calloc((size_t)  hydro_react*Nall, sizeof(int) );

    if( cell==NULL ) goto error;
    if( list==NULL ) goto error;
    if( bnd ==NULL ) goto error;
    if( x   ==NULL ) goto error;
    if( y   ==NULL ) goto error;
    if( z   ==NULL ) goto error;
    if( xabs==NULL ) goto error;
    if( yabs==NULL ) goto error;
    if( zabs==NULL ) goto error;
    if( xold==NULL ) goto error;
    if( yold==NULL ) goto error;
    if( zold==NULL ) goto error;
    if ( type==NULL) goto error;
    if ( gen==NULL) goto error;
    if ( charge==NULL) goto error;
    if ( chargevar==NULL) goto error;
    if ( chvarlist==NULL) goto error;
    if( ran0 ==NULL ) goto error;
    if( cpp0 ==NULL ) goto error;
    if( npp0 ==NULL ) goto error;
    if( enp ==NULL ) goto error;
    if( enh ==NULL ) goto error;
    if( enc ==NULL ) goto error;
    if( enc1 ==NULL ) goto error;
    if( enc2 ==NULL ) goto error;
    if( nrg==NULL ) goto error;
    if( nrg_new==NULL ) goto error;
    if( clnrg==NULL ) goto error;
    if( clnrg_new==NULL ) goto error;
    if( who==NULL ) goto error;
    if( clwho==NULL ) goto error;
    if( mycell==NULL ) goto error;
    if( cellp==NULL ) goto error;
    if( cellm==NULL ) goto error;

    return(0);
error:
    puts("*** ERROR while allocating arrays");
    exit(-1);
    return(0);
}

void init_energy()
{
    long i, j, r, a;
    double step,step1,step2,x0, x2, e, ec, alpha, rest, eqdist;
    int diag;
    FILE *fp;
    eqdist = 0.8;
    alpha = 24.0;
    rest=exp(-2.0*alpha*(range-eqdist))-2.0*exp(-alpha*(range-eqdist));
    step = range / (double)bnum1;
    for(i=1, x0=0.0; i<bnum; i++)
    {
        x0 = i*step;
        enh[i]=Estop;
        x2 = sqrt(x0);
        if( (x2>=bmin)&&(x2<=bmax) )
        {
            enh[i] = spring*pow(x2-((bmin+bmax)/2.0),2.0);
        }
        if(x2<=1.0e-10)
        {
            enp[i]=Estop;
        }
        else
        {
//        e = u0*( (double)1.0 - (2.0*pow(x2,-6.0)) + pow(x2,-12.0) );
            e = exp(-2.0*alpha*(x2-eqdist))-2.0*exp(-alpha*(x2-eqdist)) - rest;
            if ( e > Estop )
                enp[i] = Estop;
            else
                enp[i] = e;
        }
    }
    /*
    diag=4*8*8;
    lim=8/di;
    step = diag/(double)(bnum);
    NBC = bitc((long)(diag/(di*bnum)-1));
    for(i=1, x0=0.0; i<bnum; i++)
    {
        x0 = i*step;
        enc[i]=Estop;
        x2 = sqrt(x0);
        if(x2==0)
        {
            enc[i] = Estop;
        }
        else
        {
            enc[i] = 1/x2;
        }
    }

    diag=4*sysi*sysi;
    step2 = diag/(double)(bnum);
    NBC2 = bitc((long)(diag/(di*bnum)-1));
    for(i=1, x0=0.0; i<bnum; i++)
    {
        x0 = i*step2;
        enc2[i]=Estop;
        x2 = sqrt(x0);
        if(x2<bmin)
        {
            enc2[i] = Estop;
        }
        else
        {
            enc2[i] = 1/x2;
        }
    }


    enc2[0]=enc1[0]=enc[0]=enp[0] = enh[0] = Estop;
    */
    enp[0]=enh[0]=Estop;
    for(i=bnum; i<bnum+100; i++) enh[i]=Estop;
    for(i=bnum; i<bnum+100; i++) enp[i]=0;
//for(i=bnum; i<bnum+100; i++) enc[i]=enc[bnum-1];
//for(i=bnum; i<bnum+100; i++) enc1[i]=enc1[bnum-1];
//for(i=bnum; i<bnum+100; i++) enc2[i]=enc2[bnum-1];
    /*
    fp=fopen("lj.dat","w");
    for(i=0; i<bnum+100; i++)
       fprintf(fp,"%ld %lf %lf %lf %lf %lf %lf %lf %lf\n", i, enh[i], enp[i],sqrt(step*i), enc[i], sqrt(step1*i),enc1[i],sqrt(step2*i),enc2[i]);
    fclose(fp);
    */
}

void initialize(void)
{
    long errflg;
    long i;

    errflg  =  get_params( fcfg );
    if(errflg)
    {
        puts("*** ERROR in file INPUT operations. get_params");
        exit(-1);
    }
    init_mem();
    init_measure();
#ifdef IBM_PC
    printf("Core near left:%lu\n",coreleft());
    printf("Core far  left:%lu\n",farcoreleft());
#endif
    errflg |= get_monomers(fin1);
    if(errflg)
    {
        puts("*** ERROR in file INPUT operations. get_monomers");
        exit(-1);
    }
    errflg |= get_bonds(fin3);
    if(errflg)
    {
        puts("*** ERROR in file INPUT operations. get_bonds");
        exit(-1);
    }
    if(Nobst)  errflg |=  get_obstacles( fin2 );
    if(errflg)
    {
        puts("*** ERROR in file INPUT operations. get_obstacles");
        exit(-1);
    }
    if(clasterNum)  errflg |=  get_claster_cont("matrix.txt");
    if(errflg)
    {
        puts("*** ERROR in file INPUT operations. get_claster_cont");
        exit(-1);
    }
    if(uhydro)  errflg |=  get_hydro("hydro_bonds.txt");
    if(errflg)
    {
        puts("*** ERROR in file INPUT operations. get_hydro");
        exit(-1);
    }
    if(errflg)
    {
        puts("*** ERROR in file INPUT operations");
        exit(-1);
    }
    init_energy();
    /* Init cellp & cellm */
    for(i=0; i<Lcell; i++)
    {
        cellp[i] = (i+1)&(Lcell-1);
        cellm[i] = (i-1)&(Lcell-1);
    }
    /* Init List */
    for(i=0; i<Nall; i++)
    {
        list[i].next = NULL;
        list[i].pn = i;
    }
}

void finish(void)
{
    long errflg;

    errflg = put_monomers(fout1);
    errflg |= put_measurment(fout2);
    if(errflg)
    {
        puts("*** ERROR in file OUTPUT operations");
        exit(-1);
    }
    puts("*** That's all folks...");
}

/********************************
 *                              *
 *  Some   Usefull   Routines   *
 *                              *
 ********************************/

long bitc( long w )
{
    long i;
    for(i=0; w; w &= w-1, i++);
    return(i);
}

long cell_number(long i, long j, long k)
{
    long a, b;
    a = (i >> irangeL) << (cellL<<1);
    b = (j >> irangeL) << cellL;
    a |= (( k >> irangeL ) | b );
    return(a);
}

void cell_put( long ip, long nc)
{
    list[ip].next=cell[nc].next;
    cell[nc].next=&list[ip];
}

void cell_get( long ip, long nc)
{
    struct part *xp, *yp;
    yp = &list[ip];
    for(xp=&cell[nc]; (xp->next)!=NULL; xp=xp->next)
    {
        if( xp->next==yp )
        {
            xp->next = list[ip].next;
            return;
        }
    }
    printf("*** Error : particle %ld is not in cell %ld\n", ip, nc);
    exit(-1);
}

void cell_print(void)
{
    struct part *x;
    long i;
    for(i=0; i<Ncell; i++)
    {
        printf("\n cell:%ld :",i);
        for(x=cell[i].next; x!=NULL; x=x->next)
            printf(" %ld", x->pn );
    }
}

void reset_config(void)
{
    long i, j, a, ib, ci;
    int ncall[6], iball[6];
    double e, ecl;
    for(i=0; i<Nall; i++)
    {
        xabs[i]=xold[i];
        yabs[i]=yold[i];
        zabs[i]=zold[i];
        x[i] =(long)(xabs[i]);// & Lmax1;
        y[i] =(long)(yabs[i]);// & Lmax1;
        z[i] =(long)(zabs[i]);// & Lmax1;
    }

    for(i=0; i<Ncell; i++) cell[i].next=NULL;


    for(i=0; i<Nall; i++)
    {
        a = cell_number(x[i],y[i],z[i]);
        cell_put( i, a );
        mycell[i]=a;
    } /* i */
    /* Initial set of the energies */

    for(i=0; i<Nall; i++)
    {
#ifdef TESTRUN
        cell_print();
#endif
        e = 0.0;
        ecl=0.0;
        if (i<Nmono)
        {
            for (j=0; j<react; j++)
            {
                ib=bnd[react*i+j];
                iball[j]=ib;
                if (ib!=-1)
                {
                    if (ib!=-1) cell_get( ib, ncall[j]=mycell[ib] );
                    a = dist2(x[i]-x[ib], y[i]-y[ib], z[i]-z[ib]);
                    if (a>=irangesq)
                    {
                        e += Estop;
                        printf("*** Extended bond: %ld - %ld, %ld\n",i, ib, a);
                    }
                    else
                    {
                        e += enh[a>>NB];
                    }
                }
            }
        }
        cell_get( i, mycell[i] );
        change = 0;
        e += get_energy( x[i], y[i], z[i] );
        clchange = 0;
        for (ci=0; ci<Nall; ci++)
        {
            if (i!=ci) ecl += interact_coulomb(x[i],y[i],z[i],i,ci);
        }
        nrg[i] = e;
        clnrg[i] = ecl;
        if(e >= Estop) printf("*** WARNING: initialy particle %ld has BIG ENERGY...\n",i);
        if (i<Nmono)
        {
            for(j=0; j<react; j++)
            {
                ib=iball[j];
                if (ib!=-1) cell_put( ib, ncall[j]);
            }
        }
        cell_put(i,mycell[i]);
    } /* i */
}

/* #define labs(a) ( ((a)<0) ? (-a) : (a) ) */

long distance( long dx, long dy, long dz )
{
    dx = labs(dx);
    dy = labs(dy);
    dz = labs(dz);
    dx = min(dx,Lmax-dx);
    dy = min(dy,Lmax-dy);
    dz = min(dz,Lmax-dz);
    return((long)sqrt((double)((dx*dx)+(dy*dy)+(dz*dz))));
}
/*
long dist2( dx, dy, dz )
union dword dx, dy, dz;
{
unsigned long r;
r= (long)dx.word[0]*(long)dx.word[0];
r+=(long)dy.word[0]*(long)dy.word[0];
r+=(long)dz.word[0]*(long)dz.word[0];
return(r);
}
*/

/*
long dist2( long dx, long dy, long dz )
{
long r;
dx = labs(dx); dy = labs(dy); dz = labs(dz);
dx = min(dx,Lmax-dx);
dy = min(dy,Lmax-dy);
dz = min(dz,Lmax-dz);
r = (dx*dx)+(dy*dy)+(dz*dz);
return(r);
}
*/
unsigned long dist2( dx, dy, dz )
short dx, dy, dz;
{
    unsigned long r;
    r= (long)dx*dx;
    r+=(long)dy*dy;
    r+=(long)dz*dz;
    return(r);
}

