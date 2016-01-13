#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "machine.h"
#define TESTRUN
#ifdef IBM_PC
#include <alloc.h>
#endif

struct part
{
    struct part *next;
    long pn;
};


void initialize(void);
void finish(void);
long bitc( long );
void reset_config(void);
void x127( long  *m, long  *ran, long  n);
void move( long, long, long, long, long, long, long, long );
void create_all_hydro_bnd(double);
void tryjump (long, long, long);
int get_full_energy();

#define B2L 11  /* The lenght of the energy arrayis 2**(B2L) (Here 2^9=512) */
#define Estop (double)1.0e+10 /* The energy maximum > to reject a move */


long    *x,    *y,    *z, m1[128], m2[128], m3[128];
long    *ivo, *ran0, *ran,  *cpp0, *cpp, *npp0, *npp, *cellp, *cellm;
long *jmpa0, *jmpa, *jmpb0, *jmpb;
long    *who, *clwho, *mycell, change, clchange;
double  *xabs, *yabs, *zabs, *xold, *yold, *zold;
float *diag_store;
int *type, *gen, *charge, *chargevar, *chvarlist;
int *monomer_clasters, *mon_in_claster, *claster_distance, *claster_energy, *claster_cont;
double  *enp,  *enh,  *enc, *enc1, *enc2, *nrg, *nrg_new, *clnrg, *clnrg_new;

double *cmcd, *cmcdLn, *rgLn, *r2Ln;
struct  part   *list, *cell;
int *bnd, *hydro_bnd;
int ngen;
long   Lmax, Lmax1, npm, lpm, nlpm, llpm, dim, mcsmax, nrun, Nmono, Nobst, Nall, Ncell, Nch, Nvar;
long   npm1, lpm1, Nmono1, iseed, ihvar, ihvar1, Lmax2, c_Lmax2;
long   irange, irange1, c_irange1, irangeL, bnum, bnum1, c_bnum1, NB, NBC, NBC1, NBC2,lim;
long   Lcell,  LmaxL, cellL, EqL, mask=(long)2147483647L, brange, Cmsr, Lmsr, jump;
long   mask_x, mask_y, mask_z, one_x, one_y, one_z, maxtime, irangesq, maxtime1;
double sysi, temperature, range, hvar, spring, u0, bmin, bmax, di, dirange2, accp, bjerrum;
double us0,us1, uhydro, ucentr, group_cont, rmin_group;
int clasterNum;

int react, hydro_react, thydro;
char   fcfg[20], fin1[20], fin2[20], fin3[20],
       fout1[20], fout2[20], fout3[20], fout4[20], fout5[20], fout6[20], fout7[20], fout8[20], fout9[20];
long istep, irun;

double distr( rx, ry, rz )
double rx, ry, rz;
{
    double r2;
    r2 = rx*rx;
    r2 += ry*ry;
    r2 += rz*rz;
    return(r2);
}

void main(int argc, char * argv[])
{
    long   i, j, k, iflip, ijump, index, acc, mask1, mask2, mask3;
    long   newx, newy, newz, dx, dy, dz, r, p;
    char buf[80];
    clock_t t0, t1;
    int tmp;
    double dtime, rx, ry, rz, r2, work, dr1, dr2;
    FILE *fp;
    void get_ran();
    fp=fopen("en.dat","w");
    if (!argv[1])
        strcpy(fcfg, "in");
    else
        strcpy(fcfg, argv[1]);

    printf("%s",fcfg);
    initialize();
    create_wrl(fout5);
    mask1 = (2*ihvar)+1;
    mask2 = 65535L;
    get_ran(m1,m2,m3,"random");
    for(i=0; i<127; i++)
    {
        m1[i]&=mask;
        m2[i]&=mask2;
        m3[i]&=mask2;
    }

    reset_config();
//	step_measure((long)0);
    charge_distrib((long)0);
//	put_data((long)0);

    for(i=0, j=0; i<Nall; i++)//находим, какие номера у узлов
    {
        if(type[i]==1) // если мономер - узел
        {
            monomer_clasters[j]=i;
            j++;
        }
    }
    for(i=1; i<=log(clasterNum)/log(2); i++)//определяем энергию и радиус взаимодействий
    {
        claster_distance[i]=rmin_group;
//		claster_energy[i]=group_cont*(6-i);
        claster_energy[i]=group_cont/i;
    }
    for(i=1; i<=10; i++)//определяем сколько мономеров в кластере уровня и
    {
        mon_in_claster[i]=pow(2,i);
    }

    printf(" Let's start...\n");
    create_wrl(fout5);
    t0 = clock();
    for(istep=1; istep<=mcsmax*nrun; istep++)
    {
        x127( m1, ran0, Nall  );
        ran = ran0;
        x127( m2, cpp0, 3*Nall);
        cpp = cpp0;
        for(i=0; i<3*Nall; i++, cpp++)
        {
            *cpp *= mask1;
            *cpp >>= 16;
            *cpp -= ihvar;
        }
        cpp = cpp0;
        x127( m3, npp0, Nall  );
        npp = npp0;
        for( iflip=0; iflip < Nall; iflip++)
        {
            i = (npp0[iflip]*Nall)>>16;
            newx = ((dx=*cpp++) + x[i]);// & Lmax1;
            newy = ((dy=*cpp++) + y[i]);// & Lmax1;
            newz = ((dz=*cpp++) + z[i]);// & Lmax1;
            move( i, newx, newy, newz, dx, dy, dz, *ran++);
        } /* iflip */
        if(istep >= jump )
        {
            if(!(istep&(jump-1)))
            {
                x127( m1, ran0, Nall  );
                ran = ran0;
                x127( m2, jmpa0, Nvar  );
                jmpa = jmpa0;
                x127( m2, jmpb0, Nvar  );
                jmpb = jmpb0;
                for(ijump=0; ijump<Nvar; ijump++)
                {
                    i = chvarlist[(jmpa0[ijump]*Nvar)>>16];
                    j = chvarlist[(jmpb0[ijump]*Nvar)>>16];
                    if ((charge[i]==1)&&(charge[j]==0))
                    {
                        tryjump(i,j,*ran++);
//						printf ("hello: %ld %ld\n", charge[i], charge[j]);
                    }
                }
            }
        }
//		printf("%d-%d\t",istep,thydro);

        if(istep%thydro==0 && uhydro) //создаем водородные связи
        {
            printf("hydro- %f\n",(double)istep/mcsmax);
            create_all_hydro_bnd(1.0);
            /*sprintf(buf,"temp_data/view_%d.vrml",istep);
            create_wrl(buf);
            sprintf(buf,"temp_data/hydro_bonds_%d.txt",istep);
            put_hydro_bnd(buf);*/
        }

        if(istep >= Cmsr )
        {
            if(!(istep&(Cmsr-1)))
            {

                printf("%f\t",(double)istep/mcsmax);
                step_measure(istep>>Lmsr);
//				charge_distrib(istep>>Lmsr);

//				fprintf(fp, "%f\t%d\n",(double)istep/mcsmax, get_full_energy());

//				put_data(istep>>Lmsr);
//				istep=istep;
            }
        }

    }    /*  istep */
    t1 = clock();
    dtime = (double) (t1 - t0)*1000/(CLOCK_PER_SEC);
    work = Nmono*mcsmax*nrun;
    printf("\nTime:%.0lf, %.3lf upd./CPUsec\n", dtime, work/dtime );
    put_distrib(fout3);
    create_wrl(fout4);
    put_IT(fout6);
    put_strfact(fout7);
    put_local_orientaton(fout8);
    put_order_parameters(fout9);
    if(clasterNum) put_claster_cont("matrix.txt");
    if(uhydro) put_hydro_bnd("hydro_bonds.txt");
    put_pair_cont("matrix_pair.txt");
    put_chdb("chdb1_1.dat");
    finish();
    fclose(fp);
    exit(0);
}

void x127( long *m, long *r, long  n )
{
    long  k, *m0l;

    m0l=m;
    if(n<127)
    {
        for(k=0,m=m0l; k< 97; k++, m++) *m^=*(m+30);
        for(k=97; k<127; k++, m++) *m^=*(m-97);
        for(k=0; k<n; k++) r[k]=m0l[k];
    }
    else
    {
        for(k=0,m=m0l; k< 97; k++, m++) *m^=*(m+30);
        for(k=97; k<127; k++, m++) *m^=*(m-97);
        for(k=0; k<127; k++) r[k]=m0l[k];
        for(k=127; k<n; k+=16)
        {
            r[k]=r[k-97]^r[k-127];
            r[k+1] =r[k-96]^r[k-126];
            r[k+2] =r[k-95]^r[k-125];
            r[k+3] =r[k-94]^r[k-124];
            r[k+4] =r[k-93]^r[k-123];
            r[k+5] =r[k-92]^r[k-122];
            r[k+6] =r[k-91]^r[k-121];
            r[k+7] =r[k-90]^r[k-120];
            r[k+8] =r[k-89]^r[k-119];
            r[k+9] =r[k-88]^r[k-118];
            r[k+10]=r[k-87]^r[k-117];
            r[k+11]=r[k-86]^r[k-116];
            r[k+12]=r[k-85]^r[k-115];
            r[k+13]=r[k-84]^r[k-114];
            r[k+14]=r[k-83]^r[k-113];
            r[k+15]=r[k-82]^r[k-112];
        }
        for(k=n-127,m=m0l; k<n; k++,m++) *m=r[k];
    }
}

void xinit( long *m, long iseed )
{
    short  i, ib, x;
    static int flag=1;
    long  mask[32]= {0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80,
                     0x100,  0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x8000L,
                     0x10000L,0x20000L,0x40000L, 0x80000L, 0x100000L,0x200000L,0x400000L,0x800000L,
                     0x1000000L, 0x2000000L, 0x4000000L, 0x8000000L, 0x10000000L, 0x20000000L,
                     0x40000000L, 0x80000000L
                    };

    if(flag)
    {
        srand((unsigned)iseed);
        flag=0;
    }

    for(i=0; i<127; i++)
    {
        m[i]=0;
        for(ib=0; ib<=31; ib++)
        {
            x = (rand() >> 8) & 1;
            if( x ) m[i] |= mask[ib];
        }
    }

    /* Warming UP MZ 101 times */
    for(ib=0; ib<101; ib++)
    {
        for(i=0; i< 97;  i++) m[i]^=m[i+30];
        for(i=97; i<127; i++) m[i]^=m[i-97];
    }
}

void get_ran( x, y, z, name )
long *x, *y, *z;
char *name;
{
    FILE *fp;
    long i;
    fp=fopen(name,"r");
    if(fp==NULL)
    {
        puts("Error in reading Random file.");
        exit(-1);
    }
    for(i=0; i<127; i++, x++) fscanf(fp,"%ld",x);
    for(i=0; i<127; i++, y++) fscanf(fp,"%ld",y);
    for(i=0; i<127; i++, z++) fscanf(fp,"%ld",z);
    fclose(fp);
}
