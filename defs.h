#include "machine.h"

#define DDD(x,y,z)
#ifdef TEST
#define DDD(x,y,z)  printf("%s %ld %ld\n",x,y,z)
#endif

#define B2L 11  /* The lenght of the energy array is 2**(B2L) (Here 2^10=1024) */
#define Estop (double)1.0e+10 /* The energy maximum > to reject a move */
#define CLOCK_PER_SEC 1000000

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

struct part
{
    struct part *next;
    long pn;
};

extern int ngen;

extern long    *x,    *y,    *z, m1[128], m2[128], m3[128], *mycell, *who, *clwho;
extern long  *ivo, *ran0, *ran,  *cpp0, *cpp, *npp0, *npp, *cellp, *cellm, change, clchange;
extern long *jmpa0, *jmpa, *jmpb0, *jmpb;

extern double  *xabs, *yabs, *zabs, *xold, *yold, *zold;
extern int *type, *gen, *charge, *chargevar, *chvarlist;
extern int *monomer_clasters, *mon_in_claster, *claster_distance, *claster_energy, *claster_cont;
extern double *order_parameters, *order_parameters_track;
extern double  *enp,  *enh, *enc, *enc1, *enc2, *nrg,  *nrg_new, Emax;
extern double *cmcd, *cmcdLn, *rgLn, *r2Ln;
extern double *clnrg, *clnrg_new;
extern struct  part   *list, *cell;
extern int  *bnd, *hydro_bnd;
extern float *diag_store;

extern long   Lmax, Lmax1, npm, lpm, nlpm, llpm, dim, mcsmax, nrun, Nmono, Nobst, Nall, Ncell, Nch, Nvar;
extern long   npm1, lpm1, Nmono1, iseed, ihvar, ihvar1, Lmax2, c_Lmax2, Maxtry;
extern long   irange, irange1, c_irange1, irangeL, bnum, bnum1, c_bnum1, NB, NBC, NBC1, NBC2,lim;
extern long   Lcell,  LmaxL, cellL, EqL, mask, brange, irangesq, Cmsr, Lmsr, jump;
extern long   mask_x, mask_y, mask_z, one_x, one_y, one_z, maxtime, maxtime1;
extern double sysi, temperature, range, hvar, spring, u0, bmin, bmax, di, dirange2, accp, bjerrum;

extern double us0,us1, uhydro, ucentr, group_cont, rmin_group;
extern int clasterNum;
extern int react, hydro_react, thydro;
extern char   fcfg[20], fin1[20], fin2[20], fin3[20],
       fout1[20], fout2[20], fout3[20], fout4[20], fout5[20], fout6[20], fout7[20], fout8[20], fout9[20];
extern long   istep, irun;

long  get_params(char *fname);
long  get_monomers(char *);
long  get_obstacles(char *);
long  get_claster_cont(char *);
long  put_monomers(char *);
void create_wrl(char *);
long init_mem(void);
void ran_monomers(void);
void init_energy(void);
long put_measurment(char *);
long put_distrib(char *);
void init_measure(void);
void initialize(void);
void finish(void);
long bitc( long );
long cell_number( long, long, long );
void cell_put( long, long );
void cell_get( long, long );
void reset_config(void);

void ini127( long  *m, long iseed);
void i127( long  *m, long  *ran, long  n);

long distance( long, long, long );
double get_energy( long, long, long );
float interact_coulomb( long, long, long, long, long );
double e_mono( long, long, long, long );
void move( long, long, long, long, long, long, long, long );
void tryjump (long, long, long);

unsigned long dist2();
double dist_real(long, long, long);
/*
union dword{
   signed short word[2];
   };
*/
#define dist3(a,b,c) dist2( (short)a, (short)b,(short)c )


