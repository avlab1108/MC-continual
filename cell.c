#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

#ifdef IBM_PC
#include <alloc.h>
#endif


void print_matrix()
{
    int tmp1, tmp2;
    for(tmp1=0; tmp1<clasterNum; tmp1++)
    {
        for(tmp2=0; tmp2<clasterNum; tmp2++)
        {
            printf("%d ", claster_cont[tmp1*clasterNum+tmp2]);
        }
        printf("\n");
    }
    printf("\n");
}

/**
Центрально-симметричный потенциал
*/
double en_centrosymmetric(long mon_x, long mon_y, long mon_z)
{
    int i;
    static double cent[3];
    static int count=0;

    //mass center
    if(count%(10000*Nall)==0)
    {
        for(i=0, cent[0]=0, cent[1]=0, cent[2]=0; i<Nall; i++)
        {
            cent[0]+=x[i];
            cent[1]+=y[i];
            cent[2]+=z[i];
        }
        cent[0]=cent[0]/Nall;
        cent[1]=cent[1]/Nall;
        cent[2]=cent[2]/Nall;
        printf("central = %f\t%f\t%f\n",di*cent[0],di*cent[1],di*cent[2]);
    }
    count++;
    return ucentr*dist_real(mon_x-cent[0],mon_y-cent[1],mon_z-cent[2]);
}


void create_hydro_bnd(long cur_monom, long mon_x, long mon_y, long mon_z, double connect_dist)
{
    int j, free_bnd, connected, mon2, cur_free_bnd;
    //считаем количество свободных связей
    for(j=0, cur_free_bnd=0; j<hydro_react; j++)
    {
        if(hydro_bnd[cur_monom*hydro_react+j]==-1) cur_free_bnd++; //есть свободные связи
    }
    //проверяем образование связей
    for(mon2=0; mon2<Nmono; mon2++)
    {
        if(cur_free_bnd==0) return 0; //свободные связи кончились
        if(type[mon2]==1 && mon2!=cur_monom)
        {
            free_bnd=0;
            connected=0;
            for(j=0; j<hydro_react; j++)
            {
                if(hydro_bnd[mon2*hydro_react+j]==-1) free_bnd++; //есть свободные связи
                else if(hydro_bnd[mon2*hydro_react+j]==cur_monom) connected++; //уже связан с текущим звеном
            }

            if(!connected && free_bnd) //у обоих есть свободные связи и еще не связаны, проверяем образование связи
            {
                //if(mon2==10 && cur_monom==20) printf("%f|",di*sqrt(dist2(mon_x-x[mon2],mon_y-y[mon2],mon_z-z[mon2])));
                if(dist_real(mon_x-x[mon2],mon_y-y[mon2],mon_z-z[mon2]) < connect_dist) //достаточно близко, образуем связь
                {
                    for(j=0; j<hydro_react; j++)//ищем свободную связь у текушего звена
                    {
                        if(hydro_bnd[cur_monom*hydro_react+j]==-1)
                        {
                            hydro_bnd[cur_monom*hydro_react+j]=mon2;
                            j=hydro_react+1;
                        }
                    }
                    for(j=0; j<hydro_react; j++)//ищем свободную связь у текушего звена
                    {
                        if(hydro_bnd[mon2*hydro_react+j]==-1)
                        {
                            hydro_bnd[mon2*hydro_react+j]=cur_monom;
                            j=hydro_react+1;
                        }
                    }
                    //printf("create bond %ld-%d\n",cur_monom,mon2);
                    cur_free_bnd--;
                }
            }
        }
    }
}

void create_all_hydro_bnd(double connect_dist)
{
    int i;
    for(i=0; i<Nmono; i++)
    {
        create_hydro_bnd(i, x[i], y[i], z[i], 1.0);
    }
}

int break_hydro_bnd(long cur_monom, long mon_x, long mon_y, long mon_z, long rnd, double connect_dist)
{
    int j,k, mon2;
    double r;
//    printf("hydr %ld\n",cur_monom);
    //проверяем разрыв связей
    for(j=0; j<hydro_react; j++)
    {
        mon2=hydro_bnd[cur_monom*hydro_react+j];
        if(mon2>=0)
        {
            r=dist_real(mon_x-x[mon2],mon_y-y[mon2],mon_z-z[mon2]);
//	    printf("%f|",r);
            if(r>connect_dist)//дистанция увеличилась, попытка рвать связь
            {
                if(uhydro > 500) return 2; //если энергия большая то не разрешаем рвать связь и отменяем ход
                if((rnd<(long)(exp(-uhydro/temperature)*mask))) //рвем связь с вероятностью
                {
                    hydro_bnd[cur_monom*hydro_react+j]=-1; //пустая связь с возможностью образования
                    for(k=0; k<hydro_react; k++)//ищем у второго звена эту связь
                    {
                        if(hydro_bnd[mon2*hydro_react+k]==cur_monom) hydro_bnd[mon2*hydro_react+k]=-1; //пустая связь с возможностью образования
                    }
                    // printf("Drop bond %ld-%d | %f\n",cur_monom,mon2, (double)rnd/(long)(exp(-uhydro/temperature)*mask));
                }
                else
                {
                    //  printf("Don`t drop bond %ld-%d | %f\n",cur_monom,mon2, (double)rnd/(long)(exp(-uhydro/temperature)*mask));
                    return 2;
                }
            }
        }
    }
    return 0;
}

int group_conn(int cur_claster, int lev, int non_claster)//проверяем, есть ли хотя бы одна связь у данного кластера на данном уровне, кроме кластера non_claster
{
    int i, cur_m1, cur_m2;
    double r;
    cur_m2 = monomer_clasters[cur_claster];
    for(i=0; i<clasterNum; i++)
    {
        if(claster_cont[cur_claster*clasterNum+i]==lev && non_claster!=i)
        {
            cur_m1 = monomer_clasters[i];
            r=dist_real(x[cur_m2]-x[cur_m1],y[cur_m2]-y[cur_m1],z[cur_m2]-z[cur_m1]);//расстояние между мономерами
            if(r<claster_distance[lev]) return 1;
        }
    }
    return 0;
}

int en_of_groups(long cur_monom, long mon_x, long mon_y, long mon_z, long rnd)
{
// exp(-2a(r-rmin))-2exp[-a(r-rim)]

    double energ=0, r, range=1, test;
    int i, j, k, c, tmp, cur_claster, cur_monom_num, lev, maxLev, conn, mon_in_c, cur_group, in_group, mon_in_claster;
    int tmp1, tmp2;
    int flagg;
    int groups_vs[256];

    if(type[cur_monom]!=1) return 0;


    for(i=0; i<clasterNum/2; i++)
    {
        groups_vs[i]=0;
    }

    //определяем номер текущего кластера
    for(i=0; i<clasterNum; i++)
    {
        if (monomer_clasters[i]==cur_monom)
        {
            cur_claster=i;
            break;
        }
    }

    maxLev=0;
    for(i=0; i<clasterNum; i++)//проверяем, какие уровни уже есть у данного узла
    {
        if(claster_cont[cur_claster*clasterNum+i]>maxLev)
        {
            maxLev=claster_cont[cur_claster*clasterNum+i];
        }
    }

    for(lev=1; lev<=maxLev; lev++)//проверяем, разрыв связи всех уровней
    {

        conn=0;
        for(j=0; j<clasterNum; j++)//проверяем связь с j-ым коастером
        {
            if(claster_cont[cur_claster*clasterNum+j]==lev)
            {
                cur_monom_num = monomer_clasters[j];
                r=dist_real(mon_x-x[cur_monom_num],mon_y-y[cur_monom_num],mon_z-z[cur_monom_num]);//расстояние между мономерами
                if(r<claster_distance[lev])
                {
                    conn=1;
                    break;
                }
                conn+= group_conn(j, lev, cur_claster);//проверяем связь j-го кластера со всеми, кроме того, что подвинули.
                if(conn) break;

            }
        }
        if(conn==0)//на данном уровне нет ни одной связи
        {
//            test = exp(-claster_energy[maxLev]/temperature);
            if((lev==maxLev) && (rnd<(long)(exp(-claster_energy[maxLev]/temperature)*mask)))//рвем с вероятностью exp(energ/temp) связь
            {
//                printf("\n\n\nrnd>test %lf  %lf\n\n\n",rnd,test);

                for(j=0; j<clasterNum; j++)//удаляем связи уровня maxLev для j-го узла
                {
                    if(claster_cont[cur_claster*clasterNum+j]<=maxLev && claster_cont[cur_claster*clasterNum+j]>0)
                    {
                        for(k=0; k<clasterNum; k++)
                        {
                            if(claster_cont[j*clasterNum+k]==maxLev) claster_cont[j*clasterNum+k]=claster_cont[k*clasterNum+j]=0;
                        }
                    }

                }
                for(k=0; k<clasterNum; k++)//удаляем связи для выбранного узла
                {
                    if(claster_cont[cur_claster*clasterNum+k]==lev) claster_cont[cur_claster*clasterNum+k]=claster_cont[k*clasterNum+cur_claster]=0;
                }
                return 0;
            }
            else//запрещаем переход
            {
                return 1;
            }
        }
    }

    for(i=0; i<clasterNum; i++)//проверяем, образование связи maxLev+1 с кластером, в состав которого входит i
    {
        if(i!=cur_claster && claster_cont[cur_claster*clasterNum+i]==0)//не текущий кластер, и связи еще нет
        {
            flagg=0;
            for(j=0; j<clasterNum; j++)//находится ли данный кластер в кластере нужного уровня, нет ли у него еще такой связи
            {
                if(claster_cont[i*clasterNum+j]==maxLev)//находится
                {
                    flagg=1;
                }
                if(claster_cont[i*clasterNum+j]==maxLev+1)
                {
                    flagg=0;
                    break;
                }

            }
            if(flagg)//у i связи уровня max+1 еще нет, связь max есть
            {
                cur_monom_num = monomer_clasters[i];
                r=dist_real(mon_x-x[cur_monom_num],mon_y-y[cur_monom_num],mon_z-z[cur_monom_num]);//расстояние между мономерами
                //фиксированные пары
                in_group=0;
                mon_in_claster=pow(2, maxLev+1);//чило мономеров в кластере образуемого уровня
                for(cur_group=0; cur_group<clasterNum/mon_in_claster; cur_group++)//пробегаем по всем группам внутри уровня
                {
                    if(cur_claster>=mon_in_claster*cur_group && cur_claster<mon_in_claster*(cur_group+1) && i>=mon_in_claster*cur_group && i<mon_in_claster*(cur_group+1))//мономеры в одной    группе
                    {
                        //в одной группе
                        in_group=1;
                        break;
                    }
                }
                if(in_group && r<claster_distance[maxLev+1])//дистанция подходит. образуем новую связь, между кластерами уровня maxLev содержащими кластеры cur_claster и i
                {

                    if(maxLev==0)
                    {
//                        print_matrix();
                        claster_cont[cur_claster*clasterNum+i]=claster_cont[i*clasterNum+cur_claster]=1;
                        return 0;
                    }
                    else
                    {
                        groups_vs[0]=i;
                        for(k=0, mon_in_c=1; k<clasterNum; k++)//определяем узлы, всесте с которыми узел i образует нужный кластер
                        {
                            if(claster_cont[i*clasterNum+k]<=maxLev && claster_cont[i*clasterNum+k]>0)
                            {
                                groups_vs[mon_in_c]=k;
                                mon_in_c++;
                            }
                        }
                        for(c=0; c<mon_in_c; c++)//прописываем связи для cur_claster узлов
                        {
                            tmp=groups_vs[c];
                            claster_cont[cur_claster*clasterNum+tmp]=claster_cont[tmp*clasterNum+cur_claster]=maxLev+1;
                        }
                        for(k=0; k<clasterNum; k++)//прописываем связи остальный узллов кластера, в который входит cur_claster
                        {
                            if(claster_cont[cur_claster*clasterNum+k]<=maxLev && claster_cont[cur_claster*clasterNum+k]>0)
                            {
                                for(c=0; c<mon_in_c; c++)
                                {
                                    tmp=groups_vs[c];
                                    claster_cont[k*clasterNum+tmp]=claster_cont[tmp*clasterNum+k]=maxLev+1;
                                }
                            }
                        }
                        return 0;
                    }
                }
            }
        }
    }
    return 0;

}


void move( ip, newx, newy, newz, dx, dy, dz, rnd )
long ip, newx, newy, newz, dx, dy, dz, rnd;
{
    long i, oldcell, newcell, ipm, ipp, nc1, nc2, i1, i2, ci;
    unsigned long r, rm, rp;
    double e_old, ecl_old, e_new, ecl_new,  em, ep, e_old2, ecl_old2, test;
    int forbid_move;

    int counter;
    double es1;
    double cosa;
    double xlinkl, ylinkl, zlinkl;
    double bnd_vect[12];


    long j, ib, iball[6], ncall[6];
    unsigned long rb;
    double eb;

    DDD("move:",ip, mycell[ip]);
//////////////////////////////////
// no periodic boundary conditions
    if (newx*di>=sysi) return;
    if (newx<=0) return;
    if (newy*di>=sysi) return;
    if (newy<=0) return;
    if (newz*di>=sysi) return;
    if (newz<=0) return;
//////////////////////////////////

    e_new = 0.0;
    ecl_new = 0.0;
    change = 0;
    clchange = 0.0;
    xlinkl=ylinkl=zlinkl=0;


    if (ip<Nmono)
    {
        counter=0;
        for (j=0; j<react; j++)
        {
            ib=bnd[react*ip+j];
            iball[j]=ib;
            if (ib!=-1)
            {
                xlinkl=xabs[ip]+dx-xabs[ib]; //try links
                ylinkl=yabs[ip]+dy-yabs[ib];
                zlinkl=zabs[ip]+dz-zabs[ib];
                rb = dist2(newx-x[ib], newy-y[ib], newz-z[ib]);

                if (rb>irangesq) return;

                bnd_vect[4*counter]=xlinkl*di;
                bnd_vect[4*counter+1]=ylinkl*di;
                bnd_vect[4*counter+2]=zlinkl*di;
                bnd_vect[4*counter+3]=rb*di*di;

                counter++;

                eb = enh[rb>>NB];
                who[change] = ib;
                nrg_new[change++] = eb;
                e_new += eb;
            }
        }
        if (counter==2)
        {
            cosa=(bnd_vect[0]*bnd_vect[4]+bnd_vect[1]*bnd_vect[5]+bnd_vect[2]*bnd_vect[6])/(sqrt(bnd_vect[3]*bnd_vect[7]));
            if (ip<lpm) es1=us0*(1+cosa);
            else es1=us1*(1+cosa);
            e_new += es1;
        }
        if(ucentr)
        {
            e_new += en_centrosymmetric(newx, newy, newz);
        }

        if(e_new >= Estop) return;

        for(j=0; j<react; j++)
        {
            ib=iball[j];
            if (ib!=-1) cell_get( ib, ncall[j]=mycell[ib] );
        }
    }

    cell_get(ip, oldcell=mycell[ip]);



    e_new += get_energy(newx, newy, newz);

    if (charge[ip]!=0)
    {
        for (ci=0; ci<Nall; ci++)
        {
            if (charge[ci]!=0)
            {
                if (ip!=ci)
                {
                    ecl_new+=interact_coulomb(newx,newy,newz,ip,ci);
                }
            }
        }
    }

    ecl_old = clnrg[ip];
    e_old = nrg[ip];



    /* The decision for the move */
    if(e_new/temperature + ecl_new > e_old/temperature + ecl_old)
    {
        if(rnd >= (long)(exp((e_old-e_new)/temperature+(ecl_old-ecl_new))*mask))
        {
            /* Restore and exit */
            if (ip<Nmono)
            {
                for(j=0; j<react; j++)
                {
                    ib=iball[j];
                    if (ib!=-1) cell_put( ib, ncall[j]);
                }
            }
            cell_put(ip,oldcell);
            return;
        }
    }


    /* энергия кластеров и водородных связей*/
    if(type[ip]==1)
    {
        forbid_move=0;
        if(group_cont) forbid_move+=en_of_groups(ip, newx, newy, newz, rnd);
        if(uhydro) forbid_move+=break_hydro_bnd(ip, newx, newy, newz, rnd, 1.0);
//    if(forbid_move) printf("%d-%d|",forbid_move, ip);
        if(forbid_move) /* Restore and exit */
        {
            if (ip<Nmono)
            {
                for(j=0; j<react; j++)
                {
                    ib=iball[j];
                    if (ib!=-1) cell_put( ib, ncall[j]);
                }
            }
            cell_put(ip,oldcell);
            return;
        }
    }

    /* Accept the move */
    /* First correct the new energies */
    accp += 1.0;
    nrg[ip] = e_new;
    clnrg[ip] = ecl_new;
    for(i=0; i<change; i++)
    {
        nrg[who[i]] += nrg_new[i];
    }
    for(i=0; i<clchange; i++)
    {
        clnrg[clwho[i]] += clnrg_new[i];
    }
    /* Calc corections of the old energy */
    change = 0;
    clchange = 0.0;
    e_old2 = 0.0;
    ecl_old2=0.0;

    if (ip<Nmono)
    {
        counter=0;
        for (j=0; j<react; j++)
        {
            ib=iball[j];
            if (ib!=-1)
            {
                xlinkl=xabs[ip]-xabs[ib]; //try links
                ylinkl=yabs[ip]-yabs[ib];
                zlinkl=zabs[ip]-zabs[ib];
                rb = dist2(x[ip]-x[ib], y[ip]-y[ib], z[ip]-z[ib]);

                bnd_vect[4*counter]=xlinkl*di;
                bnd_vect[4*counter+1]=ylinkl*di;
                bnd_vect[4*counter+2]=zlinkl*di;
                bnd_vect[4*counter+3]=rb*di*di;

                counter++;

                who[change] = ib;
                eb = enh[rb>>NB];
                nrg_new[change++] = eb;
                e_old2 += eb;
            }
        }
        if (counter==2)
        {
            cosa=(bnd_vect[0]*bnd_vect[4]+bnd_vect[1]*bnd_vect[5]+bnd_vect[2]*bnd_vect[6])/(sqrt(bnd_vect[3]*bnd_vect[7]));
            if (ip<lpm) es1=us0*(1+cosa);
            else es1=us1*(1+cosa);
            e_old2 += es1;
        }
        if(ucentr)
        {
            e_old2 += en_centrosymmetric(x[ip],y[ip],z[ip]);
        }


    }

    e_old2 += get_energy(x[ip],y[ip],z[ip]);

    if (charge[ip]!=0)
    {
        for (ci=0; ci<Nall; ci++)
        {
            if (charge[ci]!=0)
            {
                if (ip!=ci)
                {
                    ecl_old2 +=interact_coulomb(x[ip],y[ip],z[ip],ip,ci);
                    ecl_old2=ecl_old2;
                }
            }
        }
    }

    for(i=0; i<change; i++)
    {
        nrg[who[i]] -= nrg_new[i];
        i=i;
    }
    for(i=0; i<clchange; i++)
    {
        clnrg[clwho[i]] -= clnrg_new[i];
        i=i;
    }

    mycell[ip] = newcell = cell_number(newx,newy,newz);
    cell_put( ip, newcell );
    if (ip<Nmono)
    {
        for(j=0; j<react; j++)
        {
            ib=iball[j];
            if (ib!=-1) cell_put( ib, ncall[j]);
        }
    }
    x[ip]=newx;
    y[ip]=newy;
    z[ip]=newz;
    xabs[ip]+=(double)dx;
    yabs[ip]+=(double)dy;
    zabs[ip]+=(double)dz;
    return;
}



void tryjump( ip, jp, rnd)
long ip, jp, rnd;
{
    long i, ci;
    double ecl_old, ecl_new, ecl_old2, ecl_new_i, ecl_new_j, ecl_old_i, ecl_old_j;
    int accp_jump=0, tmp;

    ecl_new = 0.0;
    clchange = 0.0;

    tmp=charge[jp];
    charge[jp]=charge[ip];
    charge[ip]=tmp;

    if (charge[ip]!=0)
    {
        for (ci=0; ci<Nall; ci++)
        {
            if (charge[ci]!=0)
            {
                if (ip!=ci)
                {
                    ecl_new_i+=interact_coulomb(x[ip],y[ip],z[ip],ip,ci);
                }
            }
        }
    }
    if (charge[jp]!=0)
    {
        for (ci=0; ci<Nall; ci++)
        {
            if (charge[ci]!=0)
            {
                if (jp!=ci)
                {
                    ecl_new_j+=interact_coulomb(x[jp],y[jp],z[jp],jp,ci);
                }
            }
        }
    }
    ecl_old_i = clnrg[ip];
    ecl_old_j = clnrg[jp];

    ecl_new=ecl_new_i+ecl_new_j;
    ecl_old=ecl_old_i+ecl_old_j;

    /* The decision for the move */
    if(ecl_new > ecl_old)
    {
        if(rnd >= (long)(exp(ecl_old-ecl_new)*mask))
        {
            /* Restore and exit */
            tmp=charge[jp];
            charge[jp]=charge[ip];
            charge[ip]=tmp;

            return;
        }
    }
    /* Accept the move */
    /* First correct the new energies */
    accp_jump += 1.0;

    clnrg[ip] = ecl_new_i;
    clnrg[jp] = ecl_new_j;
    for(i=0; i<clchange; i++)
    {
        clnrg[clwho[i]] += clnrg_new[i];
    }
    /* Calc corections of the old energy */
    clchange = 0.0;
    ecl_old2=0.0;

    if (charge[ip]!=0)
    {
        for (ci=0; ci<Nall; ci++)
        {
            if (charge[ci]!=0)
            {
                if (ip!=ci)
                {
                    ecl_old2 +=interact_coulomb(x[ip],y[ip],z[ip],ip,ci);
                }
            }
        }
    }

    if (charge[jp]!=0)
    {
        for (ci=0; ci<Nall; ci++)
        {
            if (charge[ci]!=0)
            {
                if (jp!=ci)
                {
                    ecl_old2 +=interact_coulomb(x[jp],y[jp],z[jp],jp,ci);
                }
            }
        }
    }

    for(i=0; i<clchange; i++)
    {
        clnrg[clwho[i]] -= clnrg_new[i];
    }

    return;
}



float interact_coulomb(ix,iy,iz,ip,io)
long ix, iy, iz, ip, io;
{
    float cle=0;
    long dx, dy, dz, dr;
    double inv_dr, inv_dr1;
    int ci;

    dr=di*(ix-x[io])*(ix-x[io])+di*(iy-y[io])*(iy-y[io])+di*(iz-z[io])*(iz-z[io]);
    inv_dr=1/sqrt(dr*di);
    /*    if(dr<lim)
            inv_dr=enc[dr>>NBC];
        else
            inv_dr=enc2[dr>>NBC2];
    */
    clwho[clchange]=io;
    cle = clnrg_new[clchange++] = bjerrum*charge[ip]*charge[io]*inv_dr;
    return (cle);
}


double get_energy( ix, iy, iz )
long ix, iy, iz;
{
    struct part *q;
    long nc,  nx,  ny,  nz,  sx,  sy;
    long npx, nmx, npy, nmy, npz, nmz;
    long npp, npm, nmp, nmm, nxy;
    register unsigned long dx, np;
    double e;

    e = 0.0;
    nx = ix>>irangeL;
    ny = iy>>irangeL;
    nz = iz>>irangeL;
    sx = cellL<<1;
    sy=cellL;

    nmm = (cellm[nx]<<sx)|(cellm[ny]<<sy);
    nmp = (cellm[nx]<<sx)|(cellp[ny]<<sy);
    npm = (cellp[nx]<<sx)|(cellm[ny]<<sy);
    npp = (cellp[nx]<<sx)|(cellp[ny]<<sy);
    nmy = (cellm[nx]<<sx)|(ny<<sy);
    npy = (cellp[nx]<<sx)|(ny<<sy);
    nx = nx << sx;
    nmx = nx|(cellm[ny]<<sy);
    npx = nx|(cellp[ny]<<sy);
    nxy = nx|(ny<<sy);
    nmz = cellm[nz];
    npz = cellp[nz];

    nc=nmm|nmz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nmm|nz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nmm|npz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nmy|nmz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nmy|nz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nmy|npz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nmp|nmz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nmp|nz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nmp|npz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nmx|nmz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nmx|nz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nmx|npz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nxy|nmz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nxy|nz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=nxy|npz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=npx|nmz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=npx|nz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=npx|npz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=npm|nmz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=npm|nz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=npm|npz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=npy|nmz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=npy|nz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=npy|npz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=npp|nmz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=npp|nz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    nc=npp|npz;
    for(q=cell[nc].next; q!=NULL; q=q->next)
    {
        np = q->pn;
        dx = dist2( ix-x[np], iy-y[np], iz-z[np]);
        if(dx>=irangesq) continue;
        who[change]=np;
        e += nrg_new[change++] = enp[dx>>NB];
    }
    return(e);
}

double dist_real(long mon_x,long mon_y,long mon_z )
{
    return sqrt(di*di*mon_x*mon_x+di*di*mon_y*mon_y+di*di*mon_z*mon_z);
}
