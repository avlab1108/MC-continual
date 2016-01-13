/* The RANDOM GENERATOR i127 */
#include <stdlib.h>

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

/* Declaration for Random Generators
void ini127( long  *m, long iseed);
void i127( long  *m, long  *ran, long  n);
*/

void ini127( long *m, long iseed)
{
    long  i, ib, x;
    long  mask[32]= {0x1L, 0x2L, 0x4L, 0x8L, 0x10L, 0x20L, 0x40L, 0x80L,
                     0x100L,  0x200L, 0x400L, 0x800L, 0x1000L, 0x2000L, 0x4000L, 0x8000L,
                     0x10000L,0x20000L,0x40000L,0x80000L,0x100000L,0x200000L,0x400000L,0x800000L,
                     0x1000000L, 0x2000000L, 0x4000000L, 0x8000000L, 0x10000000L, 0x20000000L,
                     0x40000000L, 0x80000000L
                    };

    srand(iseed);

    m[0] = 1;
    for(i=1; i<=127; i++)
    {
        for(ib=0; ib<=31; ib++)
        {
            x = (rand() >> 8) & 1;
            if( x )
                m[i] |= mask[ib];
            else
                m[i] &= (mask[ib] ^ ( long)0xFFFFffffL);
        }
    }
}

/* NEW FAST VERSION of 250 */
void i127( long  *m, long  *ran, long  n )
{
    long  k, kmin, kmax, l;
    long  *mp, m0 = *m;

    mp = m;
    l = n+m0-1;
huj:
    m = mp + m0;
    kmin = m0;
    kmax = min(l,97) ;
    for(k=kmin; k<=kmax; k++, ran++, m++)
        *ran= *m^=*(m+30);
    kmin = max(kmax+1,m0);
    kmax = min(l,127);
    for(k=kmin; k<=kmax; k++,ran++, m++)
        *ran=*m^=*(m-97);
    if( kmax == 127 )
    {
        m0 = 1;
        l = l - 127;
        goto huj;
    }
    *mp = kmax + 1;
}
