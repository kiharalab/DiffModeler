/* Run-time library for use with "p2c", the Pascal to C translator */

/* "p2c"  Copyright (C) 1989, 1990, 1991 Free Software Foundation.
 * By Dave Gillespie, daveg@csvax.cs.caltech.edu.  Version --VERSION--.
 * This file may be copied, modified, etc. in any way.  It is not restricted
 * by the licence agreement accompanying p2c itself.
 */



#include "p2c.h"

/* Check if at end of line (or end of entire file). */

int P_eoln(f)
FILE *f;
{
    register int ch;

    ch = getc(f);
    if (ch == EOF)
        return 1;
    ungetc(ch, f);
    return (ch == '\n');
}

int P_eof(f)
FILE *f;
{
    register int ch;

    if (feof(f))
	return 1;
    if (f == stdin)
	return 0;    /* not safe to look-ahead on the keyboard! */
    ch = getc(f);
    if (ch == EOF)
	return 1;
    ungetc(ch, f);
    return 0;
}


/* Sets are stored as an array of longs.  S[0] is the size of the set;
   S[N] is the N'th 32-bit chunk of the set.  S[0] equals the maximum
   I such that S[I] is nonzero.  S[0] is zero for an empty set.  Within
   each long, bits are packed from lsb to msb.  The first bit of the
   set is the element with ordinal value 0.  (Thus, for a "set of 5..99",
   the lowest five bits of the first long are unused and always zero.) */

/* (Sets with 32 or fewer elements are normally stored as plain longs.) */

int P_setequal(s1, s2)              /* s1 = s2 */
register long *s1, *s2;
{
    register int size = *s1++;
    if (*s2++ != size)
        return 0;
    while (--size >= 0) {
        if (*s1++ != *s2++)
            return 0;
    }
    return 1;
}



long *P_setunion(d, s1, s2)         /* d := s1 + s2 */
register long *d, *s1, *s2;
{
    long *dbase = d++;
    register int sz1 = *s1++, sz2 = *s2++;
    while (sz1 > 0 && sz2 > 0) {
        *d++ = *s1++ | *s2++;
	sz1--, sz2--;
    }
    while (--sz1 >= 0)
	*d++ = *s1++;
    while (--sz2 >= 0)
	*d++ = *s2++;
    *dbase = d - dbase - 1;
    return dbase;
}

int P_inset(val, s)                 /* val IN s */
register unsigned val;
register long *s;
{
    register int bit;
    bit = val % SETBITS;
    val /= SETBITS;
    if (val < *s++ && ((1<<bit) & s[val]))
	return 1;
    return 0;
}


long *P_addset(s, val)              /* s := s + [val] */
register long *s;
register unsigned val;
{
    register long *sbase = s;
    register int bit, size;
    bit = val % SETBITS;
    val /= SETBITS;
    size = *s;
    if (++val > size) {
        s += size;
        while (val > size)
            *++s = 0, size++;
        *sbase = size;
    } else
        s += val;
    *s |= 1<<bit;
    return sbase;
}

/* s is a "smallset", i.e., a 32-bit or less set stored
   directly in a long. */

long *P_expset(d, s)                /* d := s */
register long *d;
register long s;
{
    if (s) {
	d[1] = s;
	*d = 1;
    } else
        *d = 0;
    return d;
}


long *P_setdiff(d, s1, s2)          /* d := s1 - s2 */
register long *d, *s1, *s2;
{
    long *dbase = d++;
    register int sz1 = *s1++, sz2 = *s2++;
    while (--sz1 >= 0 && --sz2 >= 0)
        *d++ = *s1++ & ~*s2++;
    if (sz1 >= 0) {
        while (sz1-- >= 0)
            *d++ = *s1++;
    }
    while (--d > dbase && !*d) ;
    *dbase = d - dbase;
    return dbase;
}


long *P_setcpy(d, s)                /* d := s */
register long *d, *s;
{
    register long *save_d = d;

#ifdef SETCPY_MEMCPY
    memcpy(d, s, (*s + 1) * sizeof(long));
#else
    register int i = *s + 1;
    while (--i >= 0)
        *d++ = *s++;
#endif
    return save_d;
}

#ifdef NEED_MY_MEMCPY
/* The original version of my_memcpy crashed on all computers where the last
   function argument type "size_t" was unsigned, due to a bug in the while
   loop which originally read "while (--n>=0)", which is of course always
   true for unsigned values. Fixed by elmar.krieger@cmbi.kun.nl */

Anyptr my_memcpy(Anyptr d, Anyptr s, size_t n)
{
    register char *ss = (char *)s, *dd = (char *)d;
    while (n--!=0) *(dd++) = *(ss++);
    return d;
}
#endif
