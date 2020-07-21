#include "Enclave_t.h"
#include "stdio.h"
#include "stdlib.h"
#include <user_types.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
/* machine defs
   compile with -DLONG64 if a "long" is 64 bits
   else compile with no setting if "long long" is 64 bit */

/* 
 * printf: 
 *   Invokes OCALL to display the enclave buffer to the terminal.
 */
void enclave_printf(const char *fmt, ...)
{
    char buf[BUFSIZ] = {'\0'};
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, BUFSIZ, fmt, ap);
    va_end(ap);
    ocall_print(buf);
}

/* start random number generator at Nth step of stream
   routine provided by HPCC */

u64Int HPCC_starts(long n)
{
  int i, j;
  u64Int m2[64];
  u64Int temp, ran;

  while (n < 0) n += PERIOD;
  while (n > PERIOD) n -= PERIOD;
  if (n == 0) return 0x1;

  temp = 0x1;
  for (i=0; i<64; i++) {
    m2[i] = temp;
    temp = (temp << 1) ^ ((s64Int) temp < 0 ? POLY : 0);
    temp = (temp << 1) ^ ((s64Int) temp < 0 ? POLY : 0);
  }

  for (i=62; i>=0; i--)
    if ((n >> i) & 1)
      break;

  ran = 0x2;
  while (i > 0) {
    temp = 0;
    for (j=0; j<64; j++)
      if ((ran >> j) & 1)
        temp ^= m2[j];
    ran = temp;
    i -= 1;
    if ((n >> i) & 1)
      ran = (ran << 1) ^ ((s64Int) ran < 0 ? POLY : 0);
  }

  return ran;
}

int ecall_gups_loop(int nlocal, int chunk, u64Int offset, int nprocs, int logprocs, int niterate, 				int me, int logtablelocal, int *maxndata, int *maxnfinal, int *nexcess, int *nbad)
{
	u64Int *table, *data, *send;
	int i, j, iterate, index;
	u64Int nupdates;
	u64Int ran;
	int ipartner, nsend, nkeep, nrecv, ndata;
	u64Int datum, procmask;
	int nlocalm1 = nlocal -1;
	
	char *dummy = (char *) malloc(0x1c0F60); //Padding to get large page aligned
	table = (u64Int *) malloc(nlocal*sizeof(u64Int));
	data = (u64Int *) malloc(chunk*16*sizeof(u64Int));
	send = (u64Int *) malloc(chunk*16*sizeof(u64Int));

	if (!table || !data || !send) {
			enclave_printf("Table allocation failed\n");
		return 1;
	}	
	enclave_printf("Pointer addresses = %x %x %x \n", table, data, send);

	/* initialize my portion of global array
	   global array starts with table[i] = i */
	for (i = 0; i < nlocal; i++)
		table[i] = i + offset;

	/* start my random # partway thru global stream */
	nupdates = (u64Int) nprocs * chunk * niterate;
	ran = HPCC_starts(nupdates/nprocs*me);

	/* loop:
	   generate chunk random values per proc
	   communicate datums to correct processor via hypercube routing
	   use received values to update local table */

	*maxndata = 0;
	*maxnfinal = 0;
	*nexcess = 0;
	*nbad = 0;

	for (iterate = 0; iterate < niterate; iterate++)
	{
		for (i = 0; i < chunk; i++)
		{
			ran = (ran << 1) ^ ((s64Int) ran < ZERO64B ? POLY : ZERO64B);
			data[i] = ran;
		}
		ndata = chunk;

		for (j = 0; j < logprocs; j++)
		{
			nkeep = nsend = 0;
			ipartner = (1 << j) ^ me;
			procmask = ((u64Int) 1) << (logtablelocal + j);
			if (ipartner > me)
			{
				for (i = 0; i < ndata; i++)
				{
					if (data[i] & procmask) send[nsend++] = data[i];
					else data[nkeep++] = data[i];
				}
			} else {
				for (i = 0; i < ndata; i++)
				{
					if (data[i] & procmask) data[nkeep++] = data[i];
					else send[nsend++] = data[i];
				}
			}
			ndata = nkeep + nrecv;
			*maxndata = MAX(*maxndata, ndata);
		}
		*maxnfinal = MAX(*maxnfinal,ndata);
		if (ndata > chunk) 
			*nexcess += ndata - chunk;

		for (i = 0; i < ndata; i++)
		{
			datum = data[i];
			index = datum & nlocalm1;
			table[index] ^= datum;
		}

#ifdef CHECK
		procmask = ((u64Int) (nprocs-1)) << logtablelocal;
		for (i = 0; i < ndata; i++)
			if ((data[i] & procmask) >> logtablelocal != me) 
				*nbad += 1;
#endif
	}
	
	free(table);
	free(data);
	free(send);
	
	return 0;
}
