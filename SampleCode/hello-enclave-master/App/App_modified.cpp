#include <stdio.h>
#include <iostream>
#include "Enclave_u.h"
#include "sgx_urts.h"

#include <stdlib.h>
#include <mpich/mpi.h>
#include <unistd.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* machine defs
   compile with -DLONG64 if a "long" is 64 bits
   else compile with no setting if "long long" is 64 bit */

#ifdef LONG64
#define POLY 0x0000000000000007UL
#define PERIOD 1317624576693539401L
#define ZERO64B 0L
typedef long s64Int;
typedef unsigned long u64Int;
#define U64INT MPI_UNSIGNED_LONG
#else
#define POLY 0x0000000000000007ULL
#define PERIOD 1317624576693539401LL
#define ZERO64B 0LL
typedef long long s64Int;
typedef unsigned long long u64Int;
#define U64INT MPI_LONG_LONG_INT
#endif

/* Global EID shared by multiple threads */
sgx_enclave_id_t global_eid = 0;

// OCall implementations
void ocall_print(const char* str) {
    printf("%s\n", str);
}

/* Check error conditions for loading enclave */
void print_error_message(sgx_status_t ret) {
    printf("SGX error code: %d\n", ret);
}

int initialize_enclave(sgx_enclave_id_t* eid, const std::string& launch_token_path, const std::string& enclave_name) {
    const char* token_path = launch_token_path.c_str();
    sgx_launch_token_t token = {0};
    sgx_status_t ret = SGX_ERROR_UNEXPECTED;
    int updated = 0;

    /* Step 1: try to retrieve the launch token saved by last transaction
     *         if there is no token, then create a new one.
     */
    /* try to get the token saved in $HOME */
    FILE* fp = fopen(token_path, "rb");
    if (fp == NULL && (fp = fopen(token_path, "wb")) == NULL) {
        printf("Warning: Failed to create/open the launch token file \"%s\".\n", token_path);
    }

    if (fp != NULL) {
        /* read the token from saved file */
        size_t read_num = fread(token, 1, sizeof(sgx_launch_token_t), fp);
        if (read_num != 0 && read_num != sizeof(sgx_launch_token_t)) {
            /* if token is invalid, clear the buffer */
            memset(&token, 0x0, sizeof(sgx_launch_token_t));
            printf("Warning: Invalid launch token read from \"%s\".\n", token_path);
        }
    }
    /* Step 2: call sgx_create_enclave to initialize an enclave instance */
    /* Debug Support: set 2nd parameter to 1 */
    ret = sgx_create_enclave(enclave_name.c_str(), SGX_DEBUG_FLAG, &token, &updated, eid, NULL);
    if (ret != SGX_SUCCESS) {
        print_error_message(ret);
        if (fp != NULL) fclose(fp);
        return -1;
    }

    /* Step 3: save the launch token if it is updated */
    if (updated == false || fp == NULL) {
        /* if the token is not updated, or file handler is invalid, do not perform saving */
        if (fp != NULL) fclose(fp);
        return 0;
    }

    /* reopen the file with write capablity */
    fp = freopen(token_path, "wb", fp);
    if (fp == NULL) return 0;
    size_t write_num = fwrite(token, 1, sizeof(sgx_launch_token_t), fp);
    if (write_num != sizeof(sgx_launch_token_t))
        printf("Warning: Failed to save launch token to \"%s\".\n", token_path);
    fclose(fp);
    return 0;
}

/* start random number generator at Nth step of stream
   routine provided by HPCC */

u64Int HPCC_starts(s64Int n)
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

int main(int narg, char *arg[])
{
	if (initialize_enclave(&global_eid, "enclave.token", "enclave.signed.so") < 0)
	{
		std::cout << "Fail to initialize enclave." << std::endl;
		return 1;
	}
	int ptr;
	sgx_status_t ret = generate_random_number(global_eid, &ptr);
	std::cout << ret << std::endl;
	if (ret != SGX_SUCCESS)
	{
		std::cout << "noob" << std::endl;
	}
	printf("Random number: %d\n", ptr);

	sgx_destroy_enclave(global_eid);
	printf("Info: SampleEnclave successfully returned.\n");



	int me,nprocs;
	int i,j,iterate,niterate;
	int nlocal,nlocalm1,logtable,index,logtablelocal;
	int logprocs,ipartner,ndata,nsend,nkeep,nrecv,maxndata,maxnfinal,nexcess;
	int nbad,chunk,chunkbig;
	double t0,t0_all,Gups;
	u64Int *table,*data,*send;
	u64Int ran,datum,procmask,nglobal,offset,nupdates;
	u64Int ilong,nexcess_long,nbad_long;
	MPI_Status status;

	MPI_Init(&narg,&arg);
	MPI_Comm_rank(MPI_COMM_WORLD,&me);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	/* command line args = N M chunk
	   N = length of global table is 2^N
	   M = # of update sets per proc
	   chunk = # of updates in one set */

	if (narg != 4) {
		if (me == 0) printf("Syntax: gups N M chunk\n");
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	logtable = atoi(arg[1]);
	niterate = atoi(arg[2]);
	chunk = atoi(arg[3]);

	/* insure Nprocs is power of 2 */

	i = 1;
	while (i < nprocs) i *= 2;
	if (i != nprocs) {
		if (me == 0) printf("Must run on power-of-2 procs\n");
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	/* nglobal = entire table
	   nlocal = size of my portion
	   nlocalm1 = local size - 1 (for index computation)
	   logtablelocal = log of table size I store
	   offset = starting index in global table of 1st entry in local table */

	logprocs = 0;
	while (1 << logprocs < nprocs) logprocs++;

	nglobal = ((u64Int) 1) << logtable;
	nlocal = nglobal / nprocs;
	nlocalm1 = nlocal - 1;
	logtablelocal = logtable - logprocs;
	offset = (u64Int) nlocal * me;

	/* allocate local memory
	   16 factor insures space for messages that (randomly) exceed chunk size */

	chunkbig = 16*chunk;

	table = (u64Int *) malloc(nlocal*sizeof(u64Int));
	data = (u64Int *) malloc(chunkbig*sizeof(u64Int));
	send = (u64Int *) malloc(chunkbig*sizeof(u64Int));

	if (!table || !data || !send) {
		if (me == 0) printf("Table allocation failed\n");
		MPI_Abort(MPI_COMM_WORLD,1);
	}	

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

	maxndata = 0;
	maxnfinal = 0;
	nexcess = 0;
	nbad = 0;

	MPI_Barrier(MPI_COMM_WORLD);
	t0 = -MPI_Wtime();

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

			MPI_Sendrecv(send,nsend,U64INT,ipartner,0,&data[nkeep],chunkbig,U64INT,
					ipartner,0,MPI_COMM_WORLD,&status);
			MPI_Get_count(&status,U64INT,&nrecv);
			ndata = nkeep + nrecv;
			maxndata = MAX(maxndata,ndata);
		}
		maxnfinal = MAX(maxnfinal,ndata);
		if (ndata > chunk) nexcess += ndata - chunk;

		for (i = 0; i < ndata; i++)
		{
			datum = data[i];
			index = datum & nlocalm1;
			table[index] ^= datum;
		}

#ifdef CHECK
		procmask = ((u64Int) (nprocs-1)) << logtablelocal;
		for (i = 0; i < ndata; i++)
			if ((data[i] & procmask) >> logtablelocal != me) nbad++;
#endif
	}

	MPI_Barrier(MPI_COMM_WORLD);
	t0 += MPI_Wtime();

	/* stats */

	MPI_Allreduce(&t0,&t0_all,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	t0 = t0_all/nprocs;

	i = maxndata;
	MPI_Allreduce(&i,&maxndata,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
	i = maxnfinal;
	MPI_Allreduce(&i,&maxnfinal,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
	ilong = nexcess;
	MPI_Allreduce(&ilong,&nexcess_long,1,U64INT,MPI_SUM,MPI_COMM_WORLD);
	ilong = nbad;
	MPI_Allreduce(&ilong,&nbad_long,1,U64INT,MPI_SUM,MPI_COMM_WORLD);

	nupdates = (u64Int) niterate * nprocs * chunk;
	Gups = nupdates / t0 / 1.0e9;

	if (me == 0) {
		printf("Number of procs: %d\n",nprocs);
		printf("Vector size: %lld\n",nglobal);
		printf("Max datums during comm: %d\n",maxndata);
		printf("Max datums after comm: %d\n",maxnfinal);
		printf("Excess datums (frac): %lld (%g)\n",
				nexcess_long,(double) nexcess_long / nupdates);
		printf("Bad locality count: %lld\n",nbad_long);
		printf("Update time (secs): %9.3f\n",t0);
		printf("Gups: %9.6f\n",Gups);
	}

	/* clean up */

	free(table);
	free(data);
	free(send);
	MPI_Finalize();

	return 0;
}
