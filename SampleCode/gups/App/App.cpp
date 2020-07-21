#include <stdio.h>
#include <iostream>
#include "Enclave_u.h"
#include "sgx_urts.h"

#include <stdlib.h>
#include <mpich/mpi.h>
#include <unistd.h>
#include <signal.h>

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

int main(int narg, char *arg[])
{
	int me, nprocs;
	int i, niterate;
	int nglobal, nlocal, logtable, logtablelocal;
	int logprocs, maxndata, maxnfinal, nexcess;
	int nbad, chunk;
	double t0, t0_all, Gups;
	u64Int offset,nupdates;
	u64Int ilong,nexcess_long,nbad_long;
	sgx_status_t ret;
	
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
	
	if (initialize_enclave(&global_eid, "enclave.token", "enclave.signed.so") < 0)
	{
		std::cout << "Fail to initialize enclave." << std::endl;
		sgx_destroy_enclave(global_eid);
		return 1;
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

	nglobal = ((int) 1) << logtable;
	nlocal = nglobal / nprocs;
	logtablelocal = logtable - logprocs;
	offset = (u64Int) nlocal * me;

	int pid=getpid();
	int cpid=fork();

	if(cpid == 0)
	{
		//Child process
		char buf[300];
		sprintf(buf, "perf stat -e dTLB-loads:u,dTLB-stores:u,dTLB-load-misses:u,dTLB-store-misses:u,dtlb_load_misses.walk_active:u,dtlb_store_misses.walk_active:u,dtlb_load_misses.walk_completed_2m_4m:u,dtlb_load_misses.walk_completed_4k:u,dtlb_store_misses.walk_completed_2m_4m:u,dtlb_store_misses.walk_completed_4k:u -p %d",pid);
		execl("/bin/sh", "sh", "-c", buf, NULL);
	}
	else
	{
	syscall(333, getpid());
	setpgid(cpid, 0);
	t0 = -MPI_Wtime();
	ret = ecall_gups_loop(global_eid, NULL, nlocal, chunk, offset, nprocs, logprocs, niterate, me, logtablelocal, &maxndata, &maxnfinal, &nexcess, &nbad);
	syscall(333, -1);
	MPI_Barrier(MPI_COMM_WORLD);
	t0 += MPI_Wtime();
	if(ret)
	{
		printf("Error in enclave execution)\n");
		sgx_destroy_enclave(global_eid);
		return 1;
	}

	kill(-cpid, SIGINT);

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
	Gups = (double)nupdates / t0 / 1.0e9;

	if (me == 0) {
		printf("Number of procs: %d\n",nprocs);
		printf("Vector size: %d\n",nglobal);
		printf("Max datums during comm: %d\n",maxndata);
		printf("Max datums after comm: %d\n",maxnfinal);
		printf("Excess datums (frac): %lld (%g)\n",
				nexcess_long,(double) nexcess_long / nupdates);
		printf("Bad locality count: %lld\n",nbad_long);
		printf("Update time (secs): %9.3f\n",t0);
		printf("Gups: %9.6f\n",Gups);
	}

	/* clean up */
	MPI_Finalize();
	
	printf("Info: Enclave successfully returned.\n");
	scanf("%d",me);

	sgx_destroy_enclave(global_eid);
	return 0;
	}
}
