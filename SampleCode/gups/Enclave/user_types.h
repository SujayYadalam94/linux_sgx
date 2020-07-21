#ifndef _USER_TYPES_H_
#define _USER_TYPES_H_

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

#endif
