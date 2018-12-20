#include "Enclave_t.h"
#include "stdio.h"
#include "stdlib.h"

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

int generate_random_number() {
	char *p;
	p = (char *)malloc(100);
	enclave_printf("Address of p = 0x%lx\n", p);
	enclave_printf("Generating random number\n");

    return 42;
}
