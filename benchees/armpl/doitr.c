#include "bench-user.h"
#include <math.h>
#include <stdio.h>
#include <fftw3.h>
#include <string.h>

#define CONCAT(prefix, name) prefix ## name
#if defined(BENCHFFT_SINGLE)
#define FFTW(x) CONCAT(fftwf_, x)
#elif defined(BENCHFFT_LDOUBLE)
#define FFTW(x) CONCAT(fftwl_, x)
#else
#define FFTW(x) CONCAT(fftw_, x)
#endif

static const char *mkversion(void) { return FFTW(version); }
static const char *mkcc(void) { return FFTW(cc); }
static const char *mkcodelet_optim(void) { return FFTW(codelet_optim); }

BEGIN_BENCH_DOC
BENCH_DOC("name", "armpl")
BENCH_DOC("package", "ARMPL")
BENCH_DOC("year", "2017")
BENCH_DOC("author", "ARM")
BENCH_DOC("email", "support-hpc-sw@arm.com")
END_BENCH_DOC 

FFTW(plan) the_plan = 0;
unsigned the_flags = 0;

void useropt(const char *arg)
{
     if (!strcmp(arg, "patient")) the_flags |= FFTW_PATIENT;
     else if (!strcmp(arg, "estimate")) the_flags |= FFTW_ESTIMATE;
     else if (!strcmp(arg, "exhaustive")) the_flags |= FFTW_EXHAUSTIVE;
     else if (!strcmp(arg, "unaligned")) the_flags |= FFTW_UNALIGNED;

     else fprintf(stderr, "unknown user option: %s.  Ignoring.\n", arg);
}

int can_do(struct problem *p)
{
     return (p->rank == 1 && p->kind == PROBLEM_REAL);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_halfcomplex(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_halfcomplex(p, in, -1.0);
}

void setup(struct problem *p)
{
     unsigned flags = the_flags;  
     BENCH_ASSERT(can_do(p));
 
     if (p->sign == 1)
	  the_flags |= FFTW_DESTROY_INPUT;

     the_plan = FFTW(plan_r2r_1d)(p->n[0], p->in, p->out,
				  (p->sign == -1 ? FFTW_R2HC : FFTW_HC2R),
				  the_flags);

     BENCH_ASSERT(the_plan);
     
     if (verbose >= 2) {
	  FFTW(print_plan)(the_plan);
	  printf("\n");
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     FFTW(plan) plan = the_plan;

     UNUSED(p);

     for (i = 0; i < iter; ++i) 
	  FFTW(execute)(plan);
}

void done(struct problem *p)
{
     UNUSED(p);
     FFTW(destroy_plan)(the_plan);
}
