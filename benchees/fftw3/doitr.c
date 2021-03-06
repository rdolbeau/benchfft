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
BENCH_DOC("name", "fftw3-r2r")
BENCH_DOCF("version", mkversion)
BENCH_DOCF("fftw-cc", mkcc)
BENCH_DOCF("fftw-codelet-optim", mkcodelet_optim)
BENCH_DOC("package", "FFTW 3")
BENCH_DOC("year", "2003")
BENCH_DOC("author", "Matteo Frigo")
BENCH_DOC("author", "Steven G. Johnson")
BENCH_DOC("email", "fftw@fftw.org")
BENCH_DOC("url", "http://www.fftw.org")
BENCH_DOC("url-was-valid-on", "Fri Mar 28 18:46:22 EST 2003")
BENCH_DOC("language", "C")
BENCH_DOC("language", "Objective Caml")
BENCH_DOC("notes", "Using halfcomplex transform interface (R2HC/HC2R kinds).")
BENCH_DOC("copyright",
"Copyright (c) 2003 Matteo Frigo\n"
"Copyright (c) 2003 Massachusetts Institute of Technology\n"
"\n"
"This program is free software; you can redistribute it and/or modify\n"
"it under the terms of the GNU General Public License as published by\n"
"the Free Software Foundation; either version 2 of the License, or\n"
"(at your option) any later version.\n"
"\n"
"This program is distributed in the hope that it will be useful,\n"
"but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"GNU General Public License for more details.\n"
"\n"
"You should have received a copy of the GNU General Public License\n"
"along with this program; if not, write to the Free Software\n"
"Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n")
BENCH_DOC("bibitem", "M. Frigo and S. G. Johnson, FFTW: An adaptive software architecture for the FFT, Proc. ICASSP 3, 1381-1384 (1998)")
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
