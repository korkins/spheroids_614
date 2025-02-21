#include <stdio.h>  /* printf */
#include <time.h>   /* clock  */
#include "const_param_spheroids.h"
//
int optichar(double const wavel, double const refre, double const refim,
	        double const cvf, double const radf, double const sgmf,
	        double const cvc, double const radc, double const sgmc,
	        double *f11, double *f22, double *f33, double *f44, double *f12, double *f34,
			double &cext, double &ssa, double &x1, double &conc_ratio);
//
int main()
{
	FILE
		*file;
	char const
		*fname_inp = "input.dat",
		*fname_out = "aerosols.phs";
	char
		chr;
	int
		code, error_code, isca, result;
	double
		wavel, refre, refim, time_start, time_end,
		cvf, cvc, radf, radc, sgmf, sgmc, conc_ratio, cext, ssa, x1, f1;
    double
		*f11, *f22, *f33, *f44, *f12, *f34;
//
	time_start = (double)clock() /(double)CLOCKS_PER_SEC;
	error_code = -1; // default: missing file
//
	f11 = new double [nsca_fix];
	f22 = new double [nsca_fix];
	f33 = new double [nsca_fix];
	f44 = new double [nsca_fix];
	f12 = new double [nsca_fix];
	f34 = new double [nsca_fix];
	for (isca = 0; isca < nsca_fix; isca++)
	{
		f11[isca] = 0.0;
		f22[isca] = 0.0;
		f33[isca] = 0.0;
		f44[isca] = 0.0;
		f12[isca] = 0.0;
		f34[isca] = 0.0;
	} // for isca
//
	file = fopen(fname_inp, "r");
	if (file == NULL) // if file is missing
		printf("\nError opening file %s - execution terminated", fname_inp);
	else
	{
//      Step 1: read input
		printf("\nReading input from %s ... ", fname_inp);
//		line 1: skip
		chr = 'X'; while(chr != '\n') result = fscanf(file, "%c", &chr);
//		line 2: read 3 parameters and skip the rest
		result = fscanf(file, "%lf %lf %lf%c", &wavel, &refre, &refim, &chr);
		while(chr != '\n') result = fscanf(file, "%c", &chr);
//		line 3: skip
		chr = 'X'; while(chr != '\n') result = fscanf(file, "%c", &chr);
//		line 4: read 4 parameters
		result = fscanf(file, "%lf %lf %lf %lf %lf %lf%c", &cvf, &sgmf, &radf, &cvc, &sgmc, &radc, &chr);
		printf("done!");
		fclose(file);
		printf("\nFile %s is closed", fname_inp);
//
//      Step 2: run calculation
	    code = optichar(wavel, refre, refim, cvf, radf, sgmf, cvc, radc, sgmc,
					    f11, f22, f33, f44, f12, f34, cext, ssa, x1, conc_ratio);
//
//      Step 3: save to file
		file = fopen(fname_out, "w");
		printf("\nSaving output to %s ... ", fname_out);
//
		fprintf(file, "IFLAG =   %d\n", 7);
        fprintf(file, "Created using output from refactored DLS code\n");
        fprintf(file, "Line 3\n");
        fprintf(file, "LAM=  %.6f   MRR=%.4e   MRI=%.4e\n", wavel, refre, refim);
        fprintf(file, "Line 5\n");
        fprintf(file, "Line 6\n");
        fprintf(file, "R1=  %.6e   R2= %.6e\n", rmin_grid_fix, rmax_grid_fix);
        fprintf(file, "Line 8\n");
        fprintf(file, "Line 9\n");
        fprintf(file, "Line 10\n");
        fprintf(file, "CEXT=%.5e  CSCA=%.5e    W=%.5e  <COS>=%.5e\n", cext, cext*ssa, ssa, x1);
		fprintf(file, "\n");
		fprintf(file, "      <        F11        F22        F33        F44        F12        F34\n");
		for (isca = 0; isca < nsca_fix; isca++)
        {
			f1 = f11[isca];
			fprintf(file, "%6.2f %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n", sca_fix[isca],
				    f1, f1*f22[isca], f1*f33[isca], f1*f44[isca], -f1*f12[isca], f1*f34[isca]); // note '-' at f12
        } // for isca
		printf("done!");
		fclose(file);
		printf("\nFile %s is closed", fname_out);
//
		error_code = 0; // no errors
	} // end of reading file
//

	time_end = (double)clock() /(double)CLOCKS_PER_SEC;
	printf("\nruntime per run: %8.3f (sec.)", time_end - time_start);
//
	printf("\n");
	return error_code;
}
/*------------------------------------------------------------------------------------------------*/