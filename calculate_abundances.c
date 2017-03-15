/*=======================================================================

 Calculate the abundances of all species at the specified end time
 based on their initial abundances and the rates for each reaction.
 This routine calls the CVODE package to solve for the set of ODEs.
 CVODE is able to handle stiff problems, where the dynamic range of
 the rates can be very large.

-----------------------------------------------------------------------*/

/* Standard includes */

#ifdef OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

/*-----------------------------------------------------------------------*/

/* Header files with descriptions of the contents used */

#include <cvode/cvode.h>             /* CVODE functions and constants */
#include <cvode/cvode_dense.h>       /* Prototype for CVDense solver */
#include <nvector/nvector_serial.h>  /* Serial N_Vector types, functions, macros */
#include <sundials/sundials_dense.h> /* Definition of type DlsMat (dense matrix) */
#include <sundials/sundials_types.h> /* Definition of type realtype */

/*-----------------------------------------------------------------------*/

/* Functions called by the CVODE solver */

int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

int Jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J,
        void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*-----------------------------------------------------------------------*/

/* Type definition for user-supplied data passed to the solver functions */

typedef struct {
  realtype *rate, n_H, T_g, x_e;
} *User_Data;

/*-----------------------------------------------------------------------*/

/* Private function to print final abundances */

static void PrintFinalOutput(realtype *abundance, int npart, int nspec);

/* Private function to print final statistics */

static void PrintFinalStatistics(void *cvode_mem);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);

/*-----------------------------------------------------------------------*/

int calculate_abundances_(realtype *abundance, realtype *rate, realtype *density,
                          realtype *temperature, int *npart, int *nspec, int *nreac)
{
  int neq = *nspec-1; /* Number of ODEs in the system */
  realtype t0, t, tout; /* Initial, current and output times */
  realtype *data; /* Array pointer to access data stored in vectors*/
  N_Vector y; /* Vector of dependent variables that CVODE is to solve */
  User_Data user_data; /* Data to be passed to the solver routines */
  FILE *cvoderr; /* Log file for CVODE error messages */
  void *cvode_mem; /* Memory allocated to the CVODE solver */
  long int mxstep; /* Maximum number of internal steps */
  int flag, status, nfails; /* CVODE return flag, status and failure counter */

  /* Define the necessary external variables contained within the Fortran module CHEMISTRY_MODULE */
  extern realtype chemistry_module_mp_relative_abundance_tolerance_ ; /* Relative error tolerance */
  extern realtype chemistry_module_mp_absolute_abundance_tolerance_ ; /* Absolute error tolerance */

  /* Define the necessary external variables contained within the Fortran module GLOBAL_MODULE */
  extern realtype maincode_module_mp_start_time_; /* Start time for the chemical evolution (yr) */
  extern realtype maincode_module_mp_end_time_; /* End time for the chemical evolution (yr) */
  extern int global_module_mp_nelect_; /* Species index number for electron */

  realtype reltol = chemistry_module_mp_relative_abundance_tolerance_;
  realtype abstol = chemistry_module_mp_absolute_abundance_tolerance_;
  realtype start_time = maincode_module_mp_start_time_;
  realtype end_time = maincode_module_mp_end_time_;
  realtype seconds_in_year = RCONST(3.1556926e7); /* Convert from years to seconds */
  int nelect = global_module_mp_nelect_-1;

  double cpu_start, cpu_end; /* CPU times */
#ifdef OPENMP
  int nthread; /* Number of threads */
#endif

  int n, i;

  /*printf ("reltol: %E \n", reltol);
  printf ("abstol: %E \n", abstol);
  printf ("start_time: %E \n", start_time);
  printf ("end_time: %E \n", end_time);
  printf ("nelect:, %d \n", nelect);*/

  /* Open the error log files */
  /*stderr = fopen("main.log", "a");*/
  cvoderr = fopen("cvode.log", "a");

  /* Specify the maximum number of internal steps */
  mxstep = 10000000;

  /* Initialize the global status flag */
  status = 0;

#ifdef OPENMP
  cpu_start = omp_get_wtime();
#else
  cpu_start = clock();
#endif

  /*printf ("npart: %i \n",*npart);*/

#ifdef OPENMP
#pragma omp parallel for default(none) schedule(dynamic) \
  shared(npart, nspec, nreac, nelect, nthread, abundance, rate, density, temperature) \
  shared(neq, start_time, end_time, reltol, abstol, mxstep, status, cvoderr) \
  shared(seconds_in_year) \
  private(i, t0, tout, t, y, data, user_data, cvode_mem, flag, nfails)
#endif

  for (n = 0; n < *npart; n++) {

    /* Store the number of threads being used */
#ifdef OPENMP
    if (n == 0) nthread = omp_get_num_threads();
#endif
    /* Reset the integration failure counter */
    nfails = 0;

    /* Specify the start and end time of the integration (in seconds) */
    t0 = start_time*seconds_in_year;
    tout = end_time*seconds_in_year;


    /* Create a serial vector of length NEQ to contain the abundances */
    y = NULL;
    y = N_VNew_Serial(neq);
    if (check_flag((void *)y, "N_VNew_Serial", 0)) status = 1;

    /* Initialize y from the abundance array */
    data = NV_DATA_S(y);
    for (i = 0; i < neq; i++) {
#ifdef LOG_ODES
      if (abundance[n*(*nspec)+i] > 0) data[i] = log(abundance[n*(*nspec)+i]);
      else data[i] = log(abstol*abstol);
#else
      data[i] = abundance[n*(*nspec)+i];
#endif
    }
    /*for (i=212; i<215; i++) {
     printf ("n %i ,abundance( %i )= %E \n",n,i,abundance[n*(*nspec)+i]); 
    }
    printf ("\n");*/

    /* Create and allocate memory to user_data to contain the rates */
    user_data = NULL;
    user_data = (User_Data) malloc(sizeof *user_data);
    if(check_flag((void *)user_data, "malloc", 2)) status = 1;

    /* Initialize user_data with the array of reaction rate coefficients,
     * the total number density, gas temperature and electron abundance */
    user_data->rate = &rate[n*(*nreac)];
    user_data->n_H = density[n];
    user_data->T_g = temperature[n];
    user_data->x_e = abundance[n*(*nspec)+nelect];

    /* Call CVodeCreate to create the solver memory and specify the
     * use of Backward Differentiation Formula and Newton iteration */
    cvode_mem = NULL;
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) status = 1;

    /* Call CVodeSetErrFile to direct all error messages to the log file */
    flag = CVodeSetErrFile(cvode_mem, cvoderr);
    if (check_flag(&flag, "CVodeSetErrFile", 1)) status = 1;

    /* Call CVodeInit to initialize the integrator memory and specify the
     * right hand side function in ydot = f(t,y), the inital time t0, and
     * the initial dependent variable vector y */
    flag = CVodeInit(cvode_mem, f, t0, y);
    if (check_flag(&flag, "CVodeInit", 1)) status = 1;

    /* Call CVodeSStolerances to set the relative and absolute tolerances */
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag, "CVodeSStolerances", 1)) status = 1;

    /* Call CVodeSetMaxNumSteps to set the maximum number of steps */
    flag = CVodeSetMaxNumSteps(cvode_mem, mxstep);
    if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) status = 1;

    /* Specify the user-defined data to be passed to the various routines */
    flag = CVodeSetUserData(cvode_mem, user_data);
    if (check_flag(&flag, "CVodeSetUserData", 1)) status = 1;

    /* Specify that the CVDense direct dense linear solver is to be used */
    flag = CVDense(cvode_mem, neq);
    if (check_flag(&flag, "CVDense", 1)) status = 1;

    /* Specify that a user-supplied Jacobian routine (Jac) is to be used */
    /* flag = CVDlsSetDenseJacFn(cvode_mem, Jac); */
    /* if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) status = 1; */

    tout = 1.0e-4*seconds_in_year;

    do { /* Call CVode, check the return status and loop until the end time is reached */

      flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
#ifdef OPENMP
#pragma omp critical (status)
#endif
      if (flag == 0) {
        if (nfails > 0) {
          fprintf(cvoderr, "Call to CVODE successful: particle = %d, t = %8.2le yr, # steps = %ld\n\n", n, t/seconds_in_year, mxstep);
        }
      }
      else {
        nfails++;
        if (flag == -1) {
          if (nfails < 5) {
            fprintf(cvoderr, "  Particle %d: Doubling mxstep and continuing the integration.\n\n", n);
            flag = CVodeSetMaxNumSteps(cvode_mem, mxstep);
            if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) status = 1;
	  }
          else {
            /*fprintf(stderr, " WARNING! CVODE solver failed: flag = %d\n", flag);*/
            fprintf(cvoderr, "ERROR! Fifth failure - aborting integration\n\n");
            fprintf(cvoderr, "-----------------------------------------------\n");
            fprintf(cvoderr, "CVODE Parameters:\n\n");
            fprintf(cvoderr, "Particle %d\n\n", n);
            fprintf(cvoderr, "Before:\n");
            fprintf(cvoderr, "flag   = %d\n", 0);
            fprintf(cvoderr, "t      = %8.2le yr\n", start_time);
            fprintf(cvoderr, "tout   = %8.2le yr\n", end_time);
            fprintf(cvoderr, "dt     = %8.2le yr\n\n", end_time-start_time);
            fprintf(cvoderr, "After:\n");
            fprintf(cvoderr, "flag   = %d\n", flag);
            fprintf(cvoderr, "t      = %8.2le yr\n", t/seconds_in_year);
            fprintf(cvoderr, "tout   = %8.2le yr\n", tout/seconds_in_year);
            fprintf(cvoderr, "dt     = %8.2le yr\n", (tout-t)/seconds_in_year);
            fprintf(cvoderr, "-----------------------------------------------\n\n");

            status = flag;
          }
        }
        else if (flag == -4) {
          if (nfails < 5) {
            fprintf(cvoderr, "  Particle %d: Attempting to continue the integration using the same parameters.\n\n", n);
          }
          else {
            /*fprintf(stderr, " WARNING! CVODE solver failed: flag = %d\n", flag);*/
            fprintf(cvoderr, "ERROR! Fifth failure - aborting integration\n\n");
            fprintf(cvoderr, "-----------------------------------------------\n");
            fprintf(cvoderr, "CVODE Parameters:\n\n");
            fprintf(cvoderr, "Particle %d\n\n", n);
            fprintf(cvoderr, "Before:\n");
            fprintf(cvoderr, "flag   = %d\n", 0);
            fprintf(cvoderr, "t      = %8.2le yr\n", start_time);
            fprintf(cvoderr, "tout   = %8.2le yr\n", end_time);
            fprintf(cvoderr, "dt     = %8.2le yr\n\n", end_time-start_time);
            fprintf(cvoderr, "After:\n");
            fprintf(cvoderr, "flag   = %d\n", flag);
            fprintf(cvoderr, "t      = %8.2le yr\n", t/seconds_in_year);
            fprintf(cvoderr, "tout   = %8.2le yr\n", tout/seconds_in_year);
            fprintf(cvoderr, "dt     = %8.2le yr\n", (tout-t)/seconds_in_year);
            fprintf(cvoderr, "-----------------------------------------------\n\n");
            status = flag;
          }
        }
        else {
          /*fprintf(stderr, " WARNING! CVODE solver failed: flag = %d\n", flag);*/
          fprintf(cvoderr, "ERROR! CVODE failed with the following parameters:\n\n");
          fprintf(cvoderr, "-----------------------------------------------\n");
          fprintf(cvoderr, "CVODE Parameters:\n\n");
          fprintf(cvoderr, "Particle %d\n\n", n);
          fprintf(cvoderr, "Before:\n");
          fprintf(cvoderr, "flag   = %d\n", 0);
          fprintf(cvoderr, "t      = %8.2le yr\n", start_time);
          fprintf(cvoderr, "tout   = %8.2le yr\n", end_time);
          fprintf(cvoderr, "dt     = %8.2le yr\n\n", end_time-start_time);
          fprintf(cvoderr, "After:\n");
          fprintf(cvoderr, "flag   = %d\n", flag);
          fprintf(cvoderr, "t      = %8.2le yr\n", t/seconds_in_year);
          fprintf(cvoderr, "tout   = %8.2le yr\n", tout/seconds_in_year);
          fprintf(cvoderr, "dt     = %8.2le yr\n", (tout-t)/seconds_in_year);
          fprintf(cvoderr, "-----------------------------------------------\n\n");
          status = flag;
          nfails = 5;
        }
      }

      tout = tout*10;
      if (tout > end_time*seconds_in_year) tout = end_time*seconds_in_year;

    } while (t < tout && nfails < 5);

    /* Store the output values in the abundance array */
    data = NV_DATA_S(y);
    for (i = 0; i < neq; i++) {
#ifdef LOG_ODES
      abundance[n*(*nspec)+i] = exp(data[i]);
#else
      abundance[n*(*nspec)+i] = data[i];
#endif
      if(abundance[n*(*nspec)+i] < abstol) abundance[n*(*nspec)+i] = 0;
    }
    abundance[n*(*nspec)+nelect] = user_data->x_e;
   
    /* Free y vector memory allocation */
    N_VDestroy_Serial(y);

    /* Free user_data memory allocation */
    free(user_data);

    /* Free integrator memory allocation */
    CVodeFree(&cvode_mem);

  } /* End of for-loop over particles */

#ifdef OPENMP
  cpu_end = omp_get_wtime();
#else
  cpu_end = clock();
#endif

  /* Close the error log files */
  fclose(cvoderr);
  /*fclose(stderr);*/

  if (status != 0) {
    fprintf(stdout, "\n ERROR! Calculation of abundances failed with status = %d\n\n", status);
#ifdef OPENMP
    fprintf(stdout, " Threads used = %d (maximum %d)\n\n", nthread, omp_get_max_threads());
#endif
  }

  else {
#ifdef OPENMP
    fprintf(stdout, "  --> Threads used = %d (maximum %d)\n", nthread, omp_get_max_threads());
#endif
  }

  return(status);
}
/*=======================================================================*/

/*
 *--------------------------
 * Private helper functions
 *--------------------------
 */

/*=======================================================================

 Print the final abundances returned from CVODE for each particle

-----------------------------------------------------------------------*/
static void PrintFinalOutput(realtype *abundance, int npart, int nspec)
{
  int n, i;

  printf("\n");
  for (n = 0; n < npart; n++) {
    /* printf("Final y(%d) =", n+1); */
    for (i = 0; i < nspec; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      printf("Final y(%d,%2d) =%12.5Le\n", n+1, i+1, abundance[n*(nspec)+i]);
      /* printf("%12.5Le", abundance[n*(nspec)+i]); */
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      printf("Final y(%d,%2d) =%12.5le\n", n+1, i+1, abundance[n*(nspec)+i]);
      /* printf("%12.5le", abundance[n*(nspec)+i]); */
#else
      printf("Final y(%d,%2d) =%12.5e\n", n+1, i+1, abundance[n*(nspec)+i]);
      /* printf("%12.5e", abundance[n*(nspec)+i]); */
#endif
    }
    printf("\n");
  }
  return;
}
/*=======================================================================*/

/*=======================================================================

 Print the final statistics returned from the most recent call to CVODE

-----------------------------------------------------------------------*/
static void PrintFinalStatistics(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

#ifdef OPENMP
#pragma omp critical (statistics)
    printf("\nFinal Statistics (thread %d):\n", omp_get_thread_num());
#endif
    printf("nst = %-6ld nfe = %-6ld nsetups = %-6ld nfeLS = %-6ld\n",
	   nst, nfe, nsetups, nfeLS);
    printf("nje = %-6ld nni = %-6ld ncfn = %-6ld    netf = %-6ld\n",
	   nje, nni, ncfn, netf);
}
/*=======================================================================*/

/*=======================================================================

 Check SUNDIALS/CVODE function return values

    opt == 0 means the SUNDIALS function allocates memory,
             so check if a NULL pointer has been returned

    opt == 1 means the SUNDIALS function returns a flag,
             so check if the flag < 0

    opt == 2 means the called function allocates memory,
             so check if a NULL pointer has been returned

-----------------------------------------------------------------------*/
static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if the SUNDIALS function returned a NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    /*fprintf(stderr, "\n SUNDIALS ERROR! %s() failed - returned NULL pointer\n\n",
	    funcname);*/
    return(1); }

  /* Check if the function returned a flag < 0 - an error occurred */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      /*fprintf(stderr, "\n SUNDIALS ERROR! %s() failed with flag = %d\n\n",
	      funcname, *errflag);*/
      return(1); }}

  /* Check if the function returned a NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    /*fprintf(stderr, "\n MEMORY ERROR! %s() failed - returned NULL pointer\n\n",
	    funcname);*/
    return(1); }

  return(0);
}
/*=======================================================================*/
