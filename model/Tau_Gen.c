/*  NAME:         Tau_Gen.c                                       */
/*  Derived from: Gauss_2D_FITS_sal.c                             */
/*  Description:     Generates a FITS file with Gaussian  random  */
/*                   numbers with a powerspectrum and geometry    */
/*                   from an ASCII  file. Data mean is non-zero.  */
/*                   CDELT1/2 and BMAJ/MIN are written in degree  */
/*                   A time based seed build in the program.      */
/*  Developer:    Prasun Dutta, IIT (BHU), Varanasi               */
/*  Contact:      pdutta.phy@itbhu.ac.in                          */
/*  InputFIle:    Example-----------------------------------------*/
/*
# ==================== GAUSS GEN PARAMEERS =================================
# VALUE --------------------  PARAM  --------TYP---- COMMENT ---------------
#
1024             1.E-3      : NGrid, dxy    [I, D]    [,deg]
1.6E-2		 1.5E-2	    : BMaj, BMin    [D, D]    [deg, deg]
1.                          : Tau0          [D]       
*/
/* =======================================================================*/
/* ISSUES:   Normalization along with PowFits2D_dbl                       */

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <nr.h>
# include <nrutil.h>
# include <unistd.h>
# include <string.h>
# include <fitsio.h>
# include <fftw3.h>
# include <sys/stat.h>

# define EXP_NARG 7
# define SMAX  (8388608)
# define MSQR(a) (a*a)
# define DTR (M_PI/180.)

double *data;
long NGrid;
float alpha, sigma, Tau0;

float POWSP(float kval){

  float Val;
  Val = pow(kval, alpha);
  return Val;
}

int Access(char *filename){
  struct stat   buffer;   
  return (stat (filename, &buffer) == 0);
}

void printerror(int status){
  if (status){
    fits_report_error(stderr, status); // Print error report
    exit( status );    // Terminate the program, returning error status
  }
}

int Gen_Gauss_2D(long seed){
  
  int ii, jj, index;
  double *out;
  double kval, kmax, kmin, AA, BB, Norm, BeamArea, BeamPerPixel;
  double sum, power;
  fftw_complex *fmode;
  fftw_plan c2r_2D;
  char filename[100];
  FILE *fp;

  out = (double *)calloc(NGrid*(NGrid+2), sizeof(double));

  Norm = 1./sqrt(1.*NGrid*NGrid);

  fmode = (fftw_complex *)&out[0];
  c2r_2D = fftw_plan_dft_c2r_2d(NGrid, NGrid, fmode, out, FFTW_ESTIMATE);
  
  kmin = 1.;
  printf("Norm = %e\n", Norm);
  printf("kmin = %.3e\n", kmin);

  sum = 0.;
  // surface
  
  for(ii=0; ii<NGrid; ii++)
    for(jj=1; jj<NGrid/2; jj++){


	kval = kmin*sqrt(MSQR(1.*((ii>NGrid/2)?(NGrid-ii):ii)) 
			 + MSQR(1.*jj));
	index = ii*(NGrid/2+1) + jj;
	fmode[index][0] = sqrt(POWSP(kval))*gasdev(&seed);
	fmode[index][1] = sqrt(POWSP(kval))*gasdev(&seed);
	power = MSQR(1.*fmode[index][0]) + MSQR(1.*fmode[index][1]);
	sum +=  power;
	//printf("Kval1 = %e Power = %e\n", kval, power);
    }
  jj=0;
  for(ii=0; ii<NGrid/2; ii++){

    kval = kmin*sqrt(MSQR(1.*ii) + MSQR(1.*jj));

    AA=sqrt(POWSP(kval))*gasdev(&seed); 
    BB=sqrt(POWSP(kval))*gasdev(&seed); 
    
    index = ii*(NGrid/2 + 1) + jj;
    fmode[index][0] = AA; 
    fmode[index][1] = BB;
    if(ii!=0){
      power = MSQR(1.*fmode[index][0]) + MSQR(1.*fmode[index][1]);
      sum +=  power;
    }  
    index = (NGrid-ii)*(NGrid/2 + 1)+jj; 
    fmode[index][0] = AA;
    fmode[index][1] = -1.*BB;
    if(ii!=0){
      power = MSQR(1.*fmode[index][0]) + MSQR(1.*fmode[index][1]);
      sum +=  power;
    }  
    //   printf("Kval1 = %e Power = %e\n", kval, power);
  }
  fmode[0][0] = 0.; 
  fmode[0][1] = 0.;
  
  printf("total power = %.3e\n", sum);
  // finished filling Fourier coefficients

  fftw_execute(c2r_2D);
 
  double Mean, Rms;
  long NDATA;
  
  NDATA = NGrid*NGrid;
  Mean = Rms = 0.;
  
  for(ii=0; ii<NGrid; ii++)
    for(jj=0; jj<NGrid; jj++){

      index = ii*NGrid + jj;
      data[index] = (out[ii*(NGrid+2) + jj])*Norm;
      Mean += data[index];
      Rms += MSQR(data[index]);
    }
  
  Mean /= (1.*NDATA);
  Rms = sqrt(Rms/(1.*NDATA) - MSQR(Mean));
  int counter;
  counter = 0;
  for(ii=0; ii<NGrid; ii++)
    for(jj=0; jj<NGrid; jj++){

      index = ii*NGrid + jj;
      data[index] = sigma*(data[index] - Mean)/Rms + Tau0;
      if(data[index]<0.){
	data[index] = 0.;
	counter++;
      }
     }
  printf("Tau value  zero reset = %d\n", counter);
  printf("percentage zero reset = %.2f\n", (1.*counter)/(1.*NGrid*NGrid));
  return 0;
}

int main(int argc, char *argv[]){

  int ii, jj, index, NBin;
  long seed;
  FILE *fp;
  char OUTFITS[128];

  if(argc!=EXP_NARG){
    printf("\nUSAGE : %s <OUTFITS> <NGrid> <Tau0> <sigma> <alpha> <seed>\n", argv[0]);
    printf("OUTFITS: Output FITS file\n");
    printf("NGrid  : No of Grid points in image (x and y same)\n");
    printf("Tau0   : Mean Tau\n");
    printf("sigma  : standard dev of the distribution\n");
    printf("alpha  : slope of the distribution\n");
    printf("seed   : initializing random seed\n\n");
    printf("\n");
    return 1;
  }
 
  sscanf(argv[1], "%s", OUTFITS);
  NGrid = (long)atoi(argv[2]);
  Tau0  = (double)atof(argv[3]);
  sigma = (double)atof(argv[4]);
  alpha = (double)atof(argv[5]);
  seed  = (long)atoi(argv[6]);

  printf("Inputs are read\n");
  if(Access(OUTFITS)!=0){
    printf("Output File %s exists\n", OUTFITS);
    return EXIT_FAILURE;
  }
  
  alpha = (alpha>0.) ? -alpha : alpha;
  sigma = (sigma<0.) ? -sigma : sigma;
  Tau0  = (Tau0<0.)  ? -Tau0  : Tau0;

  seed = seed%SMAX;
  seed = (seed>0) ?  -seed : seed;
  /* srand(seed); */

  /* seed = seed%SMAX; */
  /* seed = (seed>0) ? - seed : seed; */
  printf("Using seed: %d\n", seed);

  fitsfile *fptro;
  int status=0, anynul=0;
  float nulval=0.;
  double DL;
  long nel, naxis, *naxes;
  char comment[FLEN_COMMENT];

  printf("NGrid = %ld\n", NGrid);
  DL = 1.;

  nel = NGrid*NGrid;
  data  = (double *)calloc(nel, sizeof(double));

  Gen_Gauss_2D(seed);
 
  status = 0;
  naxis  = 2;
  naxes = (long *)calloc(naxis, sizeof(long));
  naxes[0] = NGrid;
  naxes[1] = NGrid;

  sprintf(comment, "\\");
  if(fits_create_file(&fptro, OUTFITS, &status))
    printerror(status);
  if(fits_create_img(fptro,  DOUBLE_IMG, naxis, naxes, &status))
    printerror(status);
  if(fits_write_key(fptro, TDOUBLE, "CDELT1", &DL, comment, &status))
    printerror( status );
  if(fits_write_key(fptro, TDOUBLE, "CDELT2", &DL, comment, &status))
    printerror( status );

  if(fits_write_img(fptro, TDOUBLE, 1, nel,  data,  &status))
    printerror(status);
  printf("Data stored in datao\n");
  fits_close_file(fptro,&status);

  printf("OUTFILE = %s\n\n", OUTFITS);
  return 0;
}
