/*  NAME:         Gauss_2D_FITS_ARB.c                             */
/*  Derived from: Gauss_2D_FITS_sal.c                             */
/*  Description:     Generates a FITS file with Gaussian  random  */
/*                   numbers with a given powerspectrum specified */
/*                   in an ASCII  file. It Interpolates the values*/
/*                   at the grid point. Data mean is zero.        */
/*                   CDELT1/2 and BMAJ/MIN are written in degree  */
/*  Developer:    Prasun Dutta, IIT (BHU), Varanasi               */
/*  Contact:      pdutta.phy@itbhu.ac.in                          */
/*  InputFIle:    Example-----------------------------------------*/
/*
# ==================== GAUSS GEN PARAMEERS =================================
# VALUE --------------------  PARAM  --------TYP---- COMMENT ---------------
#
1.6E-3		 1.5E-3	    : BMaj, BMin    [D, D]    [deg, deg]
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

# define SMAX  (8388608)
# define MSQR(a) (a*a)
# define DTR (M_PI/180.)
int NPS;
double *data;
float  *U0, *P0, *P2;
float delxy;

float POWSP(float kval){

  float Val;
  splint(U0+1, P0+1, P2+1, NPS, kval, &Val);
  return Val/2.;
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

int Gen_Gauss_2D(int NGrid, float LL, long seed, float mean, float BMaj, float BMin){
  
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
  BeamArea = M_PI*BMaj*BMin/(4.*log(2.));
  BeamPerPixel = BeamArea/delxy/delxy;
  Norm =  Norm/BeamPerPixel;
  printf("delxy = %.3e\tBMaj = %.3e\tBMin = %.3e\n", delxy, BMaj, BMin);

  fmode = (fftw_complex *)&out[0];
  c2r_2D = fftw_plan_dft_c2r_2d(NGrid, NGrid, fmode, out, FFTW_ESTIMATE);
  
  kmin = 1./LL/1.e3;
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
	//    	printf("Kval1 = %e Power = %e\n", kval, power);
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
  
  for(ii=0; ii<NGrid; ii++)
    for(jj=0; jj<NGrid; jj++){

      index = ii*NGrid + jj;
      data[index] = (data[index] - Mean);
     }

  return 0;
  
}

int main(int argc, char *argv[]){

  int ii, jj, index, NBin;
  float *xx, *yy, LL, mean, kmax, kmin, tempf;
  float slope_beg, slope_end;
  double DL, BMaj, BMin;
  long seed, NGrid;
  FILE *fp;
  char CONTFILE[128], OUTFITS[128];
  char param[128];

  if(argc!=6){
    printf("\nUSAGE : %s <CONTFILE>  <OUTFITS> <seed> <sigma> <alpha>\n", argv[0]);
    printf("\nFITSFILE: Continium FITS file\n");
    printf("OUTFITS: Output FITS file\n");
    printf("seed   : seed for the random number generator\n");
    printf("sigma  : standard dev of the distribution\n");
    printf("alpha  : slope of the distribution\n");
    printf("\n");
    return 1;
  }
  sscanf(argv[1], "%s", CONTFILE);
  sscanf(argv[2], "%s", OUTFITS);
  seed = (long)atoi(argv[3]);
  sigma = (double)atof(argv[4]);
  alpha = (double)atof(argv[5]);

  if(Access(CONTFILE)==0){
    printf("Input File %s does not exists\n", CONTFILE);
    return EXIT_FAILURE;
  }
  if(Access(OUTFITS)!=0){
    printf("Output File %s exists\n", OUTFITS);
    return EXIT_FAILURE;
  }
  
  seed = seed%SMAX;
  seed = (seed>0) ? - seed : seed;

  NPS = 0;
  fp = fopen(PSPFILE, "r");
  while(fscanf(fp, "%f%*f", &tempf)!=EOF)
    NPS++;
  fclose(fp);
  printf("NPS = %d\n", NPS);

  U0 = (float *)calloc(NPS+1, sizeof(float));
  P0 = (float *)calloc(NPS+1, sizeof(float));
  P2 = (float *)calloc(NPS+1, sizeof(float));

  fp = fopen(PSPFILE, "r");
  for(ii=0; ii<NPS; ii++)
    fscanf(fp, "%e%e", &U0[ii+1], &P0[ii+1]);
  fclose(fp);
 
  slope_beg = (U0[1]-U0[2])/(P0[1]-P0[2]);
  slope_end = (U0[NPS-1]-U0[NPS])/(P0[NPS-1]-P0[NPS]);

  slope_beg = slope_end = 0.;
  printf("Slope: %e\t%e\n", slope_beg, slope_end);
  spline(U0+1, P0+1, NPS, slope_beg, slope_end, P2);

  /* Testing SPLINE_SPLINT  */
  /*
  fp = fopen("test.pow", "w");
  float UU, PU;
  for(ii=0; ii<NPS*10; ii++){

    UU = U0[1]*1.1 + (U0[NPS]*0.9-U0[1]*1.1)*(1.*ii)/(10.*NPS);
    splint(U0+1, P0+1, P2+1, NPS, UU, &PU);
    fprintf(fp, "%e\t%e\t%e\n", UU, PU, POWSP(UU)*2.);
  }
  fclose(fp);
  */
  /* Test Done */

  printf("Umin = %e\n", U0[1]);
  printf("Umax = %e\n", U0[NPS]);

  LL    = 1./(U0[1]*1.01e3);
  delxy = 1./(U0[NPS]*0.9)/2./1.e3;

  NGrid = (int)floor(LL/delxy);
  NGrid = NGrid + NGrid%2;
  
  delxy = LL/(1.*NGrid);
  delxy = delxy/DTR;
  DL    = delxy;

  fitsfile *fptro;
  int status=0, anynul=0;
  float nulval=0.;
  long nel, naxis, *naxes;
  char comment[FLEN_COMMENT];

  printf("NGrid = %ld\n", NGrid);
  printf("DL = %f (sec)\n", fabs(delxy*3600.));

  nel = NGrid*NGrid;
  data  = (double *)calloc(nel, sizeof(double));

  Gen_Gauss_2D((int)NGrid, LL, seed, 0., BMaj, BMin);
 
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
 if(fits_write_key(fptro, TDOUBLE, "BMAJ", &BMaj, comment, &status))
    printerror( status );
  if(fits_write_key(fptro, TDOUBLE, "BMIN", &BMin, comment, &status))
    printerror( status );

  if(fits_write_img(fptro, TDOUBLE, 1, nel,  data,  &status))
    printerror(status);
  printf("Data stored in datao\n");
  fits_close_file(fptro,&status);

  printf("OUTFILE = %s\n\n", OUTFITS);
  return 0;
}
