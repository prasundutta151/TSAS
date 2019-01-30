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

double *data;

# define POWSP(kval, alpha) (pow(kval, alpha)/2.)

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

int Gen_Gauss_2D(int NGrid, float LL, float Amp, float alpha, long seed, float mean, float BMaj, float BMin){
  
  int ii, jj, index;
  double *out;
  double kval, kmax, kmin, AA, BB, Norm;
  double sum, power;
  fftw_complex *fmode;
  fftw_plan c2r_2D;
  char filename[100];
  FILE *fp;

  out = (double *)calloc(NGrid*(NGrid+2), sizeof(double));

  // norm = 1./(1.*NGrid);
  Norm = 1./(sqrt(1.*NGrid*NGrid)*BMaj*BMin);
  fmode = (fftw_complex *)&out[0];
  c2r_2D = fftw_plan_dft_c2r_2d(NGrid, NGrid, fmode, out, FFTW_ESTIMATE);
  
  kmin = 1./LL/1.e3/2.;
  printf("Norm = %e\n", Norm);
  printf("kmin = %.3e\n", kmin);

  sum = 0.;
  // surface
  
  for(ii=0; ii<NGrid; ii++)
    for(jj=1; jj<NGrid/2; jj++){


	kval = kmin*sqrt(MSQR(1.*((ii>NGrid/2)?(NGrid-ii):ii)) 
			 + MSQR(1.*jj));
	index = ii*(NGrid/2+1) + jj;
	fmode[index][0] = sqrt(POWSP(kval, alpha))*gasdev(&seed);
	fmode[index][1] = sqrt(POWSP(kval, alpha))*gasdev(&seed);
	power = MSQR(1.*fmode[index][0]) + MSQR(1.*fmode[index][1]);
	sum +=  power;
    }
  jj=0;
  for(ii=0; ii<NGrid/2; ii++){

    kval = kmin*sqrt(MSQR(1.*ii) + MSQR(1.*jj));
    AA=sqrt(POWSP(kval, alpha))*gasdev(&seed); 
    BB=sqrt(POWSP(kval, alpha))*gasdev(&seed); 
    
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
      data[index] = sqrt(Amp)*(data[index] - Mean);
     }

  printf("\n\n2D MAP GENERATED WITH PWER SPECTRUM P(K) =  K ^ %.1f \n", alpha);
  printf("Mean = %e\t Rms = %e\n\n", Mean, Rms);
  return 0;
  
}

int main(int argc, char *argv[]){

  int ii, jj, index, NBin;
  float *xx, *yy, LL, DLRAD, Amp, Alp, mean, kmax, kmin;
  double DL, BMaj, BMin;
  long seed, NGrid;
  FILE *fp;
  char INPFILE[128], OUTFITS[128];
  char param[128];

  if(argc!=4){
    printf("USAGE : %s <input file>  <outfile> <seed>\n", argv[0]);
    return 1;
  }
  sscanf(argv[1], "%s", INPFILE);
  sscanf(argv[2], "%s", OUTFITS);
  seed = (long)atoi(argv[3]);

  if(Access(INPFILE)==0){
    printf("Input File %s does not exists\n", INPFILE);
    return EXIT_FAILURE;
  }
  if(Access(OUTFITS)!=0){
    printf("Output File %s exists\n", OUTFITS);
    return EXIT_FAILURE;
  }
  
  fp = fopen(INPFILE, "r");

  fgets(param,100, fp);
  fgets(param,100, fp);
  fgets(param,100, fp);
  
  fgets(param,100, fp);
  sscanf(param,"%ld%lf", &NGrid, &DL);
  
  fgets(param,100, fp);
  sscanf(param,"%f%f", &Amp, &Alp, &mean);

  fgets(param,100, fp);
  sscanf(param,"%lf%lf", &BMaj, &BMin);

  fclose(fp);

  DL = DL/3600.;

  seed = seed%SMAX;
  seed = (seed>0) ? - seed : seed;

  fitsfile *fptro;
  int status=0, anynul=0;
  float nulval=0.;
  long nel, naxis, *naxes;
  char comment[FLEN_COMMENT];

   printf("NGrid = %ld\n", NGrid);
  printf("DL = %f (sec)\n", fabs(DL*3600.));
  DLRAD = fabs(DL*M_PI/180.);
  LL = (1.*NGrid*DLRAD);
  nel = NGrid*NGrid;
  data  = (double *)calloc(nel, sizeof(double));
  BMaj *= (M_PI/180.);
  BMin *= (M_PI/180.);

  Gen_Gauss_2D((int)NGrid, LL, Amp, Alp, seed, mean, BMaj, BMin);
 
  status = 0;
  naxis  = 2;
  naxes = (long *)calloc(naxis, sizeof(long));
  naxes[0] = NGrid;
  naxes[1] = NGrid;
 
  BMaj = BMaj*(180./M_PI);
  BMin = BMin*(180./M_PI);

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
