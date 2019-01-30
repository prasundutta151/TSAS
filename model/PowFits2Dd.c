
/* ISSUES:   Normalization along with Gauss_2D_FITS_ARB                    */

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <unistd.h>
# include <nrutil.h>
# include <nr.h>
# include <unistd.h>
# include <fftw3.h>
# include <fitsio.h>
# define MYSQR(a) (a*a)
# define DTR (M_PI/180.)

 double BMaj, BMin;

int READ_FITS(char *INFILE, int npixels, double *Data){

  fitsfile *fptr;  
  int  status,  anynull,  nullval;
  long fpixel;
  double *DATAI;
  char comment[FLEN_KEYWORD], command[150];
 
  status = 0;
  if(fits_open_file(&fptr, INFILE, READONLY, &status))
    return status;
  
  fpixel   = 1;
  nullval  = 0;
  
  if(fits_read_img(fptr, TDOUBLE, fpixel, npixels, &nullval, Data, &anynull, &status))
    return status;
 
  if(fits_close_file(fptr, &status))
    return status;

  return EXIT_SUCCESS;
}


void Get_Vis_2D(double *DATA, int NGrid, double delxy){
  
  int ii, jj, index;
  fftw_plan r2c_2DI;
  fftw_complex *Vis;
  double Norm, BeamArea, BeamPerPixel;

  Vis = (fftw_complex *)&DATA[0];

  r2c_2DI = fftw_plan_dft_r2c_2d(NGrid, NGrid, DATA, Vis, FFTW_ESTIMATE);
  fftw_execute(r2c_2DI);
  fftw_destroy_plan(r2c_2DI);

  Norm = 1./sqrt(1.*NGrid*NGrid);
  BeamArea = M_PI*BMaj*BMin/(4.*log(2.));
  BeamPerPixel = BeamArea/delxy/delxy;
  Norm =  Norm*BeamPerPixel;
  printf("delxy = %.3e\tBMaj = %.3e\tBMin = %.3e\n", delxy, BMaj, BMin);

  printf("Norm = %e\n", Norm);
  for(ii=0; ii<NGrid; ii++)
    for(jj=0; jj<NGrid/2+1; jj++){
      index =  ii*(NGrid/2+1) + jj;
      Vis[index][0] = Vis[index][0]*Norm;
      Vis[index][1] = Vis[index][1]*Norm;
    }
  DATA = (double *)&Vis[0];
}
int CALPOW_2D(double *Data, int NGrid, int NBin, double delxy, char *outfile){

  int ii, jj, index, index2, pind;
  double Kval, RE, IM, buffval;
  double *DATA;
  double *kmode, *power, *sqpow, Kmax, Kmin, logalp;
  fitsfile *fptr;  
  int  status,  nfound;
  long npixels, *naxes, naxis;
  char comment[FLEN_KEYWORD];
  long *no;
  FILE *fp;
  
  printf("CALPOW_2D\n");
 
  if((DATA  = (double *)calloc(NGrid*(NGrid+2), sizeof(double)))==NULL){
    printf("An attempt to allocate %d bytes failed.\n", NGrid*(NGrid+2));
    return EXIT_FAILURE;
  }
  for(ii=0; ii<NGrid; ii++)
    for(jj=0; jj<NGrid; jj++){

      index = ii*NGrid + jj;
      index2 = ii*(NGrid+2) + jj;
      DATA[index2] = Data[index];
    }

  Get_Vis_2D(DATA, NGrid, delxy);

  if((kmode  = (double *)calloc(NBin+1, sizeof(double)))==NULL){
    printf("An attempt to allocate %d bytes failed.\n", NBin);
    return 1;
  }
  if((power  = (double *)calloc(NBin+1, sizeof(double)))==NULL){
    printf("An attempt to allocate %d bytes failed.\n", NBin);
    return 1;
  }
  if((sqpow  = (double *)calloc(NBin+1, sizeof(double)))==NULL){
    printf("An attempt to allocate %d bytes failed.\n", NBin);
    return 1;
  }
  if((no     = (long *)calloc(NBin+1, sizeof(long)))==NULL){
    printf("An attempt to allocate %d bytes failed.\n", NBin);
    return 1;
  }

  buffval = 0.1;
  Kmin = 1.e-3/(1.*NGrid*delxy)/DTR;
  Kmax = 1.e-3/(2.*delxy)/DTR;

  printf("delxy = %e\n", delxy);
  printf("Umin = %f\n", Kmin);
  printf("Kmax = %f\n", Kmax);

  logalp = NBin/log10(Kmax/Kmin);

  for(ii=0;ii<NGrid;++ii)
    for(jj=0;jj<NGrid/2;++jj){

      Kval = Kmin*sqrt(MYSQR(1.*((ii>NGrid/2)? NGrid-ii: ii)) +
		       MYSQR(1.*jj));
      index = ii*(NGrid/2 +1) + jj;      
      if((ii+jj>1)&&(Kval<Kmax)){
	
	pind = (int)floor(logalp*log10(Kval/Kmin));
	RE = DATA[index*2];
	IM = DATA[index*2+1];
	power[pind] += (RE*RE + IM*IM);
	sqpow[pind] += (RE*RE + IM*IM)*(RE*RE + IM*IM);
	kmode[pind] += Kval;
	no[pind]++; 
      }
    }
  
  fp = fopen(outfile, "w");
  for(ii=0; ii<NBin; ii++){ 
  
    if(no[ii] > 0){
     
      power[ii] /= (1.*no[ii]);
      kmode[ii] /= (1.*no[ii]);
      sqpow[ii] = sqpow[ii]/(1.*no[ii]) - MYSQR(power[ii]);
      sqpow[ii] = sqrt( sqpow[ii]/(8.*M_PI) + MYSQR(power[ii])/(1.*no[ii]/(4.*MYSQR(M_PI))));
      if(kmode[ii]<Kmax){
	fprintf(fp, "%e\t%e\t%e\t%ld\n", kmode[ii], power[ii], sqpow[ii], no[ii]);
      }
    }
  }
  fclose(fp);

  printf("PS written in file %s\n", outfile);
 
  return 0;
}


int main(int argc, char *argv[]){

  int NBin, NGrid;
  double *Data, delxy;
  fitsfile *fptr;  
  int  status,  nfound;
  long npixels, *naxes, naxis;
  char comment[FLEN_KEYWORD];

  if(argc!=4){
    printf("USAGE : %s <Infile>  <Output file> <NBin>\n", argv[0]);
    return EXIT_FAILURE;
  }
  NBin = atoi(argv[3]);
 
  if(access(argv[1], F_OK)!=0){
    printf("File %s does not exists\n", argv[1]);
    EXIT_FAILURE;
  }
  
  status = 0;
  if(fits_open_file(&fptr, argv[1], READONLY, &status))
    return status;
   if(fits_read_key_lng(fptr, "NAXIS", &naxis, comment, &status))
    return status;
  naxes = (long *)calloc(naxis, sizeof(long));
  if(fits_read_keys_lng(fptr, "NAXIS", 1, naxis, naxes, &nfound, &status))
    return status;
  if(fits_read_key_dbl(fptr, "CDELT1", &delxy, comment, &status))
    return status;
  if(fits_read_key_dbl(fptr, "BMAJ", &BMaj, comment, &status))
    return status;
  if(fits_read_key_dbl(fptr, "BMIN", &BMin, comment, &status))
    return status;

  if(fits_close_file(fptr, &status))
    return status;
  
  delxy = fabs(delxy);
  NGrid = naxes[0];
  npixels = NGrid*NGrid;

  printf("delxy = %f\n", delxy); 
  printf("NGrid = %d\n", NGrid);
  
  if((Data  = (double *)calloc(npixels, sizeof(double)))==NULL){
    printf("An attempt to allocate %ld bytes failed.\n", npixels);
    return EXIT_FAILURE;
  }
 
  printf("Reading FITS...\n");
  if ( READ_FITS(argv[1], npixels, Data) )
    return EXIT_FAILURE;

  printf("Calculating power spectra...\n");
  if ( CALPOW_2D(Data, NGrid, NBin, delxy, argv[2]) )
    return EXIT_FAILURE;

  free(Data);
  return EXIT_SUCCESS;
}
