// This program reads a single channel single stokes visibility data, grid them and write in FITS image format  [header copied from a image fits file]

// negative values in grid !! still puzzling !!

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <fitsio.h>
# include <unistd.h>

# define EXP_ARG (5)

enum{YES=0, NO=1};

void printerror(int status){
  if (status){
    fits_report_error(stderr, status); 
    exit( status );  
  }
}

int FCHECK(char filename[], int fstat){
    
  if(fstat==YES){
    if(access(filename, F_OK)!=0){
      
      fprintf(stderr, "\n ERROR: Input File %s does not exists\n\n", filename);
      return EXIT_FAILURE;
    }
  }
  else if(fstat==NO){ 
    if(access(filename, F_OK)==0){
      
     
      fprintf(stderr, "\n ERROR: Output File %s already exists\n\n", filename);
      return EXIT_FAILURE;
    }
  }
  else {
    fprintf(stderr, "\n ERROR: Wong file status: %d [options = YES, NO]\n\n", fstat);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

void HLINE(char cc){
  int ii;
  for(ii=0;ii<80;ii++)printf("%c", cc); printf("\n");
}
int main(int argc, char *argv[]){
 
  char INFITS[128], OUTFITS[128];
  double Umax, dU;
  int Nvis, Nv, Nu, Ng;
  double *Vre, *Vim, *Vsq, *Bln;
  float *data;
  fitsfile *fptr, *fptro;  
  int status, anynul;
  int iu, iv, induv, dinuv;
  long gcount, pcount, naxis, nstokes, nchan, ncmplx;
  long el1, group;

  char key_simple[FLEN_KEYWORD]="SIMPLE";
  char key_bitpix[FLEN_KEYWORD]="BITPIX";
  char key_naxis[FLEN_KEYWORD]="NAXIS";
  char key_naxis1[FLEN_KEYWORD]="NAXIS1";
  char key_naxis2[FLEN_KEYWORD]="NAXIS2";
  char key_naxis3[FLEN_KEYWORD]="NAXIS3";
  char key_naxis4[FLEN_KEYWORD]="NAXIS4";
  char key_gcount[FLEN_KEYWORD]="GCOUNT"; 
  char key_pcount[FLEN_KEYWORD]="PCOUNT";
  char key_history[FLEN_KEYWORD]="HISTORY";
  char keyvalue[FLEN_VALUE];
  char comment[FLEN_COMMENT];

  float *randpar, *cmplx, uu, vv;
  float nulval, Uval, signv, dmax, dmin;
  float Umax_data, Umin_data;
 
  el1    = 1;
  anynul = 0;
  nulval = 0.;

  if(argc!=EXP_ARG){

      printf("Usage: %s <input (uv) FITS file>  <output Grid Fits File> <Umax (Klam) > <dU (KLam)>\n", argv[0]);
      return 1;
    }
  
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[2],"%s",OUTFITS);

  Umax = (double) atof(argv[3]);
  dU   = (double) atof(argv[4]);

  if(FCHECK(INFITS, YES)!=0)
    return EXIT_FAILURE;

  if(FCHECK(OUTFITS, NO)!=0)
    return EXIT_FAILURE;

   Ng =(int) ceil(Umax/dU); 
   Nu = 2*Ng +2;            
   Nv = Ng + 2;    
   Nvis = (long) Nu*Nv; 
  
   HLINE('=');
   printf("Ng   = %d\t", Ng);
   printf("Nv   = %d\t", Nv);
   printf("Nu   = %d\t", Nu);
   printf("Nvis = %d\n", Nvis);
   HLINE('-');

   Vre = (double *)calloc(Nvis, sizeof(double));
   Vim = (double *)calloc(Nvis, sizeof(double));
   Vsq = (double *)calloc(Nvis, sizeof(double));
   Bln = (double *)calloc(Nvis, sizeof(double));
   data = (float *)calloc(Nvis, sizeof(float));
   status = 0;
   if(fits_open_file(&fptr, INFITS, READONLY, &status))
     printerror(status);

   if(fits_read_keyword(fptr,key_simple,keyvalue,comment,&status))
     printerror(status);
   if(*keyvalue!='T'){ 
     fprintf(stderr,"\n ERROR: Input File is NOT a SIMPLE FITS file\n");
     return EXIT_FAILURE;
   }
   if(fits_read_key_lng(fptr,key_naxis2,&ncmplx,comment,&status))
     printerror( status );
   printf("ncmplx = %ld\t", ncmplx);

   if(fits_read_key_lng(fptr,key_naxis3,&nstokes,comment,&status))
     printerror( status );
   printf("nstokes = %ld\t", nstokes);
   
   if(fits_read_key_lng(fptr,key_naxis4,&nchan,comment,&status))
     printerror( status );
   printf("nchan = %ld\t",nchan);
   
   if(fits_read_key_lng(fptr,key_gcount,&gcount,comment,&status))
     printerror( status );
   printf("GCOUNT = %ld\t",gcount);
   
   if(fits_read_key_lng(fptr,key_pcount,&pcount,comment,&status))
     printerror( status );
   printf("PCOUNT = %ld\n",pcount);
  
   HLINE('-');

   if(nstokes!=1){
     fprintf(stderr, "\n ERROR: NSTOKES =/= 1\n\n");
     return EXIT_FAILURE;
   }
  
  if(nchan!=1){
    fprintf(stderr, "\n ERROR: NCHAN =/= 1\n\n");
    return EXIT_FAILURE;
  }
  
  randpar = (float *)calloc(gcount, sizeof(float));
  cmplx   = (float *)calloc(ncmplx, sizeof(float));

  Umax_data = -1.e30;
  Umin_data =  1.e30;

  printf("\nGridding the visibilities\n");
  for(group=1; group<=gcount; group++){

    if(fits_read_grppar_flt(fptr, group, el1, pcount, randpar, &status))
      printerror( status );

    signv= (randpar[1]<0.) ? -1. : 1. ;
	  
    uu = signv*randpar[0]*1.e-3;
    vv = signv*randpar[1]*1.e-3;
    
    Uval = sqrt(uu*uu + vv*vv);
   
    Umax_data = (Umax_data > Uval) ? Umax_data : Uval;
    Umin_data = (Umin_data < Uval) ? Umin_data : Uval; 

    if( Uval < Umax ){
      
      iu = (int)roundf(uu/dU); 
      iu = (iu<0) ? Nu+iu : iu;
      iv = (int)roundf(vv/dU); 
      
      induv = iv*Nu + iu;

      status = 0;
      anynul = 0;
      nulval = 0.;
      el1    = 1;

      if(fits_read_img_flt(fptr, group, el1, ncmplx, nulval, cmplx, &anynul, &status))
	printerror( status );

      if(cmplx[2]>0.9){ //Checking if the wt is 1.
	
	Vre[induv] += cmplx[0];
	Vim[induv] += cmplx[1];
	Vsq[induv] += (cmplx[0]*cmplx[0] + cmplx[1]*cmplx[1]);
	Bln[induv] += 1.;
      }
    }
  }
  if(fits_close_file(fptr, &status))
    printerror ( status) ;
  printf("\nUmin = %.2e\tUmax = %.2e\n", Umin_data, Umax_data);  
  printf("\nEnsuring reality\n");

  for(iv=0; iv<Nv; iv=iv+Nv-1)
    for(iu=1; iu<=Ng; ++iu){
      
      induv = (long)iv*Nu +iu;
      dinuv = (long)iv*Nu+(Nu-iu);
      
      Vre[induv] += Vre[dinuv];
      Vim[induv] -= Vim[dinuv];
      Vsq[induv] += Vsq[dinuv];
      Bln[induv] += Bln[dinuv];
      
      Vre[induv] *= 0.5;
      Vim[induv] *= 0.5;
      Vsq[induv] *= 0.5;
      Bln[induv] *= 0.5;
      
      Vre[dinuv] = Vre[induv];
      Vim[dinuv] = -1.*Vim[induv];
      Vsq[dinuv] = Vsq[induv];
      Bln[dinuv] = Bln[induv];   
    }
  for(iu=0;iu<Nv;iu=iu+Nv-1)
	for(iv=0;iv<Nv;iv=iv+Nv-1){

	  induv=(long)iv*Nu+iu;
	  Vim[induv] = 0.;
	}
  
  
  printf("\nEstimating squares\n");
  dmax = -1.e30;
  dmin = 1.e30;  
  for(iv=0; iv<Nv; iv++)
    for(iu=0; iu<Nu; iu++){
      
      induv = (long)iv*Nu +iu;
      
      data[induv] = Vre[induv]*Vre[induv] + Vim[induv]*Vim[induv];
      data[induv] = data[induv] - Vsq[induv];
      data[induv] = (data[induv] > 0. ) ? data[induv] : 0.;
      //data[induv] = Bln[induv];
      data[induv] = (float )((Bln[induv] > 1.) ? data[induv]/(Bln[induv]*(Bln[induv]-1.)): 0.);

      dmax = (dmax>data[induv]) ? dmax : data[induv];
      dmin = (dmin<data[induv]) ? dmin : data[induv];
    }

  printf("\nComputations done\n");
  HLINE('-');

  char key_cdelt1[FLEN_KEYWORD]="CDELT1";
  char key_cdelt2[FLEN_KEYWORD]="CDELT2";
  long *naxes;

  status = 0;
  naxis  = 2;
  naxes = (long *)calloc(naxis, sizeof(long));
  naxes[0] = Nu;
  naxes[1] = Nv;
 
  sprintf(comment, "\\");
  if(fits_create_file(&fptro, OUTFITS, &status))
    printerror(status);
  if(fits_create_img(fptro, DOUBLE_IMG, naxis, naxes, &status))
    printerror(status);
  if(fits_write_key(fptro, TDOUBLE, key_cdelt1, &dU, comment, &status))
    printerror( status );
  if(fits_write_key(fptro, TDOUBLE, key_cdelt2, &dU, comment, &status))
    printerror( status );


  printf("\nWriting to fits\n");
  if(fits_write_img(fptro, TFLOAT, 1, Nvis, data, &status))
    printerror(status);
  if(fits_close_file(fptro, &status))
    printerror(status);

  free(Vre);
  free(Vim);
  free(Vsq);
  free(Bln);
  printf("\nAll done\n");
  HLINE('=');
  
  return EXIT_SUCCESS;
}
