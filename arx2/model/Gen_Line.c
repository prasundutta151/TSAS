# include <string.h>
# include <stdio.h>
# include <stdlib.h>
# include <fitsio.h>
# include <math.h>
# include <unistd.h>

# define EXP_NARG 4

void printerror(int status){
  if (status){
    fits_report_error(stderr, status); // Print error report
    exit( status );    // Terminate the program, returning error status
  }
}
int main(int argc, char *argv[]){

  int ii, npixels, nfound;
  float Multi1, Multi2;
  float *Data1, *Data2;
  fitsfile *fptr1,*fptr2, *fptro;
  long naxis1, *naxes1, naxis2, *naxes2;
  int status=0, anynul=0;
  float nulval=0.;
  char comment[128];
  if(argc!=EXP_NARG){
    printf("\nUSAGE : %s <CONFITS> <TAUFITS> <LINFITS)>\n\n", argv[0]);

    printf("OUTPUT : <outfile> = Multi1 * <FITS FILE1> + Multi2 * <FITS FILE2>\n\n");
    printf("CONFITS: Continium FITS file\n");
    printf("TAUFITS: Optical Depth FITS file\n");
    printf("LINFITS: Line Channel FITS file\n\n");
  
    return EXIT_FAILURE;
  }

  if(access(argv[1], F_OK)!=0){
    fprintf(stderr, "File %s does not exists\n", argv[1]);
    return EXIT_FAILURE;
  }
  if(access(argv[2], F_OK)!=0){
    fprintf(stderr, "File %s does not exists\n", argv[2]);
    return EXIT_FAILURE;
  }
  if(access(argv[3], F_OK)==0){
    fprintf(stderr, "File %s exists\n", argv[3]);
    return EXIT_FAILURE;
  }

  if(fits_open_file(&fptr1, argv[1], READONLY, &status))
    printerror(status);
  if(fits_open_file(&fptr2, argv[2], READONLY, &status))
    printerror(status);
  if(fits_create_file(&fptro, argv[3], &status))
    printerror(status);
  if(fits_copy_header(fptr1, fptro, &status))
    printerror(status);

  if(fits_read_key_lng(fptr1, "NAXIS", &naxis1, comment, &status))
    printerror(status);
  naxes1 = (long *)calloc(naxis1, sizeof(long));
  if(fits_read_keys_lng(fptr1, "NAXIS", 1, naxis1, naxes1, &nfound, &status))
    printerror(status);

  if(fits_read_key_lng(fptr2, "NAXIS", &naxis2, comment, &status))
    printerror(status);
  naxes2 = (long *)calloc(naxis2, sizeof(long));
  if(fits_read_keys_lng(fptr2, "NAXIS", 1, naxis2, naxes2, &nfound, &status))
    printerror(status);
  
  if(naxis1!=naxis2){
    fprintf(stderr, "Mismatch in Axes\n");
    return EXIT_FAILURE;
  }
  npixels = 1;
  for(ii=0;ii<naxis1; ii++){

    if(naxes1[ii]!=naxes2[ii]){
      fprintf(stderr, "Mismatch in Axes %d\n", ii);
      return EXIT_FAILURE;
    }
    else 
      npixels *= naxes1[ii];
  }
  
  Data1 = (float *)calloc(npixels, sizeof(float));
  Data2 = (float *)calloc(npixels, sizeof(float));
  
  if(fits_read_img(fptr1, TFLOAT, 1, npixels, &nulval, Data1, &anynul, &status))
    printerror(status);
  if(fits_read_img(fptr2, TFLOAT, 1, npixels, &nulval, Data2, &anynul, &status))
    printerror(status);
  
  for(ii=0; ii<npixels; ii++)
    Data1[ii] = Data1[ii]*exp(-Data2[ii]);
  
  if(fits_write_img(fptro, TFLOAT, 1, npixels, Data1, &status))
    printerror(status);
  if(fits_close_file(fptr1, &status))
    printerror(status);
  if(fits_close_file(fptr2, &status))
    printerror(status);
  if(fits_close_file(fptro, &status))
    printerror(status);
  return 0;
}
