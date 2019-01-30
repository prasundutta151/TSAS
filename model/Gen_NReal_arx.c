# include <string.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <omp.h>
# include <unistd.h>

# define NARG (2)

int main(int argc, char *argv[]){

  int ii, nthreads, Tid, NReal;
  float alp;
  char command[256], inputfile[128], param[128];
  FILE *fp;
  if(argc!=NARG){
    
    fprintf(stderr, "USAGE: ./argv[0] <NThreads>\n", argv[0]);
    return EXIT_FAILURE;
  }
  
  sprintf(inputfile, "tau.input"); 
  nthreads = (int)atoi(argv[1]);
  
  if(access(inputfile, F_OK)!=0){
    fprintf(stderr, "File %s does not exists\n", inputfile);
    return EXIT_FAILURE;
  }

  omp_set_num_threads(nthreads);
  
  float alp_m, alp_M, sig_m, sig_M, ta0_m, ta0_M;
  float dalp, dsig, dta0, alpha, sigma, tau0;
  int ialp, isig, ita0;
  int Nalp, Nsig, Nta0;
  
  fp=fopen(inputfile,"r");
  
  // Getting passed first seven lines
  for(ii=0; ii<7; ii++)
    fgets(param,100, fp);

  fgets(param,100, fp);
  sscanf(param,"%f%f%d",&alp_m, &alp_M, &Nalp);

  fgets(param,100, fp);
  sscanf(param,"%f%f%d",&sig_m, &sig_M, &Nsig);
  
  fgets(param,100, fp);
  sscanf(param,"%f%f%d",&ta0_m, &ta0_M, &Nta0);
 
  fclose(fp);

  dalp = (alp_M-alp_m)/(1.*Nalp);
  dsig = (sig_M-sig_m)/(1.*Nsig);
  dta0 = (ta0_M-ta0_m)/(1.*Nta0);
  
  NReal = Nalp*Nsig*Nta0;
  
#pragma omp parallel for private(Tid, ii, command, alpha, sigma, tau0, ialp, isig, ita0)
  for(ii=0; ii<NReal; ii++){

    ialp  = ii%Nalp;
    alpha = alp_m + dalp*ialp;
    isig  = (ii/Nalp)%Nsig;
    sigma = sig_m + dsig*isig;
    ita0  = (ii/Nalp)/Nsig;
    tau0  = ta0_m + dta0*ita0;

    Tid = omp_get_thread_num();
    printf("%d\t%d\t%d\t%d\n", Tid, ialp, isig, ita0);
    sprintf(command, "sh GenLinePs.sh %.2f %.2f %.2f", alpha, sigma, tau0);
    system(command);
  }

  return 0;
}
