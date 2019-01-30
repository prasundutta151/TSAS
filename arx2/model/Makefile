exec=Gen_Line Tau_Gen PowFits2Dd

LINKLIB=-lm -lnrcp -lfftw3 -lfftw3f -lcfitsio  -fopenmp
CC=gcc 

%: %.c $(DEPS)
	$(CC) -o  $@ $< ${LINKLIB}

all:
	make clean
	make ${exec}

clean:
	rm -rf ${exec} *.o *~

cleanfits:
	rm -rf contmod.fits

cleanpow:
	rm -rf *.pow

cleanall:
	make cleanfits
	make cleanpow
	rm -rf PowFits2D_dbl Gauss_2D_FITS_ARB
