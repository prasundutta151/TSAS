rm -rf *.fits
rm -rf *.pow

BMAJ=0.01153328683641
BMIN=0.01025484191047

make all

DELIF=delImod.fits
DELIP=delImod.pow
CONTD=CRABCONT.FITS
CPMOD=CRAB_CONT350_64.POW
CPMOD=ModPlawCont.POW
./Gauss_2D_FITS_ARB $CPMOD $DELIF -1232 ${BMAJ} ${BMIN}
./PowFits2D_dbl $DELIF $DELIP 32

ls -ltrh
rm -rf test.pow


ipython Gen_ContMap.py
ipython Gen_ContMap.py $CONTD $DELIF ellipse.reg contmod.fits 30. 350

ds9 contmod.fits II.fits Win.fits delI.fits &

./PowFits2D_dbl CRABCONT.FITS CRABCONT.pow 32
./PowFits2D_dbl contmod.fits contmod20.pow 32
./PowFits2D_dbl II.fits II.pow 32
./PowFits2D_dbl Win.fits Win.pow 32
./PowFits2D_dbl delI.fits delI.pow 32
