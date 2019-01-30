

SIGMA=1.6e3
ALPHA=-2.8
CONVL=10.
SCALEW=4.9

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

NBIN=64
UMIN=0.1
UMAX=50.

python Gen_powps.py $CPMOD $SIGMA $ALPHA $UMIN $UMAX $NBIN

./Gauss_2D_FITS_ARB $CPMOD $DELIF -1232 ${BMAJ} ${BMIN}
./PowFits2D_dbl $DELIF $DELIP 32

ls -ltrh
rm -rf test.pow


ipython Gen_ContMap.py
ipython Get_ContMap.py $CONTD $DELIF ellipse.reg contmod.fits $CONVL 350 $SCALEW

ds9 contmod.fits II.fits $CONTD delI.fits &

./PowFits2D_dbl CRABCONT.FITS CRABCONT.pow 32
./PowFits2D_dbl contmod.fits contmod.pow 32
./PowFits2D_dbl II.fits II.pow 32
./PowFits2D_dbl Win.fits Win.pow 32
./PowFits2D_dbl delI.fits delI.pow 32
