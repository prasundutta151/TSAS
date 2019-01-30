
ALPHA=$1
SIGMA=$2
TAU0=$3
SEED1=$4
IREAL=$5

CONFITS=$(awk 'NR==2 {print $1}' tau.input)
OBSFILE=$(awk 'NR==3 {print $1}' tau.input)
TEMPDIR=$(awk 'NR==5 {print $1}' tau.input)
UMIN=$(awk 'NR==7 {print $1}' tau.input)
UMAX=$(awk 'NR==7 {print $2}' tau.input)
NBIN=$(awk 'NR==8 {print $2}' tau.input)
NGRID=$(awk 'NR==8 {print $1}' tau.input)

#IREAL=${ALPHA}_${SIGMA}_${TAU0}_${IREAL}
OBSFILE=ILINE_60_tau_1.094.dat
TAUFITS=${TEMPDIR}/REAL/TAU_${IREAL}.FITS
LINFITS=${TEMPDIR}/REAL/LIN_${IREAL}.FITS
LINPOW=${TEMPDIR}/REAL/LIN_${IREAL}.POW


echo $SEED >> aa.seed
run_funcs(){

srname=$1
echo "*************************************************"
echo "Running $srname"
echo "$srname $code"
echo "*************************************************"
$srname $code
status=$(echo $?)
if (( $status != 0 ))
then
    echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    echo "ERROR in $srname"
    echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    echo
    exit 1
else
    echo "Run successful"
    echo "***********************************************"
fi
}
echo Contfile $CONFITS


code="$TAUFITS $NGRID $TAU0 $SIGMA $ALPHA $SEED1"
run_funcs "./Tau_Gen"
code="$CONFITS $TAUFITS $LINFITS"
run_funcs "./Gen_Line"
code="$LINFITS $LINPOW $NBIN"
run_funcs "./PowFits2Dd"

rm -rf $TAUFITS $LINFITS


