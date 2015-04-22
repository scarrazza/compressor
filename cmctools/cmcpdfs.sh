#!/bin/bash

##############################################################
#
# Script to generate a CMC-PDF set
#
##############################################################

# LHAPDF main path`
lhapdfdir=`lhapdf-config --datadir`

# Some cleanup to begin with
rm -rf $lhapdfdir/CMCPDF*

##
# Define sets in the combination
# Using most updated sets
sets=("NNPDF30" "CT14" "MMHT14" "CMCPDFcomb_")

for order in "nnlo"
do
    # Some cleanup to begin with
    echo ${sets[3]}$order
    rm -rf $lhapdfdir${sets[3]}$order
    mkdir $lhapdfdir${sets[3]}$order
    ls $lhapdfdir${sets[3]}$order

    # Copy also central replica
    # but this should never be used, only for consistency with LHAPDF6 standard
    # At the end of compression, replace with average over compressed sets
    if [ "${sets[0]}" == "NNPDF30" ]; then
	cp $lhapdfdir${sets[0]}_${order}_as_0118/${sets[0]}_${order}_as_0118_0000.dat  $lhapdfdir${sets[3]}${order}/${sets[3]}${order}_0000.dat
    fi
    
    # Set counter for replicas
    krep=1
    
    for PDFset in ${sets[0]} ${sets[1]} ${sets[2]}
    do
	
	Nrep=0
	pdfsetname="test"
	if [ "$PDFset" == "NNPDF30" ]; then
	    Nrep=100
	    if [ "$order" == "nnlo" ]; then
		pdfsetname="NNPDF30_nnlo_as_0118"
	    fi
	fi
	if [ "$PDFset" == "CT14" ]; then
	    Nrep=100
	    if [ "$order" == "nnlo" ]; then
		pdfsetname="CT14pre-1148q_rand1002"
	    fi
	fi
	if [ "$PDFset" == "MMHT14" ]; then
	    Nrep=100
	    if [ "$order" == "nnlo" ]; then
		pdfsetname="MMHT2014nnlo68cl_rand1001"
	    fi
	fi
	
	
	echo "PDFset = " $PDFset
	echo "pdfsetname = " $pdfsetname 
	echo "Number of replicas = " $Nrep
	
	for irep in `seq 1 $Nrep`;
	do
	    zeroi="0000";
	    if [ $irep -le 9 ]; then
		zeroi="000"
	    fi
	    if [ $irep -le 99 -a $irep -ge 10 ]; then
	    zeroi="00"
	fi
	if [ $irep -le 999 -a $irep -ge 100 ]; then
	    zeroi="0"
	fi
	if [ $irep -le 9999 -a $irep -ge 1000 ]; then
	    zeroi=""
	fi
	
	zerok="0000";
	if [ $krep -le 9 ]; then
	    zerok="000"
	fi
	if [ $krep -le 99 -a $krep -ge 10 ]; then
	    zerok="00"
	fi
	if [ $krep -le 999 -a $krep -ge 100 ]; then
	    zerok="0"
	fi
	if [ $krep -le 9999 -a $krep -ge 1000 ]; then
	    zerok=""
	fi
	
	#ls $pdfsetname/${pdfsetname}_${zeroi}$irep.dat
	cp $lhapdfdir$pdfsetname/${pdfsetname}_${zeroi}$irep.dat $lhapdfdir${sets[3]}${order}/${sets[3]}${order}_${zerok}$krep.dat
	#ls ${sets[3]}${order}/${sets[3]}${order}_${zerok}$krep.dat
	krep=$(($krep+1))
    done
    
    echo " -------------------------------------------------------------- "
    echo " New PDF set added to the combined set" ${sets[3]}${order}
    echo " Set used = " $pdfsetname
    echo " Number of replicas = " $Nrep
    echo " Number of replicas currently in CMCPDF = " $(($krep-1))
    echo " -------------------------------------------------------------- "
    
    done
    
done

# Finally copy info file
cp info/${sets[3]}${order}.info $lhapdfdir${sets[3]}${order}/
ls $lhapdfdir${sets[3]}${order}/${sets[3]}${order}.info
