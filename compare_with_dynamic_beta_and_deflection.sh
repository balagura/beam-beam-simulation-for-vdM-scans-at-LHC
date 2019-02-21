#!/bin/bash

file="beam_beam_MAD-X_2012.txt"

extension="${file##*.}"
filename="${file%%.*}"

kicks="precise const quadrupole const.and.quadrupole"
for kick in ${kicks}
do
    # add _kick prefix except for precise
    if [ "${kick}" == "precise" ]
    then
	dir=${filename}
    else
	dir=${filename}_${kick}
    fi
    f=${dir}.${extension}
    # create not git-controlled *_kick.txt files if they do not exist
    if [ ! -f "${f}" ]
    then
	# substitute eg. kick=precise -> kick=const etc.
	sed "/^[ \t]*kick.model[ \t]*=[ \t]*precise[ \t]*$/ s/precise/${kick}/" ${file} > ${f}
    fi
    echo "========================================"
    echo "        Simulate ${kick}"
    echo "========================================"
    make ${dir}/summary.txt
    cd comparison_with_dynamic_beta_and_deflection_in_R
    make ../${dir}/pdf/gg.cor.pdf
    cd ..
    # show final plot
    xdg-open ${dir}/pdf/gg.cor.pdf &
done
