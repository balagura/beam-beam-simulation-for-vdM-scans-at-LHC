#!/bin/bash

file="beam_beam_MAD-X_2012.txt"

extension="${file##*.}"
filename="${file%%.*}"

kicks="precise precise.minus.average average precise.minus.one.third.of.average "\
"one.third.of.average precise.minus.at.center at.center quadrupole average.and.quadrupole"

n_transional_turns="1 1000"

for kick in ${kicks}
do
    for ntrans in ${n_transional_turns}
    do
	dir=${filename}_${kick}_${ntrans}
	f=${dir}.${extension}
	# substitute
	# kick=precise -> kick=const 
	# N.transitional.turns=1000 -> N.transitional.turns=1 etc.
	sed -e "/^[ \t]*kick.model[ \t]*=[ \t]*precise[ \t]*$/ s/precise/${kick}/" \
	    -e "/^[ \t]*N.transitional.turns[ \t]*=[ \t]*1000[ \t]*$/ s/1000/${ntrans}/" \
	    ${file} > ${f}
	echo "======================================================================"
	echo "        Simulate ${kick} with ${ntrans} transitional turns"
	echo "======================================================================"
	make ${dir}/summary.txt
	cd comparison_with_dynamic_beta_and_deflection_in_R
	make ../${dir}/pdf/gg.cor.pdf
	# show final plots
	xdg-open ../${dir}/pdf/gg.cor.pdf &
	if [[ ${kick} == *precise.minus.* ]]; then
	    make ../${dir}/pdf/gg.cor.dipole.subtr.and.added.pdf
	    # show final plots
	    xdg-open ../${dir}/pdf/gg.cor.dipole.subtr.and.added.pdf &
	fi
	cd ..
    done
done
