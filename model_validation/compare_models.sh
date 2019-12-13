#!/bin/bash

file="bb_MAD-X_2012.txt"

extension="${file##*.}"
filename="${file%%.*}"

# precise.average* should be last, since their "partners" are used
# in making summary plots when processing precise.average*
kicks="precise average precise.minus.average one.third.of.average precise.minus.one.third.of.average "\
"precise.minus.at.center at.center average.and.quadrupole quadrupole"

n_transitional_turns="1 1000"

for kick in ${kicks}
do
    for ntrans in ${n_transitional_turns}
    do
	echo "======================================================================"
	echo "        Process ${kick} with ${ntrans} transitional turns"
	echo "======================================================================"
	dir=${filename}_${kick}_${ntrans}
	f=${dir}.${extension}
	# substitute
	# kick=precise -> kick=const 
	# N.transitional.turns=1000 -> N.transitional.turns=1 etc.
	# create file if it does not exist, otherwise do not touch it so that make
	# will not rerun
	test -f ${f} || \
	sed -e "/^[ \t]*kick.model[ \t]*=[ \t]*precise[ \t]*$/ s/precise/${kick}/" \
	    -e "/^[ \t]*N.transitional.turns[ \t]*=[ \t]*1000[ \t]*$/ s/1000/${ntrans}/" \
	    ${file} > ${f}
	make ${dir}/summary.txt
	test -f ${f} && rm -f ${f}
	make ${dir}/pdf/gg.cor.pdf
	# show final plots
	xdg-open ${dir}/pdf/gg.cor.pdf &
	if [[ ${kick} == *precise.minus.* ]]; then
	    # show final plots
	    xdg-open ${dir}/pdf/gg.cor.dipole.subtr.and.added.pdf &
	fi
    done
done
