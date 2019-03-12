bb=../beam_beam_
bbmad=${bb}MAD-X_2012
sum=summary.txt
subtr_pdf=pdf/gg.cor.dipole.subtr.and.added.pdf

${bb}%/pdf/gg.cor.pdf : comparison.R ${bb}%/${sum}
	Rscript --vanilla $< ${bb}$*
# $< - name of first prerequisite
# $* = stem = % in target

${bbmad}_precise.minus.%/${subtr_pdf} : precise_minus_dipole.R \
                                        ${bbmad}_precise.minus.%/${sum} \
                                        ${bbmad}_%/${sum} \
                                        ${bbmad}/${sum}
	Rscript --vanilla $< ${bbmad}_precise.minus.$*
