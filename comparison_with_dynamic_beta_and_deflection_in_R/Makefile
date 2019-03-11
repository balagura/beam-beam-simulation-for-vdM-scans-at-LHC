inputs = $(wildcard ../beam_beam_MAD-X*.txt)
dir_targets = $(patsubst ../beam_beam_%.txt, ../beam_beam_%, $(inputs))
plots = $(patsubst %, %/pdf/gg.cor.pdf, $(dir_targets))

all: $(plots) ../beam_beam_/pdf/gg.cor.orbit.subtr.and.added.pdf

../beam_beam_%/pdf/gg.cor.pdf : ../beam_beam_% comparison.R
	Rscript --vanilla comparison.R $<

../beam_beam_/pdf/gg.cor.orbit.subtr.and.added.pdf : ../beam_beam_MAD-X_2012/Rdata/results.rds precise_minus_const.R
	Rscript --vanilla precise_minus_const.R $<
