inputs = $(wildcard ../beam_beam_*.txt)
dir_targets = $(patsubst ../beam_beam_%.txt, ../beam_beam_%, $(inputs))
plots = $(patsubst %, %/pdf/gg.cor.pdf, $(dir_targets))

all: $(plots)

../beam_beam_%/pdf/gg.cor.pdf : ../beam_beam_% comparison.R
	Rscript --vanilla comparison.R $<
