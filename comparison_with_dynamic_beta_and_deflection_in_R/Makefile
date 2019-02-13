inputs = $(wildcard ../beam_beam_*.txt)
dir_targets = $(patsubst ../beam_beam_%.txt, ../beam_beam_%, $(inputs))
plots = $(patsubst %, %/pdf, $(dir_targets))

all: $(plots)

../beam_beam_%/pdf : ../beam_beam_% comparison.R
	Rscript --vanilla comparison.R $<
