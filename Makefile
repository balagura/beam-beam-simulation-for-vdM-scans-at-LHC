inputs = $(wildcard beam_beam_*.txt)
dir_targets = $(patsubst beam_beam_%.txt, beam_beam_%, $(inputs))
summary_targets = $(patsubst beam_beam_%, beam_beam_%/summary.txt, $(dir_targets))

all: beam_beam $(summary_targets)

beam_beam : config.hh config.C beam_beam.C
	g++ -O2 -std=c++14 -I. beam_beam.C config.C -o beam_beam \
          -lboost_filesystem -lboost_system -lboost_iostreams -lz

./beam_beam_%/summary.txt : beam_beam_%.txt beam_beam
	rm -rf $(@D)
	beam_beam $<
