# inputs = $(wildcard beam_beam_*.txt)
# dir_targets = $(patsubst beam_beam_%.txt, beam_beam_%, $(inputs))
# summary_targets = $(patsubst beam_beam_%, beam_beam_%/summary.txt, $(dir_targets))

# all: beam_beam $(summary_targets)

all: beam_beam \
     beam_beam_model_validation

beam_beam : beam_beam.C \
            config.hh E_field.hh multi_gaussian.hh interpolators.hh n.hh output.hh Faddeeva.hh gzstream.h \
            config.C             multi_gaussian.C  interpolators.C  n.C  output.C  Faddeeva.cc gzstream.C
	g++ -O2 -std=c++17 -I. beam_beam.C config.C multi_gaussian.C interpolators.C \
                               n.C output.C Faddeeva.cc gzstream.C -o beam_beam -lz -pthread

beam_beam_model_validation : config.hh config.C beam_beam_model_validation.C
	g++ -O2 -std=c++14 -I. beam_beam_model_validation.C config.C -o beam_beam_model_validation \
          -lboost_filesystem -lboost_system -lboost_iostreams -lz

./beam_beam_%/summary.txt : beam_beam_%.txt beam_beam
	rm -rf $(@D)
	beam_beam $<
