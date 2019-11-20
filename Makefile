# inputs = $(wildcard beam_beam_*.txt)
# dir_targets = $(patsubst beam_beam_%.txt, beam_beam_%, $(inputs))
# summary_targets = $(patsubst beam_beam_%, beam_beam_%/summary.txt, $(dir_targets))

# all: beam_beam $(summary_targets)

all: beam_beam \
     beam_beam_model_validation

beam_beam : beam_beam.o config.o multi_gaussian.o interpolators.o n.o output.o Faddeeva.o gzstream.o
	g++ -O2 beam_beam.o config.o multi_gaussian.o interpolators.o \
                            n.o output.o Faddeeva.o gzstream.o -o beam_beam -lz -pthread

beam_beam.o : beam_beam.C config.hh multi_gaussian.hh n.hh output.hh gzstream.h
	g++ -c -O2 -std=c++17 -I. $< -o $@
multi_gaussian.o : multi_gaussian.C multi_gaussian.hh interpolators.hh E_field.hh
	g++ -c -O2 -std=c++17 -I. $< -o $@
output.o : output.C output.hh config.hh n.hh gzstream.h
	g++ -c -O2 -std=c++17 -I. $< -o $@
interpolators.o : interpolators.C interpolators.hh
	g++ -c -O2 -std=c++17 -I. $< -o $@
n.o : n.C n.hh config.hh
	g++ -c -O2 -std=c++17 -I. $< -o $@
config.o : config.C config.hh
	g++ -c -O2 -std=c++17 -I. $< -o $@
Faddeeva.o : Faddeeva.cc Faddeeva.hh
	g++ -c -O2 -std=c++17 -I. $< -o $@
gzstream.o : gzstream.C gzstream.h
	g++ -c -O2 -std=c++17 -I. $< -o $@

beam_beam_model_validation.o : beam_beam_model_validation.C config.hh
	g++ -c -O2 -std=c++14 -I. $< -o $@
beam_beam_model_validation : beam_beam_model_validation.o config.o
	g++ -O2 beam_beam_model_validation.o config.o -o beam_beam_model_validation \
          -lboost_filesystem -lboost_system -lboost_iostreams -lz

./beam_beam_%/summary.txt : beam_beam_%.txt beam_beam
	rm -rf $(@D)
	beam_beam $<
