all: BxB_python.so \
     beam_beam \
     beam_beam_model_validation

beam_beam : beam_beam.o bb_config.o bb.o config.o multi_gaussian.o interpolators.o output.o Faddeeva.o gzstream.o
	g++ -O2 -fPIC $^ -o $@ -lz -pthread

BxB_python.so : bb_C.o bb.o multi_gaussian.o interpolators.o output.o Faddeeva.o gzstream.o
	g++ -O2 -fPIC -shared $^ -o $@ -lz -pthread

bb_C.o: bb_C.cpp bb_C.h bb_C_CPP.h bb.hh
	g++ -c -O2 -fPIC -shared $< -o $@
bb_config.o: bb_config.cpp bb_config.hh bb.hh bb_C_CPP.h config.hh
	g++ -c -O2 -fPIC -shared $< -o $@
bb.o: bb.cpp gzstream.hh multi_gaussian.hh interpolators.hh output.hh bb.hh bb_C_CPP.h
	g++ -c -O2 -fPIC -shared $< -o $@
beam_beam.o: beam_beam.cpp bb_config.hh bb.hh bb_C_CPP.h config.hh
	g++ -c -O2 -fPIC -shared $< -o $@
config.o: config.cpp config.hh
	g++ -c -O2 -fPIC -shared $< -o $@
Faddeeva.o: Faddeeva.cpp Faddeeva.hh
	g++ -c -O2 -fPIC -shared $< -o $@
gzstream.o: gzstream.cpp gzstream.hh
	g++ -c -O2 -fPIC -shared $< -o $@
interpolators.o: interpolators.cpp interpolators.hh
	g++ -c -O2 -fPIC -shared $< -o $@
multi_gaussian.o: multi_gaussian.cpp multi_gaussian.hh interpolators.hh \
 E_field.hh Faddeeva.hh
	g++ -c -O2 -fPIC -shared $< -o $@
output.o: output.cpp output.hh gzstream.hh bb.hh bb_C_CPP.h
	g++ -c -O2 -fPIC -shared -std=c++17 $< -o $@

beam_beam_model_validation.o: beam_beam_model_validation.C config.hh
	g++ -c -O2 -std=c++14 $< -o $@
beam_beam_model_validation : beam_beam_model_validation.o config.o
	g++ -O2 beam_beam_model_validation.o config.o -o beam_beam_model_validation \
          -lboost_filesystem -lboost_system -lboost_iostreams -lz
