#
# B*B simulation of beam-beam effects in van der Meer scans at LHC.
# Copyright (C) 2019 Vladislav Balagura (balagura@cern.ch)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------
#
# libBxB is needed for python binding and for beam_beam() calls in C++
#
# beam_beam is a standalone C++ simulation reading input parameters from
# "Config" file
#
# beam_beam_model_validation was used to validate the B*B model in the simple
# case of single IP, round single Gaussian bunches
#
#
all: libBxB.so \
     beam_beam \
     beam_beam_model_validation

beam_beam : beam_beam.o bb_config.o bb.o config.o multi_gaussian.o interpolators.o output.o Faddeeva.o gzstream.o
	g++ -O2 -fPIC $^ -o $@ -lz -pthread

libBxB.so : bb_C.o bb.o multi_gaussian.o interpolators.o output.o Faddeeva.o gzstream.o
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
