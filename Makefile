beam_beam : config.hh config.C beam_beam.C
	g++ -O2 -std=c++14 -I. beam_beam.C config.C -o beam_beam \
          -lboost_filesystem -lboost_system -lboost_iostreams -lz
