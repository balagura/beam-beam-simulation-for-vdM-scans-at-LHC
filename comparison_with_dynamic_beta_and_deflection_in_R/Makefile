beam_beam_model_validation : config.hh config.cpp beam_beam_model_validation.cpp
	g++ -O2 -std=c++14 -I. beam_beam_model_validation.cpp config.cpp -o beam_beam_model_validation \
          -lboost_filesystem -lboost_system -lboost_iostreams -lz

%/summary.txt : beam_beam_model_validation # %.txt
	rm -rf $(@D)
	beam_beam_model_validation $<
%/pdf/gg.cor.pdf : comparison.R %/summary.txt
	Rscript --vanilla $< $*
# $< - name of first prerequisite
# $* = stem = % in target


