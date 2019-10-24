compare_with_VdmBBDeflect : ../config.hh ../config.C \
                            compare_with_VdmBBDeflect.C \
                            VdmBBDeflect_Bassetti_Erskine_part.C ../E_field.hh ../Faddeeva.hh ../Faddeeva.cc
	g++ -O2 -std=c++14 -I.. -I. ../config.C \
                            compare_with_VdmBBDeflect.C \
                            VdmBBDeflect_Bassetti_Erskine_part.C ../Faddeeva.cc \
                            -o compare_with_VdmBBDeflect \
          -lboost_filesystem -lboost_system -lboost_iostreams -lz
