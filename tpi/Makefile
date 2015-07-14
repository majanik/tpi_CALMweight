CC	= g++
LD	= g++

CFLAGS	= `root-config --cflags` -O3 -g -I../SphericalHarmonics -I../therminator2_GENBOD/v12_noQSB_normal_minijets/build/include
LDFLAGS = `root-config --libs` -g -lgsl -lgslcblas

all: tpi

tpi: tpi.o PairReader.o PairWeight.o ExpCF3D.o ExpCF1D.o ExpCFSH.o ExpCFEP.o SourceMonitor.o TChainProxy.o ../SphericalHarmonics/ylm.o ../SphericalHarmonics/CorrFctnDirectYlm.o ../therminator2_GENBOD/v12_noQSB_normal_minijets/build/obj/ParticleCoor.o
	$(LD) $^ -o $@ $(LDFLAGS)

%.o: %.cxx
	$(CC) $^ -o $@ $(CFLAGS) -c 

clean: 
	rm -f *.o tpi twointeg_sh_mr test
