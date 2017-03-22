#
# Makefile for many Unix compilers using the
# "standard" command name CC
#
CC=g++
DEBUG=-g
CFLAGS=-Wall -c $(DEBUG) $(BOOST)
LFLAGS=-Wall $(DEBUG)
RCFLAGS=`root-config --cflags --libs`


OBJS += ./build/PhotonAreaDistributions.o
OBJS += ./build/PhotonFits.o
OBJS += ./build/SPEFunctions.o
OBJS += ./build/WaveformGenerator.o

all: ./AreaCorrections

./AreaCorrections: $(OBJS) ./build/AreaCorrections.o
	$(CC) $(OSFLAGS) $(RCFLAGS) $(LFLAGS) $(RLFLAGS) $^ $(LIBS) -lMinuit -o $@

./build/AreaCorrections.o: AreaCorrections.cpp $(OBJS)
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@ 

./build/PhotonAreaDistributions.o: ../PhotonTiming/PhotonAreaDistributions_3.cpp ../PhotonTiming/PhotonAreaDistributions.h ../PhotonTiming/SPEparameters.h
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $(RLFLAGS) $< -o $@

./build/PhotonFits.o: ../PhotonTiming/PhotonFits_Brian_test.cpp ../PhotonTiming/PhotonFits.h ./build/PhotonAreaDistributions.o
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $(RLFLAGS) $< -o $@

./build/SPEFunctions.o: ../PhotonTiming/SPEFunctions_Rose_test.cpp ../PhotonTiming/SPEFunctions_Rose_test.h
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@

./build/WaveformGenerator.o: ../WaveformGenerator/WaveformGenerator.cpp ../WaveformGenerator/WaveformGenerator.hh
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@


clean:
	rm -f ./build/SPEFunctions.o
	rm -f ./build/PhotonFits.o
	rm -f ./build/PhotonAreaDistributions.o
	rm -f ./build/WaveformGenerator.o
	rm -f ./build/*.o
	rm -f ./AreaCorrections

distclean:
	rm -f ./AreaCorrections
	rm -f ./build/AreaCorrections.o
	@echo Removing all object files from build directory.
	@rm -f $(OBJS)
