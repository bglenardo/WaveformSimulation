#
# Makefile for many Unix compilers using the
# "standard" command name CC
#
CC=g++
DEBUG=-g
CFLAGS=-Wall -c $(DEBUG) $(BOOST)
LFLAGS=-Wall $(DEBUG)
RCFLAGS=`root-config --cflags --libs`

OBJS += ./build/WaveformGenerator.o
OBJS += ./build/SPEFunctions.o
OBJS += ./build/WaveGeneratorTest.o

all: ./WaveGeneratorTest

./WaveGeneratorTest: $(OBJS) ./build/WaveGeneratorTest.o
	$(CC) $(OSFLAGS) $(RCFLAGS) $(LFLAGS) $(RLFLAGS) $^ $(LIBS) -lMinuit -o $@

./build/WaveGeneratorTest.o: WaveGeneratorTest.cpp $(OBJS)
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@

./build/WaveformGenerator.o: ../WaveformGenerator/WaveformGenerator.cpp ../WaveformGenerator/WaveformGenerator.hh $(OBJS)
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@

./build/SPEFunctions.o: ../WaveformGenerator/SPEFunctions.cpp ../WaveformGenerator/SPEFunctions.hh
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@

clean:
	rm -f ./build/WaveformGenerator.o
	rm -f ./build/*.o
	rm -f ./WaveGeneratorTest

distclean:
	rm -f ./WaveGeneratorTest
	rm -f ./build/WaveGeneratorTest.o
	@echo Removing all object files from build directory.
	@rm -f $(OBJS)
