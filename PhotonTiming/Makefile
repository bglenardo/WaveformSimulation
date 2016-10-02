#
# Makefile for many Unix compilers using the
# "standard" command name CC
#

include ../build/Makefile.config

OBJS  = ../build/LEvtUtilities.o
OBJS += ../build/LEvtChannel.o
OBJS += ../build/LEvtEvent.o
OBJS += ../build/LEvtFile.o
OBJS += ../build/LCvtChannel.o
OBJS += ../build/LCvtEvent.o
OBJS += ../build/LCvtFile.o
OBJS += ../build/RQFile_IO.o
OBJS += ../build/PhotonAreaDistributions.o
OBJS += ../build/PhotonFits.o
OBJS += ../build/SPEFunctions.o

all: ../bin/PhotonTiming

../bin/PhotonTiming: $(OBJS) ../build/PhotonTiming.o
	$(CC) $(OSFLAGS) $(RCFLAGS) $(LFLAGS) $(RLFLAGS) $^ $(LIBS) -lMinuit -o $@

../build/PhotonTiming.o: PhotonTiming.cpp
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@ 

../build/PhotonAreaDistributions.o: PhotonAreaDistributions_3.cpp PhotonAreaDistributions.h SPEparameters.h
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $(RLFLAGS) $< -o $@

../build/PhotonFits.o: PhotonFits_Brian_test.cpp PhotonFits.h ../build/PhotonAreaDistributions.o
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $(RLFLAGS) $< -o $@

../build/SPEFunctions.o: SPEFunctions_Rose_test.cpp SPEFunctions_Rose_test.h
	$(CC) $(OSFLAGS) $(CFLAGS) $(RCFLAGS) $< -o $@


include ../build/Makefile.objs

clean:
	rm -f ../bin/PhotonTiming
	rm -f ../build/PhotonTiming.o
	rm -f ../build/SPEFunctions.o
	rm -f ../build/PhotonFits.o
	rm -f ../build/PhotonAreaDistributions.o

distclean:
	rm -f ../bin/PhotonTiming
	rm -f ../build/PhotonTiming.o
	@echo Removing all object files from build directory.
	@rm -f $(OBJS)
