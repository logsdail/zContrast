CC=mpic++
CFLAGS=-c -O3 -Wall -Werror -pedantic -fopenmp -Wno-long-long 
LDFLAGS=
LIBRARIES=-lfftw3 -lm -fopenmp 
SOURCES=CCluster.cpp \
	CContrastImage.cpp \
	CFitting.cpp \
	CompareData.cpp \
	CompareWorkLink.cpp \
	CreateRotations.cpp \
        CVariables.cpp \
	CZContrast.cpp \
	GatherFiles.cpp \
        GlobalOptimisation.cpp \
        GlobalOptimisationMutation.cpp \
	GlobalOptimisationOffspring.cpp \
        GlobalOptimisationUtils.cpp\
	Information.cpp \
	LocalMinimisationDetails.cpp \
	Main.cpp \
        ReadVariables.cpp \
	Search.cpp \
        Utils.cpp \
        Work.cpp \

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=zContrast

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBRARIES)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ 

clean:
	rm -f ${OBJECTS} ${EXECUTABLE}
