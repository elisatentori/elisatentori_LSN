
CC=g++
CXXFLAGS= -Wall -I.
TARGET=executable
SRCLIST=$(wildcard *.cpp)
OBJLIST=$(patsubst %.cpp,%.o,$(SRCLIST))

# $@ = elementi a sinistra dei :
# $^ = elementi a destra dei :
# $< = primo elemento nella lista delle dipendenze

%.o : %.cpp

ifndef VERBOSE
	$(CC) $(CXXFLAGS) -c -o $@ $<
	
else
	@echo 'Compiling $< ...'
	@$(CC) $(CXXFLAGS) -c -o $@ $<

endif

$(TARGET) : $(OBJLIST)

ifndef VERBOSE
	$(CC) -o $(TARGET) $^
	
else 
	@echo 'Linking ....'
	@$(CC) -o $(TARGET) $^

endif

dep:
	@echo 'Making dependencies'
	$(CC) -MM *.cpp > dependencies

clean:
	rm -rf *.o

all: dep $(TARGET)

-include dependencies