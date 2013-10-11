#################################################################
########################   MakeVars   ###########################
#################################################################

CC=g++ -fopenmp

CC_OPT=-I"./include" -Wall -Wextra -g -O0

LD_LIB=

LD_OPT=

MKDIR=mkdir -p ./obj

CIBLE=life_dcd

SRC=$(wildcard ./src/*.cpp)

OBJ=$(patsubst ./src/%.cpp,./obj/%.o,$(SRC))

#################################################################
########################   Makefile   ###########################
#################################################################

all:$(CIBLE)
	@echo "Compilation Success"

$(CIBLE):Makefile

./obj/%.o:./src/%.cpp
	@$(MKDIR)
	$(CC) $(CC_OPT) -c $< -o $@

$(CIBLE):$(OBJ)
	$(CC) $(CC_OPT) $(LD_LIB) $(OBJ) -o $@ $(LD_OPT)

clean:
	rm -f $(CIBLE) ./obj/*.o
