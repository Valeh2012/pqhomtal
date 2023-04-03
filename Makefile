#This file is part of FV-NFLlib

#Copyright(C) 2016 CryptoExperts

#This program is free software; you can redistribute it and / or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software Foundation,
#Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110 - 1301 USA

#C++ Compiler
ifndef FV_CXX
ifdef CXX
FV_CXX  = $(CXX)
else 
FV_CXX 	= g++
endif
endif
CXXFLAGS    = -std=c++17 -funroll-loops -Ofast -Wall -g

#Echo function
ECHO        = echo -e

#Library and include flags
LFLAGS = 
IFLAGS = 

NFL_INCLUDE_PATH = # define nfllib include path

ifdef NFL_LIBRARY_PATH
LFLAGS      += -L$(NFL_LIBRARY_PATH)
endif

ifdef NFL_INCLUDE_PATH
IFLAGS      += -I$(NFL_INCLUDE_PATH)
endif

#Update the flags
LFLAGS      += -lnfllib -lmpfr -lgmpxx -lgmp 
IFLAGS      += -I. 
IFLAGS      += -I./lib 

#Targets
TARGETS     = Test_main Test_vericrypt


lib_fips.o:
	gcc ./lib/fips202.c -c -o lib_fips.o
	
util.o: lib_fips.o 
	$(FV_CXX) $(CXXFLAGS) $(IFLAGS) params.hpp util.hpp lib_fips.o util.cpp -c  util.o

all: $(TARGETS)

Test_main: lib_fips.o main.cpp
	$(FV_CXX) $(CXXFLAGS) $(IFLAGS) $^ -o $@ $(LFLAGS) 

Test_vericrypt: lib_fips.o vericrypt.cpp
	$(FV_CXX) $(CXXFLAGS) $(IFLAGS) $^ -o $@ $(LFLAGS)


tests: $(TARGETS)
	@$(ECHO) "\n-------------------------\n"
	@$(ECHO) "Main.cpp"
	@$(ECHO) "\n-------------------------\n"
	@./Test_main

	@$(ECHO) "\n-------------------------\n"
	@$(ECHO) "Test vericrypt"
	@$(ECHO) "\n-------------------------\n"
	@./Test_vericrypt



clean:
	-rm -rf $(TARGETS) *.dSYM
	-rm -rf *.o
	-rm -rf *.gch
