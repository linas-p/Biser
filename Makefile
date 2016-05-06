BUILD_DIR = build

all: test R_pac
test: $(BUILD_DIR)
	cd $(BUILD_DIR); cmake $(CURDIR); make

$(BUILD_DIR):
	mkdir -p $@


R_HOME := $(shell R RHOME)

## include headers and libraries for R
RCPPFLAGS := $(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS := $(shell $(R_HOME)/bin/R CMD config --ldflags)

## include headers and libraries for Rcpp interface classes
RCPPINCL := $(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := $(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)
INC_DIR = ./include
SRC_DIR = ./src
#JSON_INC =./external/json/src

%.o : %.cpp
	PKG_CPPFLAGS="-fPIC -shared -O2 $(RCPPFLAGS) $(RCPPINCL) -I$(INC_DIR)/" PKG_LIBS="$(RLDFLAGS) $(RCPPLIBS)" R CMD SHLIB $<
#PKG_CPPFLAGS="-std=c++11 -fPIC -shared $(RCPPFLAGS) $(RCPPINCL) -I$(JSON_INC) -I$(INC_DIR)/" PKG_LIBS="$(RLDFLAGS) $(RCPPLIBS)" R CMD SHLIB $<

cpp_sources := $(wildcard $(SRC_DIR)/*.cpp ./calculator_r.cpp)
cpp_sharedlibs := $(patsubst %.cpp,%.o,$(cpp_sources))

R_pac: shared_libs
	g++ -shared -o calculator.so $(wildcard $(SRC_DIR)/*.o) $(wildcard ./*.o) -L/usr/lib/R/lib -lR -L/usr/lib/R/lib -lR

shared_libs: $(cpp_sharedlibs)

clean:
	rm -rf $(BUILD_DIR) calculator calculator_r.o calculator_r.so a.out \
		$(SRC_DIR)/explicit_calculator.o $(SRC_DIR)/utils.o testing.o \
		$(SRC_DIR)/explicit_calculator.so $(SRC_DIR)/utils.o testing.so \
		calculator.so $(SRC_DIR)/utils.so ./R/calculator.so ./R/rsconnect/\
		./R/results.rda

.PHONY: all cmake clean test
