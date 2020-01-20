KOKKOS_PATH = ./kokkos-2.7.00
KOKKOS_SRC_PATH = ${KOKKOS_PATH}
# SRC = $(wildcard ${KOKKOS_SRC_PATH}/example/tutorial/07_Pi/*.cpp)
SRC = Pi_Kokkos.cpp
vpath %.cpp $(sort $(dir $(SRC)))

COMPUTING=INT32 INT64 FP32 FP64
MARSAGLIA=SHR3 CONG MWC KISS

# COMPUTING=FP32
# MARSAGLIA=MWC

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
CXXFLAGS = -O3
LINK = ${CXX}
LDFLAGS = 
EXE=cuda
else
CXX = g++
CXXFLAGS = -O3
LINK = ${CXX}
LDFLAGS =  
EXE=host
endif

DEPFLAGS = -M

OBJ = $(foreach TVAR,$(COMPUTING),$(foreach TRND,$(MARSAGLIA),$(notdir $(SRC:.cpp=))_$(TVAR)_$(TRND).o ) )
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

test: $(EXE)
	./$(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	@echo $(OBJ)
	$(foreach TVAR,$(COMPUTING),$(foreach TRND,$(MARSAGLIA), $(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(SRC:.cpp=)_$(TVAR)_$(TRND).o $(KOKKOS_LIBS) $(LIB)  -o $(SRC:.cpp=)_$(TVAR)_$(TRND).$(EXE) ; ) )

clean: kokkos-clean 
	rm -f *.o

mrproper: kokkos-clean 
	rm -f *.o *.cuda *.host

# Compilation rules

$(OBJ):$(SRC) $(KOKKOS_CPP_DEPENDS)
	$(foreach TVAR,$(COMPUTING),$(foreach TRND,$(MARSAGLIA),$(CXX) $(CFLAGS) -DT$(TVAR) -DT$(TRND) -DLONG -DTIME $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c -o $(<:.cpp=)_$(TVAR)_$(TRND).o $< ; ) )

# %.o:%.cpp $(KOKKOS_CPP_DEPENDS)
# 	@echo $(OBJ2) 
# 	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $< -o $(notdir $@)

