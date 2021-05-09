# Makefile for the correlator

CC = gcc

CMUFLAGS = -mfpmath=sse -fPIC -DVERBOSE

WFLAGS= -Wmissing-prototypes -Wall -Wextra #-Wstrict-prototypes

CFLAGS = -std=gnu11 -march=native -D_GNU_SOURCE $(WFLAGS) -D_SHARED -D_UNROLL -m64 $(CMUFLAGS)

DBGFLAGS =

LIBRARIES = -lm -lz -lpthread

INC_DIRS = include/

SOURCE = src

EXECNAME = qLD-compute

BLIS_LIB = ../../blis/blis/lib/libblis.a
BLIS_INC = ../../blis/blis/include/blis

CL_LIB = -lOpenCL
CL_INC = # In case you need to include CL directory, do it here

SRC = $(filter-out $(SOURCE)/*_gpu.c, $(wildcard $(SOURCE)/*.c))

OBJ = $(patsubst $(SOURCE)/%.c, build/%.o, $(SRC))

cblas_cpu_syrk: LDFLAGS += -DSYRK
cblas_cpu_syrk: CFLAGS += -DSYRK
cblas_cpu_syrk: cblas_cpu

cblas_cpu: INC_DIRS += $(BLIS_INC)
cblas_cpu: OBJ += $(BLIS_LIB)
cblas_cpu: LDFLAGS += -DCBLAS_USE -msse4.2 -fopenmp
cblas_cpu: CFLAGS += -DCBLAS_USE
cblas_cpu: cpu

cblas_gpu_syrk: LDFLAGS += -DSYRK
cblas_gpu_syrk: CFLAGS += -DSYRK
cblas_gpu_syrk: cblas_gpu

cblas_gpu: INC_DIRS += $(BLIS_INC)
cblas_cpu: OBJ += $(BLIS_LIB)
cblas_gpu: LDFLAGS += -DCBLAS_USE -msse4.2 -fopenmp
cblas_gpu: CFLAGS += -DCBLAS_USE
cblas_gpu: gpu

debug_cblas_cpu: INC_DIRS += $(BLIS_INC)
cblas_cpu: OBJ += $(BLIS_LIB)
debug_cblas_cpu: LDFLAGS += -DCBLAS_USE -msse4.2 -fopenmp
debug_cblas_cpu: CFLAGS += -DCBLAS_USE
debug_cblas_cpu: debug_cpu

debug_cblas_gpu: INC_DIRS += $(BLIS_INC)
cblas_cpu: OBJ += $(BLIS_LIB)
debug_cblas_gpu: LDFLAGS += -DCBLAS_USE -msse4.2 -fopenmp
debug_cblas_gpu: CFLAGS += -DCBLAS_USE
debug_cblas_gpu: debug_gpu

cstates: CFLAGS += -DCOUNTSTATES
cstates: cpu

cpu: CFLAGS += -O3 -DSORTED #-DBENCHMARK
cpu: $(EXECNAME)

gpu: CFLAGS += -DGPU -O3 #-DBENCHMARK
gpu: LIBRARIES += $(CL_LIB)
gpu: INC_DIRS += $(CL_INC)
gpu: $(EXECNAME)

debug_cstates: CFLAGS += -DCOUNTSTATES
debug_cstates: debug_cpu

debug_cpu: CFLAGS += -DVERBOSE
debug_cpu: debug

debug_gpu: CFLAGS += -DGPU -DVERBOSE
debug_gpu: LIBRARIES += $(CL_LIB)
debug_gpu: INC_DIRS += $(CL_INC) 
debug_gpu: debug

debug: DBGFLAGS += -g3
debug: CFLAGS += -O0
debug: $(EXECNAME)

${EXECNAME}: $(OBJ)
	$(CC) $(DBGFLAGS) -o $@ $^ $(LDFLAGS) $(INC_DIRS:%=-I%) $(LIBRARIES)

build/%.o: $(SOURCE)/%.c $()
	$(CC) $(DBGFLAGS) -o $@ -c $<  $(CFLAGS) $(INC_DIRS:%=-I%) $(LIBRARIES)

.PHONY: cblas_cpu_syrk cblas_gpu_syrk cblas_cpu cblas_gpu debug_cblas_cpu \
		debug_cblas_gpu cstates cpu gpu debug_cstates debug_cpu debug_gpu \
		clean mrproper

clean:
	@rm -rf build/*

mrproper: clean
	@rm -rf $(EXECNAME)*

