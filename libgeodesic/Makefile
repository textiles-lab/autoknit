TARGET = geodesic

CC = gcc
LD = gcc

INCLUDE = -I./include -I/usr/local/include -I/usr/include/suitesparse
BLAS_LIBS         = 
SUITESPARSE_LIBS  = -L/usr/local/lib -lspqr -lumfpack -lcholmod -lmetis -lcolamd -lccolamd -lcamd -lamd -lm -lsuitesparseconfig
#HSL_LIBS          = -lgfortran -lhsl_ma87 -lmetis -framework Accelerate
LIBS = $(BLAS_LIBS) $(SUITESPARSE_LIBS) $(HSL_LIBS)
CFLAGS = -O2 -Wall -Werror -pedantic $(INCLUDE) -DNDEBUG
LFLAGS = -O2 -Wall -Werror -pedantic                    

OBJS = obj/hmCholeskyFactor.o \
       obj/hmContext.o \
       obj/hmDenseMatrix.o \
       obj/hmPairSizeTDouble.o \
       obj/hmSparseMatrix.o \
       obj/hmTriDistance.o \
       obj/hmTriMesh.o \
       obj/hmVectorDouble.o \
       obj/hmVectorPairSizeTDouble.o \
       obj/hmVectorSizeT.o \
       obj/main.o

all: $(TARGET)

clean:
	rm -f $(OBJS)
	rm -f $(TARGET)
	rm -rf doc/latex doc/html

$(TARGET): $(OBJS)
	$(LD) $(LFLAGS) $(OBJS) $(LIBS) -o $(TARGET)

obj/hmCholeskyFactor.o: src/hmCholeskyFactor.c include/hmCholeskyFactor.h include/hmDenseMatrix.h include/hmConstants.h include/hmSparseMatrix.h include/hmUtility.h include/hmConstants.h
	$(CC) $(CFLAGS) -c src/hmCholeskyFactor.c -o obj/hmCholeskyFactor.o

obj/hmContext.o: src/hmContext.c include/hmContext.h include/hmConstants.h
	$(CC) $(CFLAGS) -c src/hmContext.c -o obj/hmContext.o

obj/hmDenseMatrix.o: src/hmDenseMatrix.c include/hmDenseMatrix.h include/hmConstants.h include/hmUtility.h
	$(CC) $(CFLAGS) -c src/hmDenseMatrix.c -o obj/hmDenseMatrix.o

obj/hmPairSizeTDouble.o: src/hmPairSizeTDouble.c include/hmPairSizeTDouble.h
	$(CC) $(CFLAGS) -c src/hmPairSizeTDouble.c -o obj/hmPairSizeTDouble.o

obj/hmSparseMatrix.o: src/hmSparseMatrix.c include/hmSparseMatrix.h include/hmConstants.h include/hmUtility.h
	$(CC) $(CFLAGS) -c src/hmSparseMatrix.c -o obj/hmSparseMatrix.o

obj/hmTriDistance.o: src/hmTriDistance.c include/hmTriDistance.h include/hmTriMesh.h include/hmVectorPairSizeTDouble.h include/hmPairSizeTDouble.h include/hmVectorDouble.h include/hmDenseMatrix.h include/hmConstants.h include/hmSparseMatrix.h include/hmCholeskyFactor.h include/hmVec3.h include/hmUtility.h include/hmConstants.h include/hmVectorSizeT.h
	$(CC) $(CFLAGS) -c src/hmTriDistance.c -o obj/hmTriDistance.o

obj/hmTriMesh.o: src/hmTriMesh.c include/hmTriMesh.h include/hmVectorPairSizeTDouble.h include/hmPairSizeTDouble.h include/hmVec3.h include/hmUtility.h
	$(CC) $(CFLAGS) -c src/hmTriMesh.c -o obj/hmTriMesh.o

obj/hmVectorDouble.o: src/hmVectorDouble.c include/hmVectorDouble.h include/hmUtility.h include/hmConstants.h
	$(CC) $(CFLAGS) -c src/hmVectorDouble.c -o obj/hmVectorDouble.o

obj/hmVectorPairSizeTDouble.o: src/hmVectorPairSizeTDouble.c include/hmVectorPairSizeTDouble.h include/hmPairSizeTDouble.h include/hmUtility.h include/hmConstants.h
	$(CC) $(CFLAGS) -c src/hmVectorPairSizeTDouble.c -o obj/hmVectorPairSizeTDouble.o

obj/hmVectorSizeT.o: src/hmVectorSizeT.c include/hmVectorSizeT.h include/hmUtility.h include/hmConstants.h
	$(CC) $(CFLAGS) -c src/hmVectorSizeT.c -o obj/hmVectorSizeT.o

obj/main.o: src/main.c include/hmTriDistance.h include/hmTriMesh.h include/hmVectorPairSizeTDouble.h include/hmPairSizeTDouble.h include/hmVectorDouble.h include/hmDenseMatrix.h include/hmConstants.h include/hmSparseMatrix.h include/hmCholeskyFactor.h include/hmContext.h include/hmUtility.h include/hmVectorSizeT.h
	$(CC) $(CFLAGS) -c src/main.c -o obj/main.o

