DIRS="OpenMP Pthreads MPI"

CURRENT=$PWD
ITERATIONS=1000000000

cd $CURRENT/OpenMP
for THREADS in $(seq 80) ; do for j in $(seq 10); do export OMP_NUM_THREADS=$THREADS ; echo -ne "$THREADS " ; /usr/bin/time ./Pi_OpenMP_LONG 100000000000 $THREADS >/dev/null ; done ; done > $CURRENT/PiOpenMP_$(hostname)_$(date "+%Y%m%d").log 2>&1
cd $CURRENT/Pthreads
for THREADS in $(seq 80) ; do for j in $(seq 10); do export OMP_NUM_THREADS=$THREADS ; echo -ne "$THREADS " ; /usr/bin/time ./Pi_Pthreads_LONG 100000000000 $THREADS >/dev/null ; done ; done > $CURRENT/PiPthreads_$(hostname)_$(date "+%Y%m%d").log 2>&1
cd $CURRENT/MPI
for THREADS in $(seq 80) ; do for j in $(seq 10); do export OMP_NUM_THREADS=$THREADS ; echo -ne "$THREADS " ; /usr/bin/time mpirun -np $THREADS -mca btl sm,self -x OMP_NUM_THREADS=1 ./Pi_MPI_LONG 100000000000 $THREADS >/dev/null ; done ; done > $CURRENT/PiMPI_$(hostname)_$(date "+%Y%m%d").log 2>&1
cd $CURRENT

