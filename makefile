CC = g++
C_OPTIMIZE_SWITCH = -O2 -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
LIBS = -lgsl -lgslcblas

CFLAGS = -Wall ${C_OPTIMIZE_SWITCH}  

run_secorder: run_secorder.o secorder_rec_1p.o calc_sqrtcov_rec_1p.o calc_rhos.o calc_stats_1p.o
	${CC} run_secorder.o secorder_rec_1p.o calc_sqrtcov_rec_1p.o calc_rhos.o calc_stats_1p.o -o run_secorder ${LIBS}

run_secorder_2p: run_secorder_2p.o secorder_rec_2p.o calc_sqrtcov_rec_2p.o calc_rhos.o calc_stats_2p.o
	${CC} run_secorder_2p.o secorder_rec_2p.o calc_sqrtcov_rec_2p.o calc_rhos.o calc_stats_2p.o -o run_secorder_2p ${LIBS}

run_secorder.o: secorder_rec_1p.hpp calc_stats_1p.hpp
run_secorder_2p.o: secorder_rec_2p.hpp calc_stats_2p.hpp
secorder_rec_1p.o: secorder_rec_1p.hpp calc_sqrtcov_rec_1p.hpp calc_rhos.hpp
secorder_rec_2p.o: secorder_rec_2p.hpp calc_sqrtcov_rec_2p.hpp calc_rhos.hpp calc_stats_2p.hpp
calc_sqrtcov_rec_1p.o: calc_sqrtcov_rec_1p.hpp
calc_sqrtcov_rec_2p.o: calc_sqrtcov_rec_2p.hpp
calc_rhos.o: calc_rhos.hpp
calc_stats_1p.o: calc_stats_1p.hpp
calc_stats_2p.o: calc_stats_2p.hpp


%.o : %.cpp
	${CC} -c ${CFLAGS} $<

clean:
	\rm *.o
