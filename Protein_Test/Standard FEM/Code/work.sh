make

mpirun -np 36 ./PBE -fn_mesh ../Body-fitted-Mesh/ADP.mesh -fn_pqr ../PQR/ADP-center.pqr -C_b 0.1 -options_file advanced.options
