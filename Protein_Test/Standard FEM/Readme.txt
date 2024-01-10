This program, Standard FEM-PB, is developed based on the toolbox Parallel Hierarchical Grid (PHG), whch is currently under active development at State Key Laboratory of Scientific and Engineering Computing of Chinese Academy of Sciences. Please refer to http://lsec.cc.ac.cn/phg/ to get more details about the installation and usage of PHG. Afer PHG has been successfully installed, you need to change the content of 'PHG_MAKEFILE_INC' in Makefile to the installation path of PHG on your computer. Then you can compile and run our program on a parallel computer by running 'Code/work.sh'.

You can modify 'work.sh' to adjust some parameters of our program:
-np 4 # the number of cores
-fn_pqr ../PQR/ADP-center.pqr # PQR file(in folder ../PQR)
-fn_mesh ../Body-fitted-Mesh/ADP.mesh # Body-fitted mesh file(in folder ../Body-fitted-Mesh)
-C_b #ionic concentration unit: mol/L
