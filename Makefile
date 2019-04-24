proj.out: myProject.cpp BVHData.h BVHData.o Makefile
	g++ -std=c++11 -g -o proj.out myProject.cpp BVHData.o -ldart -ldart-gui -lboost_system -lassimp -framework OpenGL -framework GLUT -I /usr/local/include/eigen3
BVHData.o: Makefile BVHData.cpp BVHData.h
	g++ -std=c++11 -g -c BVHData.cpp -I /usr/local/include/eigen3
