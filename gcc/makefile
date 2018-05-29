default: clean
	g++ -g -std=c++11 -o javatest javatest.cc hilbert.cpp -D HILBERT_NEXTA
global: clean
	g++ -g -std=c++11 -o javatest javatest.cc *.cpp -D GLOBAL
ieee: clean
	g++ -g -std=c++11 -o main  *.cpp -D TEST_IEEE
vector: clean
	g++ -g -std=c++11 -o main  *.cpp -D TEST_IEEE_VTX
point: clean
	g++ -g -std=c++11 -o main  *.cpp -D TEST_IEEE_PT
next: clean
	g++ -g -std=c++11 -o main  *.cpp -D TEST_NEXT
cmp: clean
	g++ -g -std=c++11 -o main  *.cpp -D TEST_CMP
ieeec: clean
	g++ -g -std=c++11 -o javatest javatest.cc hilbert.cpp -D TEST_IEEEC
ieeet: clean
	g++ -g -std=c++11 -o javatest javatest.cc hilbert.cpp -D IEEE_TEST
box: clean
	g++ -g -std=c++11 -o javatest javatest.cc hilbert.cpp -D HILBERT_BOX 
box2: clean
	g++ -g -std=c++11 -o javatest javatest.cc hilbert.cpp -D TEST_IEEE_VTX2 
boxp: clean
	g++ -g -std=c++11 -o javatest javatest.cc hilbert.cpp -D HILBERT_POINT
boxpe: clean
	g++ -g -std=c++11 -o javatest javatest.cc hilbert.cpp -D HILBERT_POINT_IEEE
boxpet: clean
	g++ -g -std=c++11 -o javatest javatest.cc hilbert.cpp -D HILBERT_POINT_TEMP
nextt: clean
	g++ -g -std=c++11 -o javatest javatest.cc hilbert.cpp -D HILBERT_NEXTA
clean:
	-rm javatest main
