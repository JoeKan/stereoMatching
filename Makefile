.PHONY: all

all: stereo

stereo:
	g++ -std=c++14 -w stereo.cpp Eigen.h FreeImageHelper.h FreeImageHelper.cpp -I ./libs/Eigen/ -o stereo -lfreeimage

clean:
	rm -rf stereo
