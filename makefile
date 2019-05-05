

all:
	g++ -o quantum -O3 quantum.cpp ImageUtil.cpp

run: all
	-mkdir out
	./quantum
	ffmpeg  -i out/img%03d.bmp -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4

clean: 
	-rm -r out
	-rm quantum
