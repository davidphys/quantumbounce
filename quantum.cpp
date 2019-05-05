#include <iostream>
#include <vector>
#include <cmath>
#include "ImageUtil.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;


//The space 0<x<L, 0<y<L is discretized into D*D cells.
const int D=512;
const double L=200;  
const double dx=L/D;
const double dt=dx*dx*0.4; //Stability condition is dt/dx^2<0.5 ?

//Store the real and imaginary parts of the wavefunction
vector<double> rev(D*D,0);
vector<double> imv(D*D,0);

//numerical code is in lines 112-155, the rest is just utility / plotting code.




//Some string manipulation functions for saving files. pad_int(1234,5) returns "01234". 
std::string pad_int(int arg, int padcount) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(padcount) << arg;
    return ss.str();
}

//Returns a file name in the form of "prefix00###suffix". For example "image0032.bmp"
std::string getFilename(std::string prefix, int num, int padcount, std::string suffix) {
    return prefix + pad_int(num, padcount) + suffix;
}

//From https://gist.github.com/fairlight1337/4935ae72bcbcc1ba5c72
void HSVtoRGB(float& fR, float& fG, float& fB, float& fH, float& fS, float& fV) {
  float fC = fV * fS; // Chroma
  float fHPrime = fmod(fH / 60.0, 6);
  float fX = fC * (1 - fabs(fmod(fHPrime, 2) - 1));
  float fM = fV - fC;
  
  if(0 <= fHPrime && fHPrime < 1) {
    fR = fC;
    fG = fX;
    fB = 0;
  } else if(1 <= fHPrime && fHPrime < 2) {
    fR = fX;
    fG = fC;
    fB = 0;
  } else if(2 <= fHPrime && fHPrime < 3) {
    fR = 0;
    fG = fC;
    fB = fX;
  } else if(3 <= fHPrime && fHPrime < 4) {
    fR = 0;
    fG = fX;
    fB = fC;
  } else if(4 <= fHPrime && fHPrime < 5) {
    fR = fX;
    fG = 0;
    fB = fC;
  } else if(5 <= fHPrime && fHPrime < 6) {
    fR = fC;
    fG = 0;
    fB = fX;
  } else {
    fR = 0;
    fG = 0;
    fB = 0;
  }
  
  fR += fM;
  fG += fM;
  fB += fM;
}

/* Plotting Code  */ 
void saveImage(int n){
  int imgs=D;
  DoubleImage ampx(imgs,imgs);
  DoubleImage ampy(imgs,imgs);
  DoubleImage ct(imgs,imgs);
  for(int i=0;i<D;i++){
    for(int j=0;j<D;j++){
      ampx.increase((i*imgs)/D,(j*imgs)/D,rev[D*i+j]);
      ampy.increase((i*imgs)/D,(j*imgs)/D,imv[D*i+j]);
      ct.increase((i*imgs)/D,(j*imgs)/D,1);
    }
  }
  Image img(imgs,imgs);
  for(int i=0;i<imgs;i++){
    for(int j=0;j<imgs;j++){
      if(ct.get(i,j)==0)
        continue; //possible only if imgs>D
      float x=float(ampx.get(i,j)/ct.get(i,j));
      float y=float(ampy.get(i,j)/ct.get(i,j));
      float amp=900*float(x*x+y*y);
      if(amp>1)
        amp=1;
      float phase=float(180.0/M_PI*(atan2(y,x)+M_PI));
      float saturation=0.9;
      float r,g,b;
      HSVtoRGB(r,g,b,phase,saturation,amp);
      img.put(i,j,floatToRGB(r,g,b));
    }
  }
  img.save(getFilename("out/img",n,3,".bmp"));
}

void zeroWavefunction(){
  rev= vector<double>(D*D,0);
  imv=vector<double>(D*D,0);
}
double re(int i, int j){
  if(i>=0 && i<D && j>=0 && j<D)
    return rev[D*i+j];
  return 0;
}
double im(int i, int j){
  if(i>=0 && i<D && j>=0 && j<D)
    return imv[D*i+j];
  return 0;

}
double V(double x,double y){
  return (y)*0.007;

}
void step(){
  for(int i=0;i<D;i++){
    for(int j=0;j<D;j++){
      rev[D*i+j]=rev[D*i+j]-0.5*dt*(im(i+1,j)+im(i-1,j)+im(i,j+1)+im(i,j-1)-4*imv[D*i+j])/(dx*dx)+dt*V(dx*i,dx*j)*im(i,j);
    }
  }
  for(int i=0;i<D;i++){
    for(int j=0;j<D;j++){
      imv[D*i+j]=imv[D*i+j]+0.5*dt*(re(i+1,j)+re(i-1,j)+re(i,j+1)+re(i,j-1)-4*rev[D*i+j])/(dx*dx)-dt*V(dx*i,dx*j)*re(i,j);
    }
  }
}
void addGaussian(double x0, double y0, double vx0, double vy0, double sigma){
  for(int i=0;i<D;i++){
    for(int j=0;j<D;j++){
      double x=i*dx,y=j*dx;
      double amplitude=exp(-((x-x0)*(x-x0)+(y-y0)*(y-y0))/(2*sigma*sigma))/sqrt(M_PI*sigma*sigma);
      double phase=vx0*x+vy0*y;
      rev[D*i+j]+=amplitude*cos(phase);
      imv[D*i+j]+=amplitude*sin(phase);
    }
  }
}


int main() {
  addGaussian(100,150,0.2,0.2,10);
  for(int n=0;n<300;n++){
    saveImage(n);
    for(int k=0;k<100;k++){
      step();
    }
  }
  return 0;
}

