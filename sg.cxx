#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);

void step(cmplx* psi0, cmplx* psi1, const double dt, const double dx,
	  const double lambda, const double omega, const double alpha,
	  const int Nx, const double xmin);

void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);
//-----------------------------------
int main(){
  const int Nx = 300; // Nach Aufgabenstellung
  const double xmin = -40; // Nach Aufgabenstellung
  const double xmax =  40; // Nach Aufgabenstellung
  const double Tend = 10*M_PI; // Nach Aufgabenstellung
  const double dx = (xmax-xmin)/Nx;
  const double dt = 0.01 * dx*dx; // in Anlehnung an Stabilitaetsbedigung v. Neumann
  double t = 0;
  const int Na = 10;
  int Nk = int(Tend / Na / dt + 0.5);

  const double lambda = 10;
  const double omega = 0.2; // => k = omega^2 * m = 0.2*0.2 = 0.04
  const double alpha = 0.4472136; // (0.04)^0.25 aus alpha = (m*k/hbar)^0.25

  stringstream strm;

  cmplx* psi0 = new cmplx[Nx];
  cmplx* psi1 = new cmplx[Nx];
  cmplx* haa;

  init(psi0, alpha, lambda, dx, dt, Nx, xmin);

  writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);

  for (int i = 1; i <= Na; i++) {
    for (int j = 1; j <= Nk-1; j++) {
      step(psi0,psi1,dt,dx,lambda,omega,alpha,Nx,xmin);
      
      haa = psi0;
      psi0 = psi1;
      psi1 = haa;
      
      t+=dt;
    }
    strm.str("");
    strm << "psi_" << i;
    writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
  }
  cout << "t = " << t << endl;

  return 0;
}
//-----------------------------------
void step(cmplx* f0, cmplx* f1, const double dt, const double dx,
	  const double lambda, const double omega, const double alpha,
	  const int Nx, const double xmin){
  cmplx* d = new cmplx[Nx];
  cmplx* a = new cmplx[Nx];
  double x;
  cmplx al = cmplx(0, -0.25 * dt / (dx*dx));
  
  for(int i = 0; i < Nx; i++){
    x = xmin + i * dx;
    d[i] = cmplx(1, 0.5 * dt /(dx*dx) + 0.25 * dt * x*x); // aus (1,(hbar*dt/2/m/dx/dx)+(dt/2/hbar)*V[i])
    //if((i%150) == 0){
    //  cout << "Re(d) = " << d[i].real() << "   Im(d) = " << d[i].imag() << endl; // debug
    //}
  }
  for(int i = 0; i < Nx; i++){
    a[i] = al; // aus (0, -hbar*dt/4/m/dx/dx)
  }

  a[0] = a[0] / d[0];
  f0[0] = f0[0] / d[0];
  d[0] = 1.0;
  for(int i=1;i<Nx;i++){
    // tridiagonale Matrix, mit obere und untere Nebendiagonalen a
    cmplx faktor = 1.0/(d[i] -  (-al)/d[i-1]*a[i-1]);
    d[i]  = 1;
    a[i]  = a[i] * faktor;
    f0[i] = ( f0[i] - (-al)/d[i-1]*f0[i-1] )  * faktor;
  }
  f1[Nx-1] = f0[Nx-1]/d[Nx-1]; // letzter Wert einzeln

  for (int i=Nx-2;i>=0;i--){
    f1[i] = f0[i]-a[i]*f1[i+1];
  }
  
  delete[] d; delete[] a;
}

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
  ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
  for(int i=0; i<Nx; i++){
    x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
    
    out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
  }
  out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	// const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
