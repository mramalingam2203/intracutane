using System;

class Solvers
{

// Tha analytical solution
static void First_law(out double [,] j,out  double [,] c, out double [] q, int nt, int nx,
        double h, double tmax, double k, double cv,double d) {

    j = new double[nt,nx];
    c = new double[nt,nx];
    q = new double[nt];

    double dx = h/nx;
    double dt = tmax/nt;
    double x,t,summ;
    
    //Console.WriteLine(String.Format("{0},{1},{2},{3},{4},{5},{6},{7}",
    //            nt, nx, h, tmax, k, cv, l, d));


    for (int kt =0; kt < nt; kt++){
          t = kt*dt;
        for (int kx=0; kx < nx; kx++){
            x = kx*dx;
            summ = 0.0;
        
            for (int ks=1; ks <10; ks++){
                summ = summ + (Math.Sin(ks*Math.PI*(kx+1)*x/h)
                * Math.Exp(-Math.Pow(kt*Math.PI/l,2.0)*d*t) /ks);
            };

            c[kt,kx] = k*cv*((1-x/h)- (2*summ/Math.PI));
            
            summ = 0.0;
        
            for (int ks=1; ks <10; ks++){                         
                summ  = summ + ((k*cv*d/h)*(1+(2*(Math.Pow(-1.0,ks)     
                            *Math.Exp(-Math.Pow(ks*Math.PI/h,2.0)*d*t))))); 
                };         
            
            j[kt,kx] = k*cv*(1+ 2*summ)/h;
            
        }

        summ = 0;
        for (int ks=1; ks <10; ks++){
        summ = summ +((Math.Pow(-1.0,ks)/Math.Pow(ks,2.0))
                *Math.Exp(-Math.Pow(ks*Math.PI/h,2.0)*d*t));
        }

        
        q[kt]= k*h*cv*((d*t/h*h)-0.1667-(2/Math.Pow(Math.PI,2.0)*summ));


    }
    
    for (int kt=0; kt < nt; kt+= 50){
    {
        for (int kx=0; kx < nx ; kx += 10)
            Console.Write("{0:F4} ",c[kt,kx]);
    }
       Console.WriteLine();
    }

    return;

}  

// Infinite dose - finite difference - explicit
static void IFD_explicit(out double  [,] c,out  double [] j, out double[] q, 
                int nx,int nt,double T, double X, double d, double cc)
{
    
    double dx = X/nx;
    double dt = T/nt;

    c = new double[nt,nx+1];
    j = new double[nt];
    q = new double[nt];

    double r = dt/(dx*dx);
    //Console.WriteLine(String.Format("{0},{1},{2},{3},{4},{5},{6}",nt, nx, X, T,d, cc,r*d));

    // I. C.
    j[0]=0.0;
    q[0]=0.0;
    
    for (int kx = 0; kx < nx ; kx++)
        c[0,kx]=0.0; 
    
    // B.C.
    for (int kt = 1; kt < nt; kt++){
        c[kt,0]= cc;
        c[kt,nx] = 0.0;
     }

    //infinite dose 
    for (int kt= 1; kt < nt; kt++)
        for (int kx=1; kx < nx ; kx++)
            c[kt,kx] = r*d*c[kt-1,kx-1] + (1-2*r*d)*c[kt-1,kx] + r*d*(c[kt-1,kx+1]);
        

    for (int kt=1; kt < nt; kt++){
        j[kt] = -d*(c[kt,nx]-c[kt,nx-1])/dx;
        q[kt] = q[kt-1] + j[kt]*dt ;
    }
    
//    for (int kt=1; kt < nt; kt++)
 //       for (int kx=1; kx < nx ; kx++)
 //           if(kt%100 == 0)         
 //               Console.WriteLine("{0:F4},{1:F4}",kt*dt, j[kt]);
/*
    for(int kt=1; kt < nt; kt++)
         Console.WriteLine("{0:F4}",j[kt]);
*/

}

// infinite dose - Implicit --- needs completion
static void IFD_CrankNicolson_Algorithm()                      
{                                                                              
																			   
	double[] V, L, U, ActualSolution, Z;                                       
	double FT, FX, alpha, H, K, lambda, T, X;                                  
	int N, M, M1, M2, I1, I, J;                                                
																			   
	V = new double[25]; L = new double[25];                                    
	U = new double[25]; Z = new double[25];                                    
	ActualSolution = new double[25];                                           
																			   
	FX = 1;     //These input parameters values                                
	FT = 0.5;   //are set by the user and may                                  
	alpha = 1;  //be varied accordingly                                        
	M = 10;                                                                    
	N = 50;                                                                    
																			   
	//Initialize variables with input parameter values                         
	M1 = M - 1;                                                                
	M2 = M - 2;
	H = FX / M;
	K = FT / N;
	lambda = alpha * alpha * K / (H * H);
	V[M - 1] = 0.0;
	
	//for (I = 1; I <= M1; I++) V[I - 1] = functionX(I * H);

	 /* Solve the tridiagonal linear system */
	 L[0] = 1.0 + lambda;
	 U[0] = -lambda / (2.0 * L[0]);
	
	 for (I = 2; I <= M2; I++)
	 {
		 L[I - 1] = 1.0 + lambda + lambda * U[I - 2] / 2.0;
		 U[I - 1] = -lambda / (2.0 * L[I - 1]);
	 }
	 L[M1 - 1] = 1.0 + lambda + 0.5 * lambda * U[M2 - 1];
	
	 for (J = 1; J <= N; J++)
	 {
		 T = J * K;
		 Z[0] = ((1.0 - lambda) * V[0] + lambda * V[1] / 2.0) / L[0];
		 for (I = 2; I <= M1; I++)
			 Z[I - 1] = ((1.0 - lambda) * V[I - 1] + 0.5 * lambda * (V[I] + V[I - 2] + Z[I - 2])) / L[I - 1];
		 V[M1 - 1] = Z[M1 - 1];
		 for (I1 = 1; I1 <= M2; I1++)
			 {
				 I = M2 - I1 + 1;
				 V[I - 1] = Z[I - 1] - U[I - 1] * V[I];
			 }
		 }
}


// Finite dose - seems to work fine -- but needs full check
static void FD_explicit(out double  [,] c,out  double [] j_0, out double[] q_0,
                  	out  double [] j_X, out double[] q_X,
                    int nx,int nt,double T, double X, double d,
					 double k, double cv,double A, double V)
{

	double dx = X/nx;
    double dt = T/nt;
    
     c = new double[nt,nx+1];
     j_0 = new double[nt];
     q_0 = new double[nt];
     j_X = new double[nt];
     q_X = new double[nt];
    


     double r = dt/(dx*dx);
     //Console.WriteLine(String.Format("{0},{1},{2},{3},{4},{5},{6}",nt, nx, X, T,d, cc,r*d));
    
     // I. C.
     j_0[0]=0.0;
     q_0[0]=0.0;

     j_X[0]=0.0;
     q_X[0]=0.0;
     

	for (int kx = 0; kx < nx ; kx++)
     	c[0,kx]=0.0;

    for (int kt=1; kt < nt; kt++){
        c[kt,0] = k*(cv- A*q_0[kt-1]/V);
        
        j_0[kt] =  -(d/dx)*(c[kt-1,1]-c[kt-1,0]);
        q_0[kt] = q_0[kt-1] + j_0[kt]*dt;

        { for(int kx=1; kx < nx;kx++)
            c[kt,kx] = r*d*c[kt-1,kx-1] + (1-2*r*d)*c[kt-1,kx] + r*d*(c[kt-1,kx+1]);
			}
        
        j_X[kt] =  -(d/dx)*(c[kt-1,nx-1]-c[kt-1,nx-2]);
        q_X[kt] = q_X[kt-1] + j_X[kt]*dt;		
		
    }
/*
	 for (int kt=1; kt < nt; kt++)
          for (int kx=1; kx < nx ; kx++)
              // if(kt%100 == 0)
                    Console.WriteLine("{0:F4},{1:F4},{2:F4}",kt*dt, kx*dx, c[kt,kx]);
*/

}


static void FD_implicit(double X, 
        double T, int nx,int nt, double d, double cc)
{

    double dx = X/nx;
    double dt = T/nt;
    double r = dt/(dx*dx);    
    
    double[] diag = new double[nx]; 
	double[] superdiag = new double[nx-1];
	double[] subdiag = new double[nx-1];
	double[] sol = new double[nx];
	double[] rhs = new double[nx];
	//r = 2.0;    
	//d = 0.01;

    for (int i=0; i < nx; i++){
        diag[i] = 1-2*r*d;
		sol[i] = 0.0;
		if (i==0) rhs[i] = cc;
        else rhs[i]= 0.0;
		}
    
    for (int i=0;i <  nx-1; i++){
        superdiag[i] = r;
        subdiag[i] = r;
        }

	
	tdma(nx,diag, subdiag, superdiag,rhs, out sol);
	
}


static void tdma(int n, double[] a, double[] b,double[] c,double[] d, out double[] f) {

	double[] c_star = new double[n-1];
	double[] d_star = new double[n-1];
	f = new double[n];

	c_star[0]=c[0]/b[0];
	d_star[0]=d[0]/b[0];

	for (int i=1; i<n-1; i++) {
    	double m = 1.0 / (b[i] - a[i] * c_star[i-1]);
    	c_star[i] = c[i] * m;
    	d_star[i] = (d[i] - a[i] * d_star[i-1]) * m;
//		Console.Write("{0:F4} {1:F4}",c_star[i],d_star[i]);
 	 }

	for (int i=n-1; i-- > 0; ) {
    	f[i] = d_star[i] - c_star[i] * d[i+1];
 	 }
	
	for (int i=0; i < n; i++)
		Console.Write("{0:F4}",f[i]);




}


public static void Main()
{
	
    double K = 2.1410; // partition coefficient
    double Cv = 1500; // chemical concentration in vehicle
    double L = 0.05; // skin thickness
    double D = 8.3704e-4; // diffustion coefficient
    double A = 1.77; // surface area of skin application
    double V = 1.0; // volume of vehicle
    int Ntime = 50000; // no. of time steps-- should be taken care of when using explicit method because the solution may diverge
    int Nx = 100; // no. of spatial elemeng
    
    double H = 0.05; // skin thickness -numerical
    double T = 6.0; // maximum time of computation
    double C0 = K*Cv;

    double[,] J = new double[Ntime,Nx]; // flux - analytical solution
    double[,] C = new double[Ntime,Nx]; // concentration
    double[] Q = new double[Ntime]; // concentratio
    
    double [] J1_0 = new double[Ntime]; // flux - numerical solutions
    double[] Q_0 = new double[Ntime]; 

    double [] J1_X = new double[Ntime]; // flux - numerical solutions
    double[] Q_X = new double[Ntime]; 

//    First_law(out J,out C,out Q,Ntime,Nx,H,T,K,Cv,L,D); // the analytical solver

//  IFD_explicit(out C,out J1, out Q, Nx,Ntime ,T, H, D, C0); // infinite dose - finite difference - explicit - solver 
//	FD_explicit(out C,out  J1_0, out Q_0,out J1_X, out Q_X, Nx,Ntime,T, H, D, K,Cv,A,V); // finite dose - finite difference - explicit - solver
    FD_implicit(D,H,Nx);

//    FD_implicit(H, T, Nx,Ntime, D,C0);
    double dx = H/Nx;
    double dt = T/Ntime;


    //for (int kx=0; kx < Nx+1 ; kx=kx+2)
    //   Console.Write("{0:F4} ",kx*dx);
    // Console.WriteLine();
    Console.WriteLine("TIME,         J,      Q,     J_in,    Q_in,    <------------------C along depth X--------------------------------------->");

    for (int kt=0; kt < Ntime; kt+= 50){
        Console.Write("{0:F7} {1:F4} {2:F4} {3:F4} {4:F4} ",kt*dt,J1_X[kt],Q_X[kt], J1_0[kt],Q_0[kt]);
        {
        for (int kx=0; kx < Nx ; kx += 10)
            Console.Write("{0:F4} ",C[kt,kx]);
       }
       Console.WriteLine();
    }


    }
}
    