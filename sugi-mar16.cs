using System;

class Solvers
{



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
        double T, int nx,int nt, double d, double k, double cv)
{

    double dx = X/nx;
    double dt = T/nt;
    double r = d*dt/(dx*dx);    
    
    double[] diag = new double[nx]; 
	double[] superdiag = new double[nx];
	double[] subdiag = new double[nx];
	double[] sol = new double[nx+2];
	double[] rhs = new double[nx+2];

    double[] j_X = new double[nt];
    double[] q_X = new double[nt];

    for (int i=0; i < nx+2;i++){
        if(i<nx){
        diag[i] = -1-2*r;
        superdiag[i] = r; 
        subdiag[i] = r;    }
        if (i == 0) rhs[i] = -k*cv;
        else if (i == nx-1) rhs[i] = 0;
    }

    for (int j = 1; j< nt; j++){
        tdma(nx,subdiag, diag, superdiag,rhs, out sol);
        for (int i=0; i < nx+2;i++){
            rhs[i] = -sol[i];
            Console.Write("{0:F4} ",sol[i]);
        }
        Console.WriteLine("------------------------");

     //   j_X[j] =  -(d/dx)*(c[kt-1,nx-1]-c[kt-1,nx-2]);
     //   q_X[j] = q_X[kt-1] + j_X[kt]*dt;
    }

    

    
  /*  j_X[kt] =  -(d/dx)*(c[kt-1,nx-1]-c[kt-1,nx-2]);
    q_X[kt] = q_X[kt-1] + j_X[kt]*dt;   */    
    	
}


static void tdma(int n, double[] p, double[] q,double[] r,double[] s, out double[] f) {


	double[] rcap = new double[n];
	double[] scap = new double[n];
	f = new double[n+2];

    for(int i=0; i<n; i++)
        f[i] = 0.0;


	rcap[0]=r[0]/q[0];
	scap[0]=s[0]/q[0];
    //Console.Write("{0:F4} {1:F4}\n",rcap[0],scap[0]);

    for (int i=1; i < n ; i++){
        rcap[i] = r[i]/(q[i]-p[i]*rcap[i-1]);
        scap[i] = (s[i]-p[i]*scap[i-1])/(q[i]-p[i]*rcap[i-1]);
        //Console.Write("{0:F4} {1:F4} \n", rcap[i],scap[i]);
    }
/*
   for (int i=n-1; i >= 0 ; i--)
        Console.Write("{0:F4} {1:F4} \n", rcap[i],scap[i]);
*/

    f[0] = -s[0];
    f[n] = scap[n-1];
    f[n+1] = -s[n-1];
    
    for (int i = n-1; i > 0; i--){
        f[i] = scap[i-1] - rcap[i]*f[i+1] ; 
       // Console.Write("{0:F4} {1:F4} {2:F4}\n", scap[i-1],rcap[i-1],f[i]);
    }

   /* for(int i=0; i <= n+1; i++)
        Console.Write("{0:F4}\n",f[i]);
*/

}


public static void Main()
{
	
    double K = 2.141; // partition coefficient - try a max of 5e+4
    double Cv = 1500; // chemical concentration in vehicle - try a max of 5e6
    //double L = 5.0; // skin thickness - try a max of 5.0
    double D = 8.3e-4; // diffustion coefficient- try a max of 1e-3
    double A = 1.77; // surface area of skin application
    double V = 1.0;    

    int Ntime = 200; // no. of time steps-- should be taken care of when using explicit method because the solution may diverge
    int Nx = 25; // no. of spatial elemeng
    
     
    double H = 0.05; // skin thickness -numerical
    double T = 1.0; // maximum time of computation
    double C0 = K*Cv;

    double[,] J = new double[Ntime,Nx]; // flux - analytical solution
    double[,] C = new double[Ntime,Nx]; // concentration
    double[] Q = new double[Ntime]; // concentratio
    
    double [] J1_0 = new double[Ntime]; // flux - numerical solutions
    double[] Q_0 = new double[Ntime]; 

    double [] J1_X = new double[Ntime]; // flux - numerical solutions
    double[] Q_X = new double[Ntime]; 



//  IFD_explicit(out C,out J1, out Q, Nx,Ntime ,T, H, D, C0); // infinite dose - finite difference - explicit - solver 
	//FD_explicit(out C,out  J1_0, out Q_0,out J1_X, out Q_X, Nx,Ntime,T, H, D, K,Cv,A,V); // finite dose - finite difference - explicit - solver
    //D_implicit(double X,  double T, int nx,int nt, double d, double cc)

    FD_implicit(H,T,Nx,Ntime,D,K,Cv);

    //    FD_implicit(H, T, Nx,Ntime, D,C0);
    double dx = H/Nx;
    double dt = T/Ntime;

/*
    //for (int kx=0; kx < Nx+1 ; kx=kx+2)
    //   Console.Write("{0:F4} ",kx*dx);
    // Console.WriteLine();
    Console.WriteLine("TIME,         J,      Q,     J_in,    Q_in,    <------------------C along depth X--------------------------------------->");

    for (int kt=0; kt < Ntime; kt+= 500){
        Console.Write("{0:F7} {1:F4} {2:F4} {3:F4} {4:F4} ",kt*dt,J1_X[kt],Q_X[kt], J1_0[kt],Q_0[kt]);
        {
        for (int kx=0; kx < Nx ; kx += 10)
            Console.Write("{0:F4} ",C[kt,kx]);
       }
       Console.WriteLine();
    }
*/

    }
}
    
