using System;

class Solvers
{



// Finite dose - seems to work fine -- but needs full check
static void FD_explicit(out double  [,] c,out  double [] j_0, out double[] q_0,
                  	out  double [] j_X, out double[] q_X,
                    int nx,int nt,double T, double X, double d,
					 double k, double cv,double A, double V, double a)
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
        c[kt,nx] = c[kt,nx-1]/a;	
		
    }

    /*
	 for (int kt=1; kt < nt; kt++)
          for (int kx=1; kx < nx ; kx++)
              // if(kt%100 == 0)
                    Console.WriteLine("{0:F4},{1:F4},{2:F4}",kt*dt, kx*dx, c[kt,kx]);
    */

}


public static void Main()
{

    // double K = 2.14104; // partition coefficient - try a max of 5e+4
    // double Cv = 1500; // chemical concentration in vehicle - try a max of 5e6
    // double L = 5.0; // skin thickness - try a max of 5.0
    // double D = 8.3704e-4; // diffustion coefficient- try a max of 1e-3
    // double A = 1.77; // surface area of skin application
    // double V = 1.0;    

    int Ntime = 50000; // no. of time steps-- should be taken care of when using explicit method because the solution may diverge
    int Nx = 100; // no. of spatial elemeng
        
     
    // double H = 0.05; // skin thickness -numerical
    // double T = 6.0; // maximum time of computation
    double K = 0 ;
    double Cv =0 ;
    double D = 0 ;
    double A = 0 ;
    double V = 0 ;
    double H = 0;
    double T = 0;
    double a1 = 0;

    inputs(out K, out Cv, out D, out A, out V, out H, out T,out a1);
  
    
    Console.Write("{0:F4} {0:F4} {0:F4} {0:F4} {0:F4} {0:F4} {0:F4} {0:F4}",K, Cv,D, A, V, H, T,a1);

    double C0 = K*Cv;

    double[,] J = new double[Ntime,Nx]; // flux - analytical solution
    double[,] C = new double[Ntime,Nx]; // concentration
    double[] Q = new double[Ntime]; // concentratio
    
    double [] J1_0 = new double[Ntime]; // flux - numerical solutions
    double[] Q_0 = new double[Ntime]; 

    double [] J1_X = new double[Ntime]; // flux - numerical solutions
    double[] Q_X = new double[Ntime]; 

	
    FD_explicit(out C,out  J1_0, out Q_0,out J1_X, out Q_X, Nx,Ntime,T, H, D, K,Cv,A,V, a1); 
    // finite dose - finite difference - explicit - solver
    
  
    double dx = H/Nx;
    double dt = T/Ntime;


    //for (int kx=0; kx < Nx+1 ; kx=kx+2)
    //   Console.Write("{0:F4} ",kx*dx);
    // Console.WriteLine();
    Console.WriteLine("# TIME,         J,      Q,     J_in,    Q_in,    <------------------C along depth X--------------------------------------->");

    for (int kt=0; kt < Ntime; kt+= 1000){
        Console.Write("{0:F7} {1:F4} {2:F4} {3:F4} {4:F4} ",kt*dt,J1_X[kt],Q_X[kt], J1_0[kt],Q_0[kt]);
        {
        for (int kx=0; kx < Nx ; kx += 5)
            Console.Write("{0:F4} ",C[kt,kx]);
       }
       Console.WriteLine();
    }

   

    }

 
 static  void inputs(out double k, out double cv, out double d, out double a, 
    out double v, out double h, out double t, out double a0)
 {
    
    Console.Write("Enter partition coefficient, K, ( a maximum of 50000) : ");
    k = Convert.ToDouble(Console.ReadLine());
     if (k <= 0 || k > 5e4){ 
        Console.Write("K value is out of range. Please provide input as per directions...");
        System.Environment.Exit(0);
 
    }

    Console.Write("Enter chemical concentration in vehiele, C_v, ( a maximum of 5000000) : ");
    cv = Convert.ToDouble(Console.ReadLine());
    if (cv <= 0 || cv > 5e6){ 
        Console.Write("Cv value is outside the calculable range. Please provide input as per directions...");
        System.Environment.Exit(0);
    }

    Console.Write("Enter diffusion coefficient, D, ( a maximum of 0.0003) : ");
    d = Convert.ToDouble(Console.ReadLine());
    if (d <= 0 || d > 3e-4){ 
        Console.Write("D valus is outside the calculable range. Please provide input as per directions...");
        System.Environment.Exit(0);
    }

    Console.Write("Enter skin surface area of application, A: ");
    a = Convert.ToDouble(Console.ReadLine());
    if (a <= 0 || a > 5){ 
        Console.Write("A value is outside the calculable range. Please provide input as per directions...");
        System.Environment.Exit(0);
    }
    Console.Write("Enter volume, V: ");
    v = Convert.ToDouble(Console.ReadLine());
    if (v <= 0 || v >= 5){ 
        Console.Write("V value is outside the calculable range. Please provide input as per directions...");
        System.Environment.Exit(0);
    }

    Console.Write("Enter skin thickness ( a maximum of 5.0): ");
    h = Convert.ToDouble(Console.ReadLine());
    if (h <= 0 || h > 5){ 
        Console.Write("T value is outside the calculable range. Please provide input as per directions...");
        System.Environment.Exit(0);
    }
    
    Console.Write("Enter maximum time of computation: ");
    t = Convert.ToDouble(Console.ReadLine());
    if (t <= 0 || t > 10){ 
        Console.Write("T value is outside the calculable range. Please provide input as per directions...");
        System.Environment.Exit(0);
    }

    Console.Write("Enter far bc weighting factor: ");
    a0 = Convert.ToDouble(Console.ReadLine());
    if (a0 <= 0 || a0 > 1){ 
        Console.Write("a1 value is outside the calculable range. Please provide input as per directions...");
        System.Environment.Exit(0);
    }


   
 


    


        


 }


}
    
