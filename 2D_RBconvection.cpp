#include <cstdio>
//#include <math.h>
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
 
double pi = 4.0*atan(1.0);
const int Nn = 50;
const int Nz = 101;
const int Nt = 150000;
double Pr = 0.5;
double Ra = 1.0e006;
double a = 3.0;
double dt = 3.0e-6;
double dz = 1.0 / (Nz - 1);
double dx = 1.0 / Nn;

double psi[Nz][Nn+1];
double omg[Nz][Nn+1];
double tem[Nz][Nn+1];
double v[Nz][Nn+1];
double dtemdt[Nz][Nn+1][2], domgdt[Nz][Nn+1][2];
double sub = -1.0 / (dz * dz);
double sup = -1.0 / (dz * dz);
double dia[Nz], c[Nz], d[Nz], z[Nz];
int now = 1;
int old = 0;

int main()
{
	//定义三个文件输出流，输出tem,psi,omg到文件
	ofstream f1;
	ofstream f2;
	ofstream f3;

    //初始化 
    for (int k=0 ; k < Nz ; k++ )
    {
    	for (int n=0 ; n < Nn +1 ; n++)
    	{
    			dtemdt[k][n][old]=0;
    			domgdt[k][n][old]=0;
    			dtemdt[k][n][now]=0;
    			domgdt[k][n][now]=0;
    			omg[k][n]=0;
    			psi[k][n]=0;
    			tem[k][n]=0;
		}
	}
/*
    for (int i = 0; i < Nz; i++)
    {

        z[i] = i * dz;
    }*/
    //初始 
    for (int k = 0; k < Nz-1; k++)
    {
        tem[k][0] = 1.0 - k * dz;
        tem[k][1] = 0.01 * sin( pi *  k * dz );
        //psi[k][1] = 0.01 * sin( pi *  z[k] );
        //omg[k][1] = 0.01 * pi * pi * (1+1/(a*a)) * sin(pi*k*dz);
    }


    
    for (int t = 0; t < Nt; t++)
    {
        
        for (int k = 1; k < Nz - 1 ; k++)
        {   
            //线性项 
            dtemdt[k][0][now] = ( tem[k+1][0] - 2.0 * tem[k][0] + tem[k-1][0] ) / (dz * dz);

            for (int n = 1; n < Nn +1; n++)
            {
                dtemdt[k][n][now] = (tem[k+1][n] - 2.0*tem[k][n] +tem[k-1][n]) / (dz*dz) - (n*pi/a) * (n*pi/a) * tem[k][n];
                domgdt[k][n][now] = Ra * Pr * n * pi / a * tem[k][n] + Pr*((omg[k+1][n] - 2.0*omg[k][n] + omg[k-1][n]) / (dz*dz) - (n*pi/a) * (n*pi/a) *omg[k][n]);
            }
            
            //非线性项 
            for (int n1 = 1; n1 < Nn+1; n1++)
            {
                dtemdt[k][0][now] -= pi*n1/a/2.0*((psi[k+1][n1]-psi[k-1][n1])*tem[k][n1]/2.0/dz + psi[k][n1]*(tem[k+1][n1] - tem[k-1][n1])/2.0/dz);
            }
            
            for (int n = 1; n < Nn + 1; n++)
            {
                dtemdt[k][n][now] += n*pi/a*psi[k][n]*(tem[k+1][0]-tem[k-1][0])/2.0/dz;
                for (int n1 = 1; n1 < Nn+1; n1++)
                {
                    //double cc = pi/a/dz/4.0;
                    if ((1 <= (n-n1)) && ((n-n1) <= Nn))
                    {
                        int n2=n-n1;
                        dtemdt[k][n][now] -= pi/2.0/a*((-n1)/dz/2*(psi[k+1][n2]-psi[k-1][n2])*tem[k][n1]+n2/dz/2*psi[k][n2]*(tem[k+1][n1]-tem[k-1][n1]));
                        domgdt[k][n][now] -= pi/2.0/a*((-n1)*(psi[k+1][n2]-psi[k-1][n2])/dz/2.0*omg[k][n1]+n2*psi[k][n2]*(omg[k+1][n1]-omg[k-1][n1])/dz/2.0);
                    }
                    if ((1 <= (n+n1)) && ((n+n1) <= Nn))
                    {
                        int n2=n+n1;
                        dtemdt[k][n][now] -= pi/2.0/a*(n1*(psi[k+1][n2]-psi[k-1][n2])/dz/2.0*tem[k][n1]+n2*psi[k][n2]*(tem[k+1][n1]-tem[k-1][n1])/dz/2.0);
                        domgdt[k][n][now] += pi/2.0/a*(n1*(psi[k+1][n2]-psi[k-1][n2])/dz/2.0*omg[k][n1]+n2*psi[k][n2]*(omg[k+1][n1]-omg[k-1][n1])/dz/2.0);
                    }
                    if ((1 <= (n1-n)) && ((n1-n) <= Nn))
                    {
                        int n2=n1-n;
                        dtemdt[k][n][now] -= pi/2.0/a*(n1*(psi[k+1][n2]-psi[k-1][n2])/dz/2.0*tem[k][n1]+n2*psi[k][n2]*(tem[k+1][n1]-tem[k-1][n1])/dz/2.0);
                        domgdt[k][n][now] -= pi/2.0/a*(n1*(psi[k+1][n2]-psi[k-1][n2])/dz/2.0*omg[k][n1]+n2*psi[k][n2]*(omg[k+1][n1]-omg[k-1][n1])/dz/2.0);
                    }
                }
            } 
        }
       // cout<<psi[1][1]<<" 1 "<<psi[2][1]<<" 2 "<<psi[0][1]<<" 3 "<<tem[1][1]<<" 4 "<<tem[2][1]<<" 5 "<<tem[0][1]<<endl;
         //   dtemdt[1][2][now] -= (-1.0)/dz/2.0*(psi[2][1]-psi[0][1])*tem[1][1]+1.0/dz/2.0*psi[1][1]*(tem[2][1]-tem[0][1]);//pi/2.0/3.0*((-1.0)/dz/2.0*(psi[2][1]-psi[0][1])*tem[1][1]+1.0/dz/2.0*psi[1][1]*(tem[2][1]-tem[0][1]));
		//更新	
		for (int k=0; k<Nz; k++)
        {
            for (int n=0; n<Nn+1; n++)
            {                    
				omg[k][n] += dt/2.0*(3.0*domgdt[k][n][now]-domgdt[k][n][old]);
                tem[k][n] += dt/2.0*(3.0*dtemdt[k][n][now]-dtemdt[k][n][old]);
                dtemdt[k][n][old]=dtemdt[k][n][now];
				domgdt[k][n][old]=domgdt[k][n][now];                   
            }
        }
            
	
        
		//Thomas 算法解三对角矩阵 
        for (int n = 0 ; n < Nn+1 ; n++ )
        {
            dia[Nz-1]=1;
            dia[0]=1;
				
			for (int k = 1 ; k < Nz - 1 ; k++)
            {
                dia[k]=(n*pi/a) * (n*pi/a)+2.0/(dz*dz);
            }
                
			c[0]=0;
			d[0]=0;
												
			for (int k = 1 ; k < Nz - 1 ; k++)
            {
                double m = 1.0 / (dia[k]-sub*c[k-1]);
				c[k] = sup*m;
                d[k] = (omg[k][n] - sub * d[k-1]) * m;
            }
                
            c[Nz-1]=1;
            d[Nz-1]=0;
                
			psi[Nz-1][n]=0;
				
			for (int k = Nz-2 ; k >= 0 ; k--)
            {
                psi[k][n] = d[k] - c[k] * psi[k+1][n];
            }
        }
    } 
	
	
	//检验CFL 
	for (int i=0;i<Nn+1;i++)
	{		
	    for (int k=0 ; k < Nz ; k++)
	    {
		    for ( int n = 0 ; n < Nn+1 ; n++)
		    {
			     v[k][i] += n*pi/a*cos(n*pi*i*dx)*psi[k][n];
		    }
	    }
    }
	double max = v[0][0];
	for (int k = 0; k < Nz; k++)
	{	
	    for ( int i = 0; i < Nn+1; i++)
	    {
		    if( fabs(v[k][i])> max)
		    {
		        max=abs(v[k][i]);		       
		    }
	    }  
    }
    double vz=dz/max;
	cout<<dtemdt[1][2][now]<<" "<<dt<<endl;
	 
	//出tem数据； 
	f1.open("tem.txt");
	cout<<"输出tem文件"<<endl;
    for(int i=0; i<Nz; i++)
	{
		for (int j=0; j<Nn+1; j++)
		{
            f1<<tem[i][j]<<" ";
        }
        f1<<endl;
    }
	f1.close();

	//输出psi数据输出tem数据； 
	f2.open("psi.txt");
	cout<<"输出psi文件"<<endl;
    for(int i=0; i<Nz; i++)
	{
		for (int j=0; j<Nn+1; j++)
		{
            f2<<psi[i][j]<<" ";
        }
        f2<<endl;
    }
	f2.close();

	//输出omg数据
	f3.open("omg.txt");
	cout<<"输出omg文件"<<endl;
    for(int i=0; i<Nz; i++)
	{
		for (int j=0; j<Nn+1; j++)
		{
            f3<<omg[i][j]<<" ";
        }
        f3<<endl;
    }
	f3.close();
	

	system("pause");
    return 0;
}
