#include<stdio.h>
#include<math.h>
#include<stdlib.h>
//#include<iostream.h>
#include<time.h>
#define nx 400
#define ny 400
#define pi 4*atan(1)
#define theta 15*(pi)/180.0
#define M 2.0
#define m 410
#define n 410
#define cfl 0.1
clock_t start,end;
double cpu_time_used;
int i,j,k,count=0;
double x[m][n],y[m][n],c,Vcell[m][n],L1,L2,L3,L4,Xcell[m][n],Ycell[m][n];
double nx1[m][n],ny1[m][n],nx2[m][n],ny2[m][n],nx3[m][n],ny3[m][n],nx4[m][n],ny4[m][n];
double p[m][n],d[m][n],dold[m][n],Y=1.4,R=287.0,T=300.0,a[m][n],u[m][n],v[m][n],e[m][n],t[m][n];
double U1old[m][n],U2old[m][n],U3old[m][n],U4old[m][n];

double Res,LRes,Rold,Rnew,Resd, Sum, delta=0.0;
double log10(double X);
double ML,MR,MLp,MRn,ds;
double fpm,f1p,f2p,f3p,fm,f1m,f2m,f3m;
double gpm,g1p,g2p,g3p,gm,g1m,g2m,g3m;
double tottime,dtglobalmin,dtx[m][n],dty[m][n],dtmin[m][n];
double Nx,Ny,un,ubn,Mn,Mbn,ds_x,pb[m],ub[m],vb[m],db[m],eb[m],ab[m],f,f1,f2,f3;
//double F21[m][1],F21[m][1],F31[m][1],F41[m][1];
double Fc21[m][n],Fc22[m][n],Fc23[m][n],Fc24[m][n];
double FpL21[m][n],FpL22[m][n],FpL23[m][n],FpL24[m][n];
double FpR21[m][n],FpR22[m][n],FpR23[m][n],FpR24[m][n];
double Fp21[m][n],Fp22[m][n],Fp23[m][n],Fp24[m][n];

double Fc31[m][n],Fc32[m][n],Fc33[m][n],Fc34[m][n];
double FpL31[m][n],FpL32[m][n],FpL33[m][n],FpL34[m][n];
double FpR31[m][n],FpR32[m][n],FpR33[m][n],FpR34[m][n];  
double Fp31[m][n],Fp32[m][n],Fp33[m][n],Fp34[m][n]; 
double F21[m][n],F22[m][n],F23[m][n],F24[m][n],ds_y,Msn,usn;
double F31[m][n],F32[m][n],F33[m][n],F34[m][n],ds_t,Mtn,utn;
double U1new[m][n],U2new[m][n],U3new[m][n],U4new[m][n],Mach;


int main()
{
start=clock();
FILE *fp,*fp1,*fp2;
fp=fopen("mesh.plt","w") ;
fp1=fopen("residue.plt","w");
fp2=fopen("cuputime.plt","w");

/*****GRID GENERATION*****/

	for(i=1;i<=nx+1;i++)
	{
	  for(j=1;j<=ny+1;j++)
	  {
	      x[i][j]=2.0*(i-1.0)/nx-1.0;
	      y[i][j]=2.0*(j-1.0)/ny-1.0;
	  //printf("\n %lf \t%lf",x[i][j],y[i][j]);
	  //fprintf(fp,"\n %lf \t%lf",x[i][j],y[i][j]);
	  }
	}
	
/*******VOLUME/AREA,THERE LOCATIONS AND UNIT NORMALS OF CELL*****/

	for(i=1;i<=nx;i++)
	{
	  for(j=1;j<=ny;j++)
	  {
		L1=sqrt(pow(x[i+1][j]-x[i][j],2.0)+pow(y[i+1][j]-y[i][j],2.0));
		L2=sqrt(pow(x[i+1][j+1]-x[i+1][j],2.0)+pow(y[i+1][j+1]-y[i+1][j],2.0));
		L3=sqrt(pow(x[i+1][j+1]-x[i][j+1],2.0)+pow(y[i+1][j+1]-y[i][j+1],2.0));
		L4=sqrt(pow(x[i][j+1]-x[i][j],2.0)+pow(y[i][j+1]-y[i][j],2.0));
		Vcell[i][j]=0.5*(L2+L4)*(x[i+1][j]-x[i][j]);
		Xcell[i][j]=(0.5*(x[i+1][j]+x[i][j])*L1+x[i+1][j]*L2+0.5*(x[i+1][j+1]+x[i][j+1])*L3+x[i][j]*L4)/(L1+L2+L3+L4);
		Ycell[i][j]=(0.5*(y[i+1][j]+y[i][j])*L1+0.5*(y[i+1][j+1]+y[i+1][j])*L2+0.5*(y[i+1][j+1]+y[i][j+1])*L3+0.5*(y[i][j+1]+y[i][j])*L4)/
		(L1+L2+L3+L4);
		nx1[i][j]=(y[i+1][j]-y[i][j])/L1;
		ny1[i][j]=-(x[i+1][j]-x[i][j])/L1;
		nx2[i][j]=(y[i+1][j+1]-y[i+1][j])/L2;
		ny2[i][j]= (x[i+1][j+1]-x[i+1][j])/L2;
		nx3[i][j]=-(y[i+1][j+1]-y[i][j+1])/L3;
		ny3[i][j]=(x[i+1][j+1]-x[i][j+1])/L3;
		nx4[i][j]=-(y[i][j+1]-y[i][j])/L4;
		ny4[i][j]=(x[i][j+1]-x[i][j])/L4;
		//printf("\n %lf \t%lf",L1,L2);
		L1=L2=L3=L4=0.0;
		//printf("\n %lf \t%lf \t%lf ",Vcell[i][j],Xcell[i][j],Ycell[i][j]);
		//printf("\n %lf \t%lf \t%lf \t%lf",ny1[i][j],ny2[i][j],ny3[i][j],ny4[i][j]);
		}
	}

/********INITIAL CONDITIONS**********/
for(i=nx/2+1;i<=nx;i++)
{
	for(j=ny/2+1;j<=ny;j++)
	{
		p[i][j]=1.5;
		d[i][j]=1.5;
		t[i][j]=p[i][j]/(R*d[i][j]);
		a[i][j]=sqrt((Y*p[i][j])/d[i][j]);
		u[i][j]=0.0;
		v[i][j]=0.0;
		e[i][j]=(p[i][j]/(Y-1.0))+0.5*(d[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
		//printf("\n %lf \t%lf \t%lf",d[i][j],p[i][j],e[i][j]);
	}
}
for(i=1;i<=nx/2;i++)
{
    for(j=ny/2+1;j<=ny;j++)
	{
	    p[i][j]=0.3;
		d[i][j]=0.5323;
		t[i][j]=p[i][j]/(R*d[i][j]);
		a[i][j]=sqrt((Y*p[i][j])/d[i][j]);
		u[i][j]=1.206;
		v[i][j]=0.0;
		e[i][j]=(p[i][j]/(Y-1.0))+0.5*(d[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
	
	}

}
for(i=1;i<=nx/2;i++)
{
    for(j=1;j<=ny/2;j++)
	{
	    p[i][j]=0.029;
		d[i][j]=0.1379;
		t[i][j]=p[i][j]/(R*d[i][j]);
		a[i][j]=sqrt((Y*p[i][j])/d[i][j]);
		u[i][j]=1.206;
		v[i][j]=1.206;
		e[i][j]=(p[i][j]/(Y-1.0))+0.5*(d[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
	
	}

}
for(i=nx/2+1;i<=nx;i++)
{
    for(j=1;j<=ny/2;j++)
	{
	    p[i][j]=0.3;
		d[i][j]=0.5323;
		t[i][j]=p[i][j]/(R*d[i][j]);
		a[i][j]=sqrt((Y*p[i][j])/d[i][j]);
		u[i][j]=0.0;
		v[i][j]=1.206;
		e[i][j]=(p[i][j]/(Y-1.0))+0.5*(d[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
	
	}

}

/*********INITIAL CONSERVATIVE VARIABLS*********/

for(i=1;i<=nx;i++)
{
	for(j=1;j<=ny;j++)
	{
		U1old[i][j]=d[i][j];
		U2old[i][j]=d[i][j]*u[i][j];
		U3old[i][j]=d[i][j]*v[i][j];
		U4old[i][j]=e[i][j];
	}
}

/**********TEMPORERAL DISCRETISATION********/
count=1;
Rnew=0.0;
Rold=101325.0;
Res=1.0;
tottime=0.0;
//while(Res>pow(10.0,-12.0))
//while(count<=150000)
while(tottime<=1.1)
{
	dtglobalmin=1.0;
	for(i=1;i<=nx;i++)
	{
		for(j=1;j<=ny;j++)
		{
			dtx[i][j]=(cfl*(2.0/nx))/(fabs(u[i][j])+a[i][j]);
			dty[i][j]=(cfl*(2.0/ny))/(fabs(v[i][j])+a[i][j]);
			//printf("\n%lf \t%lf",dtx[i][j],dty[i][j]);
		}
	}
	
	for(i=1;i<=nx;i++)
	{
		for(j=1;j<=ny;j++)
		{
			if(dtx[i][j]<=dty[i][j])
				dtmin[i][j]=dtx[i][j];
			else
				dtmin[i][j]=dty[i][j];
			dtglobalmin=(dtglobalmin>dtmin[i][j]? dtmin[i][j]:dtglobalmin);
		}
	}

  printf("\n %d \t%lf",count,dtglobalmin);

/*********TOP BOUNDARY CONDITION*********/

	for(i=1;i<=nx;i++)
	{
		p[i][ny+1]=p[i][ny];
		u[i][ny+1]=u[i][ny];
		v[i][ny+1]=v[i][ny];
		a[i][ny+1]=a[i][ny];
		d[i][ny+1]=d[i][ny];
		t[i][ny+1]=t[i][ny];
		e[i][ny+1]=e[i][ny];
	//	printf("\n %lf \t%lf",p[i][k+1],u[i][k+1]);
	}


	/********RIGHT BOUNDARY CONDITION**********/

	for(j=1;j<=ny;j++)
	{
		
		p[nx+1][j]=p[nx][j];
		u[nx+1][j]=u[nx][j];
		v[nx+1][j]=v[nx][j];
		a[nx+1][j]=a[nx][j];
		d[nx+1][j]=d[nx][j];
		t[nx+1][j]=t[nx][j];
		e[nx+1][j]=e[nx][j];
	//	printf("\n %lf \t%lf",p[k+1][j],u[k+1][j]);
	}

/********BOTTOM BOUNDARY CONDITION*********/

	for(i=1;i<=nx;i++)
	{
		p[i][0]=p[i][1];
		u[i][0]=u[i][1];
		v[i][0]=v[i][1];
		a[i][0]=a[i][1];
		d[i][0]=d[i][1];
		t[i][0]=t[i][1];
		e[i][0]=e[i][1];
	}
	/********LEFT BOUNDARY CONDITION**********/

	for(j=1;j<=ny;j++)
	{
		p[0][j]=p[1][j];
		u[0][j]=u[1][j];
		v[0][j]=v[1][j];
		a[0][j]=a[1][j];
		d[0][j]=d[1][j];
		t[0][j]=t[1][j];
		e[0][j]=e[1][j];
	
	}

	
    
                  /***************Fluxes at the vertical faces will be calculated at face 2 for cell no. i=1 to i=nx**********/
                    
                    for(i=0;i<=nx;i++)
                    {
                       for(j=1;j<=ny;j++)
					   {
                          						  
						  Nx = 1.0 ;
                          Ny =0.0 ;
                          ML =u[i][j]/a[i][j] ;
                          MR =u[i+1][j]/a[i+1][j] ;
                      
                          if(ML >=1.0)
                            MLp =  ML ;
                          else if(ML<=-1)
                            MLp =0.0 ;
                          else 
                            MLp =0.25*pow(ML+1.0,2) ;
                      
                          if(MR>=1.0) 
                            MRn =0.0 ;
                          else if(MR<=-1.0)
                            MRn = MR ;
                          else
                            MRn=-0.25*pow(MR-1.0,2) ;
                      
                          Mn = MLp +MRn ;
                      
                       
						if(Mn>=0.0)
                        {    
                            Fc21[i][j] = Mn*d[i][j]*a[i][j] ;  
                            Fc22[i][j] = Mn*d[i][j]*u[i][j]*a[i][j];
                            Fc23[i][j] = Mn*d[i][j]*v[i][j]*a[i][j];
                            Fc24[i][j] = Mn*d[i][j]*(Y*R*t[i][j]/(Y-1.0)+(u[i][j]*u[i][j]+v[i][j]*v[i][j])/2.0)*a[i][j];
						}
						else
						{
						    Fc21[i][j] = Mn*d[i+1][j]*a[i+1][j] ;  
                            Fc22[i][j] = Mn*d[i+1][j]*u[i+1][j]*a[i+1][j];
                            Fc23[i][j] = Mn*d[i+1][j]*v[i+1][j]*a[i+1][j];
                            Fc24[i][j] = Mn*d[i+1][j]*(Y*R*t[i+1][j]/(Y-1.0)+(u[i+1][j]*u[i+1][j]+v[i+1][j]*v[i][j])/2.0)*a[i+1][j];
						
						}
                        
						if(fabs(ML)<=1.0)
						{

							FpL21[i][j]=0;
							FpL22[i][j]=p[i][j]*pow((ML+1),2.0)*(2.0-ML)/4.0;
							FpL23[i][j]=0;
							FpL24[i][j]=0;
						}
						else
						{
                            FpL21[i][j]=0;
							FpL22[i][j]=p[i][j]*(ML+fabs(ML))/(2.0*ML);
							FpL23[i][j]=0;
							FpL24[i][j]=0;
						}

						if(fabs(MR)<=1.0)
						{
							FpR21[i][j]=0;
							FpR22[i][j]=p[i+1][j]*pow((MR-1),2.0)*(2.0+MR)/4.0;
							FpR23[i][j]=0;
							FpR24[i][j]=0;
						}
						else
						{
                            FpR21[i][j]=0;
							FpR22[i][j]=p[i+1][j]*(MR-fabs(MR))/(2.0*MR);
							FpR23[i][j]=0;
							FpR24[i][j]=0;

						}
                           
						    Fp21[i][j]=FpL21[i][j]+FpR21[i][j];
							Fp22[i][j]=FpL22[i][j]+FpR22[i][j];
							Fp23[i][j]=FpL23[i][j]+FpR23[i][j];
							Fp24[i][j]=FpL24[i][j]+FpR24[i][j];

							F21[i][j]= Fc21[i][j]+Fp21[i][j];
                            F22[i][j]= Fc22[i][j]+Fp22[i][j];
							F23[i][j]= Fc23[i][j]+Fp23[i][j];
							F24[i][j]= Fc24[i][j]+Fp24[i][j];


                            ds = fabs(y[i+1][j+1]-y[i+1][j]) ;
                            
                            F21[i][j]=F21[i][j]*ds;
				            F22[i][j]=F22[i][j]*ds;
				            F23[i][j]=F23[i][j]*ds;
			             	F24[i][j]=F24[i][j]*ds;
							
					   }
					}
                    //printf("\n %lf  \t%lf  \t%lf  \t%lf", F21[0][1],F22[0][1],F23[0][1],F24[0][1]);    
                        
             /**************Defining the fluxes at the horizontal faces 3 for cells j=0;j=ny **********/
                        
              for(i=1;i<=nx;i++)
              {
                   for(j=0;j<=ny;j++)
				   {
                      if(j==0)
                      {
                        Nx=-nx1[i][1] ;
                        Ny=-ny1[i][1] ;
					  }
                      else
                      {
                        Nx=nx3[i][j];
		                Ny=ny3[i][j];
					  }
		                un=(u[i][j]*Nx)+(v[i][j]*Ny);
	                    ML=un/a[i][j];

	               	    utn=(u[i][j+1]*Nx)+(v[i][j+1]*Ny);
		                MR=utn/a[i][j+1];
		               
		                if(ML >=1.0)
                          MLp =  ML ;
                        else if(ML<=-1)
                          MLp =0.0 ;
                        else 
                          MLp =0.25*pow(ML+1.0,2) ;
                      
                        if(MR>=1.0) 
                          MRn =0.0 ;
                        else if(MR<=-1.0)
                          MRn = MR ;
                        else
                          MRn=-0.25*pow(MR-1.0,2) ;
                      
                        Mn =MLp +MRn ;
		               
		                /***if(fabs(Mn)>=delta)
                        {    
                            Fc31[i][j] = 0.5*Mn*(d[i][j]*a[i][j]+d[i][j+1]*a[i][j+1])-0.5*fabs(Mn)*(d[i][j+1]*a[i][j+1]-d[i][j]*a[i][j]);  
                            Fc32[i][j] = 0.5*Mn*(d[i][j]*u[i][j]*a[i][j]+d[i][j+1]*u[i][j+1]*a[i][j+1])-0.5*fabs(Mn)*(d[i][j+1]*u[i][j+1]*a[i][j+1]-d[i][j]*u[i][j]*a[i][j]) ;
                            Fc33[i][j] = 0.5*Mn*(d[i][j]*v[i][j]*a[i][j]+d[i][j+1]*v[i][j+1]*a[i][j+1])-0.5*fabs(Mn)*(d[i][j+1]*v[i][j+1]*a[i][j+1]-d[i][j]*v[i][j]*a[i][j]) ;
                            Fc34[i][j] = 0.5*Mn*(d[i][j]*(Y*R*t[i][j]/(Y-1.0)+(u[i][j]*u[i][j]+v[i][j]*v[i][j])/2.0)*a[i][j]+d[i][j+1]*(Y*R*t[i][j+1]/(Y-1.0)+(u[i][j+1]*u[i][j+1]+v[i][j+1]*v[i][j+1])/2.0)*a[i][j+1])-0.5*fabs(Mn)*(d[i][j+1]*(Y*R*t[i][j+1]/(Y-1.0)+(u[i][j+1]*u[i][j+1]+v[i][j+1]*v[i][j+1])/2.0)*a[i][j+1]-d[i][j]*(Y*R*t[i][j]/(Y-1.0)+(u[i][j]*u[i][j]+v[i][j]*v[i][j])/2.0)*a[i][j]);
							
                        }
                        else 
                        {
                            Fc31[i][j] = 0.5*Mn*(d[i][j]*a[i][j]+d[i][j+1]*a[i][j+1])-0.5*((Mn*Mn+delta*delta)/(2.0*delta))*(d[i][j+1]*a[i][j+1]-d[i][j]*a[i][j]);  
                            Fc32[i][j] = 0.5*Mn*(d[i][j]*u[i][j]*a[i][j]+d[i][j+1]*u[i][j+1]*a[i][j+1])-0.5*((Mn*Mn+delta*delta)/(2.0*delta))*(d[i][j+1]*u[i][j+1]*a[i][j+1]-d[i][j]*u[i][j]*a[i][j]) ;
                            Fc33[i][j] = 0.5*Mn*(d[i][j]*v[i][j]*a[i][j]+d[i][j+1]*v[i][j+1]*a[i][j+1])-0.5*((Mn*Mn+delta*delta)/(2.0*delta))*(d[i][j+1]*v[i][j+1]*a[i][j+1]-d[i][j]*v[i][j]*a[i][j]) ;
                            Fc34[i][j] = 0.5*Mn*(d[i][j]*(Y*R*t[i][j]/(Y-1.0)+(u[i][j]*u[i][j]+v[i][j]*v[i][j])/2.0)*a[i][j]+d[i][j+1]*(Y*R*t[i][j+1]/(Y-1.0)+(u[i][j+1]*u[i][j+1]+v[i][j+1]*v[i][j+1])/2.0)*a[i][j+1])-0.5*((Mn*Mn+delta*delta)/(2.0*delta))*(d[i][j+1]*(Y*R*t[i][j+1]/(Y-1.0)+(u[i][j+1]*u[i][j+1]+v[i][j+1]*v[i][j+1])/2.0)*a[i][j+1]-d[i][j]*(Y*R*t[i][j]/(Y-1.0)+(u[i][j]*u[i][j]+v[i][j]*v[i][j])/2.0)*a[i][j]);
						}***/
						if(Mn>=0.0)
                        {    
                            Fc31[i][j] = Mn*d[i][j]*a[i][j];  
                            Fc32[i][j] = Mn*d[i][j]*u[i][j]*a[i][j] ;
                            Fc33[i][j] = Mn*d[i][j]*v[i][j]*a[i][j];
                            Fc34[i][j] = Mn*d[i][j]*(Y*R*t[i][j]/(Y-1.0)+(u[i][j]*u[i][j]+v[i][j]*v[i][j])/2.0)*a[i][j];
							
                        }
                        else 
                        {
						    Fc31[i][j] = Mn*d[i][j+1]*a[i][j+1];  
                            Fc32[i][j] = Mn*d[i][j+1]*u[i][j+1]*a[i][j+1] ;
                            Fc33[i][j] = Mn*d[i][j+1]*v[i][j+1]*a[i][j+1];
                            Fc34[i][j] = Mn*d[i][j+1]*(Y*R*t[i][j+1]/(Y-1.0)+(u[i][j+1]*u[i][j+1]+v[i][j+1]*v[i][j+1])/2.0)*a[i][j+1];
                            
						}
                        if(fabs(ML)<=1.0)
						{

							FpL31[i][j]=0;
							FpL32[i][j]=p[i][j]*pow((ML+1),2.0)*(2.0-ML)*Nx/4.0;
							FpL33[i][j]=p[i][j]*pow((ML+1),2.0)*(2.0-ML)*Ny/4.0;
							FpL34[i][j]=0;
						}
						else
						{
                            FpL31[i][j]=0;
							FpL32[i][j]=p[i][j]*Nx*(ML+fabs(ML))/(2.0*ML);
							FpL33[i][j]=p[i][j]*Ny*(ML+fabs(ML))/(2.0*ML);
							FpL34[i][j]=0;

						}

						if(fabs(MR)<=1.0)
						{

							FpR31[i][j]=0;
							FpR32[i][j]=p[i][j+1]/4.0*pow((MR-1),2.0)*(2.0+MR)*Nx;
							FpR33[i][j]=p[i][j+1]/4.0*pow((MR-1),2.0)*(2.0+MR)*Ny;
							FpR34[i][j]=0;
						}
						else
						{
                            FpR31[i][j]=0;
							FpR32[i][j]=p[i][j+1]*Nx*(MR-fabs(MR))/(2.0*MR);
							FpR33[i][j]=p[i][j+1]*Ny*(MR-fabs(MR))/(2.0*MR);
							FpR34[i][j]=0;

						}
                           
						    Fp31[i][j]=FpL31[i][j]+FpR31[i][j];
							Fp32[i][j]=FpL32[i][j]+FpR32[i][j];
							Fp33[i][j]=FpL33[i][j]+FpR33[i][j];
							Fp34[i][j]=FpL34[i][j]+FpR34[i][j];

                            F31[i][j]=Fc31[i][j]+Fp31[i][j];
							F32[i][j]=Fc32[i][j]+Fp32[i][j];
				            F33[i][j]=Fc33[i][j]+Fp33[i][j];
			             	F34[i][j]=Fc34[i][j]+Fp34[i][j];

                            ds = sqrt(pow(y[i+1][j+1] -y[i][j+1],2)+pow(x[i+1][j+1]-x[i][j+1],2))  ;
                            
                            F31[i][j]=F31[i][j]*ds;
				            F32[i][j]=F32[i][j]*ds ;
				            F33[i][j]=F33[i][j]*ds;
			             	F34[i][j]=F34[i][j]*ds;
				   }
			  }
                    
                          
             /***********For convergence***************/
                           for(i=1;i<=nx ;i++)
                           {
                               for(j=1;j<=ny;j++)
							   {
                                  dold[i][j] = d[i][j];
                               }
                           }      
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
                                             
                 
           /*************Updating the cells**********/
                            
             for(i=1;i<=nx ;i++)
             {
                   for(j=1;j<=ny;j++)
                   {
                       U1new[i][j] = U1old[i][j] - dtglobalmin*(-F31[i][j-1]+ F31[i][j] + F21[i][j]-F21[i-1][j])/Vcell[i][j] ;
                       U2new[i][j] = U2old[i][j] - dtglobalmin*(-F32[i][j-1]+ F32[i][j] + F22[i][j]-F22[i-1][j])/Vcell[i][j] ;
                       U3new[i][j] = U3old[i][j] - dtglobalmin*(-F33[i][j-1]+ F33[i][j] + F23[i][j]-F23[i-1][j])/Vcell[i][j] ;
                       U4new[i][j] = U4old[i][j] - dtglobalmin*(-F34[i][j-1]+ F34[i][j] + F24[i][j]-F24[i-1][j])/Vcell[i][j] ;
				   }
			 }
	/************CONVERGENCE*******/

	for(i=1;i<=nx;i++)
	{
		for(j=1;j<=ny;j++)
		{
			U1old[i][j]=U1new[i][j];
			U2old[i][j]=U2new[i][j];
			U3old[i][j]=U3new[i][j];
			U4old[i][j]=U4new[i][j];
		  //	printf("\n %lf \t%lf",U1old[i][j],U1new[i][j]);
		}
	}

/************NEW CONSERVATIVE VARIABLES********/

	for(i=1;i<=nx;i++)
	{
		for(j=1;j<=ny;j++)
		{
			//printf("\n%lf",U1old[i][j]);
			d[i][j]=U1old[i][j];
			u[i][j]=U2old[i][j]/d[i][j];
			v[i][j]=U3old[i][j]/d[i][j];
			e[i][j]=U4old[i][j];
			p[i][j]=(Y-1.0)*(e[i][j]-(0.5*d[i][j]*(pow(u[i][j],2.0)+pow(v[i][j],2.0))));
			a[i][j]=sqrt(Y*p[i][j]/d[i][j]);
			t[i][j]=p[i][j]/(d[i][j]*R);
	  //	printf("\n %lf \t%lf \t%lf \t%lf \t%lf \t%lf ",d[i][j],u[i][j],v[i][j],p[i][j]),a[i][j],t[i][j]);
		}
	}

	//for(i=1;i<=nx;i++)
	//{
		//for(j=1;j<=ny;j++)
		//{
			//Resd=p[i][j];
		  //printf("\n%lf",Resd) ;
			//if(Resd>Rnew)
				//{
				//Rnew=Resd;
				//}
		//}
	//}
			//Res=(Rnew-Rold)/Rnew;
			//printf("\n%lf \t%lf",Rnew,Rold);
			//Rold=Rnew;
			//printf("\n %d \t%0.12lf",count,Res);
	Sum=0.0;
	for(i=2;i<=nx;i++)
	{
		for(j=1;j<=ny;j++)
		{
			Sum = Sum+pow(((dold[i][j]-d[i][j])/dold[i][j]),2.0);
		}
	}
	
		Res=sqrt(Sum/((nx-1)*ny));
	    LRes=log10(Res);
	    
		fprintf(fp1,"\n %d \t%0.12lf",count,LRes);
		count++;
		tottime=tottime+dtglobalmin;

}

fprintf(fp,"VARIABLES=\"X\",\"Y\",\"PR\",\"DENSITY\",\"U-vel\",\"V-vel\",\"TEMP\",\"Mach no\",\nZONE T=\"BLOCK1\",I=%d,J=%d,F=POINT\n",nx/2,ny/2);

for(j=1;j<=ny/2;j++)
{
	for(i=1;i<=nx/2;i++)
	{
		Mach=sqrt(pow(u[i][j],2.0)+pow(v[i][j],2.0))/a[i][j];
		fprintf(fp,"\n %lf \t%lf \t%lf \t%0.7lf \t%lf \t%lf \t%lf \t%lf ",Xcell[i][j],Ycell[i][j],p[i][j],d[i][j],u[i][j],v[i][j],t[i][j],Mach);

	}
}

fclose(fp);
fclose(fp1);
end=clock();
cpu_time_used=((double)(end-start))/CLOCKS_PER_SEC;
fprintf(fp2,"Total CPU time used=%0.12lf",cpu_time_used);
fclose(fp2);
  return 0;
}

