//wuyong @ccnu.cn 2022.03.28 // wuyong@mails.ccnu.edu.cn
// 延迟与温度系数统计量 （此系统在温度在大于24左右会变为静息态，故温度只考虑到24） 
// 使用的是无标度网络 
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#define Pi 3.141592653589793
#define NN 200 // 定义网络的尺度 

// 自定义函数
double fv(double v,double m,double n,double h,double e);
double gn(double v,double n,double T);
double gm(double v,double m,double T);
double gh(double v,double h,double T);
double ge(double v,double e);
double alpha(double t); // 化学突触
double ran1(long *idum); 
double ran2(long *idum1);

// 运行函数 
void run_one();

// 定义文件指针 
FILE *fp1;

// 全局变量
double k=0.0; // 电磁反馈增益 
double step=0.01; // 步长 
double v_tau[NN+1][6200]={0.0}; // 数组太大无法定义局部变量，定义全局变量使用
 
int main(){
	fp1=fopen("tau-T(23.8-24.0)-R化学_D=0.1_D_xi=0.1_I=20.dat","wb");
	run_one();
	fclose(fp1);
	
	return 0;
}

void run_one(){
	double v[NN+1],m[NN+1],n[NN+1],h[NN+1],e[NN+1]; 
	double v0[NN+1],m0[NN+1],n0[NN+1],h0[NN+1],e0[NN+1],time=-200.0; // 创建变量
	int i,j,N_time=300000;   // 循环使用
	double Iext=20.0;       // 外部刺激电流 
	double T=23.8;          // 温度 
	int a_t;
	int v_i,v_j;            // 神经元循环 
	double D=0.1,Isyn,alpha_t[NN+1];      // 同步
	int tau=0;           // 时间延迟
	// 放电时间
	int flag1[NN+1]={0}; 
	double firingtime[NN+1];
	
	/////////////////// 无标度 ///////////////////
	FILE *fp_SF;
	int N_sf = NN; // 网络的节点数 
	int i_sf,j_sf,L_sf,edge[N_sf+1][N_sf+1]; // 读取数据 
//	int i, j;      // 循环使用
	fp_SF=fopen("scalefree_matrix.dat","rb");	
	for(i=1; i<N_sf+1; i++){
		for(j=1; j<N_sf+1; j++){
			fscanf(fp_SF,"%d %d %d\n",&i_sf, &j_sf, &L_sf);
//			printf("%d %d %d\n", i_sf, j_sf, L_sf);
			edge[i_sf][j_sf] = L_sf; // 网络的边 
		}
	}
	fclose(fp_SF);
	
	/////////////////// 白噪声 ///////////////////
	long *idum1, vv1,*idum2, vv2;
	double xi,D_xi=0.1;
	vv1=300;
	idum1=&vv1;
	vv2=400;
	idum2=&vv2;
	
	for(a_t=0; a_t<=1;a_t++){ 
		tau=0; // 延迟时间 
		while(tau<=6000){				
			for(v_i=1; v_i<NN+1; v_i++){ // 初始化延迟空间 
				for(v_j=0; v_j<=tau; v_j++){		
					v_tau[v_i][tau]=0.0;
				}
			}
			
			time=-200.0; // 初始化时间变量 
			
			for(v_i=1; v_i<NN+1; v_i++){ // 变量赋初值 
				v0[v_i]=0.5;
				m0[v_i]=0.5;
				n0[v_i]=1.0;
				h0[v_i]=0.6;
				e0[v_i]=0.0;
			}
				
			for(i=0; i<20000; i++){ // 去掉暂态 
				for(v_i=1; v_i<NN+1; v_i++){
					Isyn = 0.0; 
					if(i>=10000){
						for(v_j=1; v_j<NN+1; v_j++){	// 循环连接矩阵（同步）
							Isyn = Isyn+D*edge[v_i][v_j]*v_tau[v_j][0]*(0.0-v0[v_i]); // 同步(化学耦合)
						}
					}
					 
					if(i>=500) alpha_t[v_i]=alpha(time-firingtime[v_i]); // 同步alpha
					 
					xi = sqrt(-4.0*D_xi*step*log(ran1(idum1)))*cos(2.0*Pi*ran2(idum1));
					v[v_i]=v0[v_i]+(fv(v0[v_i],m0[v_i],n0[v_i],h0[v_i],e0[v_i])+Iext+Isyn)*step+xi;
					n[v_i]=n0[v_i]+gn(v0[v_i],n0[v_i],T)*step;
					m[v_i]=m0[v_i]+gm(v0[v_i],m0[v_i],T)*step;
					h[v_i]=h0[v_i]+gh(v0[v_i],h0[v_i],T)*step;
					e[v_i]=e0[v_i]+ge(v0[v_i],e0[v_i])*step;			
				}
				
				time=time+step; // 时间迭代 
			 	
			 	for(v_i=1; v_i<NN+1; v_i++){
					/*迭代变量*/
				 	v0[v_i]=v[v_i];
					n0[v_i]=n[v_i];
					m0[v_i]=m[v_i];
					h0[v_i]=h[v_i];
					e0[v_i]=e[v_i];
					
					// 放电时间 
					if(v0[v_i]>0.0&&flag1[v_i]==0){
						flag1[v_i]=1;
						firingtime[v_i]=time;
					}
					if(v0[v_i]<-20.0&&flag1[v_i]==1){
						flag1[v_i]=0;
					}
				}
				
				for(v_i=1; v_i<NN+1; v_i++){ // 存储延迟值 
					for(v_j=0; v_j<tau; v_j++){
						v_tau[v_i][v_j]=v_tau[v_i][v_j+1];
					}
					v_tau[v_i][tau]=alpha_t[v_i];
				}
			}
			
			/******* 统计量 *******/
			double R_syn=0,F_syn,R_aa=0.0,R_bb=0.0,R_cc[NN+1],R_dd[NN+1],R_c=0,R_d=0,g_syn;		
			for(v_i=1; v_i<NN+1; v_i++){
				R_cc[v_i]=0.0;
				R_dd[v_i]=0.0;
			}
				
			for(i=0; i<N_time; i++){
				F_syn=0; // 统计量
				for(v_i=1; v_i<NN+1; v_i++){
					/******* 统计量 *******/			
					F_syn+=v0[v_i]/(double)NN;
					R_cc[v_i]+=v0[v_i]*v0[v_i]/(double(N_time)+1.);
					R_dd[v_i]+=v0[v_i]/(double(N_time)+1.);
						
					Isyn = 0.0; 
					for(v_j=1; v_j<NN+1; v_j++){	// 循环连接矩阵（同步）
						Isyn = Isyn+D*edge[v_i][v_j]*v_tau[v_j][0]*(0.0-v0[v_i]); // 同步(化学耦合)
					}
					
					alpha_t[v_i]=alpha(time-firingtime[v_i]); // 同步alpha
					
					xi = sqrt(-4.0*D_xi*step*log(ran1(idum1)))*cos(2.0*Pi*ran2(idum1));
					v[v_i]=v0[v_i]+(fv(v0[v_i],m0[v_i],n0[v_i],h0[v_i],e0[v_i])+Iext+Isyn)*step+xi;
					n[v_i]=n0[v_i]+gn(v0[v_i],n0[v_i],T)*step;
					m[v_i]=m0[v_i]+gm(v0[v_i],m0[v_i],T)*step;
					h[v_i]=h0[v_i]+gh(v0[v_i],h0[v_i],T)*step;
					e[v_i]=e0[v_i]+ge(v0[v_i],e0[v_i])*step;			
				}
				/******* 统计量 *******/
				R_aa += F_syn*F_syn/(double(N_time)+1.);
				R_bb += F_syn/(double(N_time)+1.);
				
				time=time+step; // 时间迭代
			 	
			 	for(v_i=1; v_i<NN+1; v_i++){
					/*迭代变量*/
				 	v0[v_i]=v[v_i];
					n0[v_i]=n[v_i];
					m0[v_i]=m[v_i];
					h0[v_i]=h[v_i];
					e0[v_i]=e[v_i];
					
					// 放电时间 
					if(v0[v_i]>0.0&&flag1[v_i]==0) {
						flag1[v_i]=1;
						firingtime[v_i]=time;
					}
					if(v0[v_i]<-20.0&&flag1[v_i]==1) {
						flag1[v_i]=0;
					}
				}
				
				for(v_i=1; v_i<NN+1; v_i++){ // 存储延迟值 
					for(v_j=0; v_j<tau; v_j++){
						v_tau[v_i][v_j]=v_tau[v_i][v_j+1];
					}
					v_tau[v_i][tau]=alpha_t[v_i];
				}							
			}
			/******* 统计量 *******/
			F_syn=0; 		
			for(v_i=1; v_i<NN+1; v_i++){ // 加上最后一次的值 
				F_syn+=v0[v_i]/(double)NN;
				R_cc[v_i]+=v0[v_i]*v0[v_i]/(double(N_time)+1.);
				R_dd[v_i]+=v0[v_i]/(double(N_time)+1.);
			}
			R_aa += F_syn*F_syn/(double(N_time)+1.);
			R_bb += F_syn/(double(N_time)+1.);
			
			for(v_i=1; v_i<NN+1; v_i++){
				R_c+=R_cc[v_i]/(double)NN;
				R_d+=R_dd[v_i]*R_dd[v_i]/(double)NN;
			}
			
			g_syn=R_c-R_d;
			if(g_syn>-0.00000001&&g_syn<0.00000001){
				R_syn=0.0;
			} 
			else {
				R_syn=(R_aa-R_bb*R_bb)/g_syn; // synchronization factor 
			}
			
			fprintf(fp1,"%f\t%f\t%f\n", T, tau*step, R_syn); // 存储统计量
			
			if(tau%100 == 0) printf("%f\t%f\t%f\n", T, tau*step, R_syn);
			
			tau = tau+30;
		}
	T = T+0.2; // 温度	
	}
	
}

// 自定义函数
double fv(double v,double m,double n,double h,double e){
	double cm=1.0/*uF/cm2*/; 
	double result;
	result=(-36.0*n*n*n*n*(v+77.0)-120.0*m*m*m*h*(v-50.0)-0.3*(v+54.4)-k*(0.4+3.0*0.02*e*e)*v)/cm;
	return result;
}

double gn(double v,double n,double T){
	double result;
	result=0.01*(v+55.0)*pow(3,(T-6.3)/10)/(1-exp(-(v+55)/10))*(1-n)-0.125*exp(-(v+65.0)/80)*pow(3,(T-6.3)/10)*n;
	return result;
}

double gm(double v,double m,double T){
	double result;
	result=0.1*pow(3,(T-6.3)/10)*(v+40.0)/(1.0-exp(-(v+40)/10.0))*(1.0-m)-4.0*pow(3,(T-6.3)/10)*exp(-(v+65.0)/18.0)*m;
	return result;
}

double gh(double v,double h,double T){
	double result;
	result=0.07*pow(3,(T-6.3)/10)*exp(-(v+65.0)/20.0)*(1.0-h)-1.0*pow(3,(T-6.3)/10)/(1.0+exp(-(v+35.0)/10.0))*h;
	return result;
}

double ge(double v,double e){
	double result;
	double k1=0.001,k2=0.01;
	result=k1*v-k2*e;
	return result;
}

double alpha(double t){ // 化学突触 
	double result;
	double tau=2; 
	result=t*exp(-t/tau)/tau;
	return result;
}

double ran1(long *idum){
	int const IA=16807;
	long const IM=2147483647;
	double const AM=1.0/IM;
	long const IQ=127773;
	int const IR=2836;
	int const NTAB=32;
	long const NDIV=1+(IM-1)/NTAB;
	double const EPS=1.2e-7;
	double const RNMX=1.0-EPS;

	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <=0 || !iy) {
		if (-(*idum)<1) *idum=1;
		else *idum =-(*idum);
		for (j=NTAB+7; j>=0; j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum<0) *idum+=IM;
			if (j<NTAB) iv[j]=*idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum<0)*idum+=IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j]=*idum;
	if ((temp=AM*iy)>RNMX) return RNMX;
	else return temp;
}

double ran2(long *idum1){
	long const IM1=2147483563;
	long const IM2=2147483399;
	double const AM=1.0/IM1;
	long const IMM1=IM1-1;
	int const IA1=40014;
	int const IA2=40692;
	int const IQ1=53668;
	int const IQ2=52774;
	int const IR1=12211;
	int const IR2=3791;
	int const NTAB=32;
	long const NDIV=(1+IMM1/NTAB);
	double const EPS=1.2e-7;
	double const RNMX=1.0-EPS;
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum1<=0)
	{
		if (-(*idum1)<1)*idum1=1;
		else *idum1=-(*idum1);
		idum2=(*idum1);
		for (j=NTAB+7; j>=0; j--)
		{
			k=(*idum1)/IQ1;
			*idum1 =IA1*(*idum1-k*IQ1)-k*IR1;
			if (*idum1<0) *idum1+=IM1;
			if (j<NTAB) iv[j]=*idum1;
		}
		iy=iv[0];
	}

	k=(*idum1)/IQ1;
	*idum1=IA1*(*idum1-k*IQ1)-k*IR1;
	if (*idum1<0)*idum1+=IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2<0) idum2+=IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j]=*idum1;

	if (iy<1)iy+=IMM1;
	if ((temp=AM*iy)>RNMX) return RNMX;
	else return temp;
}

