//wuyong @ccnu.cn 2022.03.28 // wuyong@mails.ccnu.edu.cn
// ��ѧͻ�� ����ϵͳ���¶��ڴ���24���һ��Ϊ��Ϣ̬�����¶�ֻ���ǵ�24�� 
// ��ʱ������ͼ
// ʹ�õ����ޱ������ 
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#define Pi 3.141592653589793
#define NN 200 // ��������ĳ߶� 

// �Զ��庯��
double fv(double v,double m,double n,double h,double e);
double gn(double v,double n,double T);
double gm(double v,double m,double T);
double gh(double v,double h,double T);
double ge(double v,double e);
double alpha(double t); // ��ѧͻ��
double ran1(long *idum); 
double ran2(long *idum1);

// ���к��� 
void run_one();

// �����ļ�ָ�� 
FILE *fp1,*fp2,*fp3,*fp4;

// ȫ�ֱ���
double k=0.0; // ��ŷ������� 
double step=0.01;
double v_tau[NN+1][4100]={0.0}; // ����̫���޷�����ֲ�����������ȫ�ֱ���ʹ��
 
int main(){
//	fp1=fopen("v_t_D=0.1_D_xi=0.1__tau=23_T=5.3.dat","w");
//	fp2=fopen("�ŵ��ͼ��ѧ(2000-2100)_D=0.01_D_xi=0.1__tau=0_T=10.dat","w");
//	fp3=fopen("Ĥ��λ����ѧ��_D=0.1_D_xi=0.1__tau=30_T=11.0.dat","w");
	fp4=fopen("CV_ISI����ѧ��2_D=0.1_D_xi=0.1__tau=30_T=15.0-24.dat","w");	
	run_one();
//	fclose(fp1);
//	fclose(fp2);
//	fclose(fp3);
	fclose(fp4);
	
//	getchar();
	return 0;
}

void run_one(){
	double v[NN+1],m[NN+1],n[NN+1],h[NN+1],e[NN+1]; 
	double v0[NN+1],m0[NN+1],n0[NN+1],h0[NN+1],e0[NN+1],time=-200.0; // ��ʼ��200 // ��������
	int i,j,N_time=200000;   // ѭ��ʹ��
	double Iext=20.0;       // �ⲿ�̼����� 
	double T=15.0;          // �¶� 
	int a_t;
	int v_i,v_j;            // ��Ԫѭ�� 
	double D=0.01,Isyn,alpha_t[NN+1];  // ͬ��(��ѧ���)
	int tau=3000;           // ʱ���ӳ�
	// �ŵ�ʱ��
	int flag1[NN+1]={0}; 
	double firingtime[NN+1];
	// CVϵ������
	double max[NN+1];
	double T1[NN+1];
	double T2[NN+1];
	double sum[NN+1];
	double sum2[NN+1];
	int nn[NN+1];
	double delta[NN+1];
	double deltasquare[NN+1];
	double average[NN+1], variance[NN+1], var[NN+1], CV_SF[NN+1];
	double CV_ave; // ƽ��CVϵ�� 
	
	/////////////////// �ޱ�� ///////////////////
	FILE *fp_SF;
	int N_sf = NN; // ����Ľڵ��� 
	int i_sf,j_sf,L_sf,edge[N_sf+1][N_sf+1]; // ��ȡ���� 
//	int i, j;      // ѭ��ʹ��
	fp_SF=fopen("scalefree_matrix.dat","rb");	
	for(i=1; i<N_sf+1; i++){
		for(j=1; j<N_sf+1; j++){
			fscanf(fp_SF,"%d %d %d\n",&i_sf, &j_sf, &L_sf);
//			printf("%d %d %d\n", i_sf, j_sf, L_sf);
			edge[i_sf][j_sf] = L_sf; // ����ı� 
		}
	}
	fclose(fp_SF);
	
	/////////////////// ������ ///////////////////
	long *idum1, vv1,*idum2, vv2;
	double xi,D_xi=0.1;
	vv1=300;
	idum1=&vv1;
	vv2=400;
	idum2=&vv2;
	for(a_t=0; a_t<=45;a_t++){
		
		for(v_i=1; v_i<NN+1; v_i++){ // ��ʼ���ӳٿռ� 
			for(v_j=0; v_j<=tau; v_j++){		
				v_tau[v_i][tau]=0.0;
			}
		}
		
		time=0.0; // ��ʼ��ʱ�����
		
		// CVϵ�����㣨��ʼ���� 
		for(v_i=1; v_i<NN+1; v_i++){ // ��ʼ��
			max[v_i] = -70.0; 
			T1[v_i] = 0;
			T2[v_i] = 0;
			sum[v_i] = 0;
			sum2[v_i] = 0;
			nn[v_i] = 0;
		}
		
		for(v_i=1; v_i<NN+1; v_i++){ // ��������ֵ 
			v0[v_i]=0.5;
			m0[v_i]=0.5;
			n0[v_i]=1.0;
			h0[v_i]=0.6;
			e0[v_i]=0.0;
		}
			
		for(i=0; i<20000; i++){ // ȥ����̬ 
			
			for(v_i=1; v_i<NN+1; v_i++){
				
				Isyn = 0.0; 
				if(i>=10000){
					for(v_j=1; v_j<NN+1; v_j++){	// ѭ�����Ӿ���ͬ����					
						Isyn = Isyn+D*edge[v_i][v_j]*v_tau[v_j][0]*(0.0-v0[v_i]); // ͬ��(��ѧ���)
					}
				}
				
				if(i>=500) alpha_t[v_i]=alpha(time-firingtime[v_i]); // ͬ��alpha 
				
				xi = sqrt(-4.0*D_xi*step*log(ran1(idum1)))*cos(2.0*Pi*ran2(idum1));
				v[v_i]=v0[v_i]+(fv(v0[v_i],m0[v_i],n0[v_i],h0[v_i],e0[v_i])+Iext+Isyn)*step+xi;
				n[v_i]=n0[v_i]+gn(v0[v_i],n0[v_i],T)*step;
				m[v_i]=m0[v_i]+gm(v0[v_i],m0[v_i],T)*step;
				h[v_i]=h0[v_i]+gh(v0[v_i],h0[v_i],T)*step;
				e[v_i]=e0[v_i]+ge(v0[v_i],e0[v_i])*step;			
			}
			
			time=time+step; // ʱ����� 
		 	
		 	for(v_i=1; v_i<NN+1; v_i++){
				/*��������*/
			 	v0[v_i]=v[v_i];
				n0[v_i]=n[v_i];
				m0[v_i]=m[v_i];
				h0[v_i]=h[v_i];
				e0[v_i]=e[v_i];
				
				// �ŵ�ʱ�� 
				if(v0[v_i]>-10.0&&flag1[v_i]==0){
					flag1[v_i]=1;
					firingtime[v_i]=time;
				}
				if(v0[v_i]<-20.0&&flag1[v_i]==1){
					flag1[v_i]=0;
				}
			}
			
			for(v_i=1; v_i<NN+1; v_i++){ // �洢�ӳ�ֵ 
				for(v_j=0; v_j<tau; v_j++){
					v_tau[v_i][v_j]=v_tau[v_i][v_j+1];
				}
				v_tau[v_i][tau]=alpha_t[v_i];
			}	 
			
		}
		
//		printf("��ʼ\n");
		
		for(i=0; i<N_time; i++){
			for(v_i=1; v_i<NN+1; v_i++){
				Isyn = 0.0; 
				for(v_j=1; v_j<NN+1; v_j++){	// ѭ�����Ӿ���ͬ����
					Isyn = Isyn+D*edge[v_i][v_j]*v_tau[v_j][0]*(0.0-v0[v_i]); // ͬ��(��ѧ���)
				}
				
				alpha_t[v_i]=alpha(time-firingtime[v_i]); // ͬ��alpha 
				
				xi = sqrt(-4.0*D_xi*step*log(ran1(idum1)))*cos(2.0*Pi*ran2(idum1));
				v[v_i]=v0[v_i]+(fv(v0[v_i],m0[v_i],n0[v_i],h0[v_i],e0[v_i])+Iext+Isyn)*step+xi;
				n[v_i]=n0[v_i]+gn(v0[v_i],n0[v_i],T)*step;
				m[v_i]=m0[v_i]+gm(v0[v_i],m0[v_i],T)*step;
				h[v_i]=h0[v_i]+gh(v0[v_i],h0[v_i],T)*step;
				e[v_i]=e0[v_i]+ge(v0[v_i],e0[v_i])*step;			
			}
			
			time=time+step; // ʱ����� 
		 	
		 	for(v_i=1; v_i<NN+1; v_i++){
				/*��������*/
			 	v0[v_i]=v[v_i];
				n0[v_i]=n[v_i];
				m0[v_i]=m[v_i];
				h0[v_i]=h[v_i];
				e0[v_i]=e[v_i];
				
				// �ŵ�ʱ�� 
				if(v0[v_i]>-10.0&&flag1[v_i]==0) {
					flag1[v_i]=1;
					firingtime[v_i]=time;
				}
				if(flag1[v_i]==1&&v[v_i]>max[v_i]) {
					max[v_i]=v[v_i];
					T1[v_i]=time;
				}
				if(v0[v_i]<-20.0&&flag1[v_i]==1) {
					flag1[v_i]=0;
					nn[v_i]++;
					if(nn[v_i]>2) { 
						delta[v_i]=T1[v_i]-T2[v_i];
						deltasquare[v_i]=delta[v_i]*delta[v_i];
						sum2[v_i]=sum2[v_i]+deltasquare[v_i];
						sum[v_i]=sum[v_i]+delta[v_i];
					}
					T2[v_i]=T1[v_i];
					flag1[v_i]=0;
					max[v_i]=-70.0;
				}
				
	//			if(flag1[v_i]==1) fprintf(fp2,"%f %d\n",time+2000,v_i); // �ŵ��ͼ
			}
			
			for(v_i=1; v_i<NN+1; v_i++){ // �洢�ӳ�ֵ 
				for(v_j=0; v_j<tau; v_j++){
					v_tau[v_i][v_j]=v_tau[v_i][v_j+1];
				}
				v_tau[v_i][tau]=alpha_t[v_i];
			}
			
			// Ĥ��λ��ѡȡλ�ú�ʱ�䣩 
	//		if(i % 10 ==0 && i>=100000 && i <= 150000){ // ����������ݺܴ�ʱ��������ÿ��һ��ѭ���������磺10�Σ�50�Σ�100�Σ�...�������һ������		
	////			printf("%f\t%f\t%f\n", time, V, omega); // �������Ļ������������顢���ӳ���ͼ�����̣�
	//			fprintf(fp3,"%f %f %f\n", time, v0[2], v0[100]);	// ������ݵ������ļ�"**.dat"
	//		}
			
			// ʱ������ͼ
	//		if(i % 1 ==0){ // ����������ݺܴ�ʱ��������ÿ��һ��ѭ���������磺10�Σ�50�Σ�100�Σ�...�������һ������		
	////			printf("%f\t%f\t%f\n", time, V, omega); // �������Ļ������������顢���ӳ���ͼ�����̣�
	//			fprintf(fp1,"%f\t", time);	// ������ݵ������ļ�"**.dat"
	//			for(v_i=1; v_i<NN+1; v_i++){
	//				fprintf(fp1,"%f\t", v0[v_i]);
	//			}
	//			fprintf(fp1,"\n");
	//		}
		
		}
		
		CV_ave = 0; // ��ʼ��ƽ��CVϵ��
		for(v_i=1; v_i<NN+1; v_i++){
			average[v_i]=sum[v_i]/(double)(nn[v_i]-1);
			variance[v_i]=sum2[v_i]/(double)(nn[v_i]-1)-average[v_i]*average[v_i];
			var[v_i]=sqrt(variance[v_i]);
			CV_SF[v_i]=var[v_i]/average[v_i];
			
			CV_ave = CV_ave+CV_SF[v_i]/NN;
		}
		
		printf("%f %f %f %f %f %f %f\n", T, CV_ave, CV_SF[2], CV_SF[100], average[2], average[100], fabs(average[2]-average[100]));	
		fprintf(fp4,"%f %f %f %f %f %f %f\n", T, CV_ave, CV_SF[2], CV_SF[100], average[2], average[100], fabs(average[2]-average[100])); // �洢ͳ����
//		printf("%f %f\n", average[2], average[100]);
		
		T = T+0.2; // �¶�
	}
}

// �Զ��庯��
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

double alpha(double t){ // ��ѧͻ�� 
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

