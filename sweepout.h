#ifndef ___SWEEPOUT_H
#define ___SWEEPOUT_H

template <class TEMPLATE> inline TEMPLATE myAbs(const TEMPLATE& a){
	if(a<0)return (-1.0*a);
	else return a;
}

template <class T> void sweep_out(const T* const in,T* const out,const int n){
	T* tmp_in;
	tmp_in=new T[n*n];
	memcpy(tmp_in,in,n*n*sizeof(T));					//���ꍞ�񂾍s���ς���킯�ɂ����Ȃ��̂ŁD
	for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){								//�P�ʍs��D
		if(i==j)*(out+i*(n+1))=(T)1.0;
		else *(out+i+j*n)=(T)0.0;
	}}
	for(int i=0;i<n;i++){
		T p = *(tmp_in+i*(n+1));
		if(p==0){										//����0�̂Ƃ��D
			int j=0;
			for(j=i+1;j<n;j++){
				if(!(p = *(tmp_in+i+j*n)))break;
			}
			if(j==n)return;								//�����N�����D������-1���炢�Ԃ��悤�ɂ��ׂ���??
			else{										//�s�����D
				T* a;
				a = new T[n];
				for(int k=i;k<n;k++)a[k]=*(tmp_in+k+i*n);
				for(int k=i;k<n;k++)*(tmp_in+k+i*n)=*(tmp_in+k+j*n);
				for(int k=i;k<n;k++)*(tmp_in+k+j*n)=a[k];
				for(int k=0;k<n;k++)a[k]=*(out+k+i*n);
				for(int k=0;k<n;k++)*(out+k+i*n)=*(out+k+j*n);
				for(int k=0;k<n;k++)*(out+k+j*n)=a[k];
				delete[] a;
			}
		}
		if(p*p!=p){										//����1�ɁD
			for(int j=i;j<n;j++){*(tmp_in+j+i*n)/=p;}
			for(int j=0;j<n;j++){*(out+j+i*n)/=p;}
		}
		for(int j=0;j<n;j++){							//�e�񂩂����
			if(j==i)continue;
			if(p = *(tmp_in+i+j*n)){
				for(int k=i;k<n;k++){*(tmp_in+k+j*n)-=(p**(tmp_in+k+i*n));}
				for(int k=0;k<n;k++){*(out+k+j*n)-=(p**(out+k+i*n));}
			}
		}
	}
	delete[] tmp_in;
	return;
}

//#define N_ARRAY 4
template <class T1,class T2> void sweep_out(T1 in[N_ARRAY][N_ARRAY],T2 out[N_ARRAY][N_ARRAY]){
	T2 tmp_in[N_ARRAY][N_ARRAY];
	for(int i=0;i<N_ARRAY;i++){
	for(int j=0;j<N_ARRAY;j++){					//���ꍞ�񂾍s���ς���킯�ɂ����Ȃ��̂ŁD
		tmp_in[i][j]=in[i][j];
	}}
	
	for(int i=0;i<N_ARRAY;i++){
	for(int j=0;j<N_ARRAY;j++){								//�P�ʍs��D
		if(i==j)out[i][i]=(T2)1.0;
		else out[i][j]=0;
	}}
	for(int i=0;i<N_ARRAY;i++){
		T2 p = tmp_in[i][i];
		if(p==0){										//����0�̂Ƃ��D
			int j=0;
			for(j=i+1;j<N_ARRAY;j++){
				if(!(p = tmp_in[j][i]))break;
			}
			if(j==N_ARRAY)return;								//�����N�����D������-1���炢�Ԃ��悤�ɂ��ׂ���??
			else{										//�s�����D
				T2 a[N_ARRAY];
				for(int k=i;k<N_ARRAY;k++)a[k]=tmp_in[i][k];
				for(int k=i;k<N_ARRAY;k++)tmp_in[i][k]=tmp_in[j][k];
				for(int k=i;k<N_ARRAY;k++)tmp_in[j][k]=a[k];
				for(int k=0;k<N_ARRAY;k++)a[k]=out[i][k];
				for(int k=0;k<N_ARRAY;k++)out[i][k]=out[j][k];
				for(int k=0;k<N_ARRAY;k++)out[j][k]=a[k];
			}
		}
		if(p*p!=p){										//����1�ɁD
			for(int j=i;j<N_ARRAY;j++){tmp_in[i][j]/=p;}
			for(int j=0;j<N_ARRAY;j++){out[i][j]/=p;}
		}
		for(int j=0;j<N_ARRAY;j++){							//�e�񂩂����
			if(j==i)continue;
			if(p = tmp_in[j][i]){
				for(int k=i;k<N_ARRAY;k++){tmp_in[j][k]-=(p*tmp_in[i][k]);}
				for(int k=0;k<N_ARRAY;k++){out[j][k]-=(p*out[i][k]);}
			}
		}
	}
	return;
}
template <class T1,class T2,class T3> void lu(const T1 in[N_ARRAY][N_ARRAY],const T2 b[N_ARRAY], T3 x[N_ARRAY]){
	T2 tmp_in[N_ARRAY][N_ARRAY];
	for(int i=0;i<N_ARRAY;++i){
	for(int j=0;j<N_ARRAY;++j){					//���ꍞ�񂾍s���ς���킯�ɂ����Ȃ��̂ŁD
		tmp_in[i][j]=in[i][j];
	}}
	for(int i=0;i<N_ARRAY;++i){
		x[i]=b[i];
	}
	
	for(int i=0;i<N_ARRAY;++i){
		int d=i;
		double max = myAbs(tmp_in[i][i]);
		for(int j=i+1;j<N_ARRAY;++j){					//��Βl�ő��T���D
			if(max<myAbs(tmp_in[j][i])){
				d=j;
				max=myAbs(tmp_in[j][i]);
			}
		}
		if(d!=i){										//�s�����D
			T2 tmp=x[i];
			x[i]=x[d];
			x[d]=tmp;
			T2 a[N_ARRAY];
			for(int k=0;k<N_ARRAY;k++)a[k]=tmp_in[i][k];
			for(int k=0;k<N_ARRAY;k++)tmp_in[i][k]=tmp_in[d][k];
			for(int k=0;k<N_ARRAY;k++)tmp_in[d][k]=a[k];
		}
		if(tmp_in[i][i]){								//�{���D
			for(int j=i+1;j<N_ARRAY;++j){
				tmp_in[j][i]/=tmp_in[i][i];
				for(int k=i+1;k<N_ARRAY;++k){
					tmp_in[j][k]-=tmp_in[j][i]*tmp_in[i][k];
				}
			}
		}else {
			for(int j=0;j<N_ARRAY;++j)x[j]=0;
			return;										//�����N�����D������-1���炢�Ԃ��悤�ɂ��ׂ���??
		}
	}
	for(int i=0;i<N_ARRAY;++i){							//Uy=b
		for(int j=0;j<i;++j){
			x[i]-=tmp_in[i][j]*x[j];
		}
	}
	for(int i=N_ARRAY-1;i>=0;--i){						//Lx=��
		for(int j=N_ARRAY-1;j>i;--j){
			x[i]-=tmp_in[i][j]*x[j];
		}
		x[i]/=tmp_in[i][i];
	}
	
	return;
}
template <class T1,class T2> void lu_inv(const T1 in[N_ARRAY][N_ARRAY],T2 out[N_ARRAY][N_ARRAY]){
	T2 tmp_in[N_ARRAY][N_ARRAY];
	int index[N_ARRAY];
	for(int i=0;i<N_ARRAY;++i){
	for(int j=0;j<N_ARRAY;++j){					//���ꍞ�񂾍s���ς���킯�ɂ����Ȃ��̂ŁD
		tmp_in[i][j]=in[i][j];
	}}
	for(int i=0;i<N_ARRAY;++i){
		index[i]=i;
	}
	
	for(int i=0;i<N_ARRAY;++i){
		int d=i;
		double max = myAbs(tmp_in[i][i]);
		for(int j=i+1;j<N_ARRAY;++j){
			if(max<myAbs(tmp_in[j][i])){
				d=j;
				max=myAbs(tmp_in[j][i]);
			}
		}
		if(d!=i){
			T2 tmp=index[i];
			index[i]=index[d];
			index[d]=tmp;
			T2 a[N_ARRAY];
			for(int k=0;k<N_ARRAY;k++)a[k]=tmp_in[i][k];
			for(int k=0;k<N_ARRAY;k++)tmp_in[i][k]=tmp_in[d][k];
			for(int k=0;k<N_ARRAY;k++)tmp_in[d][k]=a[k];
		}
		if(tmp_in[i][i]){
			for(int j=i+1;j<N_ARRAY;++j){
				tmp_in[j][i]/=tmp_in[i][i];
				for(int k=i+1;k<N_ARRAY;++k){
					tmp_in[j][k]-=tmp_in[j][i]*tmp_in[i][k];
				}
			}
		}else {
			return;
		}
	}
	T2 sum;
	for (int k=0;k<N_ARRAY;k++){
		for (int i=0;i<N_ARRAY;i++){
			if(index[i]==k){sum=1;}
			else {sum=0;}
			for(int j=0;j<i;j++){
				sum-=tmp_in[i][j]*out[j][k];
			}
			out[i][k]=sum;
		}
		for(int i=N_ARRAY-1;i>=0;i--){
			sum=out[i][k];
			for(int j=i+1;j<N_ARRAY;j++){
				sum-=tmp_in[i][j]*out[j][k];
			}
			out[i][k]=sum/tmp_in[i][i];
		}
	}
	return;
}
//#undef N_ARRAY





#endif

#if 0
// �f���o���@�ɂ��t�s������߂�֐��ł��D
// �z��n*n�s��ł��邱�Ƃ�O��Ƀ��[�v�܂킵�Ă���̂�
// �����s��ȊO���ƃV���^�b�N�X�G���[��kill����܂��D���Ԃ�D
// ��͂��p��class�����̂����S�ł��ˁD
// 
// ����@��sweepout.cpp���Q��
// �P�D1�����ɗ��Ƃ����z��̐擪�|�C���^����ꍞ�ށD
// �Q�D�ÓI2�����z��̐擪�|�C���^���������ށD
// 
// �R�D�_�u���|�C���^�œ񎟌��z��𓮓I�ɘA���Ɋm�ۂ��Ă���ꍇ
// 	�A���Ȃ�ΐ擪�|�C���^����ꍞ�ށD
// 	���̂Ƃ�array,array[0],&array[0],&array[0][0]�̈Ӗ���
// 	�킩��Ȃ��Ȃ��߂Ă������ق��������D
// �S�D�_�u���|�C���^�œ񎟌��z��𓮓I�ɕs�A���Ɋm�ۂ��Ă���ꍇ
// 	����
// 
// ���ꉻ�͂�����ƕ׋����Ă��܂��D


�g�����́C
#define N_ARRAY 100
#include "sweepout.h"
#undef N_ARRAY

�ŃC���N���[�h���Ă��������D
N_ARRAY�͍s��̎����ł��D

lu�FAx=b��x�ɂ��ĉ����Ƃ�
	lu(A,b,x);
	�Ƒ�����Ă��������D
sweep_out�FA�̋t�s���A_inv�Ƃ���Ƃ�
		sweep_out(A,A_inv);
		�Ƒ�����Ă��������D
lu_inv:
		lu_inv(A,A_inv);
		�Ƒ�����Ă��������D

���ږ@�ŘA���ꎟ�������������Ȃ�΁C
�ÓI��2�����z����m�ۂł���̈�Ōv�Z����ׂ����ƁD
����ȏ�Ȃ�SOR�Ƃ��Ŏ����v�Z����ׂ����ƍl���č��܂����D
sor�͐Γc�搶�ɕ����Ă��������D

�^�C���A�^�b�N�������ʂ́C
100�s100��̃}�g���b�N�X��100����v�Z�������ʂł��D
sweep_out 13141.9     3h39m01.9 s
lu         2237.67      37m17.67s
lu_inv    12652.8     3h30m52.8 s
�P�ʂ͕b�D���I�v�e�����ɂČv�Z�D�R���p�C���I�v�V�����Ȃ��D
�v�Z�I�[�_�[�͎�����3��ɔ�Ⴕ���͂��ł��D�ڈ��܂łɁD

determinant�̍ċA�v�Z�Œ��ڋ��߂����Ƃ�����܂����C
�ƂĂ�����Ȃ��ł����v�Z�ł��Ȃ��Ǝv�����L��������̂�
���������Ă͂Ȃ��ł��D
�f���o���@�ŋt�s������߂��ɕ������������֐������ׂ���������Ȃ��ł����C
���_���LU����@�̂ق����������C���C�͂��Ȃ��ł��D

#endif

