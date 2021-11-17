#ifndef ___SWEEPOUT_H
#define ___SWEEPOUT_H

template <class TEMPLATE> inline TEMPLATE myAbs(const TEMPLATE& a){
	if(a<0)return (-1.0*a);
	else return a;
}

template <class T> void sweep_out(const T* const in,T* const out,const int n){
	T* tmp_in;
	tmp_in=new T[n*n];
	memcpy(tmp_in,in,n*n*sizeof(T));					//入れ込んだ行列を変えるわけにいかないので．
	for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){								//単位行列．
		if(i==j)*(out+i*(n+1))=(T)1.0;
		else *(out+i+j*n)=(T)0.0;
	}}
	for(int i=0;i<n;i++){
		T p = *(tmp_in+i*(n+1));
		if(p==0){										//頭が0のとき．
			int j=0;
			for(j=i+1;j<n;j++){
				if(!(p = *(tmp_in+i+j*n)))break;
			}
			if(j==n)return;								//ランク落ち．ここで-1くらい返すようにすべきか??
			else{										//行交換．
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
		if(p*p!=p){										//頭を1に．
			for(int j=i;j<n;j++){*(tmp_in+j+i*n)/=p;}
			for(int j=0;j<n;j++){*(out+j+i*n)/=p;}
		}
		for(int j=0;j<n;j++){							//各列から引く
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
	for(int j=0;j<N_ARRAY;j++){					//入れ込んだ行列を変えるわけにいかないので．
		tmp_in[i][j]=in[i][j];
	}}
	
	for(int i=0;i<N_ARRAY;i++){
	for(int j=0;j<N_ARRAY;j++){								//単位行列．
		if(i==j)out[i][i]=(T2)1.0;
		else out[i][j]=0;
	}}
	for(int i=0;i<N_ARRAY;i++){
		T2 p = tmp_in[i][i];
		if(p==0){										//頭が0のとき．
			int j=0;
			for(j=i+1;j<N_ARRAY;j++){
				if(!(p = tmp_in[j][i]))break;
			}
			if(j==N_ARRAY)return;								//ランク落ち．ここで-1くらい返すようにすべきか??
			else{										//行交換．
				T2 a[N_ARRAY];
				for(int k=i;k<N_ARRAY;k++)a[k]=tmp_in[i][k];
				for(int k=i;k<N_ARRAY;k++)tmp_in[i][k]=tmp_in[j][k];
				for(int k=i;k<N_ARRAY;k++)tmp_in[j][k]=a[k];
				for(int k=0;k<N_ARRAY;k++)a[k]=out[i][k];
				for(int k=0;k<N_ARRAY;k++)out[i][k]=out[j][k];
				for(int k=0;k<N_ARRAY;k++)out[j][k]=a[k];
			}
		}
		if(p*p!=p){										//頭を1に．
			for(int j=i;j<N_ARRAY;j++){tmp_in[i][j]/=p;}
			for(int j=0;j<N_ARRAY;j++){out[i][j]/=p;}
		}
		for(int j=0;j<N_ARRAY;j++){							//各列から引く
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
	for(int j=0;j<N_ARRAY;++j){					//入れ込んだ行列を変えるわけにいかないので．
		tmp_in[i][j]=in[i][j];
	}}
	for(int i=0;i<N_ARRAY;++i){
		x[i]=b[i];
	}
	
	for(int i=0;i<N_ARRAY;++i){
		int d=i;
		double max = myAbs(tmp_in[i][i]);
		for(int j=i+1;j<N_ARRAY;++j){					//絶対値最大を探す．
			if(max<myAbs(tmp_in[j][i])){
				d=j;
				max=myAbs(tmp_in[j][i]);
			}
		}
		if(d!=i){										//行交換．
			T2 tmp=x[i];
			x[i]=x[d];
			x[d]=tmp;
			T2 a[N_ARRAY];
			for(int k=0;k<N_ARRAY;k++)a[k]=tmp_in[i][k];
			for(int k=0;k<N_ARRAY;k++)tmp_in[i][k]=tmp_in[d][k];
			for(int k=0;k<N_ARRAY;k++)tmp_in[d][k]=a[k];
		}
		if(tmp_in[i][i]){								//本文．
			for(int j=i+1;j<N_ARRAY;++j){
				tmp_in[j][i]/=tmp_in[i][i];
				for(int k=i+1;k<N_ARRAY;++k){
					tmp_in[j][k]-=tmp_in[j][i]*tmp_in[i][k];
				}
			}
		}else {
			for(int j=0;j<N_ARRAY;++j)x[j]=0;
			return;										//ランク落ち．ここで-1くらい返すようにすべきか??
		}
	}
	for(int i=0;i<N_ARRAY;++i){							//Uy=b
		for(int j=0;j<i;++j){
			x[i]-=tmp_in[i][j]*x[j];
		}
	}
	for(int i=N_ARRAY-1;i>=0;--i){						//Lx=ｙ
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
	for(int j=0;j<N_ARRAY;++j){					//入れ込んだ行列を変えるわけにいかないので．
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
// 吐き出し法により逆行列を求める関数です．
// 配列がn*n行列であることを前提にループまわしているので
// 正方行列以外だとシンタックスエラーでkillされます．たぶん．
// やはり専用のclassを作るのが安全ですね．
// 
// 代入法はsweepout.cppを参照
// １．1次元に落とした配列の先頭ポインタを入れ込む．
// ２．静的2次元配列の先頭ポインタを押し込む．
// 
// ３．ダブルポインタで二次元配列を動的に連続に確保している場合
// 	連続ならば先頭ポインタを入れ込む．
// 	このときarray,array[0],&array[0],&array[0][0]の意味が
// 	わからないならやめておいたほうがいい．
// ４．ダブルポインタで二次元配列を動的に不連続に確保している場合
// 	無理
// 
// 特殊化はちょっと勉強してきます．


使い方は，
#define N_ARRAY 100
#include "sweepout.h"
#undef N_ARRAY

でインクルードしてください．
N_ARRAYは行列の次元です．

lu：Ax=bをxについて解くとき
	lu(A,b,x);
	と代入してください．
sweep_out：Aの逆行列をA_invとするとき
		sweep_out(A,A_inv);
		と代入してください．
lu_inv:
		lu_inv(A,A_inv);
		と代入してください．

直接法で連立一次方程式を解くならば，
静的に2次元配列を確保できる領域で計算するべきだと．
それ以上ならSORとかで収束計算するべきだと考えて作りました．
sorは石田先生に聞いてください．

タイムアタックした結果は，
100行100列のマトリックスを100万回計算した結果です．
sweep_out 13141.9     3h39m01.9 s
lu         2237.67      37m17.67s
lu_inv    12652.8     3h30m52.8 s
単位は秒．旧オプテロンにて計算．コンパイルオプションなし．
計算オーダーは次数の3乗に比例したはずです．目安までに．

determinantの再帰計算で直接求めたこともありますが，
とてもじゃないですが計算できないと思った記憶があるので
書き直してはないです．
吐き出し法で逆行列を求めずに方程式を解く関数を作るべきかもしれないですが，
理論上はLU分解法のほうが早いし，作る気力がないです．

#endif

