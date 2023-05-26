#include<algorithm> 
#include<iostream>
#include<cstring>
#include<cstdio>
#include<vector>

class Polynomial{
	public:
	long long Mod = 998244353;
	int G = 3, iG = 332748118;
	int MS;
	long long*Inv;  //逆元 
	long long Sz,*R; 
	long long InvSz;	//	NTT
	int N,M;
	long long*A1,*B1;
	
		Polynomial(int sum){
			MS=std::max(64,sum*2);
			Inv =new long long[MS];
			R=new long long[MS];
			A1=new long long[MS];
			B1=new long long[MS];
			
			for(int i=0;i<MS;i++){
				Inv[i]=0;
				R[i]=0;
				A1[i]=0;
				B1[i]=0;
			}
			InvSz=0;Sz=0;N=0;M=0;
			Init(MS);
		}
		
	void get(std::vector<int>&A,std::vector<int>&B);
	std::vector<int>& add();
	std::vector<int>& sub();
	std::vector<int>& mul();
	static int cmp(std::vector<int>&A,std::vector<int>&B); 
	
	void Clear();
	long long qPow(long long b, int e);
	void Init(int N);
	void InitFNTT(int N);
	void FNTT(long long *A, int Ty);
	void PolyMul(long long *A,long long *B, int deg);
	void PolyInv(long long *A, int N, long long *B);
	void PolyLn(long long *A, int N, long long *B);
	void PolyExp(long long *A, int N, long long *B);
	
};

class Polydiv{
	public:
	static std::vector<int>& Polymul_tool(std::vector<int> a,std::vector<int> b);	
	static void print_vector(std::vector<int>&a);
	static int ceilog(int x);
	static std::vector<int>& Inverse(std::vector<int>& A);
	static std::vector<int>& div(std::vector<int> A,std::vector<int> B,std::vector<int>&Z,int*tes);
};

class Integer{
	protected:
		std::vector<int> data;
		int flag=1;
	public:
	Integer();
	Integer(long long a);
	Integer(int flag,std::vector<int>&a);
	
	
	void check();
	friend std::istream& operator>>(std::istream& in,Integer&a); 
     friend std::ostream& operator<<(std::ostream& out,Integer&a);
	
	Integer& operator=(Integer& a);
	int operator <(Integer& a);
	int operator <=(Integer& a);
	int operator >(Integer& a);
	int operator >=(Integer& a);
 	Integer& operator+(Integer& a);
	Integer& operator-(Integer& a); 
	Integer& operator-();
	Integer& operator/(Integer& a);
	Integer& operator%(Integer& a);
	Integer& operator*(Integer& a);
	std::vector<int>& Polymul_tool(std::vector<int> a,std::vector<int> b);
	int size(){return data.size();}

};


int main(){
	Integer a,b;











	
}

Integer::Integer(){
	 	flag=1;
	 	data.push_back(0);
}

Integer::Integer(long long a){
	  	if(a<0)flag=-1;else flag=1;
	  	a=abs(a);
	  	while(a){
	  		data.push_back(a%10);
	  		a/=10;
		  }
}

Integer::Integer(int flag,std::vector<int>&a){
	data=a;
	this->flag=1;
}

int Integer::operator<(Integer& a){
	if(Polynomial::cmp(data,a.data)<0)return 1;else return 0;
	
}

int Integer::operator<=(Integer& a){
	if(Polynomial::cmp(data,a.data)>0)return 0;else return 1;
}

int Integer::operator>(Integer& a){
	if(Polynomial::cmp(data,a.data)>0)return 1;else return 0;
	
}

int Integer::operator>=(Integer& a){
	if(Polynomial::cmp(data,a.data)<0)return 0;else return 1;
}

void Integer::check(){
	for(int i=data.size()-1;i>=1;i--)if(data[i]==0)data.pop_back();else break;
	
	if(data.size()==1&&data[0]==0)flag=1;
	if(data.size()==0){
		data.push_back(0);
		flag=1;
	}
	
}

Integer& Integer:: operator=(Integer& a){
	data=a.data;
	flag=a.flag;
	return *this;
}

Integer& Integer:: operator+(Integer& a){
	if(flag==1&&a.flag==1){
		Integer*p;
 		p=new Integer;
		p->flag=1;
		Polynomial poly(data.size()+a.data.size()-1);
		poly.get(data,a.data);
		p->data=poly.add();
		p->check();
		return (*p);
	}else if(flag==-1&&a.flag==-1){
		Integer*p;
 		p=new Integer;
		p->flag=-1;
		Polynomial poly(data.size()+a.data.size()-1);
		poly.get(data,a.data);
		p->data=poly.add();
		p->check();
		return (*p);
	}else if(flag==1&&a.flag==-1){
		Integer b=a;
		b.flag=1;
		return (*this)-b;
	}else if(flag==-1&&a.flag==1){
		Integer b=(*this);
	   b.flag=1;
	   return a-b;
	}
}

Integer& Integer:: operator-(Integer& a){
	if(flag==1&&a.flag==1){
		Integer*p;
 		p=new Integer;
 		Polynomial poly(data.size()+a.data.size()-1);
		p->flag=poly.cmp(this->data,a.data);
		if(p->flag==1){
			poly.get(data,a.data);
			p->data=poly.sub();
			return (*p);
		}else{
			poly.get(a.data,data);
			p->data=poly.sub();
			return (*p);
		}
	}else if(flag==-1&&a.flag==-1){
		Integer*p;
 		p=new Integer;
 		Polynomial poly(data.size()+a.data.size()-1);
		p->flag=poly.cmp(this->data,a.data);
		if(p->flag==1){
			poly.get(data,a.data);
			p->data=poly.sub();
			p->flag*=-1;
			return (*p);
		}else{
			poly.get(a.data,data);
			p->data=poly.sub();
			p->flag*=-1;
			return (*p);
		}
	}else{
		Integer b=a;
		b.flag*=-1;
		return (*this)+b;
	}
	
}

Integer& Integer:: operator-(){
	Integer*p;
 	p=new Integer;
 	p->data.clear();
	p->flag*=-1;
	p->data=data;
	return *p;
}

std::vector<int>& Integer::Polymul_tool(std::vector<int> a,std::vector<int> b){
	
	int deg=a.size()+b.size()-1;
	
	for(int i=a.size()-1;i>=0;i--)if(a[i]==0)a.pop_back();else break;
	for(int i=b.size()-1;i>=0;i--)if(b[i]==0)b.pop_back();else break;
	
	Polynomial test(a.size()+b.size());
	test.get(a,b);
	std::vector<int>* c;
	c=&test.mul();
	
	return *c;
}

Integer& Integer:: operator*(Integer& a){
		Integer*p;
 		p=new Integer;
 		p->data.clear();
		p->flag=this->flag*a.flag;
		
		p->data=Polymul_tool(this->data,a.data);
		
	return *p;
}

Integer& Integer:: operator/(Integer& a){
	
	
	Integer*p;
 	p=new Integer;
 	p->data.clear();
	p->flag=this->flag*a.flag;
	
	
	Integer tool1;
	Integer tool2;
	Integer tool3(1);
	tool1.data=this->data;
	tool2.data=a.data; 
	
	reverse(tool1.data.begin(),tool1.data.end());
	reverse(tool2.data.begin(),tool2.data.end());
	
	std::vector<int> z;
	int tes;
	
	p->data=Polydiv::div(tool1.data,tool2.data,z,&tes);
	reverse(p->data.begin(),p->data.end());
	reverse(tool1.data.begin(),tool1.data.end());
	reverse(tool2.data.begin(),tool2.data.end());
	
	if((*p+tool3)*tool2<=tool1)*p=*p+tool3;
	return (*p);
}

Integer& Integer:: operator%(Integer& a){
	return (*this)-(*this)/a*a;
}

std::istream& operator>>(std::istream& in,Integer&a){
	a.data.clear();
	std::string date;
	std::cin>>date;
	if(date[0]=='-'){
		a.flag=-1;
		for(int i=date.size()-1;i>=1;i--){
			a.data.push_back(date[i]-48);
		}
	}else{
		a.flag=1;
		for(int i=date.size()-1;i>=0;i--){
			a.data.push_back(date[i]-48);
		}
			
	} 
	a.check();
	return in;
}

std::ostream& operator<<(std::ostream& out,Integer&a){
	a.check();
	if(a.data.size()==0)out<<'0';
	else if(a.flag==1){
		for(int i=a.data.size()-1;i>=0;i--)out<<a.data[i];
	}else {
		out<<"-";
		for(int i=a.data.size()-1;i>=0;i--)out<<a.data[i];
	}
	return out;
}


//Polydiv类的函数 
std::vector<int>& Polydiv::Polymul_tool(std::vector<int> a,std::vector<int> b){
	
	int deg=a.size()+b.size()-1;
	reverse(a.begin(),a.end());
	reverse(b.begin(),b.end());
	
	for(int i=a.size()-1;i>=0;i--)if(a[i]==0)a.pop_back();else break;
	for(int i=b.size()-1;i>=0;i--)if(b[i]==0)b.pop_back();else break;
	
	
	Polynomial test(a.size()+b.size());
	test.get(a,b);
	std::vector<int>* c;
	c=&test.mul();
	
	while(c->size()<deg)(*c).push_back(0);
	reverse(c->begin(),c->end());
	return *c;
}


void Polydiv::print_vector(std::vector<int>&a){
	for(int i=0;i<a.size();i++)std::cout<<a[i];std::cout<<std::endl;
}

int Polydiv::ceilog(int x){int ans=0;while(1<<ans<x)ans++;return ans;}

std::vector<int>& Polydiv::Inverse(std::vector<int>& A){
	std::vector<int>* B,C;
	B=new std::vector<int>;
	B->resize(2);
	(*B)[1]=100/(A[0]*10+(A.size()>1?A[1]:0));
	int x=ceilog(A.size())+1;
	
	
	
	for(register int s=1;s<=10;++s){//迭代次数 
		C.resize(1<<s),B->resize(1<<s);
		for(int i=0;i<std::min(1ll<<s,(long long)A.size());++i) C[i]=A[i];
		for(int i=std::min(1ll<<s,(long long)A.size());i<1<<s;++i) C[i]=0;
		
		C=Polymul_tool((*B),C);
		
		for(register int i=1;i<C.size();++i) C[i]=-C[i];
		C[0]=2-C[0];
		for(register int i=C.size()-1;i;--i){
			C[i-1]+=(C[i]-9)/10;
			C[i]=(C[i]+10)%10;
		}
		(*B)=Polymul_tool((*B),C);
		B->resize(1<<s);
	}
	return (*B);
}

std::vector<int>& Polydiv::div(std::vector<int> A,std::vector<int> B,std::vector<int>& Z,int* tes){
	int rest=0;
	std::vector<int>* C;
	C=new std::vector<int>;
	rest+=A.size();
	rest-=B.size();
	std::vector<int> h=Inverse(B);
	A=Polymul_tool(A,h);
	Z=A;
	(*tes)=rest;
	for(int i=0;i<=rest;++i) if(!(i==0&&A[i]==0)) (*C).push_back(A[i]);
	
	return (*C);
}


//Polynomial的函数 
void Polynomial::get(std::vector<int>&A,std::vector<int>&B){ //A B都是规范化之后的 
		N=A.size();
		for(int i=A.size()-1;i>=0;i--){
			A1[i]=A[i];
		}
		M=B.size();
		for(int i=B.size()-1;i>=0;i--){
			B1[i]=B[i];
		}
}

int Polynomial::cmp(std::vector<int>&A,std::vector<int>&B){
		if(A.size()>B.size())return 1;
		else if(A.size()<B.size())return -1;
		else{
			for(int i=A.size()-1;i>=0;i--)
			{
				if(A[i]>B[i])return 1;else if(A[i]<B[i])return -1;
			}
			return 0;
		}
		return 0;
}

std::vector<int>& Polynomial::mul(){
	std::vector<int>* P;
	P=new std::vector<int>;
	if(N==0||M==0){
		P->push_back(0);
		return (*P); 
	}
	
	PolyMul(A1,B1,N+M-2);
	long long tool=0;
	for(int i=0;i<N+M-1;i++){
		long long a=tool+(A1[i]+Mod)%Mod;
		P->push_back(a%10);
		tool=a/10;
	}
	while(tool){
		P->push_back(tool%10);
		tool/=10;
	}
	for(int i=P->size()-1;i>=1;i--){
		if((*P)[i]==0)P->pop_back();else break;
	}
	Clear();
	return (*P); 
}

std::vector<int>& Polynomial::sub(){
		std::vector<int>* P;
		P=new std::vector<int>;
		if(N==0&&M==0){
			P->push_back(0);
			return (*P); 
		}
		int tool1=0;
		int Z=std::max(N,M);
		for(int i=0;i<Z;i++){
			int tool2=A1[i]-B1[i]+tool1;
			if(tool2>=0){
				tool1=0;
				P->push_back(tool2);
			}else{
				tool1=-1;
				P->push_back(10+tool2);
			}
		}
		for(int i=P->size()-1;i>=1;i--){
			if((*P)[i]==0)P->pop_back();else break;
		}
		Clear();
		return (*P);
}

std::vector<int>& Polynomial::add(){
		std::vector<int>* P;
		P=new std::vector<int>;
		if(N==0&&M==0){
			P->push_back(0);
			return (*P); 
		}
		int tool1=0;
		int Z=std::max(N,M);
		for(int i=0;i<Z;i++){
			int tool2=A1[i]+B1[i]+tool1;
			P->push_back(tool2%10);
			tool1=tool2/10;
		}
		
		if(tool1>0){
			P->push_back(tool1);
		}
		for(int i=P->size()-1;i>=1;i--){
			if((*P)[i]==0)P->pop_back();else break;
		}
		Clear();
		return (*P);
}

void Polynomial::Clear(){
	for(int i=0;i<MS;i++){
		R[i]=0;
		A1[i]=0;
		B1[i]=0;
	}
}


long long Polynomial::qPow(long long b, int e) {   //快速幂 
	long long a = 1;
	for (; e; e >>= 1, b = b * b % Mod)
		if (e & 1) a = a * b % Mod;
	return a;
}


void Polynomial::Init(int N) { //求逆元   
	Inv[1] = 1;
	for (int i = 2; i < N; ++i) Inv[i] = -(Mod / i) * Inv[Mod % i] % Mod;
	return ;
}

void Polynomial::InitFNTT(int N) { 
	int Bt = 0;
	for (; 1 << Bt <= N; ++Bt) ;
	if (Sz == (1 << Bt)) return ;
	Sz = 1 << Bt; InvSz = -(Mod - 1) / Sz;
	for (int i = 1; i < Sz; ++i) R[i] = R[i >> 1] >> 1 | (i & 1) << (Bt - 1);
} 

void Polynomial::FNTT(long long *A, int Ty) {  
	for (int i = 0; i < Sz; ++i) if (R[i] < i) std::swap(A[R[i]], A[i]);
	for (int j = 1, j2 = 2; j < Sz; j <<= 1, j2 <<= 1) {
		long long gn = qPow(~Ty ? G : iG, (Mod - 1) / j2), g, X, Y;
		for (int i = 0, k; i < Sz; i += j2) {
			for (k = 0, g = 1; k < j; ++k, g = g * gn % Mod) {
				X = A[i + k], Y = g * A[i + j + k] % Mod;
				A[i + k] = (X + Y) % Mod, A[i + j + k] = (X - Y) % Mod;
			}
		}
	}
	if (!~Ty) for (int i = 0; i < Sz; ++i) A[i] = A[i] * InvSz % Mod; 
}

void Polynomial::PolyMul(long long *A,long long *B, int deg){   //多项式乘法  A=A*B 
	InitFNTT(deg);
	FNTT(A,1);
	FNTT(B,1);
	for(int i=0;i<Sz;i++)
    {
        A[i]=(A[i]*B[i])%Mod;
    }
	FNTT(A, -1);
	FNTT(B, -1);
}

void Polynomial::PolyInv(long long *A, int N, long long *B) {    //多项式求inv 对x^N取mod  A*B=1
	long long tA[MS], tB[MS];
	B[0] = qPow(A[0], Mod - 2);
	for (int L = 1; L < N; L <<= 1) {
		int L2 = L << 1, L4 = L << 2;
		InitFNTT(L4);
		memcpy(tA, A, 8 * L2);
		memset(tA + L2, 0, 8 * (Sz - L2));
		memcpy(tB, B, 8 * L);
		memset(tB + L, 0, 8 * (Sz - L));
		FNTT(tA, 1), FNTT(tB, 1);
		for (int i = 0; i < Sz; ++i) tB[i] = (2 - tB[i] * tA[i]) % Mod * tB[i] % Mod;
		FNTT(tB, -1); 
		for (int i = 0; i < L2; ++i) B[i] = tB[i];
	}
}


void Polynomial::PolyLn(long long *A, int N, long long *B) {  // 多项式求ln 对x^N取mod A[0]=1 B=ln(A) 
	 long long tA[MS], tB[MS];
	PolyInv(A, N - 1, tB);
	InitFNTT(N * 2 - 3);
	for (int i = 1; i < N; ++i) tA[i - 1] = i * A[i] % Mod;
	memset(tA + N - 1, 0, 8 * (Sz - N + 1));
	memset(tB + N - 1, 0, 8 * (Sz - N + 1));
	FNTT(tA, 1), FNTT(tB, 1);
	for (int i = 0; i < Sz; ++i) tA[i] = (long long)tA[i] * tB[i] % Mod; 
	FNTT(tA, -1);
	B[0] = 0;	
	
	for (int i = 1; i < N; ++i) B[i] = (long long)tA[i - 1] * Inv[i] % Mod;

}


void Polynomial::PolyExp(long long *A, int N, long long *B) {  // 多项式求exp 对x^N取mod A[0]=0  B=e^A
	long long tA[MS], tB[MS];
	B[0] = 1;
	for (int L = 1; L < N; L <<= 1) {
		int L2 = L << 1, L4 = L << 2;
		memset(B + L, 0, 8 * (L2 - L));
		PolyLn(B, L2, tB);
		InitFNTT(L4);
		memcpy(tA, B, 8 * L);
		memset(tA + L, 0, 8 * (Sz - L));
		for (int i = 0; i < L2; ++i) tB[i] = ((!i) - tB[i] + A[i]) % Mod;
		memset(tB + L2, 0, 8 * (Sz - L2));
		FNTT(tA, 1), FNTT(tB, 1);
		for (int i = 0; i < Sz; ++i) tA[i] = tA[i] * tB[i] % Mod;
		FNTT(tA, -1);
		for (int i = 0; i < L2; ++i) B[i] = tA[i];
	}
}
	
