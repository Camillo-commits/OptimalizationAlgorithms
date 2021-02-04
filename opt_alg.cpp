#include"opt_alg.h"
#include"matrix.h"
#include<math.h>
#include<fstream>
#if LAB_NO>1
double *expansion(double x0, double d, double alfa, int Nmax, matrix O)
{
	double *p = new double[2];
	//create matrix with 2 doubles
	solution X0(x0), X1(matrix(x0 + d));
	X0.fit_fun(matrix(x0));
	X1.fit_fun(matrix(x0 + d));
	
	if (X0.y == X1.y)
	{
		p[0] = x0 ;
		p[1] = x0 + d ;
		return p;
	}
	if (X0.y < X1.y )
	{
		d *= -1;
		X1.x = X0.x + d;
		X1.fit_fun(X1.x);
		if (X0.x <= X1.x )
		{
			p[0] = X1.x(0) ;
			p[1] = X0.x(0) ;
			return p;
		}
	}
	solution X2;
	int i = 1;
	while (true)
	{
		X2.x = x0 + pow(alfa, i) * d;
		X2.fit_fun() ;
		//X2.fit_fun();
		if (X2.y >= X1.y || solution::f_calls > Nmax )
			break;
		X0 = X1 ;
		X1 = X2 ;
		++i;
	}
	if (d > 0) {
		p[0] = X0.x(0);
		p[1] = X2.x(0);
	}
	else {
		p[1] = X0.x(0);
		p[0] = X2.x(0);
	
	}

	return p;
}

double logBase(double a, double base) {
	return log(a) / log(base);
}

solution fib(double a, double b, double epsilon, matrix O)
{
	int n = log(sqrt(5) * 100) / log((1 + sqrt(5)) / 2) + 1;
	double *F = new double[n] {1, 1};
	for (int i = 0; i < n; ++i)
		F[i] = F[i - 2] + F[i - 1];

	solution A(a), B(b), C, D;
	C.x = B.x(0, 0) - (F[n - 2] / F[n - 1]) * (B.x(0, 0) - A.x(0, 0));
	D.x = A.x + B.x - C.x;
	C.fit_fun();
	D.fit_fun();
	ofstream fibB_A("fibB-A.txt");
	for (int i = 0; i <= n - 3; ++i) {
		if (C.y(0, 0) < D.y(0, 0))
			B.x = D.x;
		else
			A.x = C.x;
		double tmp = A.x(0, 0) - B.x(0, 0);
		fibB_A << tmp << "\n";
		C.x = B.x - (F[n - i - 3] / F[n - i - 2]) * (B.x - A.x);
		D.x = A.x + B.x - C.x;
		C.fit_fun();
		D.fit_fun();
	}
	return C;
	/*
	int n;
	double tmp = logBase(sqrt(5) * ((b-a)/epsilon) ,((1 + sqrt(5))/2));
	n = tmp * 10;
	if (n % 10 > 0) {
		n = tmp;
		++n;
	}
	else {
		n = tmp;
	}
	int *F = new int[n];
	F[0] = 1;
	F[1] = 1;
	for (int i = 2; i < n; ++i)
		F[i] = F[i - 2] + F[i - 1];
	solution A(a), B(b);
	solution C, D;
	
	double aX = A.x(0);
	double bX = B.x(0);
	
	C.x = bX - ((F[n-2]/ F[n-1]) * (bX - aX));
	double cX = C.x(0);
	D.x = aX + bX - cX;
	//D.x = A.x(0) + B.x(0) - C.x(0);
	C.fit_fun();
	D.fit_fun();
	ofstream fibB_A("fibB-A.txt");
	for (int i = 0; i <= n - 3; ++i)
	{
		if (D.y > C.y )
			B.x = D.x(0) ;
		else
			A.x = C.x(0) ;
		aX = A.x(0);
		bX = B.x(0);
		C.x = bX - ((F[n - 1] / F[n - 2]) * (bX - aX));
		D.x = A.x(0) + B.x(0) - C.x(0);
		C.fit_fun(matrix(C.x(0)));
		D.fit_fun(matrix(D.x(0)));
		
		//saving b-a
		double streamTmp1 = aX;
		double streamTmp2 = bX;
		streamTmp2 -= streamTmp1;
		fibB_A << streamTmp2 << "\n";
	}
	fibB_A.close();
	return C;
	*/
}

solution lag(double a, double b, double epsilon, double gamma, int Nmax, matrix O)
{
	solution A(a), B(b), C((a + b) / 2), D(0);
	
	A.fit_fun(matrix(A.x(0)));
	B.fit_fun(matrix(B.x(0)));
	C.fit_fun(matrix(C.x(0)));
	double l, m;

	ofstream lagrB_A("lagrB-A.txt");
	double tmpStr = A.x(0,0) - B.x(0,0);

	while (true)
	{
		l = A.y(0)*(pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y(0)*(pow(C.x(0), 2) - pow(A.x(0), 2)) + C.y(0)*(pow(A.x(0), 2) - pow(B.x(0), 2));
		m = A.y(0)*(B.x(0) - C.x(0)) + B.y(0)*(C.x(0) - A.x(0)) + C.y(0)*(A.x(0) - B.x(0));
		if (m <= 0)
		{
			C.x = NAN;
			C.y = NAN;

			tmpStr = A.x(0, 0) - B.x(0, 0);
			lagrB_A << tmpStr << "\n";
			return C;
		}
		D.x = 0.5 * l/m ;
		D.fit_fun(matrix(D.x(0)));
		if (A.x(0) < C.x(0) && C.x(0) < D.x(0) )
		{
			//1
			if (C.y(0) > D.y(0) )
			{
				A = C ;
				C = D ;
			}//2
			else
				B = D ;

			tmpStr = A.x(0, 0) - B.x(0, 0);
			lagrB_A << tmpStr << "\n";
		}
		else if (A.x(0) < D.x(0) && D.x(0) < C.x(0) )
		{
			//4
			if (D.y(0) > C.y(0) )
			{
				B.x = C.x(0) ;
				C.x = D.x(0);
			}
			//3
			else
				A.x = D.x(0);

			tmpStr = A.x(0, 0) - B.x(0, 0);
			lagrB_A << tmpStr << "\n";
		}
		else
		{
			C.x = NAN;
			C.y = NAN;
			
			tmpStr = A.x(0, 0) - B.x(0, 0);
			lagrB_A << tmpStr << "\n";
			
			return C;
		}
		if (Nmax <= D.f_calls - 3 || fabs(D.x(0) - C.x(0)) < gamma || A.x(0) - B.x(0) < epsilon) {
			C.fit_fun(matrix(C.x(0)));
			
			tmpStr = A.x(0, 0) - B.x(0, 0);
			lagrB_A << tmpStr << "\n";
			
			return C;
		}
			
	}
}
#endif
#if LAB_NO>2
solution HJ(matrix x0, double s, double alfa, double epsilon, int Nmax, matrix O)
{
	
	double *arr = new double[2];
	arr[0] = x0(0, 0);
	arr[1] = x0(0, 1);

	solution XB(arr,2), XB_old, X;
	//solution XB, XB_old, X;
	//XB.x = x0;
	//XB.x(0, 0) = x0(0, 0);
	//XB.x(0, 1) = x0(0, 1);
	XB.fit_fun();
	int counter = 0;
	while (true)
	{
		X = HJ_trial(XB , s);
		if (X.y(0,0) < XB.y(0,0) )
		{
			while (true)
			{
				XB_old = XB;
				XB = X ;
				X.x = XB.x + (XB.x - XB_old.x);
				X.fit_fun();
				X = HJ_trial(X, s);
				//
				if (X.y(0,0) >= XB.y(0,0))
					break;
				if (counter >= Nmax)
					return XB;
				counter++;
			}
		}
		else
			s *= alfa ;
		counter++;
		if (counter >= Nmax || s < epsilon )
			return XB;
	}
}

matrix identMat(int size) {
	int nv = size;
	matrix A(nv, nv);
	for (int i = 0; i < nv; ++i)
		A(i, i)= 1;
	return A;
}

solution HJ_trial(solution XB, double s, matrix O)
{
	int *n = get_size(XB.x);
	matrix D = identMat(n[0]);
	double* arr = new double[*n];
	for (int i = 0; i < *n; ++i) {
		arr[i] = 0;
	}
	solution X(arr,*n);
	
	for (int i = 0; i < n[0]; ++i)
	{
		X.x = XB.x + s;
		X.fit_fun();
		if (X.y(0,0) < XB.y(0,0))
			XB.x(i,0) = XB.x(i,0) + s;
		else
		{
			X.x(i,0) = XB.x(i,0) - s;
			X.fit_fun();
			if (X.y(0, 0) < XB.y(0, 0))
				XB.x = XB.x - s;
		}
		//na wszelki wypadek bo zmieniamy tylko X
		XB.fit_fun();
	}


	return XB;
}


solution Rosen(matrix x0, matrix s0, double alfa, double beta, double epsilon, int Nmax, matrix O)
{
	int *n = get_size(x0);
	matrix l(n[0], 1), p(n[0], 1), s(s0);
	matrix D(n[0], n[0]);
	for (int i = 0; i < n[0]; i++)
	{
		D(i, i) = 1;
	}
	solution X, Xt;
	X.x = x0;
	X.fit_fun();
	while (true)
	{
		for (int i = 0; i < n[0]; ++i)
		{
			Xt.x = X.x + s(i) * D[i];
			Xt.fit_fun();
			if (Xt.y < X.y)
			{
				X = Xt;
				l(i) += s(i);
				s(i) *= alfa;
			}
			else
			{
				p(i) += 1;
				s(i) *= -beta;
			}
		}
		bool change = true;
		for (int i = 0; i < n[0]; ++i)
			if (l(i) == 0 || p(i) == 0)
			{
				change = false;
				break;
			}
		if (change)
		{
			matrix Q(n[0], n[0]), v(n[0], 1);
			for (int i = 0; i < n[0]; ++i)
				for (int j = 0; j <= i; ++j)
					Q(i, j) = l(i);								//macierz lambd
			Q = D * Q;
			v = Q[0] / norm(Q[0]);								
			D = set_col(D, v, 0);
			for (int i = 1; i < n[0]; ++i)
			{
				matrix temp(n[0], 1);
				for (int j = 0; j < i; ++j)
					temp = (temp + trans(Q[i]) * D[j] * D[j]);
				v = (Q[i] - temp) / norm(Q[i] - temp);
				D = set_col(D, v, i);
			}
			s = s0;
			l = matrix(n[0], 1);
			p = matrix(n[0], 1);
		}
		double max_s = abs(s(0));
		for (int i = 1; i < n[0]; ++i)
			if (max_s < abs(s(i)))
				max_s = abs(s(i));
		if (max_s < epsilon || solution::f_calls > Nmax)
			return X;
		//solution X;
	}
}

#endif
#if LAB_NO>3
//metoda funkcji kary

solution sym_NM_outside(matrix x0, double s, double alfa, double beta, double gama, double delta, double epsilon, double Nmax, matrix O)
{
	int* n = get_size(x0);
	int N = n[0] + 1;
	matrix D(n[0], n[0]);
	for (int i = 0; i < n[0]; ++i)
		D(i, i) = 1;
	solution* S = new solution[N];
	S[0].x = x0;
	S[0].fit_fun_outside(O);
	for (int i = 1; i < N; ++i) {
		S[i].x = S[0].x + s * D[i - 1];
		S[i].fit_fun_outside(O);
	}
	solution p_o, p_e, p_z;
	matrix p_sr;
	int i_min, i_max;

	while (true) {
		i_min = i_max = 0;
		for (int i = 1; i < N; ++i) {
			if (S[i_min].y > S[i].y)
				i_min = i;
			if (S[i_max].y < S[i].y)
				i_max = i;
		}
		p_sr = matrix(n[0], 1);
		for (int i = 0; i < N; ++i) {
			if (i != i_max)
				p_sr = p_sr + S[i].x;
		}
		p_sr = p_sr / (N - 1);
		p_o.x = p_sr + alfa * (p_sr - S[i_max].x);
		p_o.fit_fun_outside(O);
		if (p_o.y >= S[i_min].y && p_o.y < S[i_max].y)
			S[i_max] = p_o;
		else if (p_o.y < S[i_min].y)
		{
			p_e.x = p_sr + gama * (p_o.x - p_sr);
			p_e.fit_fun_outside(O);
			if (p_e.y < p_o.y)
				S[i_max] = p_e;
			else
				S[i_max] = p_o;
		}
		else
		{
			p_z.x = p_sr + beta * (S[i_max].x - p_sr);
			p_z.fit_fun_outside(O);
			if (p_z.y < S[i_max].y)
				S[i_max] = p_z;
			else
			{
				for (int i = 0; i < N; ++i) {
					if (i != i_min)
					{
						S[i].x = delta * (S[i].x + S[i_min].x);
						S[i].fit_fun_outside(O);
					}
				}
			}
		}
		double max_s = norm(S[0].x - S[i_min].x);
		for (int i = 1; i < N; ++i) {
			if (max_s < norm(S[i].x - S[i_min].x))
				max_s = norm(S[i].x - S[i_min].x);
		}
		if (max_s<epsilon || solution::f_calls>Nmax)
			return S[i_min];
	}
}

solution pen_outside(matrix x0, double c, double a, double epsilon, int Nmax)
{
	double alfa = 1, beta = 0.5, gama = 2, delta = 0.5, s = 0.5;
	matrix A(new double[2]{ c,a }, 2);
	solution X, X1;
	X.x = x0;
	while (true)
	{
		X1 = sym_NM_outside(X.x, s, alfa, beta, gama, delta, epsilon, Nmax, A);
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax)
			return X1;
		A(0) *= 2;		//dc=2
		X = X1;
	}
}

solution sym_NM_inside(matrix x0, double s, double alfa, double beta, double gama, double delta, double epsilon, double Nmax, matrix O)
{
	int* n = get_size(x0);
	int N = n[0] + 1;
	matrix D(n[0], n[0]);
	for (int i = 0; i < n[0]; ++i)
		D(i, i) = 1;
	solution* S = new solution[N];
	S[0].x = x0;
	S[0].fit_fun_inside(O);
	for (int i = 1; i < N; ++i) {
		S[i].x = S[0].x + s * D[i - 1];
		S[i].fit_fun_inside(O);
	}
	solution p_o, p_e, p_z;
	matrix p_sr;
	int i_min, i_max;

	while (true) {
		i_min = i_max = 0;
		for (int i = 1; i < N; ++i) {
			if (S[i_min].y > S[i].y)
				i_min = i;
			if (S[i_max].y < S[i].y)
				i_max = i;
		}
		p_sr = matrix(n[0], 1);
		for (int i = 0; i < N; ++i) {
			if (i != i_max)
				p_sr = p_sr + S[i].x;
		}
		p_sr = p_sr / (N - 1);
		p_o.x = p_sr + alfa * (p_sr - S[i_max].x);
		p_o.fit_fun_inside(O);
		if (p_o.y >= S[i_min].y && p_o.y < S[i_max].y)
			S[i_max] = p_o;
		else if (p_o.y < S[i_min].y)
		{
			p_e.x = p_sr + gama * (p_o.x - p_sr);
			p_e.fit_fun_inside(O);
			if (p_e.y < p_o.y)
				S[i_max] = p_e;
			else
				S[i_max] = p_o;
		}
		else
		{
			p_z.x = p_sr + beta * (S[i_max].x - p_sr);
			p_z.fit_fun_inside(O);
			if (p_z.y < S[i_max].y)
				S[i_max] = p_z;
			else
			{
				for (int i = 0; i < N; ++i) {
					if (i != i_min)
					{
						S[i].x = delta * (S[i].x + S[i_min].x);
						S[i].fit_fun_inside(O);
					}
				}
			}
		}
		double max_s = norm(S[0].x - S[i_min].x);
		for (int i = 1; i < N; ++i) {
			if (max_s < norm(S[i].x - S[i_min].x))
				max_s = norm(S[i].x - S[i_min].x);
		}
		if (max_s<epsilon || solution::f_calls>Nmax)
			return S[i_min];
	}
}

solution pen_inside(matrix x0, double c, double a, double epsilon, int Nmax)
{
	double alfa = 1, beta = 0.5, gama = 2, delta = 0.5, s = 0.5;
	matrix A(new double[2]{ c,a }, 2);
	solution X, X1;
	X.x = x0;
	while (true)
	{
		X1 = sym_NM_inside(X.x, s, alfa, beta, gama, delta, epsilon, Nmax, A);
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax)
			return X1;
		A(0) *= 0.1;	//dc=0,1
		X = X1;
	}
}




solution sym_NM_zew(matrix x0, double s, double alfa, double beta, double gama, double delta, double epsilon, int Nmax, matrix O);
solution pen(matrix x0, double c0, double dc, double epsilon, int Nmax, matrix O, int which)		//c0 -> pocz¹tkowa si³a funkcji kary
{																						//dc -> skok si³y funkcji kary
	double alfa = 1, beta = 0.5, gama = 2, delta = 0.5, s = 0.5;
	matrix A(new double[2]{ c0,O(0) }, 2);
	solution X, X1;
	X.x = x0;
	while (true)
	{
		if (which == 0) {
			X1 = sym_NM_zew(X.x, s, alfa, beta, gama, delta, epsilon, Nmax, A);
			A(0) *= 2;
		}
		else {
			X1 = sym_NM(X.x, s, alfa, beta, gama, delta, epsilon, Nmax, A);
			A(0) *= 0.1;
		}
		if (solution::f_calls > Nmax || norm(X.x - X1.x) < epsilon)
			return X1;
		
		X = X1;
	}
}
 
solution sym_NM_zew(matrix x0, double s, double alfa, double beta, double gama, double delta, double epsilon, int Nmax, matrix O) {
	solution::clear_calls();
	int *n = get_size(x0);
	matrix D = identMat(n[0]);
	int N = n[0] + 1;
	solution *S = new solution[N];
	S[0].x = x0;
	S[0].fit_fun_zew(O);
	for (int i = 1; i < N; ++i)
	{				//obliczanie kolejnych czêœci simpleksu
		S[i].x = S[0].x + s * D[i - 1];
		S[i].fit_fun_zew(O);
	}
	solution p_o, p_e, p_z;									//	p_o -> po odbiciu		 p_e -> po ekspansji		p_z -> po zawê¿eniu
	matrix p_sr;											//	p_sr -> punkt œrodkowy
	int i_min, i_max;
	while (true)
	{
		i_min = i_max = 0;
		for (int i = 1; i < N; ++i)							//szukamy minimalnych i maksymalnych wartoœci y w simpleksie
		{
			if (S[i].y < S[i_min].y)
				i_min = i;
			if (S[i].y > S[i_max].y)
				i_max = i;
		}
		p_sr = matrix(n[0], 1);
		for (int i = 0; i < N; ++i) {							// znajdywanie œrodka ciê¿koœci ( œrednia arytm. wszystkich y bez maksymalnego )
			if (i != i_max)
				p_sr = p_sr + S[i].x;
		}
		p_sr = p_sr / (N - 1);
		p_o.x = p_sr + alfa * (p_sr - S[i_max].x);			//1. Pierwsza operacja geometryczna (Odbicie)
		p_o.fit_fun_zew(O);										//liczenie po³o¿enia punktu odbitego
		//if (S[i_min].y <= p_o.y && p_o.y < S[i_max].y)		//		a) akceptacja odbicia	
		if (p_o.y >= S[i_min].y && p_o.y < S[i_max].y)
			S[i_max] = p_o;
		else if (p_o.y < S[i_min].y)						//		b) wykonujemy ekspansje
															//2. Ekspansja
		{
			p_e.x = p_sr + gama * (p_o.x - p_sr);
			p_e.fit_fun_zew(O);
			if (p_e.y < p_o.y)
				S[i_max] = p_e;								//krok ekspansji udany		-> nowy simpleks sk³ada siê z p_e, p_min i p, które nie by³o ani min ani max
			else
				S[i_max] = p_o;								//krok ekspansji nie udany	-> przyjmujemy simpleks odbity
		}
		else												//3. Zawê¿enie
		{
			p_z.x = p_sr + beta * (S[i_max].x - p_sr);
			p_z.fit_fun_zew(O);
			if (p_z.y < S[i_max].y)							//		a) akceptacja zawê¿enia	
				S[i_max] = p_z;
			else											//		b) wykonujemy redukcjê
															//4. Redukcja
			{
				for (int i = 0; i < N; ++i)
					if (i != i_min)
					{
						S[i].x = delta * (S[i].x + S[i_min].x);
						S[i].fit_fun_zew(O);
					}
			}
		}
		double max_s = norm(S[0].x - S[i_min].x);
		for (int i = 1; i < N; ++i)							//Warunek stopu NR 1 -> gdy odleg³oœæ ka¿dego wierzcho³ka simpleksu od najlepszego (minimalnego) jest mniejsza ni¿ epsilon
			if (max_s < norm(S[i].x - S[i_min].x))
				max_s = norm(S[i].x - S[i_min].x);
		if (max_s < epsilon || solution::f_calls > Nmax)	//Warunek stopu NR 2 -> gdy przekroczymy iloœæ iteracji
			return S[i_min];
	}
}

solution sym_NM(matrix x0, double s, double alfa, double beta, double gama, double delta, double epsilon, int Nmax, matrix O)
{
	solution::clear_calls();
	int *n = get_size(x0);
	matrix D = identMat(n[0]);
	int N = n[0] + 1; 
	solution *S = new solution[N];
	S[0].x = x0  ;
	S[0].fit_fun(O); 
	for (int i = 1; i < N; ++i)
	{				//obliczanie kolejnych czêœci simpleksu
		S[i].x = S[0].x + s * D[i - 1];
		S[i].fit_fun(O);
	}
	solution p_o, p_e, p_z;									//	p_o -> po odbiciu		 p_e -> po ekspansji		p_z -> po zawê¿eniu
	matrix p_sr;											//	p_sr -> punkt œrodkowy
	int i_min, i_max;
	while (true)
	{
		i_min = i_max = 0;
		for (int i = 1; i < N; ++i)							//szukamy minimalnych i maksymalnych wartoœci y w simpleksie
		{
			if (S[i].y < S[i_min].y)
				i_min = i;
			if (S[i].y > S[i_max].y)
				i_max = i;
		}
		p_sr = matrix(n[0], 1);
		for (int i = 0; i < N; ++i) {							// znajdywanie œrodka ciê¿koœci ( œrednia arytm. wszystkich y bez maksymalnego )
			if (i != i_max)
				p_sr = p_sr + S[i].x;
		}
		p_sr = p_sr / (N - 1);
		p_o.x = p_sr + alfa * (p_sr - S[i_max].x);			//1. Pierwsza operacja geometryczna (Odbicie)
		p_o.fit_fun(O);										//liczenie po³o¿enia punktu odbitego
		//if (S[i_min].y <= p_o.y && p_o.y < S[i_max].y)		//		a) akceptacja odbicia	
		if(p_o.y >= S[i_min].y && p_o.y < S[i_max].y)
			S[i_max] = p_o;
		else if (p_o.y < S[i_min].y)						//		b) wykonujemy ekspansje
															//2. Ekspansja
		{
			p_e.x = p_sr  + gama * (p_o.x - p_sr);
			p_e.fit_fun(O);
			if (p_e.y < p_o.y)						
				S[i_max] = p_e;								//krok ekspansji udany		-> nowy simpleks sk³ada siê z p_e, p_min i p, które nie by³o ani min ani max
			else
				S[i_max] = p_o;								//krok ekspansji nie udany	-> przyjmujemy simpleks odbity
		}
		else												//3. Zawê¿enie
		{
			p_z.x = p_sr + beta * (S[i_max].x - p_sr);
			p_z.fit_fun(O);
			if (p_z.y < S[i_max].y)							//		a) akceptacja zawê¿enia	
				S[i_max] = p_z;
			else											//		b) wykonujemy redukcjê
															//4. Redukcja
			{
				for (int i = 0; i < N; ++i)
					if (i != i_min)
					{
						S[i].x = delta * (S[i].x + S[i_min].x);
						S[i].fit_fun(O);
					}
			}
		}
		double max_s = norm(S[0].x - S[i_min].x);
		for (int i = 1; i < N; ++i)							//Warunek stopu NR 1 -> gdy odleg³oœæ ka¿dego wierzcho³ka simpleksu od najlepszego (minimalnego) jest mniejsza ni¿ epsilon
			if (max_s < norm(S[i].x - S[i_min].x))
				max_s = norm(S[i].x - S[i_min].x);
		if (max_s < epsilon || solution::f_calls > Nmax)	//Warunek stopu NR 2 -> gdy przekroczymy iloœæ iteracji
			return S[i_min];
	}
}
#endif
#if LAB_NO>4
solution SD(matrix x0, double h0, double epsilon, int Nmax, matrix O)
{ 
	solution::clear_calls();

	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2), limits = O;
	solution h;
	double b;

	ofstream x1SD_005("xSD_005.txt");
	ofstream x1SD_012("xSD_012.txt");
	ofstream x1SD_zmK("xSD_zmK.txt");

	while (true)
	{
		X.grad();
		d = -X.g;
		P = set_col(P, X.x, 0);		//macierz p s³u¿y do policzenia przesuniêcia. Zawiera: | x1	d1 |
		P = set_col(P, d, 1);		//													   | x2 d2 |
		if (h0 < 0)					//metoda zmienno-krokowa
		{
			b = compute_b(X.x, d, limits);
			h = golden(0 , b , epsilon, Nmax, P);
			X1.x = X.x + h.x * d;
		}
		else
			X1.x = X.x + h0 * d;				//metoda sta³o-krokowa\
		//wprowadzanie do pliku
		if (h0 == 0.05) {
			x1SD_005 << X1.x(0) << ";" << X1.x(1) << std::endl;
		}

		if (h0 == 0.12) {
			x1SD_012 << X1.x(0) << ";" << X1.x(1) << std::endl;
		}
		if (h0 < 0) {
			x1SD_zmK << X1.x(0) << ";" << X1.x(1) << std::endl;
		}
		///////////////////
		if (norm(X.x - X1.x) < epsilon ||
			solution::f_calls >= Nmax  ||
			solution::g_calls >= Nmax )
		{
			X1.fit_fun();
			return X1;
		}
		X = X1;
	}
}

solution CG(matrix x0, double h0, double epsilon, int Nmax, matrix O)
{
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2), limits = O;
	solution h;
	double b, beta;
	X.grad();
	d = -X.g;

	ofstream x1SD_005("xCG_005.txt");
	ofstream x1SD_012("xCG_012.txt");
	ofstream x1SD_zmK("xCG_zmK.txt");

	while (true)
	{
		P = set_col(P, X.x, 0);
		P = set_col(P, d, 1);
		if (h0 < 0)
		{
			b = compute_b(X.x, d, limits);
			h = golden(0, b, epsilon, Nmax, P);
			X1.x = X.x + h.x * d;
		}
		else
			X1.x = X.x + h0 * d;
		//wprowadzanie do pliku
		if (h0 == 0.05) {
			x1SD_005 << X1.x(0) << ";" << X1.x(1) << std::endl;
		}

		if (h0 == 0.12) {
			x1SD_012 << X1.x(0) << ";" << X1.x(1) << std::endl;
		}
		if (h0 < 0) {
			x1SD_zmK << X1.x(0) << ";" << X1.x(1) << std::endl;
		}
		///////////////////
		if (norm(X.x - X1.x) < epsilon ||
			solution::f_calls >= Nmax ||
			solution::g_calls >= Nmax)
		{
			X1.fit_fun();
			return X1;
		}
		X1.grad();
		beta = pow(norm(X1.g),2) / pow(norm(X.g),2);		//wed³ug metody Fletchera - Reevesa
		d = -X1.g + beta * d;
		X = X1;
	}
}

solution Newton(matrix x0, double h0, double epsilon, int Nmax, matrix O)
{
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2), limits = O;
	solution h;
	double b;

	ofstream x1SD_005("xN_005.txt");
	ofstream x1SD_012("xN_012.txt");
	ofstream x1SD_zmK("xN_zmK.txt");

	while (true)
	{
		X.grad();
		X.hess();
		d = - inv(X.H) * X.g;
		P = set_col(P, X.x, 0);
		P = set_col(P, d, 1);
		if (h0 < 0)
		{
			b = compute_b(X.x , d , limits);
			h = golden(0 , b , epsilon, Nmax, P);
			X1.x = X.x + h.x * d;
		}
		else
			X1.x = X.x + h0 * d;
		//wprowadzanie do pliku
		if (h0 == 0.05) {
			x1SD_005 << X1.x(0) << ";" << X1.x(1) << std::endl;
		}

		if (h0 == 0.12) {
			x1SD_012 << X1.x(0) << ";" << X1.x(1) << std::endl;
		}
		if (h0 < 0) {
			x1SD_zmK << X1.x(0) << ";" << X1.x(1) << std::endl;
		}
		///////////////////
		if (norm(X.x - X1.x) < epsilon ||
			solution::f_calls >= Nmax ||
			solution::g_calls >= Nmax ||
			solution::H_calls >= Nmax)
		{
			X1.fit_fun();
			return X1;
		}
		X = X1;
	}
}

solution golden(double a, double b, double epsilon, int Nmax, matrix O)
{
	//solution::clear_calls();
	double alfa = 1 / 1.6180340;//(a + b) / a );
	solution A, B, C, D;
	A.x = a ;
	B.x = b ;
	C.x = B.x - alfa * (B.x - A.x);
	C.fit_fun(O);
	D.x = A.x + alfa * (B.x - A.x);
	D.fit_fun(O);
	while (true)
	{
		if (C.y < D.y)
		{
			B = D ;
			D = C ;
			C.x = B.x - alfa * (B.x - A.x);
			C.fit_fun(O);
		}
		else
		{
			A = C;
			C = D;
			D.x = A.x + alfa * (B.x - A.x);
			D.fit_fun(O);
		}
		if (B.x - A.x < epsilon  || solution::f_calls >= Nmax)
		{
			A.x = (A.x + B.x) / 2.0;
			A.fit_fun(O);
			return A;
		}
	}
}

double compute_b(matrix x, matrix d, matrix limits)
{
	int *n = get_size(x);
	double b = 1e9, bi;
	for (int i = 0; i < n[0]; ++i)
	{
		if (d(i) == 0)
			bi = 1e9;
		else if (d(i) > 0)
			bi = (limits(i,1) - x(i)) / d(i);
		else
			bi = (limits(i,0) - x(i)) / d(i);
		if (b > bi)
			b = bi;
	}
	return b;
}
#endif
#if LAB_NO>5
solution Powell(matrix x0, double epsilon, int Nmax, matrix O)
{
	
	int *n = get_size(x0);
	matrix L(n[0], n[0]);
	for (int i = 0; i < n[0]; i++) {
		L(i, i) = 1;
	}
	
	matrix D = L, A(n[0], 3), limits(n[0], 2);
	limits = set_col(limits, O[0], 0);
	limits = set_col(limits, O[1], 1);
	A(0, 2) = O(0, 2);
	solution X, P, h;
	X.x = x0;
	double *ab = new double[2];

	while (true)
	{
		P = X.x; 
		
		for (int i = 0; i < n[0]; ++i)
		{
			A = set_col(A, P.x, 0);    
			A = set_col(A, D[i], 1);    
			ab = compute_ab(ab[0], ab[1], limits);
			h = golden(ab[0], ab[1], epsilon, Nmax, A); 
			P.x = P.x + h.x*D[i];
		}

		
		if (norm(X.x - P.x) < epsilon || solution::f_calls > Nmax)
		{
			P.fit_fun();
			return P; 
		}
		for (int i = 0; i < n[0] - 1; ++i)
			D = set_col(D, D[i + 1], i);    
		D = set_col(D, P.x - X.x, n[0] - 1);   
		A = set_col(A, P.x, 0);
		A = set_col(A, D[n[0] - 1], 1);
		ab = compute_ab(ab[0], ab[1], limits);
		h = golden(ab[0], ab[1], epsilon, Nmax, A);
		X.x = P.x + h.x *D[n[0] - 1];  
	}
}
/*
solution Powell(matrix x0, double epsilon, int Nmax, matrix O)
{
	int *n = get_size(x0);	//ilosc wymiarow
	matrix D(n[0],n[0]) /*= identMat(n[0])*//*, A(n[0], 3), limits(n[0], 2);
	for (int i = 0; i < n[0]; ++i) {
		D(i, i) = 1;
	}
	//matrix o_tmp(O[1]);
	limits = set_col(limits, O[0], 0);
	limits = set_col(limits, O[1], 1);
	
	A(0, 2) = O(0, 2);
	solution X, P, h;
	X.x = x0;
	double *ab;		//tak jak b
	while (true)
	{
		P = X;
		for (int i = 0; i < n[1] ; ++i)
		{
			A = set_col(A, P.x, 0);
			A = set_col(A, D[i], 1);
			ab = compute_ab(X.x , A[1] , limits);
			h = golden(ab[0] , ab[1] , epsilon, Nmax, A);
			P.x = P.x + h.x * D[i];	//ok
			std::cout << "P.x = " << P.x << std::endl;
		}
		if (solution::f_calls >= Nmax || norm(X.x - P.x) < epsilon )
		{
			P.fit_fun();
			return P;
		}
		for (int i = 0; i < n[0] - 1; ++i)
			set_col(D, D[i+1], i);
		D = set_col(D, P.x - X.x, n[0] - 1);
		A = set_col(A, P.x, 0);
		A = set_col(A, D[n[0] - 1], 1);
	
		ab = compute_ab(X.x , A[1] , limits);
		h = golden(ab[0] , ab[1] , epsilon, Nmax, A);
		X.x = X.x + h.x * D[n[0] - 1];
		std::cout << "X.x = " << X.x << std::endl;
	}
}*/

double *compute_ab(matrix x, matrix d, matrix limits)	//ok
{
	int *n = get_size(x);
	double *ab = new double[2]{ -1e9,1e9 };
	double ai, bi;
	for (int i = 0; i < n[0]; ++i)
	{
		if (d(i) == 0)
		{
			ai = -1e9 ;
			bi = 1e9 ;
		}
		//ok
		else if (d(i) > 0)
		{
			ai = (limits(i,0) - x(i)) / d(i);
			bi = (limits(i,1) - x(i)) / d(i);
		}
		//ok
		else
		{
			ai = (limits(i, 1) - x(i)) / d(i);
			bi = (limits(i, 0) - x(i)) / d(i);
		}
		if (ab[0] < ai )	
			ab[0] = ai;
		if (ab[1] > bi )	
			ab[1] = bi;
	}
	return ab;
}
#endif
#if LAB_NO>6
solution EA(int N, matrix limits, double epsilon, int Nmax, matrix O) //N -> rozmiar problemu, limits -> ograniczenia po³o¿enia x-ów, O -> (sigma) pocz wart odchyl standard, reszta standard
{
	int mi = 20;	//liczba osobników
	int lambda = 40;//generacja potomna -> z mi osobników robimy nastêpn¹ populajê lambda osobników
	solution *P = new solution[mi + lambda];	//populacja skldajaca siê z mi + lambda osobnikow, u³atwienie implementacji
	solution *Pm = new solution[mi];	//populacja pomocnicza
	random_device rd;
	default_random_engine gen;
	gen.seed(static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count()));
	normal_distribution<double> distr(0.0, 1.0);
	uniform_int_distribution<int> distribution(limits(0), limits(1));
	matrix IFF(mi, 1), temp(N, 2);	// IFF -> tablica przystosowañ
	double r, s, s_IFF;	//r, s -> liczby losowe(r nalezy do przedzia³u(0,1)), s_IFF -> suma z IFF
	double tau = 1 / sqrt(2 * N) , tau1 = 1 / sqrt(2 * sqrt(N));	//do mutacji wektora sigm -> sigma prim = sigma * exp(tau1 * r + tau * ri)
	int j_min;
	for (int i = 0; i < mi ; ++i)		//generowanie populacji pocz¹tkowej
	{									//P[i] to ka¿dy osobnik sk³adaj¹cy siê z wektora x i sigma
		P[i].x = matrix(N, 2);
		for (int j = 0; j < N; ++j)
		{
			P[i].x(j, 0) = distribution(gen);
			P[i].x(j, 1) = O(0);
		}
		P[i].fit_fun();
		if (P[i].y < epsilon)
			return P[i];
	}
	while (true)
	{
		s_IFF = 0;
		//selekcja
		for (int i = 0; i < mi ; ++i)	//gen IFF i ich suma	
		{
			IFF(i) = 1 / P[i].y(0);
			s_IFF += IFF(i);
		}
		for (int i = 0; i < lambda ; ++i)	//ko³o ruletki
		{
			r = rd() * s_IFF / rd.max();		//od 0 do s_IFF
			s = 0;
			for (int j = 0; j < mi ; ++j)
			{
				s += IFF(j);
				if (r <= s)
				{
					P[mi + i] = P[j]; //nie jestem pewien
					break;
				}
			}
		}
		//////
		//mutacja
		for (int i = 0; i < lambda ; ++i)
		{
			r = distr(gen);
			for (int j = 0; j < N; ++j)
			{
				P[mi + i].x(j, 1) *=  exp(tau1 * r + tau * distr(gen));	//mutacja sigmy
				P[mi + i].x(j, 0) += P[mi + i].x(j, 1) * distr(gen);	//mutacja wektora x
			}
		}
		////////
		//krzy¿owanie -> generowanie nowej populacji
		for (int i = 0; i < lambda ; i += 2)
		{
			r = distr(gen);
			temp = P[mi + i].x;
			P[mi + i].x = r * temp + (1 - r) * P[mi + i + 1].x;
			P[mi + i + 1].x = r * P[mi + i + 1].x + (1 - r) * temp;
		}
		////////
		//ocena rozwi¹zañ
		for (int i = 0; i < lambda ; ++i)
		{
			P[mi + i].fit_fun();
			if (P[mi + i].y < epsilon)
				return P[mi + i];
		}
		/////////
		//tworzenie nowej populacji o liczebnoœci mi
		for (int i = 0; i < mi ; ++i)
		{
			j_min = 0;
			for (int j = 1; j < lambda ; ++j)
				if (P[j_min].y>P[j].y)
					j_min = j;
			Pm[i] = P[j_min];
			P[j_min].y = 1e10;
		}
		////////
		//prze¿ucenie danych
		for (int i = 0; i < mi ; ++i)
			P[i] = Pm[i];
		///////
		//warunek stopu
		if (solution::f_calls >= Nmax)
			return P[0];
		///////
	}
}
#endif
