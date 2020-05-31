//#include <pch.h>
#include <iostream>
#include <iostream>
#include<stdlib.h>
#include <cmath>
#include <time.h>
#include <fstream>
#include <string>


using namespace std;


// greatest common div. using Euc algo.
long long gcd(long long a, long long b)
{
	// return (b=!0)? gcd(b,a%b):a;

	long long v;
	while (true)
	{
		if (a%b == 0)
			return b;
		v = b;
		b = a % b;
		a = v;
	}
	return 0;
}
long long ext_gcd(long a, long b, long long *x, long long *y)
{
	if (a == 0)
	{
		*x = 0;
		*y = 1;
		return b;
	}
	long long x1, y1;
	long long gcd = ext_gcd(b%a, a, &x1, &y1);
	*x = y1 - (b / a)*x1;
	*y = x1;
	return gcd;

	/*
	tuple<int, int, int> extended_gcd(int a, int b)
	{
		if (a == 0)
			return make_tuple(b, 0, 1);

		int gcd, x, y;

		// unpack tuple  returned by function into variables
		tie(gcd, x, y) = extended_gcd(b % a, a);

		return make_tuple(gcd, (y - (b / a) * x), x);
	}*/
}
//where b1,b2..etc is the binary rep of n.
long long binpow(long long a, long long n, long long m) {
	long long res = 1;
	a = a % m;				// we take the a mod m
	while (n > 0)			// when n=1 then result = 0
	{
		if ((n % 2) == 1)	// odd nums in binary rep starts with 1 (e.g. 1101 is odd, 1100 is even)
		{
			res = (res*a) % m;
		}
		a = a * a % m;
		n = n / 2;
	}
	return res;
}//a^(n) is b1*a^(1) + b2*a^(2) + b3*a^(4) + b4*a^(8)

long long inv(long long a, long long m)
{
	a = a % m;
	for (long long i = 1; i < m; i++)
	{
		if ((a*i) % m == 1)
			return i;
	}
	return 0;
}
int is_prime(int n)
{
	int i;
	for (i = 2; i <= n / 2; ++i) {
		if (n%i == 0)
			return 0;
	}

	return 1;
}
int is_primitive(long long a, long long m)
{
	//Prime Decomp (prime factorization)
	// it means it checks for all the primes that n is divisble by.
	int factors[100], d = 2, i = 0, j = 0;
	long long n = m - 1;

	while (n > 1)
	{
		if (n % d == 0)
		{
			n = n / d;
			if (j == 0)
			{

				factors[i] = d;
				j++;
				i++;
				cout << (n*d) << " % " << d << endl;
			}
		}
		else
		{
			//cout << "next" << endl;
			d++;
			j = 0;
		}
	}
	// check if its primitive.
	for (j = 0; j < i; j++)
	{
		if (binpow(a, (m - 1) / factors[j], m) == 1)
		{
			cout << a << " is not primitive" << endl;
			return 0;
		}
	}
	cout << a << " is primitive" << endl;
	return 1;
	/*
	// using int instead of bool.
	int primitiv = 1;
	for (j = 0; (j < i) && primitiv; j++) {
		if (binexp(a, n / factors[j], m) == 1) {
			primitiv = 0;
		}
	}
	return primitiv;
*/
}
long long findprime(long long num1, long long num2)
{
	long long p;
	do
	{
		p = num1 + rand() % (num2 - num1);
	} while (!is_prime(p));
	return p;
}
long long findprimitive(long long p)
{
	int g = 2;
	while (!is_primitive(g, p)) {
		g++;
	}
	return g;
}

void print_prime_el(long long m) {
	if (!is_prime(m))
	{
		cout << m << " is not a prime." << endl;
		return;
	}

	int j;
	for (j = 2; j < m; j++)
	{
		cout << ">> " << j << ": " << endl;
		if (is_primitive(j, m))
		{
			cout << j << " is primitive root of " << m << endl;
		}
	}
	//cout << "ERROR 404: primitive not found" << endl;
}
long long power_of(long long n, long long m)
{
	long long res = 1;
	for (int i = 0; i < m ; i++)
		res *= n;

	return res;
}


//needs a as input.
//if there is any int a where a^n mod m is not 1 or -1
//then its not a prime
//ALSO, for finding fermat's pseudo primes!!

int MRTest(long long a, long long m)
{
	int k;
	long long n = m - 1;

	for (k = 0; n % 2 == 0; k++)
		n = n / 2;

	//cout << "m prime is: " << n << endl;
	long long mr = binpow(a, n, m);
	cout << a << " ^ " << n << "(mod " << m << ") is: " << mr << endl;
	for (int i = 0; i < k; i++)
	{
		cout << "a" << i << " is: " << mr << endl;

		if (m - mr == 1)
			return 1;
		if (mr == 1)
			if (i == 0)
				return 1;
			else
				return 0;
		mr = mr * mr%m;
	}
	return 0;
}

char caesar(char x, char b)
{
	return x = x + b;
}

char caesar_dec(char x, char b)
{
	return x = x - b;
}

char affine(char x, long long a, long long b)
{
	return (char) (a * x + b) % power_of(2, 32);
}
char affine_dec(char x, long long a_inv, long long b)
{
	//long long a_inv = inv(a, power_of(2, 32));
	cout << a_inv <<endl;
	return (char) (a_inv * (x - b)) % power_of(2, 32);
}



void elgamal(long long a, long long b, long long *p, long long *g, long long *x, long long *y)
{
	*p = findprime(a, b);
	*g = findprimitive(*p);
	*x = (rand() % (*p - 1)) + 1;
	*y = power_of(*g, *x);

}
void enc_elgamal(long long m, long long p, long long g, long long y, long long *c1, long long *c2)
{
	char fname[100];
	char c;

	cout << "Filename to read from: ";
	cin >> fname;
	string sfname;
	ifstream infile(fname, ios::binary);
	//	sfname = strcat_s(fname, ".aff");
	sfname = "enc_ft.txt";
	ofstream outfile(sfname, ios::binary);
	while (!infile.eof())
	{
		infile.read(&c, 1);
		m = long long(c);

		long long k = rand() % (p - 2) + 1;
		*c1 = binpow(g, k, p);
		*c2 = (m * binpow(y, k, p)) % p;       //enc

		if (!infile.eof()) {
			c = (char)*c2;
			outfile.write(&c, 1);
		}
	}
	infile.close();
	outfile.close();

}
void dec_elgamal(long long c1, long long c2, long long p, long long g, long long x)
{
	char fname[100];
	char c;

	cout << "Filename to decrypt: ";
	cin >> fname;
	string sfname;
	ifstream infile(fname, ios::binary);
	sfname = strcat("a_",fname);
	//sfname = "dec_ft.txt";
	ofstream outfile(sfname, ios::binary);
	while (!infile.eof())
	{
		infile.read(&c, 1);
		c2 = (long long)c;
		long long b = binpow(c1, (p - x - 1), p);
		long long m = (b*c2) % p;
		//return m;

		if (!infile.eof()) {
			c = (char)m;
			outfile.write(&c, 1);
		}
	}
	infile.close();
	outfile.close();
}

void rsa(long long a, long long b, long long *p, long long *q, long long *n, long long *phin, long long *e, long long *d)
{
	*p = findprime(a, b);
	*q = findprime(a, b);
	*n = (*q) * (*p);
	*phin = ((*q) - 1)*((*p) - 1);


	// enc key.
	do
	{
		*e = (rand() % (*phin - 3) + 2);
	} while (gcd(*e, *phin) != 1);

	//dec key.
	*d = inv(*e, *phin);
}
/*
long long rsa_enc(long long m, long long e, long long n)
{
	return binpow(m, e, n);
}
long long rsa_dec(long long c, long long d, long long n)
{
	return binpow(c, d, n);
}
*/

void rsa_enc(long long e, long long n)
{
	char fname[100];
	char c;

	cout << "Filename to read from: ";
	cin >> fname;
	string sfname;
	ifstream infile(fname, ios::binary);
	//	sfname = strcat_s(fname, ".aff");
	sfname = "enc_ft.txt";
	ofstream outfile(sfname, ios::binary);
	while (!infile.eof())
	{
		infile.read(&c, 1);

		c = (char)binpow(c, e, n);         //enc

		if (!infile.eof()) {
			outfile.write(&c, 1);
		}
	}
	infile.close();
	outfile.close();

}
void rsa_dec(long long d, long long n)
{
	char fname[100];
	char c;

	cout << "Filename to decrypt: ";
	cin >> fname;
	string sfname;
	ifstream infile(fname, ios::binary);
	//	sfname = strcat(fname, ".aff");
	sfname = "dec_ft.txt";
	ofstream outfile(sfname, ios::binary);
	while (!infile.eof())
	{
		infile.read(&c, 1);
		c = binpow(long long(c), d, n);        //Dec
		//c = (c*e) % n;         //enc

		if (!infile.eof()) {
			outfile.write(&c, 1);
		}
	}
	infile.close();
	outfile.close();

}
void rsa(long long p, long long q, long long *x, long long *y) {

	long long n = q * p;
	long d, e;
	//long long a = sqrt(n) + 1;
	//long long b = sqrt(a*a - n);
	long long phin = (q - 1)*(p - 1);

	do
	{
		e = (rand() % (phin - 3) + 2);
	} while (gcd(e, phin) != 1);

	/*
	if (d < 0) {
		d = d + phin;
	}*/
	d = inv(e, phin);

	cout << phin << ", " << e << ", " << d << endl;
	*x = e;
	*y = d;
	rsa_enc(e, n);
	rsa_dec(d, n);
}


void read_file()
{
	char fname[100];
	char ch;

	cout << "Filename to read from: ";
	cin >> fname;
	string sfname;
	ifstream infile(fname, ios::binary);
	//	sfname = strcat_s(fname, ".aff");
	sfname = "enc_ft.txt";
	ofstream outfile(sfname, ios::binary);
	while (!infile.eof())
	{
		infile.read(&ch, 1);
		//enc or dec ------

		if (!infile.eof())
		{
			outfile.write(&ch, 1);
		}
	}
	infile.close();
	outfile.close();
}


int main()
{
	srand((unsigned)time(NULL));

	long long a, b, x, y;
	long long p, q, e, d, n, phin;
	long long c1, c2, c, m;
	char move;
	//int flag = 1;

	cout << power_of(2, 32)<<endl;
	cout << "AFFINE\n\nEnter the first number: ";
	cin >> a;
	cout << "Enter the second number: ";
	cin >> b;

	char fname[100];
	char ch;

	cout << "Filename to read from: ";
	cin >> fname;
	string sfname;
	ifstream infile(fname, ios::binary);
	sfname = strcat("a_",fname);
	//sfname = "dec_ft.txt";
	ofstream outfile(sfname, ios::binary);
	while (!infile.eof())
	{
		infile.read(&ch, 1);

		//long long a_inv = inv(a, power_of(2, 32));
		ch=affine(ch,a_inv,b);

		if (!infile.eof())
		{
			outfile.write(&ch, 1);
		}
	}
	infile.close();
	outfile.close();

	return 0;
}
