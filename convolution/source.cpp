#include"convolution.h"
#include<iostream>
#include<vector>
#include<ios>
#include<iomanip>
#include"bigint.h"
#include<chrono>
#include"bigmath.h"
#include"instant_timer.h"

int test2()
{
	using ubigintbin = ubigint<1 << 4>;
	ubigintbin a("0xabc22200086a8c8bd7e9fa9fa9d7cbe97a");
	a.dump();
	return 0;
}

int test3()
{
	using ubigintdec = ubigint<100000>;
	ubigintdec a("12345678123345767434");
	ubigintdec b("8342384357843445342");
	a.dump();
	b.dump();
	a *= b;
	a.dump();
	return 0;
}

int test5()
{
	using ubigintdec = ubigint<1ull << 16>;
	ubigintdec a, b;
	a.random(1ull << 26);
	b.random(1ull << 26);
	//a.dump("../results/result_a.txt", " ");
	//b.dump("../results/result_b.txt", " ");
	timer_start();
	a = a.multiply_ntt(b);
	timer_end();
	timer_print();
	//a.dump("../results/result_ab.txt", " ");
	return 0;
}

int test6()
{
	using ubigintdec = ubigint<1000000000>;
	ubigintdec a, b;
	a.random(10);
	b = -2323438937123;
	a.dump();
	b.dump();
	a.multiply_naive(b).dump();
	a.multiply_ntt(b).dump();
	//a.dump();
	return 0;
}

int test7()
{
	using bigint = ubigint<1000000000>;
	bigint a, b;
	a.random(5);
	b.random(5);
	a.dump();
	b.dump();
	if (a == a)std::cout << "a == a" << std::endl;
	if (a == b)std::cout << "a == b" << std::endl;
	if (a != b)std::cout << "a != b" << std::endl;
	if (a >= b)std::cout << "a >= b" << std::endl;
	if (a > b)std::cout << "a > b" << std::endl;
	if (a <= b)std::cout << "a <= b" << std::endl;
	if (a < b)std::cout << "a < b" << std::endl;

	return 0;
}

int test8()
{
	using bigint = ubigint<1000000000, true>;
	bigint a, b;
	a.random(100000);
	b.random(100000);
	std::vector<long long> v1 = a.dat, v2 = b.dat;
	auto v3 = bigint::multiply(v1.data(), v2.data(), v1.size(), v2.size());
	for (size_t i = 0; i < 100; ++i)
	{
		std::cout << v3[i] << ' ';
	}
	std::cout << std::endl;
	return 0;
}

int test9()
{
	using bigint = ubigint<1000000000, true>;
	bigint a, b;
	a.random(111112);
	b.random(111112);

	timer_start();
	auto c = a * b;
	timer_end();
	timer_print();

	return 0;
}

int test10()
{
	using bigint = ubigint<100000>;
	bigint a;
	a.random(128);
	a.dump();
	std::cout << std::endl;
	bigint b = a.inverse(128);
	b.dump();
	std::cout << std::endl;
	bigint c = a * b;
	c.dump(c.dat.size() - 128, c.dat.size());
	return 0;
}

int test11()
{
	using bigint = ubigint<1000000000, true>;
	bigint a, b;
	size_t len = 111112;
	a.random(len);
	timer_start();
	b = a.inverse(len);
	timer_end();
	timer_print();
	
	bigint c = a * b;
	c >>= len - 1;
	bigint d = 1;
	d <<= len;
	d -= c;
	d.dump();
	std::cout << d.dat.size() << std::endl;
	return 0;
}

int test12()
{
	using bigint = ubigint<100000>;
	bigint a = 42, b = 99999;
	a <<= 10;
	b <<= 5;

	bigint c = a / b;
	c.dump();
	c = a - c * b;
	b.dump();
	c.dump();
	if (c < b)std::cout << "ok" << std::endl;
	else std::cout << "ng" << std::endl;
	return 0;
}

int test13()
{
	using bigint = ubigint<100000>;
	bigint a, b;
	a.random(5);
	b.random(3);
	bigint c = a % b;
	a.dump();
	b.dump();
	c.dump();
	
	return 0;
}

int test14()
{
	using bigint = ubigint<1000000000, true>;
	bigint a, b;
	size_t n = 222223;
	a.random(n);
	b.random(n / 2);

	timer_start();
	auto c = a / b;
	timer_end();
	timer_print();

	if (a - b * c < b)std::cout << "ok" << std::endl;
	return 0;
}

int test15()
{
	using bigint = ubigint<1000000000>;
	timer_start();
	bigint a = bigmath::factorial<1000000000>(1000000);
	timer_end();
	timer_print();
	a.dump("../results/fact.txt", " ");
	return 0;
}

int main()
{
	return test15();
}
