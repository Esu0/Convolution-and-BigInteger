#include"convolution.h"
#include<iostream>
#include<vector>
#include<ios>
#include<iomanip>
#include"bigint.h"
#include<chrono>

std::chrono::system_clock::time_point _Start_time, _End_time;

void timer_start()
{
	_Start_time = std::chrono::system_clock::now();
}

void timer_end()
{
	_End_time = std::chrono::system_clock::now();
}

long long timer_elapsed()
{
	return std::chrono::duration_cast<std::chrono::microseconds>(_End_time - _Start_time).count();
}

void timer_print()
{
	std::cout << "execution time: " << (double)timer_elapsed() / 1000 << "ms" << std::endl;
}

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
	using ubigintdec = ubigint<10000000000>;
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
	using bigint = ubigint<10000000000>;
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
int main()
{
	return test5();
}
