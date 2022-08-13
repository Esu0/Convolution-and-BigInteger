#include"convolution.h"
#include<iostream>
#include<vector>
#include<ios>
#include<iomanip>
#include<convenience.h>
#include<fourier_trans.h>
#include"bigint.h"

using namespace convenient;
using namespace instant;

int test1()
{
	std::vector<long long> a = { 3,2,14,1,3,4,1,4 };
	number_theoretic_transform<17, 5>(a.data(), (unsigned long)a.size());
	std::cout << to_string(a) << std::endl;
	number_theoretic_transform<17, 5, true>(a.data(), (unsigned long)a.size());
	std::cout << to_string(a) << std::endl;
	return 0;
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

int test4()
{
	std::vector<long long> a = {2,4,3,7,3};
	std::vector<long long> b = { 3, 4, 4, 3 };
	std::vector<long long> c(9);
	multiply_naive_unsafe(a.data(), b.data(), c.data(), (unsigned long)a.size(), (unsigned long)b.size());
	std::cout << to_string(c) << std::endl;
	std::cout << std::hex << 0x37342ull * 0x3443 << std::endl;
	return 0;
}

int test5()
{
	using ubigintdec = ubigint<1ull << 32>;
	ubigintdec a, b;
	a.random(200000);
	b.random(200000);
	a.dump("../results/result_a.txt", " ");
	b.dump("../results/result_b.txt", " ");
	timer_start();
	a.multiply_ntt(b);
	timer_end();
	timer_print();
	a.dump("../results/result_ab.txt", " ");
	return 0;
}

int test6()
{
	using ubigintdec = ubigint<1000000000>;
	ubigintdec a, b;
	a.random(10);
	b.random(10);
	a.dump();
	b.dump();
	a *= b;
	a.dump();
	return 0;
}

int main()
{
	return test6();
}
