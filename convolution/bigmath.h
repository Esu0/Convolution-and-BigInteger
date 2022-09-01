#pragma once
#include"bigint.h"
#include"convolution.h"
#include<queue>
#include<deque>
#include"instant_timer.h"

int _count_factor(unsigned long long n, unsigned long long f)
{
	if (n == 0)return -1;
	int count = 0;
	while (n % f == 0)
	{
		++count;
		n /= f;
	}
	return count;
}

namespace bigmath
{
	template<unsigned long long num_base, bool multi_thread = false>
	ubigint<num_base, multi_thread> pow(const ubigint<num_base, multi_thread>& base, unsigned long exp)
	{
		if (exp == 0)return ubigint<num_base, multi_thread>(1);
		ubigint<num_base, multi_thread> result(base);
		unsigned long mask = msb(exp);
		mask >>= 1;
		while (mask != 0)
		{
			result *= result;
			if (mask & exp)
			{
				result *= base;
			}
			mask >>= 1;
		}
		return result;
	}

	
	template<unsigned long long num_base, bool multi_thread = false>
	ubigint<num_base, multi_thread> factorial(unsigned long x)
	{
		using ubint = ubigint<num_base, multi_thread>;
		if (x <= 1)return ubint(1);
		int* exp = new int[x];
		--exp;
		for (size_t i = 1; i <= x; ++i)exp[i] = 0;
		std::vector<unsigned long> primes;
		unsigned long sqrtx = (unsigned long)(std::sqrt)((double)x);
		size_t i = 2;
		exp[1] = -1;
		while (i <= sqrtx)
		{
			if (exp[i] >= 0)
			{
				primes.push_back((unsigned long)i);
				++(exp[i]);
				size_t j = 2;
				for (size_t j = 2; j <= x / i; j += 1)
				{
					exp[i] += _count_factor(j, i) + 1;
					exp[j * i] = -1;
				}
			}
			++i;
		}
		while (i <= x)
		{
			if (exp[i] >= 0) 
			{
				primes.push_back((unsigned long)i);
				exp[i] = (int)(x / i);
			}
			++i;
		}
		unsigned long* exp2 = new unsigned long[primes.size()];
		i = 0;
		for (size_t j = 1; j <= x; ++j)
		{
			if (exp[j] >= 0)
			{
				exp2[i] = exp[j];
				++i;
			}
		}
		exp++;
		delete[] exp;

		std::deque<ubint> que2;
		unsigned long mask = msb(exp2[0]);
		while (mask)
		{
			std::deque<ubint> que1;
			for (size_t i = 0; i < primes.size(); ++i)
			{
				if (mask & exp2[i])que1.push_back(ubint(primes[i]));
			}
			while (que1.size() >= 2)
			{
				que1.push_back(que1[0] * que1[1]);
				que1.pop_front();
				que1.pop_front();
			}
			if(!que1.empty())que2.push_back(pow<num_base, multi_thread>(que1.front(), mask));
			mask >>= 1;
		}
		while (que2.size() >= 2)
		{
			que2.push_back(que2[0] * que2[1]);
			que2.pop_front();
			que2.pop_front();
		}
		ubigint tmp(std::move(que2.front()));
		return tmp;
	}
}