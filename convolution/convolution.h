#pragma once

#ifndef _CONVOLUTION_H_
#define _CONVOLUTION_H_

#include<cstdlib>
#include<random>
#include<vector>
#include<cassert>

#ifdef _MSC_VER
#include<intrin.h>
#endif

constexpr unsigned long msb(unsigned long x)
{
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

constexpr unsigned long _log2(unsigned long x)
{
	unsigned long answer = 0;
	if (x & 0xffff0000)
	{
		x &= 0xffff0000;
		answer |= 0x00000010;
	}
	if (x & 0xff00ff00)
	{
		x &= 0xff00ff00;
		answer |= 0x00000008;
	}
	if (x & 0xf0f0f0f0)
	{
		x &= 0xf0f0f0f0;
		answer |= 0x00000004;
	}
	if (x & 0xcccccccc)
	{
		x &= 0xcccccccc;
		answer |= 0x00000002;
	}
	if (x & 0xaaaaaaaa)
	{
		x &= 0xaaaaaaaa;
		answer |= 0x00000001;
	}
	return answer;
}

constexpr unsigned long minv_R(unsigned long i, unsigned long t, unsigned long r, unsigned long mod)
{
	return (r > 0) ? ((t % 2) ? minv_R(i << 1, t >> 1, r - 1, mod) : (minv_R(i << 1, (t + mod) >> 1, r - 1, mod) | i)) : 0;
}

constexpr unsigned long bsf(unsigned long x)
{
	return _log2(x & (unsigned long)(-(int)x));
}

unsigned long bsf_fast(unsigned long x)
{
#ifdef _MSC_VER
	unsigned long i;
	_BitScanForward(&i, x);
	return i;
#else
	return bsf(x);
#endif
}

constexpr bool is_pow2(unsigned long long x)
{
	return !(x ^ (x & (unsigned long long)(-(signed long long)x)));
}


template<unsigned long mod>
class mod_productor
{
private:
	static constexpr long long _minus_mod(long long num)
	{
		return (num >= mod) ? (num - mod) : num;
	}
public:
	static_assert(mod % 2, "mod is even number.");
	static constexpr unsigned long R = msb(mod) << 1;
	static constexpr unsigned long R2 = (unsigned long)((unsigned long long)R * R % mod);
	static constexpr unsigned long mask = R - 1;
	static constexpr unsigned long shift = _log2(R);
	static constexpr unsigned long invN = minv_R(1, 0, shift, mod);

	static constexpr long long MR(long long num)
	{
		return _minus_mod((num + ((num * invN) & mask) * mod) >> shift);
	}
	static constexpr long long mprod(long long _Left, long long _Right)
	{
		return MR(MR(_Left * _Right) * R2);
	}
	//return num % mod,0 < num < mod * mod
	static constexpr long long modulo(long long num)
	{
		return MR(MR(num) * R2);
	}
};

constexpr unsigned long long mpow(unsigned long long base, unsigned long long exp, unsigned long long div)
{
	if (exp <= 0)return 1;
	if (exp == 1)return base;
	unsigned long long tmp = mpow(base, exp / 2, div);
	if (exp % 2)
		return (((tmp * tmp) % div) * base) % div;
	else
		return (tmp * tmp) % div;
}

//get primitive root
unsigned long get_prim_root(unsigned long mod)
{
	std::vector<unsigned long> prime_factor = std::vector<unsigned long>();
	--mod;
	unsigned long m = mod;
	if (m % 2 == 0)
	{
		prime_factor.push_back(2);
		m /= 2;
	}
	while (m % 2 == 0)m /= 2;
	for (unsigned long i = 3; i * i < m; i += 2)
	{
		if (m % i == 0)
		{
			prime_factor.push_back(i);
			m /= i;
		}
		while (m % i == 0)m /= i;
	}
	if (m != 1)prime_factor.push_back(m);

	static std::random_device seed;
	static std::mt19937_64 mt(seed());
	std::uniform_int_distribution<> random(2, mod - 1);
	unsigned long candidate, prod;
	bool loopend;
	while (1)
	{
		candidate = random(mt);
		loopend = true;
		for (auto tmp : prime_factor)
		{
			prod = mod / tmp;
			if (mpow(candidate, prod, mod + 1) == 1)
			{
				loopend = false;
				break;
			}
		}
		if (loopend)break;
	}
	return candidate;
}

//aとmodは互いに素
constexpr long long minv(long long a, long long mod)
{
	long long b = mod, u = 1, v = 0;
	while (b)
	{
		long long tmp = u;

		u = v;
		v = tmp - (a / b) * v;

		tmp = a;
		a = b;
		b = tmp % b;
	}
	u %= mod;
	if (u < 0)u += mod;
	return u;
}

template<unsigned long mod, unsigned long g>
class ntt
{
private:


public:
	static constexpr unsigned char rank = (unsigned char)bsf(mod - 1);

	static_assert(rank >= 2, "rankは2以上必要");

	//root[i] = 2^i乗根
	//root[i]を2^i乗すると1になる
	unsigned long root[rank + 1];
	unsigned long root3[rank - 2];

	//iroot[i] * root[i] = 1;
	unsigned long iroot[rank + 1];
	unsigned long iroot2[rank - 1];
	unsigned long iroot3[rank - 2];

	constexpr ntt()
	{
		root[rank] = (unsigned long)mpow(g, (mod - 1) >> rank, mod);
		iroot[rank] = (unsigned long)minv(root[rank], mod);
		for (unsigned char i = rank; i > 0;)
		{
			--i;
			root[i] = (unsigned long long)root[i + 1] * root[i + 1] % mod;
			iroot[i] = (unsigned long long)iroot[i + 1] * iroot[i + 1] % mod;
		}
		for (unsigned char i = 0; i < rank - 2; ++i)
		{
			root3[i] = (unsigned long long)iroot[2] * root[i + 2] % mod * root[i + 3] % mod;
			iroot3[i] = (unsigned long long)root[2] * iroot[i + 2] % mod * iroot[i + 3] % mod;
		}
		for (unsigned char i = 0; i < rank - 1; ++i)
		{
			iroot2[i] = (unsigned long long)root[1] * iroot[i + 1] % mod * iroot[i + 2] % mod;
		}
	}
};


//制約:gはmodの原始根、mod-1は素因数になるべくたくさん2を含む
//modは素数
//vは(size*8) byte以上のメモリが確保されている
//size > 2, size = 2^Nと表せる
template<unsigned mod, unsigned long g, bool inv = false>
void number_theoretic_transform(long long *v, unsigned long size)
{
	assert(is_pow2(size) && size > 2);
	static ntt<mod, g> table;
	static constexpr unsigned long long mod2 = (unsigned long long)mod * mod;
	using ull = unsigned long long;
	using ul = unsigned long;
	using stype = unsigned long;
	//逆変換
	if constexpr (inv)
	{
		stype n = size;
		stype o = bsf_fast(n);
		stype order = o;

		ul w_4 = table.iroot[2];
		if (order & 1ul)
		{
			ul w = 1;
			for (stype i = 0; i < n; i += 4)
			{
				ull a1 = v[i];
				ull a2 = v[i + 1];
				ull a3 = v[i + 2];
				ull a4 = v[i + 3];
				v[i] = (a1 + a2) % mod;
				v[i + 1] = (ull)w * (a1 + mod - a2) % mod;
				v[i + 2] = (a3 + a4) % mod;
				v[i + 3] = (ull)w_4 * w % mod * (a3 + mod - a4) % mod;
				w = (ul)((ull)w * table.iroot3[bsf_fast(~(i >> 2))] % mod);
			}
			order ^= 1ul;
		}
		while (order > 0)
		{
			stype blocks = 1 << (order - 2);
			ul w = 1;
			for (stype k = 0; k < blocks; ++k)
			{
				stype offset = k << (o - order + 2);
				stype loops = n >> order;
				ul w2 = (ul)((ull)w * w % mod);
				ul w3 = (ul)((ull)w2 * w % mod);
				for (stype i = 0; i < loops; ++i)
				{
					ul j = i + offset;
					ull a1 = v[j];
					ull a2 = v[j + loops];
					ull a3 = v[j + 2 * loops];
					ull a4 = v[j + 3 * loops];
					ull a3a4 = w_4 * (a3 + mod - a4) % mod;
					v[j] = (a1 + a2 + a3 + a4) % mod;
					v[j + loops] = (ull)w * (a1 + mod - a2 + a3a4) % mod;
					v[j + 2 * loops] = (ull)w2 * (a1 + a2 + 2 * mod - a3 - a4) % mod;
					v[j + 3 * loops] = (ull)w3 * (a1 + 2 * mod - a2 - a3a4) % mod;
				}
				w = (ul)((ull)w * table.iroot3[bsf_fast(~k)] % mod);
			}
			order -= 2;
		}
	}
	//正変換
	else
	{
		stype n = size;
		stype o = bsf_fast(n);
		stype order = o;

		ul w_4 = table.root[2];
		while (order > 1)
		{
			stype blocks = n >> order;
			unsigned long w = 1;
			for (stype k = 0; k < blocks; ++k)
			{
				stype offset = k << order;
				stype loops = 1 << (order - 2);
				ul w2 = (ul)((ull)w * w % mod);
				ul w3 = (ul)((ull)w2 * w % mod);
				for (stype i = 0; i < loops; ++i)
				{
					ul j = i + offset;
					ull a1 = v[j];
					ull a2 = w * v[j + loops];
					ull a3 = w2 * v[j + 2 * loops];
					ull a4 = w3 * v[j + 3 * loops];
					ull a2a4 = (ull)w_4 * ((a2 + mod2 - a4) % mod);
					v[j] = (a1 + a2 + a3 + a4) % mod;
					v[j + loops] = (a1 + a3 + mod2 * 2 - (a2 + a4)) % mod;
					v[j + 2 * loops] = (a1 + mod2 - a3 + a2a4) % mod;
					v[j + 3 * loops] = (a1 + 2 * mod2 - a3 - a2a4) % mod;
				}
				w = (ul)((ull)w * table.root3[bsf_fast(~k)] % mod);
			}
			order -= 2;
		}
		if (order == 1)
		{
			unsigned long w = 1;
			for (stype i = 0; i < n; i += 4)
			{
				ull a1 = v[i];
				ull a2 = (ull)w * v[i + 1];
				ull a3 = v[i + 2];
				ull a4 = (ull)w * v[i + 3] % mod * w_4;
				v[i] = (a1 + a2) % mod;
				v[i + 1] = (a1 + mod2 - a2) % mod;
				v[i + 2] = (a3 + a4) % mod;
				v[i + 3] = (a3 + mod2 - a4) % mod;
				w = (ul)((ull)w * table.root3[bsf_fast(~(i >> 2))] % mod);
			}
		}
	}
}


//num1はn1ワード、num2はn2ワード、resultはn1+n2ワードのメモリが確保されているものとする
template<unsigned long long num_base = (1ull << 32)>
void multiply_naive_unsafe(long long* num1, long long* num2, long long* result, unsigned long n1, unsigned long n2)
{
	static_assert(num_base <= 1, "進数が不正");
	if constexpr (is_pow2(num_base))
	{
		unsigned long long mask = num_base - 1;
		unsigned char shift = _log2(mask);
		unsigned long long carry = 0;
		for (unsigned long i = 0; i < n1 + n2; i++)result[i] = 0;
		for (unsigned long i = 0; i < n1; ++i)
		{
			for (unsigned long j = 0; j < n2; ++j)
			{
				unsigned long k = i + j;
				unsigned long long tmp = num1[i] * num2[j] + carry + result[k];
				carry = tmp >> shift;
				result[k] = tmp & mask;
			}
		}
	}
	else
	{
		unsigned long long carry = 0;
		for (unsigned long i = 0; i < n1 + n2; i++)result[i] = 0;
		for (unsigned long i = 0; i < n1; ++i)
		{
			for (unsigned long j = 0; j < n2; ++j)
			{
				unsigned long k = i + j;
				unsigned long long tmp = num1[i] * num2[j] + carry + result[k];
				carry = tmp / num_base;
				result[k] = tmp % num_base;
			}
		}
	}
}
#endif