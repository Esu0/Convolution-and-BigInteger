#pragma once

#ifndef _BIGINT_H_
#define _BIGINT_H_

#include"convolution.h"
#include<string>
#include<iostream>
#include<random>
#include<stdexcept>

inline unsigned char ctoi_hex(char c)
{
	switch (c)
	{
	case 'a':
	case'A':
		return 10;
		break;
	case 'b':
	case 'B':
		return 11;
		break;
	case 'c':
	case 'C':
		return 12;
		break;
	case 'd':
	case 'D':
		return 13;
		break;
	case 'e':
	case 'E':
		return 14;
		break;
	case 'f':
	case 'F':
		return 15;
		break;
	default:
		return c - '0';
		break;
	}
}

inline unsigned char ctoi_dec(char c)
{
	return c - '0';
}

template<unsigned long long num_base>
class ubigint
{
	static_assert(num_base >= 2, "進数が不正");
private:
	std::vector<long long> dat;
	long long max = num_base;

	void delete_zero()
	{
		while (dat.size() > 1 && dat.back() == 0)
		{
			dat.pop_back();
		}
	}

	void fix_carry_force()
	{
		long long carry = 0;
		long long* p = dat.data();
		for (size_t i = 0; i < dat.size(); ++i)
		{
			p[i] += carry;
			carry = p[i] / (long long)num_base;
			p[i] %= (long long)num_base;
			if (p[i] < 0)
			{
				p[i] += num_base;
				--carry;
			}
		}
		while (carry != 0 && carry != -1)
		{
			long long tmp = carry % num_base;
			carry /= num_base;
			if (tmp < 0)
			{
				tmp += num_base;
				--carry;
			}
			dat.push_back(tmp);
		}
		if (carry == -1)
		{
			dat.back() -= num_base;
		}
		delete_zero();
		max = num_base;
	}

	void fix_carry()
	{
		if (max >= LLONG_MAX >> 1)fix_carry_force();
	}

	explicit ubigint(std::vector<long long>&& v) :dat(v)
	{}
public:
	static constexpr unsigned long long max_size = 1ull << 26;

	ubigint(): dat(1,0)
	{}

	ubigint(long long num) : dat(1, num)
	{
		fix_carry_force();
	}

	explicit ubigint(const char* str)
	{
		if (str[0] == '0')
		{
			if (str[1] == 'x' || str[1] == 'X')
			{
				dat = std::move(ubigint::hex(&str[2]).dat);
				return;
			}
		}
		dat = std::move(ubigint::dec(str).dat);
	}

	//基数変換16->任意
	static ubigint<num_base> hex(const char* str)
	{
		if constexpr (is_pow2(num_base))
		{
			static constexpr unsigned char d = (unsigned char)_log2(num_base);
			unsigned long long len_str = strnlen_s(str, ubigint<num_base>::max_size * d / 4);
			size_t dat_len = len_str * 4 / d;
			std::vector<long long> v(dat_len + 1);
			v.pop_back();

			size_t j = len_str;
			size_t bit_count = 0;
			unsigned long long num = 0;
			for (size_t i = 0; i < dat_len; ++i)
			{
				while (bit_count < d)
				{
					--j;
					num |= (unsigned long long)ctoi_hex(str[j]) << bit_count;
					bit_count += 4;
				}
				v[i] = num & (num_base - 1);
				num >>= d;
				bit_count -= d;
			}
			while (j > 0)
			{
				--j;
				num |= (unsigned long long)ctoi_hex(str[j]) << bit_count;
				bit_count += 4;
			}
			if (num != 0)
			{
				v.push_back(num);
			}
			return ubigint(std::move(v));
		}
		else
		{
			return ubigint();
		}
	}

	//基数変換10->任意
	static ubigint<num_base> dec(const char* str)
	{
		if constexpr (is_pow<10>(num_base))
		{
			static constexpr unsigned char d = (unsigned char)_log10(num_base);
			size_t len_str = strnlen_s(str, max_size * d);
			size_t dat_len = len_str / d;
			ubigint result;
			std::vector<long long> v(dat_len + 1);
			v.pop_back();

			unsigned long long num;
			unsigned long long times;
			size_t j = len_str;
			for (size_t i = 0; i < dat_len; ++i)
			{
				times = 1;
				num = 0;
				for (size_t k = 0; k < d; ++k)
				{
					--j;
					num += ctoi_dec(str[j]) * times;
					times *= 10;
				}
				v[i] = num;
			}
			times = 1;
			num = 0;
			while (j > 0)
			{
				--j;
				num += ctoi_dec(str[j]) * times;
				times *= 10;
			}
			if (num != 0)v.push_back(num);
			return ubigint(std::move(v));
		}
		else
		{
			return ubigint();
		}
	}

	//コピー代入コンストラクタ
	ubigint(const ubigint& _Right): dat(_Right.dat), max(_Right.max)
	{}

	//ムーブ代入コンストラクタ
	ubigint(ubigint&& _Right)noexcept : dat(std::move(_Right.dat)), max(_Right.max)
	{}

	//コピー代入演算子
	ubigint& operator=(const ubigint& _Right)
	{
		dat = _Right.dat;
		max = _Right.max;
		return *this;
	}

	//ムーブ代入演算子
	ubigint& operator=(ubigint&& _Right)noexcept
	{
		dat = std::move(_Right.dat);
		max = _Right.max;
		return *this;
	}

	void dump()
	{
		fix_carry_force();
		for (size_t i = dat.size(); i > 0;)
		{
			--i;
			std::cout << dat[i] << ' ';
		}
		std::cout << std::endl;
	}

	void dump(const char* str)
	{
		fix_carry_force();
		for (size_t i = dat.size(); i > 0;)
		{
			--i;
			std::cout << dat[i] << str;
		}
		std::cout << std::endl;
	}
	//前置インクリメント
	ubigint& operator++()
	{
		++dat[0];
		++max;
		fix_carry();
		return *this;
	}

	//後置インクリメント
	ubigint operator++(int)
	{
		ubigint<num_base> tmp = *this;
		this->operator++();
		return tmp;
	}

	//前置デクリメント
	ubigint& operator--()
	{
		--dat[0];
		++max;
		fix_carry();
		return *this;
	}

	//後置デクリメント
	ubigint operator--(int)
	{
		ubigint<num_base> tmp = *this;
		this->operator--();
		return tmp;
	}

	//複合代入(足し算)
	ubigint& operator+=(const ubigint& _Right)
	{
		if (_Right.dat.size() > dat.size())dat.resize(_Right.dat.size(), 0);
		for (size_t i = 0; i < _Right.dat.size(); ++i)dat[i] += _Right.dat[i];
		max += _Right.max;
		fix_carry();
		return *this;
	}

	//足し算
	[[nodiscard]]
	ubigint operator+(ubigint _Right)const&
	{
		_Right.operator+=(*this);
		return _Right;
	}

	//複合代入(引き算)
	ubigint& operator-=(const ubigint& _Right)
	{
		if (_Right.dat.size() > dat.size())dat.resize(_Right.dat.size(), 0);
		for (size_t i = 0; i < _Right.dat.size(); ++i)dat[i] -= _Right.dat[i];
		max += _Right.max;
		fix_carry();
		return *this;
	}

	//Right <- Left - Right
	void sub_reverse(ubigint& _Right)const&
	{
		if (dat.size() > _Right.dat.size())_Right.dat.resize(dat.size(), 0);
		for (size_t i = 0; i < dat.size(); ++i)_Right.dat[i] = dat[i] - _Right.dat[i];
		_Right.max += max;
		_Right.fix_carry();
	}

	//引き算
	[[nodiscard]]
	ubigint operator-(ubigint _Right)const&
	{
		sub_reverse(_Right);
		return _Right;
	}

	ubigint multiply_naive(const ubigint& num1, const ubigint& num2)
	{
		std::vector<long long> resultv(num1.dat.size() + num2.dat.size());
		multiply_naive_unsafe<num_base>(num1.dat.data(), num2.dat.data(), resultv.data(), (unsigned long)num1.dat.size(), (unsigned long)num2.dat.size());
		return ubigint(std::move(resultv));
	}


	void multiply_ntt(const ubigint& num)
	{
		static constexpr unsigned long ntt_mod1 = 2013265921;
		static constexpr unsigned long ntt_mod2 = 1811939329;
		static constexpr unsigned long border1 = ((ntt_mod1 + num_base - 2) / (num_base - 1) + num_base - 2) / (num_base - 1);
		static constexpr unsigned long border2 = (unsigned long)((((unsigned long long)ntt_mod1 * ntt_mod2 + num_base - 2) / (num_base - 1) + num_base - 2) / (num_base - 1));
		unsigned long result_size = (unsigned long)(dat.size() + num.dat.size());
		unsigned long n = (unsigned long)(std::min)(dat.size(), num.dat.size());
		std::cout << border2 << std::endl;
		if (n < border1)
		{
			unsigned long new_size = (unsigned long)ceil_pow2(result_size);
			std::vector<long long> buf1(new_size);
			std::vector<long long> buf_num1(new_size);

			//コピー1回目
			for (size_t i = 0; i < dat.size(); ++i)
			{
				long long tmp = dat[i] % (long long)ntt_mod1;
				buf_num1[i] = tmp < 0 ? tmp + ntt_mod1 : tmp;
			}
			for (size_t i = 0; i < num.dat.size(); ++i)
			{
				long long tmp = num.dat[i] % (long long)ntt_mod1;
				buf1[i] = tmp < 0 ? tmp + ntt_mod1 : tmp;
			}

			//ntt_mod1でntt
			number_theoretic_transform<ntt_mod1, 137>(buf_num1.data(), new_size);
			number_theoretic_transform<ntt_mod1, 137>(buf1.data(), new_size);
			for (size_t i = 0; i < new_size; ++i)
			{
				buf1[i] = buf1[i] * buf_num1[i] % ntt_mod1;
			}
			number_theoretic_transform<ntt_mod1, 137, true>(buf1.data(), new_size);
			long long inv_size = _minv<ntt_mod1>(new_size);
			dat.resize(result_size);
			for (size_t i = 0; i < result_size; ++i)
			{
				dat[i] = buf1[i] * inv_size % ntt_mod1;
			}
			fix_carry_force();
		}
		else if (n < border2)
		{
			unsigned long new_size = (unsigned long)ceil_pow2(result_size);
			std::vector<long long> buf1(new_size);
			std::vector<long long> buf_num1(new_size);
			std::vector<long long> buf2(new_size);

			//コピー1回目
			for (size_t i = 0; i < dat.size(); ++i)
			{
				long long tmp = dat[i] % (long long)ntt_mod1;
				buf_num1[i] = tmp < 0 ? tmp + ntt_mod1 : tmp;
			}
			for (size_t i = 0; i < num.dat.size(); ++i)
			{
				long long tmp = num.dat[i] % (long long)ntt_mod1;
				buf1[i] = tmp < 0 ? tmp + ntt_mod1 : tmp;
			}

			//ntt_mod1でntt
			number_theoretic_transform<ntt_mod1, 137>(buf_num1.data(), new_size);
			number_theoretic_transform<ntt_mod1, 137>(buf1.data(), new_size);
			for (size_t i = 0; i < new_size; ++i)
			{
				buf1[i] = buf1[i] * buf_num1[i] % ntt_mod1;
			}
			number_theoretic_transform<ntt_mod1, 137, true>(buf1.data(), new_size);
			long long inv_size1 = _minv<ntt_mod1>(new_size);
			buf1.resize(result_size);
			for (size_t i = 0; i < result_size; ++i)
			{
				buf1[i] = buf1[i] * inv_size1 % ntt_mod1;
			}

			
			//コピー2回目
			size_t j;
			for (j = 0; j < dat.size(); ++j)
			{
				long long tmp = dat[j] % (long long)ntt_mod2;
				buf_num1[j] = tmp < 0 ? tmp + ntt_mod2 : tmp;
			}
			for (; j < new_size; ++j)
			{
				buf_num1[j] = 0;
			}
			for (j = 0; j < num.dat.size(); ++j)
			{
				long long tmp = num.dat[j] % (long long)ntt_mod2;
				buf2[j] = tmp < 0 ? tmp + ntt_mod2 : tmp;
			}
			//ntt_mod2でntt
			number_theoretic_transform<ntt_mod2, 136>(buf_num1.data(), new_size);
			number_theoretic_transform<ntt_mod2, 136>(buf2.data(), new_size);
			for (size_t i = 0; i < new_size; ++i)
			{
				buf2[i] = buf2[i] * buf_num1[i] % ntt_mod2;
			}
			number_theoretic_transform<ntt_mod2, 136, true>(buf2.data(), new_size);
			long long inv_size2 = _minv<ntt_mod2>(new_size);
			buf2.resize(result_size);
			for (size_t i = 0; i < result_size; ++i)
			{
				buf2[i] = buf2[i] * inv_size2 % ntt_mod2;
			}
			dat.resize(result_size);
			for (size_t i = 0; i < result_size; ++i)
			{
				dat[i] = _garner2<ntt_mod1, ntt_mod2>(buf1[i], buf2[i]);
			}
			fix_carry_force();
		}
		else
		{
			throw std::exception("border over.");
		}
	}
	ubigint& operator*=(const ubigint& _Right)
	{
		fix_carry_force();
		*this = multiply_naive(*this, _Right);
		fix_carry_force();
		return *this;
	}

	void random(size_t digits)
	{
		static std::random_device seed_gen;
		static std::mt19937_64 mt(seed_gen());
		static std::uniform_int_distribution<long long> dist(0, num_base - 1);
		dat.resize(digits);
		dat.shrink_to_fit();
		for (size_t i = 0; i < digits; ++i)
		{
			dat[i] = dist(mt);
		}
		max = num_base;
	}
};

#endif