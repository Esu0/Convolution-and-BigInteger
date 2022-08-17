#pragma once

#include"convolution.h"
#include<string>
#include<iostream>
#include<random>
#include<stdexcept>
#include<fstream>
#include<ios>
#include<iomanip>
#include<string.h>


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

template<unsigned long long num_base, bool multi_thread = false>
class ubigint
{
	static_assert(num_base >= 2, "進数が不正");
	static_assert(num_base <= 1ull << 60, "進数が大きすぎる");
private:
	std::vector<long long> dat;
	long long max = num_base;
	bool carry_fixed = true;

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
		carry_fixed = true;
	}

	void fix_carry()
	{
		if (max >= LLONG_MAX >> 1)fix_carry_force();
	}

	explicit ubigint(std::vector<long long>&& v) :dat(v)
	{}

	bool equal(const ubigint& _Right)const&
	{
		if (dat.size() != _Right.dat.size())return false;
		for (size_t i = 0; i < dat.size(); ++i)
		{
			if (dat[i] != _Right.dat[i])return false;
		}
		return true;
	}

	bool greater(const ubigint& _Right)const&
	{
		if (dat.size() > _Right.dat.size())return true;
		if (dat.size() < _Right.dat.size())return false;
		for (size_t i = dat.size(); i > 0;)
		{
			--i;
			if (dat[i] > _Right.dat[i])return true;
			if (dat < _Right.dat[i])return false;
		}
		return false;
	}

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
	static ubigint hex(const char* str)
	{
		if constexpr (is_pow2(num_base))
		{
			static constexpr unsigned char d = (unsigned char)_log2(num_base);
			unsigned long long len_str = strnlen_s(str, ubigint::max_size * d / 4);
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
	static ubigint dec(const char* str)
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
	ubigint(const ubigint& _Right): dat(_Right.dat), max(_Right.max), carry_fixed(_Right.carry_fixed)
	{}

	//ムーブ代入コンストラクタ
	ubigint(ubigint&& _Right)noexcept : dat(std::move(_Right.dat)), max(_Right.max), carry_fixed(_Right.carry_fixed)
	{}

	//コピー代入演算子
	ubigint& operator=(const ubigint& _Right)
	{
		dat = _Right.dat;
		max = _Right.max;
		carry_fixed = _Right.carry_fixed;
		return *this;
	}

	//ムーブ代入演算子
	ubigint& operator=(ubigint&& _Right)noexcept
	{
		dat = std::move(_Right.dat);
		max = _Right.max;
		carry_fixed = _Right.carry_fixed;
		return *this;
	}

	void dump()
	{
		if(!carry_fixed)fix_carry_force();
		for (size_t i = dat.size(); i > 0;)
		{
			--i;
			std::cout << dat[i] << ' ';
		}
		std::cout << std::endl;
	}

	void dump(const char* fpath, const char* str)
	{
		if(!carry_fixed)fix_carry_force();
		std::ofstream ofs(fpath);
		for (size_t i = dat.size(); i > 0;)
		{
			--i;
			ofs << dat[i] << str;
		}
		ofs << std::endl;
	}

	void dump(const char* str)
	{
		if(!carry_fixed)fix_carry_force();
		for (size_t i = dat.size(); i > 0;)
		{
			--i;
			std::cout << dat[i] << str;
		}
		std::cout << std::endl;
	}

	[[nodiscard]]
	bool operator==(const ubigint& _Right)const&
	{
		if (_Right.carry_fixed)
		{
			if (carry_fixed)
			{
				return equal(_Right);
			}
			else
			{
				ubigint tmp(*this);
				tmp.fix_carry_force();
				return tmp.equal(_Right);
			}
		}
		else
		{
			ubigint tmp(_Right);
			tmp.fix_carry_force();
			if (carry_fixed)
			{
				return equal(tmp);
			}
			else
			{
				ubigint tmp2(*this);
				tmp2.fix_carry_force();
				return tmp2.equal(tmp);
			}
		}
	}

	[[nodiscard]]
	bool operator>(const ubigint& _Right)const&
	{
		if (_Right.carry_fixed)
		{
			if (carry_fixed)
			{
				return greater(_Right);
			}
			else
			{
				ubigint tmp(*this);
				tmp.fix_carry_force();
				return tmp.greater(_Right);
			}
		}
		else
		{
			ubigint tmp(_Right);
			tmp.fix_carry_force();
			if (carry_fixed)
			{
				return greater(tmp);
			}
			else
			{
				ubigint tmp2(*this);
				tmp2.fix_carry_force();
				return tmp2.greater(tmp);
			}
		}
	}

	[[nodiscard]]
	bool operator<(const ubigint& _Right)const&
	{
		return _Right.operator>(*this);
	}

	[[nodiscard]]
	bool operator>=(const ubigint& _Right)const&
	{
		return !(operator<(_Right));
	}

	[[nodiscard]]
	bool operator<=(const ubigint& _Right)const&
	{
		return !(operator>(_Right));
	}

	[[nodiscard]]
	bool operator!=(const ubigint& _Right)const&
	{
		return !operator==(_Right);
	}

	//前置インクリメント
	ubigint& operator++()
	{
		++dat[0];
		++max;
		carry_fixed = false;
		fix_carry();
		return *this;
	}

	//後置インクリメント
	ubigint operator++(int)
	{
		ubigint tmp = *this;
		this->operator++();
		return tmp;
	}

	//前置デクリメント
	ubigint& operator--()
	{
		--dat[0];
		++max;
		carry_fixed = false;
		fix_carry();
		return *this;
	}

	//後置デクリメント
	ubigint operator--(int)
	{
		ubigint tmp = *this;
		this->operator--();
		return tmp;
	}

	//複合代入(足し算)
	ubigint& operator+=(const ubigint& _Right)
	{
		if (_Right.dat.size() > dat.size())dat.resize(_Right.dat.size(), 0);
		for (size_t i = 0; i < _Right.dat.size(); ++i)dat[i] += _Right.dat[i];
		max += _Right.max;
		carry_fixed = false;
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
		carry_fixed = false;
		fix_carry();
		return *this;
	}

	//Right <- Left - Right
	void sub_reverse(ubigint& _Right)const&
	{
		if (dat.size() > _Right.dat.size())_Right.dat.resize(dat.size(), 0);
		for (size_t i = 0; i < dat.size(); ++i)_Right.dat[i] = dat[i] - _Right.dat[i];
		_Right.max += max;
		carry_fixed = false;
		_Right.fix_carry();
	}

	//引き算
	[[nodiscard]]
	ubigint operator-(ubigint _Right)const&
	{
		sub_reverse(_Right);
		return _Right;
	}

	//普通に掛け算
	[[nodiscard]]
	ubigint multiply_naive(const ubigint& num)const&
	{
		std::vector<long long> resultv(dat.size() + num.dat.size());
		multiply_naive_unsafe<num_base>(dat.data(), num.dat.data(), resultv.data(), (unsigned long)dat.size(), (unsigned long)num.dat.size());
		return ubigint(std::move(resultv));
	}

	//数論変換を使って掛け算
	[[nodiscard]]
	ubigint multiply_ntt(const ubigint& num)const&
	{
		static constexpr unsigned long ntt_mod1 = 2013265921;
		static constexpr unsigned long ntt_mod2 = 1811939329;
		static constexpr unsigned long ntt_mod3 = 469762049;
		static constexpr unsigned long g1 = 31;
		static constexpr unsigned long g2 = 13;
		static constexpr unsigned long g3 = 3;
		static constexpr unsigned long long m1m2 = (unsigned long long)ntt_mod1 * ntt_mod2;
		static constexpr unsigned long border1 = (ntt_mod1 - 1) / (num_base - 1) / (num_base - 1) + 1;
		static constexpr unsigned long border2 = (unsigned long)((m1m2 - 1) / (num_base - 1) / (num_base - 1) + 1);
		static constexpr unsigned long long m1m2_div0 = m1m2 % (num_base - 1);
		static constexpr unsigned long long m1m2_div1 = (m1m2 / (num_base - 1)) % (num_base - 1);
		static constexpr unsigned long long m1m2_div2 = m1m2 / (num_base - 1) / (num_base - 1);
		static constexpr unsigned long long border3 = (m1m2_div2 / (num_base - 1) == 0) ? (m1m2_div2 * ntt_mod3 + m1m2_div1 * ntt_mod3 / (num_base - 1) + m1m2_div0 * ntt_mod3 / (num_base - 1) / (num_base - 1) + 1) : ULLONG_MAX;
		unsigned long result_size = (unsigned long)(dat.size() + num.dat.size());
		unsigned long n = (unsigned long)(std::min)(dat.size(), num.dat.size());
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
			convolution_mod<ntt_mod1, g1, multi_thread>(buf1.data(), buf_num1.data(), result_size, new_size);
			ubigint prod(std::move(buf1));
			prod.fix_carry_force();
			return prod;
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
			convolution_mod<ntt_mod1, g1, multi_thread>(buf1.data(), buf_num1.data(), result_size, new_size);

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
			convolution_mod<ntt_mod2, g2, multi_thread>(buf2.data(), buf_num1.data(), result_size, new_size);


			buf_num1.resize(result_size);
			for (size_t i = 0; i < result_size; ++i)
			{
				buf_num1[i] = _garner2<ntt_mod1, ntt_mod2>(buf1[i], buf2[i]);
			}
			ubigint prod(std::move(buf_num1));
			prod.fix_carry_force();
			return prod;
		}
		else if (n < border3)
		{
			unsigned long new_size = (unsigned long)ceil_pow2(result_size);
			std::vector<long long> buf1(new_size);
			std::vector<long long> buf_num1(new_size);
			std::vector<long long> buf2(new_size);
			std::vector<long long> buf3(new_size);

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
			
			//ntt_mod1で畳み込み
			convolution_mod<ntt_mod1, g1, multi_thread>(buf1.data(), buf_num1.data(), result_size, new_size);

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

			//ntt_mod2で畳み込み
			convolution_mod<ntt_mod2, g2, multi_thread>(buf2.data(), buf_num1.data(), result_size, new_size);

			//コピー3回目
			for (j = 0; j < dat.size(); ++j)
			{
				long long tmp = dat[j] % (long long)ntt_mod3;
				buf_num1[j] = tmp < 0 ? tmp + ntt_mod3 : tmp;
			}
			for (; j < new_size; ++j)
			{
				buf_num1[j] = 0;
			}
			for (j = 0; j < num.dat.size(); ++j)
			{
				long long tmp = num.dat[j] % (long long)ntt_mod3;
				buf3[j] = tmp < 0 ? tmp + ntt_mod2 : tmp;
			}

			//ntt_mod3で畳み込み
			convolution_mod<ntt_mod3, g3, multi_thread>(buf3.data(), buf_num1.data(), result_size, new_size);
			static constexpr unsigned char digit_size = _log_ceil<num_base>(1 << 25) + (unsigned char)2;
			buf_num1 = std::vector<long long>((size_t)result_size + digit_size - 1, 0);

			unsigned long long dig_1[digit_size] = {};
			for (size_t i = 0; i < result_size; ++i)
			{
				_garner3<num_base, ntt_mod1, ntt_mod2, ntt_mod3, digit_size>(buf1[i], buf2[i], buf3[i], dig_1);
				for (size_t j = 0; j < digit_size; ++j)
				{
					buf_num1[i + j] += dig_1[j];
				}
			}

			ubigint prod(std::move(buf_num1));
			prod.fix_carry_force();
			return prod;
		}
		else
		{
			throw std::out_of_range("掛け算失敗");
			return ubigint();
		}
	}

	ubigint& operator*=(ubigint& _Right)
	{
		if (!carry_fixed)fix_carry_force();
		if (!_Right.carry_fixed)_Right.fix_carry_force();
		unsigned long long ntt_size = ceil_pow2(_Right.dat.size() + dat.size());
		if (_Right.dat.size() * dat.size() < ntt_size * ((unsigned long long)_log2(ntt_size) * 3 + 1))
		{
			*this = multiply_naive(_Right);
			delete_zero();
		}
		else
		{
			*this = multiply_ntt(_Right);
		}
		fix_carry_force();
		return *this;
	}

	ubigint& operator*=(const ubigint& _Right)
	{
		if(!carry_fixed)fix_carry_force();
		if (_Right.carry_fixed)
		{
			unsigned long long ntt_size = ceil_pow2(_Right.dat.size() + dat.size());
			if (_Right.dat.size() * dat.size() < ntt_size * ((unsigned long long)_log2(ntt_size) * 3 + 1))
			{
				*this = multiply_naive(_Right);
				delete_zero();
			}
			else
			{
				*this = multiply_ntt(_Right);
			}
			fix_carry_force();
		}
		else
		{
			ubigint buf(_Right);
			operator*=(buf);
		}
		return *this;

	}

	[[nodiscard]]
	ubigint operator*(const ubigint& _Right)const&
	{
		return ubigint(*this) *= _Right;
	}

	ubigint& operator>>=(unsigned long diff)
	{
		if (dat.size() <= diff)*this = ubigint();
		else
		{
			dat.erase(dat.cbegin(), dat.cbegin() + diff);
		}
		return *this;
	}

	ubigint operator>>(unsigned long diff)const&
	{
		return ubigint(*this) >>= diff;
	}

	ubigint& operator<<=(unsigned long diff)
	{
		dat.insert(dat.cbegin(), diff, 0);
		return *this;
	}

	ubigint operator<<(unsigned long diff)const&
	{
		return ubigint(*this) <<= diff;
	}
	//逆数
	ubigint inverse(unsigned long seido)const&
	{
		ubigint x(num_base / dat.back());
		//Newton-Raphson method
		return x;
	}

	ubigint operator/(const ubigint& _Right)const&
	{
		if (_Right.dat.size() > dat.size())return ubigint();
		if (_Right.dat.size() == 1)
		{
			//simple method
		}
		else
		{
			//Newton-Raphson method
		}
		return ubigint();
	}

	ubigint& operator/=(const ubigint& _Right)
	{
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
		carry_fixed = true;
	}
};
