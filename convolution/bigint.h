#pragma once

#ifndef _BIGINT_H_
#define _BIGINT_H_

#include"convolution.h"

template<unsigned long long num_base>
class ubigint
{
private:
	std::vector<long long> dat;

	void fix_carry()
	{
		long long carry = 0;
		long long* p = dat.data();
		for (size_t i = 0; i < dat.size(); ++i)
		{
			p[i] += carry;
			carry = p[i] / num_base;
			p[i] %= num_base;
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
	}
public:
	ubigint(): dat(1,0)
	{}

	ubigint(long long num) : dat(1, num)
	{
		fix_carry();
	}

	//��ϊ�10->�C��
	static ubigint<num_base> hex(const char* str)
	{
		if constexpr (is_pow2(num_base))
		{
			return ubigint();
		}
		else
		{
			return ubigint();
		}
	}

	//��ϊ�16->�C��
	static ubigint<num_base> dec(const char* str)
	{
		return ubigint();
	}

	//�R�s�[����R���X�g���N�^
	ubigint(const ubigint& _Right): dat(_Right.dat)
	{}

	//���[�u����R���X�g���N�^
	ubigint(ubigint&& _Right)noexcept : dat(std::move(_Right.dat))
	{}

	//�R�s�[������Z�q
	ubigint& operator=(const ubigint _Right)
	{
		dat = _Right.dat;
	}

	//���[�u������Z�q
	ubigint& operator=(ubigint&& _Right)noexcept
	{
		dat = std::move(_Right.dat);
	}

	//�O�u�C���N�������g
	ubigint& operator++()
	{
		for (size_t i = 0; i < dat.size(); ++i)
		{
			++(dat[i]);
			if (dat[i] < num_base)return *this;
			else dat[i] = 0;
		}
		dat.push_back(1);
	}

	//��u�C���N�������g
	ubigint operator++(int)
	{
		ubigint<num_base> tmp = *this;
		this->operator++();
		return tmp;
	}

	//�O�u�f�N�������g
	ubigint& operator--()
	{
		if (dat.size() == 1 && dat[0] == 0)
		{
			dat[0] = LLONG_MIN;
			return *this;
		}
		for (size_t i = 0; i < dat.size(); ++i)
		{
			if (!(dat[i]))
			{
				dat[i] = LINT_BASE - 1;
			}
			else
			{
				--(dat[i]);
				return *this;
			}
		}
		return *this;
	}

	//��u�f�N�������g
	ubigint operator--(int)
	{
		ubigint<num_base> tmp = *this;
		this->operator--();
		return tmp;
	}

	//�������(�����Z)
	ubigint& operator+=(const ubigint& _Right)
	{
		if (_Right.dat.size() > dat.size())dat.resize(_Right.dat.size(), 0);
		for (size_t i = 0; i < _Right.dat.size(); ++i)dat[i] += _Right.dat[i];
		fix_carry();
	}

	//�����Z
	[[nodiscard]]
	ubigint operator+(ubigint _Right)const&
	{
		_Right.operator+=(*this);
		return _Right;
	}

	//�������(�����Z)
	ubigint& operator-=(const ubigint& _Right)
	{
		if (_Right.dat.size() > dat.size())dat.resize(_Right.dat.size(), 0);
		for (size_t i = 0; i < _Right.dat.size(); ++i)dat[i] -= _Right.dat[i];
		fix_carry();
	}

	//Right <- Left - Right
	void sub_reverse(ubigint& _Right)const&
	{
		if (dat.size() > _Right.dat.size())_Right.dat.resize(dat.size(), 0);
		for (size_t i = 0; i < dat.size(); ++i)_Right.dat[i] = dat[i] - _Right.dat[i];
		_Right.fix_carry();
	}

	//�����Z
	[[nodiscard]]
	ubigint operator-(ubigint _Right)const&
	{
		sub_reverse(_Right);
		return _Right;
	}

	ubigint& operator*=(const ubigint& _Right)
	{

	}
};

#endif