#include"instant_timer.h"
#include<chrono>
#include<iostream>

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