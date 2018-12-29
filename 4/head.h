#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>
#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;


// float || double
typedef double real;
typedef std::vector<real> vec;
typedef std::vector<vector<real>> vec2;

// Умножение на константу
inline bool operator==(const vec& a, const vec& b) {
#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
#endif
	for (int i = 0; i < a.size(); ++i)
		if (a[i] != b[i])
			return false;

	return true;
}


// Сложение векторов
inline vec operator+(const vec& a, const vec& b) {
#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
#endif
	vec result = a;
	for (int i = 0; i < b.size(); i++)
		result[i] += b[i];
	return result;
}


// Вычитание векторов
inline vec operator-(const vec& a, const vec& b) {
#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
#endif
	vec result = a;
	for (int i = 0; i < b.size(); i++)
		result[i] -= b[i];
	return result;
}



// Умножение на константу
inline vec operator*(const vec& a, double b) {
	vec result = a;
	for (int i = 0; i < result.size(); i++)
		result[i] *= b;
	return result;
}


// Умножение на константу
inline vec operator*(double b, const vec& a) {
	return operator*(a, b);
}


// Скалярное произведение
inline real operator*(const vec& a, const vec& b) {
#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
#endif
	real sum = 0;
	for (int i = 0; i < a.size(); i++)
		sum += a[i] * b[i];
	return sum;
}


// Потоковый вывод
inline std::ostream& operator<<(std::ostream& fout, const vec& v) {
	for (int i = 0; i < v.size() - 1; ++i)
		fout << v[i] << "\t";
	fout << v.back();
	return fout;
}