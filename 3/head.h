#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>

using namespace std;


// float || double
typedef double real;
typedef std::vector<real> vec;


// ��������� �� ���������
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


// �������� ��������
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


// ��������� ��������
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


// ��������� �� ���������
inline vec operator*(const vec& a, double b) {
	vec result = a;
	for (int i = 0; i < result.size(); i++)
		result[i] *= b;
	return result;
}


// ��������� �� ���������
inline vec operator*(double b, const vec& a) {
	return operator*(a, b);
}


// ��������� ������������
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


// ��������� �����
inline std::ostream& operator<<(std::ostream& out, const vec& v) {
	for (int i = 0; i < v.size() - 1; ++i)
		out << v[i] << "\t";
	out << v.back();
	return out;
}