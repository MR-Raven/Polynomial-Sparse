// Polynimial__sparse.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include<set>
#include <string>
#include<vector>

using namespace std;

template <typename T>
class Polynomial {
private:
	int degree;
	map<int, T> coef;
public:
	Polynomial(const vector<T>& coefficients) {
		degree = -1;
		for (auto num : coefficients) {
			++degree;
			if (num != T(0)) {
				coef[degree] = num;
			}
		}
		Cutting();
	}
	Polynomial(const T num = T(0)) {
		degree = -1;
		if (num != T(0)) {
			++degree;
			coef[degree] = num;
		}
	}
	template <class It>
	Polynomial(It first, It last) {
		degree = -1;
		T next_coef;
		while (first != last) {
			++degree;
			next_coef = *first;
			if (next_coef != T(0)) {
				coef[degree] = next_coef;
			}
			++first;
		}
		Cutting();
	}
	map<int, T> get_coef() const {
		return coef;
	}
	void print_coef() const {
		for (auto el : coef) {
			cout << el.second << "x^" << el.first << "  ";
		}
		cout << endl;
	}
	void Cutting() {
		while (degree > -1 && coef[degree] == T(0))
		{
			coef.erase(degree);
			--degree;
		}
		auto coef_copy = coef;
		for (auto el : coef) {
			if (el.second == T(0)) {
				coef_copy.erase(el.first);
			}
		}
		coef = coef_copy;
	}
	int Degree() const {
		return degree;
	}
	void normalization() {
		for (auto el : coef) {
			coef[el.first] = el.second / coef[degree];
		}
	}
	T module(T num) {
		if (num > T(0)) {
			return num;
		}
		else {
			return num * T(-1);
		}
	}
	Polynomial<T>& operator += (const Polynomial<T>& p) {
		degree = max(degree, p.Degree());
		for (auto el : p.get_coef()) {
			coef[el.first] += el.second;
		}
		Cutting();
		return *this;
	}
	Polynomial<T>& operator -= (const Polynomial<T>& p) {
		degree = max(degree, p.Degree());
		for (auto el : p.get_coef()) {
			coef[el.first] -= el.second;
		}
		Cutting();
		return *this;
	}
	Polynomial<T>& operator *= (const Polynomial& p) {
		vector<T> result_coef(degree + p.Degree() + 1);
		for (auto el1 : coef) {
			for (auto el2 : p.get_coef()) {
				result_coef[el1.first + el2.first] += el1.second * el2.second;
			}
		}
		Polynomial<T> result = Polynomial<T>(result_coef);
		result.Cutting();
		*this = result;
		return *this;
	}
	typename map<int, T>::const_iterator begin() const {
		return coef.begin();
	}
	typename map<int, T>::const_iterator end() const {
		return coef.end();
	}

	T operator [] (int pos) const {
		for (auto it = coef.begin(); it != coef.end(); ++it) {
			if ((*it).first == pos) {
				return (*it).second;
			}
		}
		return T(0);
	}
	T operator() (T x) const {
		T result = T(0);
		auto coef_copy = coef;
		int ind = this->Degree();
		while (ind >= 0) {
			result = result*x + coef_copy[ind];
			--ind;
		}
		return result;
	}

	bool operator == (const Polynomial<T>& p) const {
		return coef == p.get_coef();
	}
	bool operator != (const Polynomial<T>& p) const {
		return coef != p.get_coef();
	}

	friend Polynomial<T> operator + (Polynomial<T> p1, const Polynomial<T>& p2) {
		return p1 += p2;
	}
	friend Polynomial<T> operator - (Polynomial<T> p1, const Polynomial<T>& p2) {
		return p1 -= p2;
	}
	friend Polynomial<T> operator * (Polynomial<T> p1, const Polynomial<T>& p2) {
		return p1 *= p2;
	}
	friend Polynomial<T> operator & (const Polynomial<T> p1, const Polynomial<T>& p2) {
		Polynomial<T> tmp = Polynomial<T>(T(0));
		Polynomial<T> composition;
		for (int i = 0; i <= p1.Degree(); ++i) {
			tmp = Polynomial<T>(p1[i]);
			for (int j = 0; j < i; ++j) {
				tmp *= p2;
			}
			composition += tmp;
		}
		return composition;
	}
	friend Polynomial<T> operator / (const Polynomial<T>& p1, const Polynomial<T>& p2) {
		if (p2.Degree() != -1) {
			Polynomial<T> p = p1;
			int quotient_size = max(0, p.Degree() - p2.Degree() + 1);
			vector<T> quotient_coef(quotient_size);
			vector<T> monom(quotient_size, T(0));
			Polynomial<T> monomial = Polynomial<T>(monom);
			T monom_coef;
			int cur_quotient_ind;
			while (p.Degree() >= p2.Degree())
			{
				monom_coef = p[p.Degree()] / p2[p2.Degree()];
				cur_quotient_ind = p.Degree() - p2.Degree();
				quotient_coef[cur_quotient_ind] = monom_coef;
				monom[cur_quotient_ind] = monom_coef;
				monomial = Polynomial<T>(monom);
				p -= (monomial*p2);
				monom[cur_quotient_ind] = T(0);
			}
			return Polynomial<T>(quotient_coef);
		}
		else {
			return{};
		}
	}
	friend Polynomial<T> operator % (const Polynomial<T>& p1, const Polynomial<T>& p2) {
		return p1 - (p1 / p2)*p2;
	}
	friend Polynomial<T> operator , (const Polynomial<T>& p1, const Polynomial<T>& p2) {
		Polynomial<T> p1_copy = p1;
		Polynomial<T> p2_copy = p2;
		while (p1_copy.Degree() != -1 && p2_copy.Degree() != -1) {
			if (p1_copy.Degree() > p2_copy.Degree()) {
				p1_copy = p1_copy % p2_copy;
			}
			else if (p1_copy.Degree() < p2_copy.Degree()) {
				p2_copy = p2_copy % p1_copy;
			}
			else {
				if (p1_copy.module(p1_copy[p1_copy.Degree()]) > p2_copy.module(p2_copy[p2_copy.Degree()])) {
					p1_copy = p1_copy % p2_copy;
				}
				else {
					p2_copy = p2_copy % p1_copy;
				}
			}
		}
		if (p1_copy.Degree() == -1) {
			p2_copy.normalization();
			return p2_copy;
		}
		p1_copy.normalization();
		return p1_copy;
	}
};
template <typename T>
std::ostream& operator << (std::ostream& out, const Polynomial<T>& p) {
	Polynomial<T> p1 = p;
	p1.Cutting();
	if (p1.Degree() == -1) {
		out << T(0);
	}
	else {
		string plus_str;
		for (int i = p.Degree(); i >= 0; --i) {
			if (p[i] != T(0)) {
				if (p[i] > T(0)) {
					out << plus_str;;
				}
				if (i > 1) {
					if (p[i] * p[i] != T(1)) {
						out << p[i] << "*";
					}
					else if (p[i] == T(-1)) {
						out << "-";
					}
					out << "x^" << i;
				}
				else {
					if (i == 1) {
						if (p[i] * p[i] != T(1)) {
							out << p[i] << "*";
						}
						else if (p[i] == T(-1)) {
							out << "-";
						}
						out << "x";
					}
					else {
						out << p[i];
					}
				}
			}
			plus_str = "+";
		}
	}
	return out;
}


int main()
{
    return 0;
}

