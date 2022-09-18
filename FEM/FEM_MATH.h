#pragma once
#include <vector>


namespace MATH
{
	const double E = 2.71828182845904523536;
	const double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844;
	const double EPSILON = 1.0e-12;

	template <typename type, int size>    class Vec;

#define FORCEINLINE __forceinline

	template <typename type>
	FORCEINLINE type rcp(type v) {
		return type(1.0) / v;
	}

	enum class Uninitialized { Type };

	template<typename type, int size>
	class Vec {
	};

	typedef Vec <double, 2> dVec2;

	template<typename type>
	class Vec<type, 2> {

	public:

		union {
			struct {
				type x;
				type y;
			};
			struct {
				type arr[2];
			};
		};

		FORCEINLINE Vec();
		FORCEINLINE Vec(Uninitialized);
		FORCEINLINE Vec(type _x, type _y);
		FORCEINLINE Vec(type v);
		FORCEINLINE Vec(const type * v);
		FORCEINLINE type & operator[](int i);
		FORCEINLINE const type & operator[](int i) const;
		FORCEINLINE type distance(const Vec & v) const;
		FORCEINLINE type length2() const;
		FORCEINLINE type length() const;
		FORCEINLINE Vec & normalize();
		FORCEINLINE type dot(const Vec &v0) const;
		FORCEINLINE Vec operator-() const;
		FORCEINLINE Vec operator-(const Vec & v) const;
		FORCEINLINE Vec operator*(const Vec & v) const;
		FORCEINLINE Vec operator+(const Vec & v) const;
		FORCEINLINE Vec & operator*=(type v);
		FORCEINLINE Vec & operator/=(type v);
		FORCEINLINE Vec operator/(type v) const;
		FORCEINLINE Vec & operator*=(const Vec & v);
		FORCEINLINE Vec & operator/=(const Vec & v);
		FORCEINLINE Vec & operator+=(const Vec & v);
		FORCEINLINE Vec & operator-=(const Vec & v);
		FORCEINLINE Vec operator/(const Vec & v) const;
	};

	template <typename type>
	Vec<type, 2>::Vec() : x(0), y(0) {
	}

	template <typename type>
	Vec<type, 2>::Vec(Uninitialized) {
	}

	template <typename type>
	Vec<type, 2>::Vec(type _x, type _y) : x(_x), y(_y) {
	}

	template <typename type>
	Vec<type, 2>::Vec(const type * v) : x(v[0]), y(v[1]) {
	}

	template <typename type>
	Vec<type, 2>::Vec(type v) : x(v), y(v) {
	}

	template <typename type>
	type Vec<type, 2>::length2() const {
		return x * x + y * y;
	}

	template <typename type>
	type Vec<type, 2>::length() const {
		return (type)sqrt(x * x + y * y);
	}

	template <typename type>
	Vec<type, 2> & Vec<type, 2>::normalize() {
		type ilength = rsqrt(x * x + y * y);
		x *= ilength;
		y *= ilength;
		return *this;
	}

	template <typename type>
	type Vec<type, 2>::dot(const Vec<type, 2> &v0) const {
		return v0.x * x + v0.y * y;
	}

	template <typename type>
	const type & Vec<type, 2>::operator[](int i) const {
		return arr[i];
	}

	template <typename type>
	type & Vec<type, 2>::operator[](int i) {
		return arr[i];
	}

	template <typename type>
	Vec<type, 2> Vec<type, 2>::operator-() const {
		return{ -x, -y };
	}

	template <typename type>
	Vec<type, 2>    Vec<type, 2>::operator-(const Vec<type, 2> & v) const {
		return{
			x - v.x,
			y - v.y
		};
	}

	template <typename type>
	Vec<type, 2> & Vec<type, 2>::operator*=(type v) {
		x *= v;
		y *= v;
		return *this;
	}

	template <typename type>
	Vec<type, 2> & Vec<type, 2>::operator*=(const Vec<type, 2> & v) {
		x *= v.x;
		y *= v.y;
		return *this;
	}

	template <typename type>
	Vec<type, 2> Vec<type, 2>::operator*(const Vec<type, 2> & v) const {
		return{
			x * v.x,
			y * v.y
		};
	}

	template <typename type>
	Vec<type, 2> & Vec<type, 2>::operator/=(type v) {
		type iv = (type)rcp(v);
		x *= iv;
		y *= iv;
		return *this;
	}

	template <typename type>
	Vec<type, 2> & Vec<type, 2>::operator/=(const Vec<type, 2> & v) {
		x /= v.x;
		y /= v.y;
		return *this;
	}

	template <typename type>
	Vec <type, 2> Vec <type, 2>::operator/(const Vec <type, 2> & v) const {
		return{
			x / v.x,
			y / v.y
		};
	}

	template <typename type>
	Vec <type, 2> Vec <type, 2>::operator/(type v) const {
		return{
			x / v,
			y / v
		};
	}

	template <typename type>
	Vec<type, 2> & Vec<type, 2>::operator+=(const Vec<type, 2> & v) {
		x += v.x;
		y += v.y;
		return *this;
	}

	template <typename type>
	Vec<type, 2> Vec<type, 2>::operator+(const Vec<type, 2> & v) const {
		return{
			x + v.x,
			y + v.y
		};
	}

	template <typename type>
	Vec<type, 2> & Vec<type, 2>::operator-=(const Vec<type, 2> & v) {
		x -= v.x;
		y -= v.y;
		return *this;
	}

	template <typename type>
	type Vec<type, 2>::distance(const Vec<type, 2> & v) const {
		type dx = x - v.x;
		type dy = y - v.y;
		return (type)sqrt(dx * dx + dy * dy);
	}

	const dVec2 dVec2_zero(0.0);

	FORCEINLINE dVec2 normalize(const dVec2 & v) {
		dVec2 ret;
		double ilength = sqrt(v.x * v.x + v.y * v.y);
		ret.x = v.x / ilength;
		ret.y = v.y / ilength;
		return ret;
	}

	FORCEINLINE dVec2 normal(const dVec2 & v) {
		return dVec2(v.y, -v.x);
	}

	FORCEINLINE double dot(const dVec2 &a, const dVec2 &b) { return a.x * b.x + a.y * b.y; }

	FORCEINLINE double distance(const dVec2 & v, const dVec2 & v1) {
		double dx = v.x - v1.x;
		double dy = v.y - v1.y;
		return sqrt(dx * dx + dy * dy);
	}
	FORCEINLINE double cross(const dVec2 &v0, const dVec2 &v1) {
		return v0.x * v1.y - v0.y * v1.x;
	}

	FORCEINLINE double det2(double a1, double b1, double a2, double b2)
	{
		return (a1*b2) - (b1*a2);
	}
	FORCEINLINE double det3(double a1, double b1, double c1, double a2, double b2, double c2, double a3, double b3, double c3)
	{
		return a1*det2(b2, c2, b3, c3) -
			b1*det2(a2, c2, a3, c3) +
			c1*det2(a2, b2, a3, b3);
	}
}