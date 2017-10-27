#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <cstdarg>

using namespace std;

#define M_PI       3.1415926f
#define M_PI_2     1.5707963f
#define M_PI_4     0.7853981f

#define PRECISION			0.0000001f
#define SQRT(X)				sqrt((X))
#define SQR(X)				((X) * (X))
#define EQ(X, Y, EPS)		(ABS((X) - (Y)) < EPS)
#define ABS(X)				(((X) > 0.0f) ? (X) : (-(X)))
#define SWAP(type, x, y)	{ type temp = (x); (x) = (y); (y) = temp; }
#define MIN(x, y)			((x) > (y) ? (y) : (x))
#define MAX(x, y)			((x) > (y) ? (x) : (y))
#define ARR_ZERO(A, N)		memset((A), 0, sizeof(A[0]) * (N))
#define ARR_COPY(D, S, N)	memmove((D), (S), sizeof(S[0]) * (N))
#define EQ_ZERO(X, EPS)		(ABS(X) < EPS)
#define DEG2RAD(X)	((X) * (M_PI) / (180.0f))
#define RAD2DEG(X)	((X) * (180.0f) / (M_PI))

class GVector3;
class GVector;
class GMatrix;
class GLine;
class GPlane;
class GQuater;

typedef GVector3 GPoint3;


class GVector3
{
public:

	GVector3(float x = 0.0f, float y = 0.0f, float z = 0.0f);
	GVector3(const GVector3 &copy);

	// Assignment operator overloading
	GVector3& operator =(const GVector3& rhs);

	// Compound assignment operator overloading
	GVector3 &operator +=(const GVector3 &rhs);
	GVector3 &operator -=(const GVector3 &rhs);
	GVector3 &operator *=(const float &s);
	GVector3 &operator /=(const float &s);
	GVector3 &operator ^=(const GVector3 &rhs);

	
	// Equality and inequality operator overloading
	bool operator ==(const GVector3 &rhs) const;
	bool operator !=(const GVector3 &rhs) const;

	GVector3 operator +() const;
	GVector3 operator -() const;

	// The addition and subtraction of vectors
	GVector3 operator +(const GVector3& rhs) const; // u + v
	GVector3 operator -(const GVector3& rhs) const; // u - v

	// The inner product of vectors
	float operator *(const GVector3& rhs) const; // u * v

	// The outer product of vectors
	GVector3 operator ^(const GVector3& rhs) const; // u ^ v 

	// The norm and direction of vectors
	friend float norm(const GVector3& v);
	GVector3& normalize();

	GVector3& Set(const float& x, const float& y, const float& z);

	// The projection of vectors
	friend GVector3 proj(const GVector3& v, const GVector3& w); // proj 
	friend GVector3 perp(const GVector3& v, const GVector3& w); // prep

	// The output of vectors
	friend ostream& operator <<(ostream &os, const GVector3& v);

	// Subscript operator overloading
	float &operator [](const int &idx);
	const float &operator [](const int &idx) const;

	// The multiplication of vectors and scalars
	friend GVector3 operator *(const GVector3& lhs, const float& k); // k * u
	friend GVector3 operator *(const float& k, const GVector3& rhs); // u * k
	friend GVector3 operator /(const GVector3& lhs, const float& k); // u / k

	friend float dist(const GVector3 &v, const GVector3 &w);

private:
	float V[3];

};

class GVector
{
public:
	// Constructors and destructor
	GVector(int dim = 3);
	GVector(int dim, float x, ...);
	GVector(const GVector3 &copy);
	GVector(const GVector &copy);
	~GVector();

	// Assignment operator overloading
	GVector &operator =(const GVector &rhs);

	// Compound assignment operator overloading
	GVector &operator +=(const GVector &rhs);
	GVector &operator -=(const GVector &rhs);
	GVector &operator *=(const float &s);
	GVector	&operator /=(const float &s);

	// Unary operator overloading
	GVector operator +() const;
	GVector operator -() const;

	// Arithmetic operator overloading
	GVector operator +(const GVector &rhs) const;
	GVector operator -(const GVector &rhs) const;
	float  operator *(const GVector &rsh) const;
	GVector operator /(const float &s) const;

	// Equality and inequality operator overloading
	bool operator ==(const GVector &rhs) const;
	bool operator !=(const GVector &rhs) const;

	// Subscript operator overloading
	float &operator [](const int &idx);
	const float &operator [](const int &idx) const;

	// Member functions
	GVector	&Set(float x, ...);
	GVector &Set(float *pVal);
	GVector	&Normalize();
	GVector &SetZeros();
	int GetDim() const;

	friend GVector operator *(const float &s, const GVector &rhs);
	friend GVector operator *(const GVector &lhs, const float &s);
	friend GVector operator *(const GMatrix &lhs, const GVector &rhs);
	friend GMatrix operator *(const GVector &lhs, const GMatrix &rhs);
	friend float norm(const GVector &v);
	friend float dist(const GVector &v, const GVector &w);
	friend ostream &operator <<(ostream &os, const GVector &v);

	friend class GMatrix;

private:
	int N;
	float *V;
};

class GMatrix
{
public:
	// Constructors and destructor
	GMatrix(int row = 4, int col = 4, float *elem = NULL);
	GMatrix(const GMatrix &copy);
	~GMatrix();

	// Assignment operator overloading
	GMatrix &operator =(const GMatrix &rhs);

	// Compound assignment operator overloading
	GMatrix &operator +=(const GMatrix &rhs);
	GMatrix &operator -=(const GMatrix &rhs);
	GMatrix &operator *=(const GMatrix &rhs);
	GMatrix &operator *=(const float &s);
	GMatrix &operator /=(const float &s);

	// Unary operator overloading
	GMatrix operator +() const;
	GMatrix operator -() const;

	// Arithmetic operator overloading
	GMatrix operator +(const GMatrix &rhs) const;
	GMatrix operator -(const GMatrix &rhs) const;
	GMatrix operator *(const GMatrix &rhs) const;
	GMatrix operator /(const float &s) const;

	// Equality and inequality operator overloading
	bool operator ==(const GMatrix &rhs) const;
	bool operator !=(const GMatrix &rhs) const;

	// Subscript operator overloading
	float *operator [](const int idx);
	const float *operator [](const int idx) const;

	GMatrix &SetTranspose();
	GMatrix &SetIdentity();
	GMatrix &SetZeros();
	GMatrix &SetRowVec(const int idx, const GVector &v);
	GMatrix &SetColVec(const int idx, const GVector &v);
	GMatrix &ExchangeRows(const int idx0, const int idx1);
	GMatrix &ExchangeCols(const int idx0, const int idx1);

	int GetRowNum() const;
	int GetColNum() const;
	GVector GetRowVec(const int idx) const;
	GVector GetColVec(const int idx) const;

	bool IsSquare() const;

	friend GVector operator *(const GMatrix &lhs, const GVector &rhs);
	friend GMatrix operator *(const GVector& lhs, const GMatrix &rhs);
	friend GMatrix operator *(const GMatrix &lhs, const float &s);
	friend GMatrix operator *(const float &s, const GMatrix &rhs);
	friend ostream &operator <<(ostream &os, const GMatrix &m);
	friend GMatrix RowEchelonForm(const GMatrix &m);
	friend GMatrix ReducedRowEchelon(const GMatrix &m);
	friend float *form_arr(const GMatrix &m);
	friend int Rank(const GMatrix &m);
	friend int Nullity(const GMatrix &m);

private:
	// Data members
	int r;		// Dimension of row vector
	int c;		// Dimension of column vector
	float *M;	// Array of elements of a matrix
};

class GLine
{
public:
	
	GLine(const GPoint3 &_p = GPoint3(0, 0, 0), const GVector3 &_v = GVector3(0, 0, 0));
	GLine(const GLine &copy);

	GLine &operator =(const GLine &rhs);
	GPoint3 operator ()(const float t) const;

	GLine &SetPt(const GPoint3 &_p);
	GLine &SetDir(const GVector3 &_v);

	GPoint3 GetPt() const;
	GVector3 GetDir() const;

	bool IsOnLine(const GPoint3 &q) const;

	friend ostream &operator <<(ostream &os, const GLine &l);
	friend float dist(const GPoint3 &q, const GLine &l);
	friend bool intersect_line_plane(GPoint3 &p, const GLine &l, const GPlane &pi);
	friend bool intersect_line_triangle(GPoint3 &q, const GLine &l, const GPoint3 &p1, const GPoint3 &p2, const GPoint3 &p3);
	

private:

	GPoint3 p;
	GVector3 v;		// line function : l(t) = p + v * t

};

class GPlane
{
public:
	GPlane(const GVector3 &_n = GVector3(0.0, 0.0, 1.0), const GPoint3 &_p = GPoint3());
	GPlane(const GPoint3 &p1, const GPoint3 &p2, const GPoint3 &p3);
	GPlane(const float &a, const float &b, const float &c, const float &d);
	GPlane(const GPlane &copy);

	GPlane &operator =(const GPlane &rhs);

	GVector3 GetNormal() const;

	bool IsOnPlane(const GPoint3 &p) const;
	bool IsAbovePlane(const GPoint3 &p) const;
	bool IsBelowPlane(const GPoint3 &p) const;

	friend ostream &operator <<(ostream &os, const GPlane &pi);
	friend float dist(const GPlane &pi, const GPoint3 &p);
	friend bool intersect_line_plane(GPoint3 &p, const GLine &l, const GPlane &pi);
	friend bool intersect_line_triangle(GPoint3 &q, const GLine &l, const GPoint3 &p1, const GPoint3 &p2, const GPoint3 &p3);

private:
	GVector3 n;
	float d;
};

class GQuater
{
	enum EulerType
	{
		EULER_XYZ = 0,
		EULER_ZYX = 1,
	};

public:
	GQuater(float w = 1.0, float x = 0.0, float y = 0.0, float z = 0.0);
	GQuater(const GQuater &copy);
	GQuater(const float *q, const bool invOrder = false);
	GQuater(GVector3 axis, float theta, bool radian = false);
	GQuater(float theta1, float theta2, float theta3, EulerType eulerType = EULER_XYZ);
	
	GQuater &operator =(const GQuater &rhs);

	GQuater &operator +=(const GQuater &rhs);
	GQuater &operator -=(const GQuater &rhs);
	GQuater &operator *=(const GQuater &rhs);
	GQuater &operator /=(const GQuater &rhs);
	GQuater &operator *=(const float s);
	GQuater &operator /=(const float s);

	GQuater operator +() const;
	GQuater operator -() const;

	GQuater operator +(const GQuater &rhs) const;
	GQuater operator -(const GQuater &rhs) const;
	GQuater operator *(const GQuater &rhs) const;
	GQuater operator /(const GQuater &rhs) const;
	GQuater operator /(const float s) const;
	GVector3 operator *(const GVector3 &rhs) const;

	bool operator ==(const GQuater &rhs) const;
	bool operator !=(const GQuater &rhs) const;

	GQuater &Set(const float w, const float x, const float y, const float z);
	GQuater &Set(float *q, bool invOrder = false);

	GQuater &SetIdentity();
	GQuater &SetConjugate();
	GQuater &SetInverse();
	GQuater &Normalize();

	GQuater &SetFromAngleAxis(const float theta, GVector3 axis, bool radian = false);
	GQuater &SetFromEulerAngle(float theta1, float theta2, float theta3, EulerType eulerType = EULER_XYZ);
	GQuater &SetFromMatrix(const GMatrix &m);
	

	bool IsUnitQuater() const;
	bool IsIdentity() const;

	friend GQuater operator *(const GQuater &lhs, const float &s);
	friend GQuater operator *(const float &s, const GQuater &rhs);
	friend ostream &operator <<(ostream &os, const GQuater &q);
	friend float norm(const GQuater &q);
	friend GQuater inv(const GQuater &q);

private:
	
	float W;
	float X, Y, Z;
};