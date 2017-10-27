#include "GMath.h"

//--------------------------------------------------------------------------------------//
//                                                                                      //
//                                        GVector3                                      //
//                                                                                      //
//--------------------------------------------------------------------------------------//

GVector3::GVector3(float x, float y, float z)
{
	V[0] = x;
	V[1] = y;
	V[2] = z;
}

GVector3::GVector3(const GVector3 & copy)
{
	V[0] = copy.V[0];
	V[1] = copy.V[1];
	V[2] = copy.V[2];
}

GVector3 & GVector3::operator=(const GVector3 & rhs)
{
	V[0] = rhs.V[0];
	V[1] = rhs.V[1];
	V[2] = rhs.V[2];

	return *this;
}

GVector3 &GVector3::operator +=(const GVector3 &rhs)
{
	V[0] += rhs.V[0];
	V[1] += rhs.V[1];
	V[2] += rhs.V[2];
	return *this;
}

GVector3 &GVector3::operator -=(const GVector3 &rhs)
{
	V[0] -= rhs.V[0];
	V[1] -= rhs.V[1];
	V[2] -= rhs.V[2];
	return *this;
}

GVector3 &GVector3::operator *=(const float &s)
{
	V[0] *= s;
	V[1] *= s;
	V[2] *= s;
	return *this;
}

GVector3 &GVector3::operator /=(const float &s)
{
	V[0] /= s;
	V[1] /= s;
	V[2] /= s;
	return *this;
}

GVector3 &GVector3::operator ^=(const GVector3 &rhs)
{
	float x = V[0], y = V[1], z = V[2];
	V[0] = y * rhs.V[2] - z * rhs.V[1];
	V[1] = z * rhs.V[0] - x * rhs.V[2];
	V[2] = x * rhs.V[1] - y * rhs.V[0];
	return *this;
}

bool GVector3::operator ==(const GVector3 &rhs) const
{
	return !((*this) != rhs);
}

bool GVector3::operator !=(const GVector3 &rhs) const
{
	return (!EQ(V[0], rhs.V[0], PRECISION) || !EQ(V[1], rhs.V[1], PRECISION) || !EQ(V[2], rhs.V[2], PRECISION));
}

GVector3 GVector3::operator +() const
{
	return *this;
}

GVector3 GVector3::operator -() const
{
	return *this * -1;
}

GVector3 GVector3::operator+(const GVector3 & rhs) const
{
	return GVector3(V[0] + rhs.V[0], V[1] + rhs.V[1], V[2] + rhs.V[2]);
}

GVector3 GVector3::operator-(const GVector3 & rhs) const
{
	return GVector3(V[0] - rhs.V[0], V[1] - rhs.V[1], V[2] - rhs.V[2]);
}

float norm(const GVector3& v)
{
	return SQRT(v.V[0] * v.V[0] + v.V[1] * v.V[1] + v.V[2] * v.V[2]);
}

GVector3& GVector3::normalize()
{
	float len = norm(*this);
	if (len > PRECISION)
	{
		V[0] /= len;
		V[1] /= len;
		V[2] /= len;
	}
	return *this;
}

GVector3 & GVector3::Set(const float & x, const float & y, const float & z)
{
	V[0] = x;
	V[1] = y;
	V[2] = z;
	return *this;
}

float GVector3::operator*(const GVector3 & rhs) const
{
	return (V[0] * rhs.V[0] + V[1] * rhs.V[1] + V[2] * rhs.V[2]);
}

GVector3 operator*(const GVector3 & lhs, const float & k)
{
	return GVector3(lhs.V[0] * k, lhs.V[1] * k, lhs.V[2] * k);	
}

GVector3 operator*(const float & k, const GVector3 & rhs)
{
	return GVector3(rhs.V[0] * k, rhs.V[1] * k, rhs.V[2] * k);
}

GVector3 operator/(const GVector3 & lhs, const float & k)
{
	return GVector3(lhs.V[0] / k, lhs.V[1] / k, lhs.V[2] / k);
}

GVector3 GVector3::operator ^(const GVector3& rhs) const
{
	return GVector3(V[1] * rhs.V[2] - V[2] * rhs.V[1], V[2] * rhs.V[0] - V[0] * rhs.V[2], V[0] * rhs.V[1] - V[1] * rhs.V[0]);
}

GVector3 proj(const GVector3 & v, const GVector3 & w)
{
		return (v*w / SQR(norm(w)))*w;
}

GVector3 perp(const GVector3 & v, const GVector3 & w)
{
	return v - proj(v, w);
}

ostream & operator<<(ostream & os, const GVector3 & v)
{
	os << "[" << setw(10) << v.V[0] << ", " << setw(10) << v.V[1] << ", " << setw(10) << v.V[2] << "]";
	return os;
}

float &GVector3::operator [](const int &idx)
{
	return V[idx];
}

const float &GVector3::operator [](const int &idx) const
{
	return V[idx];
}

float dist(const GVector3 &v, const GVector3 &w)
{
	return SQRT(SQR(v.V[0] - w.V[0]) + SQR(v.V[1] - w.V[1]) + SQR(v.V[2] - w.V[2]));
}

//--------------------------------------------------------------------------------------//
//                                                                                      //
//                                        GVector                                       //
//                                                                                      //
//--------------------------------------------------------------------------------------//

GVector::GVector(int dim)
{
	this->N = dim;
	V = new float[N];
	ARR_ZERO(V, N);
}

GVector::GVector(int dim, float x, ...)
{
	this->N = dim;
	V = new float[N];
	va_list ap;
	va_start(ap, dim);
	for (int i = 0; i < N; i++)
		V[i] = va_arg(ap, float);
	va_end(ap);
}

GVector::GVector(const GVector3 &copy)
{
	N = 3;
	V = new float[N];
	V[0] = copy[0];
	V[1] = copy[1];
	V[2] = copy[2];
}

GVector::GVector(const GVector &copy)
{
	N = copy.N;
	V = new float[N];
	memcpy(V, copy.V, N * sizeof(float));
}

GVector::~GVector()
{
	if (V)
		delete[] V;
	V = NULL;
}

GVector &GVector::operator =(const GVector &rhs)
{
	if (V)
		delete[] V;
	N = rhs.N;
	V = new float[N];
	memcpy(V, rhs.V, N * sizeof(float));
	return *this;
}

GVector &GVector::operator +=(const GVector &rhs)
{
	assert(N == rhs.N);
	for (int i = 0; i != N; ++i)
		V[i] += rhs.V[i];
	return *this;
}

GVector &GVector::operator -=(const GVector &rhs)
{
	assert(N == rhs.N);
	for (int i = 0; i != N; ++i)
		V[i] -= rhs.V[i];
	return *this;
}

GVector &GVector::operator *=(const float &s)
{
	for (int i = 0; i != N; ++i)
		V[i] *= s;
	return *this;
}

GVector &GVector::operator /=(const float &s)
{
	assert(s != 0);
	for (int i = 0; i != N; ++i)
		V[i] /= s;
	return *this;
}

GVector GVector::operator +() const
{
	return *this;
}

GVector GVector::operator -() const
{
	return *this * -1;
}

GVector GVector::operator +(const GVector &rhs)	const
{
	assert(N == rhs.N);
	GVector ret(N);
	for (int i = 0; i != N; ++i)
		ret.V[i] = V[i] + rhs.V[i];
	return ret;
}

GVector GVector::operator -(const GVector &rhs)	const
{
	assert(N == rhs.N);
	GVector ret(N);
	for (int i = 0; i != N; ++i)
		ret.V[i] = V[i] - rhs.V[i];
	return ret;
}

float GVector::operator *(const GVector &rhs) const
{
	assert(N == rhs.N);
	float ret = 0;
	for (int i = 0; i != N; ++i)
		ret += V[i] * rhs.V[i];
	return ret;
}

GVector GVector::operator /(const float &s) const
{
	GVector ret(N);
	for (int i = 0; i != N; ++i)
		ret.V[i] = V[i] / s;
	return ret;
}

bool GVector::operator ==(const GVector &rhs) const
{
	return !((*this) != rhs);
}

bool GVector::operator !=(const GVector &rhs) const
{
	assert(N == rhs.N);
	for (int i = 0; i != N; ++i)
		if (!EQ(V[i], rhs.V[i], PRECISION))
			return true;
	return false;
}

float &GVector::operator [](const int &idx)
{
	assert(idx >= 0 && idx < N);
	return V[idx];
}

const float &GVector::operator [](const int &idx) const
{
	assert(idx >= 0 && idx < N);
	return V[idx];
}

GVector &GVector::Set(float x, ...)
{
	V[0] = x;
	va_list ap;
	va_start(ap, x);
	for (int i = 1; i != N; ++i)
		V[i] = va_arg(ap, float);
	va_end(ap);
	return *this;
}

GVector &GVector::Set(float *pVal)
{
	memcpy(V, pVal, sizeof(float) * N);
	return *this;
}

GVector &GVector::Normalize()
{
	float len = norm(*this);
	for (int i = 0; i != N; ++i)
		V[i] /= len;
	return *this;
}

GVector &GVector::SetZeros()
{
	for (int i = 0; i != N; ++i)
		V[i] = 0.0;
	return *this;
}

int GVector::GetDim() const
{
	return N;
}


GVector operator *(const GVector &lhs, const float &s)
{
	GVector ret(lhs.N);
	for (int i = 0; i != lhs.N; ++i)
		ret.V[i] = lhs.V[i] * s;
	return ret;
}

GVector operator *(const float &s, const GVector &rhs)
{
	GVector ret(rhs.N);
	for (int i = 0; i != rhs.N; ++i)
		ret.V[i] = rhs.V[i] * s;
	return ret;
}

ostream &operator <<(ostream &os, const GVector &v)
{
	for (int i = 0; i != v.N; ++i)
		if (i != v.N - 1)
			os << "| " << setw(10) << v.V[i] << " |" << endl;
		else
			os << "| " << setw(10) << v.V[i] << " |";
	return os;
}

float norm(const GVector &v)
{
	float len = 0;
	for (int i = 0; i != v.N; ++i)
		len += SQR(v.V[i]);
	len = SQRT(len);
	return len;
}

float dist(const GVector &v, const GVector &w)
{
	return norm(v - w);
}

//--------------------------------------------------------------------------------------//
//                                                                                      //
//                                        GMatrix                                       //
//                                                                                      //
//--------------------------------------------------------------------------------------//

GMatrix::GMatrix(int row, int col, float *elem)
{
	r = row;
	c = col;
	M = new float[r * c];
	if (elem)
		memcpy(M, elem, sizeof(float) * r * c);
	else
		memset(M, 0, sizeof(float) * r * c);
}

GMatrix::GMatrix(const GMatrix &copy)
{
	r = copy.r;
	c = copy.c;
	M = new float[r *c];
	memcpy(M, copy.M, sizeof(float) * r * c);
}

GMatrix::~GMatrix()
{
	if (M)
		delete[] M;
	M = NULL;
}

GMatrix &GMatrix::operator =(const GMatrix &rhs)
{
	if (M)
		delete[] M;

	r = rhs.r;
	c = rhs.c;
	M = new float[r * c];

	memcpy(M, rhs.M, sizeof(float) * r * c);
	return *this;
}

GMatrix &GMatrix::operator +=(const GMatrix &rhs)
{
	assert(r == rhs.r && c == rhs.c);
	for (int i = 0; i < r * c; ++i)
		M[i] += rhs.M[i];
	return *this;
}

GMatrix &GMatrix::operator -=(const GMatrix &rhs)
{
	assert(r == rhs.r && c == rhs.c);
	for (int i = 0; i < r * c; ++i)
		M[i] -= rhs.M[i];
	return *this;
}

GMatrix &GMatrix::operator *=(const GMatrix &rhs)
{
	assert(c == rhs.r);

	c = rhs.c;
	float *newM = new float[r * c];
	memset(newM, 0, sizeof(float) * r * c);

	for (int i = 0; i < r; ++i)
		for (int j = 0; j < c; ++j)
			for (int k = 0; k < rhs.r; ++k)
				newM[i * c + j] += M[i * rhs.r + k] * rhs.M[k * c + j];

	delete[] M;
	M = newM;

	return *this;
}

GMatrix &GMatrix::operator *=(const float &s)
{
	for (int i = 0; i < r * c; i++)
		M[i] *= s;
	return *this;
}

GMatrix &GMatrix::operator /=(const float &s)
{
	for (int i = 0; i < r * c; i++)
		M[i] /= s;
	return *this;
}

GMatrix GMatrix::operator +() const
{
	return *this;
}

GMatrix GMatrix::operator -() const
{
	return *this * -1;
}

GMatrix GMatrix::operator +(const GMatrix &rhs) const
{
	assert(r == rhs.r && c == rhs.c);
	GMatrix ret(*this);
	ret += rhs;
	return ret;
}

GMatrix GMatrix::operator -(const GMatrix &rhs) const
{
	assert(r == rhs.r && c == rhs.c);
	GMatrix ret(*this);
	ret -= rhs;
	return ret;
}

GMatrix GMatrix::operator *(const GMatrix &rhs) const
{
	assert(c == rhs.r);
	GMatrix ret(*this);
	ret *= rhs;
	return ret;
}

GMatrix GMatrix::operator /(const float &s) const
{
	GMatrix ret(*this);
	ret /= s;
	return ret;
}

bool GMatrix::operator ==(const GMatrix &rhs) const
{
	assert(r == rhs.r && c == rhs.c);
	for (int i = 0; i != r * c; ++i)
		if (abs(M[i] - rhs.M[i]) > PRECISION)
			return false;
	return true;
}

bool GMatrix::operator !=(const GMatrix &rhs) const
{
	assert(r == rhs.r && c == rhs.c);
	for (int i = 0; i != r * c; ++i)
		if (abs(M[i] - rhs.M[i]) > PRECISION)
			return true;
	return false;
}

float *GMatrix::operator [](const int idx)
{
	assert(idx >= 0 && idx < r);
	return &M[idx * c];
}

const float *GMatrix::operator [](const int idx) const
{
	assert(idx >= 0 && idx < r);
	return &M[idx * c];
}

GMatrix &GMatrix::SetTranspose()
{
	int i, j;
	if (r == c)	// Square matrix
	{
		for (i = 0; i != r; ++i)
			for (j = i + 1; j != c; ++j)
				SWAP(float, M[i * c + j], M[j * c + i]);
	}
	else		// Rectangular matrix
	{
		float *buffer = new float[r * c];
		memcpy(buffer, M, sizeof(float) * r * c);
		SWAP(int, r, c);
		for (i = 0; i != r; ++i)
			for (j = 0; j != c; ++j)
				M[i * c + j] = buffer[j * r + i];
		delete[] buffer;
	}
	return *this;
}

GMatrix &GMatrix::SetIdentity()
{
	memset(M, 0, sizeof(float) * r * c);
	int min = MIN(r, c);
	for (int i = 0; i < min; i++)
		M[i * c + i] = 1.0;
	return *this;
}

GMatrix &GMatrix::SetZeros()
{
	memset(M, 0, sizeof(float) * r * c);
	return *this;
}

GMatrix &GMatrix::SetRowVec(const int idx, const GVector &v)
{
	assert(idx < r);
	assert(v.N == c);
	for (int i = 0; i < c; i++)
		M[idx * c + i] = v.V[i];
	return *this;
}

GMatrix &GMatrix::SetColVec(const int idx, const GVector &v)
{
	assert(idx < c);
	assert(v.N == r);
	for (int i = 0; i < r; i++)
		M[i * c + idx] = v.V[i];
	return *this;
}

GMatrix &GMatrix::ExchangeRows(const int idx0, const int idx1)
{
	GVector tmp(c);
	tmp = GetRowVec(idx0);
	SetRowVec(idx0, GetRowVec(idx1));
	SetRowVec(idx1, tmp);
	return *this;
}

GMatrix &GMatrix::ExchangeCols(const int idx0, const int idx1)
{
	GVector tmp(r);
	tmp = GetColVec(idx0);
	SetColVec(idx0, GetColVec(idx1));
	SetColVec(idx1, tmp);
	return *this;
}

int GMatrix::GetRowNum() const
{
	return r;
}

int GMatrix::GetColNum() const
{
	return c;
}

GVector GMatrix::GetRowVec(const int idx) const
{
	assert(idx < r);
	GVector ret(c);
	for (int i = 0; i < c; i++)
		ret.V[i] = M[idx * c + i];
	return ret;
}

GVector GMatrix::GetColVec(const int idx) const
{
	assert(idx < c);
	GVector ret(r);
	for (int i = 0; i < r; i++)
		ret.V[i] = M[i * c + idx];
	return ret;
}

bool GMatrix::IsSquare() const
{
	return (r == c) ? true : false;
}

GVector operator *(const GMatrix &lhs, const GVector &rhs)
{
	assert(lhs.c == rhs.N);
	GVector ret(lhs.r);
	for (int i = 0; i != lhs.r; ++i)		// for each row.
		for (int j = 0; j != lhs.c; ++j)	// for each col.
			ret.V[i] += lhs.M[i * lhs.c + j] * rhs.V[j];
	return ret;
}

GMatrix operator *(const GVector &lhs, const GMatrix &rhs)
{
	assert(rhs.r == 1);
	GMatrix ret(lhs.N, rhs.c);
	for (int i = 0; i != lhs.N; ++i)		// for each row.
		for (int j = 0; j != rhs.c; ++j)	// for each col.
			ret.M[i * rhs.c + j] = lhs.V[i] * rhs.M[j];
	return ret;
}

GMatrix operator *(const GMatrix &lhs, const float &s)
{
	GMatrix ret(lhs);
	ret *= s;
	return ret;
}

GMatrix operator *(const float &s, const GMatrix &rhs)
{
	GMatrix ret(rhs);
	ret *= s;
	return ret;
}

ostream &operator <<(ostream &os, const GMatrix &m)
{
	for (int i = 0; i != m.r; i++)
	{
		os << "| ";
		for (int j = 0; j != m.c; j++)
			os << setw(10) << m.M[i * m.c + j] << " ";
		os << "|" << endl;
	}
	return os;
}

GMatrix RowEchelonForm(const GMatrix &m)
{
	int i, j, k; 
	int r = m.GetRowNum();
	int c = m.GetColNum();
	int n = MIN(r, c);
	GMatrix T(m);

	int shift = 0;
	for (i = 0; i < n; i++)
	{
		float maxi = ABS(T[i][i + shift]); // 取第1行第1列
		int pivot_idx = i; //保存当前行号
		for (j = i + 1; j < n; j++) // 找到当前列中的最大值
		{
			if (maxi < ABS(T[j][i + shift]))
			{
				maxi = ABS(T[j][i + shift]);
				pivot_idx = j; //更新行号为最大值所在行号
			}
		}

		if (EQ_ZERO(maxi, PRECISION)) // 判断最大值是否为0
		{
			shift++;
			i--;
			continue;
		}

		// 如果当前行不是包含最大值的行，则交换两行
		if (i != pivot_idx)
			T.ExchangeRows(i, pivot_idx);

		// 取出最大值
		float s = T[i][i + shift];
		for (j = i + shift; j < c; j++) // 找出leading 1（除以最大值）
			T[i][j] = T[i][j] / s;

		// 将leading 1下面的元素变为0（一行减去s倍leading 1所在行）
		for (j = i + 1; j < r; j++)
		{
			
			s = T[j][i + shift];
			for (k = i + shift; k < c; k++)
			{
				T[j][k] = T[j][k] - s * T[i][k];
			}
		}
	}

	return T;
}

GMatrix ReducedRowEchelon(const GMatrix &m)
{
	int i, j, k; 
	int r = m.GetRowNum();
	int c = m.GetColNum();
	int n = MIN(r, c);
	GMatrix T(m);

	int shift = 0;
	for (i = 0; i < n; i++)
	{
		// pivoting.
		float maxi = ABS(T[i][i + shift]);
		int pivot_idx = i;
		for (j = i + 1; j < n; j++)
		{
			if (maxi < ABS(T[j][i + shift]))
			{
				maxi = ABS(T[j][i + shift]);
				pivot_idx = j;
			}
		}

		if (EQ_ZERO(maxi, PRECISION))
		{
			shift++;
			i--;
			continue;
		}

		if (i != pivot_idx)
			T.ExchangeRows(i, pivot_idx);

		float s = T[i][i + shift];
		for (j = i + shift; j < c; j++)
			T[i][j] = T[i][j] / s;

		// 将 leading 1 上下的元素变为0
		for (j = 0; j < r; j++)
		{
			if (i == j)
				continue;

			s = T[j][i + shift];
			for (k = i + shift; k < c; k++)
			{
				T[j][k] = T[j][k] - s * T[i][k];
			}
		}
	}

	return T;
}

float *form_arr(const GMatrix &m)
{
	return m.M;
}

int Rank(const GMatrix &m)
{
	int i, r, rank = 0;
	r = m.GetRowNum();

	GMatrix T = ref(m);
	for (i = 0; i < r; i++)
	{
		GVector rVec = T.GetRowVec(i);
		if (!EQ_ZERO(norm(rVec), PRECISION))
			rank++;
	}

	return rank;
}

int Nullity(const GMatrix &m)
{
	int rnk = Rank(m);
	int c = m.GetColNum();
	return (c - rnk);
}

//--------------------------------------------------------------------------------------//
//                                                                                      //
//                                        GLine                                         //
//                                                                                      //
//--------------------------------------------------------------------------------------//

GLine::GLine(const GPoint3 &_p, const GVector3 &_v)
{
	p = _p;
	v = _v;
}


GLine::GLine(const GLine &copy) : p(copy.p), v(copy.v)
{
}

GPoint3 GLine::operator ()(const float t) const
{
	return p + t * v;
}

GLine &GLine::operator =(const GLine &rhs)
{
	this->p = rhs.p;
	this->v = rhs.v;
	return *this;
}

bool GLine::IsOnLine(const GPoint3 &q) const
{
	return EQ_ZERO( dist(q, *this),PRECISION);
}

ostream &operator <<(ostream &os, const GLine &l)
{
	os << "("
		<< l.p[0] << " + (" << l.v[0] << ") * t, "
		<< l.p[1] << " + (" << l.v[1] << ") * t, "
		<< l.p[2] << " + (" << l.v[2] << ") * t)";
	return os;
}

float dist(const GPoint3 &q, const GLine &l)
{
	return norm(proj(q - l.p, l.v) - (q - l.p));
}

GLine &GLine::SetPt(const GPoint3 &_p)
{
	p = _p;
	return *this;
}

GLine &GLine::SetDir(const GVector3 &_v)
{
	v = _v;
	return *this;
}

GPoint3 GLine::GetPt() const
{
	return p;
}

GVector3 GLine::GetDir() const
{
	return v;
}

//--------------------------------------------------------------------------------------//
//                                                                                      //
//                                        GPlane                                        //
//                                                                                      //
//--------------------------------------------------------------------------------------//

GPlane::GPlane(const GVector3 &_n, const GPoint3 &_p)
{
	n = _n;
	d = -n*_p;// -(n[0] * _p[0] + n[1] * _p[1] + n[2] * _p[2]);
}

GPlane::GPlane(const GPoint3 &p1, const GPoint3 &p2, const GPoint3 &p3)
{
	n = (p2 - p1) ^ (p3 - p1);
	d = -n*p1;
}

GPlane::GPlane(const float &a, const float &b, const float &c, const float &d)
{
	// ax + by + cz + d = 0
	// n = (a,b,x)
	this->n = GVector3(a, b, c);
	this->d = d;
}

GPlane::GPlane(const GPlane &copy)
{
	this->n = copy.n;
	this->d = copy.d;
}

GPlane &GPlane::operator =(const GPlane &rhs)
{
	this->n = rhs.n;
	this->d = rhs.d;
	return *this;
}

GVector3 GPlane::GetNormal() const
{
	return n;
}

bool GPlane::IsOnPlane(const GPoint3 &p) const
{
	float s;
	s = -n*p;
	if (EQ(s, d, PRECISION))
		return true;
	else
		return false;
}

bool GPlane::IsAbovePlane(const GPoint3 &p) const
{
	float s;
	s = -n*p;
	if (s > 0.0f)
		return true;
	else
		return false;
}

bool GPlane::IsBelowPlane(const GPoint3 &p) const
{
	float s;
	s = -n*p;
	if (s < 0.0f)
		return true;
	else
		return false;
}

ostream &operator <<(ostream &os, const GPlane &pi)
{
	os << "(" << pi.n[0] << ") * x + ("
		<< pi.n[1] << ") * y + ("
		<< pi.n[2] << ") * z + ("
		<< pi.d << ") = 0";
	return os;
}

float dist(const GPlane &pi, const GPoint3 &p)
{
	float D;
	D = (p[0] * pi.n[0] + p[1] * pi.n[1] + p[2] * pi.n[2] + pi.d) / norm(pi.n);
	return D;
}

/*!
*	\brief	求直线l和平面pi的交点
*
*	\param p 所求交点（如果相交）
*	\param l 直线
*	\param pi 平面
*
*	\return true: 相交, false: 直线和平面平行
*/
bool intersect_line_plane(GPoint3 &p, const GLine &l, const GPlane &pi)
{
	if (EQ_ZERO(l.v * pi.n, PRECISION))
	{
		return false;
	}

	float t = -(l.p[0] * pi.n[0] + l.p[1] * pi.n[1] + l.p[2] * pi.n[2] + pi.d) / (l.v * pi.n);
	p = l(t);
	return true;
}

/*!
*	\brief	求直线l和三角形p1p2p3的交点
*
*	\param q 所求交点（如果相交）
*	\param l 直线
*	\param p1 三角形顶点
*	\param p2 三角形顶点
*	\param p3 三角形顶点
*
*	\return true: 相交, false: 直线和三角形平行或交点在三角形外
*/
bool intersect_line_triangle(GPoint3 &q, const GLine &l, const GPoint3 &p1, const GPoint3 &p2, const GPoint3 &p3)
{
	GVector3 e1, e2, u, v, w;
	float det, alpha, beta, t;
	e1 = p2 - p1;
	e2 = p3 - p1; 
	u = l.v ^ e2;
	det = e1 * u;

	
	if (EQ_ZERO(det,PRECISION))
		return false;

	w = l.p - p1;
	alpha = w * u / det;
	if (alpha < 0.0 || alpha > 1.0)
		return false;

	v = w ^ e1;
	beta = l.v * v / det;
	if (beta < 0.0 || alpha + beta > 1.0)
		return false;

	t = e2 * v / det;
	

	q = l(t);
	return true;
}

//--------------------------------------------------------------------------------------//
//                                                                                      //
//                                        GQuater                                       //
//                                                                                      //
//--------------------------------------------------------------------------------------//

GQuater::GQuater(float w, float x, float y, float z)
{
	this->W = w;
	this->X = x;
	this->Y = y;
	this->Z = z;
}

GQuater::GQuater(const GQuater &copy)
{
	this->W = copy.W;
	this->X = copy.X;
	this->Y = copy.Y;
	this->Z = copy.Z;
}

GQuater::GQuater(const float *q, const bool invOrder)
{
	if (invOrder)
	{
		this->W = q[1];
		this->X = q[2];
		this->Y = q[3];
		this->Z = q[0];
	}
	else
	{
		this->W = q[0];
		this->X = q[1];
		this->Y = q[2];
		this->Z = q[3];
	}
}

GQuater::GQuater(GVector3 axis, float theta, bool radian)
{
	float rad, sn, cs;
	axis.normalize();
	if (!radian)
		rad = theta * (float)M_PI / 360.0f;

	sn = sin(rad);
	cs = cos(rad);

	sn = (abs(sn) < PRECISION) ? 0.0f : sn;
	cs = (abs(cs) < PRECISION) ? 0.0f : cs;

	W = cs;
	X = sn * axis[0];
	Y = sn * axis[1];
	Z = sn * axis[2];
}

GQuater::GQuater(float theta1, float theta2, float theta3, EulerType eulerType)
{
	float c1, c2, c3;
	float s1, s2, s3;
	theta1 = DEG2RAD(theta1);
	theta2 = DEG2RAD(theta2);
	theta3 = DEG2RAD(theta3);
	c1 = cos(theta1);
	c2 = cos(theta2);
	c3 = cos(theta3);
	s1 = sin(theta1);
	s2 = sin(theta2);
	s3 = sin(theta3);

	GMatrix mat;
	switch (eulerType)
	{
	case EULER_XYZ:
		mat[0][0] = c2 * c3;
		mat[0][1] = -c2 * s3;
		mat[0][2] = s2;
		mat[0][3] = 0.0;

		mat[1][0] = s1 * s2 * c3 + c1 * s3;
		mat[1][1] = -s1 * s2 * s3 + c1 * c3;
		mat[1][2] = -s1 * c2;
		mat[1][3] = 0.0;

		mat[2][0] = -c1 * s2 * c3 + s1 * s3;
		mat[2][1] = c1 * s2 * s3 + s1 * c3;
		mat[2][2] = c1 * c2;
		mat[2][3] = 0.0;

		mat[3][0] = mat[3][1] = mat[3][2] = 0.0;
		mat[3][3] = 1.0;
		break;
	case EULER_ZYX:
		mat[0][0] = c3 * c2;
		mat[0][1] = -s3 * c1 + c3 * s2 * s1;
		mat[0][2] = s3 * s1 + c3 * s2 * c1;
		mat[0][3] = 0.0;

		mat[1][0] = s3 * c2;
		mat[1][1] = c3 * c1 + s3 * s2 * s1;
		mat[1][2] = -c3 * s1 + s3 * s2 * c1;
		mat[1][3] = 0.0;

		mat[2][0] = -s2;
		mat[2][1] = c2 * s1;
		mat[2][2] = c2 * c1;
		mat[2][3] = 0.0;

		mat[3][0] = mat[3][1] = mat[3][2] = 0.0;
		mat[3][3] = 1.0;
		break;
	}
	SetFromMatrix(mat);
}

GQuater &GQuater::operator =(const GQuater &rhs)
{
	W = rhs.W;
	X = rhs.X;
	Y = rhs.Y;
	Z = rhs.Z;
	return *this;
}

GQuater &GQuater::operator +=(const GQuater &rhs)
{
	W += rhs.W;
	X += rhs.X;
	Y += rhs.Y;
	Z += rhs.Z;
	return *this;
}

GQuater &GQuater::operator -=(const GQuater &rhs)
{
	W -= rhs.W;
	X -= rhs.X;
	Y -= rhs.Y;
	Z -= rhs.Z;
	return *this;
}

GQuater &GQuater::operator *=(const GQuater &rhs)
{
	float w = W, x = X, y = Y, z = Z;
	this->W = w * rhs.W - x * rhs.X - y * rhs.Y - z * rhs.Z;
	this->X = w * rhs.X + rhs.W * x + y * rhs.Z - z * rhs.Y;
	this->Y = w * rhs.Y + rhs.W * y + z * rhs.X - x * rhs.Z;
	this->Z = w * rhs.Z + rhs.W * z + x * rhs.Y - y * rhs.X;
	return *this;
}

GQuater &GQuater::operator /=(const GQuater &rhs)
{
	(*this) = inv(rhs) * (*this);
	Normalize();
	return *this;
}

GQuater &GQuater::operator *=(const float s)
{
	W *= s;
	X *= s;
	Y *= s;
	Z *= s;
	return *this;
}

GQuater &GQuater::operator /=(const float s)
{
	W /= s;
	X /= s;
	Y /= s;
	Z /= s;
	return *this;
}

GQuater GQuater::operator -() const
{
	return *this * -1;
}

GQuater GQuater::operator +() const
{
	return *this;
}

GQuater GQuater::operator +(const GQuater &rhs) const
{
	GQuater ret(*this);
	ret += rhs;
	return ret;
}

GQuater GQuater::operator -(const GQuater &rhs) const
{
	GQuater ret(*this);
	ret -= rhs;
	return ret;
}

GQuater GQuater::operator *(const GQuater &rhs) const
{
	GQuater ret(*this);
	ret *= rhs;
	return ret;
}

GQuater GQuater::operator /(const GQuater &rhs) const
{
	GQuater ret(*this);
	ret /= rhs;
	return ret;
}

GQuater GQuater::operator /(const float s) const
{
	GQuater ret(*this);
	ret /= s;
	return ret;
}

GVector3 GQuater::operator *(const GVector3 &rhs) const
{
	assert(IsUnitQuater());
	GVector3 ret;
	GQuater v(0.0, rhs[0], rhs[1], rhs[2]);
	GQuater rq = (*this) * v * inv(*this);
	ret.Set(rq.X, rq.Y, rq.Z);
	return ret;
}

bool GQuater::operator ==(const GQuater &rhs) const
{
	return !((*this) != rhs);
}

bool GQuater::operator !=(const GQuater &rhs) const
{
	return (!EQ(W, rhs.W, PRECISION) || !EQ(X, rhs.X, PRECISION) ||
		!EQ(Y, rhs.Y, PRECISION) || !EQ(Z, rhs.Z, PRECISION));
}

GQuater &GQuater::Set(const float w, const float x, const float y, const float z)
{
	this->W = w;
	this->X = x;
	this->Y = y;
	this->Z = z;
	return *this;
}

GQuater &GQuater::Set(float *q, bool invOrder)
{
	if (invOrder)
	{
		this->W = q[1];
		this->X = q[2];
		this->Y = q[3];
		this->Z = q[0];
	}
	else
	{
		this->W = q[0];
		this->X = q[1];
		this->Y = q[2];
		this->Z = q[3];
	}
	return *this;
}

GQuater &GQuater::SetFromAngleAxis(const float theta, GVector3 axis, bool radian)
{
	float rad, sn, cs;
	axis.normalize();
	if (!radian)
		rad = theta * M_PI / 360.0f;
	else
		rad = theta / 2.0f;

	sn = sin(rad);
	cs = cos(rad);

	sn = (abs(sn) < PRECISION) ? 0.0f : sn;
	cs = (abs(cs) < PRECISION) ? 0.0f : cs;

	W = cs;
	X = sn * axis[0];
	Y = sn * axis[1];
	Z = sn * axis[2];

	return *this;
}

GQuater &GQuater::SetFromEulerAngle(float theta1, float theta2, float theta3, EulerType eulerType)
{
	float c1, c2, c3;
	float s1, s2, s3;
	theta1 = DEG2RAD(theta1);
	theta2 = DEG2RAD(theta2);
	theta3 = DEG2RAD(theta3);
	c1 = cos(theta1);
	c2 = cos(theta2);
	c3 = cos(theta3);
	s1 = sin(theta1);
	s2 = sin(theta2);
	s3 = sin(theta3);

	GMatrix mat;
	switch (eulerType)
	{
	case EULER_XYZ:
		mat[0][0] = c2 * c3;
		mat[0][1] = -c2 * s3;
		mat[0][2] = s2;
		mat[0][3] = 0.0;

		mat[1][0] = s1 * s2 * c3 + c1 * s3;
		mat[1][1] = -s1 * s2 * s3 + c1 * c3;
		mat[1][2] = -s1 * c2;
		mat[1][3] = 0.0;

		mat[2][0] = -c1 * s2 * c3 + s1 * s3;
		mat[2][1] = c1 * s2 * s3 + s1 * c3;
		mat[2][2] = c1 * c2;
		mat[2][3] = 0.0;

		mat[3][0] = mat[3][1] = mat[3][2] = 0.0;
		mat[3][3] = 1.0;
		break;
	case EULER_ZYX:
		mat[0][0] = c3 * c2;
		mat[0][1] = -s3 * c1 + c3 * s2 * s1;
		mat[0][2] = s3 * s1 + c3 * s2 * c1;
		mat[0][3] = 0.0;

		mat[1][0] = s3 * c2;
		mat[1][1] = c3 * c1 + s3 * s2 * s1;
		mat[1][2] = -c3 * s1 + s3 * s2 * c1;
		mat[1][3] = 0.0;

		mat[2][0] = -s2;
		mat[2][1] = c2 * s1;
		mat[2][2] = c2 * c1;
		mat[2][3] = 0.0;

		mat[3][0] = mat[3][1] = mat[3][2] = 0.0;
		mat[3][3] = 1.0;
		break;
	}
	return SetFromMatrix(mat);
}

GQuater &GQuater::SetFromMatrix(const GMatrix &m)
{
	int r, c;
	r = m.GetRowNum();
	c = m.GetColNum();
	assert(r == 3 || r == 4);
	assert(c == 3 || c == 4);

	float q[4];
	float tr, s;
	int i, j, k;
	int nxt[3] = { 1, 2, 0 };

	tr = m[0][0] + m[1][1] + m[2][2];

	if (tr > 0.0f)
	{
		s = SQRT(tr + 1.0f);
		W = s * 0.5f;
		s = 0.5f / s;
		X = (m[2][1] - m[1][2]) * s;
		Y = (m[0][2] - m[2][0]) * s;
		Z = (m[1][0] - m[0][1]) * s;
	}
	else
	{
		i = 0;
		if (m[1][1] > m[0][0])
			i = 1;
		if (m[2][2] > m[i][i])
			i = 2;
		j = nxt[i];
		k = nxt[j];
		s = SQRT((m[i][i] - (m[j][j] + m[k][k])) + 1.0f);
		q[i] = s * 0.5f;
		s = 0.5f / s;
		W = (m[k][j] - m[j][k]) * s;
		q[j] = (m[j][i] + m[i][j]) * s;
		q[k] = (m[k][i] + m[i][k]) * s;
		X = q[0];
		Y = q[1];
		Z = q[2];
	}
	Normalize();
	return *this;
}

GQuater &GQuater::SetIdentity()
{
	this->W = 1.0;
	this->X = 0.0;
	this->Y = 0.0;
	this->Z = 0.0;
	return *this;
}

GQuater& GQuater::SetConjugate()
{
	this->X *= -1.0;
	this->Y *= -1.0;
	this->Z *= -1.0;
	return *this;
}

GQuater &GQuater::SetInverse()
{
	if (!IsUnitQuater())
	{
		float norm_sqr = SQR(W) + SQR(X) + SQR(Y) + SQR(Z);
		*this /= norm_sqr;
	}
	SetConjugate();
	return *this;
}

GQuater &GQuater::Normalize()
{
	float len = norm(*this);
	this->W /= len;
	this->X /= len;
	this->Y /= len;
	this->Z /= len;
	return *this;
}

bool GQuater::IsUnitQuater() const
{
	float norm = SQR(W) + SQR(X) + SQR(Y) + SQR(Z);
	return EQ(norm, 1.0, PRECISION) ? true : false;
}

bool GQuater::IsIdentity() const
{
	return (EQ(W, 1.0, 1.0e-5) && EQ(X, 0.0, 1.0e-5) &&
		EQ(Y, 0.0, 1.0e-5) && EQ(Z, 0.0, 1.0e-5));
}

GQuater operator *(const GQuater &lhs, const float &s)
{
	GQuater ret(lhs);
	ret *= s;
	return ret;
}

GQuater operator *(const float &s, const GQuater &rhs)
{
	GQuater ret(rhs);
	ret *= s;
	return ret;
}

float norm(const GQuater &q)
{
	return SQRT(SQR(q.W) + SQR(q.X) + SQR(q.Y) + SQR(q.Z));
}

GQuater inv(const GQuater &q)
{
	GQuater ret(q);
	if (!ret.IsUnitQuater())
	{
		float norm_sqr = SQR(ret.W) + SQR(ret.X) + SQR(ret.Y) + SQR(ret.Z);
		ret /= norm_sqr;
	}
	ret.Set(ret.W, -ret.X, -ret.Y, -ret.Z);
	return ret;
}

ostream &operator <<(ostream &os, const GQuater &q)
{
	os << "(" << q.W << ") + " << "(" << q.X << ")i + " << "(" << q.Y << ")j + " << "(" << q.Z << ")k" << endl;
	return os;
}