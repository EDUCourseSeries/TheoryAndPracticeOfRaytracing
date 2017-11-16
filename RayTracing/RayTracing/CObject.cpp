#include "CObject.h"


// ---------------------------------- CObject ----------------------------------

CObject::CObject(const GVector3 & _Ka, const GVector3 & _Kd, const GVector3 & _Ks, const float & _shininess, const float & _reflectivity)
{
	Ka = _Ka;
	Kd = _Kd;
	Ks = _Ks;
	Shininess = _shininess;
	Reflectivity = _reflectivity;
}

CObject::~CObject()
{
}


// ---------------------------------- CPlane ----------------------------------

CPlane::CPlane(const GVector3 & _Ka, const GVector3 & _Kd, const GVector3 & _Ks, const float & _shininess, const float & _reflectivity, const GVector3 & _N, const GPoint3 & _P)
{
	Ka = _Ka;
	Kd = _Kd;
	Ks = _Ks;
	Shininess = _shininess;
	Reflectivity = _reflectivity;
	N = _N;
	P = _P;
	D = -N*P;
}

GVector3 CPlane::getN(GVector3 _p)
{
	return N;
}

bool CPlane::isIntersect(CRay _ray)
{
	double v = N*(_ray.getV());

	return EQ_ZERO(v, 0.000001f) ? false : true;
}

GVector3 CPlane::getIntersection(CRay _ray)
{
	GVector3 p = _ray.getP();
	GVector3 v = (_ray.getV()).normalize();

	float t = -(p[0] * N[0] + p[1] * N[1] + p[2] * N[2] + D) / (v * N);

	return _ray.getPoint(t);
}

// ---------------------------------- CSphere ----------------------------------

CSphere::CSphere(const GVector3 & _Ka, const GVector3 & _Kd, const GVector3 & _Ks, const float & _shininess, const float & _reflectivity, const GVector3 & _C, const float & _R)
{
	Ka = _Ka;
	Kd = _Kd;
	Ks = _Ks;
	Shininess = _shininess;
	Reflectivity = _reflectivity;
	C = _C;
	R = _R;
}

GVector3 CSphere::getN(GVector3 _p)
{
	return (_p - C).normalize();
}

bool CSphere::isIntersect(CRay _ray)
{
	GVector3 u = C - _ray.getP();
	GVector3 d = perp(u, _ray.getV());

	if (norm(d) - R <= .0f)
		return true;
	else
		return false;
}

GVector3 CSphere::getIntersection(CRay _ray)
{
	GVector3 v = _ray.getP() - C;

	float B = 2.0f * (_ray.getV() * v);
	float C = v*v - R*R;

	float discriminant = (B*B) - 4.0f*C;

	discriminant = sqrt(discriminant);

	float s0 = (-B + discriminant) / 2.0f;
	float s1 = (-B - discriminant) / 2.0f;

	return _ray.getPoint(s0<s1 ? s0 : s1);
}
