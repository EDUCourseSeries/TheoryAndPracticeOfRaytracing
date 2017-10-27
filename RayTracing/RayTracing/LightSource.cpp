#include "LightSource.h"



CLightSource::CLightSource()
{
}


CLightSource::~CLightSource()
{
}


CDirectionalLight::CDirectionalLight(const GVector3 & _position, const GVector3 & _direction, const GVector3 & _Ka, const GVector3 & _Kd, const GVector3 & _Ks)
{
	Position = _position;
	Direction = _direction;
	Ka = _Ka;
	Kd = _Kd;
	Ks = _Ks;
}

CDirectionalLight::~CDirectionalLight()
{
}

GVector3 CDirectionalLight::EvalAmbient(const GVector3 & _Ka)
{
	return GVector3(Ka[0] * _Ka[0], Ka[1] * _Ka[1], Ka[2] * _Ka[2]);
	
}

GVector3 CDirectionalLight::EvalDiffuse(const GVector3 & _n, const GVector3 & _l, const GVector3 & _Kd)
{
	GVector3 IdKd = GVector3(Kd[0] * _Kd[0], Kd[1] * _Kd[1], Kd[2] * _Kd[2]);

	float NdotL = MAX(_n*_l, 0.0f);

	return IdKd*NdotL;
	
}

GVector3 CDirectionalLight::EvalSpecluar(const GVector3 & _n, const GVector3 & _l, const GVector3 & _v, const GVector3 & _Ks, const float & _shininess)
{
	GVector3 IsKs = GVector3(Ks[0] * _Ks[0], Ks[1] * _Ks[1], Ks[2] * _Ks[2]);

	GVector3 H = (_l + _v).normalize();

	float NdotL = MAX(_n*_l, 0.0f);
	float NdotH = pow(MAX(_n*H, 0.0f), _shininess);

	if (NdotL <= 0.0)
		NdotH = 0.0;

	return IsKs*NdotH;
}
