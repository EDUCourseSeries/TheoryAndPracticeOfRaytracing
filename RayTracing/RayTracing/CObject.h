#pragma once
#include "CRay.h"

class CObject
{
public:

	CObject(const GVector3& _Ka = GVector3(1.0f, 0.0f, 0.0f),
		    const GVector3& _Kd = GVector3(1.0f, 0.0f, 0.0f),
		    const GVector3& _Ks = GVector3(1.0f, 0.0f, 0.0f),
		    const float& _shininess = 64.0f,
			const float& _reflectivity=0.0f);

	virtual ~CObject();

	// Getters and Setters
	inline void setKa(const GVector3& _Ka) { Ka = _Ka; }
	inline void setKd(const GVector3& _Kd) { Kd = _Kd; }
	inline void setKs(const GVector3& _Ks) { Ks = _Ks; }
	inline void setShininess(const float& _shininess) { Shininess = _shininess; }
	inline void setReflectivity(const float& _reflectivity) { Reflectivity = _reflectivity; }

	inline GVector3 getKa() { return Ka; }
	inline GVector3 getKd() { return Kd; }
	inline GVector3 getKs() { return Ks; }
	inline float getShininess() { return Shininess; }
	inline float getReflectivity() { return Reflectivity; }

	// 根据不类型的物体获取物体上一点的法线
	virtual GVector3 getN(GVector3 _p = GVector3(0.0, 0.0, 0.0)) = 0;

	// 判断射线是否和物体有交点
	virtual bool isIntersect(CRay _ray) = 0;
	virtual GVector3 getIntersection(CRay _ray) = 0;

protected:

	GVector3 Ka;	// 环境反射系数
	GVector3 Kd;	// 漫反射系数
	GVector3 Ks;	// 镜面反射系数
	float Shininess;	// 光滑程度
	float Reflectivity;	// 环境反射强度
};

////////////////////////////////////////////////////////////////////////////////////////////////////////

class CPlane : public CObject
{
public:
	CPlane(const GVector3& _Ka = GVector3(1.0f, 0.0f, 0.0f),
		const GVector3& _Kd = GVector3(1.0f, 0.0f, 0.0f),
		const GVector3& _Ks = GVector3(1.0f, 0.0f, 0.0f),
		const float& _shininess = 64.0f,
		const float& _reflectivity = 0.0f,
		const GVector3& _N = GVector3(),
		const GPoint3& _P = GPoint3());

	inline void setN(GVector3 _N) { N = _N.normalize(); }
	inline void setP(const GPoint3& _P) { P = _P; }

	GVector3 getN(GVector3 _p = GVector3(0.0, 0.0, 0.0)); // 获取平面上_P点的法线
	bool isIntersect(CRay _ray);
	GVector3 getIntersection(CRay _ray);



private:
	GVector3 N; // 法线
	GPoint3  P; // 平面上一点
	float D; // 方程 Ax + By + cZ + D = 0 中的常数D

};

////////////////////////////////////////////////////////////////////////////////////////////////////////

class CSphere : public CObject
{
public:
	CSphere(const GVector3& _Ka = GVector3(1.0f, 0.0f, 0.0f),
		const GVector3& _Kd = GVector3(1.0f, 0.0f, 0.0f),
		const GVector3& _Ks = GVector3(1.0f, 0.0f, 0.0f),
		const float& _shininess = 64.0f,
		const float& _reflectivity = 0.0f,
		const GVector3& _C = GVector3(),
		const float& _R = 1.0f);

	// 设置球心位置
	inline void setC(const GVector3& _C) { C = _C; }

	// 设置球半径
	inline void setR(const float& _R) { R = _R; }

	// 获取球心
	inline GVector3 getC() { return C; }

	// 获取球半径
	inline float getR() { return R; }
	
	GVector3 getN(GVector3 _p = GVector3(0.0, 0.0, 0.0)); // 获取球面上_P点的法线
	bool isIntersect(CRay _ray);
	GVector3 getIntersection(CRay _ray);

private:
	GVector3 C; // 球心
	float R;	// 半径
};