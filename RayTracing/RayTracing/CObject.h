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

	// ���ݲ����͵������ȡ������һ��ķ���
	virtual GVector3 getN(GVector3 _p = GVector3(0.0, 0.0, 0.0)) = 0;

	// �ж������Ƿ�������н���
	virtual bool isIntersect(CRay _ray) = 0;
	virtual GVector3 getIntersection(CRay _ray) = 0;

protected:

	GVector3 Ka;	// ��������ϵ��
	GVector3 Kd;	// ������ϵ��
	GVector3 Ks;	// ���淴��ϵ��
	float Shininess;	// �⻬�̶�
	float Reflectivity;	// ��������ǿ��
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

	GVector3 getN(GVector3 _p = GVector3(0.0, 0.0, 0.0)); // ��ȡƽ����_P��ķ���
	bool isIntersect(CRay _ray);
	GVector3 getIntersection(CRay _ray);



private:
	GVector3 N; // ����
	GPoint3  P; // ƽ����һ��
	float D; // ���� Ax + By + cZ + D = 0 �еĳ���D

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

	// ��������λ��
	inline void setC(const GVector3& _C) { C = _C; }

	// ������뾶
	inline void setR(const float& _R) { R = _R; }

	// ��ȡ����
	inline GVector3 getC() { return C; }

	// ��ȡ��뾶
	inline float getR() { return R; }
	
	GVector3 getN(GVector3 _p = GVector3(0.0, 0.0, 0.0)); // ��ȡ������_P��ķ���
	bool isIntersect(CRay _ray);
	GVector3 getIntersection(CRay _ray);

private:
	GVector3 C; // ����
	float R;	// �뾶
};