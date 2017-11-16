#pragma once
#include "GMath.h"

class CRay
{
public:
	CRay(const GPoint3& _P = GPoint3(), const GVector3& _V = GVector3());

	inline void setP(const GPoint3& _P) { P = _P; }
	inline void setV(const GVector3& _V) { V = _V; }
	
	// ��ȡ�������
	inline GPoint3  getP() { return P; }
	
	// ��ȡ���߷���
	inline GVector3 getV() { return V; }

	// ��ȡ����tʱ�ĵ�
	inline GPoint3 getPoint(const float& t) { return P + V.normalize()*t; }


private:

	GPoint3 P;	// ���
	GVector3 V; // ����
};