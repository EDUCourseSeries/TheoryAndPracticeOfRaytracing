#pragma once
#include "GMath.h"

class CRay
{
public:
	CRay(const GPoint3& _P = GPoint3(), const GVector3& _V = GVector3());

	inline void setP(const GPoint3& _P) { P = _P; }
	inline void setV(const GVector3& _V) { V = _V; }
	
	// 获取射线起点
	inline GPoint3  getP() { return P; }
	
	// 获取射线方向
	inline GVector3 getV() { return V; }

	// 获取参数t时的点
	inline GPoint3 getPoint(const float& t) { return P + V.normalize()*t; }


private:

	GPoint3 P;	// 起点
	GVector3 V; // 方向
};