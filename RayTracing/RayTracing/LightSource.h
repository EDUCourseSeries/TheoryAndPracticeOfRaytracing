#pragma once
#include "GMath.h"

//------------------------------------------ CLightSource ------------------------------------------//

//所有光源类的基类
class CLightSource
{
public:
	CLightSource();
	~CLightSource();

protected:

	// 光源位置
	GVector3 Position;

	// 光源环境光系数 Ambinet
	GVector3 Ka;

	// 光源漫反射系数 Diffuse
	GVector3 Kd;

	// 光源镜面反射系数 Specular
	GVector3 Ks;

public:
	
	inline void setPosition(const GVector3& _position) { Position = _position; }
	inline void setKa(const GVector3& _Ka) { Ka = _Ka; }
	inline void setKd(const GVector3& _Ka) { Ka = _Ka; }
	inline void setKs(const GVector3& _Ka) { Ka = _Ka; }

	inline GVector3 getPosition() { return Position; }
	inline GVector3 getKa() { return Ka; }
	inline GVector3 getKd() { return Kd; }
	inline GVector3 getKs() { return Ks; }


	/*
	*	\brief	计算该光源在物体表面任意一点产生的环境反射
	*
	*	\param _Ka 物体表面的环境反射系数
	*/
	virtual GVector3 EvalAmbient(const GVector3& _Ka) = 0;

	/*
	*	\brief	计算该光源在物体表面任意一点产生的漫反射
	*
	*	\param _n 该点的法线
	*	\param _l 光源的方向
	*	\param _Kd 物体表面的漫反射系数
	*/
	virtual GVector3 EvalDiffuse(const GVector3& _n, const GVector3& _l, const GVector3& _Kd) = 0;

	/*
	*	\brief	计算该光源在物体表面任意一点产生的镜面反射
	*
	*	\param _n 该点的法线
	*	\param _l 光源的方向
	*	\param _v 观察方向
	*	\param _Ka 物体表面的环境反射系数
	*	\param _shininess 该点的光滑程度
	*/
	virtual GVector3 EvalSpecluar(const GVector3& _n, const GVector3& _l, const GVector3& _v, const GVector3& _Ks, const float& _shininess) = 0;
};



//------------------------------------------ CDirectionalLight ------------------------------------------//

// 平行光源类
class CDirectionalLight : public CLightSource
{
public:

	/*
	*	\brief	创建一个平行光源
	*
	*	\param _position  光源的位置
	*	\param _direction 光源的方向
	*	\param _Ka 光源环境光系数
	*	\param _Kd 光源漫反射系数
	*	\param _Ks 光源镜面反射系数
	*/
	CDirectionalLight(const GVector3& _position, const GVector3& _direction, const GVector3& _Ka = GVector3(1.0f, 1.0f, 1.0f), const GVector3& _Kd = GVector3(1.0f, 1.0f, 1.0f), const GVector3& _Ks = GVector3(1.0f, 1.0f, 1.0f));

	~CDirectionalLight();

private:
	GVector3 Direction;

public:

	inline void setDirection(const GVector3& _direction) { Direction = _direction; }
	inline GVector3 getDirection() { return Direction; }

public:
	GVector3 EvalAmbient(const GVector3& _Ka);
	GVector3 EvalDiffuse(const GVector3& _n, const GVector3& _l, const GVector3& _Kd);
	GVector3 EvalSpecluar(const GVector3& _n, const GVector3& _l, const GVector3& _v, const GVector3& _Ks, const float& _shininess);
};