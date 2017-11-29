#pragma once
#include "GMath.h"

//------------------------------------------ CLightSource ------------------------------------------//

//���й�Դ��Ļ���
class CLightSource
{
public:
	CLightSource();
	~CLightSource();

protected:

	// ��Դλ��
	GVector3 Position;

	// ��Դ������ϵ�� Ambinet
	GVector3 Ka;

	// ��Դ������ϵ�� Diffuse
	GVector3 Kd;

	// ��Դ���淴��ϵ�� Specular
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
	*	\brief	����ù�Դ�������������һ������Ļ�������
	*
	*	\param _Ka �������Ļ�������ϵ��
	*/
	virtual GVector3 EvalAmbient(const GVector3& _Ka) = 0;

	/*
	*	\brief	����ù�Դ�������������һ�������������
	*
	*	\param _n �õ�ķ���
	*	\param _l ��Դ�ķ���
	*	\param _Kd ��������������ϵ��
	*/
	virtual GVector3 EvalDiffuse(const GVector3& _n, const GVector3& _l, const GVector3& _Kd) = 0;

	/*
	*	\brief	����ù�Դ�������������һ������ľ��淴��
	*
	*	\param _n �õ�ķ���
	*	\param _l ��Դ�ķ���
	*	\param _v �۲췽��
	*	\param _Ka �������Ļ�������ϵ��
	*	\param _shininess �õ�Ĺ⻬�̶�
	*/
	virtual GVector3 EvalSpecluar(const GVector3& _n, const GVector3& _l, const GVector3& _v, const GVector3& _Ks, const float& _shininess) = 0;
};



//------------------------------------------ CDirectionalLight ------------------------------------------//

// ƽ�й�Դ��
class CDirectionalLight : public CLightSource
{
public:

	/*
	*	\brief	����һ��ƽ�й�Դ
	*
	*	\param _position  ��Դ��λ��
	*	\param _direction ��Դ�ķ���
	*	\param _Ka ��Դ������ϵ��
	*	\param _Kd ��Դ������ϵ��
	*	\param _Ks ��Դ���淴��ϵ��
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