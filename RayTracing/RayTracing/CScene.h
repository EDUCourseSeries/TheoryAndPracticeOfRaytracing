#pragma once
#include <vector>
#include "LightSource.h"
#include "CObject.h"

class CScene
{
public:
	CScene();

	inline void AddLight(CLightSource* _ptr_light) { Lights.push_back(_ptr_light); }
	inline void AddObject(CObject* _ptr_obj) { Objects.push_back(_ptr_obj); }

	inline CLightSource* GetLight(const int& idx) 
	{ 
		assert(idx >= 0);
		return Lights[idx]; 
	}

	inline CObject* GetObject(const int& idx)
	{
		assert(idx >= 0);
		return Objects[idx];
	}

	inline void SetCamPos(const GVector3& _pos) { CamPos = _pos; }
	inline void SetCamDir(const GVector3& _dir) { CamDir = _dir; }
	inline GVector3 GetCamPos() { return CamPos; }
	inline GVector3 GetCamDir() { return CamDir; }
	inline int GetObjectsNum() { return Objects.size(); }
	inline int GetLightsNum() { return Lights.size(); }


private:

	vector<CLightSource*> Lights;
	vector<CObject*> Objects;
	
	GVector3 CamPos;
	GVector3 CamDir;

};