#pragma once
#include "CRay.h"
#include "CScene.h"



class CTracer
{
public:
	static GVector3 Trace(CScene* _scene, CRay R, int _traceDepth);
	static int TraceDepth;
};