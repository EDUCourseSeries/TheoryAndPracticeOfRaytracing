#include "Tracer.h"

int CTracer::TraceDepth = 1;

GVector3 CTracer::Trace(CScene * _scene, CRay R, int _traceDepth)
{
	GVector3 color;
	float shade = 1.0f;

	// 遍历每一个物体
	for (int k = 0; k < _scene->GetObjectsNum(); k++)
	{
		if (_scene->GetObject(k)->isIntersect(R))
		{
			// 计算交点
			GVector3 P = _scene->GetObject(k)->getIntersection(R);
			
			// 计算法线
			GVector3 N = _scene->GetObject(k)->getN(P);
			N.normalize();

			// 遍历每一个光源
			for (int m = 0; m < _scene->GetLightsNum(); m++)
			{
				GVector3 L = _scene->GetLight(m)->getPosition() - P;
				L.normalize();

				// 判断P点是否在阴影中
				for (int n = k + 1; n < _scene->GetObjectsNum(); n++)
				{
					if (_scene->GetObject(n)->isIntersect(CRay(_scene->GetLight(m)->getPosition(), -L)));
					{
						shade = .0f;
						break;
					}
				}

				// 计算phong光照模型
				if (shade > .0f)
				{
					GVector3 V = _scene->GetCamPos() - P;
					V.normalize();

					GVector3 ambient = _scene->GetLight(m)->EvalAmbient(_scene->GetObject(k)->getKa());
					GVector3 diffuse = _scene->GetLight(m)->EvalDiffuse(N,L,_scene->GetObject(k)->getKd());
					GVector3 specular = _scene->GetLight(m)->EvalSpecluar(N, L,V, _scene->GetObject(k)->getKs(), _scene->GetObject(k)->getShininess());

					color += ambient + diffuse + specular;
				}
			}
			
			if (shade > .0f)
			{
				if (TraceDepth != _traceDepth)
				{
					CRay Reflect(P, R.getV() - 2.0f*(R.getV()*N)*N);
					GVector3 c = Trace(_scene, Reflect, ++_traceDepth);
					color += c;
				}
			}
		}
	}

	return color*shade;

	
}
