#include "Tracer.h"

int main(void)
{
	CScene *scene = new CScene();

	CDirectionalLight *DL1 = new CDirectionalLight(GVector3(), GVector3(.0f, .0f, 1.0f));
	CPlane *Pi1 = new CPlane();

	scene->AddLight(DL1);
	scene->AddObject(Pi1);

	CRay r;

	CTracer::Trace(scene, r, 1);
	return 1;

}