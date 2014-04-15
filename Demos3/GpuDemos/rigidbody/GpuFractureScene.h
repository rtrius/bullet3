#ifndef GPU_FRACTURE_SCENE_H
#define GPU_FRACTURE_SCENE_H

#include "GpuConvexScene.h"

class GpuFractureScene : public GpuBoxPlaneScene
{
public:
	GpuFractureScene(){}
	virtual ~GpuFractureScene(){}
	virtual const char* getName()
	{
		return "Rigid Fracture";
	}

	static GpuDemo* MyCreateFunc()
	{
		GpuDemo* demo = new GpuFractureScene;
		return demo;
	}

	virtual void setupScene(const ConstructionInfo& ci);
	
	virtual void renderScene();
	
	virtual void clientMoveAndDisplay();
	
	virtual bool mouseMoveCallback(float x, float y);
	virtual bool mouseButtonCallback(int button, int state, float x, float y);
};

#endif //GPU_FRACTURE_SCENE_H
