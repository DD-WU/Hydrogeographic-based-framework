#pragma once
class Abstract_Feature_Interaction_Hub
{
public:
	virtual void beforeInit() ;
	virtual void init() = 0;
	virtual void afterInit();
	virtual void beforeUpdate();
	virtual void update() = 0;
	virtual void afterUpdate();
	virtual void beforeFinalize();
	virtual void finalize() = 0;
	virtual void afterFinalize();
};






