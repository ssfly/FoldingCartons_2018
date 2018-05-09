#pragma once
#include "Vertex.h"
class LineSegment
{
public:
	LineSegment();
	LineSegment(Vertex v1, Vertex v2);
	~LineSegment();

	Vertex GetV1();
	void SetV1(Vertex v1);
	void SetV1(float v_x, float v_y, float v_z);
	Vertex GetV2();
	void SetV2(Vertex v2);
	void SetV2(float v_x, float v_y, float v_z);

	void SetV1SeparateId(int f);
	int GetV1SeparateId();

	void SetV2SeparateId(int f);
	int GetV2SeparateId();

	int GetId();
	int SetId(int id);

	int GetIsDash();
	int SetIsDash(int isdash);

	int GetInWhichLoop();
	int SetInWhichLoop(int id);

	int Floatcmp(float a, float b);
	float Dot(float x1, float y1, float x2, float y2);
	float Cross(float x1, float y1, float x2, float y2);
	float NormalizedCross(float x1, float y1, float x2, float y2);
	int PointOnLine(Vertex &a, Vertex &b, Vertex &c);
	float ABCrossAC(Vertex &a, Vertex &b, Vertex &c);
	int IntersectionByTwoLines(LineSegment &l, Vertex &v);

	//判断两个线段是否为同一条
	bool operator ==(const LineSegment& rhs) const
	{
		if ((_V1 == rhs._V1 && _V2 == rhs._V2) || (_V1 == rhs._V2 && _V2 == rhs._V1))
			return true;
		else
			return false;
	}

	double P2DLength();

private:
	Vertex _V1;
	Vertex _V2;

	int _InWhichLoop;
	int _Id;
	int _IsDash;
};

