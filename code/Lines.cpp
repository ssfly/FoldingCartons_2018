#include "Lines.h"
#include <math.h>


LineSegment::LineSegment()
{
	_InWhichLoop = -1;
}

LineSegment::LineSegment(Vertex v1, Vertex v2)
{
	_V1 = v1;
	_V2 = v2;
}


LineSegment::~LineSegment()
{
}

Vertex LineSegment::GetV1()
{
	return _V1;
}

void LineSegment::SetV1(Vertex v1)
{
	_V1 = v1;
}

void LineSegment::SetV1(float v_x, float v_y, float v_z)
{
	_V1.SetX(v_x);
	_V1.SetY(v_y);
	_V1.SetZ(v_z);
}


Vertex LineSegment::GetV2()
{
	return _V2;
}

void LineSegment::SetV2(Vertex v2)
{
	_V2 = v2;
}

void LineSegment::SetV2(float v_x, float v_y, float v_z)
{
	_V2.SetX(v_x);
	_V2.SetY(v_y);
	_V2.SetZ(v_z);
}

void LineSegment::SetV1SeparateId(int f)
{
	_V1.SetIsSeparateIn3D(f);
}

int LineSegment::GetV1SeparateId()
{
	return _V1.GetIsSeparateIn3D();
}

void LineSegment::SetV2SeparateId(int f)
{
	_V2.SetIsSeparateIn3D(f);
}

int LineSegment::GetV2SeparateId()
{
	return _V2.GetIsSeparateIn3D();;
}




int LineSegment::GetId()
{
	return _Id;
}

int LineSegment::SetId(int id)
{
	_Id = id;
	return _Id;
}

int LineSegment::GetIsDash()
{
	return _IsDash;
}

int LineSegment::SetIsDash(int isdash)
{
	_IsDash = isdash;
	return _IsDash;
}

int LineSegment::GetInWhichLoop()
{
	return _InWhichLoop;
}

int LineSegment::SetInWhichLoop(int id)
{
	_InWhichLoop = id;
	return _InWhichLoop;
}

int LineSegment::Floatcmp(float a, float b)
{
	if (fabs(a - b) <= g_parameters.MyEpsilon) return 0;
	if (a > b) return 1;
	else return -1;
}

float LineSegment::Dot(float x1, float y1, float x2, float y2)
{
	//点积判断点是否在线段上
	return x1*x2 + y1*y2;
}

float LineSegment::Cross(float x1, float y1, float x2, float y2)
{

	return (x1*y2 - x2*y1);
}

float LineSegment::NormalizedCross(float x1, float y1, float x2, float y2)
{
	return(x1*y2 - x2*y1) / (sqrt(x1*x1 + y1*y1)*sqrt(x2*x2 + y2*y2));
}

int LineSegment::PointOnLine(Vertex &a, Vertex &b, Vertex &c)
{
	//求a,b,c是否共线，>0不在，=0与端点重合，<0在。
	float a_x = a.GetX();
	float a_y = a.GetY();
	float b_x = b.GetX();
	float b_y = b.GetY();
	float c_x = c.GetX();
	float c_y = c.GetY();
	return Floatcmp(Dot(b_x - a_x, b_y - a_y, c_x - a_x, c_y - a_y), 0);
}

float LineSegment::ABCrossAC(Vertex &a, Vertex &b, Vertex &c)
{
	//ab与ac的叉积
	float a_x = a.GetX();
	float a_y = a.GetY();
	float b_x = b.GetX();
	float b_y = b.GetY();
	float c_x = c.GetX();
	float c_y = c.GetY();
	return NormalizedCross(b_x - a_x, b_y - a_y, c_x - a_x, c_y - a_y);
}

int LineSegment::IntersectionByTwoLines(LineSegment &l, Vertex &v)
{
	//Z = 0平面
	//求L1是否与L2相交，交点为p。1规范相交，2,3,4,5交点是一线段的端点，-1不相交
	float s1, s2;
	int d1, d2, d3, d4;
	float a_x = _V1.GetX();
	float a_y = _V1.GetY();
	float b_x = _V2.GetX();
	float b_y = _V2.GetY();
	float c_x = l.GetV1().GetX();
	float c_y = l.GetV1().GetY();
	float d_x = l.GetV2().GetX();
	float d_y = l.GetV2().GetY();

	d1 = Floatcmp(ABCrossAC(_V1, _V2, l.GetV1()), 0);
	d2 = Floatcmp(ABCrossAC(_V1, _V2, l.GetV2()), 0);
	d3 = Floatcmp(ABCrossAC(l.GetV1(), l.GetV2(), _V1), 0);
	d4 = Floatcmp(ABCrossAC(l.GetV1(), l.GetV2(), _V2), 0);

	//如果规范相交则求交点
	if ((d1^d2) == -2 && (d3^d4) == -2)
	{
		s1 = Cross(b_x - a_x, b_y - a_y, c_x - a_x, c_y - a_y);
		s2 = Cross(b_x - a_x, b_y - a_y, d_x - a_x, d_y - a_y);
		v.SetX((c_x*s2 - d_x*s1) / (s2 - s1));
		v.SetY((c_y*s2 - d_y*s1) / (s2 - s1));
		return 1;
	}

	//如果不规范相交
	if (d1 == 0 && PointOnLine(l.GetV1(), _V1, _V2) <= 0)
	{
		v = l.GetV1();
		return 2;
	}
	if (d2 == 0 && PointOnLine(l.GetV2(), _V1, _V2) <= 0)
	{
		v = l.GetV2();
		return 3;
	}
	if (d3 == 0 && PointOnLine(_V1, l.GetV1(), l.GetV2()) <= 0)
	{
		v = _V1;
		return 4;
	}
	if (d4 == 0 && PointOnLine(_V2, l.GetV1(), l.GetV2()) <= 0)
	{
		v = _V2;
		return 5;
	}

	//如果不相交
	return -1;
}

double LineSegment::P2DLength()
{
	double dis;
	dis = (_V1.GetX() - _V2.GetX())*(_V1.GetX() - _V2.GetX()) + (_V1.GetY() - _V2.GetY())*(_V1.GetY() - _V2.GetY());
	return sqrt(dis);
}