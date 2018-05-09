#include "Vertex.h"
#include"Lines.h"
#include"halfedge.h"


Vertex::Vertex()
{
	_X = 0;
	_Y = 0;
	_Z = 0;
	 
	isSeparateIn3D = 1 ;  //初始化为都需要分离，后续遇到虚线，其两端的顶点值改为0
}


Vertex::~Vertex()
{
}

Vertex::Vertex(float x, float y, float z)
{
	_X = x;
	_Y = y;
	_Z = z;
}

float Vertex::GetX()
{
	return _X;
}

void Vertex::SetX(float x)
{
	_X = x;
}

float Vertex::GetY()
{
	return _Y;
}

void Vertex::SetY(float y)
{
	_Y = y;
}

float Vertex::GetZ()
{
	return _Z;
}

void Vertex::SetZ(float z)
{
	_Z = z;
}

void Vertex::SetIsSeparateIn3D(int f)
{
	isSeparateIn3D = f;
}

int Vertex::GetIsSeparateIn3D()
{
	return isSeparateIn3D;
}

Vertex Vertex::Cross(Vertex a)
{
	Vertex cc;
	float Xa = a.GetX();
	float Ya = a.GetY();
	float Za = a.GetZ();

	float u = _Y * Za - Ya * _Z;
	float v = _Z * Xa - Za * _X;
	float w = _X * Ya - Xa * _Y;

	float s = sqrtf(u*u + v*v + w*w);

	cc.SetX(u / s);
	cc.SetY(v / s);
	cc.SetZ(w / s);

	return cc;


}


float Vertex::DistanceBetweenTwoVertices(Vertex a)
{
	float ax = a.GetX();
	float ay = a.GetY();
	float az = a.GetZ();

	float bx = _X;
	float by = _Y;
	float bz = _Z;

	float dis = (ax - bx) * (ax - bx) + (ay - by) * (ay - by) + (az - bz) * (az - bz);
	dis = sqrtf(dis);
	return dis;
}

bool Vertex::IsThreePointsInALine(Vertex a, Vertex b)
{
	// using distance
	float dis1 = this->DistanceBetweenTwoVertices(a);
	float dis2 = this->DistanceBetweenTwoVertices(b);
	float dis3 = a.DistanceBetweenTwoVertices(b);

	if (abs(dis1 + dis2 - dis3) < g_parameters.MyEpsilon || abs(dis1 + dis3 - dis2) < g_parameters.MyEpsilon || abs(dis2 + dis3 - dis1) < g_parameters.MyEpsilon)
		return true;
	else
		return false;
}

void Vertex::sort_o_halfedges()
{
	std::sort(Ohalfedges.begin(), Ohalfedges.end(), [](Halfedge *a, Halfedge *b) {
		return ((a->angle()) < (b->angle()));
	});
	for (int i = 0, N = Ohalfedges.size(); i<N; ++i)
	{
		Halfedge *prevh = Ohalfedges[i]->opposite();
		Halfedge *nexth = Ohalfedges[(i + 1) % N];
		prevh->set_next(nexth);
		nexth->set_prev(prevh);
	}
}
