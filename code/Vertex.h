#ifndef VERTEX_H
#define VERTEX_H

#include "Parameters.h"
#include <math.h>
#include "../glut/glut.h"
#include"MyMesh.h"

class Halfedge;
class LineSegment;

class Vertex
{
public:
	Vertex();
	~Vertex();

	Vertex(float x, float y, float z);

	float GetX();
	void SetX(float x);
	float GetY();
	void SetY(float y);
	float GetZ();
	void SetZ(float z);

	void SetIsSeparateIn3D(int f);
	int GetIsSeparateIn3D();


	bool operator ==(const Vertex& rhs) const
	{
		if (fabs(_X-rhs._X)< g_parameters.ThresholdForSeparateLines && fabs(_Y-rhs._Y)< g_parameters.ThresholdForSeparateLines && fabs(_Z-rhs._Z)< g_parameters.ThresholdForSeparateLines)
			return true;
		else
			return false;
	}

	//void SetParameters(Parameters para);
	
	Vertex Cross(Vertex a);

	float DistanceBetweenTwoVertices(Vertex a);
	bool IsThreePointsInALine(Vertex a, Vertex b);

	void add_o_halfedge(Halfedge *h) { Ohalfedges.push_back(h); }
	void add_i_halfedge(Halfedge *h) { Ihalfedges.push_back(h); }
	vector<Halfedge*> get_o_halfedges() { return Ohalfedges; }

	void sort_o_halfedges();
	void clear_halfedges() { Ohalfedges.clear(); Ihalfedges.clear(); }

	void set_idx(int id) { m_id = id; };
	int idx() { return m_id; };

private:
	float _X;
	float _Y;
	float _Z;

	int m_id;

	int isSeparateIn3D;

	vector<Halfedge *> Ohalfedges;
	vector<Halfedge *> Ihalfedges;
};

#endif