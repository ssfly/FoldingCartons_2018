#pragma once
#include"Lines.h"
#include"Vertex.h"
#include"halfedge.h"
#include<unordered_set>

using namespace std;

class StereoMesh
{
public:
	StereoMesh();
	StereoMesh(MyMesh sm);
	~StereoMesh();

	float computeRotateAngle(float angle, MyMesh::EdgeHandle edge, int f1, int f2);
	void RotateArbitraryAxis(float angle, MyMesh::EdgeHandle edge, float rtmatrix[4][4], int start_p);
	void changeVertexCoordinate(float **rtmatrix, int f_id);
	void changeVertexCoordinate(float ** rtmatrix, MyMesh::VertexHandle vh);
	unordered_set<MyMesh::VertexHandle> vertexOnUnvisitedAdjface(int fid, vector<int> visited);
	MyMesh folding();
	MyMesh folding(int hasvisited, float finalAngle);

private:
	vector<int> m_f1_idx;
	vector<int> m_f2_idx;
	vector<int> m_v1_idx;
	vector<int> m_v2_idx;
	MyMesh m_sm;
	MyMesh m_tm;
	vector<MyMesh::FaceHandle> m_fhlist;
	vector<MyMesh::VertexHandle> m_vhlist;

	float ***m_RTMatrix;
	float m_initialAngle;  //folding后所有邻接面的角度值

	double calFaceArea(MyMesh sm, MyMesh::FaceHandle fh);
	int MaxFaceIdx(MyMesh sm);

	//Animation
	MyMesh amesh;
	int age_;
};

