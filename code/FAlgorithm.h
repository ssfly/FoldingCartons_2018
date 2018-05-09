#pragma once
#include"Lines.h"
#include"Vertex.h"
#include"halfedge.h"
#include<Eigen\Dense>
#include<Eigen\Eigen>
#include "Solver.h"
#include "Constraint.h"

using namespace std;

class FAlgorithm
{
public:
	FAlgorithm();
	~FAlgorithm();

	void LoadParameters();

	float computeDistance(MyMesh m, MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2);

	void InitializeShapeConstrain(MyMesh m);
	MyMesh ShapeOptimization(MyMesh mesh, vector<MyMesh::VertexHandle> vhlist);
	MyMesh ShapeOptimization(MyMesh mesh, vector<MyMesh::FaceHandle> fhlist);
	MyMesh ShapeOptimization(MyMesh mesh, vector<MyMesh::VertexHandle> vhlist, vector<MyMesh::VertexHandle> vhlist2);

	vector<MyMesh::VertexHandle> symmetricVertex(MyMesh mesh, MyMesh::VertexHandle &vh);
	//bool compare(MyMesh::VertexHandle v1h, MyMesh::VertexHandle v2h);
	vector<vector<MyMesh::VertexHandle>> symmetricVertexPair(MyMesh mesh, MyMesh::VertexHandle &v1h, MyMesh::VertexHandle &v2h);


	vector<vector<MyMesh::VertexHandle>> symmetryDetection(MyMesh m, vector<MyMesh::VertexHandle> vhlist);
	vector<vector<MyMesh::VertexHandle>> modifiedSymmetryDetection(MyMesh m, vector<MyMesh::VertexHandle> vhlist);

	vector<vector<MyMesh::VertexHandle>> autoDetectionbySameVertex(MyMesh m);
	vector<vector<MyMesh::VertexHandle>> autoDetection(MyMesh m);
	vector<vector<MyMesh::VertexHandle>> autoDetection2(MyMesh m);
	vector<vector<MyMesh::FaceHandle>> autoDectionofFace(MyMesh m);

	vector<int> GetClassForSameVertexin3D(MyMesh mesh);
	vector<int> GetVertexClass();

	MyMesh LayoutOptimization(MyMesh FMesh, MyMesh SMesh, MyMesh::VertexHandle & v1h, MyMesh::VertexHandle & v2h);
	MyMesh LayoutOptimization(MyMesh FMesh, MyMesh SMesh, vector<MyMesh::VertexHandle> vhs);

private:
	vector<MyMesh::VertexHandle> m_vertices;
	ShapeOp::Solver s;
	ShapeOp::Matrix3X p;
	vector<int> labelForSameVertexin3D; //相同的标识意味着在3D MODEL中应该被当做同一个点

	

	bool isInSameFace(MyMesh m, MyMesh::VertexHandle v1h, MyMesh::VertexHandle v2h);
	MyMesh JointMove(MyMesh SMesh, MyMesh::VertexHandle & v1h, MyMesh::VertexHandle & v2h);
};

