#pragma once
#include"Lines.h"
#include"Vertex.h"
#include"halfedge.h"
#include<vector>

using namespace std;
static GLfloat colorSpace[][3] = {
	{ 241. / 255., 109. / 255., 122. / 255. },{ 241. / 255., 204. / 255.,184. / 255. },{ 219. / 255., 136. / 255., 120. / 255. },
	{ 183. / 255., 210. / 255., 141. / 255. },{ 184. / 255., 241. / 255., 237. / 255. },{ 241. / 255., 241. / 255., 184. / 255. },
	{ 231. / 255., 218. / 255., 201. / 255. },{ 241. / 255., 184. / 255., 228. / 255. },{ 243. / 255., 214. / 255., 78. / 255. },
	{ 254. / 255., 249. / 255., 69. / 255. } ,{ 227. / 255., 124. / 255., 91. / 255. },{ 255. / 255., 130. / 255., 64. / 255. },
	{ 207. / 255., 136. / 255., 136. / 255. },{ 226. / 255., 158. / 255., 75. / 255. },{ 183. / 255., 210. / 255., 141. / 255. }
};
static int colorNum = 15;
class FlatMesh
{
public:
	FlatMesh();
	~FlatMesh();

	void createMesh(string filepath, MyMesh &m);
	



private:

	void UpdateMaxMinValue(const float x, float& minx, float& maxx);
	void ReadLinesFromTxt(string filepath);
	//切割线段
	void CutOffLines();

	void generateVertices(Vertex &v);
	int findVertexId(Vertex &v);
	void generateHalfedgePair(int src, int dst);
	void initializeData();

	bool HasSeparateVertex(vector<OpenMesh::VertexHandle> &f1);
	bool HasSameSeparateVertices(MyMesh &m, vector<OpenMesh::VertexHandle> &f1, vector<OpenMesh::VertexHandle> &f2);

	void clear_halfedges();

	vector<LineSegment> m_lineList;
	vector<Vertex> m_vertices;
	vector<Halfedge*> m_halfedges;
	vector<int> m_whichLineofHalfedge;  //记录每条halfedge对应哪条line
	MyMesh m_flatMesh;

};

