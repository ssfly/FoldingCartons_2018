#pragma once

#include <qgl.h>
#include <gl\GLU.h>
#include <QMouseEvent>
#include <QGLViewer/qglviewer.h>
#include <vector>
#include <qobject.h>
#include <fstream>
#include <cv.h>
#include "mymesh.h"

#define E 1

using namespace std;
using namespace cv;

static GLfloat colorSpace_[][3] = {
	{ 241. / 255., 109. / 255., 122. / 255. },{ 241. / 255., 204. / 255.,184. / 255. },{ 219. / 255., 136. / 255., 120. / 255. },
	{ 183. / 255., 210. / 255., 141. / 255. },{ 184. / 255., 241. / 255., 237. / 255. },{ 241. / 255., 241. / 255., 184. / 255. },
	{ 231. / 255., 218. / 255., 201. / 255. },{ 241. / 255., 184. / 255., 228. / 255. },{ 243. / 255., 214. / 255., 78. / 255. },
	{ 254. / 255., 249. / 255., 69. / 255. } ,{ 227. / 255., 124. / 255., 91. / 255. },{ 255. / 255., 130. / 255., 64. / 255. },
	{ 207. / 255., 136. / 255., 136. / 255. },{ 226. / 255., 158. / 255., 75. / 255. },{ 220. / 255., 255. / 255., 147. / 255. },
	{ 242. / 255., 222. / 255., 189. / 255. },{ 207. / 255., 136. / 255., 136. / 255. },{ 254. / 255., 207. / 255., 69. / 255. },
	{ 184. / 255., 241. / 255., 204. / 255. },{ 254. / 255., 151. / 255., 120. / 255. },{ 255. / 255., 229. / 255., 67. / 255. }
};

static int colorNum_ = 21;

class Viewer : public QGLViewer, public QObject
{
	Q_OBJECT
public:
	Viewer() : selectionmode_(NONE), type_(NONE_), v_data(NULL) { isCandidate = false; isUsed = false; needUpdateColor = true; showMergeVertices = false;  };
	//~Viewer();

	enum SelectionMode
	{
		NONE,
		ADD,
		REMOVE,
		EDIT
	};
	enum Type
	{
		NONE_,
		VERTEX,
		FACE,
		EDGE
	};
	enum EditMode
	{
		FREE,
		ONEFACE,
		TWOFACE
	};
	virtual void draw();
	virtual void init();
	virtual void animate();
	virtual void drawWithNames();
	virtual void endSelection(const QPoint &);

	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void mouseReleaseEvent(QMouseEvent *e);
	virtual void mouseDoubleClickEvent(QMouseEvent *e);

	void clearViewer();
	void setColorIdxAs(vector<int> coloridx);
	vector<int> getColorIdx();

	void getEdgeVertex();

	MyMesh mesh;
	SelectionMode selectionmode_;
	Type type_;
	QList<int> selected_vertices;
	QList<int> selected_edges;
	QList<int> selected_faces;
	bool isCandidate;
	vector<int> mergeVertices;
	vector<int> symVertices;

	vector<int> mergeFaces;
	int viewIdx;
	bool isUsed;
	bool needUpdateColor;
	bool showMergeVertices;

	vector<int> getEdgeVertex0() { return edgeVertex0; };
	vector<int> getEdgeVertex1() { return edgeVertex1; };

	vector<int> m_label;   ////相同的标识意味着在3D MODEL中应该被当做同一个点

	bool isAnimation = false;  //用于动画
	MyMesh flat_mesh;
	int hasvisited ; //visited edge
	float curAngle ; //当前访问的边的角度

signals:
	void isDoubleClicked(int);
	void updatelayout();
	void computeLengthError();

private:
	vector<int> edgeVertex0;
	vector<int> edgeVertex1;

	void drawSelectionRectangle() const;
	void drawMoveLine() const;

	void addIdToSelection(int id);
	void removeIdFromSelection(int id);

	void drawVertexWithName();
	void drawFaceWithName();
	void drawEdgeWithName();

	QRect rectangle_;
	QLine line_;

	void drawVertex(int id);
	void drawEdge(int id);
	void drawFace(int id);

	void transferPos(int x, int y, qreal *pos0, qreal *pos1, qreal *pos2);
	vector<float> CalPoint(vector<float> planeVector, vector<float>  planePoint, vector<float>  lineVector, vector<float>  linePoint);
	vector<int> findFace(MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2);
	vector<float> calFaceNormal(MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2);
	void updateEdit();
	EditMode editmode_;

	void updateMeshColor();
	vector<int> faceColorIdx;

	GLUtesselator * tesser();
	vector<vector<vector<GLdouble>>> v_data;
	vector<int> fvSize;
	void updateData();

	

};

