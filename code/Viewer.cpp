#include "Viewer.h"
#include "StereoMesh.h"
#include <qdebug.h>


void Viewer::getEdgeVertex()
{
	int n = mesh.n_edges();
	vector<int> vtmp0(n, -1);
	vector<int> vtmp1(n, -1);
	//initVector(edgeVertex0, n);
	//initVector(edgeVertex1, n);
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		int vidx = (*v_it).idx();
		for (MyMesh::VertexEdgeIter ve_it = mesh.ve_iter(*v_it); ve_it.is_valid(); ve_it++)
		{
			int eidx = (*ve_it).idx();
			if (vtmp0[eidx] == -1)
				vtmp0[eidx] = vidx;
			else vtmp1[eidx] = vidx;
		}
	}
	edgeVertex0 = vtmp0;
	edgeVertex1 = vtmp1;

	//for (int i = 0; i < fvSize.size(); i++)
	//{
	//	for (int j = 0; j < fvSize[i]; j++)
	//		//if (v_data[i][j] != NULL)
	//			//delete[]v_data[i][j];
	//	delete[]v_data[i];
	//}
	//ofstream out("out.txt");
	//for (int i = 0; i < fvSize.size(); i++)
	//{
	//	for (int j = 0; j < fvSize[i]; j++)
	//	{
	//		out << v_data[i][j][0] << " " << v_data[i][j][1] << " " << v_data[i][j][2] << endl;
	//	}
	//}
	//out.close();
	/*if (v_data !=NULL) delete[]v_data;	*/
	
	//vector<vector<GLdouble>> tmpList;
	//vector<GLdouble> tmp(3, 0);

	//for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	//{
	//	int count = 0;
	//	for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
	//		count++;
	//	fvSize.push_back(count);
	//}

	//for (int i = 0; i < fvSize.size(); i++)
	//{
	//	for (int j = 0; j < fvSize[i]; j++)
	//		tmpList.push_back(tmp);

	//	v_data.push_back(tmpList);
	//	tmpList.clear();
	//}

	//for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	//{
	//	int count = 0;
	//	for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
	//	{
	//		v_data[(*f_it).idx()][count][0] = mesh.point(*fv_it)[0];
	//		v_data[(*f_it).idx()][count][1] = mesh.point(*fv_it)[1];
	//		v_data[(*f_it).idx()][count][2] = mesh.point(*fv_it)[2];
	//		count++;
	//	}
	//}
	updateData();
	faceColorIdx.resize(mesh.n_faces());
	for (int i = 0; i < mesh.n_faces(); i++)
		faceColorIdx[i] = i;
	if (needUpdateColor) updateMeshColor();
}

void Viewer::updateData()
{
	vector<vector<GLdouble>> tmpList;
	vector<GLdouble> tmp(3, 0);

	fvSize.clear();
	for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		int count = 0;
		for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
			count++;
		fvSize.push_back(count);
	}
	v_data.clear();
	for (int i = 0; i < fvSize.size(); i++)
	{
		for (int j = 0; j < fvSize[i]; j++)
			tmpList.push_back(tmp);

		v_data.push_back(tmpList);
		tmpList.clear();
	}

	for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		int count = 0;
		for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			v_data[(*f_it).idx()][count][0] = mesh.point(*fv_it)[0];
			v_data[(*f_it).idx()][count][1] = mesh.point(*fv_it)[1];
			v_data[(*f_it).idx()][count][2] = mesh.point(*fv_it)[2];
			count++;
		}
	}
}
void Viewer::updateMeshColor()
{
	mesh.request_face_normals();
	for (MyMesh::FaceIter f_it0 = mesh.faces_begin(); (f_it0 + 1) != mesh.faces_end(); f_it0++)
	{
		MyMesh::FaceVertexIter fv0 = mesh.fv_iter(*f_it0);
		for (MyMesh::FaceIter f_it1 = mesh.faces_begin() + 1; f_it1 != mesh.faces_end(); f_it1++)
		{
			bool areSame = true;
			for (MyMesh::FaceVertexIter fv = mesh.fv_iter(*f_it1); fv.is_valid(); fv++)
			{
				float dist_x = (mesh.point(*fv)[0] - mesh.point(*fv0)[0]) * mesh.normal(*f_it0)[0];
				float dist_y = (mesh.point(*fv)[1] - mesh.point(*fv0)[1]) * mesh.normal(*f_it0)[1];
				float dist_z = (mesh.point(*fv)[2] - mesh.point(*fv0)[2]) * mesh.normal(*f_it0)[2];
				float dist = dist_x + dist_y + dist_z;
				if ( abs(dist) > E )
				{
					areSame = false;
					continue;
				}
			}
			//if (areSame && faceColorIdx[(*f_it1).idx()] != faceColorIdx[(*f_it0).idx()])
			//{
			//	//mesh.set_color((*f_it1), MyMesh::Color(mesh.color(*f_it0)[0], mesh.color(*f_it0)[1], mesh.color(*f_it0)[2]));
			//	faceColorIdx[(*f_it1).idx()] = faceColorIdx[(*f_it0).idx()];
			//}
		}
	}
}
void Viewer::drawVertex(int id)
{
	MyMesh::VertexHandle vh(id);
	Vec3d pos = { mesh.point(vh)[0],  mesh.point(vh)[1] ,mesh.point(vh)[2] };
	glPointSize(15);
	glBegin(GL_POINTS);
	glVertex3f(pos[0], pos[1], pos[2]);
	glEnd();
	//glPushMatrix();
	//glTranslatef(pos[0], pos[1], pos[2]);
	//GLUquadric *quad = gluNewQuadric();
	//gluSphere(quad, 4., 4., 4.);
	//glPopMatrix();
}

void Viewer::drawEdge(int id)
{
	MyMesh::EdgeHandle eh(id);
	//MyMesh::HalfedgeHandle ehh = mesh.halfedge_handle(eh, 0);
	//MyMesh::VertexHandle vh1 = mesh.from_vertex_handle(ehh);
	//MyMesh::VertexHandle vh2 = mesh.to_vertex_handle(ehh);
	MyMesh::VertexHandle vh1(edgeVertex0[id]);
	MyMesh::VertexHandle vh2(edgeVertex1[id]);
	glLineWidth(5);
	glBegin(GL_LINE_STRIP);
	glVertex3f(mesh.point(vh1)[0], mesh.point(vh1)[1], mesh.point(vh1)[2]);
	glVertex3f(mesh.point(vh2)[0], mesh.point(vh2)[1], mesh.point(vh2)[2]);
	glEnd();
}

void Viewer::drawFace(int id)
{
	MyMesh::FaceHandle fh(id);
	//mesh.set_color(fh, MyMesh::Color(0,0,0));

	glBegin(GL_POLYGON);
	for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it.is_valid(); fv_it++)
	{
		glVertex3f(mesh.point(*fv_it)[0], mesh.point(*fv_it)[1], mesh.point(*fv_it)[2] );
	}
	glEnd();
}
void Viewer::clearViewer()
{
	mesh.clear();
	edgeVertex0.clear();
	edgeVertex1.clear();
	mergeVertices.clear();
	//selectionmode_ = NONE;
	//type_ = NONE_;
	selected_edges.clear();
	selected_vertices.clear();
	selected_faces.clear();
	isCandidate = false;
	isUsed = false;
	/*if (!fvSize.empty())
	{
		for (int i = 0; i < fvSize.size(); i++)
		{
			for (int j = 0; j < fvSize[i]; j++)
				delete[]v_data[i][j];
			delete[]v_data[i];
		}
		delete[]v_data;
	}*/
	v_data.clear();
	fvSize.clear();
}

GLUtesselator * Viewer::tesser()
{
	GLUtesselator * tess;
	tess = gluNewTess();
	gluTessCallback(tess, GLU_TESS_BEGIN, (void(_stdcall *)()) &glBegin);
	gluTessCallback(tess, GLU_TESS_VERTEX, (void(_stdcall *)()) &glVertex3dv);
	gluTessCallback(tess, GLU_TESS_END, (void(_stdcall *)()) &glEnd);
	return tess;
}

vector<int> Viewer::getColorIdx()
{
	return faceColorIdx;
}

void Viewer::setColorIdxAs(vector<int> colorIdx)
{
	faceColorIdx.clear();
	faceColorIdx = colorIdx;
}

vector<float> Viewer::CalPoint(vector<float> planeVector, vector<float>  planePoint, vector<float>  lineVector, vector<float>  linePoint)
{
	vector<float> res(3);
	float vpt = lineVector[0] * planeVector[0] + lineVector[1] * planeVector[1] + lineVector[2] * planeVector[2];
	if (abs(vpt) < E) return res;
	else
	{
		float t = ((planePoint[0] - linePoint[0]) * planeVector[0] + (planePoint[1] - linePoint[1]) * planeVector[1] + (planePoint[2] - linePoint[2]) * planeVector[2]) / vpt;

		res[0] = linePoint[0] + lineVector[0] * t;
		res[1] = linePoint[1] + lineVector[1] * t;
		res[2] = linePoint[2] + lineVector[2] * t;
	}
	return res;
}

vector<int> Viewer::findFace(MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2)
{
	vector<int> res(2, -1);
	int vh1NeighborCount = 0;
	for (MyMesh::VertexFaceIter vf_it1 = mesh.vf_iter(vh1); vf_it1.is_valid(); vf_it1++)
	{
		bool isface = true;
		vh1NeighborCount++;
		for (MyMesh::VertexFaceIter vf_it2 = mesh.vf_iter(vh2); vf_it2.is_valid(); vf_it2++)
		{
			if (*vf_it1 == *vf_it2)
			{
				isface = false;
				continue;
			}
		}
		if (isface)
		{
			res[0] = (*vf_it1).idx();
			//break;
		}
	}
	int vh2NerghborCount = 0;
	for (MyMesh::VertexFaceIter vf_it2 = mesh.vf_iter(vh2); vf_it2.is_valid(); vf_it2++)
	{
		bool isface = true;
		vh2NerghborCount++;
		for (MyMesh::VertexFaceIter vf_it1 = mesh.vf_iter(vh1); vf_it1.is_valid(); vf_it1++)
		{
			if (*vf_it1 == *vf_it2)
			{
				isface = false;
				continue;
			}
		}
		if (isface)
		{
			res[1] = (*vf_it2).idx();
			//break;
		}
	}

	if ((vh1NeighborCount == 1 && vh2NerghborCount == 1)) editmode_ = FREE;
	else if ((vh1NeighborCount == 1 && vh2NerghborCount == 2) || (vh1NeighborCount == 2 && vh2NerghborCount == 1)) editmode_ = ONEFACE;
	else editmode_ = TWOFACE;

	return res;
}
void Viewer::draw()
{

	glDisable(GL_LIGHTING);
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glClear(GL_STENCIL_BUFFER_BIT);

	GLUtesselator* tess = tesser();
	if (!tess) return;
	//mesh.request_face_colors();
	for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{

		//glColor3f(colorSpace[colorIdx % colorNum][0], colorSpace[colorIdx % colorNum][1], colorSpace[colorIdx % colorNum][2]);
		//colorIdx++;
		float r = mesh.color(*f_it)[0];
		//glColor3f(mesh.color(*f_it)[0] / 255., mesh.color(*f_it)[1] / 255., mesh.color(*f_it)[2] / 255.);
		glColor3f(colorSpace_[faceColorIdx[(*f_it).idx()] % colorNum_][0], colorSpace_[faceColorIdx[(*f_it).idx()] % colorNum_][1], colorSpace_[faceColorIdx[(*f_it).idx()] % colorNum_][2]);

		gluTessBeginPolygon(tess, NULL);
		gluTessBeginContour(tess);

		for (int i = 0; i<fvSize[(*f_it).idx()]; i++)
			gluTessVertex(tess, &v_data[(*f_it).idx()][i][0], &v_data[(*f_it).idx()][i][0]);


		gluTessEndContour(tess);
		gluTessEndPolygon(tess);
	}

	int n = mesh.n_edges();
	for (MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
	{
		//drawEdge((*e_it).idx());
		MyMesh::VertexHandle vh1(edgeVertex0[(*e_it).idx()]);
		MyMesh::VertexHandle vh2(edgeVertex1[(*e_it).idx()]);
		glLineWidth(2);
		glBegin(GL_LINE_STRIP);
		glColor3f(0.f, 0.f, 0.f);
		glVertex3f(mesh.point(vh1)[0], mesh.point(vh1)[1], mesh.point(vh1)[2]);
		glVertex3f(mesh.point(vh2)[0], mesh.point(vh2)[1], mesh.point(vh2)[2]);
		glEnd();
	}

	if (showMergeVertices)
	{
		glColor3f(1.f, 0.f, 0.f);
		for (int i = 0; i < mergeVertices.size(); i++)
		{
			drawVertex(mergeVertices[i]);
		}
		glColor3f(0.f, 0.f, 1.f);
		for (int i = 0; i < symVertices.size(); i++)
		{
			drawVertex(symVertices[i]);
		}
	}

	if (isCandidate)
	{
		if (!mergeVertices.empty())
		{
			glColor3f(1.f, 0.f, 0.f);
			for (int i = 0; i < mergeVertices.size() ; i++)
			{
				drawVertex(mergeVertices[i]);
			}
			//MyMesh::VertexHandle vh1(mergeVertices[0]);
			//for (int i = 0; i < mergeVertices.size(); i++)
			//{
			//	drawVertex(mergeVertices[i]);
			//	glColor3f(0.f, 1.f, 0.f);
			//	MyMesh::VertexHandle vh2(mergeVertices[i]);
			//	glLineStipple(2, 0x5555);
			//	glLineWidth(3);
			//	glEnable(GL_LINE_STIPPLE);
			//	glBegin(GL_LINE_STRIP);
			//	glVertex3f(mesh.point(vh1)[0], mesh.point(vh1)[1], mesh.point(vh1)[2]);
			//	glVertex3f(mesh.point(vh2)[0], mesh.point(vh2)[1], mesh.point(vh2)[2]);
			//	glEnd();
			//	glDisable(GL_LINE_STIPPLE);
			//}
		}
		if (!mergeFaces.empty())
		{
			glColor3f(1.f, 0.f, 0.f);
			for (int i = 0; i < mergeFaces.size(); i++)
			{
				MyMesh::FaceHandle fh(mergeFaces[i]);
				for (MyMesh::FaceEdgeIter fe_it = mesh.fe_iter(fh); fe_it.is_valid(); fe_it++)
				{
					drawEdge((*fe_it).idx());
				}
			}
		}
	}
	switch (type_)
	{
	case VERTEX:
	{
		glColor3f(1.f, 0.f, 0.f);
		//for (QList<int>::const_iterator it = selected_vertices.begin(), end = selected_vertices.end(); it != end; it++)
		//{
		//	drawVertex(selected_vertices.at(*it));
		//}
		for (int i = 0; i < selected_vertices.size(); i++)
		{
			drawVertex(selected_vertices.at(i));
		}
		glColor3f(0.f, 1.f, 0.f);
		for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
		{
			if (!selected_vertices.contains((*v_it).idx()))
				drawVertex((*v_it).idx());
		}
		if (selectionmode_ != NONE) drawSelectionRectangle();
		break;
	}
	case EDGE:
	{
		glColor3f(1.f, 0.f, 0.f);
		for (int i = 0; i < selected_edges.size(); i++)
		{
			drawEdge(selected_edges.at(i));
		}
		glColor3f(0.f, 1.f, 0.f);
		int n = mesh.n_edges();
		for (MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
		{
			if (!selected_edges.contains((*e_it).idx()))
			{
				drawEdge((*e_it).idx());
			}
		}
		if (selectionmode_ != NONE && selectionmode_ != EDIT) drawSelectionRectangle();
		else if (selectionmode_ == EDIT)
		{
			//drawSelectionRectangle();
			if (selected_edges.size() == 1)
			{
				//drawMoveLine();
				updateEdit();
				line_.setP1(line_.p2());
				updateData();
				emit updatelayout();
				//updateEdit();
			}
		}
		break;
	}
	case FACE:
		glColor4f(0.f, 0.f, 0.f, 1.f);

		for (int i = 0; i < selected_faces.size(); i++)
		{
			MyMesh::FaceHandle fh(selected_faces.at(i));
			glColor3f(1.f, 0.f, 0.f);
			for (MyMesh::FaceEdgeIter fe_it = mesh.fe_iter(fh); fe_it.is_valid(); fe_it++)
			{
				drawEdge((*fe_it).idx());
			}
			//mesh.set_color(fh, MyMesh::Color(0, 0, 0));
			//drawFace(selected_faces.at(i));
		}
		if (selectionmode_ != NONE) drawSelectionRectangle();
		break;
	default: break;
	}

}

void Viewer::animate()
{
	//angle 180->90
	int angleInterval = 5;
	int n_needtoFold = flat_mesh.n_faces() - 1;
	if(isAnimation)
	{
		StereoMesh tmp_m = StereoMesh(flat_mesh);
		if (curAngle <= 90 && hasvisited >= n_needtoFold)
		{
			curAngle = 180;
			hasvisited = 1;
		}
		else if (curAngle > 90 && hasvisited <= n_needtoFold)
			curAngle = curAngle - angleInterval >= 90 ? curAngle - angleInterval : 90;
		else if (curAngle = 90 && hasvisited < n_needtoFold)
		{
			curAngle = 180;
			hasvisited++;
		}
		else
		//还有什么情况没有考虑到么？
		{
			curAngle = 180;
			hasvisited = 1;
		}
		qDebug() << hasvisited << "  " << curAngle << endl;
		mesh = tmp_m.folding(hasvisited, curAngle);
		updateData();
	}
}

vector<float> Viewer::calFaceNormal(MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2)
{
	vector<int> vidx;
	for (MyMesh::VertexVertexIter vv_it = mesh.vv_iter(vh1); vv_it.is_valid(); vv_it++)
	{
		if ((*vv_it) != vh2)
		{
			vidx.push_back((*vv_it).idx());
		}
	}
	MyMesh::VertexHandle v0(vidx[0]);
	MyMesh::VertexHandle v1(vidx[1]);
	float x1 = mesh.point(vh1)[0] - mesh.point(v0)[0];
	float y1 = mesh.point(vh1)[1] - mesh.point(v0)[1];
	float z1 = mesh.point(vh1)[2] - mesh.point(v0)[2];

	float x2 = mesh.point(vh1)[0] - mesh.point(v1)[0];
	float y2 = mesh.point(vh1)[1] - mesh.point(v1)[1];
	float z2 = mesh.point(vh1)[2] - mesh.point(v1)[2];

	vector<float> res(3);
	res[0] = y1 * z1 - z1 * y2;
	res[1] = x2 * z1 - x1 * z2;
	res[2] = x1 * y2 - y1 * x2;
	return res;
}

void Viewer::updateEdit()
{
	//qreal start[3], end[3];
	//transferPos(line_.x1(), line_.y1(), start, start + 1, start + 2);
 //	transferPos(line_.x2(), line_.y2(), end, end + 1, end + 2);
	//qreal dist[3] = { (end[0] - start[0]) / 100 , (end[1] - start[1]) /100 , (end[2] - start[2]) / 100 };

	////MyMesh::EdgeHandle eh(selected_edges.at(0));
	//MyMesh::VertexHandle vh1(edgeVertex0[selected_edges.at(0)]);
	//MyMesh::VertexHandle vh2(edgeVertex1[selected_edges.at(0)]);
	//vector<float> lineV(3);
	////float *lineV = new float(3);
	//lineV[0] = mesh.point(vh1)[0] - mesh.point(vh2)[0];
	//lineV[1] = mesh.point(vh1)[1] - mesh.point(vh2)[1];
	//lineV[2] = mesh.point(vh1)[2] - mesh.point(vh2)[2];
	//vector <float> linePoint(3);
	////float *linePoint = new float(3);

	////haha,gluUnProject error, test by *10
	//linePoint[0] = mesh.point(vh1)[0] + dist[0] ;
	//linePoint[1] = mesh.point(vh1)[1] + dist[1] ;
	//linePoint[2] = mesh.point(vh1)[2] + dist[2] ;

	//vector<int> faceidx;
	//faceidx = findFace(vh1, vh2);
	//MyMesh::FaceHandle fh1(faceidx[0]), fh2(faceidx[1]);
	////float tmp = mesh.normal(fh1)[0];
	//vector<float> planeV1(3);
	////float *planeV1 = new float(3);
	//planeV1[0] = mesh.normal(fh1)[0];
	//planeV1[1] = mesh.normal(fh1)[1];
	//planeV1[2] = mesh.normal(fh1)[2];
	//vector<float> planeP1(3);
	////float *planeP1 = new float(3);
	//planeP1[0] = mesh.point(vh1)[0];
	//planeP1[1] = mesh.point(vh1)[1];
	//planeP1[2] = mesh.point(vh1)[2];
	//vector<float> planeV2(3);
	////float *planeV2 = new float(3);
	//planeV2[0] = mesh.normal(fh2)[0];
	//planeV2[1] = mesh.normal(fh2)[1];
	//planeV2[2] = mesh.normal(fh2)[2];
	//vector<float> planeP2(3);
	////float *planeP2 = new float(3);
	//planeP2[0] = mesh.point(vh2)[0];
	//planeP2[1] = mesh.point(vh2)[1];
	//planeP2[2] = mesh.point(vh2)[2];

	//vector<float> p1, p2;
	//p1 = CalPoint(planeV1, planeP1, lineV, linePoint);
	//p2 = CalPoint(planeV2, planeP2, lineV, linePoint);

	//
	//mesh.point(vh1)[0] = p1[0];
	//mesh.point(vh1)[1] = p1[1];
	//mesh.point(vh1)[2] = p1[2];
	//mesh.point(vh2)[0] = p2[0];
	//mesh.point(vh2)[1] = p2[1];
	//mesh.point(vh2)[2] = p2[2];
	qreal start[3], end[3];
	transferPos(line_.x1(), line_.y1(), start, start + 1, start + 2);
	transferPos(line_.x2(), line_.y2(), end, end + 1, end + 2);
	//qreal dist[3] = { end[0] - start[0] , end[1] - start[1] , end[2] - start[2] };
	qreal dist[3] = { (end[0] - start[0])  , (end[1] - start[1])  , (end[2] - start[2]) };

	//MyMesh::EdgeHandle eh(selected_edges.at (0));
	MyMesh::VertexHandle vh1(edgeVertex0[selected_edges.at(0)]);
	MyMesh::VertexHandle vh2(edgeVertex1[selected_edges.at(0)]);
	vector<float> lineV(3);
	//float *lineV = new float(3);
	lineV[0] = mesh.point(vh1)[0] - mesh.point(vh2)[0];
	lineV[1] = mesh.point(vh1)[1] - mesh.point(vh2)[1];
	lineV[2] = mesh.point(vh1)[2] - mesh.point(vh2)[2];
	vector <float> linePoint(3);
	//float *linePoint = new float(3);
	linePoint[0] = mesh.point(vh1)[0] + dist[0];
	linePoint[1] = mesh.point(vh1)[1] + dist[1];
	linePoint[2] = mesh.point(vh1)[2] + dist[2];

	vector<int> faceidx;

	faceidx = findFace(vh1, vh2);
	if (editmode_ == FREE)
	{
		MyMesh::VertexFaceIter vf_it = mesh.vf_iter(vh1);
		vector<float> direction(3);
		float x1 = mesh.normal(*vf_it)[0];
		float y1 = mesh.normal(*vf_it)[1];
		float z1 = mesh.normal(*vf_it)[2];

		float x2 = mesh.point(vh1)[0] - mesh.point(vh2)[0];
		float y2 = mesh.point(vh1)[1] - mesh.point(vh2)[1];
		float z2 = mesh.point(vh1)[2] - mesh.point(vh2)[2];

		direction[0] = y1 * z1 - z1 * y2;
		direction[1] = x2 * z1 - x1 * z2;
		direction[2] = x1 * y2 - y1 * x2;

		float newdist = direction[0] * dist[0] + direction[1] + dist[1] + direction[2] * dist[2] ;
		float tmp = direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2];

		vector<float> newd(3);
		newd[0] = newdist * direction[0] / tmp;
		newd[1] = newdist * direction[1] / tmp;
		newd[2] = newdist * direction[2] / tmp;

		mesh.point(vh1)[0] = mesh.point(vh1)[0] + newd[0];
		mesh.point(vh1)[1] = mesh.point(vh1)[1] + newd[1];
		mesh.point(vh1)[2] = mesh.point(vh1)[2] + newd[2];

		mesh.point(vh2)[0] = mesh.point(vh2)[0] + newd[0];
		mesh.point(vh2)[1] = mesh.point(vh2)[1] + newd[1];
		mesh.point(vh2)[2] = mesh.point(vh2)[2] + newd[2];
	}
	else if (editmode_ == ONEFACE)
	{
		int fidx;
		MyMesh::VertexHandle vh, vh_;
		int flag;
		if (faceidx[0] == -1)
		{
			fidx = faceidx[1];
			vh = vh2;
			vh_ = vh1;
			flag = 1;
		}
		else
		{
			fidx = faceidx[0];
			vh = vh1;
			vh_ = vh2;
			flag = -1;
		}
		MyMesh::FaceHandle fh(fidx);
		vector<float> planeV(3);
		planeV[0] = mesh.normal(fh)[0];
		planeV[1] = mesh.normal(fh)[1];
		planeV[2] = mesh.normal(fh)[2];

		vector<float> planeP(3);
		planeP[0] = mesh.point(vh)[0];
		planeP[1] = mesh.point(vh)[1];
		planeP[2] = mesh.point(vh)[2];

		vector<float> p;
		p = CalPoint(planeV, planeP, lineV, linePoint);
		mesh.point(vh)[0] = p[0];
		mesh.point(vh)[1] = p[1];
		mesh.point(vh)[2] = p[2];

		mesh.point(vh_)[0] = p[0] + flag * lineV[0];
		mesh.point(vh_)[1] = p[1] + flag * lineV[1];
		mesh.point(vh_)[2] = p[2] + flag * lineV[2];
	}
	else if (editmode_ == TWOFACE)
	{
		vector<float> planeV1(3);
		vector<float> planeP1(3);
		vector<float> planeV2(3);
		vector<float> planeP2(3);

		planeP2[0] = mesh.point(vh2)[0];
		planeP2[1] = mesh.point(vh2)[1];
		planeP2[2] = mesh.point(vh2)[2];

		planeP1[0] = mesh.point(vh1)[0];
		planeP1[1] = mesh.point(vh1)[1];
		planeP1[2] = mesh.point(vh1)[2];

		if (faceidx[0] != -1 && faceidx[1] != -1)
		{
			MyMesh::FaceHandle fh1(faceidx[0]), fh2(faceidx[1]);
			//float tmp = mesh.normal(fh1)[0];

			//float *planeV1 = new float(3);
			planeV1[0] = mesh.normal(fh1)[0];
			planeV1[1] = mesh.normal(fh1)[1];
			planeV1[2] = mesh.normal(fh1)[2];


			planeV2[0] = mesh.normal(fh2)[0];
			planeV2[1] = mesh.normal(fh2)[1];
			planeV2[2] = mesh.normal(fh2)[2];

		}
		else
		{
			planeV1 = calFaceNormal(vh1, vh2);
			planeV2 = calFaceNormal(vh2, vh1);

		}
		vector<float> p1, p2;
		p1 = CalPoint(planeV1, planeP1, lineV, linePoint);
		p2 = CalPoint(planeV2, planeP2, lineV, linePoint);


		mesh.point(vh1)[0] = p1[0];
		mesh.point(vh1)[1] = p1[1];
		mesh.point(vh1)[2] = p1[2];
		mesh.point(vh2)[0] = p2[0];
		mesh.point(vh2)[1] = p2[1];
		mesh.point(vh2)[2] = p2[2];
	}
	
	//joint move
	//In order to let all merging vertexes move with given vertex
	//change the vertex coordinate 
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		if ( (*v_it)!=vh1 && m_label[(*v_it).idx()] ==m_label[vh1.idx()])
			mesh.set_point(*v_it, mesh.point(vh1));

		else if ((*v_it) != vh2 && m_label[(*v_it).idx()] == m_label[vh2.idx()])
			mesh.set_point(*v_it, mesh.point(vh2));

	}

}

void Viewer::transferPos(int x, int y, qreal* pos0, qreal* pos1, qreal* pos2)
{
	bool found;
	qglviewer::Vec pos = this->camera()->pointUnderPixel( QPoint(x, y), found);
	*pos0 = pos[0];
	*pos1 = pos[1];
	*pos2 = pos[2];
	//GLint viewport[4];
	//GLdouble modelview[16];
	//GLdouble projection[16];
	//GLfloat winX, winY, winZ;

	//winX = x;

	//glPushMatrix();

	//glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	//glGetDoublev(GL_PROJECTION_MATRIX, projection);
	//glGetIntegerv(GL_VIEWPORT, viewport);

	////winY = this->height() - y;
	//winY = viewport[3] - y;
	//glPopMatrix();


	//glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
	//gluUnProject(winX, winY, winZ, modelview, projection, viewport, pos0, pos1, pos2);

}

void Viewer::mousePressEvent(QMouseEvent *e)
{
	if (type_ != NONE_)
	{
		rectangle_ = QRect(e->pos(), e->pos());
		line_ = QLine(e->pos(), e->pos());
		if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::ShiftModifier))
		{
			selectionmode_ = ADD;
		}
		else
			if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::AltModifier))
			{
				selectionmode_ = REMOVE;
			}
			else if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::ControlModifier))
			{
				selectionmode_ = EDIT;
			}
			else
			{
				QGLViewer::mousePressEvent(e);
			}
	}
	else QGLViewer::mousePressEvent(e);
}

void Viewer::mouseMoveEvent(QMouseEvent *e)
{
	if ((type_ != NONE_) && (selectionmode_ != NONE) && (selectionmode_ != EDIT))
	{
		rectangle_.setBottomRight(e->pos());
		update();
	}
	else if ((type_ != NONE_) && (selectionmode_ == EDIT))
	{
		//rectangle_.setBottomRight(e->pos());
		line_.setP2(e->pos());
		//若偏移太大，会进行限制
		/*int restrict_translate = 10;
		if (line_.dx() > restrict_translate)
		{
			QPoint a(line_.x1() + restrict_translate, e->pos().y);
			line_.setP2(a);
		}			
		if (line_.dy() > restrict_translate)
		{
			QPoint b(e->pos().x, line_.y1() + restrict_translate);
			line_.setP2(b);
		}*/
			
		update();
	}
	else
		QGLViewer::mouseMoveEvent(e);
}

void Viewer::mouseReleaseEvent(QMouseEvent *e)
{
	if (selectionmode_ != NONE && selectionmode_ != EDIT)
	{
		rectangle_ = rectangle_.normalized();
		setSelectRegionWidth(rectangle_.width());
		setSelectRegionHeight(rectangle_.height());
		select(rectangle_.center());
		update();
	}
	else if (selectionmode_ == EDIT)
	{
		//updateEdit();
		//updateData();
		update();
		selectionmode_ = NONE;

		emit computeLengthError();
	}
	else QGLViewer::mouseReleaseEvent(e);
}

void Viewer::mouseDoubleClickEvent(QMouseEvent *e)
{
	if (isCandidate)
	{
		emit isDoubleClicked(viewIdx);
	}
	else QGLViewer::mouseDoubleClickEvent(e);
}

void Viewer::drawWithNames()
{
	switch (type_)
	{
	case VERTEX:
		drawVertexWithName();
		break;
	case FACE:
		drawFaceWithName();
		break;
	case EDGE:
		drawEdgeWithName();
		break;
	default:
		break;
	}
}

void Viewer::drawVertexWithName()
{
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		glPushName(i);
		drawVertex(i);
		glPopName();
	}
}
void Viewer::drawFaceWithName()
{
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		glPushName(i);
		drawFace(i);
		glPopName();
	}
}

void Viewer::drawEdgeWithName()
{
	for (int i = 0; i < mesh.n_edges(); i++)
	{
		glPushName(i);
		drawEdge(i);
		glPopName();
	}
}

void Viewer::endSelection(const QPoint &)
{
	glFlush();
	GLint nbHits = glRenderMode(GL_RENDER);
	if (nbHits > 0)
	{
		for (int i = 0; i < nbHits; i++)
			switch (selectionmode_)
			{
			case ADD:
				addIdToSelection((selectBuffer())[4 * i + 3]);
				break;
			case REMOVE:
				removeIdFromSelection((selectBuffer())[4 * i + 3]);
				break;
			default:
				break;
			}
	}
	selectionmode_ = NONE;
}

void Viewer::addIdToSelection(int id)
{
	switch (type_)
	{
	case VERTEX:
	{
		if (!selected_vertices.contains(id))
			selected_vertices.push_back(id);
		break;
	}
	case EDGE:
	{
		if (!selected_edges.contains(id))
			selected_edges.push_back(id);
		break;
	}
	case FACE:
	{
		if (!selected_faces.contains(id))
			selected_faces.push_back(id);
		break;
	}
	default:
		break;
	}

}

void Viewer::removeIdFromSelection(int id)
{
	switch (type_)
	{
	case Viewer::NONE_:
		break;
	case Viewer::VERTEX:
		selected_vertices.removeAll(id);
		break;
	case Viewer::FACE:
	{
		selected_faces.removeAll(id);
		//MyMesh::FaceHandle fh(id);
		//mesh.set_color(fh, MyMesh::Color(colorSpace_[id % colorNum_][0] * 255, colorSpace_[id % colorNum_][1] * 255, colorSpace_[id % colorNum_][2] * 255));
		break;
	}
	case Viewer::EDGE:
		selected_edges.removeAll(id);
		break;
	default:
		break;
	}

}

void Viewer::drawMoveLine() const
{

	startScreenCoordinatesSystem();
	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);

	glLineWidth(2.0);
	glColor4f(0.4f, 0.4f, 0.5f, 0.5f);
	glBegin(GL_LINES);
	glVertex2i(line_.x1(), line_.y1());
	glVertex2i(line_.x2(), line_.y2());

	glEnd();

	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	stopScreenCoordinatesSystem();
}

void Viewer::drawSelectionRectangle() const
{
	startScreenCoordinatesSystem();
	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);

	glColor4f(0.0, 0.0, 0.3f, 0.3f);
	glBegin(GL_QUADS);
	glVertex2i(rectangle_.left(), rectangle_.top());
	glVertex2i(rectangle_.right(), rectangle_.top());
	glVertex2i(rectangle_.right(), rectangle_.bottom());
	glVertex2i(rectangle_.left(), rectangle_.bottom());
	glEnd();

	glLineWidth(2.0);
	glColor4f(0.4f, 0.4f, 0.5f, 0.5f);
	glBegin(GL_LINE_LOOP);
	glVertex2i(rectangle_.left(), rectangle_.top());
	glVertex2i(rectangle_.right(), rectangle_.top());
	glVertex2i(rectangle_.right(), rectangle_.bottom());
	glVertex2i(rectangle_.left(), rectangle_.bottom());
	glEnd();

	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	stopScreenCoordinatesSystem();
}

void Viewer::init()
{
	// Restore previous viewer state.
	setBackgroundColor(Qt::white);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glEnable(GL_BLEND);
	//glBlendFunc(GL_ONE, GL_ONE);
	//glEnable(GL_LIGHTING);
	//glEnable(GL_POINT_SMOOTH);
	//glEnable(GL_LINE_SMOOTH);
	//glEnable(GL_POLYGON_SMOOTH);
	//glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	//glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	//glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	camera()->setFieldOfView((15. / 180.)*3.1415926);
	restoreStateFromFile();
	//  help();

	hasvisited = 1;
	curAngle = 180;

	startAnimation();
}

