#include "StereoMesh.h"
#include<queue>
#include<qDebug>
#include <gl\GLU.h>

StereoMesh::StereoMesh()
{
	
}

StereoMesh::StereoMesh(MyMesh sm)
{
	m_sm = sm;
	m_tm = sm;


	// 计算出每条边的相邻两个面的id
	for (int i = 0; i < m_sm.n_edges(); i++)
	{
		m_f1_idx.push_back(-1);
		m_f2_idx.push_back(-1);
	}

	for (MyMesh::FaceIter f_it = m_sm.faces_begin(); f_it != m_sm.faces_end(); ++f_it)
	{
		m_fhlist.push_back(*f_it);
		for (MyMesh::FaceEdgeIter e_it = m_sm.fe_begin(*f_it); e_it != m_sm.fe_end(*f_it); ++e_it)
		{
			if (-1 == m_f1_idx[(*e_it).idx()])
				m_f1_idx[(*e_it).idx()] = (*f_it).idx();
			else if (((*f_it).idx() > m_f1_idx[(*e_it).idx()]) && (-1 == m_f2_idx[(*e_it).idx()]))
				m_f2_idx[(*e_it).idx()] = (*f_it).idx();
			else if (((*f_it).idx() < m_f1_idx[(*e_it).idx()]) && (-1 == m_f2_idx[(*e_it).idx()]))
			{
				m_f2_idx[(*e_it).idx()] = m_f1_idx[(*e_it).idx()];
				m_f1_idx[(*e_it).idx()] = (*f_it).idx();
			}
			else
				std::cerr << "More than two faces for one edge ? some thing is wrong!" << endl;
		}
	}

	//计算每条边的起点和终点id
	for (int i = 0; i < m_sm.n_edges(); i++)
	{
		m_v1_idx.push_back(-1);
		m_v2_idx.push_back(-1);
	}

	for (MyMesh::VertexIter v_it = m_sm.vertices_begin(); v_it != m_sm.vertices_end(); ++v_it)
	{
		m_vhlist.push_back(*v_it);
		for (MyMesh::VertexEdgeIter e_it = m_sm.ve_begin(*v_it); e_it != m_sm.ve_end(*v_it); ++e_it)
		{
			if (-1 == m_v1_idx[(*e_it).idx()])
				m_v1_idx[(*e_it).idx()] = (*v_it).idx();
			else if (((*v_it).idx() > m_v1_idx[(*e_it).idx()]) && (-1 == m_v2_idx[(*e_it).idx()]))
				m_v2_idx[(*e_it).idx()] = (*v_it).idx();
			else if (((*v_it).idx() < m_v1_idx[(*e_it).idx()]) && (-1 == m_v2_idx[(*e_it).idx()]))
			{
				m_v2_idx[(*e_it).idx()] = m_v1_idx[(*e_it).idx()];
				m_v1_idx[(*e_it).idx()] = (*v_it).idx();
			}
			else
				std::cerr << "More than two vertices for one edge ? some thing is wrong!" << endl;
		}
	}

	//创建旋转矩阵
	m_RTMatrix = new float **[m_sm.n_faces()];
	for (int i = 0; i < m_sm.n_faces(); i++)
	{
		m_RTMatrix[i] = new float *[4];
		for (int j = 0; j < 4; j++)
			m_RTMatrix[i][j] = new float[4];
	}

	//初始化旋转矩阵
	for (int i = 0; i < m_sm.n_faces(); i++)
		for (int j = 0; j < 4; j++)
			for (int k = 0; k < 4; k++)
			{
				if ((j == 0 && k == 0)
					|| (j == 1 && k == 1)
					|| (j == 2 && k == 2)
					|| (j == 3 && k == 3))
					m_RTMatrix[i][j][k] = 1;
				else
					m_RTMatrix[i][j][k] = 0;
			}

	m_initialAngle = 90;
	
}


StereoMesh::~StereoMesh()
{
	for (int i = 0; i < m_sm.n_faces(); i++)
	{
		for (int j = 0; j < 4; j++)
			delete m_RTMatrix[i][j];
		delete m_RTMatrix[i];
	}

	delete m_RTMatrix;
}

float StereoMesh::computeRotateAngle(float angle, MyMesh::EdgeHandle edge, int f1, int f2)
{
	m_sm.request_face_normals();
	MyMesh::Normal n1 = m_sm.normal(m_fhlist[f1]);
	MyMesh::Normal n2 = m_sm.normal(m_fhlist[f2]);

	OpenMesh::Vec3f cross = OpenMesh::cross(n1, n2);

	MyMesh::VertexHandle v1 = m_vhlist[m_v1_idx[edge.idx()]];
	MyMesh::VertexHandle v2 = m_vhlist[m_v2_idx[edge.idx()]];

	float a = m_sm.point(v1)[0];
	float b = m_sm.point(v1)[1];
	float c = m_sm.point(v1)[2];

	float u = m_sm.point(v2)[0] - a;
	float v = m_sm.point(v2)[1] - b;
	float w = m_sm.point(v2)[2] - c;

	float nor = sqrtf(u*u + v*v + w*w);

	if (nor != 0)
	{
		u = u / nor;
		v = v / nor;
		w = w / nor;
	}

	float temp = u*cross[0] + v*cross[1] + w*cross[2];

	if (temp < -0.1)
		angle = angle - 180;
	else
		angle = 180 - angle;
	return angle;
}

void StereoMesh::RotateArbitraryAxis(float angle, MyMesh::EdgeHandle edge, float rtmatrix[4][4], int start_p)
{
	// when end_plane rotate around start_plane, the axis's direction will be the same as the edge within start_plane
	// use halfedge data structure
	m_tm.request_face_normals();
	MyMesh::VertexHandle v1 = m_vhlist[m_v1_idx[edge.idx()]];
	MyMesh::VertexHandle v2 = m_vhlist[m_v2_idx[edge.idx()]];
	for (MyMesh::FaceHalfedgeIter he_it = m_sm.fh_begin(m_fhlist[start_p]); he_it != m_sm.fh_end(m_fhlist[start_p]); ++he_it)
	{
		auto from_v = m_sm.from_vertex_handle(*he_it);
		auto to_v = m_sm.to_vertex_handle(*he_it);
		if (v1 == from_v && v2 == to_v)
			break;
		else if (v1 == to_v && v2 == from_v)
		{
			auto temp = v1;
			v1 = v2;
			v2 = temp;
			break;
		}
	}




	float a = m_sm.point(v1)[0];
	float b = m_sm.point(v1)[1];
	float c = m_sm.point(v1)[2];

	float u = m_sm.point(v2)[0] - a;
	float v = m_sm.point(v2)[1] - b;
	float w = m_sm.point(v2)[2] - c;

	float nor = sqrtf(u*u + v*v + w*w);

	if (nor != 0)
	{
		u = u / nor;
		v = v / nor;
		w = w / nor;
	}

	float uu = u * u;
	float uv = u * v;
	float uw = u * w;
	float vv = v * v;
	float vw = v * w;
	float ww = w * w;
	float au = a * u;
	float av = a * v;
	float aw = a * w;
	float bu = b * u;
	float bv = b * v;
	float bw = b * w;
	float cu = c * u;
	float cv = c * v;
	float cw = c * w;

	float costheta = cosf(angle*g_parameters.MyPi / 180);
	float sintheta = sinf(angle*g_parameters.MyPi / 180);

	rtmatrix[0][0] = uu + (vv + ww) * costheta;
	rtmatrix[1][0] = uv * (1 - costheta) + w * sintheta;
	rtmatrix[2][0] = uw * (1 - costheta) - v * sintheta;
	rtmatrix[3][0] = 0;

	rtmatrix[0][1] = uv * (1 - costheta) - w * sintheta;
	rtmatrix[1][1] = vv + (uu + ww) * costheta;
	rtmatrix[2][1] = vw * (1 - costheta) + u * sintheta;
	rtmatrix[3][1] = 0;

	rtmatrix[0][2] = uw * (1 - costheta) + v * sintheta;
	rtmatrix[1][2] = vw * (1 - costheta) - u * sintheta;
	rtmatrix[2][2] = ww + (uu + vv) * costheta;
	rtmatrix[3][2] = 0;

	rtmatrix[0][3] = (a * (vv + ww) - u * (bv + cw)) * (1 - costheta) + (bw - cv) * sintheta;
	rtmatrix[1][3] = (b * (uu + ww) - v * (au + cw)) * (1 - costheta) + (cu - aw) * sintheta;
	rtmatrix[2][3] = (c * (uu + vv) - w * (au + bv)) * (1 - costheta) + (av - bu) * sintheta;
	rtmatrix[3][3] = 1;
}

void StereoMesh::changeVertexCoordinate(float ** rtmatrix, MyMesh::VertexHandle vh)
{
	float x = rtmatrix[0][0] * m_sm.point(vh)[0] + rtmatrix[0][1] * m_sm.point(vh)[1] + rtmatrix[0][2] * m_sm.point(vh)[2] + rtmatrix[0][3];
	float y = rtmatrix[1][0] * m_sm.point(vh)[0] + rtmatrix[1][1] * m_sm.point(vh)[1] + rtmatrix[1][2] * m_sm.point(vh)[2] + rtmatrix[1][3];
	float z = rtmatrix[2][0] * m_sm.point(vh)[0] + rtmatrix[2][1] * m_sm.point(vh)[1] + rtmatrix[2][2] * m_sm.point(vh)[2] + rtmatrix[2][3];
	m_tm.set_point(vh, MyMesh::Point(x, y, z));
}

void StereoMesh::changeVertexCoordinate(float **rtmatrix, int f_id)
{
	for (MyMesh::FaceVertexIter v_it = m_sm.fv_begin(m_fhlist[f_id]); v_it != m_sm.fv_end(m_fhlist[f_id]); ++v_it)
	{
		float x = rtmatrix[0][0] * m_sm.point(*v_it)[0] + rtmatrix[0][1] * m_sm.point(*v_it)[1] + rtmatrix[0][2] * m_sm.point(*v_it)[2] + rtmatrix[0][3];
		float y = rtmatrix[1][0] * m_sm.point(*v_it)[0] + rtmatrix[1][1] * m_sm.point(*v_it)[1] + rtmatrix[1][2] * m_sm.point(*v_it)[2] + rtmatrix[1][3];
		float z = rtmatrix[2][0] * m_sm.point(*v_it)[0] + rtmatrix[2][1] * m_sm.point(*v_it)[1] + rtmatrix[2][2] * m_sm.point(*v_it)[2] + rtmatrix[2][3];
		m_tm.set_point(*v_it, MyMesh::Point(x, y, z));
	}
}

double StereoMesh::calFaceArea(MyMesh sm, MyMesh::FaceHandle fh)
{
	//要求Z坐标均为0
	vector<double> fvx;
	vector<double> fvy;
	for (MyMesh::FaceVertexIter fv_it = sm.fv_begin(fh); fv_it != sm.fv_end(fh); fv_it++)
	{
		if (sm.point(*fv_it)[2] != 0)
			qDebug() << "Z coordinate is not zero!!!" << endl;
		else
		{
			fvx.push_back(sm.point(*fv_it)[0]);
			fvy.push_back(sm.point(*fv_it)[1]);
		}
	}

	double s = 0;
	for (int i = 0; i < fvx.size(); i++)
		s += fvx[i] * fvy[(i + 1) % fvx.size()] - fvy[i] * fvx[(i + 1) % fvx.size()];

	return fabs(s / 2);
}

int StereoMesh::MaxFaceIdx(MyMesh sm)
{
	int res;
	double max_s = 0;
	for (MyMesh::FaceIter f_it = sm.faces_begin(); f_it != sm.faces_end(); f_it++)
	{
		double tmp = calFaceArea(sm, *f_it);
		if (tmp > max_s)
		{
			max_s = tmp;
			res = (*f_it).idx();
		}
	}
	return res;
}

MyMesh StereoMesh::folding()
{
	vector<int> visited(m_sm.n_faces(), 0);
	int n_visited = 0; // sum of visited number

	//push the biggest face into queue
	int max_face = MaxFaceIdx(m_sm);
	visited[max_face] = 1;


	//bfs
	//visited[0] = 1;
	n_visited += 1;

	queue<int> my_q;
	int f_id;
	float t_angle;
	float rtmatrix[4][4], t_rtm[4][4];
	for (MyMesh::FaceEdgeIter e_it = m_sm.fe_begin(m_fhlist[max_face]); e_it != m_sm.fe_end(m_fhlist[max_face]); ++e_it)
	{
		if (m_fhlist[max_face].idx() == m_f1_idx[(*e_it).idx()])
			f_id = m_f2_idx[(*e_it).idx()];
		else
			f_id = m_f1_idx[(*e_it).idx()];

		if (f_id != -1 && visited[f_id] == 0)
		{
			my_q.push(f_id);
			t_angle = computeRotateAngle(m_initialAngle, *e_it, f_id, m_fhlist[max_face].idx());

			RotateArbitraryAxis(t_angle, *e_it, rtmatrix, m_fhlist[max_face].idx());

			//这个全局矩阵左乘的顺序要注意！
			for (int p = 0; p < 4; p++)
				for (int q = 0; q < 4; q++)
				{
					t_rtm[p][q] = 0;
					for (int r = 0; r < 4; r++)
						t_rtm[p][q] += m_RTMatrix[m_fhlist[max_face].idx()][p][r] * rtmatrix[r][q];
				}

			for (int p = 0; p < 4; p++)
				for (int q = 0; q < 4; q++)
				{
					m_RTMatrix[f_id][p][q] = t_rtm[p][q];
				}

			changeVertexCoordinate(m_RTMatrix[f_id], f_id);
			n_visited += 1;
			visited[f_id] = 1;
		}
	}

	while (!my_q.empty())
	{
		int t_q = my_q.front();
		my_q.pop();
		for (MyMesh::FaceEdgeIter e_it = m_sm.fe_begin(m_fhlist[t_q]); e_it != m_sm.fe_end(m_fhlist[t_q]); ++e_it)
		{
			if (m_fhlist[t_q].idx() == m_f1_idx[(*e_it).idx()])
				f_id = m_f2_idx[(*e_it).idx()];
			else
				f_id = m_f1_idx[(*e_it).idx()];

			if (f_id != -1 && visited[f_id] == 0)
			{
				my_q.push(f_id);
				t_angle = computeRotateAngle(m_initialAngle, *e_it, f_id, m_fhlist[t_q].idx());

				RotateArbitraryAxis(t_angle, *e_it, rtmatrix, t_q);

				//这个全局矩阵左乘的顺序要注意！
				for (int p = 0; p < 4; p++)
					for (int q = 0; q < 4; q++)
					{
						t_rtm[p][q] = 0;
						for (int r = 0; r < 4; r++)
							t_rtm[p][q] += m_RTMatrix[m_fhlist[t_q].idx()][p][r] * rtmatrix[r][q];
					}

				for (int p = 0; p < 4; p++)
					for (int q = 0; q < 4; q++)
					{
						m_RTMatrix[f_id][p][q] = t_rtm[p][q];
					}
				changeVertexCoordinate(m_RTMatrix[f_id], f_id);
				n_visited += 1;
				visited[f_id] = 1;
			}

		}
	}

	m_tm.update_face_normals();

	return m_tm;
}

//for animation

unordered_set<MyMesh::VertexHandle> StereoMesh::vertexOnUnvisitedAdjface(int fid, vector<int> visited)
{
	unordered_set<MyMesh::VertexHandle> res;
	for (MyMesh::FaceVertexIter fv_it = m_sm.fv_begin(m_fhlist[fid]); fv_it != m_sm.fv_end(m_fhlist[fid]); ++fv_it)
		res.insert(*fv_it);

	queue<int> my_q;
	my_q.push(fid);

	int f_id;

	while (!my_q.empty())
	{
		int t_q = my_q.front();
		my_q.pop();
		for (MyMesh::FaceEdgeIter e_it = m_sm.fe_begin(m_fhlist[t_q]); e_it != m_sm.fe_end(m_fhlist[t_q]); ++e_it)
		{
			if (m_fhlist[t_q].idx() == m_f1_idx[(*e_it).idx()])
				f_id = m_f2_idx[(*e_it).idx()];
			else
				f_id = m_f1_idx[(*e_it).idx()];
			if (f_id != -1 && visited[f_id] == 0)
			{
				my_q.push(f_id);
				for (MyMesh::FaceVertexIter fv_it = m_sm.fv_begin(m_fhlist[f_id]); fv_it != m_sm.fv_end(m_fhlist[f_id]); ++fv_it)
				{
					if (res.find(*fv_it) == res.end())
						res.insert(*fv_it);
				}
				visited[f_id] = 1;
			}
		}
	}

	return res;
}


MyMesh StereoMesh::folding(int hasvisited, float finalAngle)
{
	vector<int> visited(m_sm.n_faces(), 0);
	int n_visited = 0; // sum of visited number

					   //push the biggest face into queue
	int max_face = MaxFaceIdx(m_sm);
	visited[max_face] = 1;


	//bfs
	//visited[0] = 1;
	n_visited += 1;

	if (n_visited <= hasvisited)
	{
		queue<int> my_q;
		int f_id;
		float t_angle;
		float rtmatrix[4][4], t_rtm[4][4];
		for (MyMesh::FaceEdgeIter e_it = m_sm.fe_begin(m_fhlist[max_face]); e_it != m_sm.fe_end(m_fhlist[max_face]); ++e_it)
		{
			if (n_visited <= hasvisited)
			{
				if (m_fhlist[max_face].idx() == m_f1_idx[(*e_it).idx()])
					f_id = m_f2_idx[(*e_it).idx()];
				else
					f_id = m_f1_idx[(*e_it).idx()];

				if (f_id != -1 && visited[f_id] == 0)
				{
					my_q.push(f_id);

					if (n_visited == hasvisited)
						t_angle = computeRotateAngle(finalAngle, *e_it, f_id, m_fhlist[max_face].idx());
					else
						t_angle = computeRotateAngle(m_initialAngle, *e_it, f_id, m_fhlist[max_face].idx());

					RotateArbitraryAxis(t_angle, *e_it, rtmatrix, m_fhlist[max_face].idx());

					//这个全局矩阵左乘的顺序要注意！
					for (int p = 0; p < 4; p++)
						for (int q = 0; q < 4; q++)
						{
							t_rtm[p][q] = 0;
							for (int r = 0; r < 4; r++)
								t_rtm[p][q] += m_RTMatrix[m_fhlist[max_face].idx()][p][r] * rtmatrix[r][q];
						}

					for (int p = 0; p < 4; p++)
						for (int q = 0; q < 4; q++)
						{
							m_RTMatrix[f_id][p][q] = t_rtm[p][q];
						}

					unordered_set<MyMesh::VertexHandle> tmp = vertexOnUnvisitedAdjface(f_id, visited);
					for (unordered_set<MyMesh::VertexHandle>::iterator it = tmp.begin(); it != tmp.end(); it++)
						changeVertexCoordinate(m_RTMatrix[f_id], *it);

					//changeVertexCoordinate(m_RTMatrix[f_id], f_id);
					n_visited += 1;
					visited[f_id] = 1;
				}
			}
			else
			{
				m_tm.update_face_normals();

				return m_tm;
			}
		}

		while (!my_q.empty() && n_visited <= hasvisited)
		{
			int t_q = my_q.front();
			my_q.pop();
			for (MyMesh::FaceEdgeIter e_it = m_sm.fe_begin(m_fhlist[t_q]); e_it != m_sm.fe_end(m_fhlist[t_q]); ++e_it)
			{
				if (n_visited <= hasvisited)
				{
					if (m_fhlist[t_q].idx() == m_f1_idx[(*e_it).idx()])
						f_id = m_f2_idx[(*e_it).idx()];
					else
						f_id = m_f1_idx[(*e_it).idx()];

					if (f_id != -1 && visited[f_id] == 0)
					{
						my_q.push(f_id);

						if (n_visited == hasvisited)
							t_angle = computeRotateAngle(finalAngle, *e_it, f_id, m_fhlist[max_face].idx());
						else
							t_angle = computeRotateAngle(m_initialAngle, *e_it, f_id, m_fhlist[max_face].idx());

						RotateArbitraryAxis(t_angle, *e_it, rtmatrix, t_q);

						//这个全局矩阵左乘的顺序要注意！
						for (int p = 0; p < 4; p++)
							for (int q = 0; q < 4; q++)
							{
								t_rtm[p][q] = 0;
								for (int r = 0; r < 4; r++)
									t_rtm[p][q] += m_RTMatrix[m_fhlist[t_q].idx()][p][r] * rtmatrix[r][q];
							}

						for (int p = 0; p < 4; p++)
							for (int q = 0; q < 4; q++)
							{
								m_RTMatrix[f_id][p][q] = t_rtm[p][q];
							}
						//changeVertexCoordinate(m_RTMatrix[f_id], f_id);

						unordered_set<MyMesh::VertexHandle> tmp = vertexOnUnvisitedAdjface(f_id, visited);
						for (unordered_set<MyMesh::VertexHandle>::iterator it = tmp.begin(); it != tmp.end(); it++)
							changeVertexCoordinate(m_RTMatrix[f_id], *it);

						n_visited += 1;
						visited[f_id] = 1;
					}
				}
				else
				{
					m_tm.update_face_normals();

					return m_tm;
				}
				

			}
		}
	}
	

	m_tm.update_face_normals();

	return m_tm;
}