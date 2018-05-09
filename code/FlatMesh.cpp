#include "FlatMesh.h"
#include<unordered_set>


FlatMesh::FlatMesh()
{
}


FlatMesh::~FlatMesh()
{
}

void FlatMesh::UpdateMaxMinValue(const float x, float& minx, float& maxx)
{
	minx = x < minx ? x : minx;
	maxx = x > maxx ? x : maxx;
}

void FlatMesh::ReadLinesFromTxt(string filepath)
{
	cerr << "Reading Line Segments From TXT..." << endl;
	LineSegment temp_l;
	Vertex temp_v;

	float t_value1 = 0;
	float t_value2 = 0;
	float t_value3 = 0;
	float t_value4 = 0;
	int id = 1;

	string oneLine;
	ifstream txt_in(filepath);
	string typeofLine;

	//需要知道二维坐标点的大概范围，以便使用Flood方法
	//查找二维点坐标中XY的最大最小值
	float t_max_x = 0;
	float t_max_y = 0;
	float t_min_x = 100000;
	float t_min_y = 100000;

	int isdash;

	//m_lineList.clear();


	//逐行读取
	while (txt_in)
	{
		txt_in >> typeofLine;
		if (typeofLine == "ISDASH")
			txt_in >> isdash;
		if (typeofLine == "LINESEGMENT" || typeofLine == "BEZIER")
		{
			txt_in >> t_value1 >> t_value2 >> t_value3 >> t_value4;

			//防止线段两个端点一样
			if (fabs(t_value1 - t_value3) < g_parameters.MyEpsilon && fabs(t_value2 - t_value4) < g_parameters.MyEpsilon)
				continue;

			//查找二维点坐标中XY的最大最小值
			UpdateMaxMinValue(t_value1, t_min_x, t_max_x);
			UpdateMaxMinValue(t_value3, t_min_x, t_max_x);
			UpdateMaxMinValue(t_value2, t_min_y, t_max_y);
			UpdateMaxMinValue(t_value4, t_min_y, t_max_y);

			/*
			if (t_value1 > t_value3)
			{
			if (t_max_x < t_value1)
			t_max_x = t_value1;
			if (t_min_x > t_value3)
			t_min_x = t_value3;
			}
			else
			{
			if (t_max_x < t_value3)
			t_max_x = t_value3;
			if (t_min_x > t_value1)
			t_min_x = t_value1;
			}

			if (t_value2 > t_value4)
			{
			if (t_max_y < t_value2)
			t_max_y = t_value2;
			if (t_min_y > t_value4)
			t_min_y = t_value4;
			}
			else
			{
			if (t_max_y < t_value4)
			t_max_y = t_value4;
			if (t_min_y > t_value2)
			t_min_y = t_value2;
			}
			*/

			temp_l.SetV1(t_value1, t_value2, 0);
			temp_l.SetV2(t_value3, t_value4, 0);
			temp_l.SetIsDash(isdash);
			m_lineList.push_back(temp_l);
			id++;
		}

	}
	txt_in.close();

	vector<LineSegment>::iterator it, it1;

	for (it = m_lineList.begin(); it != m_lineList.end(); )
	{
		//防止线段起始点相同
		if ((*it).GetV1() == (*it).GetV2())
		{
			it = m_lineList.erase(it);
		}
		else
			it++;
	}

	//防止线段重复读取
	for (it = m_lineList.begin(); it != m_lineList.end(); it++)
	{
		for (it1 = it + 1; it1 != m_lineList.end(); )
		{
			if (*it1 == *it)
			{
				if ((*it1).GetIsDash() != (*it).GetIsDash())
					(*it).SetIsDash(0);
				it1 = m_lineList.erase(it1);
			}
			else
				++it1;
		}
	}

	//测试List里是否存入正确数值
	/*for (int i = 0; i < 6; i++)
	{
	cout << _LineList[i].GetP1().GetX() << "    " << _LineList[i].GetP1().GetY() << "     "
	<< _LineList[i].GetP2().GetX() << "     " << _LineList[i].GetP2().GetY() << endl;
	}
	cout << _LineList.size() << endl;*/

	//切割线段
	CutOffLines();

	/*ofstream file1("./Lines_600.txt");
	for (int i = 0; i < _LineList.size(); i++)
	{
	file1 << _LineList[i].GetP1().GetX() << " " << _LineList[i].GetP1().GetY() << " " << _LineList[i].GetP2().GetX() << " " << _LineList[i].GetP2().GetY() << endl;
	}
	file1.close();*/
}

void FlatMesh::CutOffLines()
{
	Vertex v;
	int flag = 1;

	//若线段长度小于0.5，则不进行细分
	int Min_Length = 0.5;

	while (flag == 1)
	{
		flag = 0;
		for (int i = 0; i < m_lineList.size() - 1; i++)
		{
			if (m_lineList[i].P2DLength() > Min_Length)
			{
				for (int j = i + 1; j < m_lineList.size(); j++)
				{
					if (m_lineList[j].P2DLength() > Min_Length)
					{
						bool judge1 = (m_lineList[i].GetV1() == m_lineList[j].GetV2()) || (m_lineList[i].GetV2() == m_lineList[j].GetV1()) || (m_lineList[i].GetV1() == m_lineList[j].GetV1()) || (m_lineList[i].GetV2() == m_lineList[j].GetV2());   //两个顶点重叠
						bool judge2 = ((m_lineList[i].GetV1() == m_lineList[j].GetV1()) && (m_lineList[i].GetV2() == m_lineList[j].GetV2())) || ((m_lineList[i].GetV1() == m_lineList[j].GetV2()) && (m_lineList[i].GetV2() == m_lineList[j].GetV1()));   //同一条线段
						if (m_lineList[i].IntersectionByTwoLines(m_lineList[j], v) > 0 && !judge1)
						{
							flag = 1;
							switch (m_lineList[i].IntersectionByTwoLines(m_lineList[j], v))
							{
							case 1:
							{
								m_lineList.push_back(m_lineList[i]);
								m_lineList[i].SetV2(v);
								m_lineList[m_lineList.size() - 1].SetV1(v);
								m_lineList.push_back(m_lineList[j]);
								m_lineList[j].SetV2(v);
								m_lineList[m_lineList.size() - 1].SetV1(v);
								break;
							}
							case 2:
							{
								if (!judge1 && !judge2)
								{
									m_lineList.push_back(m_lineList[i]);
									m_lineList[i].SetV2(v);
									m_lineList[m_lineList.size() - 1].SetV1(v);
								}
								break;
							}
							case 3:
							{
								if (!judge1 && !judge2)
								{
									m_lineList.push_back(m_lineList[i]);
									m_lineList[i].SetV2(v);
									m_lineList[m_lineList.size() - 1].SetV1(v);
								}
								break;
							}
							case 4:
							{
								if (!judge1 && !judge2)
								{
									m_lineList.push_back(m_lineList[j]);
									m_lineList[j].SetV2(v);
									m_lineList[m_lineList.size() - 1].SetV1(v);
								}
								break;
							}
							case 5:
							{
								if (!judge1 && !judge2)
								{
									m_lineList.push_back(m_lineList[j]);
									m_lineList[j].SetV2(v);
									m_lineList[m_lineList.size() - 1].SetV1(v);
								}
								break;
							}                                                                
							default:
								break;
							}
						}
					}

				}
			}

		}

	}

	//防止线段重复读取
	vector<LineSegment>::iterator it, it1;
	for (it = m_lineList.begin(); it != m_lineList.end(); it++)
		for (it1 = it + 1; it1 != m_lineList.end(); )                                                                                  
		{
			if (*it1 == *it)
			{
				//cout << (*it1).GetP1().GetX() << "    " << (*it1).GetP1().GetY() << endl;
				if ((*it1).GetIsDash() != (*it).GetIsDash())
					(*it).SetIsDash(0);
				it1 = m_lineList.erase(it1);
			}
			else
				++it1;
		}

	//set id
	for (int i = 0; i < m_lineList.size(); i++)
		m_lineList[i].SetId(i + 1);

	//set vertices separate
	for (int i = 0; i < m_lineList.size(); i++)
	{
		if (m_lineList[i].GetIsDash() == 1)
		{
			m_lineList[i].SetV1SeparateId(0);
			m_lineList[i].SetV2SeparateId(0);
		}
	}
}

void FlatMesh::generateVertices(Vertex & v)
{
	//没有找到则生成
	//考虑该点在三维上是否需要分离的情况
	int found = findVertexId(v);
	if (-1 == found)
	{
		m_vertices.push_back(v);
		//if(v.GetIsSeparateIn3D() == 1)
		//	m_vertices.push_back(v); //如果需要分离，则生成同样坐标的两个顶点，加入list
	}
	else
	{
		//找到顶点的话，若IsSeparateIn3D的值不同，则改为0
		if (m_vertices[found].GetIsSeparateIn3D() == 1 && v.GetIsSeparateIn3D() == 0)
			m_vertices[found].SetIsSeparateIn3D(0);
	}
}

int FlatMesh::findVertexId(Vertex & v)
{
	//找到每个顶点v在vertices数组中的序号
	//没有则为-1
	if (0 == m_vertices.size())
		return -1;

	for (int i = 0; i < m_vertices.size(); i++)
	{
		if (m_vertices[i] == v)
			return i;
	}

	//没有找到则返回-1
	return -1;
}

void FlatMesh::generateHalfedgePair(int src, int dst)
{
	Halfedge *h_from, *h_to;
	//生成halfedge
	h_from = new Halfedge(&m_vertices[src], &m_vertices[dst]);
	h_to = new Halfedge(&m_vertices[dst], &m_vertices[src]);
	h_from->set_opposite(h_to);
	h_to->set_opposite(h_from);

	h_from->set_idx(m_halfedges.size());
	m_halfedges.push_back(h_from);

	h_to->set_idx(m_halfedges.size());
	m_halfedges.push_back(h_to);
}

void FlatMesh::initializeData()
{
	//生成顶点序列
	for (int i = 0; i < m_lineList.size(); i++)
	{
		generateVertices(m_lineList[i].GetV1());
		generateVertices(m_lineList[i].GetV2());
	}
	//initialize halfedges

	for (int i = 0; i < m_lineList.size(); i++)
	{
		//note:不能一边新建vertices一边新建halfedges，因为vertices的地址会变，导致之前建立的halfedge有误
		int src, dst;

		src = findVertexId(m_lineList[i].GetV1());
		dst = findVertexId(m_lineList[i].GetV2());
		generateHalfedgePair(src, dst);
		//记录halfedge对应哪条line
		m_whichLineofHalfedge.push_back(i);
		m_whichLineofHalfedge.push_back(i);

		////考虑vertex在三维上需要分离的情况
		//if (m_lines[i].GetV1SeparateId() == 1 && m_lines[i].GetV2SeparateId() == 1)			
		//	generateHalfedgePair(src + 1, dst + 1);
		//else if (m_lines[i].GetV1SeparateId() == 1 && m_lines[i].GetV2SeparateId() == 0)
		//	generateHalfedgePair(src + 1, dst);
		//else if(m_lines[i].GetV1SeparateId() == 0 && m_lines[i].GetV2SeparateId() == 1)
		//	generateHalfedgePair(src, dst+1);		

	}

	//for each vertex, sort o_halfedges
	for (Vertex v : m_vertices)
		v.sort_o_halfedges();
}

bool FlatMesh::HasSeparateVertex(vector<OpenMesh::VertexHandle>& f1)
{
	for (int i = 0; i < f1.size(); i++)
	{
		if (m_vertices[f1[i].idx()].GetIsSeparateIn3D() == 1)
			return true;
	}
	return false;
}

bool FlatMesh::HasSameSeparateVertices(MyMesh & m, vector<OpenMesh::VertexHandle>& f1, vector<OpenMesh::VertexHandle>& f2)
{
	if (HasSeparateVertex(f1) && HasSeparateVertex(f2))
	{
		unordered_set<int> s;
		for (int i = 0; i < f1.size(); i++)
			if (m_vertices[f1[i].idx()].GetIsSeparateIn3D() == 1)
				s.insert(f1[i].idx());

		int flag = 0; //判断是否有相同的需分离的顶点
		for (int i = 0; i < f2.size(); i++)
		{
			if (s.find(f2[i].idx()) != s.end())
			{
				flag = 1;
				Vertex *v = &m_vertices[f2[i].idx()];
				OpenMesh::Vec3f p(v->GetX(), v->GetY(), v->GetZ());
				OpenMesh::VertexHandle vh = m.add_vertex(p);
				v->set_idx(vh.idx());
				f2[i] = vh;
				m_vertices.push_back(*v);
			}
		}
		if (flag == 1)
			return true;
	}

	return false;
}

void FlatMesh::createMesh(string filepath, MyMesh & m)
{
	//m_filePath.copy()
	ReadLinesFromTxt(filepath);

	initializeData();

	m.request_face_colors();
	m.request_face_normals();

	// add vertices
	for (int i = 0; i<m_vertices.size(); i++)
	{
		Vertex *v = &m_vertices[i];
		OpenMesh::Vec3f p(v->GetX(), v->GetY(), v->GetZ());
		OpenMesh::VertexHandle vh = m.add_vertex(p);
		v->set_idx(vh.idx());    //为什么在chengcheng的repo中这个函数是在halfedge的类里？
	}

	//add faces
	vector< vector<Halfedge*> > faces;
	for (Halfedge *h : m_halfedges)
	{
		if (!h->has_face())
		{
			Halfedge* face_halfedge = h;
			vector<Halfedge*> face;
			do {
				face.push_back(face_halfedge);
				face_halfedge->set_has_face();
				face_halfedge = face_halfedge->next();
			} while ((face_halfedge != NULL) && (!face_halfedge->has_face()));
			faces.push_back(face);
		}
	}

	sort(faces.begin(), faces.end(), [](vector<Halfedge*> a, vector<Halfedge*> b) {
		return (a.size() < b.size());
	});
	faces.pop_back();     

	vector<vector<OpenMesh::VertexHandle>> fflist;
	for (vector<Halfedge*> f : faces)
	{
		//不添加由实线包裹的面
		int flag = 0;
		for (Halfedge * h : f)
		{
			if (m_lineList[m_whichLineofHalfedge[h->idx()]].GetIsDash() == 1)  //若平面内有一条线段为虚线，则添加平面
			{
				flag = 1;
				break;
			}
		}

		if (flag == 1)
		{
			vector<OpenMesh::VertexHandle> fvs;
			for (Halfedge * h : f)
				fvs.push_back(OpenMesh::VertexHandle((h->from_vertex())->idx()));

			//对face进行删减，不添加degenerated faces和size = 2的face
			if (fvs.size() > 2)
			{
				std::unordered_set<OpenMesh::VertexHandle> s(fvs.begin(), fvs.end());
				if (s.size() == fvs.size())  //若vector中没有重复顶点
											 //m.add_face(fvs);
					fflist.push_back(fvs);
			}
		}

	}

	for (int i = 0; i<fflist.size() - 1; i++)
		for (int j = i + 1; j < fflist.size(); j++)
		{
			HasSameSeparateVertices(m, fflist[i], fflist[j]);
		}

	for (int i = 0; i < fflist.size(); i++)
		m.add_face(fflist[i]);

	
	//获取预先定义属性
	if (!m.has_vertex_status())
		m.request_vertex_status();
	if (!m.has_face_status())
		m.request_face_status();
	if (!m.has_edge_status())
		m.request_edge_status();

	//delete isolated vertices and face
	m.delete_isolated_vertices();

	vector<MyMesh::FaceHandle> needDelete;
	for (MyMesh::FaceIter f_it = m.faces_begin(); f_it != m.faces_end(); f_it++)
	{
		if (m.ff_begin(*f_it) == m.ff_end(*f_it))
			needDelete.push_back(*f_it);
	}

	for (int i = 0; i < needDelete.size(); i++)
	m.delete_face(needDelete[i],true);


	//执行删除操作
	m.garbage_collection();

	//释放预先定义属性
	if (m.has_vertex_status())
		m.release_vertex_status();/**/
	if (m.has_face_status())
		m.release_face_status();
	if (m.has_edge_status())
		m.release_edge_status();

	//color
	int colorIdx = 0;
	for (MyMesh::FaceIter f_it = m.faces_begin(); f_it != m.faces_end(); ++f_it)
	{
		// set plane color
		float color_r = colorSpace[colorIdx % colorNum][0] * 255;
		float color_g = colorSpace[colorIdx % colorNum][1] * 255;
		float color_b = colorSpace[colorIdx % colorNum][2] * 255;

		m.set_color(*f_it, MyMesh::Color(color_r, color_g, color_b));
		colorIdx++;
	}
}

void FlatMesh::clear_halfedges()
{
	for (auto v : m_vertices)
	{
		vector<Halfedge *> halfedges = v.get_o_halfedges();
		for (Halfedge *h : halfedges)
			delete h;
	}
	for (auto v : m_vertices) v.clear_halfedges();
}
