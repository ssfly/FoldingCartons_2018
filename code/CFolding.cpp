#include "CFolding.h"
#include <qdebug.h>
#include<fstream>

CFolding::CFolding(QWidget *parent)
	: QMainWindow(parent), m_steroMesh(NULL)
{
	ui.setupUi(this);

	m_flatMesh = new FlatMesh();
	m_algorithm = new FAlgorithm();
	//m_steroMesh = new StereoMesh();

	m_view2d = new Viewer();
	m_view3d = new Viewer();
	ui.show2dLayout->addWidget(m_view2d);
	ui.show3dLayout->addWidget(m_view3d);

	m_candiArea = new QHBoxLayout();
	ui.scrollAreaWidgetContents->setLayout(m_candiArea);
	canWidth = ui.scrollAreaWidgetContents->height();

	m_faceCandi = new QHBoxLayout();
	ui.scrollAreaWidgetContents_2->setLayout(m_faceCandi);

	connect(ui.actionLoad_2, SIGNAL(triggered()), this, SLOT(loadButton()));
	connect(ui.actionSave, SIGNAL(triggered()), this, SLOT(saveButton()));
	connect(ui.actionSave2D, SIGNAL(triggered()), this, SLOT(save2DButton()));
	connect(ui.actionInitialize, SIGNAL(triggered()), this, SLOT(initializeButton()));
	connect(ui.actionOptimize, SIGNAL(triggered()), this, SLOT(optimizeButton()));
	connect(ui.MoveTool, SIGNAL(clicked()), this, SLOT(moveButton()));
	connect(ui.VertexTool, SIGNAL(clicked()), this, SLOT(vertexButton()));
	connect(ui.EdgeTool, SIGNAL(clicked()), this, SLOT(edgeButton()));
	connect(ui.faceTool, SIGNAL(clicked()), this, SLOT(faceButton()));
	connect(ui.SymmetryTool, SIGNAL(clicked()), this, SLOT(symmetryButton()));
	connect(ui.opVersaTool, SIGNAL(clicked()), this, SLOT(optimizationVersa()));
	//connect(ui.actionSave, SIGNAL(triggered()), this, SLOT(saveButton()));
	connect(m_view3d, SIGNAL(updatelayout()), this, SLOT(optimizationVersa()));
	connect(m_view3d, SIGNAL(computeLengthError()), this, SLOT(computeLengthDiff()));

	connect(ui.actionAnimation, SIGNAL(triggered()), this, SLOT(AnimationButton()));
	connect(ui.actionUndo, SIGNAL(triggered()), this, SLOT(UndoButton()));
	connect(ui.actionLoadtest, SIGNAL(triggered()), this, SLOT(test()));
	//connect(ui, SIGNAL(triggered()), this, SLOT());

	
}

void CFolding::test()
{
	QFileDialog *fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle("Open 3D Model");
	fileDialog->setDirectory(".");
	fileDialog->setFilter(QDir::Files | QDir::NoSymLinks);
	if (fileDialog->exec() == QDialog::Accepted)
	{
		modelPath = fileDialog->selectedFiles()[0];
		QMessageBox::information(NULL, tr("Path"), tr("You selected ") + modelPath);
	}
	else
	{
		QMessageBox::information(NULL, tr("Path"), tr("You didn't select any model."));
	}
	m_view3d->mesh.request_vertex_normals();
	m_view3d->mesh.request_face_colors();
	if (!m_view3d->mesh.has_vertex_normals())
	{
		QMessageBox::information(NULL, tr("ERROR"), tr("ERROR: error mesh."));
	}
	OpenMesh::IO::Options opt;
	string pathtmp = modelPath.toStdString();
	if (!OpenMesh::IO::read_mesh(m_view3d->mesh, pathtmp, opt))
	{
		QMessageBox::information(NULL, tr("Path"), tr("Error: Cannot read mesh from") + modelPath);
	}
	if (!opt.check(OpenMesh::IO::Options::VertexNormal))
	{
		m_view3d->mesh.request_face_normals();
		m_view3d->mesh.update_normals();
		m_view3d->mesh.release_face_normals();
	}
	for (MyMesh::FaceIter f_it = m_view3d->mesh.faces_begin(); f_it != m_view3d->mesh.faces_end(); f_it++)
	{
		float r = rand() % 255;
		float g = rand() % 255;
		float b = rand() % 255;
		m_view3d->mesh.set_color(*f_it, MyMesh::Color(r, g, b));
	}
	//showmesh = view3d->mesh;
	double xmin = 100000000.;
	double xmax = -100000000.;
	double ymin = 100000000.;
	double ymax = -100000000.;
	double zmin = 100000000.;
	double zmax = -100000000.;

	for (MyMesh::VertexIter v_it = m_view3d->mesh.vertices_begin(); v_it != m_view3d->mesh.vertices_end(); ++v_it)
	{
		xmin = m_view3d->mesh.point(*v_it).data()[0] < xmin ? m_view3d->mesh.point(*v_it).data()[0] : xmin;
		ymin = m_view3d->mesh.point(*v_it).data()[1] < ymin ? m_view3d->mesh.point(*v_it).data()[1] : ymin;
		zmin = m_view3d->mesh.point(*v_it).data()[2] < zmin ? m_view3d->mesh.point(*v_it).data()[2] : zmin;
		xmax = m_view3d->mesh.point(*v_it).data()[0] > xmax ? m_view3d->mesh.point(*v_it).data()[0] : xmax;
		ymax = m_view3d->mesh.point(*v_it).data()[1] > ymax ? m_view3d->mesh.point(*v_it).data()[1] : ymax;
		zmax = m_view3d->mesh.point(*v_it).data()[2] > zmax ? m_view3d->mesh.point(*v_it).data()[2] : zmax;
	}
	//bbMax = qglviewer::Vec(xmax, ymax, zmax);
	//bbMin = qglviewer::Vec(xmin, ymin, zmin);
	m_view3d->getEdgeVertex();
	m_view3d->camera()->setRevolveAroundPoint(qglviewer::Vec(.5*(xmin + xmax), .5*(ymin + ymax), .5*(zmin + zmax)));
	m_view3d->camera()->setSceneBoundingBox(qglviewer::Vec(xmin, ymin, zmin), qglviewer::Vec(xmax, ymax, zmax));
	m_view3d->camera()->showEntireScene();
	m_view3d->update();
}

void CFolding::focusOnObj(MyMesh &mesh, Viewer *view)
{
	double xmin = 100000000.;
	double xmax = -100000000.;
	double ymin = 100000000.;
	double ymax = -100000000.;
	double zmin = 100000000.;
	double zmax = -100000000.;
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		xmin = mesh.point(*v_it).data()[0] < xmin ? mesh.point(*v_it).data()[0] : xmin;
		ymin = mesh.point(*v_it).data()[1] < ymin ? mesh.point(*v_it).data()[1] : ymin;
		zmin = mesh.point(*v_it).data()[2] < zmin ? mesh.point(*v_it).data()[2] : zmin;
		xmax = mesh.point(*v_it).data()[0] > xmax ? mesh.point(*v_it).data()[0] : xmax;
		ymax = mesh.point(*v_it).data()[1] > ymax ? mesh.point(*v_it).data()[1] : ymax;
		zmax = mesh.point(*v_it).data()[2] > zmax ? mesh.point(*v_it).data()[2] : zmax;
	}
	bbMax = qglviewer::Vec(xmax, ymax, zmax);
	bbMin = qglviewer::Vec(xmin, ymin, zmin);
	view->camera()->setRevolveAroundPoint(qglviewer::Vec(.5 * (bbMin[0] + bbMax[0]), .5 * (bbMin[1]+ bbMax[1]), .5 * (bbMin[2] + bbMax[2])));
	view->camera()->setSceneBoundingBox(bbMin, bbMax);
	view->camera()->showEntireScene();
}
void CFolding::loadButton()
{
	QFileDialog *fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle("Open 2D layout");
	fileDialog->setDirectory(".");
	fileDialog->setFilter(QDir::Files | QDir::NoSymLinks);
	if (fileDialog->exec() == QDialog::Accepted)
	{
		layoutPath = fileDialog->selectedFiles()[0];
		QMessageBox::information(NULL, tr("Path"), tr("You selected ") + layoutPath);
	}
	else
	{
		QMessageBox::information(NULL, tr("Path"), tr("You didn't select any layout. "));
		return;
	}
	if (m_view2d->isUsed)
	{
		m_view2d->clearViewer();
		m_view3d->clearViewer();
		m_flatMesh = new FlatMesh();
		m_algorithm = new FAlgorithm();
		for (int i = 0; i < viewCandi.size(); i++)
		{
			m_candiArea->removeWidget(viewCandi[i]);
			delete viewCandi[i];
		}
		viewCandi.clear();
		m_candiVertices.clear();
		if (m_steroMesh != nullptr) delete m_steroMesh;
	}
	m_flatMesh->createMesh(layoutPath.toStdString(), m_view2d->mesh);
	m_view2d->isUsed = true;
	//////////////////////////////////////////
	//flatMesh->readLineListFromTxt(layoutPath.toStdString);
	//Get the 2d mesh;
	//m_view2d->mesh = ...... ////send the mesh to the class view2d;
	focusOnObj(m_view2d->mesh, m_view2d);
	m_view2d->needUpdateColor = false;
	m_view2d->getEdgeVertex();
	m_view2d->update();
	//////////////////////////////////////////
}

void CFolding::updateCandidate(MyMesh &mesh)
{
	if (!viewCandi.empty())
	{
		for (int i = 0; i < viewCandi.size(); i++)
		{
			m_candiArea->removeWidget(viewCandi[i]);
			delete viewCandi[i];
		}
		viewCandi.clear();
		m_candiVertices.clear();
	}
	if (!faceCandi.empty())
	{
		for (int i = 0; i < faceCandi.size(); i++)
		{
			m_faceCandi->removeWidget(faceCandi[i]);
			delete faceCandi[i];
		}
		faceCandi.clear();
		m_candiFaces.clear();
	}
 	m_candiVertices = m_algorithm->autoDetection2(mesh);
	int n = m_candiVertices.size();
	//vector<MyMesh::VertexHandle> candiTmp;
	ui.scrollAreaWidgetContents->setFixedSize(n * canWidth, canWidth);
	viewCandi.resize(n);

	m_candiFaces = m_algorithm->autoDectionofFace(mesh);
	int m = m_candiFaces.size();
	ui.scrollAreaWidgetContents_2->setFixedSize(m * canWidth, canWidth);
	faceCandi.resize(m);

	for (int i = 0; i < n; i++)
	{
		viewCandi[i] = new Viewer();
		viewCandi[i]->viewIdx = i;
		viewCandi[i]->resize(canWidth, canWidth);
		viewCandi[i]->mesh = m_view3d->mesh;
		viewCandi[i]->isCandidate = true;
		m_candiArea->addWidget(viewCandi[i]);
		for (int j = 0; j < m_candiVertices[i].size(); j++)
		{
			viewCandi[i]->mergeVertices.push_back(m_candiVertices[i][j].idx());
		}
		focusOnObj(viewCandi[i]->mesh, viewCandi[i]);
		viewCandi[i]->getEdgeVertex();
		connect(viewCandi[i], SIGNAL(isDoubleClicked(int)), this, SLOT(chooseCandidate(int)));
		viewCandi[i]->update();
	}

	for (int i = 0; i < m; i++)
	{
		faceCandi[i] = new Viewer();
		faceCandi[i]->viewIdx = i;
		faceCandi[i]->resize(canWidth, canWidth);
		faceCandi[i]->mesh = m_view3d->mesh;
		faceCandi[i]->isCandidate = true;
		m_faceCandi->addWidget(faceCandi[i]);
		for (int j = 0; j < m_candiFaces[i].size(); j++)
		{
			faceCandi[i]->mergeFaces.push_back(m_candiFaces[i][j].idx());
		}
		focusOnObj(faceCandi[i]->mesh, faceCandi[i]);
		faceCandi[i]->getEdgeVertex();
		connect(faceCandi[i], SIGNAL(isDoubleClicked(int)), this, SLOT(chooseCandidate_face(int)));
		faceCandi[i]->update();
	}
}

void CFolding::computeLengthDiff()
{
	float diff = 0;
	MyMesh::EdgeIter me_it = m_view3d->mesh.edges_begin();
	for (MyMesh::EdgeIter e_it = m_view2d->mesh.edges_begin(); e_it != m_view2d->mesh.edges_end(); e_it++)
	{
		diff += fabs(m_view2d->mesh.calc_edge_length(*e_it) - m_view3d->mesh.calc_edge_length(*me_it));
		me_it++;
	}
	diff /= m_view2d->mesh.n_edges();

	ofstream out("lengthDiff.txt",ios::app);
	out << layoutPath.data() << "    " << diff << endl;
	out.close();
}

void CFolding::chooseCandidate_face(int idx)
{
	vector<MyMesh::FaceHandle> merge_f;
	for (int i = 0; i < faceCandi[idx]->mergeFaces.size(); i++)
	{
		MyMesh::FaceHandle htmp(faceCandi[idx]->mergeFaces[i]);
		merge_f.push_back(htmp);
	}
	m_view3d->mesh = m_algorithm->ShapeOptimization(m_view3d->mesh, merge_f);

	m_view3d->getEdgeVertex();
	m_view3d->m_label = m_algorithm->GetClassForSameVertexin3D(m_view3d->mesh);

	focusOnObj(m_view3d->mesh, m_view3d);
	m_view3d->update();
	updateCandidate(m_view3d->mesh);
	m_view2d->setColorIdxAs(m_view3d->getColorIdx());
	m_view2d->update();

}
void CFolding::saveButton()
{
	m_view3d->saveSnapshot(false, false);
}

void CFolding::save2DButton()
{
	m_view2d->saveSnapshot(false, false);
}

void CFolding::initializeButton()
{
	m_steroMesh = new StereoMesh(m_view2d->mesh);
	m_view3d->mesh = m_steroMesh->folding();
	m_view3d->getEdgeVertex();
	m_view3d->m_label = m_algorithm->GetClassForSameVertexin3D(m_view3d->mesh);
	focusOnObj(m_view3d->mesh, m_view3d);
	m_view3d->update();

	m_view2d->setColorIdxAs(m_view3d->getColorIdx());
	m_view2d->update();
	qDebug() << m_view3d->getColorIdx();
	updateCandidate(m_view3d->mesh);
 //	m_candiVertices = m_algorithm->autoDetection2(m_view3d->mesh);
	//int n = m_candiVertices.size();
	//vector<MyMesh::VertexHandle> candiTmp;
	//ui.scrollAreaWidgetContents->setFixedSize( n * canWidth, canWidth);
	//viewCandi.resize(n);

	//for (int i = 0; i < n; i++)
	//{
	//	viewCandi[i] = new Viewer();
	//	viewCandi[i]->viewIdx = i;
	//	viewCandi[i]->resize(canWidth, canWidth);
	//	viewCandi[i]->mesh = m_view3d->mesh;
	//	viewCandi[i]->isCandidate = true;
	//	m_candiArea->addWidget(viewCandi[i]);
	//	for (int j = 0; j < m_candiVertices[i].size(); j++)
	//	{
	//		viewCandi[i]->mergeVertices.push_back(m_candiVertices[i][j].idx());
	//	}
	//	focusOnObj(viewCandi[i]->mesh, viewCandi[i]);
	//	viewCandi[i]->getEdgeVertex();
	//	connect(viewCandi[i], SIGNAL(isDoubleClicked(int)), this, SLOT(chooseCandidate(int)));
	//	viewCandi[i]->update();
	//}

	//initialize solver, add basic constrains
 	m_algorithm->InitializeShapeConstrain(m_view3d->mesh);

	pre_2D = m_view2d->mesh;
	pre_3D = m_view3d->mesh;

	//m_algorithm->GetClassForSameVertexin3D(m_view3d->mesh);

	///////////////////////////////////////////////
	//get the initial model
	//m_view3d->mesh = .......
	//m_view3d->update;
	//get the mergeVerticesList
	//vector<Viewer*> viewCandi(n);
	///////////////////////////////////////////////
}
void CFolding::chooseCandidate(int idx)
{

	vector<MyMesh::VertexHandle> merge_v;
	for (int i = 0; i < viewCandi[idx]->mergeVertices.size(); i++)
	{
		MyMesh::VertexHandle htmp(viewCandi[idx]->mergeVertices[i]);
		merge_v.push_back(htmp);
	}
	vector<vector<MyMesh::VertexHandle>> tmpSD = m_algorithm->modifiedSymmetryDetection(m_view3d->mesh, merge_v);
	if (!tmpSD.empty())
	{
		windows = new QWidget();
		windows->setAttribute(Qt::WA_DeleteOnClose);
		showSelected = new Viewer();
		showSelected->mesh = m_view3d->mesh;
		showSelected->mergeVertices = viewCandi[idx]->mergeVertices;

		for (int i = 0; i < tmpSD[0].size(); i++)
		{
			showSelected->symVertices.push_back(tmpSD[0][i].idx());
		}
		showSelected->showMergeVertices = true;

		QVBoxLayout *layout = new QVBoxLayout();
		QPushButton *yes = new QPushButton();
		yes->setText("YES");
		QPushButton *no = new QPushButton();
		no->setText("NO");
		QHBoxLayout *button = new QHBoxLayout();
		QPushButton *cancel = new QPushButton();
		cancel->setText("CANCEL");
		button->addWidget(yes);
		button->addWidget(no);
		button->addWidget(cancel);
		button->setSpacing(100);
		button->setContentsMargins(QMargins(150, 5, 150, 5));
		//button->setMargin(100);

		//QVBoxLayout *llay = new QVBoxLayout();
		QLabel *label = new QLabel();
		QFont ft("Microsoft YaHei", 15, 10);
		//ft.setPointSize(15);
		label->setFont(ft);
		label->setText("Symetric vertex set is detected as shown by blue points.<br/>Do you want to optimize the carton with these vertices merging?");
		label->setAlignment(Qt::AlignHCenter);
		QPalette pa;
		pa.setColor(QPalette::Background, QColor(255, 255, 255));
		label->setAutoFillBackground(true);
		label->setPalette(pa);
		// ay->addWidget(label);
		// ay->setContentsMargins(QMargins(10, 0, 10, 0));

		showSelected->getEdgeVertex();
		focusOnObj(showSelected->mesh, showSelected);

		layout->addWidget(showSelected);
		layout->addWidget(label);
		//layout->addLayout(llay);
		layout->addLayout(button);
		layout->setStretch(0, 10);
		layout->setStretch(1, 1);
		layout->setStretch(2, 2);
		layout->setSpacing(10);

		windows->setLayout(layout);
		windows->show();
		connect(yes, SIGNAL(clicked()), this, SLOT(canYes()));
		connect(no, SIGNAL(clicked()), this, SLOT(canNo()));
		connect(cancel, SIGNAL(clicked()), this, SLOT(canCancel()));

	}
	else
	{
		
		pre_2D = m_view2d->mesh;
		pre_3D = m_view3d->mesh;

		m_view3d->mesh = m_algorithm->ShapeOptimization(m_view3d->mesh, merge_v);
		m_view3d->getEdgeVertex();
		m_view3d->m_label = m_algorithm->GetClassForSameVertexin3D(m_view3d->mesh);
		focusOnObj(m_view3d->mesh, m_view3d);
		m_view3d->update();
		updateCandidate(m_view3d->mesh);

		m_view2d->mesh = m_algorithm->LayoutOptimization(m_view2d->mesh, m_view3d->mesh, merge_v);
		m_view2d->getEdgeVertex();
		m_view2d->setColorIdxAs(m_view3d->getColorIdx());
		focusOnObj(m_view2d->mesh, m_view2d);
		m_view2d->update();

	}
	// optimize on chosen candidate
	//m_view3d->mesh = m_algorithm->ShapeOptimization(m_view3d->mesh, merge_v);
	//m_view3d->getEdgeVertex();
	//focusOnObj(m_view3d->mesh, m_view3d);
	//m_view3d->update();
}
void CFolding::canYes()
{
	vector<int> mergeVertices = showSelected->mergeVertices;
	vector<int> symVertics = showSelected->symVertices;
	windows->close();
	vector<MyMesh::VertexHandle> vhs1;
	for (int i = 0; i < mergeVertices.size(); i++)
	{
		MyMesh::VertexHandle vh(mergeVertices[i]);
		vhs1.push_back(vh);
	}

	vector<MyMesh::VertexHandle> vhs2;
	for (int i = 0; i < symVertics.size(); i++)
	{
		MyMesh::VertexHandle vh(symVertics[i]);
		vhs2.push_back(vh);
	}

	pre_2D = m_view2d->mesh;
	pre_3D = m_view3d->mesh;

	m_view3d->mesh  = m_algorithm->ShapeOptimization(m_view3d->mesh, vhs1, vhs2);
	m_view3d->getEdgeVertex();
	m_view3d->m_label = m_algorithm->GetClassForSameVertexin3D(m_view3d->mesh);
	focusOnObj(m_view3d->mesh, m_view3d);
	m_view3d->update();
	updateCandidate(m_view3d->mesh);

	vector<MyMesh::VertexHandle> vsall = vhs1;
	vsall.insert(vsall.end(), vhs2.begin(), vhs2.end());
	m_view2d->mesh = m_algorithm->LayoutOptimization(m_view2d->mesh, m_view3d->mesh, vsall);
	m_view2d->getEdgeVertex();
	m_view2d->setColorIdxAs(m_view3d->getColorIdx());
	focusOnObj(m_view2d->mesh, m_view2d);
	m_view2d->update();

}

void CFolding::canNo()
{
	vector<int> mergeVertices = showSelected->mergeVertices;
	windows->close();

	vector<MyMesh::VertexHandle> vhs1;
	for (int i = 0; i < mergeVertices.size(); i++)
	{
		MyMesh::VertexHandle vh(mergeVertices[i]);
		vhs1.push_back(vh);
	}

	pre_2D = m_view2d->mesh;
	pre_3D = m_view3d->mesh;

	m_view3d->mesh = m_algorithm->ShapeOptimization(m_view3d->mesh, vhs1);
	m_view3d->getEdgeVertex();
	m_view3d->m_label = m_algorithm->GetClassForSameVertexin3D(m_view3d->mesh);
	focusOnObj(m_view3d->mesh, m_view3d);
	m_view3d->update();
	updateCandidate(m_view3d->mesh);

	m_view2d->mesh = m_algorithm->LayoutOptimization(m_view2d->mesh, m_view3d->mesh, vhs1);
	m_view2d->getEdgeVertex();
	m_view2d->setColorIdxAs(m_view3d->getColorIdx());
	focusOnObj(m_view2d->mesh, m_view2d);
	m_view2d->update();
}

void CFolding::canCancel()
{
	windows->close();
}
void CFolding::AnimationButton()
{
	windows = new QWidget();
	windows->setAttribute(Qt::WA_DeleteOnClose);
	showAnimation = new Viewer();
	showAnimation->mesh = m_view2d->mesh;
	showAnimation->flat_mesh = m_view2d->mesh;
	showAnimation->isAnimation = true;

	QVBoxLayout *layout = new QVBoxLayout();

	showAnimation->getEdgeVertex();
	showAnimation->setColorIdxAs(m_view3d->getColorIdx());
	focusOnObj(showAnimation->mesh, showAnimation);

	layout->addWidget(showAnimation);

	windows->setLayout(layout);
	windows->show();

}

void CFolding::UndoButton()
{
	m_view3d->mesh = pre_3D;
	m_view3d->getEdgeVertex();
	m_view3d->m_label = m_algorithm->GetClassForSameVertexin3D(m_view3d->mesh);
	focusOnObj(m_view3d->mesh, m_view3d);
	m_view3d->update();
	updateCandidate(m_view3d->mesh);

	m_view2d->mesh = pre_2D;
	m_view2d->getEdgeVertex();
	m_view2d->setColorIdxAs(m_view3d->getColorIdx());
	focusOnObj(m_view2d->mesh, m_view2d);
	m_view2d->update();
}

void CFolding::optimizeButton()
{
	////////////////////////////////////////////
	//optimize()
	//m_view3d->mesh = ...
	//m_view3d->update();
	////////////////////////////////////////////
	

	if (m_view3d->selected_vertices.size() > 0)
	{
		vector<MyMesh::VertexHandle> vhs;
		for (int i = 0; i < m_view3d->selected_vertices.size(); i++)
		{
			MyMesh::VertexHandle vh((m_view3d->selected_vertices).at(i));
			vhs.push_back(vh);
		}
		vector<vector<MyMesh::VertexHandle>> tmpSD = m_algorithm->modifiedSymmetryDetection(m_view3d->mesh, vhs);
		if (!tmpSD.empty())
		{
			windows = new QWidget();
			windows->setAttribute(Qt::WA_DeleteOnClose);
			showSelected = new Viewer();
			showSelected->mesh = m_view3d->mesh;
			for (int i = 0; i < m_view3d->selected_vertices.size(); i++)
			{
				showSelected->mergeVertices.push_back(m_view3d->selected_vertices.at(i));
			}
			//showSelected->mergeVertices = m_view3d->selected_vertices;

			for (int i = 0; i < tmpSD[0].size(); i++)
			{
				showSelected->symVertices.push_back(tmpSD[0][i].idx());
			}
			showSelected->showMergeVertices = true;

			QVBoxLayout *layout = new QVBoxLayout();
			QPushButton *yes = new QPushButton();
			yes->setText("YES");
			QPushButton *no = new QPushButton();
			no->setText("NO");
			QHBoxLayout *button = new QHBoxLayout();
			QPushButton *cancel = new QPushButton();
			cancel->setText("CANCEL");
			button->addWidget(yes);
			button->addWidget(no);
			button->addWidget(cancel);
			button->setSpacing(100);
			button->setContentsMargins(QMargins(150, 5, 150, 5));
			//button->setMargin(100);

			//QVBoxLayout *llay = new QVBoxLayout();
			QLabel *label = new QLabel();
			QFont ft("Microsoft YaHei", 15, 10);
			//ft.setPointSize(15);
			label->setFont(ft);
			label->setText("Symetric vertex set is detected as shown by blue points.<br/>Do you want to optimize the carton with these vertices merging?");
			label->setAlignment(Qt::AlignHCenter);
			QPalette pa;
			pa.setColor(QPalette::Background, QColor(255, 255, 255));
			label->setAutoFillBackground(true);
			label->setPalette(pa);
			// ay->addWidget(label);
			// ay->setContentsMargins(QMargins(10, 0, 10, 0));

			showSelected->getEdgeVertex();
			focusOnObj(showSelected->mesh, showSelected);

			layout->addWidget(showSelected);
			layout->addWidget(label);
			//layout->addLayout(llay);
			layout->addLayout(button);
			layout->setStretch(0, 10);
			layout->setStretch(1, 1);
			layout->setStretch(2, 2);
			layout->setSpacing(10);

			windows->setLayout(layout);
			windows->show();
			connect(yes, SIGNAL(clicked()), this, SLOT(canYes()));
			connect(no, SIGNAL(clicked()), this, SLOT(canNo()));
			connect(cancel, SIGNAL(clicked()), this, SLOT(canCancel()));
		}
		else
		{

			pre_2D = m_view2d->mesh;
			pre_3D = m_view3d->mesh;

			m_view3d->mesh = m_algorithm->ShapeOptimization(m_view3d->mesh, vhs);
			m_view3d->m_label = m_algorithm->GetClassForSameVertexin3D(m_view3d->mesh);
		}
	}

	if (m_view3d->selected_faces.size() > 0)
	{
		vector<MyMesh::FaceHandle> fhs;
		for (int i = 0; i < m_view3d->selected_faces.size(); i++)
		{
			MyMesh::FaceHandle fh((m_view3d->selected_faces).at(i));
			fhs.push_back(fh);
		}

		pre_2D = m_view2d->mesh;
		pre_3D = m_view3d->mesh;

		m_view3d->mesh = m_algorithm->ShapeOptimization(m_view3d->mesh, fhs);


	}
	

	m_view3d->getEdgeVertex();
	m_view3d->m_label = m_algorithm->GetClassForSameVertexin3D(m_view3d->mesh);
	m_view2d->setColorIdxAs(m_view3d->getColorIdx());
	m_view2d->update();
	focusOnObj(m_view3d->mesh, m_view3d);
	m_view3d->update();
	updateCandidate(m_view3d->mesh);

}

void CFolding::optimizationVersa()
{

	pre_2D = m_view2d->mesh;
	pre_3D = m_view3d->mesh;

	if (m_view3d->selected_edges.size() == 1)
	{
		MyMesh::VertexHandle vh1(m_view3d->getEdgeVertex0()[m_view3d->selected_edges.at(0)]);
		MyMesh::VertexHandle vh2(m_view3d->getEdgeVertex1()[m_view3d->selected_edges.at(0)]);
		m_view2d->mesh = m_algorithm->LayoutOptimization(m_view2d->mesh, m_view3d->mesh, vh1, vh2);
	}

	m_view2d->getEdgeVertex();
	focusOnObj(m_view2d->mesh, m_view2d);
	m_view2d->setColorIdxAs(m_view3d->getColorIdx());
	m_view2d->update();

	//OpenMesh::IO::Options wopt;
	//wopt += OpenMesh::IO::Options::FaceColor;
	////output 3D
	//try
	//{
	//	if (!OpenMesh::IO::write_mesh(m_view3d->mesh,"editModel.off", wopt)) //"./data/ShapeOp/op-" + g_parameters.InputFileName
	//	{
	//		std::cerr <<"  Initial Mesh Cannot be writen" << std::endl;
	//		return;
	//	}
	//}
	//catch (std::exception& x)
	//{
	//	std::cerr << x.what() << std::endl;
	//	return;
	//}

	//try
	//{
	//	if (!OpenMesh::IO::write_mesh(m_view2d->mesh, "editLayout.off", wopt)) //"./data/ShapeOp/op-" + g_parameters.InputFileName
	//	{
	//		std::cerr << "  Initial Mesh Cannot be writen" << std::endl;
	//		return;
	//	}
	//}
	//catch (std::exception& x)
	//{
	//	std::cerr << x.what() << std::endl;
	//	return;
	//}
}

void CFolding::moveButton()
{
	m_view3d->type_ = Viewer::Type::NONE_;
	m_view3d->selected_edges.clear();
	m_view3d->selected_faces.clear();
	m_view3d->selected_vertices.clear();
}

void CFolding::vertexButton()
{
	//set the choose mode as VERTEX
	m_view3d->type_ = Viewer::Type::VERTEX;  
	m_view3d->selected_edges.clear();
	m_view3d->selected_faces.clear();
	m_view3d->selected_vertices.clear();
}

void CFolding::edgeButton()
{
	//set the choose mode as EDGE
	m_view3d->type_ = Viewer::Type::EDGE;
	m_view3d->selected_edges.clear();
	m_view3d->selected_faces.clear();
	m_view3d->selected_vertices.clear();
}

void CFolding::faceButton()
{
	m_view3d->type_ = Viewer::Type::FACE;
	m_view3d->selected_edges.clear();
	m_view3d->selected_faces.clear();
	m_view3d->selected_vertices.clear();
}
void CFolding::symmetryButton()
{
	//m_view3d->selected_vertices
	//test for symmetry vertex pair
	
	for (int i = 0; i < viewCandi.size(); i++)
		m_candiArea->removeWidget(viewCandi[i]);

	/*MyMesh::VertexHandle vh1((m_view3d->selected_vertices).at(0));
	MyMesh::VertexHandle vh2((m_view3d->selected_vertices).at(1));

	m_candiVertices = m_algorithm->symmetricVertexPair(m_view3d->mesh, vh1 , vh2);*/

	vector<MyMesh::VertexHandle> vhs;
	for (int i = 0; i < m_view3d->selected_vertices.size(); i++)
	{
		MyMesh::VertexHandle vh((m_view3d->selected_vertices).at(i));
		vhs.push_back(vh);
	}
	//m_candiVertices = m_algorithm->symmetryDetection(m_view3d->mesh, vhs);
	m_candiVertices = m_algorithm->modifiedSymmetryDetection(m_view3d->mesh, vhs);

	int n = m_candiVertices.size();
	vector<MyMesh::VertexHandle> candiTmp;
	ui.scrollAreaWidgetContents->setFixedSize(n * canWidth, canWidth);
	viewCandi.resize(n);

	for (int i = 0; i < n; i++)
	{
		viewCandi[i] = new Viewer();
		viewCandi[i]->resize(canWidth, canWidth);
		viewCandi[i]->mesh = m_view3d->mesh;
		viewCandi[i]->isCandidate = true;
		m_candiArea->addWidget(viewCandi[i]);
		for (int j = 0; j < m_candiVertices[i].size(); j++)
		{
			viewCandi[i]->mergeVertices.push_back(m_candiVertices[i][j].idx());
		}
		viewCandi[i]->getEdgeVertex();
		focusOnObj(viewCandi[i]->mesh, viewCandi[i]);
		viewCandi[i]->update();
	}
}
