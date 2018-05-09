#pragma once

#include <QtWidgets/QMainWindow>
#include <QtGui>
#include <qlabel.h>
#include <qlayout.h>
#include <qtoolbutton.h>
#include <fstream>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <math.h>
#include <qscrollarea.h>
#include <QGraphicsScene>
#include <QTimer>
#include <vector>
#include <string>
#include <iostream>
#include <qpushbutton.h>
#include "ui_CFolding.h"
#include "Viewer.h"
#include "FlatMesh.h"
#include "StereoMesh.h"
#include "FAlgorithm.h"

class CFolding : public QMainWindow
{
	Q_OBJECT

public:
	CFolding(QWidget *parent = Q_NULLPTR);

private:
	Ui::CFoldingClass ui;
	Viewer *m_view3d;
	Viewer *m_view2d;

	FlatMesh *m_flatMesh;
	StereoMesh *m_steroMesh;
	FAlgorithm *m_algorithm;

	QHBoxLayout *m_candiArea;
	QHBoxLayout *m_faceCandi;

	QString modelPath;
	QString layoutPath;

	vector < vector < MyMesh::VertexHandle >> m_candiVertices;
	vector<vector<MyMesh::FaceHandle>> m_candiFaces; 
	vector<Viewer*> viewCandi;
	vector<Viewer*> faceCandi;
	int canWidth, canHeight;

	qglviewer::Vec bbMax;
	qglviewer::Vec bbMin;

	QWidget *windows;
	Viewer *showSelected;
	Viewer *showAnimation;
	void focusOnObj(MyMesh &mesh, Viewer *view);
	void updateCandidate( MyMesh &mesh);	

	MyMesh pre_2D;
	MyMesh pre_3D;

public slots:
	void loadButton();
	void saveButton();
	void save2DButton();
	void initializeButton();
	void optimizeButton();
	void moveButton();
	void vertexButton();
	void edgeButton();
	void faceButton();
	void symmetryButton();
	void chooseCandidate(int idx);
	void chooseCandidate_face(int idx);
	void optimizationVersa();
	void computeLengthDiff();
	void test();
	void canYes();
	void canNo();
	void canCancel();

	void AnimationButton();
	void UndoButton();
	//void updateL();
};
