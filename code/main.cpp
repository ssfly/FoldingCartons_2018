#include "CFolding.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	CFolding w;
	w.show();
	return a.exec();
}