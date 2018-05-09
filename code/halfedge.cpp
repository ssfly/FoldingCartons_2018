#include "halfedge.h"
#include "Vertex.h"
#include <math.h>

using namespace std;

Halfedge::Halfedge()
{
}

Halfedge::Halfedge(Vertex *from, Vertex *to)
{
	from_v = from;
	to_v = to;
	from->add_o_halfedge(this);
	to->add_i_halfedge(this);
	fidx = -1;
	next_halfedge = 0;
	opposite_halfedge = 0;
	prev_halfedge = 0;
	//hidx = -1;
	hasface = false;
}

double Halfedge::angle()
{
	//计算矢量与x轴正方向的逆时针夹角
	float x1 = to_v->GetX() - from_v->GetX();
	float y1 = to_v->GetY() - from_v->GetY();
	return atan2(y1,x1);

}

