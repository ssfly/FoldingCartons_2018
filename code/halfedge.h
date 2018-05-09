#ifndef HALFEDGE_H
#define HALFEDGE_H

class Vertex;

class Halfedge
{
public:
	Halfedge();
	Halfedge(Vertex* , Vertex* );
	Halfedge *next()    const { return next_halfedge; }
	Halfedge *opposite()const { return opposite_halfedge; }
	Halfedge *prev()    const { return prev_halfedge; }
	void set_next(Halfedge *h) { next_halfedge = h; }
	void set_opposite(Halfedge *h) { opposite_halfedge = h; }
	void set_prev(Halfedge *h) { prev_halfedge = h; }
	Vertex* from_vertex() { return from_v; }
	Vertex* to_vertex() { return to_v; }
	int face_idx() { return fidx; }
	void set_face_idx(int f) { fidx = f; }
	void set_idx(int h) { hidx = h; }
	int idx() { return hidx; }
	void set_has_face(bool b = true) { hasface = b; }
	bool has_face() { return hasface; }
	void set_spline(bool b = true) { isspline = b; }
	bool is_spline() { return isspline; }
	double angle();

private:
	Halfedge *next_halfedge,*prev_halfedge;
	Halfedge *opposite_halfedge;
	Vertex *from_v, *to_v;
	int fidx;
	int hidx;
	bool hasface;
	bool isspline;
};

#endif // HALFEDGE_H
