//
//  sofa-bounds.cpp
//  SofaBounds Version 1.0
//
//  This source file is part of the SofaBounds software package
//  SofaBounds is a program that proves upper bounds in the moving sofa problem,
//  as described in the paper "Improved upper bounds in the moving sofa problem",
//  by Yoav Kallus and Dan Romik
//
//  For more information go to the project web page: *** TODO: add URL ***
//
//  Copyright © 2017 Yoav Kallus and Dan Romik
//

#include "sofa-bounds.hpp"

std::string rat_str(ExactRational x){
    //convert rational to string
    char mystring[500];
    mpq_get_str(mystring,10,x.mpq());
    std::string the_string(mystring);
    return the_string;
}

void split_interval(struct interval b, struct interval *b1, struct interval *b2){
    //split an interval in two
    b1->left = b.left; b1->right = (b.right + b.left) / 2;
    b2->left = (b.right + b.left) / 2; b2->right = b.right;
}

struct interval midpoint_of_interval(struct interval b){
    //return a zero-length interval at center of given interval
    struct interval r;
    r.left = r.right = (b.left + b.right)/2;
    return r;
}

void get_point_coordinates(Point p, ExactRational *x, ExactRational *y){
    //compute a cartesian representation of a CGAL point which might be in non-cartesian representation (e.g. homogeneous coordinates)
    *x = ExactRational(p.hx(),p.hw());
    *y = ExactRational(p.hy(),p.hw());;
}

struct box a_priori_bounds(struct bb_thread_params slopes, unsigned int i){
    //compute a priori bounds on the translation variables based on the angle parameters
    //the optimization problem is guaranteed to obtain a maximum within this bounded region
    //if beta = pi/2, i is the index of the corridor one coordinate of which will be fixed to remove irrelevant degree of freedom
    struct box newbox;
    ExactRational x1,x2,s,c;
    unsigned int j;

    newbox.coord_bound_intervals.resize(2*slopes.num_intermediate);
    //first restrict the intersection to a finite width of the horizontal strip, [x1,x2]X[0,1]
    if (!slopes.has_final) {
	//if beta = pi/2, this is done by fixing the u (second) coordinate for one angle at zero,
	//and restricting the v (first) coordinate so that the intersection of the corridor with H
	//is neither contained within another such intersection nor disconnected.
        newbox.coord_bound_intervals[2*i+0].left = 0;
        newbox.coord_bound_intervals[2*i+0].right = ExactRational(slopes.intermediate[i].c,slopes.intermediate[i].b);

        newbox.coord_bound_intervals[2*i+1].left = 0;
        newbox.coord_bound_intervals[2*i+1].right = 0;

	x2 = ExactRational(slopes.intermediate[i].c,slopes.intermediate[i].b);
	x1 = -x2 * ExactRational(slopes.intermediate[i].c + slopes.intermediate[i].b,slopes.intermediate[i].a);
    } else {
	//if beta < pi/2, this is done by finding smallest width that contains the intersection of H and the butterfly set
	//set i to value such that no angle will be skipped in next loop
	i=slopes.num_intermediate;
	x1 = ExactRational(-slopes.final_max.a,slopes.final_max.b);
	x2 = ExactRational(slopes.final_max.c,slopes.final_max.b);
    }

    for (j=0;j<slopes.num_intermediate;j++) {
	//with x1,x2 determined, all corridors (except i, which has already been restricted)
	//can be restricted to translations such that the intersection with [x1,x2]X[0,1]
	//is not contained within another such intersection.
	if (j==i) continue;
	s = ExactRational(slopes.intermediate[j].a,slopes.intermediate[j].c);
	c = ExactRational(slopes.intermediate[j].b,slopes.intermediate[j].c);

        newbox.coord_bound_intervals[2*j+0].left = -s*x2 + (1-s) * s/c;
        newbox.coord_bound_intervals[2*j+0].right = -s*x1 + c - 1;

        newbox.coord_bound_intervals[2*j+1].left = c*x1 + (1-c) * c/s;
        newbox.coord_bound_intervals[2*j+1].right = c*x2 + s - 1;
    }
    newbox.depth = 0;

    return newbox;
}

ExactRational initial_lower_bound(struct slope m){
    //beta dependent initial lower bound
    //ExactRational a1(15,8), a2(19,12);
    //ExactRational s(m.a,m.c), c(m.b,m.c);

    return ExactRational(0,1);
    /*if (10*m.a < 9*m.c) return ExactRational(3,4);
    return a1*( s - (c-1) ) - a2;*/
}

ExactRational area(Nef_polygon N, std::vector<ExactRational> *polygon) {
    // Calculates the area of a NEF polygon
    // Optionally, if polygon != NULL, store the vertex coordinates in this vector
    ExactRational a,x1,x2,y1,y2;
    Explorer::Face_const_iterator fi;
    Explorer::Halfedge_const_handle hh,hh0;
    Explorer::Vertex_const_handle vh;
    Point p;
    int j;

    Explorer E = N.explorer();
    
    a = 0;
    //iterate over faces
    for (fi=E.faces_begin(); fi!=E.faces_end(); ++fi){
        //skip face if excluded
        if ( ((int) E.mark(fi)) == 0 ) continue;
        //get halfedge to start iterating on
        hh = E.halfedge(fi);
        //if the face is not bounded (even by the box at infinity), skip it
        if (hh == Explorer::Halfedge_const_handle()) continue;
        //iterate over edges in cyclic order
        hh0=hh; j=0;

        vh = E.source(hh);
        // if edge source is not a point at infinity
        if (E.is_standard(vh)) {
            //get affine embedding of the point
            p = E.point(vh);
            get_point_coordinates(p,&x2,&y2);
            //store current point as previous point
            x1 = x2; y1 = y2;
        }
        else {
            //if any point is at infinity, area is infinite
            a = -1;
            return a;
        }
        hh = E.next(hh);
        j++;

	if (polygon!=NULL) polygon->clear();

        do {
            vh = E.source(hh);
            // if edge source is not a point at infinity
            if (E.is_standard(vh)) {
                //get affine embedding of the point
                p = E.point(vh);
		get_point_coordinates(p,&x2,&y2);
		if (polygon!=NULL) {
		    polygon->push_back(x2);
		    polygon->push_back(y2);
		}
                //add twice oriented area of the triangle formed by origin, this point, and previous point
                a += x1*y2 - x2*y1;
                //uncomment for debug output
                //fprintf(stderr,"%lf (%lf,%lf)\n",CGAL::to_double(a),CGAL::to_double(x1),CGAL::to_double(y1));
                //store current point as previous point
                x1 = x2; y1 = y2;
            }
            else {
                //if any point is at infinity, area is infinite
                a = -1;
                return a;
            }
            hh = E.next(hh);
            j++;
        } while (j==1 || hh != E.next(hh0));
    }
    //divide by two to get actual area
    return a/2;
}

inline Nef_polygon ell_side(struct slope myslope, struct interval myinterval, int updown, int inout, int lowhi, int overunder) {
    // return a half plane defining a sweep of rotated ells
    // b = interval of inner corner coord
    // updown = 0 ~ up (left) arm of the ell; 1 ~ down (right)
    // inout = 0 ~ line through inner (0) or outer corner
    // lowhi = lowest or highest line (corresponds to left and right ends of interval b)
    // overunder = the halfplane over or under this line

    ExactInteger ma,mb,mc;
    ExactInteger ai,bi,ci;
    ExactRational x;

    ma = myslope.a;
    mb = myslope.b;
    mc = myslope.c;

    //determine the equation of the line before translation ai*X + bi*Y = ci
    if (updown == 0) {
        ai = -ma; bi = mb;
    } else {
        ai = mb; bi = ma;
    }

    if (inout == 1)  ci = -mc;
    else ci = 0;

    //now, transform ai,bi,ci to represent translation ci -> ci - x
    //because x is rational, need to multiply out by its denominator
    if (lowhi == 0) x = myinterval.left;
    else x = myinterval.right;

    ci *= x.denominator();

    ci -= x.numerator()*mc;
    
    ai*=x.denominator();
    bi*=x.denominator();

    //reverse the orientation of the line if we want the halfplane under it
    if (overunder == 1) {
        ai*=-1; bi*=-1; ci*=-1;
    }

    //construct nef polygon representing this halfplane
    Line l(ai,bi,ci);
    Nef_polygon N(l, Nef_polygon::INCLUDED);

    return N;
}

Nef_polygon rotated_ell(struct slope myslope, struct interval xb, struct interval yb, int dim) {
    // construct the union of unit-width ells with given slope and inner corner at all coordinates
    // (along corridor wing directions) in the cartesian product of the intervals xb and yb
    //
    // if dim = 2, return union of such ells
    // if dim = 1, return union of the boundaries of such ells
    // if dim = 0, return union of the vertices of such ells
    //
    // if dim = 3, return version of (dim = 1) with only left side
    // if dim = 4, return version of (dim = 1) with only right side
    
    Nef_polygon N1, N2, N3, N4;

    switch(dim) {
        case 1:
        N1 = ell_side(myslope,xb,0,0,0,0).intersection(ell_side(myslope,xb,0,0,1,1)).intersection(ell_side(myslope,yb,1,0,1,1)); //left arm bottom
        N2 = ell_side(myslope,xb,0,1,0,0).intersection(ell_side(myslope,xb,0,1,1,1)).intersection(ell_side(myslope,yb,1,1,1,1)); //left arm top
        N3 = ell_side(myslope,yb,1,0,0,0).intersection(ell_side(myslope,yb,1,0,1,1)).intersection(ell_side(myslope,xb,0,0,1,1)); //right arm bottom
        N4 = ell_side(myslope,yb,1,1,0,0).intersection(ell_side(myslope,yb,1,1,1,1)).intersection(ell_side(myslope,xb,0,1,1,1)); //right arm top
        return N1.join(N2).join(N3).join(N4);

        case 0:
        N1 = ell_side(myslope,xb,0,0,0,0).intersection(ell_side(myslope,xb,0,0,1,1)).intersection(ell_side(myslope,yb,1,0,0,0)).intersection(ell_side(myslope,yb,1,0,1,1)); //bottom
        N2 = ell_side(myslope,xb,0,1,0,0).intersection(ell_side(myslope,xb,0,1,1,1)).intersection(ell_side(myslope,yb,1,1,0,0)).intersection(ell_side(myslope,yb,1,1,1,1)); //top
        return N1.join(N2);

        case 3:
        N1 = ell_side(myslope,xb,0,0,0,0).intersection(ell_side(myslope,xb,0,0,1,1)).intersection(ell_side(myslope,yb,1,0,1,1)); //left arm bottom
        N2 = ell_side(myslope,xb,0,1,0,0).intersection(ell_side(myslope,xb,0,1,1,1)).intersection(ell_side(myslope,yb,1,1,1,1)); //left arm top
        return N1.join(N2);

        case 4:
        N3 = ell_side(myslope,yb,1,0,0,0).intersection(ell_side(myslope,yb,1,0,1,1)).intersection(ell_side(myslope,xb,0,0,1,1)); //right arm bottom
        N4 = ell_side(myslope,yb,1,1,0,0).intersection(ell_side(myslope,yb,1,1,1,1)).intersection(ell_side(myslope,xb,0,1,1,1)); //right arm top
        return N3.join(N4);


	default:
        N1 = ell_side(myslope,xb,0,0,0,0).intersection(ell_side(myslope,xb,0,1,1,1)).intersection(ell_side(myslope,yb,1,1,1,1)); //left arm
        N2 = ell_side(myslope,yb,1,0,0,0).intersection(ell_side(myslope,yb,1,1,1,1)).intersection(ell_side(myslope,xb,0,1,1,1)); //right arm
        return N1.join(N2);

    }

}

Nef_polygon the_union(struct box mybox, struct bb_thread_params slopes){
    //return the intersection H^U_1^U_2^...^U_k^B
    //where U_i = union of L_i(u,v) as (u,v) varies over I_2i-1,I_2i
    //mybox = I_1 X I_2 X ... X I_2k
    unsigned int i;

    struct interval z;
    z.left=0; z.right=0;
    
    //First construct H (we allow a different initial angle than zero)
    Nef_polygon Nx = ell_side(slopes.initial,z,0,0,0,0);
    Nx = Nx.intersection(ell_side(slopes.initial,z,0,1,1,1));
    //If beta<pi/2, construct B
    if (slopes.has_final) {
        Nef_polygon Nfinal = ell_side(slopes.final_max,z,1,0,0,0).intersection(ell_side(slopes.final_min,z,1,1,1,1));
        Nfinal = Nfinal.join(ell_side(slopes.final_min,z,1,0,0,0).intersection(ell_side(slopes.final_max,z,1,1,1,1)));
        Nx = Nx.intersection(Nfinal);
    }

    for (i=0;i<slopes.num_intermediate;i++) {
	//construct U_i
        Nx = Nx.intersection(rotated_ell(slopes.intermediate[i],mybox.coord_bound_intervals[2*i],mybox.coord_bound_intervals[2*i+1],2));
    }

    return Nx;
}

void polygon_of_union(struct box mybox, struct bb_thread_params slopes, std::vector<ExactRational> *polygon){
    //stores the vertex coordinates of the polygon constructed by "the_union" in the vector referred to by polygon
    area(the_union(mybox,slopes),polygon);
    return;
}

ExactRational area_of_union(struct box mybox, struct bb_thread_params slopes){
    //return the area of the polygon constructed by "the_union".
    //This is \mathcal{G} in the paper
    return area(the_union(mybox,slopes),NULL);
}

ExactRational lower_bound_on_max_in_box(struct box mybox, struct bb_thread_params slopes){
    //calculates the value of g in the reference point (midpoint) of the box
    unsigned int i;
    struct box midpoint;
    
    midpoint = mybox;
    for (i=0;i<2*slopes.num_intermediate;i++) {
        midpoint.coord_bound_intervals[i] = midpoint_of_interval(mybox.coord_bound_intervals[i]);
    }

    return area_of_union(midpoint,slopes);
}

unsigned int index_of_coordinate_to_split(struct box mybox, struct bb_thread_params slopes){
    //determines which coordinate to split the given box along
    unsigned int i,j;
    struct interval z;
    ExactRational a,b;
    z.left=0; z.right=0;
    Nef_polygon N1,Nx;

    //construct the polygon whose area is \mathcal{G}
    Nx = the_union(mybox,slopes);

    //For each coordinate, determine the difference between the union varying over it and the intersection
    //The coordinate given the largest intersection of this difference with Nx is the one to be split
    a = 0; j=0;
    for (i=0;i<slopes.num_intermediate;i++) {
        N1 = rotated_ell(slopes.intermediate[i],mybox.coord_bound_intervals[2*i],mybox.coord_bound_intervals[2*i+1],3);
        b = area(Nx.intersection(N1),NULL);
        if (b>a) {a = b; j = 2*i + 0;}

        N1 = rotated_ell(slopes.intermediate[i],mybox.coord_bound_intervals[2*i],mybox.coord_bound_intervals[2*i+1],4);
        b = area(Nx.intersection(N1),NULL);
        if (b>a) {a = b; j = 2*i + 1;}
    }

    return j;
}


