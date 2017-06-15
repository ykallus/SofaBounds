//
//  sofa-bounds.hpp
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

#include <CGAL/Exact_integer.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <thread>
#include <chrono>
#include <math.h>
#include <gmp.h>

typedef CGAL::Exact_integer ExactInteger;
typedef CGAL::Exact_rational ExactRational;
typedef CGAL::Filtered_extended_homogeneous<ExactInteger> Extended_kernel;
typedef CGAL::Nef_polyhedron_2<Extended_kernel> Nef_polygon;
typedef Nef_polygon::Point Point;
typedef Nef_polygon::Line  Line;
typedef Nef_polygon::Explorer Explorer;

struct slope {
    ExactInteger a; //rise
    ExactInteger b; //run
    ExactInteger c; //hypothenuse
};

struct interval {
    ExactRational left; //left endpoint
    ExactRational right; //right endpoint
};

struct box {
    std::vector<struct interval> coord_bound_intervals;
    int depth;
    ExactRational upper_bound_on_max_in_box;
};

struct bb_thread_params {
    //problem specification
    std::vector<struct slope> intermediate; //alpha_1, ..., alpha_k
    struct slope initial; //should always be zero for a_priori_bounds to be correct
    struct slope final_min; //beta_1
    struct slope final_max; //beta_2
    bool has_final; //should be true if beta_2 < pi/2
    unsigned int num_intermediate; //k

    //runtime variables
    ExactRational lower_bound; //best lower bound so far
    ExactRational upper_bound; //best upper bound so far
    std::vector<ExactRational> lower_bound_witness;
    std::vector<ExactRational> lower_bound_polygon;
    std::vector<ExactRational> upper_bound_polygon;
    bool stop_flag; //set to true to make the branch-and-bound loop stop
    bool report_flag; //set to true to make branch-and-bound loop output upper and lower bound polygons to the file fpoly_name
    bool is_thread_running; //true while running
    std::thread *mythread; //thread running the branch-and-bound calculation
    std::string fpoly_name;
    unsigned long iterations;
    std::chrono::high_resolution_clock::time_point t_start;
    std::chrono::high_resolution_clock::time_point t_stop;

    bool reporteverysec_on, reporteveryiter_on, reporteveryjump_on;
    std::chrono::high_resolution_clock::time_point reporteverysec_last;
    unsigned int reporteverysec_inc;
    unsigned long reporteveryiter_last,reporteveryiter_inc;
    double reporteveryjump_last,reporteveryjump_inc;
};

struct CompareBoxes {
    //this struct is defined to provide the priority queue template access to an operator comparing box structs
    bool operator() (const struct box &op1, const struct box &op2) const { return (op1.upper_bound_on_max_in_box < op2.upper_bound_on_max_in_box); }
};

//these functions are implemented in sofa-bounds.cpp
std::string rat_str(ExactRational x);
void split_interval(struct interval b, struct interval *b1, struct interval *b2);
struct interval midpoint_of_interval(struct interval b);
ExactRational area(Nef_polygon N, std::vector<ExactRational> *polygon);
inline Nef_polygon ell_side(struct slope myslope, struct interval myinterval, int updown, int inout, int lowhi, int overunder);
Nef_polygon rotated_ell(struct slope myslope, struct interval xb, struct interval yb, int dim);
Nef_polygon the_union(struct box mybox, struct bb_thread_params slopes);
void polygon_of_union(struct box mybox, struct bb_thread_params slopes, std::vector<ExactRational> *polygon);
ExactRational area_of_union(struct box mybox, struct bb_thread_params slopes);
ExactRational lower_bound_on_max_in_box(struct box mybox, struct bb_thread_params slopes);
unsigned int index_of_coordinate_to_split(struct box mybox, struct bb_thread_params slopes);
struct box a_priori_bounds(struct bb_thread_params slopes, unsigned int i);
ExactRational initial_lower_bound(struct slope m);

//this function is implemented in branch-and-bound.cpp
int branch_and_bound(struct bb_thread_params *the_bb_thread_params);

//this function is implemented in frontend.cpp
void short_inspect(struct bb_thread_params *my_bb_thread_params);
