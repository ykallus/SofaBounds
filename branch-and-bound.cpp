//
//  branch-and-bound.cpp
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

#include <queue>
#include <vector>
#include "sofa-bounds.hpp"

int branch_and_bound(struct bb_thread_params *my_bb_thread_params) {
    unsigned int i,j;
    std::priority_queue<struct box, std::vector<struct box>, struct CompareBoxes> boxqueue;
    struct box currentbox;
    struct box newbox[2];
    ExactRational newbound[2];
    ExactRational a,b,c,d;
    ExactRational bestyet;

    //initialize lower bound
    bestyet = my_bb_thread_params->lower_bound;

    //initailize the queue
    newbox[0] = a_priori_bounds(*my_bb_thread_params,my_bb_thread_params->num_intermediate/2);
    newbox[0].upper_bound_on_max_in_box = area_of_union(newbox[0],*my_bb_thread_params);
    boxqueue.push(newbox[0]);
    
    //initialize interation counter and vector recording the lower bound witness
    my_bb_thread_params->iterations = 0UL;
    my_bb_thread_params->lower_bound_witness.resize(2*my_bb_thread_params->num_intermediate);

    //main branch-and-bound loop
    while (!boxqueue.empty() && !my_bb_thread_params->stop_flag) {
	my_bb_thread_params->iterations++;

        //pop top box from the stack, and declare its upper bound the universal upper bound
        currentbox = boxqueue.top();
        boxqueue.pop();
	my_bb_thread_params->upper_bound = currentbox.upper_bound_on_max_in_box;

	//if requested, construct polygons corresponding to the current best lower bound and upper bound
	//and print their vertex coordinates to a file
	if (my_bb_thread_params->report_flag) {
	    //upper bound polygon
	    //first line is number of vertices x 2
	    polygon_of_union(currentbox,*my_bb_thread_params,&my_bb_thread_params->upper_bound_polygon);
	    std::ofstream fpoly(my_bb_thread_params->fpoly_name);
	    fpoly << my_bb_thread_params->upper_bound_polygon.size() << std::endl;
	    for (i=0; i<my_bb_thread_params->upper_bound_polygon.size();i++)
		//each subsequent line alternates between x and y coordinates of successive vertices
		fpoly << rat_str(my_bb_thread_params->upper_bound_polygon[i]) << std::endl;
	    //lower bound polygon
	    fpoly << my_bb_thread_params->lower_bound_polygon.size() << std::endl;
	    for (i=0; i<my_bb_thread_params->lower_bound_polygon.size();i++)
		fpoly << rat_str(my_bb_thread_params->lower_bound_polygon[i]) << std::endl;
	    fpoly.close();
	    my_bb_thread_params->report_flag = false;
	}

        //if local lower bound on max in current box is better than current global bound, update global bound
        a = lower_bound_on_max_in_box(currentbox,*my_bb_thread_params);
        if (a > bestyet) {
	    bestyet = a;
	    my_bb_thread_params->lower_bound = bestyet;
	    newbox[0] = currentbox;
	    //also record the witness variables and the associated polygon
	    for (i=0;i<2*my_bb_thread_params->num_intermediate;i++) {
		newbox[0].coord_bound_intervals[i] = midpoint_of_interval(currentbox.coord_bound_intervals[i]);
		my_bb_thread_params->lower_bound_witness[i] = newbox[0].coord_bound_intervals[i].left;
	    }
	    polygon_of_union(newbox[0],*my_bb_thread_params,&my_bb_thread_params->lower_bound_polygon);
	}

	if (my_bb_thread_params->reporteverysec_on) {
	    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
	    //long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end-my_bb_thread_params->reporteverysec_last).count();
	    //if (ms >= 1000*my_bb_thread_params->reporteverysec_inc) {
	    if (t_end - my_bb_thread_params->reporteverysec_last >= std::chrono::seconds(my_bb_thread_params->reporteverysec_inc)) {
		std::cout << std::endl;
		short_inspect(my_bb_thread_params);
		std::cout << std::flush;
		//my_bb_thread_params->reporteverysec_last = t_end;
	        while (t_end - my_bb_thread_params->reporteverysec_last >= std::chrono::seconds(my_bb_thread_params->reporteverysec_inc)) my_bb_thread_params->reporteverysec_last += std::chrono::seconds(my_bb_thread_params->reporteverysec_inc);

	    }
	}
	if (my_bb_thread_params->reporteveryiter_on) {
	    if (my_bb_thread_params->iterations - my_bb_thread_params->reporteveryiter_last >= my_bb_thread_params->reporteveryiter_inc) {
		std::cout << std::endl;
		short_inspect(my_bb_thread_params);
		std::cout << std::flush;
		my_bb_thread_params->reporteveryiter_last = (my_bb_thread_params->iterations/my_bb_thread_params->reporteveryiter_inc)*my_bb_thread_params->reporteveryiter_inc;
	    }
	}
	if (my_bb_thread_params->reporteveryjump_on) {
	    if (my_bb_thread_params->reporteveryjump_last - CGAL::to_double(my_bb_thread_params->upper_bound) >= my_bb_thread_params->reporteveryjump_inc) {
		std::cout << std::endl;
		short_inspect(my_bb_thread_params);
		std::cout << std::flush;
		my_bb_thread_params->reporteveryjump_last = round(CGAL::to_double(my_bb_thread_params->upper_bound)/my_bb_thread_params->reporteveryjump_inc)*my_bb_thread_params->reporteveryjump_inc;
	    }
	}

        //uncomment for verbose output
        //fprintf(stdout,"%12lu (%6d, %12lu): %18.15lf %18.15lf %18.15lf\n",iter,currentbox.depth,boxqueue.size(),CGAL::to_double(a),CGAL::to_double(currentbox.upper_bound_on_max_in_box),CGAL::to_double(bestyet));
        
        //split the box into two descendents
        j = index_of_coordinate_to_split(currentbox,*my_bb_thread_params);
        newbox[0] = currentbox;
        newbox[1] = currentbox;
        split_interval(currentbox.coord_bound_intervals[j],&(newbox[0].coord_bound_intervals[j]),&(newbox[1].coord_bound_intervals[j]));
	//calculate the upper bound for the descendents
        for (i=0;i<2;i++) {
            newbound[i] = area_of_union(newbox[i],*my_bb_thread_params);
        }

        for (i=0;i<2;i++) {
	    //if descendent cannot be discarded
            if ( (newbound[i] < 0) || (newbound[i] > bestyet) ) {
		//update its depth and upper bound and push it to the queue
                newbox[i].depth++;
                newbox[i].upper_bound_on_max_in_box = newbound[i];
                boxqueue.push(newbox[i]);
            }
	    //if it can be discarded, just do nothing
        }

    }

    //this point should only be reach if the initial "lower bound" was greater than the actual supremum,
    //otherwise the loop should continue ad infinitum giving progressively better and better bounds.
    std::cerr << "stopped after " << my_bb_thread_params->iterations << " iterations" << std::endl;

    return 0;
}

