//
//  frontend.cpp
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


typedef struct {
    unsigned int a;
    unsigned int b;
    unsigned int c;
} PythagoreanTriple;

#define MAX_CORRIDORS       20
unsigned int num_corridors = 1;
PythagoreanTriple final_slope_min = { 1, 0, 1};
PythagoreanTriple final_slope_max = { 1, 0, 1};
PythagoreanTriple corridor_slopes[MAX_CORRIDORS] = { 119, 120, 169};
std::thread mythread;

const std::string welcome_string = "SofaBounds version 1.0\n\nType \"help\" for instructions.\n\n";
const std::string info_string = "SofaBounds version 1.0\n
Created by Yoav Kallus and Dan Romik\n
You are free to use and/or modify the software and source code. Enjoy!\n";
const std::string help_string = "\nSofaBounds version 1.0\n\n\
Valid commands are: \n\
help:\t\t\tprints this help message\n\
info:\t\t\tprints version and license information\n
settings:\t\tprints current list of settings\n\
quit:\t\t\tquit the program\n\
setcorr [num]:\tset the number of corridors to [num]\n\
setslope [ind] [a] [b] [c]:\tset the Pythagorean triple associated with the [ind]-th slope to ([a],[b],[c])\n\
setfinalmin [a] [b] [c]:\t\tset the Pythagorean triple associated with the minimum final slope to ([a],[b],[c])\n\
setfinalmax [a] [b] [c]:\t\tset the Pythagorean triple associated with the maximum final slope to ([a],[b],[c])\n\
reportevery [x] [sec|iter|jump]:\twhile running, print short progress report every x seconds, x iterations, or imporvement of x in the upper bound. A value of 0 cancels this reporting mode\n\
run:\t\t\trun thread with current settings\n\
stop:\t\t\tstop currently running thread\n\
inspect:\t\treport progress of current thread\n\
load [filename]:\t\tload settings from file\n\
save [filename]:\t\tsave settings and results to file\n\
savepoly [filename]:\t\tsave polygons associated with the current upper and lower bounds to file\n\
\n";


void inspect(struct bb_thread_params *my_bb_thread_params){
    //print progress report for currently running branch and bound calculation
    ExactRational x;
    std::chrono::high_resolution_clock::time_point t_end;
    if (my_bb_thread_params->is_thread_running) t_end = std::chrono::high_resolution_clock::now();
    else t_end = my_bb_thread_params->t_stop;
    x = my_bb_thread_params->lower_bound;
    std::cout << "lower bound " << std::fixed << std::setfill(' ') << std::fixed << std::setw(18) << std::setprecision(15) << CGAL::to_double(x) << " (" << rat_str(x) << ")" << std::endl;
    x = my_bb_thread_params->upper_bound;
    std::cout << "upper bound " << std::fixed << std::setfill(' ') << std::setw(18) << std::setprecision(15) << CGAL::to_double(x) << " (" << rat_str(x) << ")" << std::endl;
    std::cout << "iterations: " << my_bb_thread_params->iterations << std::endl;
    
    long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end-my_bb_thread_params->t_start).count();
    std::cout << "time: " << ms/(3600000L) <<":";
    std::cout << std::fixed << std::setfill('0') << std::setw(2) << ms%(3600000L)/60000L << ":";
    std::cout << std::fixed << std::setfill('0') << std::setw(2) << ms%(60000L)/1000L << ".";
    std::cout << std::fixed << std::setfill('0') << std::setw(3) << ms%1000L << std::endl;
} 

void short_inspect(struct bb_thread_params *my_bb_thread_params){
    std::chrono::high_resolution_clock::time_point t_end;
    if (my_bb_thread_params->is_thread_running) t_end = std::chrono::high_resolution_clock::now();
    else t_end = my_bb_thread_params->t_stop;
    long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end-my_bb_thread_params->t_start).count();
    std::cout << "<iterations=" << my_bb_thread_params->iterations;
    if (my_bb_thread_params->iterations == 0) {
	std::cout << "> ";
	return;
    }
    std::cout << " | upper bound=" << std::fixed << std::setprecision(3) << CGAL::to_double(my_bb_thread_params->upper_bound);
    std::cout << " | time="<< ms/(3600000L) <<":";
    std::cout << std::fixed << std::setfill('0') << std::setw(2) << ms%(3600000L)/60000L << ":";
    std::cout << std::fixed << std::setfill('0') << std::setw(2) << ms%(60000L)/1000L << "> "; //<< ".";
    //std::cout << std::fixed << std::setfill('0') << std::setw(3) << ms%1000L << "> ";
} 

std::thread run(struct bb_thread_params *my_bb_thread_params){
    //construct a bb_thread_params object from current settings and launch thread executing the branch-and-bound calculation
    std::thread mythread;

    my_bb_thread_params->iterations = 0;
    my_bb_thread_params->stop_flag = false;
    my_bb_thread_params->report_flag = false;
    my_bb_thread_params->num_intermediate = num_corridors;
    my_bb_thread_params->intermediate.resize(num_corridors);
    for (unsigned int i = 0; i < num_corridors; i++) {
	my_bb_thread_params->intermediate[i].a = (unsigned long) corridor_slopes[i].a;
	my_bb_thread_params->intermediate[i].b = (unsigned long) corridor_slopes[i].b;
	my_bb_thread_params->intermediate[i].c = (unsigned long) corridor_slopes[i].c;
    }
    my_bb_thread_params->final_min.a = (unsigned long) final_slope_min.a;
    my_bb_thread_params->final_min.b = (unsigned long) final_slope_min.b;
    my_bb_thread_params->final_min.c = (unsigned long) final_slope_min.c;
    my_bb_thread_params->final_max.a = (unsigned long) final_slope_max.a;
    my_bb_thread_params->final_max.b = (unsigned long) final_slope_max.b;
    my_bb_thread_params->final_max.c = (unsigned long) final_slope_max.c;
    my_bb_thread_params->initial.a = 0UL;
    my_bb_thread_params->initial.b = 1UL;
    my_bb_thread_params->initial.c = 1UL;
    if (final_slope_max.b > 0) {
	//if beta_2<pi/2, set flag which tells the relevant functions to construct the butterfly set and intersect with it
	my_bb_thread_params->has_final = true;
	//also, the initial lower bound 2.2 cannot be used, instead use a beta-dependent one (currently just 0 for all beta)
	my_bb_thread_params->lower_bound = initial_lower_bound(my_bb_thread_params->final_max);
    }
    else {
	//if beta_2=pi/2 H^B = H, so no need to construct the butterfly set
	my_bb_thread_params->has_final = false;
	//2.2 is a valid lower bound on G >= mu(pi/2) >= mu_gerver > 2.2
	my_bb_thread_params->lower_bound = ExactRational(11,5);
    }
    //record start time in order to later report elapsed time
    my_bb_thread_params->t_start = std::chrono::high_resolution_clock::now();
    my_bb_thread_params->reporteverysec_last = std::chrono::high_resolution_clock::now();;
    my_bb_thread_params->reporteveryiter_last = 0;
    my_bb_thread_params->reporteveryjump_last = 1e6;
    //launch thread
    mythread = std::thread(branch_and_bound, my_bb_thread_params);

    return mythread;
}

int process_command(struct bb_thread_params *my_bb_thread_params, std::string input_line) {
    std::string input_command = "";
    std::thread timerthread;
    ExactRational x;

    std::stringstream input_stream(input_line);
    input_stream >> input_command;

    if (input_command == "quit") {
	std::cout << "Quitting... ";
	if (my_bb_thread_params->is_thread_running) {
	    my_bb_thread_params->stop_flag = true;
	    my_bb_thread_params->mythread->join();
	    my_bb_thread_params->t_stop = std::chrono::high_resolution_clock::now();
	}
	std::cout << "Goodbye!\n\n";
	return 1;
    }
    else if (input_command == "help") {
	std::cout << help_string;
    }
    else if (input_command == "info") {
	std::cout << info_string;
    }
    else if (input_command == "settings") {
	std::cout << "\nNumber of corridors: " << num_corridors << "\n";
	std::cout << "\n";
	for (unsigned int i = 0; i < num_corridors; i++)
	    std::cout << "Slope " << i+1 << ":\t\t" << corridor_slopes[i].a << "\t" << corridor_slopes[i].b << "\t" << corridor_slopes[i].c << "\t\t(angle: " <<
		    180.0 / M_PI * asin((double)corridor_slopes[i].a / (double)corridor_slopes[i].c) << " deg)\n";
	std::cout << "Minimum final slope:\t" << final_slope_min.a << "\t" << final_slope_min.b << "\t" << final_slope_min.c << "\t\t(angle: " <<
		    180.0 / M_PI * asin((double)final_slope_min.a / (double)final_slope_min.c) << " deg)\n";
	std::cout << "Maximum final slope:\t" << final_slope_max.a << "\t" << final_slope_max.b << "\t" << final_slope_max.c << "\t\t(angle: " <<
		    180.0 / M_PI * asin((double)final_slope_max.a / (double)final_slope_max.c) << " deg)\n";
	if (my_bb_thread_params->reporteverysec_on || my_bb_thread_params->reporteveryiter_on || my_bb_thread_params->reporteveryjump_on)
	    std::cout << "\nReporting progress every:\t";
        if (my_bb_thread_params->reporteverysec_on) std::cout << "\t" << my_bb_thread_params->reporteverysec_inc << " seconds" << std::endl << "\t\t\t\t";
        if (my_bb_thread_params->reporteveryiter_on) std::cout << "\t" << my_bb_thread_params->reporteveryiter_inc << " iterations" << std::endl << "\t\t\t\t";
        if (my_bb_thread_params->reporteveryjump_on) std::cout << "\t" << my_bb_thread_params->reporteveryjump_inc << " decrease in upper bound" << std::endl << "\t\t\t\t";
	if (my_bb_thread_params->reporteverysec_on || my_bb_thread_params->reporteveryiter_on || my_bb_thread_params->reporteveryjump_on)
	    std::cout << std::endl;
    }
    else if (input_command == "setcorr") {
	unsigned int new_num_corridors;
	input_stream >> new_num_corridors;
	if (new_num_corridors > 0 && new_num_corridors <= MAX_CORRIDORS) {
	    num_corridors = new_num_corridors;
	}
	else {
	    std::cout << "Error: number of corridors must be between 1 and " << MAX_CORRIDORS << "\n";
	}
    }
    else if (input_command == "setslope") {
	unsigned int slope_index;
	PythagoreanTriple new_slope;
	input_stream >> slope_index >> new_slope.a >> new_slope.b >> new_slope.c;
	if (slope_index > 0 && slope_index <= num_corridors) {
	    if (new_slope.a*new_slope.a + new_slope.b*new_slope.b == new_slope.c*new_slope.c)
		corridor_slopes[slope_index-1] = new_slope;
	    else
		std::cout << "Error: the triple (" << new_slope.a << "," << new_slope.b << "," << new_slope.c << ") is not Pythagorean.\n";
	}
	else
	    std::cout << "Error: the index " << slope_index << " is not between 1 and " << num_corridors << "\n";
    }
    else if (input_command == "setfinalmin") {
	PythagoreanTriple new_slope;
	input_stream >> new_slope.a >> new_slope.b >> new_slope.c;
	if (new_slope.a*new_slope.a + new_slope.b*new_slope.b == new_slope.c*new_slope.c)
	    final_slope_min = new_slope;
	else
	    std::cout << "Error: the triple (" << new_slope.a << "," << new_slope.b << "," << new_slope.c << ") is not Pythagorean.\n";
    }
    else if (input_command == "setfinalmax") {
	PythagoreanTriple new_slope;
	input_stream >> new_slope.a >> new_slope.b >> new_slope.c;
	if (new_slope.a*new_slope.a + new_slope.b*new_slope.b == new_slope.c*new_slope.c)
	    final_slope_max = new_slope;
	else
	    std::cout << "Error: the triple (" << new_slope.a << "," << new_slope.b << "," << new_slope.c << ") is not Pythagorean.\n";
    }
    else if (input_command == "run") {
	if (!my_bb_thread_params->is_thread_running) {
	    mythread = run(my_bb_thread_params);
	    my_bb_thread_params->mythread = &mythread;
	    my_bb_thread_params->is_thread_running = true;
	} else {
	    std::cout << "Only one thread at a time.\n";
	}
    }
    else if (input_command == "stop") {
	if (my_bb_thread_params->is_thread_running) {
	    my_bb_thread_params->stop_flag = true;
	    my_bb_thread_params->mythread->join();
	    my_bb_thread_params->t_stop = std::chrono::high_resolution_clock::now();
	    inspect(my_bb_thread_params);
	    my_bb_thread_params->is_thread_running = false;
	} else {
	    std::cout << "No thread running.\n";
	}
    }
    else if (input_command == "inspect") {
	inspect(my_bb_thread_params);
    }
    else if (input_command == "witness") {
	for (unsigned int i = 0; i < 2*my_bb_thread_params->num_intermediate; i++){
	    x = my_bb_thread_params->lower_bound_witness[i];
	    std::cout << CGAL::to_double(x) << " (" << rat_str(x) << ")";
	    if (i%2 == 0) std::cout << ", ";
	    else std::cout << std::endl;
	}
    }
    else if (input_command == "load") {
	std::string filename, line;
	input_stream >> filename;
	std::ifstream infile(filename);
	bool success = 0;
	while (std::getline(infile, line)) {
	    process_command(my_bb_thread_params,line);
	    success = 1;
        }
	if (success)
	    std::cout << "File '" << filename << "' loaded successfully.\n";
	else
	    std::cout << "File '" << filename << "' could not be loaded.\n";
	}
    else if (input_command == "save") {
	std::string filename;
        input_stream >> filename;
	std::ofstream outfile(filename);

	outfile << "#" << std::endl;
	outfile << "# SofaBounds profile file" << std::endl;
	outfile << "#" << std::endl;
	outfile << "# Number of corridors: " << num_corridors << std::endl;
	outfile << "#" << std::endl;
	for (unsigned int i = 0; i < num_corridors; i++)
	    outfile << "# Slope " << i+1 << ":\t\t" << corridor_slopes[i].a << "\t" << corridor_slopes[i].b << "\t" << corridor_slopes[i].c << "\t\t(angle: " <<
		    180.0 / M_PI * asin((double)corridor_slopes[i].a / (double)corridor_slopes[i].c) << " deg)" << std::endl;
	outfile << "# Minimum final slope:\t" << final_slope_min.a << "\t" << final_slope_min.b << "\t" << final_slope_min.c << "\t\t(angle: " <<
		    180.0 / M_PI * asin((double)final_slope_min.a / (double)final_slope_min.c) << " deg)" << std::endl;
	outfile << "# Maximum final slope:\t" << final_slope_max.a << "\t" << final_slope_max.b << "\t" << final_slope_max.c << "\t\t(angle: " <<
		    180.0 / M_PI * asin((double)final_slope_max.a / (double)final_slope_max.c) << " deg)" << std::endl;
	outfile << std::endl;
	outfile << "setcorr " << num_corridors << std::endl;
	for (unsigned int i=0;i<num_corridors;i++)
	    outfile << "setslope " << i+1 << " " << corridor_slopes[i].a << " " << corridor_slopes[i].b << " " << corridor_slopes[i].c << std::endl;
	outfile << "setfinalmin " << final_slope_min.a << " " << final_slope_min.b << " " << final_slope_min.c << std::endl;
	outfile << "setfinalmax " << final_slope_max.a << " " << final_slope_max.b << " " << final_slope_max.c << std::endl;
	if (my_bb_thread_params->reporteverysec_on) outfile << "reportevery " << my_bb_thread_params->reporteverysec_inc << " sec" << std::endl;
	else outfile << "reportevery 0 sec" << std::endl;
	if (my_bb_thread_params->reporteveryiter_on) outfile << "reportevery " << my_bb_thread_params->reporteveryiter_inc << " iter" << std::endl;
	else outfile << "reportevery 0 iter" << std::endl;
	if (my_bb_thread_params->reporteveryjump_on) outfile << "reportevery " << my_bb_thread_params->reporteveryjump_inc << " jump" << std::endl;
	else outfile << "reportevery 0 jump" << std::endl;
	outfile << std::endl;

	if (my_bb_thread_params->iterations > 0) {
	    std::chrono::high_resolution_clock::time_point t_end;
	    if (my_bb_thread_params->is_thread_running) t_end = std::chrono::high_resolution_clock::now();
	    else t_end = my_bb_thread_params->t_stop;
	    outfile << "# iterations: " << my_bb_thread_params->iterations << std::endl;
	    long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end-my_bb_thread_params->t_start).count();
	    outfile << "# time: " << ms/(3600000L) <<":";
	    outfile << std::fixed << std::setfill('0') << std::setw(2) << ms%(3600000L)/60000L << ":";
	    outfile << std::fixed << std::setfill('0') << std::setw(2) << ms%(60000L)/1000L << ".";
	    outfile << std::fixed << std::setfill('0') << std::setw(3) << ms%1000L << std::endl;
	    outfile << "# lower bound: " << rat_str(my_bb_thread_params->lower_bound) << " = " << CGAL::to_double(my_bb_thread_params->lower_bound) << std::endl;
	    outfile << "# upper bound: " << rat_str(my_bb_thread_params->upper_bound) << " = " << CGAL::to_double(my_bb_thread_params->upper_bound) << std::endl;
	    outfile << "# lower bound witness arguments: " << std::endl;
	    outfile << "# ";
	    for (unsigned int i = 0; i<2*my_bb_thread_params->num_intermediate; i++){
		outfile << rat_str(my_bb_thread_params->lower_bound_witness[i]);
		if (i%2 == 0) outfile << ", ";
		else outfile << std::endl << "# ";
	    }
	    outfile << std::endl;
	}
	outfile.close();
    }
    else if (input_command == "savepoly") {
	std::string filename;
        input_stream >> my_bb_thread_params->fpoly_name;
	my_bb_thread_params->report_flag = true;
    }
    else if (input_command == "reportevery") {
	std::string str_inc,type;
	input_stream >> str_inc >> type;
	if (type == "sec") {
	    unsigned int inc_sec;
	    std::stringstream(str_inc) >> inc_sec;
	    if (inc_sec == 0) my_bb_thread_params->reporteverysec_on = false;
	    else {
		my_bb_thread_params->reporteverysec_inc = inc_sec;
		my_bb_thread_params->reporteverysec_on = true;
	    }
	} else if (type == "iter") {
	    unsigned long inc_iter;
	    std::stringstream(str_inc) >> inc_iter;
	    if (inc_iter == 0) my_bb_thread_params->reporteveryiter_on = false;
	    else {
		my_bb_thread_params->reporteveryiter_inc = inc_iter;
		my_bb_thread_params->reporteveryiter_on = true;
	    }
	} else if (type == "jump") {
	    double inc_jump;
	    std::stringstream(str_inc) >> inc_jump;
	    if (inc_jump == 0) my_bb_thread_params->reporteveryjump_on = false;
	    else {
		my_bb_thread_params->reporteveryjump_inc = inc_jump;
		my_bb_thread_params->reporteveryjump_on = true;
	    }
	} else {
	    std::cout << "Unknown increment type for reportevery" << std::endl;
	}
    }
    else if (input_command != "" && input_command[0] != '#' /* to allow inputting comments */ ) {
	std::cout << "Unknown command \"" << input_command << "\". Type \"help\" for instructions.\n";
    }
    return 0;
}

int main(int argc, const char * argv[]) {
    std::cout << welcome_string;
    std::string input_line;
    struct bb_thread_params my_bb_thread_params;

    my_bb_thread_params.iterations = 0;
    my_bb_thread_params.is_thread_running = false;
    my_bb_thread_params.reporteverysec_on = true;
    my_bb_thread_params.reporteverysec_inc = 10;
    my_bb_thread_params.reporteveryiter_on = false;
    my_bb_thread_params.reporteveryjump_on = false;
    while (true) {
        if (my_bb_thread_params.is_thread_running)
	    short_inspect(&my_bb_thread_params);
        else
	    std::cout << "> ";
        getline(std::cin, input_line); 
	if (std::cin.eof()) input_line = "quit\n";
	if (process_command(&my_bb_thread_params,input_line)==1) return 0;
    }

    return 0;
}
