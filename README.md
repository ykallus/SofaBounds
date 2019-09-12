# sofa-bounds
Calculates bounds on the sofa moving problem

The program depends on CGAL (The Computational Geometry Algorithms Library) and GMP (GNU Multiple Precision Arithmetic Library).
Tested to work with CGAL 4.12 and Boost 1.67.

On a typical setup, you will compile the program with the following command:

    g++ -o SofaBounds sofa-bounds.cpp branch-and-bound.cpp frontend.cpp -lCGAL -lgmp -lmpfr -O3 -frounding-math -finline-functions -lpthread

Running the executable will present you with a prompt.

    $ ./SofaBounds 
    SofaBounds version 1.0

    Type "help" for instructions.

    > 

To set up a branch-and-bound calculation, you need to set the number of intermediate angles, the slopes of all the intermediate angles, and the slopes of the minimum and maximum final angles, represented as a Pythagorean triple.


    > setcorr 3
    > setslope 1 7 24 25
    > setslope 2 33 56 65
    > setslope 3 119 120 169
    > setfinalmin 56 33 65
    > setfinalmax 24 7 25

When you issue the command "run", the program will start executing the branch-and-bound calculation in the background.

    > run
    <iterations=0> 

The prompt is now changed to give a brief report on the progress of the calculation. To keep track of the calculation you can retrigger the prompt, by hitting the carriage return, or issue the command "inspect" for a more detailed report.

    <iterations=0>
    <iterations=794 | upper bound=2.228 | time=0:00:59.616>
    <iterations=1179 | upper bound=2.179 | time=0:01:28.943> inspect
    lower bound  1.878819407167702 (101556597033993982997/54053410693201920000)
    upper bound  2.161169767609245 (272643940385/126155725696)
    iterations: 1422
    time: 0:01:47.520
    <iterations=1422 | upper bound=2.161 | time=0:01:47.520>

You can also tell the program to automatically trigger the brief reports every certain number of seconds, every certain number of iterations, or every time the upper bound has been improved by a certain amount.

    > reportevery 10 sec
    > reportevery 1000 iter
    > reportevery 0.1 jump

To cancel any of these automatic triggers, reissue the command with a parameter of 0.

    > reportevery 0 sec
    > reportevery 0 iter
    > reportevery 0 jump

You can stop a currently running calculation with the command "stop", which will also print a final report.

    <iterations=1789 | upper bound=2.136 | time=0:02:15.776> stop
    stopped after 1820 iterations
    lower bound  1.893117951427195 (6283362327082082057/3319054854635520000)
    upper bound  2.133771423128615 (3712203056657044111/1739737919638118400)
    iterations: 1820
    time: 0:02:18.202

You can print out the current settings with the command "settings", or save them to a file, together with the results of the calculation, if it has already started with the command "save [filename]". You can also load previously saved settings with "load [filename]". Finally, the command "savepoly [filename]" saves the polygons associated with the current upper and lower bounds to a file (this command only works if there is a thread currently running).
