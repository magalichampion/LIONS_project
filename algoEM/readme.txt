Version : 0.2
(may 15, 2015)


Unpack the archive containing the jar files in a directory of your choice.
The deregulation detection tool is ddt.jar, make it executable or give it as an argument to the java utility with the option -jar.

Required parameters :
-network <path>  : A sample network is given in this archive
-expression <path> : A sample expression file is given in this archive
-scores <path> : where the posterior deregulation probabilities will be written

Optional parameters :
-convergence <amount> : default is 0.0001. Iteration stops when epsilon varies by less than this amount during one EM iteration.
-iterations <n> : even if the abov convergence criterion is not reached, stop after <n> iterations. Default is 1000.
-bp <n> : number of belief propagation passes done in step E. Default is 10.
-plot <path> : where to output the evolution of the parameters with each iteration.
-log <path> : where to output log messages about the program execution instead of standard output.

Choosing the EM initial parameter set. This is usually not necessary, the only recommended use is to set appropriate initial means if the expression values are not centered on zero. Give three real numbers separated by colons. e.g. -imu -4:-2:0.5
-imu <Minus:Zero:Plus>
-isigma <Minus:Zero:Plus>
-ialpha <Minus:Zero:Plus>
-iepsilon <value>

