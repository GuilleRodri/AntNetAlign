--- COMPILATION ---

A Makefile is provided. 
IF you plan to compile the Negative Learning version (AntNetAlign_Neg), the CPLEXDIR and CONCERTDIR paths need to be substituted by the directory where your CPLEX program is.

LINUX DISTRIBUTIONS:
To compile both versions, type 'make' into the command line. 
If you only want to compile one of the versions, type 'make AntNetAlign' for the base algorithm and 'make AntNetAlign_Neg' for the Negative Learning version.


--- PARAMETERS ---

-g1 [file]: Source network file in edgelist format. [MANDATORY]
-g2 [file]: Target network file in edgelist format. [MANDATORY]
-similarity [file]:  A similarity file with an V1 x V2 matrix. If not provided, the topological similarity from NETAL will be used instead.
-max_constructions [integer]: Maximum number of constructed solutions. Positive integer. Default: 1000.
-n_ants [integer]: Number of solution constructions per iteration (i.e., number of ants). Positive integer. Default: 10.
-pos_l_rate [float]: Positive learning rate. Float. 0 < l_rate <= 0.5. Default: 0.3.
-d_rate_select [float]: Determinism rate for selecting the next node to align. Float. 0 < d_rate_select <= 1. Default: 0.8.
-d_rate_align [float]: Determinism rate for candidate selection. Float. 0 < d_rate_align<= 0.9. Default: 0.9.
-score [S3/EC/ICS]: Score to maximize. Default: S3.
-max_t [float]: Maximum running time in seconds. Float. Default: 3600.
-tuning: If this flag is used, only the obtained final score will be displayed on the standard output. Used for tuning using the irace tool.
-save: If this flag is used, the obtained results will be saved in files.
-output_file [name]: Name that will be used for the generated output files (only if the -save flag is used). Default: "output".
-output_dir [path]: Path of the directory where the output files will be stored (only if the -save flag is used). Default: ".".

EXTRA PARAMETERS FOR THE NEGATIVE LEARNING VERSION:
-neg_l_rate [float]: Negative learning rate. Float. 0 < l_rate <= 0.5. Default: 0.1.
-t_sub [float]: Maximum time for subinstance solving in seconds. Float. Default: 60.
-emphasis [0/1/2/3/4]: MIP emphasis parameter of CPLEX. Default: 1.

-- TEST ---
Two sample networks (B_taurus.el and O_sativa.el) are provided, along with an example of a similarity file (B_taurus-O_sativa.similarity)

Command line example:
./AntNetAlign -g1 ./Test_Files/B_taurus.el -g2 ./Test_Files/O_sativa.el -similarity ./Test_Files/B_taurus-O_sativa.similarity


