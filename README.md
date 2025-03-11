# compnet_ms_analyses
Analyses accompanying Rosenblad (2025) compnet: Bayesian dyadic regression reveals hidden signals of interspecific competition in observational data

There are five subdirectories that each correspond to a distinct analysis.

Within each subdirectory, scripts may have numbers at the beginnings of the file names. If they do, this means scripts with earlier numbers need to be run before those with later numbers. e.g., 1_myscript.R must be run before 2_myotherscript.R.

compnet_ms_example runs the toy example shown in Box 1 with compnet's built-in example data.

compnet_ms_fia runs the empirical case study of trees in Payette National Forest.

compnet_ms_metacom_sim_uncorr runs the metacommunity simulations, and accompanying data analyses, with uncorrelated species traits. This generates the spatiotemporal landscapes for the correlated scenario as well, so run this one first.

compnet_ms_metacom_sim_corr runs the metacommunity simulations, and accompanying data analyses, with correlated species traits.

compnet_ms_metacom_sim_summary produces the summary figure encompassing both simulation scenarios. The script in this directory must be run after the above two analyses have run.

The two metacommunity simulation analyses require a lot of computing resources. I recommend running them on a cluster. You may need to adjust the number of processor cores if you don't have enough RAM and/or enough cores (see comments near top of each script regarding how to do this).

Please contact me at kyle_rosenblad@berkeley.edu with any questions.