simparams = soslasso_sim_setup(1);
% modify simparams structures

[simparams, X] = soslasso_sim_setup(1, simparams);

simdata = soslasso_sim_makesimdata(X,simparams);