function [simparams,X_truth,Y,ActiveVoxels] = soslasso_sim_setup(varargin)
%  SOSLASSO_SIM_SETUP Setup a simulation experiment and create a data
%  template.
%
%  USAGE:
%    simparams = SOSLASSO_SIM_SETUP() Returns default simparams structure, 
%    which can be modified and passed back into the function.
%
%    [simparams, X_truth, G] = SOSLASSO_SIM_SETUP(simparams) Using the information 
%    in simparams, generate the data template and group assignments, and log 
%    additional information in simparams.
%    
%    [...] = SOSLASSO_SIM_SETUP(simparams,'verbose') Same as above, but 
%    print plots and other useful information to the screen.
%
%    Key to Data Types: (N.B. Groups can be ``Modular'' by logical flag)
%    1 Same Sparse Groups
%    2 Shifted Sparse Groups
%    3 Different Sparse Groups
%    4 Identical No Groups
%    5 Different No Groups
%
%    If Modular==true, then active members within groups will all be 
%    responsive to the same category.
%
%    See also:
%    SOSLASSO_SIM_HELP>VARIABLES

	if nargin == 0
		simparams.nsubjects = uint32(16);	% number of subjects
		simparams.nvoxels = uint32(1024);   % number of voxels
		simparams.groupsize = uint32(64);	% group size
		simparams.groupshift = uint32(32);  % group shift (i.e. distance from G1(1)
											%     to G2(1) in voxel space.) 
		simparams.ntrials = uint32(64);     % number of trials
		simparams.nactgroups = uint32(4);   % number of active groups
		simparams.nactvoxels = uint32(8);   % number of active voxels per trial
		simparams.DataTypeInd = uint32(1);  % See help.
		simparams.Modular = false;          % See help.
		return
	end
	% N.B. the group size and group shift define the ground truth---how the
	% underlying signal is actually going to be created. The group size and shift
	% considered by SOSlasso can be defined separately. These are separate because
	% in practice the true groups are probably unknown, and we would like to know
	% how robust SOSlasso is to searching over groups that don't exactly correspond
	% to the truth.
	if nargin > 1
		VERBOSE = true;
	else
		VERBOSE = false;
	end

	simparams = varargin{1};
	if simparams.DataTypeInd == 1 
		simparams.DataType = 'Same Sparse Groups';
	elseif simparams.DataTypeInd == 2 
		simparams.DataType = 'Shifted Sparse Groups';
	elseif simparams.DataTypeInd == 3 
		simparams.DataType = 'Different Sparse Groups';
	elseif simparams.DataTypeInd == 4 
		simparams.DataType = 'Identical No Groups';
	elseif simparams.DataTypeInd == 5 
		simparams.DataType = 'Different No Groups';
	end

	%% Define Data
	[X_truth,Y] = define_data(simparams); % private function
	if VERBOSE
        for i=1:6
            subplot(2,3,i);
            imagesc(X_truth{i})
        end
	end
    
    %% Identify Active Voxels by Subject
    ActiveVoxels = cell2mat(cellfun(@any,X_truth,'Unif',0));
end


%% SUB-FUNCTIONS
function [X_truth,Y] = define_data(simparams)
% DEFINE_DATA Sub-function to soslasso_sim_setup
%
% See also:
% SOSLASSO_SIM_HELP>DATA_TYPES

	X_truth = cell(simparams.nsubjects,1);
    [X_truth{:}] = deal(zeros(simparams.ntrials,simparams.nvoxels));
    ani_trial_end = idivide(simparams.ntrials,2,'floor');
    art_trial_beg = ani_trial_end + 1;
    
    %% Create target indexes
	Y = cell(simparams.nsubjects,1);
	y = [ones(idivide(simparams.ntrials,2,'floor'),1);-ones(idivide(simparams.ntrials,2,'ceil'),1)];
	Y(:) = deal({y});
    
    a = 1:simparams.groupshift:(simparams.nvoxels-simparams.groupsize+1); % group start ind
    b = (simparams.groupsize):simparams.groupshift:simparams.nvoxels; % group end ind
    ngroups = length(a);
    G = cell(ngroups,1);
    for i=1:ngroups
        G{i} = a(i):b(i);
    end
    clear a b;
    
	switch simparams.DataType 
		case 'Same Sparse Groups'
			active_groups = uint32(randperm(ngroups,simparams.nactgroups));
		
            if simparams.Modular
				ani_groups = active_groups(1:idivide(simparams.nactgroups,2,'floor'));
				art_groups = active_groups(idivide(simparams.nactgroups,2,'ceil'):end);
				v_ani = unique(cell2mat(G(ani_groups))); % N.B. it is possible for ani and art
				v_art = unique(cell2mat(G(art_groups))); % voxels to overlap, whenever groups  
														 % overlap.
			else
				v = unique(cell2mat(G(active_groups)));
				N = uint32(length(v));
				v_randperm = v(randperm(N));
				v_ani = v_randperm(1:idivide(N,2,'floor')); % N.B. ani and art
				v_art = v_randperm(idivide(N,2,'ceil'):N);  % voxels will never overlap 
            end
		
			N_ani = uint32(length(v_ani));
			N_art = uint32(length(v_art));
		
            for i = 1:simparams.nsubjects
				% Animal Loop (first half of trials)
                for j=1:ani_trial_end;
					temp = uint32(randperm(N_ani,simparams.nactvoxels));
					active_voxels = v_ani(temp);
					X_truth{i}(j,active_voxels) = 1;
                end
				% Artifact Loop (second half of trials)
                for j=art_trial_beg:simparams.ntrials;
					temp = uint32(randperm(N_art,simparams.nactvoxels));
					active_voxels = v_art(temp);
					X_truth{i}(j,active_voxels) = 1;
                end
            end
            
		case 'Shifted Sparse Groups' 
			active_groups = uint32(randperm(ngroups,simparams.nactgroups));
		
            for i = 1:simparams.nsubjects
				step_size = 2;
				max_steps = idivide(simparams.groupsize,simparams.groupshift*step_size,'floor');
				n_steps = randi(max_steps);
				direction = sign(randn(1));
				shift = n_steps * step_size * direction;
				active_groups = mod(active_groups + shift,ngroups) + 1;
					 
                if simparams.Modular
					ani_groups_end  = idivide(simparams.nactgroups,2,'floor');
					art_groups_beg = idivide(simparams.nactgroups,2,'ceil');
					ani_groups = active_groups(1:ani_groups_end);
					art_groups = active_groups(art_groups_beg:end);
					v_ani = unique(cell2mat(G(ani_groups))); % N.B. it is possible for ani and art
					v_art = unique(cell2mat(G(art_groups))); % voxels to overlap, whenever groups  
															 % overlap.
                else
					v = unique(cell2mat(G(active_groups)));
					N = uint32(length(v));
					v_randperm = v(randperm(N));
					v_ani = v_randperm(1:idivide(N,2,'floor')); % N.B. ani and art
					v_art = v_randperm(idivide(N,2,'ceil'):N);  % voxels will never overlap 
                end
					
				N_ani = uint32(length(v_ani));
				N_art = uint32(length(v_art));

				% Make X_truth
       			% Animal Loop (first half of trials)
                for j=1:ani_trial_end;
					temp = uint32(randperm(N_ani,simparams.nactvoxels));
					active_voxels = g_ani(temp);
					X_truth{i}(j,active_voxels) = 1;
                end
				% Artifact Loop (second half of trials)
                for j=art_trial_beg:simparams.ntrials;
					temp = uint32(randperm(N_art,simparams.nactvoxels));
					active_voxels = g_art(temp);
					X_truth{i}(j,active_voxels) = 1;
                end
            end
 
		case 'Different Sparse Groups'
            for i = 1:simparams.nsubjects
				active_groups = uint32(randperm(ngroups,simparams.nactgroups));
					 
				if simparams.Modular
					ani_groups_end = idivide(simparams.nactgroups,2,'floor');
					art_groups_beg = idivide(simparams.nactgroups,2,'ceil');
					ani_groups = active_groups(1:ani_groups_end);
					art_groups = active_groups(art_groups_beg:end);
					v_ani = unique(cell2mat(G(ani_groups))); % N.B. it is possible for ani and art
					v_art = unique(cell2mat(G(art_groups))); % voxels to overlap, whenever groups  
															 % overlap.
				else
					v = unique(cell2mat(G(active_groups)));
					N = uint32(length(v));
					v_randperm = v(randperm(N));
					v_ani = v_randperm(1:idivide(N,2,'floor')); % N.B. ani and art
					v_art = v_randperm(idivide(N,2,'ceil'):N);  % voxels will never overlap 
				end

				N_ani = uint32(length(v_ani));
				N_art = uint32(length(v_art));
		
				% Make X_truth
				% Animal Loop (first half of trials)
                for j=1:ani_trial_end;
					temp = uint32(randperm(N_ani,simparams.nactvoxels));
					active_voxels = g_ani(temp);
					X_truth{i}(j,active_voxels) = 1;
                end
				% Artifact Loop (second half of trials)
                for j=art_trial_start:simparams.ntrials;
					temp = uint32(randperm(N_art,simparams.nactvoxels));
					active_voxels = g_art(temp);
					X_truth{i}(j,active_voxels) = 1;
                end
            end
	
		case 'Identical No Groups' % Cannot be modular
			active_groups = uint32(randperm(ngroups,simparams.nactgroups));
			v = unique(cell2mat(G(active_groups)));
			N = uint32(length(v));
			v_randperm = randperm(simparams.nvoxels,N);
			v_ani = v_randperm(1:idivide(N,2,'floor')); % N.B. ani and art
			v_art = v_randperm(idivide(N,2,'ceil'):N);  % voxels will never overlap 

			N_ani = uint32(length(v_ani));
			N_art = uint32(length(v_art));
            

            for i = 1:simparams.nsubjects
				% Make X_truth
				% Animal Loop (first half of trials)
                for j=1:ani_trial_end;
					temp = uint32(randperm(N_ani,simparams.nactvoxels));
					active_voxels = v_ani(temp);
					X_truth{i}(j,active_voxels) = 1;
                end
				% Artifact Loop (second half of trials)
                for j=art_trial_start:simparams.ntrials;
					temp = uint32(randperm(N_art,simparams.nactvoxels));
					active_voxels = v_art(temp);
					X_truth{i}(j,active_voxels) = 1;
                end
            end
	
		case 'Different No Groups'
			active_groups = uint32(randperm(ngroups,simparams.nactgroups));
			v = unique(cell2mat(G(active_groups)));
			N = uint32(length(v));
			v_randperm = randperm(simparams.nvoxels,N);
		
			N_ani = idivide(N,2,'floor');
			N_art = N - N_ani;

            for i = 1:simparams.nsubjects
				v_ani = v_randperm(1:idivide(N,2,'floor')); % N.B. ani and art
				v_art = v_randperm(idivide(N,2,'ceil'):N);  % voxels will never overlap 

				% Make X_truth
				X_truth{i} = zeros(simparams.ntrials,simparams.nvoxels);
				% Animal Loop (first half of trials)
                for j=1:ani_trial_end;
					temp = uint32(randperm(N_ani,simparams.nactvoxels));
					active_voxels = v_ani(temp);
					X_truth{i}(j,active_voxels) = 1;
                end
				% Artifact Loop (second half of trials)
                for j=art_trial_start:simparams.ntrials;
					temp = uint32(randperm(N_art,simparams.nactvoxels));
					active_voxels = v_art(temp);
					X_truth{i}(j,active_voxels) = 1;
                end
            end
	end
end
