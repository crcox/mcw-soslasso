function soslasso_sim_help()
% SOSLASSO_SIM_HELP Help with running simulations.
%
% See also:
% SOSLASSO_SIM_HELP>VARIABLES    Key to all variables.
% SOSLASSO_SIM_HELP>DATA_TYPES   Descriptions of different data types.
    help soslasso_sim_help
end
function variables()
    % KEY TO VARIABLES:
    %
    %                       * PARAMETERS *
    % NSUBJECTS The number of unique matrices that will be generated and
    % considered during optimization.  The relation of these matrices to
    % one another is controlled by the "data type". Different subjects will
    % always have different random, I.I.D. Gaussian noise masking their
    % signal, but the signal itself may be either exactly the same,
    % similar, or arbitrarily different in different subjects.  See
    % "soslasso_sim_help>data_types".
    %
    % NVOXELS Each subject's matrix has two dimensions: trials and voxels.
    % Voxels are the columns. The optimization routine will discover a set
    % of weights that are applied to the columns of the matrices to predict
    % the category of each trial. This parameter just defines how many
    % voxels there are in the whole brain. We assume all subjects have
    % identically sized brains that are perfectly aligned.
    %
    % NTRIALS See NVOXELS. Trials conceptually correspond to stimulus
    % presentations. Each row in the matrix is thus the simulated "evoked
    % response" to a given stimulus. This parameter changes how many trials
    % there are. The code is written to split this total number of trials
    % and create two categories. Conceptually:
    %    A = 1:floor(NTRIALS/2);
    %    B = ceil(NTRIALS/2):NTRIALS;
    %
    % NACTVOXELS The number of voxels that will be active on each trial.
    % Each trial, conceptually, evokes some pattern of activity that is
    % represented in the rows of each subject's matrix.  This parameter
    % determines exactly how many voxels will participate in any given
    % representation. The code is written so that every trail has exactly
    % the same number of voxels that are truly part of the representation.
    % The sparsity of this true representation can affect different methods
    % differently. See "soslasso_sim_help>data_types" for more information
    % on how these active voxels can be distributed.
    %
    % GROUPSIZE It may be the case that voxels that are closer to each
    % other in space will participate in representing similar things.  To
    % implement this, the 1-D voxel space is divided into groups, or
    % uniform size, and which may overlap. Depending on the "data type"
    % this will play an important role in defining how signal is
    % distributed. SET A GROUP SIZE NO MATTER WHAT. Even if a datatype does
    % not conceptually involve groups, setting a group size will make sure
    % that the number of voxels over which signal can possibly be
    % represented is the same across simulations experimenting with with
    % different types.
    %
    % GROUPSHIFT This defines the amount of overlap between groups, in
    % terms of the distance between the first voxel index of Group_n and
    % Group_n+1 (as if you were "shifting" a moving-window of size
    % "GROUPSIZE" along the voxel indexes). In principle, GROUPSHIFT >
    % GROUPSIZE is possible, but some voxels will not belong to any group
    % (which is probably not desirable). SET A GROUP SHIFT NO MATTER WHAT.
    %
    % NACTGROUPS: Similar to NACTVOXELS, but at the level of groups. When
    % deciding where to distribute the true representations, first a set of
    % active groups are chosen, and then a set of active voxels are
    % selected from the subset of voxels that belong to the active groups.
    % If NACTGROUPS is small and NACTVOXELS is large, the representations
    % will be densely localized, and vice versa (if the "data type"
    % considers groups at all).
    %
    % DATATYPEIND Generating data of different types is what makes these
    % simulations interesting.  DATATYPEIND is a simple index to select
    % which datatype to use. See "soslasso_sim_help>data_types".
    %    1 Same Sparse Groups
    %    2 Shifted Sparse Groups
    %    3 Different Sparse Groups
    %    4 Identical No Groups
    %    5 Different No Groups
    %
    %                      * VARIABLES *
    % V: Number of (unique) voxels across all selected groups. When
    % constructing simulated data, first a set of active groups are
    % selected (even if they will not be used), and V is computed. Then,
    % voxels are selected to compose the representation---but only within
    % from a subset of the total voxels, constrained by V in all cases.
    % Conceptually, this controls for the proportion of the brain that
    % contributes to "conceptual knowledge representation", while allowing
    % us to manipulate how this proportion of the brain is organized.
    % Otherwise, there would be no control over the proportion of the brain
    % utilized in data types that do not respect groups.
    %
    % X_truth: Data is created, at a high gloss, in two steps. The first
    % involves embedding the true signal. This involves inserting a 1 at
    % every voxel that is responding on a particular trial, for all trials.
    % The second adds in I.I.D. Gaussian noise to every point in the
    % matrix. X_truth is the data before adding noise.
    %
    % X: The data after adding noise (see X_truth).
    %
    % To be continued...
    
    help soslasso_sim_help>variables
end
function data_types()
    % Explanation of data types:
    %
    % SAME SPARSE GROUPS: A set of active groups is defined one time,
    % and used when constructing the data for all subject.  These means
    % that the same subset of voxels will contain the true representations
    % in all subjects, but the particular representations will involve
    % different permutations of active voxels within that subset.
    %    Expectation: SOSLASSO > UNIVARIATE > LASSO
    %
    % SHIFTED SPARSE GROUPS: Same as above, except that the group indexes
    % used for each subject are slightly jittered around the set of active
    % groups defined at the outset.  This means that the subset of voxels
    % utilized in each subject will differ to some degree, but will overlap
    % quite a bit.  There will be similar spatial structure in all
    % subjects, but not identical spatial structure.
    %    Expectation: SOSLASSO > UNIVARIATE > LASSO
    %
    % DIFFERENT SPARSE GROUPS: A more extreme perturbation---the set of
    % active groups is chosen anew for each subject. So, there is group
    % structure within each subject, but the groups are unlikely to
    % correspond across subjects.  
    %    Expectation: LASSO > SOSLASSO ~= UNIVARIATE
    %
    % IDENTICAL NO GROUPS: Active voxels are selected from exactly the same
    % subset of voxels in all subjects, but there is no spatial group
    % structure.  
    %    Expectation: UNIVARIATE > LASSO > SOSLASSO
    %
    % DIFFERENT NO GROUPS: Chaos. There is no spatial group structure, and
    % no correspondance across subjects. 
    %    Expectation: UNIVARIATE ~= LASSO ~= SOSLASSO == FAILURE

end