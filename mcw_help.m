function mcw_help()
% MCW_HELP Index to help on using scripts/running experiments effectively.
%
% See also:
% MCW_HELP>HOW_TO_EXP Protocol for experimenting
% MCW_HELP>ON_SPACE   Notes on the different representational spaces
    help mcw_help
end

function how_to_exp()
% How to run experiments: 
% 1. Create a new informatively named directory.
% 2. Copy this script into the new directory.
% 3. Rename informatively.
% 4. Use the outline below to lay out the experiment's parameters and
% metadata.
% 5. After inspecting, save params, metadata, and GroupInfo to .mat files
% in the directory.
% 6. If params.Save == true, then mcw_prep_data() saves prepped versions of
% the X matrices into the directory. This is the data returned in cell 
% array X. [Highly recommended]. 
% 7. If there is anything noteworthy in the setup, create a README.txt in
% the directory.
% 8. Finally, execute classify_CRC() to run the experiment.
%
% Do not reuse this directory for other experiments. If a second experiment
% will use exactly the same X matrices, you might set params.Save = false,
% and create symbolic links to the original data (using ln -s
% ../[old_exp_dir]/*.[data_type].mat ./ from in the current experiment
% directory. This should also be noted in the README.
%
% Nikhil Rao and Chris Cox | University of Wisconsin-Madison | 2013-08-06
end

function on_space()
% ON_SPACE Musings on use of space
%
% NATIVE vs. TALAIRACH: This has to do with the XYZ coordinates
% assigned to each voxel in the space. XYZ coordinates in Talairach space
% have been warped so that and a particular XYZ coordinate in any subject
% refers to a similar physical point in the brain.  These coordinates are
% stored in XYZ_tlrc. Native space refers to the XYZ coordinates before
% warping to the common Talairach space.  In this space, a particular
% coordinate can refer to very different points in the brains of different
% people, due to a number of factors such as different brain size,
% different brain morphology, different positioning within the scanner,
% etc.  The XYZ coordinate system is in millimeters, with origin at the
% center of the brain.  Whether left or right hemisphere gets positive or
% negative numbers associated with it varies by standard. In the RAI
% (Right, Anterior, Inferior) standard, values increase as you move Right, 
% Anterior, or inferior (i.e. the smallest values will be in the left most,
% front most, and bottom most positions) RAI is also called DICOM or 
% "Radiologist perspective". Our data is RAI oriented.
%
% DIFFERENT INDEX SPACES: Most of the software, however, does not operate
% using XYZ information. Each XYZ coordinate maps to a particular column in
% the full data matrix. This full matrix is never interacted with, however.
% The full matrix would have columns for voxels which lie outside of the
% brain, or in our case, which don't belong to the cortex. So, we create an
% index relative to this submatrix of meaningful, cortical voxels that
% relates to the XYZ coordinates. 
%
%     However, for multitask learning/relating across subjects, we need to
% create matrices of the same size across all subjects.  The way we do this
% is to identify the min([x y z]) and max([x y z]) across all subjects, and
% define a new space---our "bounding box" for our dataset.  Representing
% each subject's data in terms of this full space is very inefficient---
% because this space is on a 1mm grid, and subjects functional data was 
% originally collected on ~3mm grids. So, there are many columns in these 
% full matrixes that are completely zero across subjects.
%
%     HOWEVER, we cannot simply drop these coordinates right away. They are
% informative for creating group assignments (see <a href="matlab: help formFMRIgroups">formFMRIgroups</a>). First 
% groups are formed uniformly. After the expanded group space matrices have
% been created (group space in both an index and Talairach sense... make
% sure you understand what this means) and group assignments for indexes in
% this space have been made, then the data matrices can be reduced and the
% group indexes can be reduced and remapped accordingly.
% 
%     By eliminating zero voxels, we change the "index space". The index
% space changes any time a column is added or removed, but there may be
% need to relate between these two.  These complications require careful
% consideration and documentation of which space is being worked with.
%
% See also:
% DEFINE_REP_SPACE
%
% Nikhil Rao and Chris Cox | University of Wisconsin-Madison | 2013-07-25
end