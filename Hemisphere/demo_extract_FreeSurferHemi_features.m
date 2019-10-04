%% DEMO FILE
% We show how to load controls and patients into separate tables and then
% merge both tables and join them with a corresponding metadata file.
% Note that in general there is no need to do it separately, but here we do
% to demonstrate the merging process.

%The lib folder and subfolders need to be added to path to work. Note that
%we included the freesurfer matlab functions again in the lib for
%completeness. Please refer to the freesurfer manual for further
%information on this.

%You can see some example outputs in the subfolder demo (as .mat files 
%containing tables).

addpath(genpath('../lib'))

%% load controls
% subjdir contains the FreeSurfer files, one directory per subject
subjdir='/data/yujiang/data/NeuroImaging/TLEall/';

% IDs are in the form 'C001', 'C002', ...
% you have to select which subjects you want to calculate by string IDs!
ids = num2str([1:30]', '%03u'); % formating ID (leading zeros)
ids = string([repmat('C', [length(ids),1]), ids]);

% pattern for foldername of each subject
format='*_T1.nii.gz'; % '*' will be replaced by the ID - can leave empty if you directly provide the directory names

% How to handle hemispheres? 'avg' for averaging, 'sum' for summing
% (caution: might make sense for area/volume but not for e.g. thickness),
% 'left' or 'right' for including only one of them respectively,
% or 'both' for having separate rows for left and right hemispheres.
% The column 'Hemisphere' indicates 'avg', 'sum' or 'left'/'right'.
hemi = "both"; %recommended

% run extraction - this is the main function doing the data extraction!
controls = extract_FreeSurferHemi_features(subjdir, ids, [], 'format', format, 'hemi', hemi);

%% load patients and merge with controls
ids = ["930"];

% merging is done by specifying the old table to add to/replacing old by new ids
[tbl, corrupt] = extract_FreeSurferHemi_features(subjdir, ids, tbl, ...
                                            'format', format, 'hemi', hemi, ...
                                            'replace', true, ...
                                            'verbose', true, ...
                                            'libdir', '../lib' ...
                                           ...'saveas', 'tle_controls'
                                            );
%% join with metadata

% load metadata table
load('./demo/metadata.mat'); % will load variable 'metadata'

% join table with metadata table by SubjectID (common column in the tables)
tle_controls = join(tle_controls, metadata, 'Keys', 'SubjectID');