%% Instructions
% Copy this script to your running directory and run it before running
% kldivMain.m


%% Main input variables:

% Main data directory for perturbed and reference directories
settings.mydir = '/PATH/TO/TEST/DIR';
settings.refdir =  '/PATH/TO/REF/DIR';
settings.mainName = 'TEST';
settings.refName = 'REF';

% Path for the database used by md2path
settings.databasePath = '/PATH/TO/DATABASE';
settings.systemName = 'NAME'; % Points to the residue table in the database

% Name of the trajectory file in each run
settings.xtcName = 'traj.xtc';

% Number of frames to remove at the start of the simulation
settings.frames2skip = 500;

% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
settings.chains = 'A-B';

% Reference PDB code that this simulation is based on
settings.pdbCode = '3SN6';
settings.pdbChains = 'RA-';

% Inactive state reference PDB, used for dihedral comparison
settings.pdbCodeInactive = '3D4S';
settings.pdbInactiveChains = 'A--';
% Number of runs
settings.numRuns = 7;

% How many frames to skip when saving xtcs to dcds
settings.stride = 5;

settings.helices = [];

% KLDiv options

% number of blocks to divide the simulation for signficance testing,
% recommended number is the number of independent simulations or less.
% MUST BE an EVEN number
settings.nBlocks = 6; 

% Histogram finite size effect correction
settings.Grassberger = false ;

% Significance threshold for KLDiv, recommended: 0.1 for nBlocks = 4
% 0.05 for nBlocks >= 6
settings.st = 0.05;

% Specific residues to be highlighted during KL1 visualziation and the
% corresponding text
settings.highlightRes =  [ 1 2 3 ];
settings.highlightText = 'HIGHLIGHT';

% Activates GPCR specfic options
settings.isGPCR = true;
settings.includeLigPathwayCalc = false;