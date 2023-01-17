%% Instructions
% Copy this script to your running directory and run it before running
% md2path.m

%% Main input variables:

% Main data directory
settings.mydir = 'Path\To\demo';

% Path for the database used by md2path
settings.databasePath = 'Path\To\database';
settings.systemName = 'CXCR4'; % Points to the residue table in the database

% Name of the trajectory file in each run
settings.xtcName = 'traj.xtc';

% Number of frames to remove at the start of the simulation
settings.frames2skip = 20;

% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
settings.chains = 'A-B';

% Reference PDB code that this simulation is based on
settings.pdbCode = '4XT1';
settings.pdbChains = 'A-B';

% Inactive state reference PDB, used for dihedral comparison
settings.pdbCodeInactive = '4RWS';
settings.pdbInactiveChains = 'A--';

% Extra PDBs for reference (in conformational landscape plots for example)
settings.pdbCodeExtra = {'6LFO','7F1R'};
settings.pdbChainsExtra = {'RAD','RA-'};

% Number of runs
settings.numRuns = 3;

% How many frames to skip when saving xtcs to dcds
settings.stride = 20;

% Activates GPCR specfic options (calculating 7 TM helices and
% GPCR allosteric path calculation, as well as planned GPCR order parameters)
settings.isGPCR = true;

% Residues to be ignored in the pathway calculation, usually you want to
% ignore floppy termini or long loops that introduce noise into the system.
% Defaults to first and last 4 residues of the receptor chain. Assumes that
% protein-renum, true to its name, is renumbered starting from 1.
settings.excludedResidues = [     1     4
    275 277
   287   288];

% The residues belonging to the 7 helical turns of a GPCR, If
% empty, the code will attempt to find them by calling VMD and extracting a
% secondary structure variable.
settings.helices = [5 39
46 74
80 114
127 149
168 202
207 242
248 276];

% If binding site/ligand interacting residues are not defined,
% the code will calculate contact maps of the protein and ligand.
settings.bsResidues = [];

% If intracellular binding partner is present, check the contacts from the
% original PDB structure and NOT from the simulation, since G-protein is
% usually truncated in simulations and some contacts are lost.
settings.gpiResidues  = [61, 62, 63, 64, 65, 92, 93, 94, 95, 96, 177, 178, 179, 223, 224, 225, 226, 272]; 


% Same as gpi but for beta arrestin, keep empty unless you want to do a beta
% arrestin bias calculation.
settings.baiResidues = [];

% Do principle component analysis on the ligand binding pose as
% well as the receptor conformations
settings.doPCA = true;

% Number of clusters for PCA analysis of the ligand, if empty,
% the code will calculate it automatically, BUT it will take time!!!
settings.kClusters = [];

% Maximum number of clusters to consider during PCA clustering,
% only give this option if kClusters is empty!!!
settings.kmax = 15;

% Include ligand (peptide only for now) in pathway calculations
% NOTE: take extra care when using this option, as the pathway part of
% the code does not deal well with different chains, to avoid problems, do
% the following:
% 1- Make sure there are no residue numbers repeated in the peptide and
% receptor
% 2- Renumber both receptor and peptide in the same order in which they
% appear in the pdb file
settings.includeLigPathwayCalc = true;


% Allosteric pathway run options (Options controlled from md2path script now)

% Distance cutoff for considering "allosteric" signal, default is 10 A
settings.disCutoff = 10;

% Minimal fraction of MI to be conserved when filtering pathways
settings.miFractionCutoff = 0.85;

% Cutoff for pathways to be considered near to cluster them into one pipeline
settings.nearCutoff = 7.5;

% Amount of overlap required for pathways to be considered near
settings.overlapCutoff = 0.75;
