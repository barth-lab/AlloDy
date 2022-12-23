%% Instructions
% Copy this script to your running directory and run it before running
% md2path.m

% Preparations needed: PDB file with separate chains for receptor, ligand,
% and G-protein: 1st chain is receptor, 2nd is ligand, 3rd is Gp
% Input: directory, number of runs to analyze, and stride of trajectory
% NOTE: all input residues should be in the same pdb numbering as the input
% PDB


%% Main input variables:

% Main data directory
settings.mydir = '/PATH/TO/DIR';

% Path for the database used by md2path
settings.databasePath = '/PATH/TO/DATABASE';
settings.systemName = 'NAME'; % Points to the residue table in the database

% Name of the trajectory file in each run
settings.xtcName = 'traj.xtc';

% Number of frames to remove at the start of the simulation
settings.frames2skip = 500;

% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
settings.chains = 'ACB';

% Reference PDB code that this simulation is based on
settings.pdbCode = '3SN6';
settings.pdbChains = 'RA-';

% Inactive state reference PDB, used for dihedral comparison
settings.pdbCodeInactive = '3D4S';
settings.pdbInactiveChains = 'A--';

% Extra PDBs for reference (in conformational landscape plots for example)
settings.pdbCodeExtra = {'7DHR','7F1R'};
settings.pdbChainsExtra = {'RA-','RA-'};

% Number of runs
settings.numRuns = 5;

% How many frames to skip when saving xtcs to dcds
settings.stride = 5;

% Actiavtes GPCR specfic options (calculating 7 TM helices and
% GPCR allosteric path calculation, as well as planned GPCR order parameters)
settings.isGPCR = true;

% Residues to be ignored in the pathway calculation, usually you want to
% ignore floppy termini or long loops that introduce noise into the system.
% Defaults to first and last 4 residues of the receptor chain. Assumes that
% protein-renum, true to its name, is renumbered starting from 1.
settings.excludedResidues = [];

% The residues belonging to the 7 helical turns of a GPCR, If
% empty, the code will attempt to find them by calling VMD and extracting a
% secondary structure variable.
settings.helices = [];

% If binding site/ligand interacting residues are not defined,
% the code will calculate contact maps of the protein and ligand.
settings.bsResidues = [];

% If intracellular binding partner is present, check the contacts from the
% original PDB structure and NOT from the simulation, since G-protein is
% usually truncated in simulations and some contacts are lost.
settings.gpiResidues = [];

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
settings.includeLigPathwayCalc = false;


% Allosteric pathway run options (Options controlled from md2path script now)

% Distance cutoff for considering "allosteric" signal, default is 10 A
settings.disCutoff = 10;

% Minimal fraction of MI to be conserved when filtering pathways
settings.miFractionCutoff = 0.85;

% Cutoff for pathways to be considered near to cluster them into one pipeline
settings.nearCutoff = 7.5;

% Amount of overlap required for pathways to be considered near
settings.overlapCutoff = 0.75;
