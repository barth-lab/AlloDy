% Features to add:
% Add ligand binding residues to RMSF plot
% Add Ligand atoms instead of atom numbers on RMSF ligand plot

%% This set of scripts/functions aims to take the output of an MD simulation
% and calculate allosteric pathways, going thru several steps of analysis

% Preparations needed: PDB file with separate chains for receptor, ligand,
% and G-protein: chains can be set up in the input_md2path script
% NOTE: all input residues should be in the same pdb numbering as the input
% PDB
% DO NOT FORGET TO SETUP AND RUN input_md2path before running this script

% Small roadmap:
% 1- transform xtc to dcd
% 2- load trajs
% 3- Calc RMSD and RMSF
% 4- Calc contact map with ligand and/or Gp (if they exist)
% 5- Calc GPCR order param for activation/state
% 6- Make a PCA of ligand binding poses and receptor coordinates
% 7- Calc dihedrals
% 8- Calc entropies
% 9- Evalutate convergence of entropies
% 10- Prep input for pathway calculation
% 11- Run pathway calculation


%% Check that the data directory exists

assert(logical(exist(settings.mydir, 'dir')), "Directory does not exist!!! Please point to an existing directory!");


%% Load, fetch and align PDB files

database = Database(settings.databasePath);

% Chains: [receptor, G protein, ligand]
database.read(fullfile(settings.mydir, "prot.pdb"), settings.chains, "Protein");

if ~isempty(settings.pdbCode)
    database.fetch(settings.pdbCode, settings.pdbChains);

    if ~isempty(settings.pdbCodeInactive)
        database.fetch(settings.pdbCodeInactive, settings.pdbInactiveChains);
    end
    for i = 1:length(settings.pdbCodeExtra)
        database.fetch(settings.pdbCodeExtra{i}, settings.pdbChainsExtra{i});
    end
end

% Align sequences and produce the database.residues table
database.align();

% Align structures to the first database entry
database.alignStructures();

% Store names for each chain
Chains.receptor = 1;
Chains.gprotein = 2;
Chains.ligand = 3;

% Keep references to entries and chains used frequently
mainEntry = database.entries{1};
mainChain = mainEntry.chains{Chains.receptor};

if mainEntry.hasChain(Chains.ligand)
    ligandChain = mainEntry.chains{Chains.ligand};
end

% Temp section until Mathworks decides to fix how copyfile works on Linux
% systems
if isunix
    slash = '/';
    copy = 'cp';
elseif ispc
    slash = '\';
    copy = 'copy';
end

%% Initialize the directory where output and figures will be saved

% TODO: Remove? Only used once

[~,name,~] = fileparts(settings.mydir); % gives path and current folder name
name(name=='_')='-';

if length(name) > 3 && name(2) == '-' % just my directory naming convention
    name = name(3:end);
elseif length(name) > 3 && name(3) == '-'
    name = name(4:end);
end

mainEntry.name = name;


md2pathdir = fullfile(settings.mydir, "md2path");
saveVarName = fullfile(md2pathdir, "workspace.mat"); % Where to save important variables

% Create md2pathdir if it does not exist
if ~exist(md2pathdir, 'dir')
    mkdir(md2pathdir);
    add2log(md2pathdir, "Creating md2path directory in " + settings.mydir);
end
if ~exist(saveVarName,'file')
    save(saveVarName,'md2pathdir');
else
    save(saveVarName,'md2pathdir','-append');
end

% Possibly copy some of the PDB to md2pathdir
% e.g. copyfile(database.entries{1}.path, fullfile(md2pathdir, "foo.pdb"));


%% Add labels

% Add labels with Ballesteros-Weinstein notation, aligned with the 2nd
% database entry
% Exported from GPCRdb (https://gpcrdb.org/residue/residuetable)
if length(database.entries) > 1 && isfield(settings, 'systemName')
    residueTablePath = fullfile(database.dir, settings.systemName + "_residue_table.xlsx");
    database.label(3, residueTablePath);
end

%% Find mutations compared to a downloaded reference sequence:

if exist(fullfile(settings.databasePath, settings.systemName + ".fasta"),'file') % Quick and dirty if statement, change later
    fastaPath = fullfile(settings.databasePath, settings.systemName + ".fasta");
    mutations = database.findMut(1, Chains.receptor, fastaPath);
    mutations
end
%% Transform .xtc files into .dcd

% Are the DCD files there? If yes, skip this step
% areThereDCDs = dir([settings.mydir '/run*/traj.dcd']);
areThereDCDs = dir(fullfile(settings.mydir, "run*", "traj.dcd"));

if length(areThereDCDs) < settings.numRuns
    % run VMD from command line:  vmd -dispdev text -e
    pathToScript = fullfile(pwd(), "load_save.tcl");  % assumes script is in current directory
    cmdStr       = "vmd -dispdev text -e " + pathToScript + " -args " + settings.mydir + " " + num2str(settings.numRuns) + " " + num2str(settings.stride) + " " + database.entries{1}.path + " " + settings.xtcName;
    system(cmdStr);
    add2log(md2pathdir, "Transformed .xtc files to .dcd files with the following command: """ + cmdStr + """");
end


%% Import trajectories, decenter and align them using CAs to first frame

mainSim = mainEntry.addSimulation(fullfile(settings.mydir, "run*"),'align2chain',settings.chains(Chains.receptor));


%% Calculate RMSD and RMSF of CAs and ligand

% Receptor:
% Calculate RMSD and RMSF of the CA atoms over the course of the sims

receptorAtoms = mainChain.getAtoms();
RMSD = mainSim.computeRmsd(receptorAtoms);
RMSF = mainSim.computeRmsf(receptorAtoms, 'StartFrame', settings.frames2skip + 1);
save(saveVarName,'RMSD','RMSF','-append');

bfactor = zeros(1, mainEntry.atomCount);
bfactor(1, receptorAtoms) = mean(RMSF, 2);

if mainEntry.hasChain(Chains.ligand)
    ligandAtoms = ligandChain.getLigandAtoms();
    RMSD_L = mainSim.computeRmsd(ligandAtoms);
    RMSF_L = mainSim.computeRmsf(ligandAtoms, 'StartFrame', settings.frames2skip + 1);
    bfactor(1, ligandAtoms) = mean(RMSF_L, 2);
    save(saveVarName,'RMSD_L','RMSF_L','-append');
end

% Save the RMSF as the B-factor
writepdb(fullfile(md2pathdir, "prot_rmsf.pdb"), mainEntry.pdb, [], 'default', bfactor);


% Options:
errStrideRMSF = 3; % how far are the error bars on the RMSF plot

plot_save_rmsd_rmsf;
% Plots RMSD/RMSF plots for receptor (and ligand) and saves them to md2pathdir


%% Calculate contact map of receptor with ligand
% If ligand/gp exists in simulation!

if isempty(settings.bsResidues) && mainEntry.hasChain(Chains.ligand)
    add2log(md2pathdir, 'Calculating contact map of receptor with ligand');
    % Dual cutoff contact scheme: two atoms come into contact if within
    % rcut1 but only disconnect when they're further than rcut2
    rcut1 = 3.5;
    rcut2 = 5;

    [contactsPerRun, chainBAtomIndices, contactRes] = mainSim.computeContacts(mainChain, ligandChain, 'StartFrame', settings.frames2skip + 1,'RCut1', rcut1, 'RCut2', rcut2);
    receptorLigandResIds = plot_save_contact_map(contactsPerRun, mainChain, ligandChain, chainBAtomIndices, ...
        'ImportantCutoff', 0.55, 'RCut', rcut2, 'SaveDir', md2pathdir); %FIX MUMBER SHOWN ON TITLE OF HEATMAP

    save(saveVarName,'contactsPerRun','chainBAtomIndices','receptorLigandResIds','-append');
elseif ~isempty(settings.bsResidues)
    receptorLigandResIds = settings.bsResidues;
end


%% Calculate contact map of receptor with G protein

if isempty(settings.gpiResidues) && (mainEntry.hasChain(Chains.gprotein) || (length(database.entries) > 1))
    add2log(md2pathdir, "Calculating Gp interacting residues from downloaded PDB: " + settings.pdbCode);

    if length(database.entries) > 1
        refEntry = database.entries{2};
        [~, receptorGpResIds] = calc_plot_gp_contact_map(mainEntry, refEntry, Chains, 'SaveDir', md2pathdir, 'StartFrame', settings.frames2skip + 1);
    else
        [receptorGpResIds, ~] = calc_plot_gp_contact_map(mainEntry, [], Chains, 'SaveDir', md2pathdir, 'StartFrame', settings.frames2skip + 1);
    end

    % Write out gp contacting residues to file
    writematrix(receptorGpResIds', fullfile(md2pathdir, "GPI_residues.txt"), 'Delimiter', 'space');

    % Save figure as pdf and as .fig
    % figPath = fullfile(md2path, "contact_map" + name + "_G-protein");
    % savefig(figPath);
    % print2pdf(figPath, 1);

    add2log(md2pathdir, ["G-protein contacting residues written to GPI_residues.txt", "Contact maps saved to pdf and fig", ""]);
elseif ~isempty(settings.gpiResidues)
    % Set the variable used by the script for Gp interacting residues
    receptorGpResIds = settings.gpiResidues;
end

%% PCA section:

if settings.doPCA
    add2log(md2pathdir, "PCA analysis commenced");
    pcadir = fullfile(md2pathdir, "pca");

    if ~exist(pcadir, 'dir')
        mkdir(pcadir);
        add2log(md2pathdir, "Creating pca directory in " + md2pathdir);
    end

    %% Do a PCA analysis on the ligand pose
    if  mainEntry.hasChain(Chains.ligand)
        % Ligand: we do PCA on heavy atoms in case of small ligand and CA in
        % case of peptide
         rng(1); % For reproducibility

        atomIndices = ligandChain.getLigandAtoms();
        traj = mainSim.concatRuns('Atoms', atomIndices, 'StartFrame', settings.frames2skip + 1);

        kPrinComp = 2; % Number of principal components to take into consideration
        [indexOfCluster_pca, centroid_pca, p, ind_centers] = cluster_traj_pca(traj, kPrinComp, settings.kClusters, [], settings.kmax);
        save(saveVarName,'indexOfCluster_pca','centroid_pca','p','ind_centers','-append');
        clear traj;

        figPath = fullfile(pcadir, "pca_scatter_" + name);
        savefig( figPath + ".fig");
        print2pdf(figPath);

        add2log(md2pathdir, "PCA clustering of ligand pose performed!");

       %% Run the ligand PCA module
        pca_1_ligand;
        % output:
        % FILE: pcaCenterFrame_run.txt
        %  ligandReceptorpca.pdb
        %  ligandpca.pdb
        %
        % * pca_frame_centers: frames and runs of PCA cluster centers
        %
        save(saveVarName,'rmsdLigandpcaCenter','-append');
    end

    % Run the receptor PCA module
%     pca_2_prot;
    % output: same as pca_ligand but for receptor

end


%% Calculate GPCR order parameters to assess state of the system

if settings.isGPCR
   params = calcPlotOrderParameters(database, database.entries{3}, 'SavePath', md2pathdir,'frames2skip',settings.frames2skip + 1);
   if settings.doPCA
   calcPlotOrderParametersClusters(database, database.entries{3},indexOfCluster_pca , 'SavePath', md2pathdir,'frames2skip',settings.frames2skip + 1);
   end

   save(saveVarName,'params','-append');

end


%% Calculate dihedrals from trajs

if settings.includeLigPathwayCalc
    receptorResIds = mainChain.concatResIds(ligandChain);
else
    receptorResIds = mainChain.resIds;
end

mainSim.computeDihedrals( ...
    'Path', fullfile(md2pathdir, "dihedrals.mat"), ...
    'ReSortPath', fullfile(md2pathdir, "reSort.txt"), ...
    'ResIds', receptorResIds, ...
    'StartFrame', settings.frames2skip + 1 ...
);


%% Calculate MI!

mainSim.computeMI('Path', fullfile(md2pathdir, "MI.mat"));


%% Assess convergence of entropies

assessEntropyConvergence(mainSim, 'SavePath', fullfile(md2pathdir, "convergence"));
add2log(md2pathdir, "S2 convergence calculated and plotted");


%% Prepare input files for path calculation!

pathCalcdir = prepareAlloPathCalc( ...
    mainSim, ...
    mainChain, ...
    settings, ...
    'Dir', md2pathdir, ...
    'LogPath', md2pathdir, ...
    'ReceptorLigandResIds', receptorLigandResIds, ...
    'ReceptorGpResIds', receptorGpResIds, ...
    'ReceptorResIds', receptorResIds ...
);


%% Finally run path calculation!!

% These options are now specified in the input script
% % Distance cutoff for considering "allosteric" signal, default is 10 A
% disCutOff = 10;
%
% % cutoff for pathways to be considered near, and amount of overlap required
% % to cluster them into one pipeline, default to 7.5 and 0.75
% nearcutoff = 7.5;
% overlapcutoff = 0.75;

load(fullfile(pathCalcdir,"workspace.mat"))

prepMI;
MIStatsResLevel;
graphanalysis;
ClusterMIpathways;
analyzeClusters;
writeChannels;
analyzePathDomains;
MIanalysisBS2Effector;
if settings.isGPCR
    pathwayBiasAnalysis;
end
[pHere] = visualizeClsGraph(PDB,pathstruc,Gmatmajor,1,'MIFractionCutoff',MIFractionCutoff);
% If protein is inversed along z, use this command to make it upright again
% set(gca,'zdir','reverse')

% If the "chirality" is off, you can use this:
% set(gca,'ydir','reverse')
