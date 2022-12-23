%% Prepare input files: common input files for all "meta" scripts
% Every "meta" script needs a different set of input, as shown below:

databasePath = '/PATH/TO/DATABASE'; 
% Father directory for all the runs
metadir = '/PATH/TO/DIR';
foldersToStudy = {'SUBDIR1','SUBDIR2'};

thisSysLabel = {'NAME1','NAME2'};

name = 'NAME';


%% Parameters for metaPathCalcAnalysis.m
abc = char(65:65+25); % Can you sing the abc?

% Analysis parameters:
nPipes = 10; % number of pipelines to consider in analysis
topHubs = 25; % Number of top scoring global hubs to consider
isGPCR = true;

%% Parameters for metaContactMap.m
% Generic GPCR numbering from https://gpcrdb.org/residue/residuetable
% You need to DL the table manually from the website for every receptor you
% want to study and save it as:
%
% gpcrdbRefName_residue_table.xlsx
%
% Replace gpcrdbRefName with the system name defined below, of course!
gpcrdbRefName = 'NAME';

%% Parameters for metapca_4_multisys and metaDihedralAnalysisProtorype.m

numRunsSys = [5 5]; % Number of runs available for each system studied
frames2skip = 500; % remove first 50 ns of the simulation
nCentersLig = 10; % Number of points closest to the mean point of a system 
% consider for analysis and pdb saving
pdbName = 'prot.pdb';
inputmodelName = 'model.pdb'; % USeful for reference in PCA plots