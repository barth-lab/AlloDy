% metaKLDiv: meta-analysis of KLDiv from several directories

% Input section
settings.databasePath = '/PATH/TO/DATABASE'; 
% Father directory for all the runs
settings.metadir = '/PATH/TO/DIR';

settings.foldersToStudy = {'SUBDIR1','SUBDIR2'};

settings.thisSysLabel = {'NAME1','NAME2'};

settings.name = 'NAME';

settings.gpcrdbRefName = 'NAME';

% Reference system for the KLDiv calculations
settings.refdir =   '/PATH/TO/REF/DIR';
settings.refName = 'REF';

% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
settings.chains = 'A--';

% Reference PDB code that this simulation is based on
settings.pdbCode = '3SN6';
settings.pdbChains = 'RA-';

% Inactive state reference PDB, used for dihedral comparison
settings.pdbCodeInactive = '3D4S';
settings.pdbInactiveChains = 'A--';

settings.isGPCR = true;
settings.includeLigPathwayCalc = false;

% Which types of dihedrals to be considered for KL1 calculation:
% phi = 1; psi = 2; chi1 = 3; chi2 = 4; ... chi5 = 7
% settings.dihedralTypes = [1 2 3 4 5 6 7]; 
settings.dihedralTypes = 3;
%% Make the meta-direcory
if isunix
    slash = '/';
    copy = 'cp';
elseif ispc
    slash = '\';
    copy = 'copy';
end
md2pathdirMeta = [settings.metadir slash 'md2pathMeta' slash];
% md2pathdir = [metadir slash 'md2pathMeta' slash name slash];
md2pathdir = [settings.metadir slash 'md2pathMeta' slash];

% Create meta directory
if ~exist(md2pathdirMeta, 'dir')
    mkdir(md2pathdirMeta)
    add2log(md2pathdirMeta, ['Creating md2path directory in ' settings.metadir]);
end

% Create subdirectories for every bundle of systems studied
if ~exist(md2pathdir, 'dir')
    mkdir(md2pathdir)
    add2log(md2pathdir, ['Creating md2path directory in ' md2pathdirMeta]);
end

%% Load, fetch and align PDB files and trajectories
database = Database(settings.databasePath);

% Read PDB of every system to be studied:
for thisSys = 1:length(settings.foldersToStudy)
    % Get local directory
    mydir = ([settings.metadir slash settings.foldersToStudy{thisSys}]); 
    database.read(fullfile(mydir, "prot.pdb"), settings.chains, settings.thisSysLabel{thisSys});

end

% Load the KLDiv reference structure:
database.read(fullfile(settings.refdir, "prot.pdb"), settings.chains, settings.refName );
refEntry = database.entries{(length(settings.foldersToStudy)+1)}; % This entry spot usually reserved for active ref pdb

% Store names for each chain
Chains.receptor = 1;
Chains.gprotein = 2;
Chains.ligand = 3;

refChain = refEntry.chains{Chains.receptor};

if refEntry.hasChain(Chains.ligand)
    ligandChain = refEntry.chains{Chains.ligand};
end

% Get list of residues
if settings.includeLigPathwayCalc
    receptorResIds = refChain.concatResIds(ligandChain);
else
    receptorResIds = refChain.resIds;
end

% Read reference PDB structures
if ~isempty(settings.pdbCode)
    database.fetch(settings.pdbCode, settings.pdbChains);

    if ~isempty(settings.pdbCodeInactive)
        database.fetch(settings.pdbCodeInactive, settings.pdbInactiveChains);
    end
end

% Align sequences and produce the database.residues table
database.align();

% Align structures to the first database entry
database.alignStructures();

% Add labels with Ballesteros-Weinstein notation, aligned with the ref.
% database entry from the PDB
% Exported from GPCRdb (https://gpcrdb.org/residue/residuetable)
if length(database.entries) > (length(settings.foldersToStudy)+1) && isfield(settings, 'gpcrdbRefName') && settings.isGPCR
    residueTablePath = fullfile(database.dir, settings.gpcrdbRefName + "_residue_table.xlsx");
    database.label(length(settings.foldersToStudy)+2, residueTablePath); 
end

%% Read KLDiv from different systems:

kl1 = cell(length(settings.foldersToStudy),1);
kl1Res = cell(length(settings.foldersToStudy),1); % Per residue 
reSortCommon = cell(length(settings.foldersToStudy),1);

for thisSys = 1:length(settings.foldersToStudy)
    mydir = ([settings.metadir slash settings.foldersToStudy{thisSys} slash 'md2path']); 
%     mydir = ([metadir slash foldersToStudy{thisSys} slash]); 
    temp =  importdata([mydir slash 'klDiv.mat']);
%     reSortCommon{thisSys} = importdata([mydir slash 'reSortKLDiv.mat']);
     tempStruct = load([mydir slash 'reSortKLDiv.mat'],'reSortCommon');
    reSortCommon{thisSys} = tempStruct.reSortCommon;
    if isstruct(temp)
        kl1{thisSys} = temp.kl1;
    else
        kl1{thisSys} = temp;
    end
    
    % Get custom residue wise KL1: (according to specified dihedral types)
    kl1Res{thisSys} = zeros(length(receptorResIds),1);

    for i = 1:length(receptorResIds)
        resHere = receptorResIds(i);
        ndxHere = find(reSortCommon{thisSys}(:,1) == resHere); % Index for dihedrals of this residue
%         reSortHere = (reSortCommon{thisSys}(ndxHere,2));
        
%         finalNdx = [];
%         % Dih Types:
%         for dihType = settings.dihedralTypes
%             finalNdx = [finalNdx ];
% 
%         end

        kl1Res{thisSys}(i) = sum(kl1{thisSys}(ndxHere(1:min(length(ndxHere) ,settings.dihedralTypes ))));
    end
end

