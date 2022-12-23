%% Make the meta-direcory
if isunix
    slash = '/';
    copy = 'cp';
elseif ispc
    slash = '\';
    copy = 'copy';
end
md2pathdirMeta = [metadir slash 'md2pathMeta' slash];
md2pathdir = [metadir slash 'md2pathMeta' slash name slash];

% Create meta directory
if ~exist(md2pathdirMeta, 'dir')
    mkdir(md2pathdirMeta)
    add2log(md2pathdirMeta, ['Creating md2path directory in ' metadir]);
end

% Create subdirectories for every bundle of systems studied
if ~exist(md2pathdir, 'dir')
    mkdir(md2pathdir)
    add2log(md2pathdir, ['Creating md2path directory in ' md2pathdirMeta]);
end

%% Load, fetch and align PDB files

metaLoadAlign

%% Load contact maps from separate folders

% I take the maximum contact from every contact map, which is not STRICTLY
% CORRECT, since it may miss some overlapping contacts between different
% atoms from the ligand to the same protein residue 

contactCutMeta = 0.3; % Cutoff for signficant contacts, ALL rows should 
% show contact less than cutoff or else it won't zero them out

meta_contacts = [];
for thisSys = 1:length(foldersToStudy)
    L = load([metadir slash foldersToStudy{thisSys} slash 'md2path' slash  'workspace.mat']);
    meta_contacts = [meta_contacts ; max(mean(L.contactsPerRun,3))]; 
end


% meta_contacts(meta_contacts<contactCutMeta)=0; % Zero non-significant contacts
toCut = sum(meta_contacts<contactCutMeta); % If ALL variants have less than 
% cutoff contact frequency, then I'll zero the column

for i = 1:length(toCut)
   if toCut(i) == length(foldersToStudy)
       meta_contacts(:,i) = zeros(length(foldersToStudy),1);
   end
end

[row, col] = find(meta_contacts); % Find nonzero contacts over all runs (averaged)
contact_map = meta_contacts(unique(row),unique(col));

liProtRes = receptorResIds(unique(col)); % Contacting receptor residues

%% Plot the figure:

% Grab BW numbering of ligand binding residues from database/chain
resname_num = refChain.formatResidues(liProtRes,'SecondaryEntry',database.entries{length(foldersToStudy)+1},'BWonly', true);
X = round(contact_map,2,"significant");
figure
h = heatmap(resname_num,thisSysLabel,X,'Colormap',parula)

% xlabel('Protein residue', 'FontSize', fontsz)
h.Title = 'Ligand contact frequency from simulations';
h.XLabel = 'Protein residue';
h.YLabel = 'Ligand';
set(gca,'FontSize',12)

%% Save!!

savefig([md2pathdir 'contact_map_' name]);
print2pdf([md2pathdir 'contact_map_' name]);


%% Why not do RMSD too:

legend_count = 1;
errStrideRMSD = 50;
figure
for thisSys = 1:length(foldersToStudy)
    L = load([metadir slash foldersToStudy{thisSys} slash 'md2path' slash  'workspace.mat']);
    
    RMSD_mean_L= mean(L.RMSD_L,2);
    RMSD_std_L= std(L.RMSD_L,0,2);
    plot(RMSD_mean_L,'LineWidth',1);
    hold on
    errorbar(1:errStrideRMSD:length(RMSD_mean_L),RMSD_mean_L(1:errStrideRMSD:length(RMSD_mean_L)),...
        RMSD_std_L(1:errStrideRMSD:length(RMSD_mean_L))/sqrt( size(L.RMSD,2)),'o','MarkerSize',5)
    legend_entries{legend_count} = thisSysLabel{thisSys};
    legend_count = legend_count + 1;
    legend_entries{legend_count} = '';
    legend_count = legend_count + 1;
end

% Add annotations
xlabel('Frame'); ylabel('RMSD [Angstrom]')
legend(legend_entries,'location','best');
legend boxoff
title('RMSD of non-H ligand atoms')
set(gca,'FontSize',16)


%% Save!!

savefig([md2pathdir 'ligand_RMSD_' name]);
print2pdf([md2pathdir 'ligand_RMSD_' name]);