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

%% Load dihedrals and reSorts

dihedrals_mat = cell(length(foldersToStudy),1);
reSortCell = cell(length(foldersToStudy),1);
for thisSys = 1:length(foldersToStudy)
    
    mydir = ([metadir slash foldersToStudy{thisSys}]); 
    
    % Load dihedrals of this system:
    load([mydir slash 'md2path' slash 'dihedrals.mat']); % loads dihedrals 
    % and reSort
    reSortCell{thisSys} = reSort;
    
    % make a matrix of all dihedrals
    dihedrals_mat{thisSys} = [];
    for dih = 1:length(dihedrals)
        X = dihedrals{dih};
        if dih == 1 && isnan(X(1,1))
            X(:,1) = []; % Remove the NaNs from Nterm column 1
        elseif dih == length(dihedrals) && isnan(X(1,2))
            X(:,2) = []; % Remove the NaNs from Cterm column 2
        end
        dihedrals_mat{thisSys} = [dihedrals_mat{thisSys} X];
    end
    dihedrals_mat{thisSys} = real(dihedrals_mat{thisSys});
    [Nframes, Ndih] = size(dihedrals_mat{thisSys});
end

%% Fetch active/inactive state rotamers from available experimental structures

pdbCodes = {'3D4S','3SN6'};
pdbStates = [0 2]; % 0 for inactive, 1 for intermediate, 2 for active
pdbsRef = cell(length(pdbCodes),1);
crdsRef = cell(length(pdbCodes),1);
reSortCellRef = cell(length(pdbCodes),1);
dihedralsRef_mat = cell(length(pdbCodes),1);

% resQueryRef = [122 382];
 resQueryRef = [382];
higherOrder = 'all';

for i = 1:length(pdbCodes) % I'm ignoring 1 since 6CM4 has missing atoms on I122
    refFileName = [databasePath pdbCodes{i} '.pdb'];
    [pdbsRef{i}, crdsRef{i}] = readpdb(refFileName);
    
    % if active state, take receptor only:
    if pdbStates(i) == 2
        indexChainR = selectname(pdbsRef{i}.chainid, 'R');
        pdbTemp = substruct(pdbsRef{i}, indexChainR);
        crdTemp = crdsRef{i}(to3(indexChainR));
    else
        indexChainR = selectname(pdbsRef{i}.chainid, 'A');
        pdbTemp = substruct(pdbsRef{i}, indexChainR);
        crdTemp = crdsRef{i}(to3(indexChainR));
%         pdbTemp = pdbsRef{i};
%         crdTemp = crdsRef{i};
    end
    [dihedralsRef,~,reSortRef] = calcalldihedralsfromtrajs(pdbTemp,crdTemp,resQueryRef,1,higherOrder);
    reSortCellRef{i} = reSortRef;
    
    
    dihedralsRef_mat{i} = [];
    for dih = 1:length(dihedralsRef)
        X = dihedralsRef{dih};
        if dih == 1 && isnan(X(1,1))
            X(:,1) = []; % Remove the NaNs from Nterm column 1
        elseif dih == length(dihedralsRef) && isnan(X(1,2))
            X(:,2) = []; % Remove the NaNs from Cterm column 2
        end
        dihedralsRef_mat{i} = [dihedralsRef_mat{i} X];
    end
    dihedralsRef_mat{i} = real(dihedralsRef_mat{i});
end



%% Plot specific dihedrals of queried residues


% resQuery = [91 217]; % Residues to plot dihedrals for
resQuery = [217]; % Residues to plot dihedrals for

mySys = [ 1 6]; % System to compare (avoids cluttering)
dihedralToPlot = 3; % 1 2 BB dihedrals, 3 -> 5 SC dihedrals

titleName = name; % Clean name so it appears nice in the plot title
titleName(titleName=='_')= ' ';

for resi = 1:length(resQuery)
    figure
    thisRes = resQuery(resi);
   for thisSys = mySys
        thisReSort = reSortCell{thisSys}(reSortCell{thisSys}(:,1)==thisRes,:);
        temp = find(reSortCell{thisSys}(:,1)==thisRes);
        
        % plot chi1 (third element in temp)
        scatter(1:size(dihedrals_mat{thisSys},1),dihedrals_mat{thisSys}(:,temp(dihedralToPlot)),5,'filled');
        hold on
   end
   
   % Add the reference plots:
   refIsUsed = zeros(length(pdbCodes),1); % Keeps track whethere a reference structure was used (for legend)
   for i = 1:length(pdbCodes)
       if  ~isempty(dihedralsRef_mat{i})
           refIsUsed(i) = 1;
           tempRef = find(reSortCellRef{i}(:,1) == resi);
           xValues = dihedralsRef_mat{i}(tempRef(dihedralToPlot))*ones(size(dihedrals_mat{thisSys},1),1);
           plot(xValues,'--','LineWidth',2.0)
           hold on
       end
   end
   
   title([titleName ', Residue ' num2str(thisRes)])
   legend([thisSysLabel(mySys)  pdbCodes(find(refIsUsed))])
   legend boxoff
   
   savefig([md2pathdir 'chi1_res' num2str(thisRes) '_' name]);
   print2pdf([md2pathdir 'chi1_res' num2str(thisRes) '_' name]);
end


%% Try doing RMSD to inactive/active structures instead, as the difference
% is NOT appearing with dihedrals sadly
% Something similar was done in: Dror et al. doi/full/10.1073/pnas.1110499108

