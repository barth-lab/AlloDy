% Load trajectories of different systems, and then perform a common PCA to all
% of the common trajs

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
%% Do PCA on more than 1 system

% Grab union of all ligand binding residues
thisList = [];
for thisSys = 1:length(foldersToStudy)
    L = importdata([metadir slash foldersToStudy{thisSys} slash 'md2path' slash 'BS_residues.txt']);
    thisList = [thisList L];

end
liRes = unique(thisList); % Make a list of common ligand binding residues from all systems

tempProtTraj = [];
tempLigTraj = [];
RMSDCell =  cell(length(foldersToStudy),1);
pdbModel = cell(length(foldersToStudy),1);
crdModel = cell(length(foldersToStudy),1);
nFramesEff = cell(length(foldersToStudy),1); % I will use this to extract
% the pdb of the center
trajCeption = cell(length(thisSysLabel),1);

CSys = []; % Contains labels as: [ system run frameNdx ]
for thisSys = 1:length(foldersToStudy)

    mydir = ([metadir slash foldersToStudy{thisSys}]);
    numRuns = numRunsSys(thisSys);
    % load the trajetories:
    import_align_traj;

    % Load input model and align it: (if it exists)
    if exist([mydir slash inputmodelName],'file')
        [pdbModel{thisSys}, temp] = readpdb([mydir slash inputmodelName]);
        [~, crdModel{thisSys}] = superimpose(crd, temp, find(index_CA_all));

    end

    index_CB_all = selectname(pdb.name, 'CB');
    index_CB = index_chainR & index_CB_all; % CBs of receptor only

    % Select ligand binding residues CA and CB only
    ligBindSele = selectid(pdb.resseq,liRes) & (index_CA_all | index_CB_all);

%     sele = index_CB | index_CA; % Do PCA on CA and CB
    sele = index_CA;
    if isLigSmall
        seleLig = index_L; % Non-H atoms
    else
        seleLig =  index_chainL & (index_CA_all | index_CB_all) ; % Same for ligand
    end
    % sele = index_chainR & noH;
    % sele = index_CA; % Selection for doing PCA
    % seleRMSD = index_chainR & noH;% Selection for doing RMSD calcs between centers


    nFramesEff{thisSys} = zeros(numRuns,1);
    RMSDCell{thisSys} = zeros(minTraj,numRuns);
%     startingFrame = 0;
    for runi = 1:numRuns
       tempProtTraj = [tempProtTraj ; traj{runi}(frames2skip:end,to3(sele))]; % Make one fat trajectory
       tempLigTraj = [tempLigTraj ;  traj{runi}(frames2skip:end,to3(seleLig))]; % Make one fat trajectory
       % Calc RMSD along the way
       RMSDCell{thisSys}(:,runi) = calcrmsd(traj{runi}(1,:), traj{runi}(1:minTraj,:),ligBindSele);
%        if runi ~= 1 % Keep startingFrame = 0 for first run,
%          startingFrame = startingFrame + size(traj{runi-1},1);
%        end
       framesNdx = (frames2skip:size(traj{runi},1))';

       nFramesEff{thisSys}(runi) = size(traj{runi}(frames2skip:end,to3(sele)),1);
       Ctemp = [thisSys*ones(nFramesEff{thisSys}(runi),1)  ...
           runi*ones(nFramesEff{thisSys}(runi),1) framesNdx];
       CSys=[CSys ;Ctemp]; % Used for coloring and for labeling: [ system run frameNdx]
    end

    % Add input models if available:
    if exist([mydir '\step1_pdbreader.pdb'],'file')
        tempProtTraj = [tempProtTraj ; crdModel{thisSys}(to3(sele))];
        tempLigTraj = [tempLigTraj ;  crdModel{thisSys}(to3(seleLig))];
        Ctemp = [thisSys 0 0]; % 0 for input models
        CSys=[CSys ;Ctemp]; % Used for coloring and for labeling: [ system run frameNdx]
    end
end

[modelNdx,~] = find(CSys==0);

%% Side quest: get RMSD of ligand binding residues:

% Plot RMSD
figure
errStrideRMSD = ceil(nFrames/60);
for thisSys = 1:length(foldersToStudy)

    RMSD_mean= mean(RMSDCell{thisSys},2);
    RMSD_std= std(RMSDCell{thisSys},0,2);


    plot(RMSD_mean,'LineWidth',1);
    hold on
    errorbar(1:errStrideRMSD:length(RMSD_mean),RMSD_mean(1:errStrideRMSD:length(RMSD_mean)), ...
        RMSD_std(1:errStrideRMSD:length(RMSD_mean))/sqrt(numRuns),'o','MarkerSize',5)
    hold on
    legend_entries{2*thisSys-1} = thisSysLabel{thisSys};
    legend_entries{2*thisSys} = '';
end

xlabel('Frame'); ylabel('RMSD [Angstrom]')
legend(legend_entries,'location','best');
legend boxoff
title('RMSD of ligand binding res, C\alpha and C\beta')
set(gca,'FontSize',16)

savefig([md2pathdir 'RMSD_ligand_binding']);
print2pdf([md2pathdir 'RMSD_ligand_binding']);
%% Do the actual PCA now:

% Protein

kClusters = [];
kmax = 15;
kPrinComp = 2; % Number of principal components to take into consideration

%{[indexOfCluster_pcaProt, centroid_pcaProt, pProt, ind_centersProt] = cluster_traj_pca( tempProtTraj,kPrinComp, kClusters,[],kmax);
%
% savefig([md2pathdir 'pcaProt_cluster_' name]);
% print2pdf([md2pathdir 'pcaProt_cluster_' name]);
% title(['Clustering with ' num2str(kPrinComp) ' receptor PCs']);
%}
%     clear tempProtTraj

% Ligand:

[indexOfCluster_pca, centroid_pca, p, ind_centers] = cluster_traj_pca( tempLigTraj,kPrinComp, kClusters,[],kmax, colors);

title(['Clustering with ' num2str(kPrinComp) ' ligand PCs']);
savefig([md2pathdir 'pcaLig_cluster_' name]);
print2pdf([md2pathdir 'pcaLig_cluster_' name]);
%% Plot with runs:
pcaRange = range(pProt,'all');
pCentersProtLig = zeros(length(foldersToStudy),2);
pInputModelsProtLig = zeros(length(foldersToStudy),2);

figure
% subplot(3,1,1)
% scatter(pProt(:,1),pProt(:,2),5,CSys(:,1),'filled')
% hold on
% subplot(3,1,2)
% scatter(p(:,1),p(:,2),5,CSys(:,1),'filled')
% hold on
% subplot(3,1,3)
scatter(pProt(:,1),p(:,1),5,CSys(:,1),'filled')
hold on
  xlabel('Receptor PC 1', 'fontsize', 25);
  ylabel('Ligand PC 1', 'fontsize', 25);
%
for thisSys = 1:length(foldersToStudy)
    numRuns = numRunsSys(thisSys);

%     (find(CSys(:,1)==thisSys))

    p1 = mean(pProt(find(CSys(:,1)==thisSys),1));
    p2 = mean(p(find(CSys(:,1)==thisSys),1));
    pCentersProtLig(thisSys,:) = [p1 p2];
%     scatter(p1,p2,60,'MarkerEdgeColor',[0 .5 .5],...
%                   'MarkerFaceColor',[0 .7 .7],...
%                   'LineWidth',1.5)
    colorHere = colors(max(round((size(parula,1)/(length(thisSysLabel)-1)*(thisSys-1))),1),:);

    scatter(p1,p2,60,'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',colorHere,...
                  'LineWidth',1.5)
    hold on
    text(p1+pcaRange/50,p2+pcaRange/50,thisSysLabel{thisSys},'FontSize',14)

    % If input model is available, plot it!
    if ~isempty(crdModel{thisSys})
        p1 = pProt(modelNdx(thisSys),1);
        p2 = p(modelNdx(thisSys),1);
        scatter(p1,p2,60,'MarkerEdgeColor',[0.5 .5 0],...
                      'MarkerFaceColor',[0.7 .7 0],...
                      'LineWidth',1.5)
        hold on
        text(p1+pcaRange/50,p2+pcaRange/50,[thisSysLabel{thisSys} 'M'],'FontSize',14)
        pInputModelsProtLig(thisSys,:) = [p1 p2];
    end


%     for runi = 1:numRuns
%         % Find centroids of the scatter of each run and label it
%         if runi==1
%            lower_lim = 1;
%         else
%            lower_lim = sum(nFramesEff(1:runi-1))+1;
%         end
%         upper_lim = sum(nFramesEff(1:runi));
%         p1 = mean(pProt(lower_lim:upper_lim,1));
%         p2 = mean(pProt(lower_lim:upper_lim,2));
%
%         scatter(p1,p2,60,'MarkerEdgeColor',[0 .5 .5],...
%                       'MarkerFaceColor',[0 .7 .7],...
%                       'LineWidth',1.5)
%                   hold on
%         % Add text label for runs with proper offset
%         text(p1+pcaRange/50,p2+pcaRange/50,['Run ' num2str(runi)],'FontSize',14)
%     end
end

%   xlabel('PC 1', 'fontsize', 25);
%   ylabel('PC 2', 'fontsize', 25);
%   title(['Scatter plot of ' num2str(kPrinComp) ' receptor PCs, colored by run'])
    legend('Data colored by run','Centroids of runs','Input models')
  legend boxoff

  savefig([md2pathdir 'pcaProtLig_system_scatter_' name]);
print2pdf([md2pathdir 'pcaProtLig_system_scatter_' name]);
writematrix(pCentersProtLig,[md2pathdir 'pCentersProtLigCoordPC1.txt'],'Delimiter','space')

  %% Ligand only PCAs:
  pcaRange = range(p,'all');
  pCentersLig = zeros(length(foldersToStudy),2);
  pInputModelsLig = zeros(length(foldersToStudy),2);
  ind_centersRunsLig = zeros(length(foldersToStudy),nCentersLig);
  run_FrameNdx = cell(length(foldersToStudy),1);

  figure

colors = colormap('parula');

% Plot the PCs first
for thisSys = 1:length(foldersToStudy)
    colorHere = colors(max(round((size(parula,1)/(length(thisSysLabel)-1)*(thisSys-1))),1),:);
    scatter(p(CSys(:,1)==thisSys,1),p(CSys(:,1)==thisSys,2),5,colorHere,'filled');
    hold on
end

% Now plot centers and add text labels
for thisSys = 1:length(foldersToStudy)
    numRuns = numRunsSys(thisSys);

%     (find(CSys(:,1)==thisSys))

    p1 = mean(p(find(CSys(:,1)==thisSys),1));
    p2 = mean(p(find(CSys(:,1)==thisSys),2));
    pCentersLig(thisSys,:) = [p1 p2];

    colorHere = colors(max(round((size(parula,1)/(length(thisSysLabel)-1)*(thisSys-1))),1),:);

    scatter(p1,p2,60,'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',colorHere,...
                  'LineWidth',1.5)
    hold on
    text(p1+pcaRange/50,p2+pcaRange/50,thisSysLabel{thisSys},'FontSize',14)

    % If input model is available, plot it!
    if ~isempty(crdModel{thisSys})
        p1 = p(modelNdx(thisSys),1);
        p2 = p(modelNdx(thisSys),2);
        scatter(p1,p2,60,'MarkerEdgeColor',[0.5 .5 0],...
                      'MarkerFaceColor',[0.7 .7 0],...
                      'LineWidth',1.5)
        hold on
        text(p1+pcaRange/50,p2+pcaRange/50,[thisSysLabel{thisSys} 'M'],'FontSize',14)
        pInputModelsLig(thisSys,:) = [p1 p2];
    end

    % Save the pdbs of the cluster centers:
    % Only consider pcs from this system:
    pTemp = p(CSys(:,1)==thisSys,1:kPrinComp);
    % The closest point to the mean point
   [~,ind_centersRunsLig(thisSys)] = min(vecnorm(pTemp-pCentersLig(thisSys,:),2,2));

   % The X closest points to the mean point
   [B,I] = sort(vecnorm(pTemp-pCentersLig(thisSys,:),2,2));
   ind_centersRunsLig(thisSys,:) = I(1:nCentersLig);

   % ind_centersRunsLig is the index of the center of all the runs of a
   % system, the number is the frame number of given sys with @frames2skip
   % removed

    % Now use frame index to grab the real frame and extract it from traj:
    CSysTemp = CSys(CSys(:,1)==thisSys,:); % Csys for this system only

    run_FrameNdx{thisSys} = CSysTemp(ind_centersRunsLig(thisSys,:),2:3); % Got 'em!

end
  xlabel('PC 1', 'fontsize', 25);
  ylabel('PC 2', 'fontsize', 25);
%     legend('Data colored by run','Centroids of runs','Input models')
  legend(thisSysLabel)
  legend boxoff
  set(gca,'FontSize',16)

   title(['Scatter plot of ' num2str(kPrinComp) ' ligand PCs, colored by run'])

   % Save the plot
savefig([md2pathdir 'pcaLig_system_scatter_' name]);
print2pdf([md2pathdir 'pcaLig_system_scatter_' name]);
% Write the PCA values of the centers and the models
writematrix(pCentersLig,[md2pathdir 'pCentersLigCoordPC1_2.txt'],'Delimiter','space')
writematrix(p(CSys(:,3)==0,1:2),[md2pathdir 'pModelCentersLigCoordPC1_2.txt'],'Delimiter','space')

%% Caclulate density plots of every system and save highest density cluster
% centers

  pHighestDen = zeros(length(foldersToStudy),2);
  ind_HighestDenLig = zeros(length(foldersToStudy),nCentersLig);
  run_FrameHighestDenNdx = cell(length(foldersToStudy),1);
  figure
  nPlots = ceil(sqrt(length(foldersToStudy)));

for thisSys = 1:length(foldersToStudy)
    subplot(nPlots,nPlots,thisSys);
     % Only consider pcs from this system:
    pTemp = p(CSys(:,1)==thisSys,1:kPrinComp);

    [values, centers] = hist3([pTemp(:, 1) pTemp(:, 2)],[51 51]); % bins may need modification
    imagesc(centers{:},values.')
    colorbar
    axis xy
    xlabel('PCA 1', 'fontsize', 25);
    ylabel('PCA 2', 'fontsize', 25);
    title(['Density plot of the first ' num2str(kPrinComp) ' PCs, ' thisSysLabel{thisSys}])

    % Find highest density point/s

    % How to decide how many maxes to take?
    maxDen = max(values,[],'all');
    [a,b] = find(values==maxDen);
    pHighestDen(thisSys,:) = [centers{1}(a) centers{2}(b)]; % Highest density PC

    % Find closest PC point to our highest density center
   [~,ind_HighestDenLig(thisSys)] = min(vecnorm(pTemp-pHighestDen(thisSys,:),2,2));
   % The X closest points to the mean point
   [B,I] = sort(vecnorm(pTemp-pHighestDen(thisSys,:),2,2));
   ind_HighestDenLig(thisSys,:) = I(1:nCentersLig);

   CSysTemp = CSys(CSys(:,1)==thisSys,:); % Csys for this system only

   run_FrameHighestDenNdx{thisSys} = CSysTemp(ind_HighestDenLig(thisSys,:),2:3); % Got 'em!
    % extract structure et voila!

end
savefig([md2pathdir 'pcaLig_density_' name]);
print2pdf([md2pathdir 'pcaLig_density_' name]);
% Write the PCA values of the highest density points
writematrix(pHighestDen,[md2pathdir 'pHighestDenLigCoordPC1_2.txt'],'Delimiter','space')

%% Finally save the pdbs of the centers of the clusters of separate systems:

% CSys and PCs of the ligand cluster centers
CSysClusterCenters = [CSys(ind_centers,:)  (1:length(ind_centers))'];
pClusterCenters = p(ind_centers,1:kPrinComp);
pdbCell = cell(length(foldersToStudy),1);
indexLigCell = cell(length(foldersToStudy),1);

 pcaTraj = cell(length(foldersToStudy),1);
 pcaTrajHighestDen = cell(length(foldersToStudy),1);
 pcaClusterTraj = cell(length(ind_centers),1);
  for thisSys = 1:length(foldersToStudy)

    mydir = ([metadir slash foldersToStudy{thisSys}]);
    numRuns = numRunsSys(thisSys);
    % load the trajetories:
    import_align_traj;
    pdbCell{thisSys} = pdb;
    indexLigCell{thisSys} = index_chainL;
    % Save pdbs of centers of systems
    for modelN = 1:nCentersLig
        pcaTraj{thisSys} = [ pcaTraj{thisSys} ;traj{run_FrameNdx{thisSys}(modelN,1)}(run_FrameNdx{thisSys}(modelN,2),:)];
        pcaTrajHighestDen{thisSys} = [ pcaTrajHighestDen{thisSys} ;traj{run_FrameHighestDenNdx{thisSys}(modelN,1)}(run_FrameHighestDenNdx{thisSys}(modelN,2),:)];
    end
    writepdbndx([md2pathdir 'pcaLigand_' thisSysLabel{thisSys} '.pdb'], pdb, [], 0, pcaTraj{thisSys});

    % Save pdbs of highest density points of every system
    writepdbndx([md2pathdir 'pcaLigandHighestDen_' thisSysLabel{thisSys} '.pdb'], pdb, [], 0, pcaTrajHighestDen{thisSys});

    % Save pdbs of PCA cluster centers:

    for thisCenter = 1:length(ind_centers)
        if CSysClusterCenters(thisCenter,1) == thisSys % Cluster is in this system, save PDB!
            pcaClusterTraj{thisCenter} = traj{CSysClusterCenters(thisCenter,2)}(CSysClusterCenters(thisCenter,3),:);
            writepdbndx([md2pathdir 'pcaClusterLigand_' thisSysLabel{thisSys} '_C' num2str(thisCenter) '.pdb'], pdb, [], 0, pcaClusterTraj{thisCenter});
        end
    end

  end

