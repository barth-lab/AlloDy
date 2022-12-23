%% Do PCA for the receptor coordinates

atomIndices = mainChain.getAtoms();
traj = mainSim.concatRuns('Atoms', atomIndices, 'StartFrame', settings.frames2skip + 1);

kPrinComp = 2; % Number of principal components to take into consideration
[indexOfCluster_pcaProt, centroid_pcaProt, pProt, ind_centersProt] = cluster_traj_pca(traj, kPrinComp, settings.kClusters);
clear traj;

figPath = fullfile(pcadir, "prot_pca_scatter" + mainEntry.name);
savefig(figPath);
print2pdf(figPath);

add2log(md2pathdir, "PCA clustering of receptor performed!");


%% Plot PCA1/2 colored by run and by cluster

figure;
pcaRange = range(pProt, 'all');

scatter(pProt(:,1),pProt(:,2),5,C(:,1),'filled');
hold on;

for i = 1:mainSim.runCount
    % Find centroids of the scatter of each run and label it
    if i==1
       lower_lim = 1;
    else
       lower_lim = sum(nFramesEff(1:i-1))+1;
    end
    upper_lim = sum(nFramesEff(1:i));
    p1 = mean(pProt(lower_lim:upper_lim,1));
    p2 = mean(pProt(lower_lim:upper_lim,2));

    scatter(p1,p2,60,'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',[0 .7 .7],...
                  'LineWidth',1.5)
              hold on
    % Add text label for runs with proper offset
    text(p1+pcaRange/50,p2+pcaRange/50,['Run ' num2str(i)],'FontSize',14)
end
  xlabel('PC 1', 'fontsize', 25);
  ylabel('PC 2', 'fontsize', 25);
  title(['Scatter plot of ' num2str(kPrinComp) ' receptor PCs, colored by run'])
    legend('Data colored by run','Centroids of runs')
  legend boxoff

% Save the plot
figPath = fullfile(pcadir, "prot_pca_run_scatter" + mainEntry.name);
savefig(figPath);
print2pdf(figPath);


%% Now get the representative structures:
% ind_centers contains the frame numbers of the cluster centers

pca_frame_centers = zeros(length(ind_centersProt),3);
protpcaTraj = zeros(length(ind_centersProt),mainEntry.atomCount * 3);
for thisCenter = 1:length(ind_centersProt)
    % Grab the run that the center belongs to:
    runCenter = C(ind_centersProt(thisCenter)); %Misuse color vector to grab the run
    % Make a workaround for finding lower limit with trajectories of
    % unequal length
    if runCenter==1
       lower_lim_center = 1;
    else
       lower_lim_center = sum(nFramesEff(1:runCenter-1))+1;
    end
    pca_frame_centers(thisCenter,1) = ind_centersProt(thisCenter) - lower_lim_center + settings.frames2skip + 1;
    pca_frame_centers(thisCenter,2) = runCenter;
    pca_frame_centers(thisCenter,3) = length(find(indexOfCluster_pcaProt==thisCenter));
    protpcaTraj(thisCenter,:) = mainSim.traj{runCenter}(pca_frame_centers(thisCenter,1),:);
end

% Save frame center data and write out pdbs:
writematrix(pca_frame_centers, fullfile(pcadir, "prot_pcaCenterFrame_run_Nelements.txt"), 'Delimiter','space');
% writepdbndx([pcadir 'ligandpca.pdb'], pdb, index_chainL, 0, ligandpcaTraj);
writepdbndx(fullfile(pcadir, sprintf("prot_pcaLigandReceptor_%s.pdb", mainEntry.name)), mainEntry.pdb, [], 0, protpcaTraj);

add2log(md2pathdir,{'frame_centers saved to prot_pcaCenterFrame_run.txt','pdbs saved to prot_pcaLigandReceptor.pdb'});


%% Calculate RMSD matrix between centers

rmsdLigandpcaCenter = zeros(length(ind_centersProt));
seleRMSD = mainEntry.getAtoms('Chain', mainChain.index, 'NoName', 'H*');

for i =1:length(ind_centersProt)-1
    for j = i+1:length(ind_centersProt)
        rmsdLigandpcaCenter(i,j) = calcrmsd(protpcaTraj(i,to3(seleRMSD)),protpcaTraj(j,to3(seleRMSD)));
    end
end

figure;
hm = heatmap(rmsdLigandpcaCenter,'Colormap',parula);
hm.Title = ['Receptor RMSD [' char(197) '] matrix PCA centers '];
hm.XLabel = 'PCA Cluster';
hm.YLabel = 'PCA Cluster';

figPath = fullfile(pcadir, "prot_pcarmsdMatrix");
savefig(figPath);
print2pdf(figPath);
