function params = calcPlotOrderParametersClusters(database, refEntry, clusterNdx, options)
  arguments
    database
    refEntry
    clusterNdx % For example, PCA cluster index: indexOfCluster_pca
    options.SaveName = "%s"
    options.SavePath
    options.frames2skip = 1 % 1 does not skip any frames, ironically
    options.nClusters =  min(max(clusterNdx),3) % Defaults to plotting 3 clusters, or less if there are less than 3
  end

  res3 = database.findResidue("3.50");
  res6 = database.findResidue("6.30");
  res7 = database.findResidue("7.53");


  % Parameters

  params.npxxy.data = database.calcRmsd(database.findFeature("7.53", 4, 0), refEntry);
  params.npxxy.label = ["RMSD of NPxxY with respect", "to inactive (" + refEntry.name + ") (Å)"];

  params.tm36.data = database.calcDistance(res3, res6);
  params.tm36.label = "TM3-6 distance (Å)";

  params.tm37.data = database.calcDistance(res3, res7);
  params.tm37.label = "TM3-7 distance (Å)";


  % Plots
  plots = cell(2, 1);
  plots{1}.x = params.tm36;
  plots{1}.y = params.npxxy;
  plots{1}.title = "TM3-6 vs NPxxY";
  plots{1}.name = "param_tm36_npxxy_clusters";

  plots{2}.x = params.tm36;
  plots{2}.y = params.tm37;
  plots{2}.title = "TM3-6 vs TM3-7";
  plots{2}.name = "param_tm36_tm37_clusters";

  for currentPlotIndex = 1:size(plots)
    currentPlot = plots{currentPlotIndex};

    tempX = cell(size(currentPlot.x.data{1},1),1);
    tempY = cell(size(currentPlot.x.data{1},1),1);
    for thisRun=1:size(currentPlot.x.data{1},1)
        tempX{thisRun} = currentPlot.x.data{1}{thisRun}((options.frames2skip):end);
        tempY{thisRun} = currentPlot.y.data{1}{thisRun}((options.frames2skip):end);
    end
    mainX = cell2mat(tempX(:));
    mainY = cell2mat(tempY(:));
%     mainX = [];
%     mainY = [];
%     for thisRun = 1:size(tempX,1) % Make mainX into a linear vector
%         mainX = [ mainX tempX{thisRun}];
%         mainY = [ mainY tempY{thisRun}];
%     end

    xl = minmaxAll([mainX; cell2mat(currentPlot.x.data(2:end))']);
    yl = minmaxAll([mainY; cell2mat(currentPlot.y.data(2:end))']);

    figure('Name', currentPlot.title, 'Position', [10 10 1500 600]);
    tiledlayout(1, 3); % Plot the top 3 clusters

    xi = xl(1):0.05:xl(2);
    yi = yl(1):0.02:yl(2);

    gaussianWidth = [0.1 0.1];
    
    for thisCluster = 1:options.nClusters
    nexttile;
%      scatter(mainX, mainY, 10, clusterNdx, 'o','filled');
    pmf = calcpmf2d([mainX(clusterNdx==thisCluster), mainY(clusterNdx==thisCluster)], xi, yi, gaussianWidth);
    % contourLevels = linspace(min(pmf, [], 'all'), max(pmf, [], 'all'), 10);
    contourLevels = 0:0.25:3.5;
    
  
    

    landscape(xi, yi, pmf, contourLevels);
    hold on;
    
    if thisCluster ==3  
        colorbar;    
    end
    
    handles = zeros(1, length(database.entries) - 1);

    for entryIndex = 2:length(database.entries)
      h = scatter(currentPlot.x.data{entryIndex}, currentPlot.y.data{entryIndex}, 100, 'o', 'DisplayName', database.entries{entryIndex}.name);

      set(h, 'MarkerFaceColor', get(h, 'MarkerEdgeColor'));
      set(h, 'MarkerEdgeColor', 'w');

      handles(entryIndex - 1) = h;
    end

    xlabel(currentPlot.x.label);
    if thisCluster==1
        ylabel(currentPlot.y.label);
    
        legend(handles, 'location','best');
        legend boxoff;
    end
    title(['Cluster C' num2str(thisCluster)])
    xlim(xl);
    ylim(yl);
    formatplot2;

    end


    sgtitle(currentPlot.title, 'FontName','Helvetica', 'FontSize', 25);

    if isfield(options, 'SavePath')
      figPath = fullfile(options.SavePath, sprintf(options.SaveName, currentPlot.name));
      savefig(figPath);
      print2pdf(figPath);
    end
  end
end


%% Functions

function output = minmaxAll(arr)
  output = [min(arr, [], 'all'), max(arr, [], 'all')];
end
