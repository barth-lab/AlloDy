function params = calcPlotOrderParameters(database, refEntry, options)
  arguments
    database
    refEntry
    options.SaveName = "%s"
    options.SavePath
    options.frames2skip = 1 % 1 does not skip any frames, ironically
    options.ColorMap % For example, color by PCA cluster index: indexOfCluster_pca
    
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
  plots{1}.name = "param_tm36_npxxy";

  plots{2}.x = params.tm36;
  plots{2}.y = params.tm37;
  plots{2}.title = "TM3-6 vs TM3-7";
  plots{2}.name = "param_tm36_tm37";

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

    figure('Name', currentPlot.title, 'Position', [10 10 1200 600]);
    tiledlayout(1, 2);

    nexttile;
   if exist('options.ColorMap','var')
        scatter(mainX, mainY, 5, options.ColorMap, '.');
   else
       for thisRun = 1:size(tempX,1)
            scatter(tempX{thisRun}', tempY{thisRun}', 5, '.');
            hold on
       end
   end
    hold on;

    formatplot2;

    xlabel(currentPlot.x.label);
    ylabel(currentPlot.y.label);

    xlim(xl);
    ylim(yl);

    nexttile;

    xi = xl(1):0.05:xl(2);
    yi = yl(1):0.02:yl(2);

    gaussianWidth = [0.1 0.1];
    pmf = calcpmf2d([mainX, mainY], xi, yi, gaussianWidth);
    % contourLevels = linspace(min(pmf, [], 'all'), max(pmf, [], 'all'), 10);
    contourLevels = 0:0.25:3.5;

    landscape(xi, yi, pmf, contourLevels);
    hold on;

    handles = zeros(1, length(database.entries) - 1);

    for entryIndex = 2:length(database.entries)
      h = scatter(currentPlot.x.data{entryIndex}, currentPlot.y.data{entryIndex}, 100, 'o', 'DisplayName', database.entries{entryIndex}.name);

      set(h, 'MarkerFaceColor', get(h, 'MarkerEdgeColor'));
      set(h, 'MarkerEdgeColor', 'w');

      handles(entryIndex - 1) = h;
    end

    legend(handles);
    legend boxoff;
    colorbar;
    xlabel(currentPlot.x.label);

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
