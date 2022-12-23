clear conectrecord

strongMIlimit = 0.0; % This is obsolete with the proper filtering step we have now
MIavgnear = mean(MIres(dismat<disCutOff & MIres>0));
MIavgfar = mean(MIres(dismat>disCutOff & MIres>0));
Gmat = zeros(Nres,Nres);
Gmatmajor = zeros(Nres,Nres);
MImax = max(reshape(MIres,Nres*Nres,1));
MIavg = mean(MIres(MIres>0));

% connect residues closer than disCutOff in the graph using their MI values
% (inverted so that min distance becomes max MI)
for i=1:Nres-1
    for j=i+1:Nres
        if (dismat(i,j)<=disCutOff)
            Gmat(i,j) = 1.1*MImax-MIres(i,j);
            Gmat(j,i) = Gmat(i,j);
            Gmatmajor(i,j) = Gmat(i,j);
            Gmatmajor(j,i) = Gmat(j,i);
        end
    end
end
% graphallshortestpaths is to be removed from future Matlab releases
% G = sparse((Gmatmajor));
% Gdis = graphallshortestpaths(G, 'Directed', 'false');
G = graph(Gmatmajor);
Gdis = distances(G);

% filter residue pairs that are farther than disCutOff (defaults to 10A)
% Also filter any residue pairs that do not have allosteric paths between
% them

%MIpathsum = zeros(Nres,Nres);
pathstruc = struct();
count = 1;
Npath = 0;
for i=1:Nres-1
    for j=i+1:Nres
        % eliminate res pairs in direct contact or weak MI
        if( (MIres(i,j)<=strongMIlimit*MIavg) || (dismat(i,j)<=disCutOff) )
            Gdis(i,j) = 0;
            Gdis(j,i) = 0;
            MIallo = 0;
        elseif isinf(Gdis(i,j)) % eliminate res pairs with no direct path
            MIallo = 0;
        else
            MIallo = MIres(i,j);
        end

        % Calculate "shortest path" for residues further than disCutOff and
        % having significant MI
        if MIallo 

%             graphshortestpath is going to be removed in future Matlab
%             releases
%             [dum1,path,dum2] = graphshortestpath(G,i,j,'Directed', 'false');
            path = shortestpath(G,i,j);
            if length(path)>=3
                pathstruc(count).path = path;
                pathstruc(count).Npath = length(path);
                sum1 = 0;
                for k=1:length(path)-1
                    sum1 = sum1+MIres(path(k),path(k+1));
                end
                MIpathsum(i,j) = sum1/(length(path)-1);
                MIpathsum(j,i) = sum1/(length(path)-1);
                pathstruc(count).MI = MIallo;
                pathstruc(count).meanpathMI = sum1/(length(path)-1);
                %pathstruc(count).Gdis = Gdis(i,j)/(length(path)-1);
                %[sum1 Gdis(i,j)]
                Npath = Npath+1;
            else
                pathstruc(count).path = [];
                pathstruc(count).Npath = 0;
                pathstruc(count).MI = 0;
                pathstruc(count).meanpathMI = 0;
            end
        else
            pathstruc(count).path = [];
            pathstruc(count).Npath = 0;
            pathstruc(count).MI = 0;
            pathstruc(count).meanpathMI = 0;
            %pathstruc(count).Gdis = 0;
        end
        count = count+1;
    end
end

% sort allosteric pairs by MI


%[Y,I] = sortrows(reshape(MIallo,Nres*Nres,1),-1);
[Y,I] = sortrows([pathstruc(:).MI]',-1);
rankMI(I) = (1:length(I));
%[X,J] = sortrows(reshape(Gdis,Nres*Nres,1),1);
[X,J] = sortrows([pathstruc(:).meanpathMI]',-1);
rankmeanpathMI(J) = (1:length(J));
[dum,K] = sortrows((rankMI.*rankmeanpathMI)',1);


cumMI = cumsum([pathstruc(I).MI]);
Npath = find(cumMI > cumMI(end) * MIFractionCutoff, 1);

% Moved to options in auto script
% nearcutoff = 7.5;
% overlapcutoff = 0.75;


%% Plot pair MI

figure;
hold on;

cumMIfiltered = normalizeMiCumsum(MI);

significantCount = find(cumMIfiltered == 100, 1);
a = plot(normalizeMiCumsum(MI));
xline(significantCount, 'Color', a.Color, 'LineStyle', '--');

plot(normalizeMiCumsum(MIraw));

xlabel("Ranked dihedral pairs");
ylabel("MI normalized cumulative sum");
ytickformat('percentage');

title("Pair MI");
legend("Filtered MI", sprintf("Top %.1f%% dihedrals", significantCount / length(cumMIfiltered) * 100), "Raw MI", 'Location', 'best');
legend boxoff;


%% Plot pathway MI

figure;
hold on;
plot(cumMI / cumMI(end) * 100, 'LineWidth', 1);

xline(Npath, 'Color', 'red', 'LineStyle', '--');
yline(cumMI(Npath) / cumMI(end) * 100, 'Color', 'red', 'LineStyle', '--');

xlabel("Ranked pathways");
ylabel("MI normalized cumulative sum");
ytickformat('percentage');

title("Pathway MI");
legend("Filtered MI", ("Top " + num2str(Npath / length(pathstruc) * 100, 3) + "% pathways, " + num2str(cumMI(Npath) / cumMI(end) * 100, 3) + "% MI"), 'Location', 'best');
legend boxoff;

savefig(fullfile(pathCalcdir,"pathway_MI_percentage"));
print2pdf(fullfile(pathCalcdir,"pathway_MI_percentage"));

%% Utility functions

function output = flattenTriu(input)
    output = nonzeros(input - tril(input) + triu(ones(size(input)), 1)) - 1;
end

function output = normalizeMiCumsum(mi)
    output = cumsum(sort(flattenTriu(mi), 'descend'));
    output = output / output(end) * 100;
end
