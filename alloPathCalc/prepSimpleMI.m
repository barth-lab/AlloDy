% prepSimpleMI: a barebones prepMI function:

randMIlim = 0; % This is obsolete with the new filterCorrectMI function (Thanks Simon!)
reSortRenum = load(fullfile(pathCalcdir,"reSort.txt"));
resind = reSortRenum(:,1);
dihtype = reSortRenum(:,2);
dihedral = reSortRenum(:,4:5);
excluded = load(fullfile(pathCalcdir,"excluded_res.txt"));

if isGPCR % Add more general secondary structure?
    % read helices definitions
    helices = load(fullfile(pathCalcdir,"helices.txt"));
end

excluded_res = [];
for i=1:size(excluded,1)
    excluded_res = [excluded_res excluded(i,1):excluded(i,2)];
end

MI = load(miName,'I');
MI = MI.I;
MI(MI<0) = 0;
MI(isnan(MI)) = 0;
MIraw = MI; % Save unfiltered MI for comparison with filtered


%% Filter MI according to a statistical significance criterion:

load(dihedralsPath);
MI = filterCorrectMI(dihedralsMat, MIraw, 'SaveDir', pathCalcdir,'SignificanceThreshold',0.05);

%% get residue list
PDB = pdbread((fullfile(pathCalcdir,"protRenum.pdb")));
atomno = [PDB.Model.Atom(:).AtomSerNo];
resname = cellstr(char(PDB.Model.Atom(:).resName));
resno = [PDB.Model.Atom(:).resSeq];
atomname = char(PDB.Model.Atom(:).AtomName);
CAind = find( strncmp(cellstr(atomname),'CA',3) );
reslist = unique(resno);
atomnameStr = cellstr(atomname);

%%

Nres = length(reslist);
MIres = zeros(Nres,Nres);
% compute residue summed MI

for i=1:Nres-1
    for j=i+1:Nres
        temp = MI(resind==(i),resind==(j));
        temp = temp(:);
        MIres(i,j) = sum(temp);
        MIres(j,i) = MIres(i,j);
    end
end


% Clear excluded residues by zeroing out their MI 
for i=1:Nres
    for j=1:Nres
        if ismember(i,excluded_res) || ismember(j,excluded_res) 
            MIres(i,j) = 0;
        end
    end
end

save(fullfile(pathCalcdir, 'MIres'),'MIres');

figure;
imagesc(MIres)
colorbar
xlabel('Residue')
ylabel('Residue')
title('Residue-wise filtered MI')
