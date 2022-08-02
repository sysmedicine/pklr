%% Figure 2A
content = Cheng_readTXT('tpm_name.txt');
sampleIDs = content(1,1:end-1);
genes = content(2:end,1);
expMat = str2double(content(2:end,2:end));

metaData = Cheng_readTXT('metadata.txt');
metaData(cellfun(@isempty,metaData(:,1)),:) = [];

ChowCTRLind = ismember(sampleIDs,metaData(ismember(metaData(:,3),'CTRL')&ismember(metaData(:,4),'Liver')&ismember(metaData(:,2),'Chow')&ismember(metaData(:,6),'8w')&~contains(metaData(:,2),'Both'),1));
ChowKOind = ismember(sampleIDs,metaData(ismember(metaData(:,3),'KO')&ismember(metaData(:,4),'Liver')&ismember(metaData(:,2),'Chow')&ismember(metaData(:,6),'8w')&~contains(metaData(:,2),'Both'),1));
HSDCTRLind = ismember(sampleIDs,metaData(ismember(metaData(:,3),'CTRL')&ismember(metaData(:,4),'Liver')&ismember(metaData(:,2),'HSD')&ismember(metaData(:,6),'8w')&~contains(metaData(:,2),'Both'),1));
HSDKOind = ismember(sampleIDs,metaData(ismember(metaData(:,3),'KO')&ismember(metaData(:,4),'Liver')&ismember(metaData(:,2),'HSD')&ismember(metaData(:,6),'8w')&~contains(metaData(:,2),'Both'),1));

selGene = 'Fasn';
selGeneInd = ismember(genes,selGene);
boxplot([expMat(selGeneInd,ChowCTRLind)';expMat(selGeneInd,ChowKOind)';expMat(selGeneInd,HSDCTRLind)';expMat(selGeneInd,HSDKOind)';],[repmat({'Chow-CTRL'},sum(ChowCTRLind),1);repmat({'Chow-KO'},sum(ChowKOind),1);repmat({'HSD-CTRL'},sum(HSDCTRLind),1);repmat({'HSD-KO'},sum(HSDKOind),1);]);
xtickangle(50)
title(selGene)
ylabel('TPM')

%% Figure 2B & 3B & S5A & S5B
tissue = 'Liver';
% tissue = 'WAT';
% tissue = 'Muscle';
% tissue = 'Heart';
content = Cheng_readExcel('Table S6.xlsx',['HSDCTRL' tissue 'vsChowCTRL' tissue]);
HSD_DE.genes = strrep(content(2:end,1),'"','');
HSD_DE.lgFC = str2double(content(2:end,2));
HSD_DE.padj = str2double(content(2:end,4));

content = Cheng_readExcel('Table S6.xlsx',['HSDKO' tissue 'vsHSDCTRL' tissue]);
% content = Cheng_readExcel('C:\work\my papers\Under review\PKL KO project\For manuscript\For submission\Table S9.xlsx',['ChowKO' tissue 'vsChowCTRL' tissue]);
HSDKO_DE.genes = strrep(content(2:end,1),'"','');
HSDKO_DE.lgFC = str2double(content(2:end,2));
HSDKO_DE.padj = str2double(content(2:end,4));


HSD_UP = HSD_DE.genes(HSD_DE.lgFC>0&HSD_DE.padj<0.05);
HSD_DN = HSD_DE.genes(HSD_DE.lgFC<0&HSD_DE.padj<0.05);
HSDKO_UP = HSDKO_DE.genes(HSDKO_DE.lgFC>0&HSDKO_DE.padj<0.05);
HSDKO_DN = HSDKO_DE.genes(HSDKO_DE.lgFC<0&HSDKO_DE.padj<0.05);



overlap = length(intersect(HSD_UP,HSDKO_DN))+length(intersect(HSD_DN,HSDKO_UP));
HSDspecific = length(HSD_UP)+length(HSD_DN)-overlap;
HSDKOspecific = length(HSDKO_DN)+length(HSDKO_UP)-overlap;
Cheng_hyperGeometricTest(union(intersect(HSD_UP,HSDKO_DN),intersect(HSD_DN,HSDKO_UP)),intersect(HSD_DE.genes,HSDKO_DE.genes),union(HSD_UP,HSD_DN))

Cheng_printCellArrayToFile([[HSD_UP;HSD_DN],[repmat({'red'},size(HSD_UP));repmat({'blue'},size(HSD_DN))]],'HSDdiff_KEGG_liver.txt')
Cheng_printCellArrayToFile([[HSDKO_UP;HSDKO_DN],[repmat({'red'},size(HSDKO_UP));repmat({'blue'},size(HSDKO_DN))]],'HSDKOdiff_KEGG_liver.txt')
Cheng_printCellArrayToFile([intersect(HSD_UP,HSDKO_DN),repmat({'blue'},size(intersect(HSD_UP,HSDKO_DN)));intersect(HSD_DN,HSDKO_UP),repmat({'red'},size(intersect(HSD_DN,HSDKO_UP)))],'HSDoverlap_KEGG_liver_58.txt')

%% Figure 2D & 3B & S5C & S5D
tissue = 'Liver';
% tissue = 'WAT';
% tissue = 'Muscle';
% tissue = 'Heart';
content = Cheng_readExcel('Table S6.xlsx',['HSDCTRL' tissue 'vsChowCTRL' tissue]);
HSD_DE.genes = strrep(content(2:end,1),'"','');
HSD_DE.lgFC = str2double(content(2:end,2));
HSD_DE.padj = str2double(content(2:end,4));

content = Cheng_readExcel('Table S6.xlsx',['HSDKO' tissue 'vsHSDCTRL' tissue]);
HSDKO_DE.genes = strrep(content(2:end,1),'"','');
HSDKO_DE.lgFC = str2double(content(2:end,2));
HSDKO_DE.padj = str2double(content(2:end,4));


HSD_UP = HSD_DE.genes(HSD_DE.lgFC>0&HSD_DE.padj<0.05);
HSD_DN = HSD_DE.genes(HSD_DE.lgFC<0&HSD_DE.padj<0.05);
HSDKO_UP = HSDKO_DE.genes(HSDKO_DE.lgFC>0&HSDKO_DE.padj<0.05);
HSDKO_DN = HSDKO_DE.genes(HSDKO_DE.lgFC<0&HSDKO_DE.padj<0.05);


[A,B] = ismember(HSD_DE.genes,HSDKO_DE.genes);
lgFC1 = HSD_DE.lgFC(A);
lgFC2 = HSDKO_DE.lgFC(B(A));
inds = ~isnan(lgFC1)&~isnan(lgFC2);
padj1 = HSD_DE.padj(A);
padj2 = HSDKO_DE.padj(B(A));
% Only overlap
sz = 10*ones(size(padj1))+50*(padj1<0.05&padj2<0.05);
c = [0*(padj2<0.05)+1*(padj1<0.05&padj2<0.05),zeros(size(padj1)),zeros(size(padj1))];

rmInds = lgFC1>10|lgFC1<-10|lgFC2>10|lgFC2<-10;
inds = inds&~rmInds;
rho = corr(lgFC1(inds),lgFC2(inds),'type','spearman');
x = lgFC1(inds);
y = lgFC2(inds);
sz = sz(inds);
c = c(inds,:);
[~,inds] = sort(sz,'ascend');
scatter(x(inds),y(inds),sz(inds),c(inds,:),'filled');


title([tissue ' Spearman corr. = ' num2str(rho)])
xlabel(['lgFC-HSDvsChow-' tissue])
ylabel(['lgFC-HSDKOvsHSD-' tissue])
