%clustergram code
%dat_ref contains the refined profile data, class_label contains class
%labels for cases
%%z-score continuous analytes
dat_z = zscore(round(table2array(dat_ref), 2))';
dat_all = [dat_z];

%%make clustergram of people + metabolites
cgzdat = clustergram(dat_all, 'LabelsWithMarkers', true, 'Cluster', 3);
cgzdat.RowLabels = dat_ref.Properties.VariableNames;
cgzdat.ColumnLabels = class_label;
cgzdat.RowPDist = 'spearman';
cgzdat.ColumnPDist = 'spearman';
%cgzdat.LogTrans = 0;
cgzdat.Colormap = redgreencmap(200);
cgzdat.DisplayRange = [3];
cgzdat.Dendrogram = [0.78 1.1];

%label columns with a color
colLabels = class_label';
colColors = cellfun(@str2num, cellstr(renamecats(categorical(class_label), {'Control', 'WTC-LI'}, {'[0 0 178]/255', '[229 26 26]/255'})), 'UniformOutput', false)';
cgzdat.LabelsWithMarkers = 1;
cgzdat.ColumnLabelsColor = struct('Labels', colLabels, 'Colors', colColors);

markzdat = struct('GroupNumber', {26, 23, 22, 25, 27}, 'Annotation', {'1', '2', '3', '4', '5'}, 'Color', {[217 136 128]/255, [195 155 211]/255, [127 179 213]/255, [118 215 196]/255, [247 220 111]/255});
cgzdat.RowGroupMarker = markzdat;

markcas = struct('GroupNumber', {28, 26}, 'Annotation', {'Control', 'WTC-LI'}, 'Color', {[0 0 178]/255, [229 26 26]/255});
cgzdat.ColumnGroupMarker = markcas;

cgzdat.RowLabels = dat_ref.Properties.VariableNames;