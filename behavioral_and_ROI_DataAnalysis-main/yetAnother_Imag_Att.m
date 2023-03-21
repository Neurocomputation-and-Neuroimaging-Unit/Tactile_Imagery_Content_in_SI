% yet another skript for plots and statistics

path = 'C:\Users\saraw\Desktop\revision\goodImag\controlCusters\ConjConj_FWE_clust\';
cd(path)
data = dir('*.mat')
rois = 5;
MeansVTIQ = [2.75; 4.375; 1.25; 2.1875; 1.9375; 3.4375; 2.625; 3.375; 3.625; 3.375; 3.8125; 4.25; 3.0625; 3.5625; 3.625; 3.1875; 2.75; 3.875; 2.625; 3.375; 3.4375; 4.5625; 3.8125; 3.3125; 3.5; 3.4375; 3.875; 3.8125];
MeansVTIQ2 = MeansVTIQ';
MeansVTIQ2([3 4 13 19 21 24 27]) = []; % 21 % 3 13 19 26 % 3 4 13 19 21 24 27

exclude = []

all = [];

for r = 1


    % load
    
    thisData = load([path data(r).name]);
    thisData = thisData.gdata{1,1};
    for e = 1:length(exclude)
        del = exclude(e)-(e-1)
        thisData(del,:) = [];
    end
    size(thisData)

    % calculate & plot

%     meanData = [mean(thisData(:,1:3), 2) mean(thisData(:,4:6), 2) mean(thisData(:,9:10), 2)];
    meanData = [mean(thisData(:,1:3), 2) mean(thisData(:,4:6), 2)];
%     meanData = [thisData(:,1) thisData(:,2)];
%     meanData = [mean(thisData(:,1:3), 2) mean(thisData(:,4:6), 2) mean(thisData(:,7), 2)];

%     all = [all meanData];
%     stim = mean(thisData(:,1:3), 2);
%     imag = mean(thisData(:,4:6), 2);
%     att = mean(thisData(:,7), 2);
%     att = mean(thisData(:,9:10), 2);
    
%     diffData1 = meanData(:,1)-meanData(:,3);
%     diffData2 = meanData(:,2)-meanData(:,3);
    data_sem = std(meanData)/sqrt(27);

    meany = mean(meanData, 1);
% % 
%     f = figure
%     b = bar(meany,'FaceColor', 'flat');
%     hold on
%     errorbar(1:3, meany, data_sem, 'LineStyle', 'none')
%     set(gca, 'xtick', 1:3, 'xticklabel', {'Stim', 'Imag', 'Att'})
%     set(gca, 'ylim', [-1 7]);

    % statistics
% 
%     both = [MeansVTIQ2' meanData(:,2) meanData(:,3) diffData2];
%     corrplot(both)
%     hfig = gcf;
%     haxes = findobj(hfig, 'Type', 'Axes');
%     arrayfun(@(ax) ylim(ax, [-3 5.5]), haxes);

    both = [MeansVTIQ2' meanData(:,1), meanData(:,2)];
    corrplot(both)
    hfig = gcf;
    haxes = findobj(hfig, 'Type', 'Axes');

%     arrayfun(@(ax) ylim(ax, [-3 5.5]), haxes);
% % 
%     [rho1, pval1] = corrcoef(MeansVTIQ2, diffData1)
%     [rho2, pval2] = corrcoef(MeansVTIQ2, diffData2)

    [rho1, pval1] = corrcoef(MeansVTIQ2, meanData(:,1));
    stim = [rho1(1,2) pval1(1,2)]
    [rho2, pval2] = corrcoef(MeansVTIQ2, meanData(:,2));
    imag = [rho2(1,2) pval2(1,2)]
%     [rho1, pval1] = corrcoef(MeansVTIQ2, diffData2)

%     'Stim'
%     [h, p, ci, stats] = ttest(diffData1, 0, 'Tail', 'right')
%     'Imag'
%     [h, p, ci, stats] = ttest(diffData2, 0, 'Tail', 'right')
%     'Stim'
%     [h, p, ci, stats] = ttest(diffData1, 0, 'Tail', 'right')
%     'Imag'
%     [h, p, ci, stats] = ttest(diffData2, 0, 'Tail', 'right')

end

mean(all, 1)
