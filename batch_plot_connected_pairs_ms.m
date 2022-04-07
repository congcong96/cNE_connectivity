addpath(genpath('/home/conghu/MatlabCodes/cNE_connectivity'))
addpath(genpath('/home/conghu/MatlabCodes/support/figureplot'))
data_path = '/data/congcong/rat_MGB_A1/3_singleunit/dmr';
figfolder = '/data/congcong/rat_MGB_A1/figure/singleunit/connect/dmr_spon';
%% find connected pairs
pairfiles = dir(fullfile(data_path, '*-connected_pairs.mat'));
binsize = 0.5;
hw = 100;
centers = -hw:binsize:hw;
colors = cbrewer('qual', 'Paired', 12);
for ii = 1:length(pairfiles)
    
    load(fullfile(data_path, pairfiles(ii).name), 'connected_pairs')
    if isempty(connected_pairs)
        continue
    end
    
    %% plot connected pairs

    for jj = 1:length(connected_pairs)
        idx_MGB = connected_pairs(jj).idx_MGB;
        idx_A1 = connected_pairs(jj).idx_A1;
        
        figure('Renderer', 'painters', 'Position', [10 50 800 400], 'visible', 'off');
        
        % dmr raw
        subplot(221)
        ccg = flip(connected_pairs(jj).ccg_stim);
        baseline = flip(connected_pairs(jj).baseline_stim);
        bar(centers, ccg, 1, 'k')
        hold on
        centers2 = -4.5:0.5:-1;
        bar(centers2, ccg(192:199), 1, 'edgecolor', 'r', 'facecolor', 'r')
        bar(centers2, baseline(192:199), 1, 'k')
        plot([0 0],[min(ylim) max(ylim)],'r--')
        thresh = flip(connected_pairs(jj).high_bound_stim);
        plot(centers, baseline, 'b', 'linewidth', 1.5)
        plot(centers, thresh, 'b:', 'linewidth', 1)
        ylabel('# of MGB spikes')
        title(sprintf('stim sig-%d', connected_pairs(jj).sig_stim))
        set(gca,'color','none','box','off')

        
        subplot(223)
        hold on
        bar(centers, ccg(:)-baseline(:), 1, 'k')
        plot([0 0],[min(ylim) max(ylim)],'r--')
        ylabel('# of MGB spikes')
        xlabel('time from A1 spike (ms)')
        set(gca,'color','none','box','off')

        
        % spon raw
        subplot(222)
        hold on
        ccg = flip(connected_pairs(jj).ccg_spon);
        baseline = flip(connected_pairs(jj).baseline_spon);
        bar(centers, ccg, 1, 'k')
        bar(centers2, ccg(192:199), 1, 'edgecolor', 'r', 'facecolor', 'r')
        bar(centers2, baseline(192:199), 1, 'k')
        plot([0 0],[min(ylim) max(ylim)],'r--')
        thresh = flip(connected_pairs(jj).high_bound_spon);
        plot(centers, baseline, 'b', 'linewidth', 1.5)
        plot(centers, thresh, 'b:', 'linewidth', 1)
        ylabel('# of MGB spikes')
        title(sprintf('spon sig-%d', connected_pairs(jj).sig_spon))
        set(gca,'color','none','box','off')

        subplot(224)
        hold on
        bar(centers, ccg(:)-baseline(:), 1, 'k')
        plot([0 0],[min(ylim) max(ylim)],'r--')
        ylabel('# of MGB spikes')
        xlabel('time from A1 spike (ms)')
        set(gca,'color','none','box','off')

        saveas(gcf, fullfile(figfolder, ...
            sprintf('%s-tha%d-crtx%d.jpg', ...
            pairfiles(ii).name(1:13), idx_MGB, idx_A1)))
        close all
    end
end