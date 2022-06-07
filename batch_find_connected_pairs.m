addpath(genpath('/home/conghu/MatlabCodes/cNE_connectivity'))
addpath(genpath('/home/conghu/MatlabCodes/support/figureplot'))

data_path = '/data/congcong/rat_MGB_A1/3_singleunit/dmr';
figfolder = '/data/congcong/rat_MGB_A1/figure/singleunit/connect/dmr_spon';
if ~isfolder(figfolder)
    mkdir(figfolder)
end

crtxfiles = dir(fullfile(data_path, '*H22x32*split.mat'));
crtx_exp = cellfun(@(x) x(1:13), {crtxfiles.name}, 'UniformOutput', false);
binsize = 0.5;
bandwidth = [75, Inf];
hw = 100;
nfiles = length(crtxfiles);

%% find connected pairs
for ii = 1:length(crtxfiles)
    
    exp = crtxfiles(ii).name(1:13);
    site = regexp(crtxfiles(ii).name, '(?<=site)\d{1,2}', 'match', 'once');
    outfile = sprintf('%s-site%s-connected_pairs.mat', ...
        exp,  site);
    
    if exist(fullfile(data_path,outfile), 'file')
        fprintf('(%d/%d)Already processed %s\n', ii, nfiles, crtxfiles(ii).name)
        continue
    end
    fprintf('(%d/%d)processing %s\n', ii, nfiles, crtxfiles(ii).name)
    
    data = load(fullfile(data_path, crtxfiles(ii).name), ...
        'spk', 'trigger', 'spktrain', 'spktrain_spon');
    spk_A1 = data.spk;
    fs = spk_A1.fs;
    trigger_A1 = data.trigger;
    locators_A1_dmr = data.spktrain;
    locators_A1_spon = data.spktrain_spon;
    
    MGBfile = dir(fullfile(data_path, ...
        sprintf('%s*-H31x64-*split.mat', spk_A1.exp)));
    data = load(fullfile(MGBfile.folder, MGBfile.name), ...
        'spk', 'trigger', 'spktrain', 'spktrain_spon');
    spk_MGB = data.spk;
    trigger_MGB = data.trigger;
    locators_MGB_dmr = data.spktrain;
    locators_MGB_spon = data.spktrain_spon;
    
    if any(trigger_MGB ~= trigger_A1)
        fprintf('%s', MGBfile.name)
        warning('Triggers don''t match. Something''s wrong!')
        continue
    end
    
    % make the spon spktrain match length by appending zeros at the end
    if size(locators_MGB_spon, 2) ~= size(locators_A1_spon, 2)
        lt = length(locators_MGB_spon);
        lc = length(locators_A1_spon);
        l = max(lt, lc);
        locators_MGB_spon = [locators_MGB_spon ...
            zeros(size(locators_MGB_spon, 1), l-lt)];
        locators_A1_spon = [locators_A1_spon ...
            zeros(size(locators_A1_spon, 1), l-lc)];
        spktrain_spon = locators_A1_spon;
        parsave(crtxfiles(ii).name, spktrain_spon, 'spktrain_spon')
        spktrain_spon = locators_MGB_spon;
        parsave(MGBfile.name, spktrain_spon, 'spktrain_spon')
    end
    
    connected_pairs = find_connected_pairs(...
        locators_MGB_spon, locators_A1_spon, ...
        locators_MGB_dmr, locators_A1_dmr, ...
        hw, binsize, [spk_MGB.spk.unit], [spk_A1.spk.unit], bandwidth, ...
        {'Gaussian_convolve'});
    parsave(fullfile(data_path,outfile), connected_pairs, 'connected_pairs')
    
end
% plot connected pairs ccg
cd(data_path)
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
        ccg=ccg(:);
        baseline = flip(connected_pairs(jj).baseline_stim);
        baseline = baseline(:);
        bar(centers, ccg, 1, 'k')
        hold on
        centers2 = -4.5:0.5:-1;
        bar(centers2, ccg(192:199), 1, 'edgecolor', 'r', 'facecolor', 'r')
        bar(centers2, min([ccg(192:199),baseline(192:199)],[],2), 1, 'k')
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
        plot([5 5],[min(ylim) max(ylim)],'b--')
        ylabel('# of MGB spikes')
        xlabel('time from A1 spike (ms)')
        set(gca,'color','none','box','off')

        
        % spon raw
        subplot(222)
        hold on
        ccg = flip(connected_pairs(jj).ccg_spon);
        ccg = ccg(:);
        baseline = flip(connected_pairs(jj).baseline_spon);
        baseline = baseline(:);
        bar(centers, ccg, 1, 'k')
        bar(centers2, ccg(192:199), 1, 'edgecolor', 'r', 'facecolor', 'r')
        bar(centers2, min([ccg(192:199),baseline(192:199)], [],2), 1, 'k')
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
        plot([5 5],[min(ylim) max(ylim)],'b--')
        ylabel('# of MGB spikes')
        xlabel('time from A1 spike (ms)')
        set(gca,'color','none','box','off')

        saveas(gcf, fullfile(figfolder, ...
            sprintf('%s-tha%d-crtx%d.jpg', ...
            pairfiles(ii).name(1:13), idx_MGB, idx_A1)))
        close all
    end
end