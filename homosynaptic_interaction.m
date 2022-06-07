addpath(genpath('/home/conghu/MatlabCodes/cNE_connectivity'))
addpath(genpath('/home/conghu/MatlabCodes/support/figureplot'))
data_path = '/data/congcong/rat_MGB_A1/3_singleunit/dmr';
pair_files = dir(fullfile(data_path, '*pairs.mat'));
fig_folder = '/data/congcong/rat_MGB_A1/figure/singleunit/connect/homosynaptic';
cd(data_path)

%% get efficacy of connected pairs
plot_opt = 1;
binsize = 0.5;
centers = -100:binsize:100;
DT = 50;
isi_min = 5;
for ii = 1:length(pair_files)
    
    filename = fullfile(data_path, pair_files(ii).name);
    load(filename, 'connected_pairs')
    if isempty(connected_pairs)
        continue
    end
%     if isfield(connected_pairs, 'efficacy_spon')
%         continue
%     end
    % load spk
    exp = pair_files(ii).name(1:13);
    file_MGB = dir(fullfile(data_path, [exp '*H31x64*split.mat']));
    data = load(fullfile(data_path, file_MGB.name), 'spk');
    spk_MGB = data.spk.spk;
    file_A1 = dir(fullfile(data_path, [exp '*H22x32*split.mat']));
    data = load(fullfile(data_path, file_A1.name), 'spk');
    spk_A1 = data.spk.spk;
    
    for jj = 1:length(connected_pairs)
        
            
        idx_MGB = connected_pairs(jj).idx_MGB;
        idx_A1 = connected_pairs(jj).idx_A1;
        % for spontaneous response------------------------------
        connected_pairs(jj).nspk_MGB_spon = length(spk_MGB(idx_MGB).spiketimes_spon);
        connected_pairs(jj).nspk_A1_spon = length(spk_A1(idx_A1).spiketimes_spon);
        
        spiketimes = spk_A1(idx_A1).spiketimes_spon;
        spiketimes_all = spiketimes;
        % get raster of all spikes
        reftimes = spk_MGB(idx_MGB).spiketimes_spon;
        reftimes_all = reftimes;
        dur_spon = max(max(spiketimes), max(reftimes));
        r_all = raster(spiketimes, reftimes, 40);
        connected_pairs(jj).raster_spon = r_all;
        
        % efficacy
        ccg = connected_pairs(jj).ccg_spon;
        baseline = connected_pairs(jj).baseline_spon;
        idx = 203:210;
        causal_spikes = (ccg(idx) -  baseline(idx)');
        efficacy = sum(causal_spikes(causal_spikes>0)) / length(reftimes);
        connected_pairs(jj).efficacy_spon = efficacy*100;

        % get raster of spike pairs
        [fspk, sspk] = spike_pairs(reftimes, 30, DT);
        [spk_pairs, r_pair, efficacy] = pair_efficacy(fspk, sspk, spiketimes, isi_min);
        connected_pairs(jj).spk_pair_spon = spk_pairs;
        connected_pairs(jj).raster_pair_spon = r_pair;
        connected_pairs(jj).efficacy_pair_spon = efficacy;
        
        % get isi
        isi_A1_spon = diff(spiketimes);
        isi_MGB_spon = diff(reftimes);
        
        % for dmr response ------------------------------------
        connected_pairs(jj).nspk_MGB_stim = length(spk_MGB(idx_MGB).spiketimes_dmr);
        connected_pairs(jj).nspk_A1_stim = length(spk_A1(idx_A1).spiketimes_dmr);
        
        spiketimes = spk_A1(idx_A1).spiketimes_dmr;
        spiketimes_all = [spiketimes_all; spiketimes+dur_spon+100]; % add 100ms between spon and stim response
        % get raster of all spikes
        reftimes = spk_MGB(idx_MGB).spiketimes_dmr;
        reftimes_all = [reftimes_all; reftimes+dur_spon+100];
        r_all = raster(spiketimes, reftimes, 40);
        connected_pairs(jj).raster_stim = r_all;
        
         % efficacy
        ccg = connected_pairs(jj).ccg_stim;
        baseline = connected_pairs(jj).baseline_stim;
        idx = 203:210;
        causal_spikes = (ccg(idx) -  baseline(idx)');
        efficacy = sum(causal_spikes(causal_spikes>0)) / length(reftimes);
        connected_pairs(jj).efficacy_stim = efficacy*100;
        
        % get raster of spike pairs
        [fspk, sspk] = spike_pairs(reftimes, 30, DT);
        [spk_pairs, r_pair, efficacy] = pair_efficacy(fspk, sspk, spiketimes, isi_min);
        connected_pairs(jj).spk_pair_stim = spk_pairs;
        connected_pairs(jj).raster_pair_stim = r_pair;
        connected_pairs(jj).efficacy_pair_stim = efficacy;
        
        % get isi
        isi_A1_stim = diff(spiketimes);
        isi_MGB_stim = diff(reftimes);
        
        % ------ for combined response ------------------------------------
        [fspk, sspk] = spike_pairs(reftimes_all, 30, DT);
        [spk_pairs, r_pair, efficacy] = pair_efficacy(fspk, sspk, spiketimes_all, isi_min);
        connected_pairs(jj).spk_pair_all = spk_pairs;
        connected_pairs(jj).raster_pair_all = r_pair;
        connected_pairs(jj).efficacy_pair_all = efficacy;
        if plot_opt
            figure('Renderer', 'painters', 'Position', [10 50 1200 800],'visible', 'on')
            
            % -------raster and ccg for all MGB spikes ---------------------------------
            xl = 'Time After MGB Spike\n (ms)';
            yl = 'MGB Spike #';
            t = 'All MGB Spikes';
            % raster plot (spon)
            subplot(2, 6, 1)
            r = connected_pairs(jj).raster_spon;
            nspk = connected_pairs(jj).nspk_MGB_spon;
            raster_format(r, nspk, t, xl, yl)
            % ccg
            subplot(4, 6, 13)
            ccg = connected_pairs(jj).ccg_spon;
            baseline = connected_pairs(jj).baseline_spon;
            efficacy = connected_pairs(jj).efficacy_spon;
            ccg_format(centers, ccg, baseline, xl, efficacy)
            ylabel(sprintf('A1 Neuron Response\n(spikes)'))
            
            % raster plot (stim)
            subplot(2, 6, 3)
            r = connected_pairs(jj).raster_stim;
            nspk = connected_pairs(jj).nspk_MGB_stim;
            raster_format(r, nspk, t, xl, yl)
            % ccg
            subplot(4, 6, 15)
            ccg = connected_pairs(jj).ccg_stim;
            baseline = connected_pairs(jj).baseline_stim;
            efficacy = connected_pairs(jj).efficacy_stim;
            ccg_format(centers, ccg, baseline, xl, efficacy)
            
             % raster plot (all)
            subplot(2, 6, 5)
            tmp = connected_pairs(jj).raster_stim;
            nspk_spon = connected_pairs(jj).nspk_MGB_spon;
            tmp(:,2) = tmp(:,2) + nspk_spon;
            r = [connected_pairs(jj).raster_spon; tmp];
            nspk = connected_pairs(jj).nspk_MGB_stim + nspk_spon;
            raster_format(r, nspk, t, xl, yl)
            plot([-40 40], connected_pairs(jj).nspk_MGB_spon*[1 1], 'k', 'linewidth', 1)
            % ccg
            subplot(4, 6, 17)
            ccg = connected_pairs(jj).ccg_stim ...
                + connected_pairs(jj).ccg_spon;
            baseline = connected_pairs(jj).baseline_stim ...
                + connected_pairs(jj).baseline_spon;
            causal_spikes = (ccg(idx) -  baseline(idx)');
            nspk = connected_pairs(jj).nspk_MGB_spon ...
                + connected_pairs(jj).nspk_MGB_stim;
            efficacy = sum(causal_spikes(causal_spikes>0)) / nspk * 100;
            connected_pairs(jj).efficacy_all = efficacy;
            ccg_format(centers, ccg, baseline, xl, efficacy)
            
            %----------raster and ccg for spike pairs -------------------------------------------
            yl = 'pair ISI (ms)';
            t = 'MGB Spike Pairs';
            % spon
            subplot(2, 6, 2)
            r = connected_pairs(jj).raster_pair_spon;
            nspk = length(connected_pairs(jj).spk_pair_spon);
            raster_format(r, nspk, t, xl, yl)
            raster_format2(r, nspk)
            %ccg
            subplot(4, 6, 14)
            ccg_format2(r)
            
            % dmr
            subplot(2, 6, 4)
            r = connected_pairs(jj).raster_pair_stim;
            nspk = length(connected_pairs(jj).spk_pair_stim);
            raster_format(r, nspk, t, xl, yl)
            raster_format2(r, nspk)
            %ccg
            subplot(4, 6, 16)
            ccg_format2(r)
            
            % all
            subplot(2, 6, 6)
            r = connected_pairs(jj).raster_pair_all;
            nspk = length(connected_pairs(jj).spk_pair_all);
            raster_format(r, nspk, t, xl, yl)
            raster_format2(r, nspk)
            
            %ccg
            subplot(4, 6, 18)
            ccg_format2(r)
            
            %----------efficacy vs ISI -------------------------------------------
            %spon
            subplot(4, 6, 19)
            efficacy_pair = connected_pairs(jj).efficacy_pair_spon;
            efficacy = connected_pairs(jj).efficacy_spon;
            plot_efficacy_isi(efficacy_pair, efficacy, 'spon');
            
            %stim
            subplot(4, 6, 20)
            efficacy_pair = connected_pairs(jj).efficacy_pair_stim;
            efficacy = connected_pairs(jj).efficacy_stim;
            plot_efficacy_isi(efficacy_pair, efficacy, 'dmr');
            
            %all
            subplot(4, 6, 21)
            efficacy_pair = connected_pairs(jj).efficacy_pair_all;
            efficacy = connected_pairs(jj).efficacy_all;
            plot_efficacy_isi(efficacy_pair, efficacy, 'all');
            
            %---------- ISI -------------------------------------------
            subplot(4, 6, 22)
            edges = 0:0.5:40;
            isi_spon = histcounts(isi_MGB_spon, edges);
            isi_stim = histcounts(isi_MGB_stim, edges);
            c =  (edges(1:end-1) + edges(2:end))/2;
            bar(c, isi_spon, 'facecolor', 0.4*[1 1 1], 'edgecolor', 0.4*[1 1 1])
            hold on
            plot(c, isi_stim, 'k')
            ymax = max([isi_spon, isi_stim])*1.1;
            ylim([0 ymax])
            text(10, ymax*0.8, 'MGB ISI')
            
            subplot(4, 6, 23)
            edges = 0:0.5:20;
            isi_spon = histcounts(isi_A1_spon, edges);
            isi_stim = histcounts(isi_A1_stim, edges);
            c =  (edges(1:end-1) + edges(2:end))/2;
            bar(c, isi_spon, 'facecolor', 0.4*[1 1 1], 'edgecolor', 0.4*[1 1 1])
            hold on
            plot(c, isi_stim, 'k')
            ymax = max([isi_spon, isi_stim])*1.1;
            ylim([0 ymax])
            text(10, ymax*0.8, 'A1 ISI')
            
            subplot(4, 6, 24)
            text(0.1, 0.9, sprintf('A1 nspk (spon): %d', length(isi_A1_spon)+1))
            text(0.1, 0.8, sprintf('A1 nspk (dmr): %d', length(isi_A1_stim)+1))
            text(0.1, 0.7, sprintf('MGB nspk (spon): %d', length(isi_MGB_spon)+1))
            text(0.1, 0.6, sprintf('MGB nspk (dmr): %d', length(isi_MGB_stim)+1))
            text(0.1, 0.5, sprintf('MGB npairs (spon): %d', ...
                size(connected_pairs(jj).spk_pair_spon, 1)))
            text(0.1, 0.4, sprintf('MGB npairs (dmr): %d', ...
                size(connected_pairs(jj).spk_pair_stim, 1)))
            xlim([0, 1])
            ylim([0 1])
            axis off
            box off

            
            saveas(gcf, fullfile(fig_folder, sprintf('%s-tha%d-crtx%d.jpg', exp, idx_MGB, idx_A1)))
            close
        end
        
    end
    
    save(filename, 'connected_pairs')
end
%% ccg aligned to first & second spike
t = 103:110;
centers = -50:0.5:50;
xl = 'Time After MGB Spike';

for ii = 1:length(pair_files)
    filename = fullfile(data_path, pair_files(ii).name);
    load(filename, 'connected_pairs')
       

    if isempty(connected_pairs)
        continue
    end

    % load spk
    exp = pair_files(ii).name(1:13);
    file_A1 = dir(fullfile(data_path, [exp '*H22x32*split.mat']));
    data = load(fullfile(data_path, file_A1.name), 'spk');
    spk_A1 = data.spk.spk;
    
    for jj = 1:length(connected_pairs)
        
        figure('Renderer', 'painters', 'Position', [10 50 400 600],'visible', 'on')
        
        idx_MGB = connected_pairs(jj).idx_MGB;
        idx_A1 = connected_pairs(jj).idx_A1;
        
        efficacy_all = zeros(1, 6);
        ccg_all = zeros(4, 201);
        baseline_all =  zeros(4, 201);
        nspk_all = 0;
        
        conditions = {'spon', 'stim'};
        c = 1;
        for kk = 1:length(conditions)
            condition = conditions{kk};
            
            spk_pairs = eval(sprintf('connected_pairs(jj).spk_pair_%s', condition));
            if strcmp(condition, 'stim')
                spiketimes = spk_A1(idx_A1).spiketimes_dmr;
            else
            spiketimes = eval(sprintf('spk_A1(idx_A1).spiketimes_%s', condition));
            end
            nspk = length(spk_pairs);
            nspk_all = nspk_all+nspk;
            % first spike
            subplot(3,2,kk*2-1)
            first_spikes = spk_pairs(:, 1);
            [efficacy, ccg, baseline] = ccg_format3(spiketimes, first_spikes, nspk);
            efficacy_all(c) = efficacy;
            ccg_all(c,:) = ccg;
            baseline_all(c,:) = baseline;
            title(sprintf('%s 1st', condition))
            c = c+1;
            
            % second spike
            subplot(3,2,kk*2)
            second_spikes = spk_pairs(:, 2);
            [efficacy, ccg, baseline] = ccg_format3(spiketimes, second_spikes, nspk);
            efficacy_all(c) = efficacy;
            ccg_all(c,:) = ccg;
            baseline_all(c,:) = baseline;
            title(sprintf('%s 2nd', condition))
            c = c+1;
        end
        
        % combined
        %first spike
        subplot(3, 2, 5)
        ccg = sum(ccg_all([1,3], :));
        baseline = sum(baseline_all([1,3], :));
        ccg_baseline = ccg - baseline;
        causal_spikes = ccg_baseline(t);
        efficacy_all(5) = sum(causal_spikes(causal_spikes > 0))/ nspk_all*100;
        ccg_format(centers, ccg, baseline, xl, efficacy_all(5))
        title('all 1st')
        
        subplot(3, 2, 6)
        ccg = sum(ccg_all([2,4], :));
        baseline = sum(baseline_all([2,4], :));
        ccg_baseline = ccg - baseline;
        causal_spikes = ccg_baseline(t);
        efficacy_all(6) = sum(causal_spikes(causal_spikes > 0))/ nspk_all*100;
        ccg_format(centers, ccg, baseline, xl, efficacy_all(6))
        title('all 2nd')
        
        % save results
        connected_pairs(jj).ccg_pair = ccg_all;
        connected_pairs(jj).baseline_pair = baseline_all;
        connected_pairs(jj).efficacy_pair = efficacy_all;
        
        saveas(gcf, fullfile(fig_folder, sprintf('ccg-%s-tha%d-crtx%d.jpg', exp, idx_MGB, idx_A1)))
        close all
    end
    save(filename, 'connected_pairs')
end
%% summary plot for efficacy
efficacy_all = [];
efficacy1 = [];
efficacy2 = [];
for ii = 1:length(pair_files)
    filename = fullfile(data_path, pair_files(ii).name);
    load(filename, 'connected_pairs')
    for jj = 1:length(connected_pairs)
        efficacy1 = [efficacy1, connected_pairs(jj).efficacy_pair(5)];
        efficacy2 = [efficacy2, connected_pairs(jj).efficacy_pair(6)];
        efficacy_all = [efficacy_all, connected_pairs(jj).efficacy_all];
    end
end

figure('Position', [10 50 250 250])
scatter(efficacy_all, efficacy1, 'k')
hold on
plot([0 20], [0, 20], 'k')
xlim([0 20])
ylim([0 20])
xlabel('efficacy of all spikes')

figure('Position', [10 50 250 250])
scatter(efficacy_all, efficacy2, 'k')
hold on
plot([0 20], [0, 20], 'k')
xlim([0 20])
ylim([0 20])
xlabel('efficacy of all spikes')
ylabel('efficacy of 2nd spike in spike pairs')

figure('Position', [10 50 250 250])
scatter(efficacy1, efficacy2, 'k')
hold on
plot([0 20], [0, 20], 'k')
xlim([0 20])
ylim([0 20])
xlabel('efficacy of 1st spike in spike pairs')
ylabel('efficacy of 2nd spike in spike pairs')
%% plot of # of A1 spikes after each MGB spike
cd('/data/congcong/rat_MGB_A1/3_singleunit/dmr')
files = dir('*pairs.mat');
isi = 1:30;
nA1spk1 = nan(1, 30);
nA1spk2 = nan(1, 30);
nMGBspk = nan(1, 30);
for ii = 1:length(files)
    load(files(ii).name, 'connected_pairs')
    exp = files(ii).name(1:13);
    for jj = 1:length(connected_pairs)
        
        idx_MGB = connected_pairs(jj).idx_MGB;
        idx_A1 = connected_pairs(jj).idx_A1;
        spk_pairs = connected_pairs(jj).raster_pair_all;
        
        for kk = 1:length(isi) 
            idx = spk_pairs(:, 3) > kk-1 & spk_pairs(:, 3) <= kk;
            if sum(idx) == 0
                continue
            end
            spk = spk_pairs(idx, :);
            spktimes2 = spk(:, 1) - spk(:, 3);
            nA1spk1(kk) = sum(spk(:,1) > 0 & spk(:,1) <= 5);
            nA1spk2(kk) = sum(spktimes2 > 0 & spktimes2 <= 5);
            nMGBspk(kk) = length(spk);
            
        end
        figure('visible', 'off')
        plot(1:30, nA1spk1./nMGBspk)
        hold on
        plot(1:30, nA1spk2./nMGBspk)
        xlabel('ISI (ms)')
        ylabel('# of A1 spikes / # of MGB spike pairs')
        legend({'1st spike', '2nd spike'})
        saveas(gcf, fullfile(fig_folder, sprintf('nspk-%s-tha%d-crtx%d.jpg', exp, idx_MGB, idx_A1)))
        close

    end
end


%% helper functions
function plot_efficacy_isi(efficacy_pair, efficacy, l)
plot(1:30, efficacy_pair(:,1), 'k--')
hold on
plot(1:30, efficacy_pair(:,2), 'k')
plot([1, 30], efficacy*[1, 1], 'k:')
ymax = max([max(max(efficacy_pair)), efficacy])*1.1;
ylim([0 ymax])
xlabel(sprintf('ISI (MGB2-MGB1) \n (ms)'))
ylabel('Efficacy %')
text(20, ymax*0.8, l)
end

function raster_format(raster, ymax, t, xl, yl)
scatter(raster(:,1), raster(:,2), 1.5, ...
                'markeredgecolor', 0.4*[1 1 1])
hold on
plot([0 0], [1 ymax], 'b--', 'linewidth', 1)
box on
xlim([-40 40])
ylim([1, ymax])
yticks(ylim)
ylabel(yl)
xlabel(sprintf(xl))
yh = get(gca,'ylabel');
p = get(yh,'position');
p(1) = p(1)*0.8;
set(yh, 'position', p)
title(t)
end

function raster_format2(r, nspk)
isi = unique(ceil(r(:,3)));
yt = zeros(1, length(isi));
for kk = 1:length(isi)
    idx = find(r(:,3) <= isi(kk), 1, 'last');
    plot([-40, 40], r(idx, 2)*[1, 1], 'k')
    yt(kk) = r(idx, 2);
end
idx = find(diff(yt) == 0);
yt(idx) = [];
isi(idx) = [];
yt(end) = nspk;
yticks(yt)
yticklabels(isi)
yt = [0 yt];
for kk = 1:length(yt)-1
    plot(isi(kk)*[1 1], [yt(kk), yt(kk+1)], 'b')
end

end

function ccg_format(centers, ccg, baseline, xl, efficacy)

bar(centers, ccg, 1, 'k')
hold on
ymax = max(ccg)*1.1;
plot([0, 0], [0, ymax], 'b--', 'linewidth', 1)
plot(centers, baseline, 'b')
xlim([-40 40])
ylim([0, ymax])
xlabel(sprintf(xl))
text(15, 0.8*ymax, sprintf('%.2f', efficacy))
end

function ccg_format2(r)
edges =  -40:0.5:40;
ccg = histcounts(r(:,1), edges);
c = (edges(1:end-1) + edges(2:end))/2;
bar(c, ccg, 1, 'k')
hold on
ymax = max(ccg)*1.1;
plot([0, 0], [0, ymax], 'b--', 'linewidth', 1)
ylim([0 ymax])
xlabel(sprintf('Time After 1st MGB Spike\n (ms)'))
end

function [efficacy, ccg, baseline] = ccg_format3(spiketimes, reftimes, nspk)
t = 103:110;
centers = -50:0.5:50;
xl = 'Time After MGB Spike';

[baseline, ccg] = get_baseline(spiketimes, reftimes, 0.5);
ccg_baseline = ccg - baseline;
causal_spikes = ccg_baseline(t);
efficacy = sum(causal_spikes(causal_spikes > 0))/ nspk*100;
ccg_format(centers, ccg, baseline, xl, efficacy)
end

function [spk_pairs, r_pair, efficacy] = pair_efficacy(fspk, sspk, spiketimes,t)
spk_pairs = [fspk, sspk]; % 1-first spike; 2-second spike
spk_pairs(:,3) = spk_pairs(:,2) - spk_pairs(:,1); % 3-isi
spk_pairs(spk_pairs(:,3) < t-1, :) = [];
spk_pairs(:,4) = 1:length(spk_pairs);% 4-current spike order (based on first spike time)
spk_pairs = sortrows(spk_pairs, 3); 
spk_pairs(:,5) = 1:length(spk_pairs);% 5-order based isi
spk_pairs = sortrows(spk_pairs, 4); % reorder spike times
r_pair = raster(spiketimes, spk_pairs(:,1), 40);
% 1-cortical spike relative time; 
% 2-corresponding MGB spike label as in column4
r_pair(:, 3) = spk_pairs(r_pair(:, 2),3); % 3-get isi for the spike pair
r_pair(:, 4) = spk_pairs(r_pair(:, 2),1); % 4- first spike
r_pair(:, 5) = spk_pairs(r_pair(:, 2),2); % 5 - second spike
r_pair(:, 2) = spk_pairs(r_pair(:, 2),5); % 2 - order based on isi in spk_pair

r_pair = sortrows(r_pair, 3);

efficacy = nan(30, 2);
for ii = t: 30
    idx = r_pair(:, 3) > ii-1 & r_pair(:, 3) <= ii;
    first_spikes = r_pair(idx, 4);
    
    if length(first_spikes) > 10
        
        nspk = length(first_spikes);
        second_spikes = r_pair(idx, 5);
        
        [baseline, ccg] = get_baseline(spiketimes, first_spikes, 0.5);
        ccg_baseline = ccg - baseline;
        t = 103:110;
        causal_spikes = ccg_baseline(t);
        efficacy(ii,1) = sum(causal_spikes(causal_spikes > 0))/ nspk*100;
        
        [baseline, ccg] = get_baseline(spiketimes, second_spikes, 0.5);
        ccg_baseline = ccg - baseline;
        causal_spikes = ccg_baseline(t);
        efficacy(ii,2) = sum(causal_spikes(causal_spikes > 0))/ nspk*100;
        
    end
end

end