data_path = '/data/congcong/rat_MGB_A1/3_singleunit/dmr';
pair_files = dir(fullfile(data_path, '*pairs.mat'));
fig_folder = '/data/congcong/rat_MGB_A1/figure/singleunit/connect/homosynaptic';
cd(data_path)

%% get efficacy of connected pairs
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
        
        spktimes_spon_MGB = spk_MGB(idx_MGB).spiketimes_spon;
        spktimes_stim_MGB = spk_MGB(idx_MGB).spiketimes_dmr;
        spktimes_spon_A1 = spk_A1(idx_A1).spiketimes_spon;
        spktimes_stim_A1 = spk_A1(idx_A1).spiketimes_dmr;
        
        spk_pairs_spon = connected_pairs(jj).spk_pair_spon;
        spk_pairs_stim = connected_pairs(jj).spk_pair_stim;
        
        figure('Renderer', 'painters', 'Position', [50 50 1000 750],'visible', 'off')
        
        % -------raster for single spikes and spike pairs ---------------------------------
        xl = 'Time After MGB Spike\n (ms)';
        
        % raster plot - single spikes
        subplot('position', [.05, .08, .18, .85])
        r_spon = connected_pairs(jj).raster_spon;
        [r_spon, nspk_spon, reftimes_spon] = get_single_spk(...
            r_spon, spktimes_spon_MGB, DT);
        r_stim = connected_pairs(jj).raster_stim;
        [r_stim, nspk_stim, reftimes_stim] = get_single_spk(...
            r_stim, spktimes_stim_MGB, DT);
        connected_pairs(jj).raster_spon_single = r_spon;
        connected_pairs(jj).raster_stim_single = r_stim;
        connected_pairs(jj).nspk_MGB_spon_single = nspk_spon;
        connected_pairs(jj).nspk_MGB_stim_single = nspk_stim;
        r_stim(:, 2) = r_stim(:, 2) + nspk_spon;
        r = [r_spon; r_stim];
        nspk = nspk_stim + nspk_spon;
        raster_format(r, nspk, 'MGB single spikes', xl, 'Spike #')
        plot([-40 40], nspk_spon*[1 1], 'k', 'linewidth', 1)
        
        % raster plot - spike pairs
        fig =  subplot('position', [.27, .08, .18, .85]);
        r = connected_pairs(jj).raster_pair_all;
        nspk = length(connected_pairs(jj).spk_pair_all);
        raster_format(r, nspk, 'MGB spike pairs', xl, 'Spike Pair #')
        raster_format2(r, nspk, fig)

        % -------ccg for single spikes and spike pairs ---------------------------------
        subplot('position', [.55, .68, .18, .24]);
        dur = max([spktimes_spon_MGB;spktimes_spon_A1]);
        reftimes_single = [reftimes_spon; reftimes_stim + dur + 100];
        spiketimes_all = [spktimes_spon_A1; spktimes_stim_A1 + dur + 100];
        nspk = length(reftimes_single);
        [efficacy, ccg, baseline] = ccg_format3(spiketimes_all, reftimes_single, nspk);
        connected_pairs(jj).efficacy_single_all = efficacy;
        connected_pairs(jj).ccg_single_all = ccg;
        connected_pairs(jj).baseline_single_all = baseline;
        xlabel('')
        title('Single Spike')
        
        subplot('position', [.55, .38, .18, .24]);
        reftimes_1st = [spk_pairs_spon(:,1); ...
            spk_pairs_stim(:,1) + dur + 100];
        nspk = length(reftimes_1st);
        [efficacy, ccg, baseline] = ccg_format3(...
            spiketimes_all, reftimes_1st, nspk);
        connected_pairs(jj).efficacy_1st_all = efficacy;
        connected_pairs(jj).ccg_1st_all = ccg;
        connected_pairs(jj).baseline_1st_all = baseline;
        ylabel('A1 Neuron Response (spikes)')
        title('1st Spike')
        xlabel('')
        
        subplot('position', [.55, .08, .18, .24]);
        reftimes_2nd = [spk_pairs_spon(:,2); ...
            spk_pairs_stim(:,2) + dur + 100];
        nspk = length(reftimes_2nd);
        [efficacy, ccg, baseline] = ccg_format3(...
            spiketimes_all, reftimes_2nd, nspk);
        connected_pairs(jj).efficacy_2nd_all = efficacy;
        connected_pairs(jj).ccg_2nd_all = ccg;
        connected_pairs(jj).baseline_2nd_all = baseline;
        xlabel('Time after MGB spike (ms)')
        title('2nd Spike')
        
        % -------# of A1 spikes folowing spike pairs-----------------------
        isi = 1:30;
        nA1spk1 = nan(1, 30);
        nA1spk2 = nan(1, 30);
        nMGBspkpair = nan(1, 30);
        spk_pairs = connected_pairs(jj).raster_pair_all;
        spk_single = [connected_pairs(jj).raster_spon_single;...
            connected_pairs(jj).raster_stim_single];
        for kk = 1:length(isi)
            
            idx = spk_pairs(:, 3) > kk-1 & spk_pairs(:, 3) <= kk;
            if sum(idx) ~= 0
               spk = spk_pairs(idx, :);
            spktimes2 = spk(:, 1) - spk(:, 3);
            nA1spk1(kk) = sum(spk(:,1) > 0 & spk(:,1) <= 5);
            nA1spk2(kk) = sum(spktimes2 > 0 & spktimes2 <= 5);
            nMGBspkpair(kk) = length(spk);
            end
        end
        nA1spk = sum(spk_single(:,1) > 0 & spk_single(:,1) <= 5);
        nMGBsingle = connected_pairs(jj).nspk_MGB_spon_single...
            + connected_pairs(jj).nspk_MGB_stim_single;
        
        subplot('position', [.8, .68, .18, .24]);
        nA1spk1 =  nA1spk1./nMGBspkpair;
        nA1spk1 = smooth(nA1spk1, 3, 'moving');
        nA1spk2 =  nA1spk2./nMGBspkpair;
        nA1spk2 = smooth(nA1spk2, 3, 'moving');
        nA1spk =  nA1spk/nMGBsingle;
        plot(1:30, nA1spk1, 'k--')
        hold on
        plot(1:30, nA1spk2, 'k-')
        xlabel('ISI (ms)')
        ylabel(sprintf('# of A1 spikes \n/ # of MGB spikes'))
        plot([1 30], nA1spk * [1 1], 'k:')
        legend({'1st', '2nd', 'single'})
        ymax = max([nA1spk1; nA1spk2; nA1spk]);
        ylim([0, ymax*1.2])
        
        %---------- ISI -------------------------------------------
        isi_MGB_spon = diff(spktimes_spon_MGB);
        isi_MGB_stim = diff(spktimes_stim_MGB);
        isi_A1_spon = diff(spktimes_spon_A1);
        isi_A1_stim = diff(spktimes_stim_A1);
        
        subplot('position', [.8, .38, .18, .24]);
        edges = 0:0.5:40;
        isi_spon = histcounts(isi_MGB_spon, edges);
        isi_stim = histcounts(isi_MGB_stim, edges);
        c =  (edges(1:end-1) + edges(2:end))/2;
        bar(c, isi_spon, 'facecolor', 0.4*[1 1 1], 'edgecolor', 0.4*[1 1 1])
        hold on
        plot(c, isi_stim, 'k')
        ymax = max([isi_spon, isi_stim])*1.1;
        ylim([0 ymax])
        xlabel('MGB ISI (ms)')
        ylabel('# of spikes')
        
        subplot('position', [.8, .08, .18, .24]);
        edges = 0:0.5:40;
        isi_spon = histcounts(isi_A1_spon, edges);
        isi_stim = histcounts(isi_A1_stim, edges);
        c =  (edges(1:end-1) + edges(2:end))/2;
        bar(c, isi_spon, 'facecolor', 0.4*[1 1 1], 'edgecolor', 0.4*[1 1 1])
        hold on
        plot(c, isi_stim, 'k')
        ymax = max([isi_spon, isi_stim])*1.1;
        ylim([0 ymax])
        xlabel('A1 ISI (ms)')
        ylabel('# of spikes')

        saveas(gcf, fullfile(fig_folder, sprintf('homosynaptic-%s-tha%d-crtx%d.jpg', exp, idx_MGB, idx_A1)))
        close
        
    end
    
    save(filename, 'connected_pairs')
end
%% summary plot of efficacy
efficacy_single = [];
efficacy1 = [];
efficacy2 = [];
for ii = 1:length(pair_files)
    filename = fullfile(data_path, pair_files(ii).name);
    load(filename, 'connected_pairs')
    for jj = 1:length(connected_pairs)
        efficacy1 = [efficacy1, connected_pairs(jj).efficacy_1st_all];
        efficacy2 = [efficacy2, connected_pairs(jj).efficacy_2nd_all];
        efficacy_single = [efficacy_single, ...
            connected_pairs(jj).efficacy_single_all];
    end
end

figure('Position', [10 50 250 250])
scatter(efficacy_single, efficacy1, 'k')
hold on
plot([0 20], [0, 20], 'k')
xlim([0 20])
ylim([0 20])
xlabel('efficacy of single spikes')
ylabel('efficacy of 1st spikes')
p = signrank(efficacy1, efficacy_single);
if p< 0.05
    title('p=%.3f', p)
end
saveas(gcf, fullfile(fig_folder, 'efficacy_1st_vs_single.jpg'))
close

figure('Position', [10 50 250 250])
scatter(efficacy_single, efficacy2, 'k')
hold on
plot([0 20], [0, 20], 'k')
xlim([0 20])
ylim([0 20])
xlabel('efficacy of single spikes')
ylabel('efficacy of 2nd spikes')
p = signrank(efficacy_single, efficacy2);
if p< 0.05
    title('p=%.3f', p)
end
saveas(gcf, fullfile(fig_folder, 'efficacy_2nd_vs_single.jpg'))
close

figure('Position', [10 50 250 250])
scatter(efficacy2, efficacy1, 'k')
hold on
plot([0 20], [0, 20], 'k')
xlim([0 20])
ylim([0 20])
p = signrank(efficacy1, efficacy2);
if p< 0.05
    title('p=%.3f', p)
end
xlabel('efficacy of 1st spikes')
ylabel('efficacy of 2nd spikes')
saveas(gcf, fullfile(fig_folder, 'efficacy_2nd_vs_1st.jpg'))
close 

idx = efficacy1 > efficacy_single;
figure('Position', [10 50 250 250])
scatter(efficacy_single(idx), efficacy1(idx))
hold on
scatter(efficacy_single(~idx), efficacy1(~idx))
plot([0 20], [0, 20], 'k')
xlim([0 20])
ylim([0 20])
xlabel('efficacy of single spikes')
ylabel('efficacy of 1st spikes')
saveas(gcf, fullfile(fig_folder, 'efficacy_1st_vs_single_group.jpg'))
close

figure('Position', [10 50 250 250])
scatter(efficacy_single(idx), efficacy2(idx))
hold on
scatter(efficacy_single(~idx), efficacy2(~idx))
plot([0 20], [0, 20], 'k')
xlim([0 20])
ylim([0 20])
xlabel('efficacy of single spikes')
ylabel('efficacy of 2nd spikes')
saveas(gcf, fullfile(fig_folder, 'efficacy_2nd_vs_single_group.jpg'))
close

figure('Position', [10 50 250 250])
scatter(efficacy1(idx), efficacy2(idx))
hold on
scatter(efficacy1(~idx), efficacy2(~idx))
plot([0 20], [0, 20], 'k')
xlim([0 20])
ylim([0 20])
xlabel('efficacy of 1st spikes')
ylabel('efficacy of 2nd spikes')
saveas(gcf, fullfile(fig_folder, 'efficacy_2nd_vs_1st_group.jpg'))
close 
%% helper functions
function [r, nspk, sspktimes] = get_single_spk(r, spiketimes, DT)
isi = diff(spiketimes);
single_spk = find(isi(1:end-1) > DT & isi(2:end) > DT)+1;
sspktimes = spiketimes(single_spk);
nspk = length(single_spk);
r = r(ismember(r(:,2), single_spk), :);
idx_map = zeros(1, max(single_spk));
single_spk = [single_spk, (1:length(single_spk))'];
idx_map(single_spk(:,1)) = single_spk(:,2);
r(:,2) = idx_map(r(:,2));
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

function raster_format2(r, nspk, fig)
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
yt = [0 yt];
for kk = 1:length(yt)-1
    plot(isi(kk)*[1 1], [yt(kk), yt(kk+1)], 'b-')
end
yyaxis right
set(fig, 'YColor', [0 0 0])
ylim([1 nspk])
yt(1) = [];
yticks(yt)
yticklabels(isi)
ylabel('ISI (ms)')
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