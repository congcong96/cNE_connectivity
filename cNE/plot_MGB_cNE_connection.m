% run MGB_cNE_member_A1_connection first
addpath(genpath('/home/conghu/MatlabCodes/cNE_MGB/connection'))
addpath(genpath('/home/conghu/MatlabCodes/MultiUnit_IntanAnalysis'))
addpath(genpath('/home/conghu/MatlabCodes/SingleUnit'))
addpath(genpath('/home/conghu/MatlabCodes/support'))
addpath(genpath('/home/conghu/MatlabCodes/ks2_to_clust3r'))
spkfolder = '/data/congcong/rat_MGB_A1/3_singleunit/dmr';
nefolder = '/data/congcong/rat_MGB_A1/4_cNE/MGB';
figfolder = '/data/congcong/rat_MGB_A1/figure/cNE/MGB';
cd(nefolder)
nefiles = dir(fullfile(nefolder, '*H31x64*20dft.mat'));
%%
connect_files = dir(fullfile(spkfolder, '*pairs.mat'));
pval = 0.02;
pspacex = 0.074;
pwidth = 0.052;
colors = cbrewer('qual', 'Set1', 9);
for ii = 1:length(connect_files)
    connected_file = fullfile(connect_files(ii).folder, connect_files(ii).name);
    vars = whos('-file', connected_file);
    if ~ismember('cNE_connected', {vars.name})
        fprintf('(%d/%d) No cNE with connected cortical neuron in %s\n',...
            ii, length(connect_files), connect_files(ii).name)
        continue
    end
    fprintf('(%d/%d) Plotting cNE connection to A1 %s\n',...
        ii, length(connect_files), connect_files(ii).name)
    % load connected pair information
    load(connected_file, 'cNE_connected');
    % load nefile
    nefile = dir(fullfile(nefolder, [nefiles(ii).name(1:13) '*H31x64*20dft.mat']));
    load(fullfile(nefile.folder, nefile.name), 'exp_site_nedata');
    nedata = exp_site_nedata.nedata;
    df = nedata.df;
    % load MGB neuron spk data
    spkfile = dir(fullfile(spkfolder, [nefiles(ii).name(1:13) '*H31x64*split.mat']));
    data = load(fullfile(spkfile.folder, spkfile.name), 'spk', 'waveform', 'strf', 'rtfparams', 'crh', 'spktrain');
    spk = data.spk;
    wavestruct = data.waveform;
    strf = data.strf;
    rtf = data.rtfparams;
    crh = data.crh;
    spktrain = data.spktrain;
    % load A1 neuron spk data
    spkfile_A1 = dir(fullfile(spkfolder, [nefiles(ii).name(1:13) '*H22x32*split.mat']));
    data = load(fullfile(spkfile_A1.folder, spkfile_A1.name), 'spk', 'waveform', 'strf', 'rtfparams', 'crh');
    spk_A1 = data.spk;
    wavestruct_A1 = data.waveform;
    strf_A1 = data.strf;
    rtf_A1 = data.rtfparams;
    crh_A1 = data.crh;
    
    % downsample spktrain to 1ms for xcorr
    nbins = 2;
    nneuron = size(spktrain,1);
    extra = mod(size(spktrain,2), nbins);
    spktrain_tmp = spktrain(:,1:end-extra);
    spktrain_tmp = reshape(spktrain_tmp, [nneuron*nbins, size(spktrain_tmp,2)/nbins]);
    spktrain_downsampled = zeros(nneuron,  size(spktrain_tmp,2));
    for kk = 1:nbins
        spktrain_downsampled = spktrain_downsampled + spktrain_tmp(nneuron*(kk-1)+1:kk*nneuron,:);
    end
    spktrain_downsampled = spktrain_downsampled';
    
    % get penetration probe information
    [probinfo] = neuronexus_prob({'H31x64'});
    probe_x = probinfo.posi_x;
    probe_y = spk.depth - probinfo.posi_depth;
    position = reshape(cell2mat(nedata.position), [2, nneuron])';
    [~, idx_depth] = sortrows(position, 'ascend');
    [~, idx_depth] = sort(idx_depth);
    for jj = 1:length(cNE_connected)
        members = cNE_connected(jj).members;
        cNE = cNE_connected(jj).cNE;
        idx_A1 = cNE_connected(jj).idx_A1;
        nmembers = length(members);
        if nmembers < 2
            continue
        end
        figure('visible', 'off'); figuresetup2seeonthescreen;
        % plot Patterns
        subplot('position', [0.04 3/8 1/20 4/7]) %'H31x64'
        hold on
        stem(idx_depth, nedata.Patterns(:,cNE), 'filled', 'color', 0.5*[1 1 1])
        stem(idx_depth(members), nedata.Patterns(members,cNE), 'filled', 'color', 'r')
        plot([0 nneuron+1], nedata.CI(1)*[1 1], 'b')
        plot([0 nneuron+1], nedata.CI(2)*[1 1], 'b')
        xlim([0 nneuron+1])
        view([90 90])
        ylabel('ICweight')
        xlabel('Neuron #')
        box off
        for kk = 1:nmembers
            text(idx_depth(members(kk)), nedata.Patterns(members(kk), cNE)+.1, sprintf('u%03d', spk.spk(members(kk)).unit))
        end
        % get probe and location infromation for member neurons
        subplot('position', [0.13 3/8 1/20 4/7]) %'H31x64'
        scatter(probe_x, probe_y, '.','MarkerEdgeColor', 0.5*[1 1 1])
        ax = gca;
        ax.XAxis.Visible = 'off';
        ax.YAxis.Direction = 'reverse';
        title('depth (um)')
        hold on
        for kk = 1:nmembers
            position = spk.spk(members(kk)).position;
            scatter(position(1), position(2), 'ro', 'filled')
            text(position(1), position(2), sprintf(' u%0d',spk.spk(members(kk)).unit), 'fontsize', 8)
        end
        if nmembers>3
            nrows = nmembers + 2;
        else
            nrows = 5;
        end
        pheight = 1/nrows*0.5;
        pspacey = 1/nrows*0.65;
        
        % plot correlation matrix for the members neurons
        subplot('position', [0.04 .2 .1 .12])
        corrmat = corr(nedata.spktrain(members,:)');
        corrmat = corrmat - eye(size(corrmat));
        cmap = cschemes('rdbu', 9);
        imagesc(corrmat),colormap(gca, cmap), colorbar, caxis([-0.4, 0.4])
        unit = [spk.spk(members).unit];
        xticks(1:nmembers);
        xticklabels(unit)
        xtickangle(45)
        yticks(1:nmembers)
        yticklabels(unit)
        
        % get cross correlation of members and memberv.s. non-members
        member_pair = nchoosek(members, 2);
        nonmembers = setdiff(1:nneuron, members);
        [A, B] = meshgrid(members, nonmembers);
        member_nonmember = [A(:), B(:)];
        % get corss correlation of neuronal pairs - member
        nbins = 50;
        xcorr_member = zeros(nbins*2+1, size(member_pair, 1));
        for kk = 1:size(member_pair, 1)
            xcorr_member(:,kk) = zscore(...
                xcorr(...
                spktrain_downsampled(:, member_pair(kk,1)), ...
                spktrain_downsampled(:, member_pair(kk,2)), ...
                nbins));
        end
        % get corss correlation of neuronal pairs : member-nonmember
        xcorr_nonmember = zeros(nbins*2+1, size(member_nonmember, 1));
        for kk = 1:size(member_nonmember, 1)
            xcorr_nonmember(:,kk) = zscore(...
                xcorr(...
                spktrain_downsampled(:, member_nonmember(kk,1)), ...
                spktrain_downsampled(:, member_nonmember(kk,2)), ...
                nbins));
        end
        subplot('position', [0.04 .04 .06 .1])
        hold on
        curve1 = (mean(xcorr_member,2) + std(xcorr_member,[], 2))';
        curve2 = (mean(xcorr_member,2) - std(xcorr_member,[], 2))';
        taxis = -50:50;
        a = fill([taxis, fliplr(taxis)], [curve1, fliplr(curve2)], 'k');
        a.FaceColor = 0.5*[1 1 1];
        a.EdgeColor = 0.5*[1 1 1];
        plot(taxis, mean(xcorr_member,2), 'k', 'linewidth', 2)
        y1 = ylim;
        xlabel('time lag (ms)')
        ylabel('c.c.(z-scored)')
        title('member pairs')
        
        subplot('position', [0.12 .04 .06 .1])
        hold on
        curve1 = (mean(xcorr_nonmember,2) + std(xcorr_nonmember,[], 2))';
        curve2 = (mean(xcorr_nonmember,2) - std(xcorr_nonmember,[], 2))';
        taxis = -50:50;
        a = fill([taxis, fliplr(taxis)], [curve1, fliplr(curve2)], 'k');
        a.FaceColor = 0.5*[1 1 1];
        a.EdgeColor = 0.5*[1 1 1];
        plot(taxis, mean(xcorr_nonmember,2), 'k', 'linewidth', 2)
        y2 = ylim;
        xlabel('time lag (ms)')
        title('member pairs')
        ylim([min([y1, y2]), max([y1, y2])])
        y = ylim;
        subplot('position', [0.04 .04 .06 .1])
        ylim(y);
        % raster plot for cNE and members
        subplot('position', [.22, .85, 0.75, .1])
        hold on
            % find the 5s window with most NE activity
        netrain = nedata.NEtrain(cNE,:);
        extrabins = mod(length(netrain), 1e4);
        netrain = reshape(netrain(1:end-extrabins), [1e4, floor(length(netrain)/1e4)]);
        netrain = sum(netrain);
        [~, idx_5s] = max(netrain);
        idx_5ms = ((idx_5s-1)*1e4 + 1): idx_5s*1e4;
        idx_spk = find(nedata.NEtrain(cNE,idx_5ms) > 0);
        x = (idx_5ms(1) + idx_spk)*0.5/1e3; % cNE event times (s)
        plot([x; x], [zeros(1, length(x)); 0.9*ones(1, length(x))], 'r')
        for kk = 1:nmembers
            idx_spk = find(spktrain(members(kk),idx_5ms) > 0);
            x = (idx_5ms(1) + idx_spk)*0.5/1e3; % spk times (s)
            plot([x; x], [kk*ones(1, length(x));(kk+0.9) * ones(1, length(x))], 'k')
        end
        set(gca, 'YDir', 'reverse')
        box off
        xticks([])
        ylim([0 nmembers+1])
        yticks(0.5:(nmembers+1))
        xlim([idx_5ms(1) idx_5ms(end)]*0.5/1e3)
            % set ylabels
        yl = cell(1, nmembers+1);
        for kk = 2:nmembers+1
            yl{kk} = sprintf('u%03d', spk.spk(members(kk-1)).unit);
        end
        yl{1} = sprintf('cNE%d', cNE);
        yticklabels(yl)
        set(gca, 'TickLength', [2e-3 2e-3])
        
        % plot for NE activity
        subplot('position', [.22, .75, 0.75, .08])
        hold on        
        idx_df = unique(ceil(idx_5ms/df));
        taxis = idx_df*df/2/1e3;
        plot(taxis, nedata.Activities(cNE,idx_df), 'k', 'LineWidth', 1.8)
        thresh = nedata.NEthresh_2018(1,cNE);
        if thresh > 0
            xl = xlim;
            plot(xl, thresh*[1 1], 'color', colors(5,:))
        end
        yl = ylim;
        text(taxis(end)-0.5, yl(1)+(yl(2)-yl(1))*0.8, sprintf('cNE thresh: %.2f', thresh), 'color', 'r')
        text(taxis(end)-0.5, yl(1)+(yl(2)-yl(1))*0.6, '(99 percentile)', 'color', 'r')
        xlabel('Time From Stimulus Onset (s)')
        ylabel('Activity (a.u.)')
        box off
        %% cortical neuron infromation
        subplot('position', [.2, .7-pspacey, pwidth, pheight]) 
        addTextToAxis(gca, 10,...
            sprintf('unit ID: A1-u%d', spk_A1.spk(idx_A1).unit),...
            sprintf('FR(Hz): %.2f', spk_A1.spk(idx_A1).fr_dmr))
        axis off
        
        % plot waveform
        subplot('position', [.2+pspacex, .7-pspacey, pwidth, pheight])
        plot_waveform(wavestruct_A1(idx_A1), 0.5)
        axis off
        
        % plot strf
        subplot('position', [.2+2*pspacex, .7-pspacey, pwidth, pheight])
        dur = spk.spk(1).dur_sec_dmr;
        rf = strf_A1(idx_A1);
        plot_strf(rf.rfcontra, 0, pval, rf.mdb, dur, rf.taxis, rf.faxis,...
            'timelabels', 0:25:75, ...
            'freqlabels', 2.^(-1:5), 'contur', 'off');
        xlabel('')
        ylabel('')
        
        % plot rtf
        subplot('position', [.2+3*pspacex, .7-pspacey, pwidth, pheight])
        plot_rtf(rtf_A1(idx_A1))
        xlabel('')
        title('RTF')
        
        % plot crh
        subplot('position', [.2+4*pspacex, .7-pspacey, pwidth, pheight])
        plot_CRH(crh_A1(idx_A1).mtfhist, crh_A1(idx_A1).tmfaxis, crh_A1(idx_A1).smfaxis, 'spectral')
        xlabel('')
        title('CRH')
        
        subplot('position', [.2+5*pspacex, .7-pspacey, pwidth, pheight])
        title('CCG spon')
        axis off
        subplot('position', [.2+6*pspacex, .7-pspacey, pwidth, pheight])
        title('ne CCG spon')
        axis off
        subplot('position', [.2+7*pspacex, .7-pspacey, pwidth, pheight])
        title('nonne CCG spon')
        axis off
        subplot('position', [.2+8*pspacex, .7-pspacey, pwidth, pheight])
        title('CCG dmr')
        axis off
        subplot('position', [.2+9*pspacex, .7-pspacey, pwidth, pheight])
        title('ne CCG dmr')
        axis off
        subplot('position', [.2+10*pspacex, .7-pspacey, pwidth, pheight])
        title('nonne CCG dmr')
        axis off
        %% cNE information
        subplot('position', [.2, .7-2*pspacey, pwidth, pheight]) 
        dur = spk.spk(1).dur_sec_dmr;
        n0 = sum(nedata.NEtrain(cNE,:));
        fr = n0/dur;
        addTextToAxis(gca, 10,...
            sprintf('unit ID: cNE%d', cNE),...
            sprintf('#ne spk: %d', sum(nedata.NEtrain(cNE,:))),...
            sprintf('FR(Hz): %.2f', fr))
        axis off
        
        % plot strf
        subplot('position', [.2+2*pspacex, .7-2*pspacey, pwidth, pheight])
        rf = nedata.NEstrf.ne{cNE};
        strf_taxis = nedata.NEstrf.taxis;
        strf_faxis = nedata.NEstrf.faxis;
        plot_strf(rf, n0, pval, 40, dur, ...
            strf_taxis, strf_faxis, ...
            'freqlabels', 2.^(-1:5), ...
            'contur', 'off');
        xlabel('')
        ylabel('')
        
        % plot rtf
        subplot('position', [.2+3*pspacex, .7-2*pspacey, pwidth, pheight])
        plot_rtf(nedata.NErtf.ne{cNE})
        xlabel('')
        
        % plot crh
        subplot('position', [.2+4*pspacex, .7-2*pspacey, pwidth, pheight])
        plot_CRH(nedata.NEcrh.ne(:,cNE), nedata.NEcrh.taxis, nedata.NEcrh.faxis, 'spectral')
        xlabel('')
                
        %% members neuron properties
        % get ymax
        conditions = {'spon', 'dmr'};
        for kk = 1:length(conditions)
            eval(sprintf(...
                ['ccg = [cNE_connected(jj).ccg_%s,'...
                'cNE_connected(jj).ccg_ne_%s,'...
                'cNE_connected(jj).ccg_nonne_%s];'], ...
                conditions{kk}, conditions{kk}, conditions{kk}))
            
            eval(sprintf(...
                ['dur = [cNE_connected(jj).n0_%s,'...
                'cNE_connected(jj).n0_ne_%s,'...
                'cNE_connected(jj).n0_nonne_%s]' ...
                ' * .5 * 1e-3;'], ...
                conditions{kk}, conditions{kk}, conditions{kk}))
            
            fr = ccg./dur;
            eval(sprintf('ymax_%s = max(fr(:))*1.1;', conditions{kk}));
        end
        
        for kk = 1:length(members)
            % information
            subplot('position', [.2, .7-(kk+2)*pspacey, pwidth, pheight])
            depth_frac = (spk.spk(members(kk)).position(2)-min(probe_y))/(max(probe_y)-min(probe_y));
            tpd = wavestruct(members(kk)).tpd;
            addTextToAxis(gca, 10, ...
                sprintf('unit ID: unit%03d', spk.spk(members(kk)).unit),...
                sprintf('#ne spk: %d', sum(nedata.ne_spikes{cNE}(kk,:))),...
                sprintf('FR(Hz): %.2f', spk.spk(members(kk)).fr_dmr))
            axis off
            
            % plot waveform
            subplot('position', [.2+1*pspacex, .7-(kk+2)*pspacey, pwidth, pheight])
            plot_waveform(wavestruct(members(kk)), 0.35)
            axis off
            
            % plot strf
            subplot('position', [.2+2*pspacex, .7-(kk+2)*pspacey, pwidth, pheight])
            rf = strf(members(kk));
            plot_strf(rf.rfcontra, 0, pval, rf.mdb, dur, rf.taxis, rf.faxis,...
                'timelabels', 0:25:75, ...
                'freqlabels', 2.^(-1:5), 'contur', 'off');
            if kk < length(members)
                xlabel('')
                ylabel('')
            end
            % plot rtf
            subplot('position', [.2+3*pspacex, .7-(kk+2)*pspacey, pwidth, pheight])
            plot_rtf(rtf(members(kk)))
            if kk < length(members)
                xlabel('')
            end
            % plot crh
            subplot('position', [.2+4*pspacex, .7-(kk+2)*pspacey, pwidth, pheight])
            plot_CRH(crh(members(kk)).mtfhist, crh(members(kk)).tmfaxis, crh(members(kk)).smfaxis, 'spectral')
            if kk < length(members)
                xlabel('')
            end
            
            % plot ccg spon
            subplot('position', [.2+5*pspacex, .7-(kk+2)*pspacey, pwidth, pheight])
            dur = cNE_connected(jj).n0_spon(kk) *.5 * 1e-3;
            plot_ccg(cNE_connected(jj).ccg_taxis, ...
                cNE_connected(jj).ccg_spon(:,kk)/dur,...
                cNE_connected(jj).baseline_spon(:,kk)/dur,...
                '', '', ...
                cNE_connected(jj).ccg_efficacy_spon(kk), ...
                'r',...
                'ymax', ymax_spon)
            if kk == length(members)
                ylabel('spikes per sec')
                xlabel('time preceding MGB spike (ms)')
            end
            
            % plot ne ccg spon
            subplot('position', [.2+6*pspacex, .7-(kk+2)*pspacey, pwidth, pheight])
            dur = cNE_connected(jj).n0_ne_spon(kk) *.5 * 1e-3;
            plot_ccg(cNE_connected(jj).ccg_taxis, ...
                cNE_connected(jj).ccg_ne_spon(:,kk)/dur,...
                cNE_connected(jj).baseline_ne_spon(:,kk)/dur,...
                '', '', ...
                cNE_connected(jj).ccg_efficacy_ne_spon(kk), ...
                'r', ...
                'ymax', ymax_spon)
           
            % plot non-ne ccg
            subplot('position', [.2+7*pspacex, .7-(kk+2)*pspacey, pwidth, pheight])
            dur = cNE_connected(jj).n0_nonne_spon(kk) *.5 * 1e-3;
            plot_ccg(cNE_connected(jj).ccg_taxis, ...
                cNE_connected(jj).ccg_nonne_spon(:,kk)/dur,...
                cNE_connected(jj).baseline_nonne_spon(:,kk)/dur,...
                '', '', ...
                cNE_connected(jj).ccg_efficacy_nonne_spon(kk), ...
                'r', ...
                'ymax', ymax_spon)
            
            % plot ccg dmr
            subplot('position', [.2+8*pspacex, .7-(kk+2)*pspacey, pwidth, pheight])
            dur = cNE_connected(jj).n0_dmr(kk) *.5 * 1e-3;
            plot_ccg(cNE_connected(jj).ccg_taxis, ...
                cNE_connected(jj).ccg_dmr(:,kk)/dur,...
                cNE_connected(jj).baseline_dmr(:,kk)/dur,...
                '', '', ...
                cNE_connected(jj).ccg_efficacy_dmr(kk), ...
                'r', ...
                'ymax', ymax_dmr)
            
            % plot ne ccg dmr
            subplot('position', [.2+9*pspacex, .7-(kk+2)*pspacey, pwidth, pheight])
            dur = cNE_connected(jj).n0_ne_dmr(kk) *.5 * 1e-3;
            plot_ccg(cNE_connected(jj).ccg_taxis, ...
                cNE_connected(jj).ccg_ne_dmr(:,kk)/dur,...
                cNE_connected(jj).baseline_ne_dmr(:,kk)/dur,...
                '', '', ...
                cNE_connected(jj).ccg_efficacy_ne_dmr(kk), ...
                'r', ...
                'ymax', ymax_dmr)

            % plot nonne ccg dmr
            subplot('position', [.2+10*pspacex, .7-(kk+2)*pspacey, pwidth, pheight])
            dur = cNE_connected(jj).n0_nonne_dmr(kk) *.5 * 1e-3;
            plot_ccg(cNE_connected(jj).ccg_taxis, ...
                cNE_connected(jj).ccg_nonne_dmr(:,kk)/dur,...
                cNE_connected(jj).baseline_nonne_dmr(:,kk)/dur,...
                '', '', ...
                cNE_connected(jj).ccg_efficacy_nonne_dmr(kk), ...
                'r', ...
                'ymax', ymax_dmr)
        end
        
%         printPDFandPSC(gcf, fullfile(figfolder, ...
%             sprintf('cNE_member_connection-%s-cNE%d-MGB_unit%d-A1_unit%d',spk.exp, cNE, cNE_connected(jj).unit_MGB, cNE_connected(jj).unit_A1)));
        saveas(gcf, fullfile(figfolder, ...
            sprintf('cNE_member_connection-%s-cNE%d-MGB_unit%d-A1_unit%d.jpg',spk.exp, cNE, cNE_connected(jj).unit_MGB, cNE_connected(jj).unit_A1)));
        close
    end
end