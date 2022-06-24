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
pval = 0.02;
pspacex = 0.074;
pwidth = 0.052;
colors = cbrewer('qual', 'Set1', 9);
for ii = 1:length(nefiles)
    load(fullfile(nefiles(ii).folder, nefiles(ii).name), 'exp_site_nedata');
    fprintf('(%d/%d) Plotting ne member properties for %s\n', ii, length(nefiles), nefiles(ii).name)
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
    for jj = 1:size(nedata.Patterns,2)
        members = nedata.NEmembers_2018{jj};
        nmembers = length(members);
        if nmembers < 2
            continue
        end
        figure('visible', 'off'); figuresetup2seeonthescreen;
        % plot Pattern
        subplot('position', [0.04 3/8 1/20 4/7]) %'H31x64'
        hold on
        stem(idx_depth, nedata.Patterns(:,jj), 'filled', 'color', 0.5*[1 1 1])
        stem(idx_depth(members), nedata.Patterns(members,jj), 'filled', 'color', 'r')
        plot([0 nneuron+1], nedata.CI(1)*[1 1], 'b')
        plot([0 nneuron+1], nedata.CI(2)*[1 1], 'b')
        xlim([0 nneuron+1])
        view([90 90])
        ylabel('ICweight')
        xlabel('Neuron #')
        box off
        for kk = 1:nmembers
            text(idx_depth(members(kk)), nedata.Patterns(members(kk), jj)+.1, sprintf('u%03d', spk.spk(members(kk)).unit))
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
        if nmembers>4
            nrows = nmembers + 1;
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
        title('non-member pairs')
        ylim([min([y1, y2]), max([y1, y2])])
        y = ylim;
        subplot('position', [0.04 .04 .06 .1])
        ylim(y);
        % raster plot for cNE and members
        subplot('position', [.22, .85, 0.75, .1])
        hold on
        % find the 5s window with most NE activity
        netrain = nedata.NEtrain(jj,:);
        extrabins = mod(length(netrain), 1e4);
        netrain = reshape(netrain(1:end-extrabins), [1e4, floor(length(netrain)/1e4)]);
        netrain = sum(netrain);
        [~, idx_5s] = max(netrain);
        idx_5ms = ((idx_5s-1)*1e4 + 1): idx_5s*1e4;
        idx_spk = find(nedata.NEtrain(jj,idx_5ms) > 0);
        x = (idx_5ms(1) + idx_spk-1)*0.5/1e3; % cNE event times (s)
        plot([x; x], [zeros(1, length(x)); 0.9*ones(1, length(x))], 'r')
        for kk = 1:nmembers
            idx_spk = find(spktrain(members(kk),idx_5ms) > 0);
            x = (idx_5ms(1) + idx_spk-1)*0.5/1e3; % spk times (s)
            plot([x; x], [kk*ones(1, length(x));(kk+0.9) * ones(1, length(x))], 'k')
        end
        set(gca, 'YDir', 'reverse')
        box off
        xticks([])
        ylim([0 nmembers+1])
        yticks(0.5:(nmembers+1))
            % set ylabels
        yl = cell(1, nmembers+1);
        for kk = 2:nmembers+1
            yl{kk} = sprintf('u%03d', spk.spk(members(kk-1)).unit);
        end
        yl{1} = sprintf('cNE%d', jj);
        yticklabels(yl)
        set(gca, 'TickLength', [2e-3 2e-3])
        xlim([idx_5ms(1), idx_5ms(end)]*0.5/1e3);
        
        % plot for NE activity
        subplot('position', [.22, .75, 0.75, .08])
        hold on        
        idx_df = unique(ceil(idx_5ms/df));
        taxis = idx_df*df/2e3;
        plot(taxis, nedata.Activities(jj,idx_df), 'k', 'LineWidth', 1.8)
        thresh = nedata.NEthresh_2018(1,jj);
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
        
        % cNE information
        subplot('position', [.2, .7-pspacey, pwidth, pheight]) 
        dur = spk.spk(1).dur_sec_dmr;
        n0 = sum(nedata.NEtrain(jj,:));
        fr = n0/dur;
        addTextToAxis(gca, 10,...
            sprintf('unit ID: cNE%d', jj),...
            'type: ne',...
            sprintf('#ne spk: %d', sum(nedata.NEtrain(jj,:))),...
            sprintf('FR(Hz): %.2f', fr))
        axis off
        
        subplot('position', [.2+2*pspacex, .7-pspacey, pwidth, pheight])
        title('STRF(all)')
        axis off
        subplot('position', [.2+3*pspacex, .7-pspacey, pwidth, pheight])
        title('RTF(all)')
        axis off
        subplot('position', [.2+4*pspacex, .7-pspacey, pwidth, pheight])
        title('CRH(all)')
        axis off
        
        
        % plot strf
        subplot('position', [.2+5*pspacex, .7-pspacey, pwidth, pheight])
        rf = nedata.NEstrf.ne{jj};
        strf_taxis = nedata.NEstrf.taxis;
        strf_faxis = nedata.NEstrf.faxis;
        plot_strf(rf, n0, pval, 40, dur, ...
            strf_taxis, strf_faxis, ...
            'freqlabels', 2.^(-1:5), ...
            'contur', 'off');
        xlabel('')
        title('STRF(ne)')
        
        % plot rtf
        subplot('position', [.2+6*pspacex, .7-pspacey, pwidth, pheight])
        plot_rtf(nedata.NErtf.ne{jj})
        xlabel('')
        title('RTF(ne)')
        
        % plot crh
        subplot('position', [.2+7*pspacex, .7-pspacey, pwidth, pheight])
        plot_CRH(nedata.NEcrh.ne(:,jj), nedata.NEcrh.taxis, nedata.NEcrh.faxis, 'spectral')
        xlabel('')
        title('CRH(ne)')
        
        subplot('position', [.2+8*pspacex, .7-pspacey, pwidth, pheight])
        title('STRF(non-ne)')
        axis off
        subplot('position', [.2+9*pspacex, .7-pspacey, pwidth, pheight])
        title('RTF(non-ne)')
        axis off
        subplot('position', [.2+10*pspacex, .7-pspacey, pwidth, pheight])
        title('CRH(non-ne)')
        axis off
        
        for kk = 1:length(members)
            % plot properties of the A1 neuron
            subplot('position', [.2, .7-(kk+1)*pspacey, pwidth, pheight])
            depth_frac = (spk.spk(members(kk)).position(2)-min(probe_y))/(max(probe_y)-min(probe_y));
            tpd = wavestruct(members(kk)).tpd;
            if tpd <= 0.35
                type = 'NS';
            else
                type = 'BS';
            end
            addTextToAxis(gca, 10, ...
                sprintf('unit ID: unit%03d', spk.spk(members(kk)).unit),...
                sprintf('type: %s', type),...
                sprintf('depth(um): %d', spk.spk(members(kk)).position(2)),...
                sprintf('#ne spk: %d', sum(nedata.ne_spikes{jj}(kk,:))),...
                sprintf('FR(Hz): %.2f', spk.spk(members(kk)).fr_dmr))
            axis off
            
            % plot waveform
            subplot('position', [.2+1*pspacex, .7-(kk+1)*pspacey, pwidth, pheight])
            plot_waveform(wavestruct(members(kk)), 0.35)
            xlim([-1.5 1.5])
            axis off
            % plot strf
            subplot('position', [.2+2*pspacex, .7-(kk+1)*pspacey, pwidth, pheight])
            rf = strf(members(kk));
            plot_strf(rf.rfcontra, 0, pval, rf.mdb, dur, rf.taxis, rf.faxis,...
                'timelabels', 0:25:75, ...
                'freqlabels', 2.^(-1:5), 'contur', 'off');
            if kk < length(members)
                xlabel('')
            end
            % plot rtf
            subplot('position', [.2+3*pspacex, .7-(kk+1)*pspacey, pwidth, pheight])
            plot_rtf(rtf(members(kk)))
            if kk < length(members)
                xlabel('')
            end
            % plot crh
            subplot('position', [.2+4*pspacex, .7-(kk+1)*pspacey, pwidth, pheight])
            plot_CRH(crh(members(kk)).mtfhist, crh(members(kk)).tmfaxis, crh(members(kk)).smfaxis, 'spectral')
            if kk < length(members)
                xlabel('')
            end
            
            % plot ne strf
            subplot('position', [.2+5*pspacex, .7-(kk+1)*pspacey, pwidth, pheight])
            plot_strf(nedata.NEstrf.member{jj}(:,kk), 0, pval, rf.mdb, dur, strf_taxis, strf_faxis, 'freqlabels', 2.^(-1:5), 'contur', 'off');
            if kk < length(members)
                xlabel('')
            end
            
            % plot ne rtf
            subplot('position', [.2+6*pspacex, .7-(kk+1)*pspacey, pwidth, pheight])
            plot_rtf(nedata.NErtf.member{jj}(kk))
            if kk < length(members)
                xlabel('')
            end
            
            % plot ne crh
            subplot('position', [.2+7*pspacex, .7-(kk+1)*pspacey, pwidth, pheight])
            plot_CRH(nedata.NEcrh.member{jj}(:,kk), nedata.NEcrh.taxis, nedata.NEcrh.faxis, 'spectral')
            if kk < length(members)
                xlabel('')
            end
            
            % plot nonne strf
            subplot('position', [.2+8*pspacex, .7-(kk+1)*pspacey, pwidth, pheight])
            plot_strf(nedata.NEstrf.member_nonne{jj}(:,kk), 0, pval, rf.mdb, dur, strf_taxis, strf_faxis, 'freqlabels', 2.^(-1:5), 'contur', 'off');
            if kk < length(members)
                xlabel('')
            end
            % plot nonne rtf
            subplot('position', [.2+9*pspacex, .7-(kk+1)*pspacey, pwidth, pheight])
            plot_rtf(nedata.NErtf.member_nonne{jj}(kk))
            if kk < length(members)
                xlabel('')
            end
            
            % plot nonne crh
            subplot('position', [.2+10*pspacex, .7-(kk+1)*pspacey, pwidth, pheight])
            plot_CRH(nedata.NEcrh.member_nonne{jj}(:,kk), nedata.NEcrh.taxis, nedata.NEcrh.faxis, 'spectral')
            if kk < length(members)
                xlabel('')
            end
            
        end
        
        saveas(gcf, fullfile(figfolder, ...
            sprintf('unit_properties-%s-cNE%d.jpg',spk.exp, jj)));
        close
    end
end