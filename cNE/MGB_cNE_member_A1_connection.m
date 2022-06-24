addpath(genpath('/home/conghu/MatlabCodes/cNE_connectivity'))

nefolder = '/data/congcong/rat_MGB_A1/4_cNE/MGB';
spkfolder = '/data/congcong/rat_MGB_A1/3_singleunit/dmr';
cd(spkfolder)
%% get efficacy and causal spikes
connect_files = dir(fullfile(spkfolder, '*pairs.mat'));
ccg_taxis = -50:0.5:50;
causal_window = 1:0.5:4.5;
for ii = 1:length(connect_files)
    data = load(fullfile(connect_files(ii).folder, connect_files(ii).name), 'connected_pairs');
    cp = data.connected_pairs;
    if isempty(cp)
        fprintf('(%d/%d) No connected pairs in %s...\n',  ii, length(connect_files), connect_files(ii).name)
        continue
    end
    
    fprintf('(%d/%d) Processing %s...\n', ii, length(connect_files), connect_files(ii).name)
    c = 1;
    % load ne file
    nefile = dir(fullfile(nefolder, [connect_files(ii).name(1:13) '*H31x64*20dft.mat']));
    load(fullfile(nefolder, nefile.name), 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    nNE = size(nedata.Patterns,2);
    % load MGB neuron spk data
    spkfile_MGB = dir(fullfile(spkfolder, [connect_files(ii).name(1:13) '*H31x64*split.mat']));
    data = load(fullfile(spkfile_MGB.folder, spkfile_MGB.name), 'spk', 'spktrain', 'spktrain_spon', 'trigger');
    spk_MGB = data.spk;
    spktrain_dmr_MGB = data.spktrain;
    spktrain_spon_MGB = data.spktrain_spon;
    trigger = data.trigger;
    % load A1 neuron spk data
    spkfile_A1 = dir(fullfile(spkfolder, [connect_files(ii).name(1:13) '*H22x32*split.mat']));
    data = load(fullfile(spkfile_A1.folder, spkfile_A1.name), 'spk', 'spktrain', 'spktrain_spon');
    spk_A1 = data.spk;
    spktrain_dmr_A1 = data.spktrain;
    spktrain_spon_A1 = data.spktrain_spon;
    
    for jj = 1:length(cp)
        idx_MGB = cp(jj).idx_MGB;
        idx_A1= cp(jj).idx_A1;
        for kk = 1:nNE
            members = nedata.NEmembers_2018{kk};
            nmembers = length(members);
            if ~sum(members == idx_MGB)
                continue
            end
            cNE_connected(c).idx_MGB = cp(jj).idx_MGB;
            cNE_connected(c).idx_A1 = cp(jj).idx_A1;
            cNE_connected(c).unit_MGB = cp(jj).unit_MGB;
            cNE_connected(c).unit_A1 = cp(jj).unit_A1;
            cNE_connected(c).cNE = kk;
            cNE_connected(c).members = members;
            %% get ccg and efficacy under all conditions
            conditions = {'spon', 'dmr'};
            for cc = 1:length(conditions)
                condition = conditions{cc};
                % -------------neuron------------------------
                ccg = nan(201, nmembers);
                baseline = nan(201, nmembers);
                efficacy = nan(1, nmembers);
                n0 = nan(1, nmembers);
                % get cortical spiketimes
                spiketimes = eval(sprintf('spktrain_%s_A1(idx_A1,:)', condition));
                spiketimes = find(spiketimes > 0) * .5;
                for mm = 1:nmembers
                    reftimes = eval(sprintf('spktrain_%s_MGB(members(mm),:)', condition));
                    reftimes = find(reftimes > 0) * .5;
                    nspk_MGB = length(reftimes);
                    n0(mm) = nspk_MGB;
                    [baseline(:,mm), ccg(:, mm)] = get_baseline(spiketimes, reftimes, 0.5);
                    efficacy(mm) = get_efficacy(ccg(:,mm), baseline(:,mm), ccg_taxis, causal_window, nspk_MGB);
                end
                
                cNE_connected(c).ccg_taxis = ccg_taxis;
                cNE_connected(c).(sprintf('n0_%s', condition)) = n0;
                cNE_connected(c).(sprintf('baseline_%s', condition)) = baseline;
                cNE_connected(c).(sprintf('ccg_%s', condition)) = ccg;
                cNE_connected(c).(sprintf('ccg_efficacy_%s', condition)) = efficacy;
                
                % NE activity
                if strcmp(condition, 'spon')
                    reftimes = 0.5*(find(nedata.NEtrain_spon(kk,:)>0));
                else
                    reftimes = 0.5*(find(nedata.NEtrain(kk,:)>0));  
                end
                [baseline, ccg] = get_baseline(spiketimes, reftimes, 0.5);
                efficacy = get_efficacy(ccg, baseline, ccg_taxis, causal_window, length(reftimes));
                cNE_connected(c).(sprintf('baseline_cNE_%s', condition)) = baseline;
                cNE_connected(c).(sprintf('ccg_cNE_%s', condition)) = ccg;
                cNE_connected(c).(sprintf('ccg_efficacy_cNE_%s', condition)) = efficacy;
                
                % ne spikes
                ccg = nan(201, nmembers);
                baseline = nan(201, nmembers);
                efficacy = nan(1, nmembers);
                n0 = nan(1, nmembers);
                for mm = 1:nmembers
                    if strcmp(condition, 'dmr')
                        reftimes = nedata.ne_spikes{kk}(mm,:);
                    elseif strcmp(condition, 'spon')
                        reftimes = nedata.ne_spikes_spon{kk}(mm,:);
                    end
                    reftimes = find(reftimes > 0) * .5;
                    n0(mm) = length(reftimes);
                    [baseline(:,mm), ccg(:, mm)] = get_baseline(spiketimes, reftimes, 0.5);
                    efficacy(mm) = get_efficacy(ccg(:, mm), baseline(:, mm), ccg_taxis, causal_window, n0(mm));
                end
                cNE_connected(c).(sprintf('n0_ne_%s', condition)) = n0;
                cNE_connected(c).(sprintf('baseline_ne_%s', condition)) = baseline;
                cNE_connected(c).(sprintf('ccg_ne_%s', condition)) = ccg;
                cNE_connected(c).(sprintf('ccg_efficacy_ne_%s', condition)) = efficacy;
                
                % ccg for nonne spikes
                ccg = nan(201, nmembers);
                baseline = nan(201, nmembers);
                efficacy = nan(1, nmembers);
                n0 = nan(1, nmembers);
                spktrain = eval(sprintf('spktrain_%s_MGB(members,:)', condition));
                if strcmp(condition, 'spon')
                    ne_spikes = nedata.ne_spikes_spon{kk};
                elseif strcmp(condition, 'dmr')
                    ne_spikes = nedata.ne_spikes{kk};
                end
                extra = length(spktrain) - length(ne_spikes);
                spktrain = spktrain(:, 1:end-extra);
                nonne_spikes = spktrain - ne_spikes;
                for mm = 1:nmembers
                    reftimes = find(nonne_spikes(mm,:) > 0) * .5;
                    n0(mm) = length(reftimes);
                    [baseline(:,mm), ccg(:, mm)] = get_baseline(spiketimes, reftimes, 0.5);
                    efficacy(mm) = get_efficacy(ccg(:, mm), baseline(:, mm), ccg_taxis, causal_window, n0(mm));
                end
                cNE_connected(c).(sprintf('n0_nonne_%s', condition)) = n0;
                cNE_connected(c).(sprintf('baseline_nonne_%s', condition)) = baseline;
                cNE_connected(c).(sprintf('ccg_nonne_%s', condition)) = ccg;
                cNE_connected(c).(sprintf('ccg_efficacy_nonne_%s', condition)) = efficacy;
                
            end
            c = c+1;
        end
    end
    if exist('cNE_connected', 'var')
        save(fullfile(connect_files(ii).folder, connect_files(ii).name), 'cNE_connected', '-append')
        clear cNE_connected
    end
end