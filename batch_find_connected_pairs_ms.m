addpath(genpath('/home/conghu/MatlabCodes/cNE_connectivity'))
data_path = '/data/congcong/rat_MGB_A1/3_singleunit/ms_dmr/6std';
%% find connected pairs
crtxfiles = dir(fullfile(data_path, '*H22x32*split.mat'));
binsize = 0.5;
bandwidth = [25, Inf];
hw = 100;
parfor ii = 1:length(crtxfiles)
    
    
    exp = crtxfiles(ii).name(1:13);
    site = regexp(crtxfiles(ii).name, '(?<=site)\d{1,2}', 'match', 'once');
    outfile = sprintf('%s-site%s-connected_pairs_test_05ms_25Hz.mat', ...
        exp,  site);
    
    if exist(fullfile(data_path,outfile), 'file')
        fprintf('(%d/%d)Already processed %s\n', ii, length(crtxfiles), crtxfiles(ii).name)
        continue
    end
    fprintf('(%d/%d)processing %s\n', ii, length(crtxfiles), crtxfiles(ii).name)
    
    
    
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
    
    % downsample spktrains according to binsize
    if binsize > 0.5
        downsample_factor = round(binsize*2);
        % dmr response
        extrabin = mod(size(locators_A1_dmr, 2), downsample_factor);
        locators_A1_dmr = locators_A1_dmr(:,1:end-extrabin);
        locators_MGB_dmr = locators_MGB_dmr(:,1:end-extrabin);
        nsamples = size(locators_A1_dmr, 2)/downsample_factor;
        locators_A1_dmr = reshape(locators_A1_dmr, ...
            [size(locators_A1_dmr, 1), downsample_factor, nsamples]);
        locators_A1_dmr = squeeze(sum(locators_A1_dmr, 2));
        locators_MGB_dmr = reshape(locators_MGB_dmr, ...
            [size(locators_MGB_dmr, 1), downsample_factor, nsamples]);
        locators_MGB_dmr = squeeze(sum(locators_MGB_dmr, 2));
        % spontaneous response
        extrabin = mod(size(locators_A1_spon, 2), downsample_factor);
        locators_A1_spon = locators_A1_spon(:,1:end-extrabin);
        locators_MGB_spon = locators_MGB_spon(:,1:end-extrabin);
        nsamples = size(locators_A1_spon, 2)/downsample_factor;
        locators_A1_spon = reshape(locators_A1_spon, ...
            [size(locators_A1_spon, 1), downsample_factor, nsamples]);
        locators_A1_spon = squeeze(sum(locators_A1_spon, 2));
        locators_MGB_spon = reshape(locators_MGB_spon, ...
            [size(locators_MGB_spon, 1), downsample_factor, nsamples]);
        locators_MGB_spon = squeeze(sum(locators_MGB_spon, 2));
    end
    
    
    
    connected_pairs = find_connected_pairs(...
        locators_MGB_spon, locators_A1_spon, ...
        locators_MGB_dmr, locators_A1_dmr, ...
        hw, binsize, [spk_MGB.spk.unit], [spk_A1.spk.unit], bandwidth, ...
        {'Normal_weak_both', 'Normal_strong', 'Normal_strong_both'...
        'Poisson_chi_square', 'Poisson_normal', ...
        'Poisson_chi_square_both', 'Poisson_normal_both'});
    parsave(fullfile(data_path,outfile), connected_pairs, 'connected_pairs')
    
end