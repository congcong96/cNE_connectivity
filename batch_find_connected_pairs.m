addpath(genpath('/home/conghu/MatlabCodes/cNE_connectivity'))
data_path = '/data/congcong/rat_MGB_A1/3_singleunit/dmr/ks25_ks2';
crtxfiles = dir(fullfile(data_path, '*H22x32*split.mat'));
crtx_exp = cellfun(@(x) x(1:13), {crtxfiles.name}, 'UniformOutput', false);
binsize = 0.5;
bandwidth = [75, Inf];
hw = 100;
nfiles = length(crtxfiles);

%% find connected pairs
parfor ii = 1:length(crtxfiles)
    
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