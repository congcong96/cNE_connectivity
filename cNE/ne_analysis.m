addpath(genpath('/home/conghu/MatlabCodes/cNE_MGB_A1'));
addpath(genpath('/home/conghu/MatlabCodes/SqMo_cNE'));
addpath(genpath('/home/conghu/MatlabCodes/support')) 

spkfolder = '/data/congcong/rat_MGB_A1/3_singleunit/dmr';
save_path = '/data/congcong/rat_MGB_A1/4_cNE';
trig_path = '/data/congcong/rat_MGB_A1/trigger';
stimfolder = '/data/congcong/stimulus/thalamus';
stimfile = fullfile(stimfolder, 'rn1-500flo-40000fhi-0-4SM-0-40TM-40db-96khz-48DF-15min-seed190506_DFt1_DFf5_stim.mat');
cd(spkfolder)

%% PART1: get cNEs
files = dir('*H31x64*-split.mat');
files = {files.name};
binsize = 10; %ms
nefolder = '/data/congcong/rat_MGB_A1/4_cNE/MGB';
ne_rn_data_processing2(files, 2*binsize, 1, nefolder, {stimfile})

%% PART2-a: get NE activity threshold and NE spikes
cd(nefolder)
df = 20;
nefiles = dir('*H31x64*20dft.mat');
for ii = 1:length(nefiles)
    spkfile = dir(fullfile(spkfolder, [nefiles(ii).name(1:end-13) '*split.mat']));
    load(fullfile(spkfile.folder, spkfile.name), 'spktrain_spon')
    
    load(nefiles(ii).name, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    Activities = assembly_activity_v2(nedata.Patterns, nedata.spktrain, nedata.NEmembers_2018);
    exp_site_nedata.nedata.Activities = Activities;
    
    extra_bins = mod(length(spktrain_spon), df);
    spktrain_spon = spktrain_spon(:, 1:end-extra_bins);
    spktrain_spon(spktrain_spon > 0) = 1;
    spktrain_spon = reshape(spktrain_spon, [size(spktrain_spon, 1), df, size(spktrain_spon, 2)/df]);
    spktrain_spon = squeeze(sum(spktrain_spon, 2));
    Activities = assembly_activity_v2(nedata.Patterns, spktrain_spon, nedata.NEmembers_2018);
    exp_site_nedata.nedata.Activities_spon = Activities;
    
    [thresh, alpha] = ne_calc_NE_act_thresholds(exp_site_nedata,'circular', 50, 99:0.1:99.9);
    exp_site_nedata.nedata.NEthresh_2018 = thresh;
    exp_site_nedata.nedata.NEthresh_alpha_2018 = alpha;
    
    save(nefiles(ii).name, 'exp_site_nedata')
end
%%
NEtrain_method = 'last'; % sliding
if strcmp(NEtrain_method, 'last')
    % 0.5 ms binned NEtrain, where the last spike is
    % also checks the complexity of NE events
    ne_get_NEtrain_ne_spikes(nefiles, spkfolder)
    
elseif strcmp(NEtrain_method, 'sliding')
    for ii = 1:length(nefiles)
        % get 0.5 ms binned NEtrain (sliding window) and ne spikes
        % ne event times are when there are multiple member firing in the
        % sliding window
        % mark all 0.5ms bins within range of multiple NE member spikes as
        % event
        nefile = fullfile(nefiles(ii).folder, nefiles(ii).name);
        load(nefile, 'exp_site_nedata');
        nedata = exp_site_nedata.nedata;
        spkfile = dir(fullfile(spkfolder, [nefiles(ii).name(1:end-13) '*split.mat']));
        load(fullfile(spkfile.folder, spkfile.name), 'spktrain')
        nedata = ne_get_NEtrain_sliding_window(nedata, spktrain, 'NEmembers_2018');
        nedata = ne_get_ne_spikes(nedata, spktrain, 'NEmembers_2018');
        exp_site_nedata.nedata = nedata;
        save(nefile, 'exp_site_nedata', '-append')
    end
end
%% PART5: get strf crh and rtf for ne and neurons
ne_calc_ne_neuron_strf_rtf_crh(nefiles, spkfolder, stimfile)