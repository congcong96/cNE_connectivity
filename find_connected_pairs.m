function ccg_strct = find_connected_pairs(...
    spktrain_ref_spon, spktrain_target_spon, ...
    spktrain_ref_stim, spktrain_target_stim, ...
    hw, binsize, unit_ref, unit_target, bandwidth, methods)
% 
% Input:
%   spktrain_ref_spon: m x n matrix, m is the number of neurons, n is the
%   number of time bins
%       spontaneous reference spike trians binned at binsize
%   spktrain_target_spon: m x n matrix
%       spontaneous target spike trians binned at binsize
%   spktrain_ref_stim: m x n matrix
%       reference spike trians in response to stimulus binned at binsize
%   spktrain_target_stim: m x n matrix
%       target spike trians in response to stimulus binned at binsize
%   hw: integer
%       half window size, in ms
%   binsize: in ms
%   Unit_ref: vector of length m
%       unit label for thalamic neurons
%   Unit_target: vector of length m
%       unit label for cortical neurons
%   bandwith: bandwith to filter ccg
%   method: cell with string indicating method to test significance of ccg
%       Normal_weak_both
%               The significance level was set at a probability of
%               0.2%, assuming a normal distribution in the baseline
%               amplitudes after filtering. mean + 3.1std
%       Normal_strong
%               at least one condition (stim or spon) has ccg peak > 5std
%       Normal_strong_both
%               ccg peak > 5std under both conditions
%       Poisson_chi_square
%               any condition, ccg peak > 0.5*chi2inv(0.99, 2*(mean+1))
%       Poisson_normal
%               any condition ccg peak > mean + 3.1*sqrt(mean)
%       Poisson_chi_square_both
%               both condition, ccg peak > 0.5*chi2inv(0.99, 2*(mean+1))
%       Poisson_normal_both
%               both condition ccg peak > mean + 3.1*sqrt(mean)
%       Gaussian_convolve
%               cross-correlograms binned in 0.5-ms windows were convolved 
%               with a 7-ms standard deviation Gaussian window resulting 
%               in a predictor of the baseline rate
ccg_strct = struct(...
    'binsize', {}, ...
    'bandwidth', {}, ...
    'idx_MGB',{},...
    'idx_A1',{},...
    'unit_MGB',{},...
    'unit_A1',{},...
    'ccg_spon',{},...
    'ccg_filtered_spon', {}, ...
    'lag_spon', {}, ...
    'thresh_spon', {}, ...
    'ccg_stim',{},...
    'ccg_filtered_stim', {}, ...
    'lag_stim', {}, ...
    'thresh_stim', {}, ...
    'sig', {}, ...
    'methods', {});

nhwbins = hw / binsize;% number of bins in half window
c = 0;
for ii = 1:size(spktrain_ref_spon, 1)
    %fprintf('searching pairs for thalamic neuron #%d/%d\n', ii, size(spktrain_ref_spon, 1))
    % get shifted ref spike trains from -hw to hw
    shifted_spk_ref_spon = zeros(2*nhwbins + 1, size(spktrain_ref_spon, 2));
    shifted_spk_ref_stim = zeros(2*nhwbins + 1, size(spktrain_ref_stim, 2));
    for jj = -nhwbins:nhwbins
        n = jj + nhwbins + 1;
        shifted_spk_ref_spon(n,:) = circshift(spktrain_ref_spon(ii,:), jj);
        shifted_spk_ref_stim(n,:) = circshift(spktrain_ref_stim(ii,:), jj);
    end
    % exclude the edge of recordings
    shifted_spk_ref_spon = shifted_spk_ref_spon(:, nhwbins + 1 : end - nhwbins); 
    shifted_spk_ref_stim = shifted_spk_ref_stim(:, nhwbins + 1 : end - nhwbins);
    ccg_spon = shifted_spk_ref_spon * spktrain_target_spon(:, nhwbins + 1 : end - nhwbins)';
    ccg_stim = shifted_spk_ref_stim * spktrain_target_stim(:, nhwbins + 1 : end - nhwbins)';
    
    for jj = 1:size(ccg_spon,2)
        
        ccg_spon_tmp = ccg_spon(:,jj);
        ccg_stim_tmp = ccg_stim(:,jj);
        
        % if the maximum value exists outside of 1-5ms before cortical
        % spikes skip the pair
        [m, ~] = max(ccg_spon_tmp);
        if sum(ccg_spon_tmp([1:(nhwbins+1/binsize+1) (nhwbins+5/binsize+1):2*nhwbins+1]) == m) > 0
            continue
        end
        [m, ~] = max(ccg_stim_tmp);
        if sum(ccg_stim_tmp([1:(nhwbins+1/binsize+1) (nhwbins+5/binsize+1):2*nhwbins+1]) == m) > 0
            continue
        end
        
        % filter ccg with peak value within range
        if bandwidth(2) == Inf
            ccg_filtered_spon = highpass(ccg_spon_tmp, bandwidth(1),1000/binsize);
            ccg_filtered_stim = highpass(ccg_stim_tmp, bandwidth(1),1000/binsize);
        else
            ccg_filtered_spon = bandpass(ccg_spon_tmp, bandwidth,1000/binsize);
            ccg_filtered_stim = bandpass(ccg_stim_tmp, bandwidth,1000/binsize);
        end
        
        % get location of lag and significance of ccg
        [M_spon, Midx_spon] = max(ccg_filtered_spon);
        lag_spon = (Midx_spon - nhwbins -1) * binsize;
        [M_stim, Midx_stim] = max(ccg_filtered_stim);
        lag_stim = (Midx_stim - nhwbins -1) * binsize;
        
        if (lag_stim >= 1 && lag_stim <= 5) && (lag_spon >= 1 && lag_spon <= 5)
            % when lag of the peak is between 1-5ms
            sig = zeros(1, length(methods));
            thresh_stim = zeros(1, length(methods));
            thresh_spon = zeros(1, length(methods));
            for kk = 1:length(methods)
                method = methods{kk};
                switch method
                    case 'Normal_weak_both'
                        % The significance level was set at a probability of
                        % 0.2%, assuming a normal distribution in the baseline
                        % amplitudes after filtering. mean + 3.1std
                        thresh_spon(kk) = 3.1 * std( ccg_filtered_spon);
                        thresh_stim(kk) = 3.1 * std( ccg_filtered_stim);
                        sig(kk) = determine_sig(M_spon, thresh_spon(kk), ...
                            M_stim, thresh_stim(kk), ...
                            ccg_filtered_spon, ccg_filtered_stim, 'both');
                        
                    case 'Normal_strong'
                        % one condition (stim or spon) has ccg peak > 5std
                        thresh_spon(kk) = 5 * std( ccg_filtered_spon);
                        thresh_stim(kk) = 5 * std( ccg_filtered_stim);
                        sig(kk) = determine_sig(M_spon, thresh_spon(kk), ...
                            M_stim, thresh_stim(kk), ...
                            ccg_filtered_spon, ccg_filtered_stim, 'any');
                        
                    case 'Normal_strong_both'
                        % ccg peak > 5std under both conditions
                        thresh_spon(kk) = 5 * std( ccg_filtered_spon);
                        thresh_stim(kk) = 5 * std( ccg_filtered_stim);
                        sig(kk) = determine_sig(M_spon, thresh_spon(kk), ...
                            M_stim, thresh_stim(kk), ...
                            ccg_filtered_spon, ccg_filtered_stim, 'both');
                        
                    case 'Poisson_chi_square'
                        % one condition ccg peak > 0.5*chi2inv(0.99, 2*(mean+1))
                        thresh_spon(kk) = 0.5*chi2inv(0.99, 2*(mean(ccg_spon_tmp)+1));
                        thresh_stim(kk) = 0.5*chi2inv(0.99, 2*(mean(ccg_stim_tmp)+1));
                        sig(kk) = determine_sig(M_spon, thresh_spon(kk), ...
                            M_stim, thresh_stim(kk), ...
                            ccg_filtered_spon, ccg_filtered_stim, 'any');
                        
                    case 'Poisson_normal'
                        % one condition ccg peak > 3.1*sqrt(mean)
                        thresh_spon(kk) = 3.1*sqrt(mean(ccg_spon_tmp));
                        thresh_stim(kk) = 3.1*sqrt(mean(ccg_stim_tmp));
                        sig(kk) = determine_sig(M_spon, thresh_spon(kk), ...
                            M_stim, thresh_stim(kk), ...
                            ccg_filtered_spon, ccg_filtered_stim, 'any');
                        
                    case 'Poisson_chi_square_both'
                        % both condition ccg peak > 0.5*chi2inv(0.99, 2*(mean+1))
                        thresh_spon(kk) = 0.5*chi2inv(0.99, 2*(mean(ccg_spon_tmp)+1));
                        thresh_stim(kk) = 0.5*chi2inv(0.99, 2*(mean(ccg_stim_tmp)+1));
                        sig(kk) = determine_sig(M_spon, thresh_spon(kk), ...
                            M_stim, thresh_stim(kk), ...
                            ccg_filtered_spon, ccg_filtered_stim, 'both');
                        
                    case 'Poisson_normal_both'
                        % both condition ccg peak > 3.1*sqrt(mean)
                        thresh_spon(kk) = 3.1*sqrt(mean(ccg_spon_tmp));
                        thresh_stim(kk) = 3.1*sqrt(mean(ccg_stim_tmp));
                        sig(kk) = determine_sig(M_spon, thresh_spon(kk), ...
                            M_stim, thresh_stim(kk), ...
                            ccg_filtered_spon, ccg_filtered_stim, 'both');
                        
                    case 'Gaussian_convolve'
                        % significnact after subtracting baseline
                        ccg_data_spon = get_ccg_data( spktrain_ref_spon(ii,:), spktrain_target_spon(jj,:), binsize, hw, 0 );
                        ccg_data_stim = get_ccg_data( spktrain_ref_stim(ii,:), spktrain_target_stim(jj,:), binsize, hw, 0);
                        sig(kk) = ccg_data_spon.ab_comod_adj_peak_sig == 1 ...
                            || ccg_data_stim.ab_comod_adj_peak_sig == 1;
                end
            end
            
            if sum(sig) > 0
                c = c + 1;

                ccg_strct(c).binsize = binsize;
                ccg_strct(c).bandwidth = bandwidth;
                ccg_strct(c).idx_MGB = ii;
                ccg_strct(c).idx_A1 = jj;
                ccg_strct(c).unit_MGB = unit_ref(ii);
                ccg_strct(c).unit_A1 = unit_target(jj);
                ccg_strct(c).ccg_spon = ccg_spon_tmp;
                ccg_strct(c).ccg_filtered_spon = ccg_filtered_spon;
                ccg_strct(c).lag_spon = lag_spon;
                ccg_strct(c).thresh_spon = thresh_spon;
                ccg_strct(c).ccg_stim = ccg_stim_tmp;
                ccg_strct(c).ccg_filtered_stim = ccg_filtered_stim;
                ccg_strct(c).lag_stim = lag_stim;
                ccg_strct(c).thresh_stim = thresh_stim;
                ccg_strct(c).sig = sig;
                ccg_strct(c).methods = methods;
                
                if contains(method, 'Gaussian_convolve')
                    
                    ccg_strct(c).baseline_spon = ccg_data_spon.r_ab_smooth;
                    ccg_strct(c).ccov_spon = ccg_data_spon.q_ab;
                    ccg_strct(c).ccov_baseline_spon = ccg_data_spon.q_ab_smooth;
                    ccg_strct(c).conflimit_spon = ccg_data_spon.conf_limit;
                    ccg_strct(c).sig_spon = ccg_data_spon.ab_comod_adj_peak_sig;
                    ccg_strct(c).lag_spon2 = ccg_data_spon.ab_comod_adj_peak_delay;
                    ccg_strct(c).hw_spon = ccg_data_spon.ab_comod_adj_peak_hw;
                    
                    ccg_strct(c).baseline_stim = ccg_data_stim.r_ab_smooth;
                    ccg_strct(c).ccov_stim = ccg_data_stim.q_ab;
                    ccg_strct(c).ccov_baseline_stim = ccg_data_stim.q_ab_smooth;
                    ccg_strct(c).conflimit_stim = ccg_data_stim.conf_limit;
                    ccg_strct(c).sig_stim = ccg_data_stim.ab_comod_adj_peak_sig;
                    ccg_strct(c).lag_stim2 = ccg_data_stim.ab_comod_adj_peak_delay;
                    ccg_strct(c).hw_stim = ccg_data_stim.ab_comod_adj_peak_hw;
                end
                
            end
        end
    end
end

function sig = determine_sig(M_spon, thresh_spon, M_stim, thresh_stim, ...
    ccg_filtered_spon, ccg_filtered_stim, tag)

if strcmp(tag, 'both')
    condition1 = (M_spon > thresh_spon) && (M_stim > thresh_stim);
elseif strcmp(tag, 'any')
    condition1 = (M_spon > thresh_spon) || (M_stim > thresh_stim);
end

if condition1
    sig_spon = ccg_filtered_spon > thresh_spon(kk);
    sig_stim = ccg_filtered_stim > thresh_stim(kk);
    if sig_spon(Midx_spon-1) || sig_spon(Midx_spon+1) ...
            || sig_stim(Midx_stim-1) || sig_stim(Midx_stim+1)
        sig = 1;
    else
        sig = 0;
    end
else
    sig = 0;
end
