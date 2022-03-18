function ccg_strct = find_connected_pairs(...
    spktrain_MGB_spon, spktrain_A1_spon, ...
    spktrain_MGB_dmr, spktrain_A1_dmr, ...
    hw, binsize, unit_MGB, unit_A1, bandwidth, methods)
% 
% Input:
%   spktrain_MGB_spon: m x n matrix, m is the number of neurons, n is the
%   number of time bins
%       spontaneous thalamic spike trians binned at binsize
%   spktrain_A1_spon: m x n matrix
%       spontaneous cortical spike trians binned at binsize
%   spktrain_MGB_dmr: m x n matrix
%       thalamic spike trians in response to dmr binned at binsize
%   spktrain_MGB_dmr: m x n matrix
%       cortical spike trians in response to dmr binned at binsize
%   hw: integer
%       half window size, in ms
%   binsize: in ms
%   Unit_MGB: vector of length m
%       unit label for thalamic neurons
%   Unit_A1: vector of length m
%       unit label for cortical neurons
%   bandwith: bandwith to filter ccg
%   method: cell with string indicating method to test significance of ccg
%       Normal_weak_both
%               The significance level was set at a probability of
%               0.2%, assuming a normal distribution in the baseline
%               amplitudes after filtering. mean + 3.1std
%       Normal_strong
%               at least one condition (dmr or spon) has ccg peak > 5std
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
%
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
    'ccg_dmr',{},...
    'ccg_filtered_dmr', {}, ...
    'lag_dmr', {}, ...
    'thresh_dmr', {}, ...
    'sig', {}, ...
    'methods', {});

nhwbins = hw / binsize;% number of bins in half window
c = 0;
for ii = 1:size(spktrain_MGB_spon, 1)
    %fprintf('searching pairs for thalamic neuron #%d/%d\n', ii, size(spktrain_MGB_spon, 1))
    % get shifted MGB spike trains from -hw to hw
    shifted_spk_MGB_spon = zeros(2*nhwbins + 1, size(spktrain_MGB_spon, 2));
    shifted_spk_MGB_dmr = zeros(2*nhwbins + 1, size(spktrain_MGB_dmr, 2));
    for jj = -nhwbins:nhwbins
        n = jj + nhwbins + 1;
        shifted_spk_MGB_spon(n,:) = circshift(spktrain_MGB_spon(ii,:), jj);
        shifted_spk_MGB_dmr(n,:) = circshift(spktrain_MGB_dmr(ii,:), jj);
    end
    % exclude the edge of recordings
    shifted_spk_MGB_spon = shifted_spk_MGB_spon(:, nhwbins + 1 : end - nhwbins); 
    shifted_spk_MGB_dmr = shifted_spk_MGB_dmr(:, nhwbins + 1 : end - nhwbins);
    ccg_spon = shifted_spk_MGB_spon * spktrain_A1_spon(:, nhwbins + 1 : end - nhwbins)';
    ccg_dmr = shifted_spk_MGB_dmr * spktrain_A1_dmr(:, nhwbins + 1 : end - nhwbins)';
    
    for jj = 1:size(ccg_spon,2)
        
        ccg_spon_tmp = ccg_spon(:,jj);
        ccg_dmr_tmp = ccg_dmr(:,jj);
        
        % if the maximum value exists outside of 1-5ms before cortical
        % spikes skip the pair
        [m, ~] = max(ccg_spon_tmp);
        if sum(ccg_spon_tmp([1:(nhwbins+1/binsize+1) (nhwbins+5/binsize+1):2*nhwbins+1]) == m) > 0
            continue
        end
        [m, ~] = max(ccg_dmr_tmp);
        if sum(ccg_dmr_tmp([1:(nhwbins+1/binsize+1) (nhwbins+5/binsize+1):2*nhwbins+1]) == m) > 0
            continue
        end
        
        % filter ccg with peak value within range
        if bandwidth(2) == Inf
            ccg_filtered_spon = highpass(ccg_spon_tmp, bandwidth(1),1000/binsize);
            ccg_filtered_dmr = highpass(ccg_dmr_tmp, bandwidth(1),1000/binsize);
        else
            ccg_filtered_spon = bandpass(ccg_spon_tmp, bandwidth,1000/binsize);
            ccg_filtered_dmr = bandpass(ccg_dmr_tmp, bandwidth,1000/binsize);
        end
        
        % get location of lag and significance of ccg
        [M_spon, Midx_spon] = max(ccg_filtered_spon);
        lag_spon = (Midx_spon - nhwbins -1) * binsize;
        [M_dmr, Midx_dmr] = max(ccg_filtered_dmr);
        lag_dmr = (Midx_dmr - nhwbins -1) * binsize;
        
        if (lag_dmr >= 1 && lag_dmr <= 5) && (lag_spon >= 1 && lag_spon <= 5)
            % when lag of the peak is between 1-5ms
            sig = zeros(1, length(methods));
            thresh_dmr = zeros(1, length(methods));
            thresh_spon = zeros(1, length(methods));
            for kk = 1:length(methods)
                method = methods{kk};
                switch method
                    case 'Normal_weak_both'
                        % The significance level was set at a probability of
                        % 0.2%, assuming a normal distribution in the baseline
                        % amplitudes after filtering. mean + 3.1std
                        thresh_spon(kk) = 3.1 * std( ccg_filtered_spon);
                        thresh_dmr(kk) = 3.1 * std( ccg_filtered_dmr);
                        sig(kk) = (M_spon > thresh_spon(kk)) ...
                            && (M_dmr > thresh_dmr(kk));
                    case 'Normal_strong'
                        % one condition (dmr or spon) has ccg peak > 5std
                        thresh_spon(kk) = 5 * std( ccg_filtered_spon);
                        thresh_dmr(kk) = 5 * std( ccg_filtered_dmr);
                        sig(kk) = (M_spon > thresh_spon(kk)) ...
                                || (M_dmr > thresh_dmr(kk));
                    case 'Normal_strong_both'
                        % ccg peak > 5std under both conditions
                        thresh_spon(kk) = 5 * std( ccg_filtered_spon);
                        thresh_dmr(kk) = 5 * std( ccg_filtered_dmr);
                        sig(kk) = (M_spon > thresh_spon(kk)) ...
                            && (M_dmr > thresh_dmr(kk));
                    case 'Poisson_chi_square'
                        % one condition ccg peak > 0.5*chi2inv(0.99, 2*(mean+1))
                        thresh_spon(kk) = 0.5*chi2inv(0.99, 2*(mean(ccg_spon_tmp)+1));
                        thresh_dmr(kk) = 0.5*chi2inv(0.99, 2*(mean(ccg_dmr_tmp)+1));
                        sig(kk) = ((M_spon + mean(ccg_spon_tmp)) > thresh_spon(kk)) ...
                            || ((M_dmr + mean(ccg_dmr_tmp)) > thresh_dmr(kk));
                    case 'Poisson_normal'
                        % one condition ccg peak > 3.1*sqrt(mean)
                        thresh_spon(kk) = 3.1*sqrt(mean(ccg_spon_tmp));
                        thresh_dmr(kk) = 3.1*sqrt(mean(ccg_dmr_tmp));
                        sig(kk) = (M_spon  > thresh_spon(kk)) ...
                            || (M_dmr  > thresh_dmr(kk));
                    case 'Poisson_chi_square_both'
                        % both condition ccg peak > 0.5*chi2inv(0.99, 2*(mean+1))
                        thresh_spon(kk) = 0.5*chi2inv(0.99, 2*(mean(ccg_spon_tmp)+1));
                        thresh_dmr(kk) = 0.5*chi2inv(0.99, 2*(mean(ccg_dmr_tmp)+1));
                        sig(kk) = ((M_spon + mean(ccg_spon_tmp)) > thresh_spon(kk)) ...
                            && ((M_dmr + mean(ccg_dmr_tmp)) > thresh_dmr(kk));
                    case 'Poisson_normal_both'
                        % both condition ccg peak > 3.1*sqrt(mean)
                        thresh_spon(kk) = 3.1*sqrt(mean(ccg_spon_tmp));
                        thresh_dmr(kk) = 3.1*sqrt(mean(ccg_dmr_tmp));
                        sig(kk) = (M_spon  > thresh_spon(kk)) ...
                            && (M_dmr  > thresh_dmr(kk));
                end
            end
            
            if sum(sig) > 0
                c = c + 1;

                ccg_strct(c).binsize = binsize;
                ccg_strct(c).bandwidth = bandwidth;
                ccg_strct(c).idx_MGB = ii;
                ccg_strct(c).idx_A1 = jj;
                ccg_strct(c).unit_MGB = unit_MGB(ii);
                ccg_strct(c).unit_A1 = unit_A1(jj);
                ccg_strct(c).ccg_spon = ccg_spon_tmp;
                ccg_strct(c).ccg_filtered_spon = ccg_filtered_spon;
                ccg_strct(c).lag_spon = lag_spon;
                ccg_strct(c).thresh_spon = thresh_spon;
                ccg_strct(c).ccg_dmr = ccg_dmr_tmp;
                ccg_strct(c).ccg_filtered_dmr = ccg_filtered_dmr;
                ccg_strct(c).lag_dmr = lag_dmr;
                ccg_strct(c).thresh_dmr = thresh_dmr;
                ccg_strct(c).sig = sig;
                ccg_strct(c).methods = methods;
            end
        end
    end
end
