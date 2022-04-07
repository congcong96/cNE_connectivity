function ccg_params = get_ccg_params( delay, r_ab, r_ab_smooth, high_bound )
%
% ccg_params = get_ccg_params( delay, q_ab, q_ab_smooth, conf_limit )
%
%    INPUT --------------
%    delay : lag of cross-covariance function, in ms
%    r_ab : cross-correlation function
%    high_bound : high bound to assess significance
%
%   OUTPUT --------------
%   ccg_params is an output structure with the following param fields.
%       ab_peak_sig: True if max of q_ab is greater than higher_bound
%       ab_peak_val: Value of ab_peak_sig
%       ab_peak_delay: Latency of ab_peak_sig
%       ab_peak_hw: Half-width of ab_peak_sig

%% Initialize vars
binwidth = delay(2)-delay(1); % Assuming equal bins
max_delay = 5;
cutoff_val = 0.5; % For defining peak/trough halfwidth

%% A->B params -------------------------------------------
idx = find( delay >= 1 & delay < max_delay );
r_tmp = r_ab(idx);
conf_limit = high_bound(idx);
[~, Midx] = max(r_ab-r_ab_smooth);
if ~sum(idx == Midx)
    % if the peak after correction for baseline is not in the time window
    ab_peak_sig = false;
    ab_peak_val = NaN;
    ab_peak_delay = NaN;
    ab_peak_hw = NaN;
else
    
    ab_peak_sig = true;
    r_tmp0 = r_tmp;
    sig_bins = get_sig_bins( r_tmp > conf_limit ); % Restrict peak calc to 2 adjacent sig bins
    if sum(sig_bins) == 0
        ab_peak_sig = false;
        ab_peak_val = NaN;
        ab_peak_delay = NaN;
        ab_peak_hw = NaN;
    else
        r_tmp0(~sig_bins) = 0;
        [ M, I ] = max( r_tmp0 );
        ab_peak_val = M;
        ab_peak_delay = I * binwidth + 1;
        idx_peak = I + numel(find( delay < 0 )) + 2;
        x1 = find( r_ab(1:idx_peak) < M*cutoff_val , 1, 'last' );
        x2 = idx_peak + find( r_ab(idx_peak+1:end) < M*cutoff_val , 1 );
        ab_peak_hw = (x2-x1-1)*binwidth;
    end
   
    
end

%% Package output params into structure
ccg_params.ab_peak_sig = ab_peak_sig;
ccg_params.ab_peak_val = ab_peak_val;
ccg_params.ab_peak_delay = ab_peak_delay;
ccg_params.ab_peak_hw = ab_peak_hw;

end

% Helper funciton --------------------------------------
function sig_bins = get_sig_bins(vec)
% There must be a vectorized implementation??
sig_bins = zeros(size(vec));
for qq = 1:numel(vec)
    if qq == 1
        sig_bin = vec(qq) && vec(qq+1);
    elseif qq == numel(vec)
        sig_bin = vec(qq) && vec(qq-1);
    else
        sig_bin = vec(qq) && (vec(qq-1) || vec(qq+1));
    end
    sig_bins(qq) = sig_bin;
end

end
