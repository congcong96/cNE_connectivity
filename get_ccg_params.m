function ccg_params = get_ccg_params( delay, q_ab, q_ab_smooth, conf_limit )
%
% ccg_params = get_ccg_params( delay, q_ab, q_ab_smooth, conf_limit )
%
%    INPUT --------------
%    delay : lag of cross-covariance function, in ms
%    q_ab : cross-covariance function
%    q_ab_smooth : smoothed cross-covariance function (estimate of slow co-modulation)
%    conf_limit : confidence limit used to assess significance
%
%   OUTPUT --------------
%   ccg_params is an output structure with >2 dozen param fields.
%   The following keys simplify the organization of these fields:
%   In #1-3, ab_ signifies influence of neuron_a -> neuron_b (reference -> target).
%   #1: ab_ccg_sig: True if any feature of q_ab exceeds q_ab +/- conf_limit
%   #2: ab_peak_sig: True if max of q_ab is greater than q_ab + conf_limit
%       ab_peak_val: Value of ab_peak_sig
%       ab_peak_delay: Latency of ab_peak_sig
%       ab_peak_hw: Half-width of ab_peak_sig
%   #3: Substitute trough for peak in #2 for params if q_ab is less than q_ab - conf_limit
%   #4: Substitute ab_comod_adj for ab_ in keys #2-3 for params reflecting q_ab_smooth, i.e., q_ab adjusted for slow co-modulation
%   #5: Substitute ba_ for ab_ in keys in #1-4 for params reflecting influence of neuron_b -> neuron_a (target -> reference)
%   #6: centroid and asymmetry_index reflect full q_ab as estimate bias of influence between neuron_a <-> neuron_b
%
% Example of full ccg_params =
%   struct with fields:
%                    ab_ccg_sig: 1
%                   ab_peak_sig: 1
%                   ab_peak_val: 4.6377e-04
%                 ab_peak_delay: 2
%                    ab_peak_hw: 1.5000
%                 ab_trough_sig: 0
%                 ab_trough_val: NaN
%               ab_trough_delay: NaN
%                  ab_trough_hw: NaN
%         ab_comod_adj_peak_sig: 1
%         ab_comod_adj_peak_val: 4.2817e-04
%       ab_comod_adj_peak_delay: 2
%          ab_comod_adj_peak_hw: 1.5000
%       ab_comod_adj_trough_sig: 0
%       ab_comod_adj_trough_val: NaN
%     ab_comod_adj_trough_delay: NaN
%        ab_comod_adj_trough_hw: NaN
%                    ba_ccg_sig: 1
%                   ba_peak_sig: 0
%                   ba_peak_val: NaN
%                 ba_peak_delay: NaN
%                    ba_peak_hw: NaN
%                 ba_trough_sig: 1
%                 ba_trough_val: -7.5096e-05
%               ba_trough_delay: 2
%                  ba_trough_hw: 2
%         ba_comod_adj_peak_sig: 0
%         ba_comod_adj_peak_val: NaN
%       ba_comod_adj_peak_delay: NaN
%          ba_comod_adj_peak_hw: NaN
%       ba_comod_adj_trough_sig: 1
%       ba_comod_adj_trough_val: -1.0395e-04
%     ba_comod_adj_trough_delay: 2
%        ba_comod_adj_trough_hw: 3.5000
%                      centroid: 4.1608
%               asymmetry_index: 0.4822

%% Initialize vars
q_ab_comod_adj = q_ab - q_ab_smooth; % xcov: adjusted to remove influence of slow comodulation
binwidth = delay(2)-delay(1); % Assuming equal bins
ccg_center_ms = 5;
cutoff_val = 0.5; % For defining peak/trough halfwidth

%% A->B params -------------------------------------------
idx = find( delay >= 1 & delay < ccg_center_ms );
q_tmp = q_ab(idx);
q_tmp_comod_adj = q_ab_comod_adj(idx);

% 1/5 Sig A->B CCG? - - - - - - - - - - - - - - - - - - - -
[~, midx] = max(q_ab);
if ~sum(idx == midx) || isempty( find( diff( find( q_tmp > conf_limit | q_tmp < -conf_limit ) ) == 1, 1 ) )
    ab_ccg_sig = false;
else
    ab_ccg_sig = true;
end

% 2/5 Sig A->B peak? - - - - - - - - - - - - - - - - - - - -
q_full = q_ab;
if isempty( find( diff( find( q_tmp > conf_limit ) ) == 1, 1 ) )
    ab_peak_sig = false;
    ab_peak_val = NaN;
    ab_peak_delay = NaN;
    ab_peak_hw = NaN;
else
    ab_peak_sig = true;
    q_tmp0 = q_tmp; 
    sig_bins = get_sig_bins( q_tmp > conf_limit ); % Restrict peak calc to 2 adjacent sig bins
    q_tmp0(~sig_bins) = 0; 
    [ M, I ] = max( q_tmp0 );
    ab_peak_val = M;
    ab_peak_delay = I * binwidth + 1;
    idx_peak = I + numel(find( delay < 0 )) + 2;
    x1 = find( q_full(1:idx_peak) < M*cutoff_val , 1, 'last' );
    x2 = idx_peak + find( q_full(idx_peak+1:end) < M*cutoff_val , 1 );
    ab_peak_hw = (x2-x1-1)*binwidth;
end

% 3/5 Sig A->B trough? - - - - - - - - - - - - - - - - - - - -
q_full = q_ab;
if isempty( find( diff( find( q_tmp < -conf_limit ) ) == 1, 1 ) )
    ab_trough_sig = false;
    ab_trough_val = NaN;
    ab_trough_delay = NaN;
    ab_trough_hw = NaN;
else
    ab_trough_sig = true;
    q_tmp0 = q_tmp; 
    sig_bins = get_sig_bins( q_tmp < -conf_limit ); % Restrict trough calc to 2 adjacent sig bins
    q_tmp0(~sig_bins) = 0; 
    [ M, I ] = min( q_tmp0 );
    ab_trough_val = M;
    ab_trough_delay = I * binwidth + 1;
    idx_peak = I + numel(find( delay < 0 )) + 2;
    x1 = find( q_full(1:idx_peak) > M*cutoff_val, 1, 'last');
    x2 = idx_peak + find( q_full(idx_peak+1:end) > M*cutoff_val, 1, 'first' );
    ab_trough_hw = (x2-x1-1)*binwidth;
end

% 4/5 Sig A->B peak after slow comodulation adjustment? - - - - - - - - - -
q_full = q_ab_comod_adj;
[~, midx] = max(q_ab_comod_adj);
if ~sum(idx == midx) || isempty( find( diff( find( q_tmp_comod_adj > conf_limit ) ) == 1, 1 ) )
    ab_comod_adj_peak_sig = false;
    ab_comod_adj_peak_val = NaN;
    ab_comod_adj_peak_delay = NaN;
    ab_comod_adj_peak_hw = NaN;
else
    q_tmp_comod_adj0 = q_tmp_comod_adj; 
    sig_bins = get_sig_bins( q_tmp_comod_adj > conf_limit ); % Restrict peak calc to 2 adjacent sig bins
    q_tmp_comod_adj0(~sig_bins) = 0; 
    [ M, I ] = max( q_tmp_comod_adj0 );
    idx_peak = I + numel(find( delay < 0 ))+2;
    if idx_peak ~= midx
        ab_comod_adj_peak_sig = false;
        ab_comod_adj_peak_val = NaN;
        ab_comod_adj_peak_delay = NaN;
        ab_comod_adj_peak_hw = NaN;
    else
        ab_comod_adj_peak_sig = true;
        ab_comod_adj_peak_val = M;
        ab_comod_adj_peak_delay = I * binwidth + 1;
        x1 = find( q_full(1:idx_peak) < M*cutoff_val, 1, 'last' ) ;
        x2 = idx_peak + find( q_full(idx_peak+1:end) < M*cutoff_val, 1, 'first' );
        ab_comod_adj_peak_hw = (x2-x1)*binwidth;
    end
end

% 5/5 Sig A->B trough after slow comodulation adjustment? - - - - - - - - - -
q_full = q_ab_comod_adj;
if isempty( find( diff( find( q_tmp_comod_adj < -conf_limit ) ) == 1, 1) )
    ab_comod_adj_trough_sig = false;
    ab_comod_adj_trough_val = NaN;
    ab_comod_adj_trough_delay = NaN;
    ab_comod_adj_trough_hw = NaN;
else
    ab_comod_adj_trough_sig = true;
    q_tmp_comod_adj0 = q_tmp_comod_adj; 
    sig_bins = get_sig_bins( q_tmp_comod_adj < -conf_limit ); % Restrict trough calc to 2 adjacent sig bins
    q_tmp_comod_adj0(~sig_bins) = 0; 
    [ M, I ] = min( q_tmp_comod_adj0 );
    ab_comod_adj_trough_val = M;
    ab_comod_adj_trough_delay = I * binwidth;
    idx_peak = I + numel(find( delay < 0 ));
    x1 = find( q_full(1:idx_peak) > M*cutoff_val , 1, 'last' );
    x2 = idx_peak + find( q_full(idx_peak+1:end) > M*cutoff_val , 1 );
    ab_comod_adj_trough_hw = (x2-x1)*binwidth;
end

%% B->A params -------------------------------------------
idx = find( delay < 0 & delay >= -ccg_center_ms );
q_tmp = fliplr( q_ab(idx) );
q_tmp_comod_adj = fliplr( q_ab_comod_adj(idx) );

% 1/5 Sig B->A CCG? - - - - - - - - - - - - - - - - - - - -
if isempty( find( diff( find( q_tmp >= conf_limit | q_tmp <= -conf_limit ) ) == 1, 1 ) )
    ba_ccg_sig = false;
else
    ba_ccg_sig = true;
end

% 2/5 Sig B->A peak? - - - - - - - - - - - - - - - - - - - -
q_full = fliplr( q_ab );
if isempty( find( diff( find( q_tmp >= conf_limit ) ) == 1, 1 ) )
    ba_peak_sig = false;
    ba_peak_val = NaN;
    ba_peak_delay = NaN;
    ba_peak_hw = NaN;
else
    ba_peak_sig = true;
    q_tmp0 = q_tmp; 
    sig_bins = get_sig_bins( q_tmp >= conf_limit ); % Restrict peak calc to 2 adjacent sig bins
    q_tmp0(~sig_bins) = 0; 
    [ M, I ] = max( q_tmp0 );
    ba_peak_val = M;
    ba_peak_delay = I * binwidth;
    idx_peak = I + numel(find( delay < 0 ));
    x1 = find( q_full(1:idx_peak) < M*cutoff_val , 1, 'last' );
    x2 = idx_peak + find( q_full(idx_peak+1:end) < M*cutoff_val , 1 );
    ba_peak_hw = (x2-x1)*binwidth;
end

% 3/5 Sig B->A trough? - - - - - - - - - - - - - - - - - - - -
q_full = fliplr( q_ab );
if isempty( find( diff( find( q_tmp <= -conf_limit ) ) == 1, 1 ) )
    ba_trough_sig = false;
    ba_trough_val = NaN;
    ba_trough_delay = NaN;
    ba_trough_hw = NaN;
else
    ba_trough_sig = true;
    q_tmp0 = q_tmp; 
    sig_bins = get_sig_bins( q_tmp <= -conf_limit ); % Restrict peak calc to 2 adjacent sig bins
    q_tmp0(~sig_bins) = 0; 
    [ M, I ] = min( q_tmp0 );
    ba_trough_val = M;
    ba_trough_delay = I * binwidth;
    idx_peak = I + numel(find( delay < 0 ));
    x1 = find( q_full(1:idx_peak) > M*cutoff_val , 1, 'last' );
    x2 = idx_peak + find( q_full(idx_peak+1:end) > M*cutoff_val , 1 );
    ba_trough_hw = (x2-x1)*binwidth;
end

% 4/5 Sig B->A peak after slow comodulation adjustment? - - - - - - - - - -
q_full = fliplr( q_ab_comod_adj );
if isempty( find( diff( find( q_tmp_comod_adj >= conf_limit ) ) == 1, 1 ) )
    ba_comod_adj_peak_sig = false;
    ba_comod_adj_peak_val = NaN;
    ba_comod_adj_peak_delay = NaN;
    ba_comod_adj_peak_hw = NaN;
else
    ba_comod_adj_peak_sig = true;
    q_tmp_comod_adj0 = q_tmp_comod_adj; 
    sig_bins = get_sig_bins( q_tmp_comod_adj >= conf_limit ); % Restrict peak calc to 2 adjacent sig bins
    q_tmp_comod_adj0(~sig_bins) = 0; 
    [ M, I ] = max( q_tmp_comod_adj0 );
    ba_comod_adj_peak_val = M;
    ba_comod_adj_peak_delay = I * binwidth;
    idx_peak = I + numel(find( delay < 0 ));
    x1 = find( q_full(1:idx_peak) < M*cutoff_val , 1, 'last' );
    x2 = idx_peak + find( q_full(idx_peak+1:end) < M*cutoff_val , 1 );
    ba_comod_adj_peak_hw = (x2-x1)*binwidth;
end

% 5/5 Sig B->A trough after slow comodulation adjustment? - - - - - - - - - -
q_full = fliplr( q_ab_comod_adj );
if isempty( find( diff( find( q_tmp_comod_adj <= -conf_limit ) ) == 1, 1 ) )
    ba_comod_adj_trough_sig = false;
    ba_comod_adj_trough_val = NaN;
    ba_comod_adj_trough_delay = NaN;
    ba_comod_adj_trough_hw = NaN;
else
    ba_comod_adj_trough_sig = true;
    q_tmp_comod_adj0 = q_tmp_comod_adj; 
    sig_bins = get_sig_bins( q_tmp_comod_adj <= -conf_limit ); % Restrict trough calc to 2 adjacent sig bins
    q_tmp_comod_adj0(~sig_bins) = 0; 
    [ M, I ] = min( q_tmp_comod_adj0 );
    ba_comod_adj_trough_val = M;
    ba_comod_adj_trough_delay = I * binwidth;
    idx_peak = I + numel(find( delay < 0 ));
    x1 = find( q_full(1:idx_peak) > M*cutoff_val , 1, 'last' );
    x2 = idx_peak + find( q_full(idx_peak+1:end) > M*cutoff_val , 1 );
    ba_comod_adj_trough_hw = (x2-x1)*binwidth;
end


%% Centroid & Asymmetry Index

ccg_center_window = round( ccg_center_ms / binwidth );

idx_0 = find(delay == 0);
idx_neg = find( delay <= -ccg_center_window, 1, 'last' );
idx_pos = find( delay >= ccg_center_window, 1, 'first' );

idx_center = idx_neg:idx_pos;
qab_center = abs( q_ab(idx_center) );
delay_center = delay(idx_center);

% Centroid of central portion of cross-covariance function (abs value to account for inhibited spikes)
centroid = sum( qab_center .* delay_center ) / sum(qab_center);

% Asymmetry Index from cross-covariance function (abs value to account for inhibited spikes)
vals_right = abs( q_ab(idx_0+1:idx_pos) );
vals_left = abs( q_ab(idx_neg:idx_0-1) );
asymmetry_index = ( sum(vals_right) - sum(vals_left) ) / ( sum(vals_right) + sum(vals_left) );

%% Package output params into structure

% A->B params
ccg_params.ab_ccg_sig = ab_ccg_sig;

ccg_params.ab_peak_sig = ab_peak_sig;
ccg_params.ab_peak_val = ab_peak_val;
ccg_params.ab_peak_delay = ab_peak_delay;
ccg_params.ab_peak_hw = ab_peak_hw;

ccg_params.ab_trough_sig = ab_trough_sig;
ccg_params.ab_trough_val = ab_trough_val;
ccg_params.ab_trough_delay = ab_trough_delay;
ccg_params.ab_trough_hw = ab_trough_hw;

ccg_params.ab_comod_adj_peak_sig = ab_comod_adj_peak_sig;
ccg_params.ab_comod_adj_peak_val = ab_comod_adj_peak_val;
ccg_params.ab_comod_adj_peak_delay = ab_comod_adj_peak_delay;
ccg_params.ab_comod_adj_peak_hw = ab_comod_adj_peak_hw;

ccg_params.ab_comod_adj_trough_sig = ab_comod_adj_trough_sig;
ccg_params.ab_comod_adj_trough_val = ab_comod_adj_trough_val;
ccg_params.ab_comod_adj_trough_delay = ab_comod_adj_trough_delay;
ccg_params.ab_comod_adj_trough_hw = ab_comod_adj_trough_hw;

% B->A params
ccg_params.ba_ccg_sig = ba_ccg_sig;

ccg_params.ba_peak_sig = ba_peak_sig;
ccg_params.ba_peak_val = ba_peak_val;
ccg_params.ba_peak_delay = ba_peak_delay;
ccg_params.ba_peak_hw = ba_peak_hw;

ccg_params.ba_trough_sig = ba_trough_sig;
ccg_params.ba_trough_val = ba_trough_val;
ccg_params.ba_trough_delay = ba_trough_delay;
ccg_params.ba_trough_hw = ba_trough_hw;

ccg_params.ba_comod_adj_peak_sig = ba_comod_adj_peak_sig;
ccg_params.ba_comod_adj_peak_val = ba_comod_adj_peak_val;
ccg_params.ba_comod_adj_peak_delay = ba_comod_adj_peak_delay;
ccg_params.ba_comod_adj_peak_hw = ba_comod_adj_peak_hw;

ccg_params.ba_comod_adj_trough_sig = ba_comod_adj_trough_sig;
ccg_params.ba_comod_adj_trough_val = ba_comod_adj_trough_val;
ccg_params.ba_comod_adj_trough_delay = ba_comod_adj_trough_delay;
ccg_params.ba_comod_adj_trough_hw = ba_comod_adj_trough_hw;

% Centroid & Assymetry Index
ccg_params.centroid = centroid;
ccg_params.asymmetry_index = asymmetry_index;

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

end
