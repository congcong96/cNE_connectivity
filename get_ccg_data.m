function ccg_data = get_ccg_data( spktrain_a, spktrain_b, binsize, maxlag, plotopt )

% INPUT --------------
% spktrain_a: reference unit, time-binned spike count vector
% spktrain_b: target unit, time-binned spike count vector
% binsize: time binwidth for spike trains in ms
% maxlag: CCG analysis window in ms
% plotopt: plot CCG data if True.
% 
% OUTPUT --------------
% ccg_data: Structure with all relevant CCG data, including param estimates from get_ccg_params 

% Note: Why not use Matlab xcov function? Because the 'baseline' or null
% estimate is shifted by non-zero firing for pairs with sig CCG (increased by excitation)
% [c,lags] = xcov( spktrain_b(:)', spktrain_a(:)', maxlag_ms/dt, 'biased' );
% plot(lags,c)


% Define filter for slow co-modulation estimate - - - - - - - - - - - - - -
% "cross-correlograms binned in 0.5-ms windows were convolved with a 7-ms standard deviation Gaussian window 
% resulting in a predictor of the baseline rate"
% Senzai et al. (2019) Layer-Specific Physiological Features and Interlaminar Interactions in the Primary Visual Cortex of the Mouse
sigma_ms = 7; 
hsize = [(sigma_ms^2)/binsize 1];
sigma = sigma_ms/binsize; 
h = fspecial('gaussian',hsize,sigma);

spktrain_a_smooth = conv( spktrain_a, h, 'same' );
spktrain_b_smooth = conv( spktrain_b, h, 'same' );

[r_ab, delay] = xcorr( spktrain_b(:)', spktrain_a(:)', maxlag/binsize );
r_ab = real(r_ab);
r_ab(r_ab < 0) = 0;
delay = delay .* binsize; % delay in ms, not bin number

r_ab_smooth  = xcorr( spktrain_b_smooth(:)', spktrain_a_smooth(:)', maxlag/binsize );
r_ab_smooth = real(r_ab_smooth);
r_ab_smooth(r_ab_smooth < 0) = 0;

% upper bound
high_bound = poissinv(0.99, r_ab_smooth);
ccg_params = get_ccg_params( delay, r_ab, r_ab_smooth, high_bound);

% Store output data in struct
ccg_data.dt = binsize;
ccg_data.delay = delay;
ccg_data.r_ab = r_ab;
ccg_data.r_ab_smooth = r_ab_smooth;
ccg_data.high_bound = high_bound;
%...ccg_params:
ccg_data.ab_peak_sig = ccg_params.ab_peak_sig;
ccg_data.ab_peak_delay = ccg_params.ab_peak_delay;
ccg_data.ab_peak_hw = ccg_params.ab_peak_hw;
end