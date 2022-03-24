function ccg_data = get_ccg_data( spktrain_a, spktrain_b, binwidth_ms, maxlag_ms, plotopt )

% INPUT --------------
% spktrain_a: reference unit, time-binned spike count vector
% spktrain_b: target unit, time-binned spike count vector
% binwidth_ms: time binwidth for spike trains in ms
% maxlag_ms: CCG analysis window in ms
% plotopt: plot CCG data if True.
% 
% OUTPUT --------------
% ccg_data: Structure with all relevant CCG data, including param estimates from get_ccg_params 

% Note: Why not use Matlab xcov function? Because the 'baseline' or null
% estimate is shifted by non-zero firing for pairs with sig CCG (increased by excitation)
% [c,lags] = xcov( spktrain_b(:)', spktrain_a(:)', maxlag_ms/dt, 'biased' );
% plot(lags,c)

if ~exist('plotopt','var')
    plotopt = 0;
end

% Define filter for slow co-modulation estimate - - - - - - - - - - - - - -
% "cross-correlograms binned in 0.5-ms windows were convolved with a 7-ms standard deviation Gaussian window 
% resulting in a predictor of the baseline rate"
% Senzai et al. (2019) Layer-Specific Physiological Features and Interlaminar Interactions in the Primary Visual Cortex of the Mouse
sigma_ms = 7; 
hsize = [(sigma_ms^2)/binwidth_ms 1];
sigma = sigma_ms/binwidth_ms; 
h = fspecial('gaussian',hsize,sigma);

spktrain_a_smooth = conv( spktrain_a, h, 'same' );
spktrain_b_smooth = conv( spktrain_b, h, 'same' );

stimdur = numel(spktrain_a) * binwidth_ms;

nspk_a = sum(spktrain_a);
nspk_b = sum(spktrain_b);

[r_ab, delay] = xcorr( spktrain_b(:)', spktrain_a(:)', maxlag_ms/binwidth_ms );
r_ab = real(r_ab);
r_ab(r_ab < 0) = 0;
delay = delay .* binwidth_ms; % delay in ms, not bin number

r_ab_smooth  = xcorr( spktrain_b_smooth(:)', spktrain_a_smooth(:)', maxlag_ms/binwidth_ms );
r_ab_smooth = real(r_ab_smooth);
r_ab_smooth(r_ab_smooth < 0) = 0;

[q_ab, conf_limit, stimdur] = get_cross_covariance(r_ab, nspk_a, nspk_b, binwidth_ms, stimdur );
[q_ab_smooth, conf_limit_smooth ] = get_cross_covariance(r_ab_smooth, nspk_a, nspk_b, binwidth_ms, stimdur );

rho = corrstrength(delay, r_ab, nspk_a, nspk_b, stimdur); % correlation coefficient

ccg_params = get_ccg_params( delay, q_ab, q_ab_smooth, conf_limit );

% Store output data in struct
ccg_data.dt = binwidth_ms;
ccg_data.nspk_a = nspk_a;
ccg_data.nspk_a = nspk_b;
ccg_data.delay = delay;
ccg_data.r_ab = r_ab;
ccg_data.q_ab = q_ab;
ccg_data.r_ab_smooth = r_ab_smooth;
ccg_data.q_ab_smooth = q_ab_smooth;
ccg_data.conf_limit = conf_limit;
ccg_data.rho = rho;
%...ccg_params:
ccg_data.ab_ccg_sig = ccg_params.ab_ccg_sig;
ccg_data.ab_peak_sig = ccg_params.ab_peak_sig;
ccg_data.ab_peak_delay = ccg_params.ab_peak_delay;
ccg_data.ab_peak_hw = ccg_params.ab_peak_hw;
ccg_data.ab_trough_sig = ccg_params.ab_trough_sig;
ccg_data.ab_trough_delay = ccg_params.ab_trough_delay;
ccg_data.ab_trough_hw = ccg_params.ab_trough_hw;
ccg_data.ab_comod_adj_peak_sig = ccg_params.ab_comod_adj_peak_sig;
ccg_data.ab_comod_adj_peak_delay = ccg_params.ab_comod_adj_peak_delay;
ccg_data.ab_comod_adj_peak_hw = ccg_params.ab_comod_adj_peak_hw;
ccg_data.ab_comod_adj_trough_sig = ccg_params.ab_comod_adj_trough_sig;
ccg_data.ab_comod_adj_trough_delay = ccg_params.ab_comod_adj_trough_delay;
ccg_data.ab_comod_adj_trough_hw = ccg_params.ab_comod_adj_trough_hw;
ccg_data.ba_ccg_sig = ccg_params.ba_ccg_sig;
ccg_data.ba_peak_sig = ccg_params.ba_peak_sig;
ccg_data.ba_peak_delay = ccg_params.ba_peak_delay;
ccg_data.ba_peak_hw = ccg_params.ba_peak_hw;
ccg_data.ba_trough_sig = ccg_params.ba_trough_sig;
ccg_data.ba_trough_delay = ccg_params.ba_trough_delay;
ccg_data.ba_trough_hw = ccg_params.ba_trough_hw;
ccg_data.ba_comod_adj_peak_sig = ccg_params.ba_comod_adj_peak_sig;
ccg_data.ba_comod_adj_peak_delay = ccg_params.ba_comod_adj_peak_delay;
ccg_data.ba_comod_adj_peak_hw = ccg_params.ba_comod_adj_peak_hw;
ccg_data.ba_comod_adj_trough_sig = ccg_params.ba_comod_adj_trough_sig;
ccg_data.ba_comod_adj_trough_delay = ccg_params.ba_comod_adj_trough_delay;
ccg_data.ba_comod_adj_trough_hw = ccg_params.ba_comod_adj_trough_hw;
ccg_data.centroid = ccg_params.centroid;
ccg_data.asymmetry_index = ccg_params.asymmetry_index;

if plotopt == 1
    
    %%
    close all
    figure('units','inches','position',[.1 .1 12 8])
    subplot(221);
    bar(delay, r_ab, 'k', 'BarWidth', 1);
    line([0 0],[min(ylim) max(ylim)],'color','r','linestyle',':')
    xlabel('Time (ms)')
    ylabel('Counts')
    set(gca,'color','none','box','off')
    
    subplot(222);
    bar(delay, r_ab, 'k', 'BarWidth', 1);
    line([0 0],[min(ylim) max(ylim)],'color','r','linestyle',':')
    xlabel('Time (ms)')
    ylabel('Counts')
    xlim([-10 10])
    set(gca,'color','none','box','off')
    
    subplot(223);
    hold on;
    bar(delay, q_ab, 'k', 'BarWidth', 1);
    plot(delay, q_ab_smooth, 'b','linewidth',1.5);
    plot(delay, q_ab_smooth+conf_limit_smooth, 'b','linewidth',1,'linestyle',':');
    plot(delay, q_ab_smooth-conf_limit_smooth, 'b','linewidth',1,'linestyle',':');
    xlabel('Time (ms)')
    ylabel('Cross Covariance')
    ylm = ylim;
    set(gca,'ylim',ylm);
    line([0 0],[min(ylim) max(ylim)],'color','r','linestyle',':')
    line([min(delay) max(delay)], [0 0], 'color','k');
    line([min(delay) max(delay)], [conf_limit conf_limit], 'color','r');
    line([min(delay) max(delay)], [-conf_limit -conf_limit], 'color','r');
    set(gca,'color','none','box','off')
    text(.05,.95,'Null interval','units','normalized','color','r')
    text(.05,.9,'Slow Co-Mod','units','normalized','color','b')

    
    subplot(224);
    hold on;
    bar(delay, q_ab, 'k', 'BarWidth', 1);
    plot(delay, q_ab_smooth, 'b','linewidth',1.5);
    plot(delay, q_ab_smooth+conf_limit_smooth, 'b','linewidth',1,'linestyle',':');
    plot(delay, q_ab_smooth-conf_limit_smooth, 'b','linewidth',1,'linestyle',':');
    xlabel('Time (ms)')
    ylabel('Cross Covariance')
    set(gca,'ylim',ylm); 
    line([0 0],[min(ylim) max(ylim)],'color','r','linestyle',':')
    line([min(delay) max(delay)], [0 0], 'color','k');
    line([min(delay) max(delay)], [conf_limit conf_limit], 'color','r');
    line([min(delay) max(delay)], [-conf_limit -conf_limit], 'color','r');
    set(gca,'color','none','box','off')
    
    lwp = 1.5; 
    if ccg_data.ab_comod_adj_trough_sig
        line([ccg_data.ab_comod_adj_trough_delay ccg_data.ab_comod_adj_trough_delay],[min(ylim) max(ylim)],'color','g','linewidth',lwp)
    end
    if ccg_data.ab_comod_adj_peak_sig
        line([ccg_data.ab_comod_adj_peak_delay ccg_data.ab_comod_adj_peak_delay],[min(ylim) max(ylim)],'color','m','linewidth',lwp)
    end
    if ccg_data.ba_comod_adj_trough_sig
        line([-ccg_data.ba_comod_adj_trough_delay -ccg_data.ba_comod_adj_trough_delay],[min(ylim) max(ylim)],'color','g','linewidth',lwp)
    end
    if ccg_data.ba_comod_adj_peak_sig
        line([-ccg_data.ba_comod_adj_peak_delay -ccg_data.ba_comod_adj_peak_delay],[min(ylim) max(ylim)],'color','m','linewidth',lwp)
    end
    if ccg_data.ab_comod_adj_peak_sig || ccg_data.ba_comod_adj_peak_sig
        text(.05,.95,'Peak','units','normalized','color','m')
    end
    if ccg_data.ab_comod_adj_trough_sig || ccg_data.ba_comod_adj_trough_sig
        text(.05,.9,'Trough','units','normalized','color','g')
    end
    xlim([-10 10])
    
    %%
    
end