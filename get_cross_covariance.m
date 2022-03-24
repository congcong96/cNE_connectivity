function [qab, conf_limit, stimdur] = get_cross_covariance( correlogram, na, nb, dt, stimdur)
% get_cross_covariance Cross-covariance function from cross-correlation function
% 
%    [qab, conf_limit, stimdur] = get_cross_covariance(correlogram, na, nb, dt, stimdur)
%    correlogram : cross-correlation function for two spike trains
%    na : number of spikes in train A
%    nb : number of spikes in train B
%    dt : bin size, in ms, of spike trains
%    stimdur : duration of stimulus, in ms
% 
%    qab : cross-covariance function
%    conf_limit : 99% confidence limits for cross-covariance function
%    stimdur : stimulus duration used for analysis, in ms

narginchk(4, 5)

if ( nargin == 4 )
   stimdur = 900000;
end

if ( dt < 0.05 )
   dt = dt * 1000;
end

pa = na / stimdur;
pb = nb / stimdur;

pab = correlogram / dt / stimdur;

% pabasymptote = pa*pb; 
% Use empirical mixing value. This will approximate pa*pb
pabasymptote = mean([pab(1:5) pab(end-5:end)]); 

qab = pab - pabasymptote;
% qab = pab - pa*pb; % theoretical calculation; we use experimental value
qabmean = 0;
qabstd = (pa * pb / dt / stimdur)^(1/2);
conf_limit = qabmean + 4 * qabstd;