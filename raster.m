function r = raster(spiketimes, reftimes, hw)

% get raster of spktimes aligned with reftime
%
% INPUT
% ----------------
% spiketimes: vector
%       times pf spikes, in ms
% reftimes: vector
%       times of reference, in ms
% hw: double
%       half-window to get raster, in ms;
%       40 means get raster -40 to 40ms around reference time points
%
% OUTPUT
% ----------------
% r: matrix (number of reference point by 2)
%       raster data
%       the first column are spike times related to reference 
%       the second column are indices of reference 

r = cell(length(reftimes), 1);
for ii = 1:length(reftimes)
    
    r_tmp = spiketimes(spiketimes > reftimes(ii)-hw ...
        & spiketimes < reftimes(ii)+hw) - reftimes(ii);
    r{ii} = [r_tmp, ii*ones(size(r_tmp))];
end
r = cell2mat(r);