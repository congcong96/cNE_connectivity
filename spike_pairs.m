function [first_spikes, second_spikes] = spike_pairs(spiketimes, isi, DT, spiketimes2)
% get time of spikes in spike pairs
%
% INPUT
% ----------------
% spiketimes: vector
%       times of spikes, in ms
% isi: double
%       isi of spike pairs, get all pairs with interval less than isi 
% DT: double
%       dead time before the first spike
% spiketimes2: vector
%       times of second spike train, in ms
%       for heterosynaptic interaction

%
% OUTPUT
% ----------------
% first_spikes: vector
%       times of first spikes, in ms
% second_spikes: vector
%
% reference:
% "Synaptic interactions between thalamic inputs to 
% simple cells in cat visual cortex" (Usrey et al., 2000)
%
% created by Congcong, 2022-04-07

% interaction of homosynaptic spikes
if nargin == 3
    d = diff(spiketimes);
    idx = find(d(1:end-1) > DT & d(2:end) <= isi);
    first_spikes = spiketimes(idx+1);
    second_spikes = spiketimes(idx+2);
end