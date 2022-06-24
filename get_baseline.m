function [baseline, varargout] = get_baseline(spiketimes1, spiketimes2, binsize)

if nargin == 2
    binsize = 0.5;
end
maxlag = 50;
maxtime = max([spiketimes1(:); spiketimes2(:)]);
spktrain1 = histcounts(spiketimes1, 0:binsize:maxtime);
spktrain2 = histcounts(spiketimes2, 0:binsize:maxtime);

sigma_ms = 7; 
hsize = [(sigma_ms^2)/binsize 1];
sigma = sigma_ms/binsize; 
h = fspecial('gaussian',hsize,sigma);

spktrain1_smooth = conv( spktrain1, h, 'same' );
spktrain2_smooth = conv( spktrain2, h, 'same' );
baseline  = xcorr( spktrain1_smooth, spktrain2_smooth, maxlag/binsize );

varargout{1} = xcorr( spktrain1, spktrain2, maxlag/binsize );
