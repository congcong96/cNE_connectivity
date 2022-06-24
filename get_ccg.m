function [CCG, taxis] = get_ccg(train1, train2, tw, binsize)
% old name: calculateCCG
% Input:
%   train1: binned spktrain
%   train2: binned spktrain
%   tw: time window to calculate CCG with
%   binsizeL bin size for train1 and train2
% tw and binsize in ms

%Based on (Coordinated Neuronal Activity Enhances Corticocortical Communication, Amin
%Zandvakili, 2015)

Tbinnum = length(train1);
T = Tbinnum*binsize/1000;
CCG = zeros(size(train2,1),tw/binsize*2+1);

ii = 1;
for tau = -tw:binsize:tw
    idxtau = tau/binsize;
    idxspk1tau = find(train1 > 0)+idxtau;
    idxspk1tau(idxspk1tau<=0 | idxspk1tau >=Tbinnum) = [];
    CCG(:,ii) = sum(train2(:,idxspk1tau),2)/(length(idxspk1tau)*binsize/1000);
    ii = ii + 1;
end
taxis = -tw:binsize:tw;