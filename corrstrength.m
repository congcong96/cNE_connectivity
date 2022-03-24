function rho = corrstrength(delay, r12, n1, n2, stimdur)
% corrstrength  Connection strength from correlogram
% 
%    rho = corrstrength(delay, r12, n1, n2, stimdur)
% 
%    Inputs:
%    ------------------------------
%    delay : lag, in ms
%    r12 : correlogram
%    n1 : number of spikes in train 1
%    n2 : number of spikes in train 2
%    stimdur : duration of stimulus, in number of bins
% 
%    Outputs:
%    ------------------------------
%    rho : correlation coefficient
% 
%          rho = ( r12max - n1*n2/stimdur ) / ...
%             ( (n1-n1^2/stimdur) * (n2-n2^2/stimdur) )^(1/2);
%
%


% indices
indneg = find(delay<=-20);
indneg = max(indneg);
indpos = find(delay>=20);
indpos = min(indpos);

indcenter = indneg:indpos;

r12center = r12(indcenter);
delaycenter = delay(indcenter);


% Maximum value in correlogram
[temp, indmax] = max(r12center);
indpd = max(indmax);
pd = delaycenter(indpd);
r12max = r12center(indpd);

rho = ( r12max - n1*n2/stimdur ) / ( (n1-n1^2/stimdur) * (n2-n2^2/stimdur) )^(1/2);

return;






