function eff = get_efficacy(ccg, baseline, time, interval, nspk)

[~, idx] = ismember(interval, time);
causal_spikes = ccg(idx) - baseline(idx);
eff = sum(causal_spikes(causal_spikes>0))/nspk;