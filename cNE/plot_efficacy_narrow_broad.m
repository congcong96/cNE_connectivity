narrow_spiking_criteria = 0.5;
cmap = cbrewer('qual', 'Dark2', 8);

%%
spkfolder = '/data/congcong/rat_MGB_A1/3_singleunit/dmr';
figfolder = '/data/congcong/rat_MGB_A1/figure/singleunit/connect';
cd(spkfolder)
connect_files = dir(fullfile(spkfolder, '*pairs.mat'));
%%
efficacy = cell(1, length(connect_files));
efficacy_ne =  cell(1, length(connect_files));
cell_type =  cell(1, length(connect_files)); % 1 for narrow spiking
ne_label = cell(1, length(connect_files));

for ii = 1:length(connect_files)
    v = whos('-file', connect_files(ii).name);
    if ismember('cNE_connected', {v.name})
        load(connect_files(ii).name, 'cNE_connected');
        spkfile = dir(sprintf('%s*H22x32*.mat',connect_files(ii).name(1:13) ));
        load(spkfile.name, 'waveform')
        n_pair = length(cNE_connected);
        pair = cell(0, 1);
        for jj = 1:n_pair
            pair_tmp = [cNE_connected(jj).cNE, ...
                cNE_connected(jj).idx_A1];
            if sum(cellfun(@(x) isequal(x, pair_tmp), pair)) == 0
                pair{jj} = pair_tmp;
                
                % copy efficacy of the whole spike train or ne train
                efficacy_ne{ii}{jj} = [cNE_connected(jj).ccg_efficacy_ne_spon];
                efficacy{ii}{jj} = [cNE_connected(jj).ccg_efficacy_spon];
                % file , NE and cortical neuron label
                ne_label{ii}{jj} = repmat([ii, pair_tmp]', size(efficacy{ii}{jj}));
                % label cell type of cortical neuron
                if waveform(pair_tmp(2)).tpd < narrow_spiking_criteria
                    cell_type{ii}{jj} = ones(size(efficacy{ii}{jj}));
                else
                    cell_type{ii}{jj} = zeros(size(efficacy{ii}{jj}));
                end
            end
        end
    end
end

%%
efficacy = cellfun(@(x) cell2mat(x), efficacy,  'UniformOutput', false);
efficacy_ne = cellfun(@(x) cell2mat(x), efficacy_ne,  'UniformOutput', false);
cell_type = cellfun(@(x) cell2mat(x), cell_type,  'UniformOutput', false);
ne_label = cellfun(@(x) cell2mat(x), ne_label,  'UniformOutput', false);
efficacy = cell2mat(efficacy);
efficacy_ne = cell2mat(efficacy_ne);
cell_type = cell2mat(cell_type);
ne_label = cell2mat(ne_label);
idx = efficacy > .02;% only caculate gain for pairs with efficacy > 0.02
efficacy = efficacy(idx);
efficacy_ne = efficacy_ne(idx);
cell_type = cell_type(idx);
ne_label = ne_label(:, idx);

gain = efficacy_ne./efficacy;

figure
bar(1, mean(gain(cell_type==1)),'FaceColor', cmap(1,:))
hold on
bar(2, mean(gain(cell_type==0)),'FaceColor', cmap(2,:))
xticks([1 2])
xticklabels({'narrow', 'broad'})
errorbar(1, mean(gain(cell_type==1)), std(gain(cell_type==1)), 'k')
errorbar(2, mean(gain(cell_type==0)), std(gain(cell_type==0)), 'k')
ylabel('efficacy gain')
saveas(gcf, 'narrow_vs_broad_efficacy_gain.jpg')
saveas(gcf, 'narrow_vs_broad_efficacy_gain.fig')
printPDFandPSC(gcf, 'narrow_vs_broad_efficacy_gain');

figure
bar(1, mean(efficacy(cell_type==1)),'FaceColor', cmap(1,:))
hold on
bar(2, mean(efficacy(cell_type==0)),'FaceColor', cmap(2,:))
xticks([1 2])
xticklabels({'narrow', 'broad'})
errorbar(1, mean(efficacy(cell_type==1)), std(efficacy(cell_type==1)), 'k')
errorbar(2, mean(efficacy(cell_type==0)), std(efficacy(cell_type==0)), 'k')
ylabel('efficacy')
saveas(gcf, 'narrow_vs_broad_efficacy.jpg')
saveas(gcf, 'narrow_vs_broad_efficacy.fig')
printPDFandPSC(gcf, 'narrow_vs_broad_efficacy');

figure
bar(1, mean(efficacy_ne(cell_type==1)),'FaceColor', cmap(1,:))
hold on
bar(2, mean(efficacy_ne(cell_type==0)),'FaceColor', cmap(2,:))
xticks([1 2])
xticklabels({'narrow', 'broad'})
errorbar(1, mean(efficacy_ne(cell_type==1)), ...
    std(efficacy_ne(cell_type==1)), 'k')
errorbar(2, mean(efficacy_ne(cell_type==0)), ...
    std(efficacy_ne(cell_type==0)), 'k')
ylabel('efficacy')
saveas(gcf, 'narrow_vs_broad_efficacy_ne.jpg')
saveas(gcf, 'narrow_vs_broad_efficacy_ne.fig')
printPDFandPSC(gcf, 'narrow_vs_broad_efficacy_ne');

filenumbers = unique(ne_label(1,:));
n_A1_narrow = 0;
n_A1_broad = 0;
for ii = 1:length(filenumbers)
    filenumber = filenumbers(ii);
    A1_numbers = ne_label(3,...
        ne_label(1,:)==filenumber ...
        & cell_type == 1);
    n_A1_narrow = n_A1_narrow + length(unique(A1_numbers));
    A1_numbers = ne_label(3,...
        ne_label(1,:)==filenumber ...
        & cell_type == 0);
    n_A1_broad = n_A1_broad + length(unique(A1_numbers));
end
fprintf('Number of A1 neuron (narrow spiking): %d\n', n_A1_narrow)
fprintf('Number of A1 neuron (broad spiking): %d\n', n_A1_broad)
fprintf('Number of cNE member - A1 neuron pairs (narrow spiking): %d\n', ...
    sum(cell_type==1))
fprintf('Number of cNE member - A1 neuron pairs (broad spiking): %d\n', ...
    sum(cell_type==0))