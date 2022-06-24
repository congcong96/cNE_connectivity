function plot_ccg(centers, ccg, baseline, xl, yl, efficacy, c, varargin)

p = inputParser;
addParameter(p, 'ymax', [])
parse(p,varargin{:});
ymax = p.Results.ymax;

if nargin == 6
    c = 'b';
end
bar(centers, ccg, 1, 'k')
hold on
if isempty(ymax)
    ymax = max([max(ccg)*1.1, 1]);
end
plot([0, 0], [0, ymax], [c,'--'], 'linewidth', 1)
plot(centers, baseline, c)
xlim([-40 40])
ylim([0, ymax])
xlabel(sprintf(xl))
ylabel(sprintf(yl))
text(10, 0.8*ymax, sprintf('%.2f', efficacy), 'FontSize', 8)
end