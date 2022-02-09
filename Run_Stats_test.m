function [p_kw,tbl,stats, p_anova] = run_stats_tests(data, groups)
[p_kw,tbl,stats] = kruskalwallis(data,groups);
title('Kruskal-Wallis') 
ylabel('magnitude (mV)')

% Trying to connect median lines to see the trend better

lines = findobj(gcf, 'type', 'line', 'Tag','Median');
xMed = mean(vertcat(lines.XData),2);
yMed = vertcat(lines.YData);
hold on
% plot(xMed, yMed, 'r.-') %this one has small dots to mark median
plot(xMed, yMed, 'ro-') %this one has a circle to mark median
% plot(xMed, yMed, 'r+-') %this one has a plus sign to mark median


% ANOVA TEST
% p_anova = anova1(data,groups);
% title('ANOVA1') 
% ylabel('magnitude (mV)')
end 