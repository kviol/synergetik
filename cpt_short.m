function [row_mean,row_var]=cpt_short(A,max_change)
% CPs with respect to the shifts of the mean
TF_mean=ischange(A,'mean','MaxNumChanges',max_change); %1: entlang der Spalten
row_mean=NaN(1,max_change);
CPs_mean=find(TF_mean==1);  % gibt Indizes alle Zeilen und Spalten mit 1
for i=1:length(CPs_mean)
    row_mean(i)=CPs_mean(i);
end

% CPs with respect to the shifts of the variance
TF_var=ischange(A,'variance','MaxNumChanges',max_change); %1: entlang der Spalten
row_var=NaN(1,max_change); % make sure that NaNs are stated when no change point or less than the max are found
CPs_var=find(TF_var==1);  % gibt Indizes alle Zeilen und Spalten mit 1
for i=1:length(CPs_var)
    row_var(i)=CPs_var(i);
end
end % function


