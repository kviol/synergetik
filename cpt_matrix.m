% phase transition of recurrence plots (CP-Analysis)
function [cp_mean_M, cp_var_M] = cpt_matrix(data,max,L_TS)
L=length(data(:,1)); % number of lines (assuming time is columns)
cpts_mean=NaN(L,max);
cpts_var=NaN(L,max);
for t=1:L
    [cp_mean,cp_var]=cpt_short(data(t,:),max);
    cpts_mean(t,:)=cp_mean;
    cpts_var(t,:)=cp_var;
    clear cp_mean cp_var
end
outl_mean=isoutlier(cpts_mean);  % define outliers
outl_var=isoutlier(cpts_var);
cpts_mean(outl_mean==1)=NaN;   % remove outliers
cpts_var(outl_var==1)=NaN;
for i=1:max
    cpt_mean(i,:)=stats(cpts_mean(:,i),length(data),L_TS); % keep significant cps only
    cpt_var(i,:)=stats(cpts_var(:,i),length(data),L_TS); % keep significant cps only

    cp_mean_M(i)=round(nanmean(cpt_mean(i,:)));
    cp_var_M(i)=round(nanmean(cpt_var(i,:)));
end
end