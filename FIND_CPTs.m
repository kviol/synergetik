function [RES]=FIND_CPTs(TimeSer,max_CP,output)
% input: vector or matrix (important: length of the time series has to be
% larger than the number of time series)
% output: [RES] structure with results

if nargin<3  % no output (default)
    output = 0;
end

if size(TimeSer,2)>size(TimeSer,1)  % make sure that lines = time points 
    TimeSer=TimeSer';
end

L = size(TimeSer,1);
Cols = size(TimeSer,2);

for c=1:Cols
    TS=TimeSer(:,c);
    %% on raw time series
    [cp_mean,cp_var]=cpt_short(TS,max_CP);
    CP.raw.mean=cp_mean;
    CP.raw.var=cp_var;
    Data.raw=TS;
    
    
    %% on recurrence plots
    dim=3; tau=1; % can be changed by user
    y = phasespace(TS,dim,tau); % by Hui Yang
    recplot = cerecurr_y(y); % by Hui Yang
    [cp_mean,cp_var]=cpt_matrix(recplot,max_CP,L);
    CP.rp.mean=cp_mean;
    CP.rp.var=cp_var;
    Data.rp=recplot;
    
    %% on Time-Frequency-Distributions
    [st,t,f]=TFD_fun(TS); % by Glen Stockwell
    z=imcomplement(abs(st));
    [cp_mean,cp_var]=cpt_matrix(z,max_CP,L);
    CP.tfd.mean=cp_mean;
    CP.tfd.var=cp_var;
    Data.tfd.z=st;
    Data.tfd.f=f;
    Data.tfd.t=t;
    
    %% on dynamic complexity
    win=7;
    DyK=DynK(TS,win,max(TS)-min(TS)); % by Günter Schiepek
    [cp_mean,cp_var]=cpt_short(DyK,max_CP);
    CP.dc.mean=cp_mean;
    CP.dc.var=cp_var;
    Data.dc=DyK; 
   
    %% remove outliers
    RES(c).CP=CP; % save individual CPs to structure
    RES(c).Data=Data;   
    
    CPs=[CP.raw.mean; CP.raw.var; CP.rp.mean; CP.rp.var; CP.tfd.mean; CP.tfd.var; CP.dc.mean; CP.dc.var];
    
    if max_CP<2
        outl=isoutlier(CPs);
        CPs(outl==1)=NaN;
    end
    
    %% cross-validate
    %set NaN if no result from RAW or RP (to reduce false-positives)
    if isnan(CPs(1)) && isnan(CPs(2)) && isnan(CPs(3)) && isnan(CPs(4))
        RES(c).CP.tfd.mean="NaN";
        RES(c).CP.tfd.var="NaN";
        RES(c).CP.dc.mean="NaN";
        RES(c).CP.dc.var="NaN";
        CPs(1:end)="NaN";
    end
    
    %% TP (cluster)
    if max_CP>1
        x=CPs(:);
        [idx,TP]=kmeans(x,max_CP);
        RES(c).TP=TP;        
        for i=1:max_CP
            RES(c).cluster(i).CPs=x(idx==i);
        end
        clear x
    else
        RES(c).TP=nanmean(CPs);
        RES(c).cluster(1).CPs=CPs;
    end
 
    %% stats
    for i=1:max_CP
        y=RES(c).cluster(i).CPs;
        [x,p,CI,x1_IQR_mean,x2_IQR] = stats(y,length(TS));
        CPs=x;
        clear x y
 
        RES(c).cluster(i).CPs=CPs;
        RES(c).Stats(i).p=p;
        RES(c).Stats(i).CI=CI;
        RES(c).Stats(i).x_rand=x1_IQR_mean;
        RES(c).Stats(i).x=x2_IQR;
    end    
    
    
    %% output
    if output
        figure
        plot(TS,'k')
        hold on
        for i=1:max_CP % mean for number of change points
            if isnan(nanmean(CPs(:,i)))==0
                xline(nanmean(CPs(:,i)),'LineWidth',5);
            end
        end
        for k=1:length(CPs(:,1)) % for several methods
            for i=1:max_CP
                if isnan(CPs(k,i))==0
                    plot(CPs(k,i),TS(CPs(k,i)),'.b','MarkerSize',18)
                    hold on
                end
            end
        end
        clear k i
            
        % stats output
        for i=1:max_CP
            if RES(c).Stats(i).p=="significant"
                dec=") is ";
                st_out2="The change point is significant.";
            else
                dec=") is NOT ";
                st_out2="The change point is NOT significant.";
            end
            disp(strcat('Change point #',num2str(i),':'))
            if isnan(nanmean(CPs(:,i)))==0
                disp(strcat('Change Point found at: ',num2str(nanmean(CPs(:,i)))))
                st_out=strcat('The IQR of this change point (',num2str(RES(c).Stats(i).x),dec,'smaller than the lower 95% confidence interval value (',num2str(RES(c).Stats(i).CI(1)),') of random change points.');
                disp(st_out)
                disp(st_out2)
            else
                disp("No change point found")
            end
        end
    end
end
end
