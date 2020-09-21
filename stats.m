function [CPs,p,CI,x1_IQR_mean,x2_IQR] = stats(CPs,L,L_TS)
if nargin < 3
    L_TS=L;
end
%% Test if distribution is significantly different from random 
CPs_clean=CPs;
CPs_clean(sum(isnan(CPs_clean), 2) == 1, :) = [];  % delete NaNs

% TS with 1 for cp and 0 for none
CPs_TS=zeros(L,1);
for j=1:length(CPs_clean)
    if ~isnan(CPs_clean(j))
        CPs_TS(CPs_clean(j))=1;
    end
end

% compare IQRs
rng('default')
N=100; % number of random tests
numCPs=length(CPs_clean);

x_rand=unidrnd(L,N,numCPs);   % creates N random data from equal distribution

x2_IQR=iqr(CPs_clean);% IQR für real data berechnen

for n=1:N % für alle Zufallsziehungen
    x1_IQR(n)=iqr(x_rand(n,:));
end

CI=CIs(x1_IQR,0.05,2);
x1_IQR_mean=mean(x1_IQR);

if x2_IQR < CI(1)
    p="significant";
    
else
    p="not significant";
    CPs=NaN(length(CPs),1);
end
%% remove all mean CPs at the extremes
% if the mean is within the extremes, set all to NaN
if nanmean(CPs)<=0.1*L_TS || nanmean(CPs)>=0.9*L_TS
    CPs=NaN(length(CPs),1);
end

end

%% function CI (confidence intervals)
function [yCI] = CIs(x,a,tail)
if nargin<3
    tail=2;
end
if nargin<2
    tail=2;
    a=0.05;
end
if size(x,2)>size(x,1) % convert Spalten to Spaltenvektor
    x=x';
end
% x: Variable to be assessed
% a: p-value
% tail: 1 for one-sided, 2 for two-sided

% NaNs entfernen
x_clean=x(any(~isnan(x),2));

b=1-a/tail;
N = length(x_clean);
ME=mean(x_clean);
SEM = std(x_clean)/sqrt(N);
T = tinv(b, N-1);
CI=T*SEM;
yCI = [ME-CI, ME+CI];
end
