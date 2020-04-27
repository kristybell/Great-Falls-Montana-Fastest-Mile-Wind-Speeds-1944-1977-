% At Great Falls, Montana, the largest yearly fastest-mile wind speeds at
% 10m above ground during the period 1944-1977 (sample size n = 34) were:
% 57, 65, 62, 58, 64, 65, 59, 65, 59, 60, 64, 65, 73, 60, 67, 50, 74, 60,
% 66, 55, 51, 60, 55, 60, 51, 51, 62, 51, 54, 52, 59, 56, 52, 49
% (mph). The sample mean and the sample standard deviation for these data
% are X = 59 mph and s = 6.41 mph. For N = 50 years and N= 1000years,
% velocity(50) ~ 76mph    SD(velocity(50)) ~ 3.7mph
% velocity(100) ~ 91mph    SD(velocity(100)) ~ 6.4mph

clear
close all
clc

% INPUT data (variable: U_fmile_data)
% Fastest-mile wind speed data
U_fmile_data=[57;65;62;58;64;65;59;65;59;60;64;65;73;60;67;50;...
 74;60;66;55;51;60;55;60;51;51;62;51;54;52;59;56;52;49];

%% DATA ANALYSIS AND GUMBEL FIT


% Plotting position
Ufm_sorted=sort(U_fmile_data);      % Sorting the data in ascending order
N=length(Ufm_sorted);               % Sample size

Iv=[1:1:N]/(N+1);                   % Plotting position

% (Linear) Regression on Gumbel 'chart'
Y=-log(-log(Iv));                   % Reduced variate
Pev=polyfit(Y',Ufm_sorted,1)

mu=Pev(2);                          % Location
sigma=Pev(1);                       % Dispersion

figure(1)
plot(Y,Ufm_sorted,'bd',Y,Y*sigma+mu,'r--')
title('Reduced variate vs. velocity (mph)')
xlabel('Reduced variate')

%% PLOTTING THE THEORETICAL CDF FUNCTION & COMPARE RESULTS AGAINST DATA

% The 'evcdf' and 'evpdf' standard functions of Matlab are based on
% EV-1 Gumbel of MINIMA. In order to calculate the
% Gumbel CDF of the MAXIMA see information below
% y = evpdf(-x, -1*location parameter, scale parameter)
% Fy = 1-evcdf(-x,-1*location parameter, scale parameter)

Fv=1-evcdf(-Ufm_sorted,-mu,sigma);

figure(2)
plot(Ufm_sorted,Iv,'k.',Ufm_sorted,Fv)
title('CDF of velocity annual maxima')
xlabel('Annual velocity maxima, mph')
legend('Data','Gumbel CDF','location','northwest')


%% ESTIMATE WIND SPEEDS AT VARIOUS RECURRENCE INTERVALS
N_bar=[50,100,500,1000];            % Recurrence interval (years)
S_hat=std(U_fmile_data);            % Standard deviation of the sample
n=N;                                % Sample size

% Standard deviation of the wind speed estimator
SDVN=0.78*[1.64+1.46.*log(N_bar)-0.577+1.1.*...
 (log(N_bar)-0.577).^2].^.5*S_hat/sqrt(n);

% Find 95% confidence-interval wind speeds at various recurrence intervals
p=1-1./N_bar;                       % probability of NON-exceedance
VNe=-log(-log(p))*sigma+mu;         % estimate the VNe value via Gumbel CDF
VN_95p=VNe+3*SDVN;                  % adjust VNe using SDVN (confid. int.)
