clc; clear
close all
% Created by Ze Jiang on 22/07/2024

%% synthetic data generation
% initialize random number generator
iseed = 101;    % seed number for random number generator
rng(iseed,'twister')

N = 500;

% case 1
%X = randn(N,1); % random predictors

% case 2
fs = 50;
dt = 1/fs;
t = 0:dt:dt*(N-1); 
%X = (sin(2*pi*t)+0.1*randn(1,N))';
X = (sin(2*pi*t+randn(1,1))+0.1*randn(1,N))';
plot(t,X)

% wavelet with N vanishing moments 
n_vanish = 1; % db1 is equivalent to haar
% name of wavelet filter
switch 1
    case 1
        %Daubechies wavelet with N vanishing moments, where N is a positive integer from 1 to 45.
        wname = ['db' num2str(n_vanish)] 
    case 2
        %Symlets wavelet with N vanishing moments, where N is a positive integer from 2 to 45.
        wname = ['sym' num2str(n_vanish)] 
    case 3 
        %Coiflets wavelet with N vanishing moments, where N is a positive integer from 1 to 5.
        wname = ['coif' num2str(n_vanish)] 
end

% maximum decomposition level: floor(log2(size(X,1)))
lev = floor(log2(N))-1; 

%% scale to frequency
% Set sampling period DELTA (dt) and wavelet name.
delta = 1/12; 
% Define scales.
scales = 2.^(1:lev);

% Compute associated pseudo-frequencies.
% scal2frq: returns the pseudo-frequencies corresponding to the 
% scales  given by A and the wavelet function 'wname' 
% and the sampling period DELTA.
f = scal2frq(scales,wname,delta); % f = centfrq(wname)./(a.*dt); 

% Compute associated pseudo-periods.
per = 1./f

%% wavelet transforms
X_tmp = X(:,1); % choose any variable, predictor or response

disp('DWT MRA')
X_DWT_MRA = dwtmra(X_tmp, wname, lev);
disp(['Additive:' num2str(sum(abs(sum(X_DWT_MRA,2)-X_tmp)))])
disp(['Variance:' num2str(sum(var(X_DWT_MRA))-var(X_tmp))])

disp('MODWT')
X_MODWT = (modwt(X_tmp, wname, lev))'; 
disp(['Additive:' num2str(sum(abs(sum(X_MODWT,2)-X_tmp)))])
disp(['Variance:' num2str(sum(var(X_MODWT))-var(X_tmp))])

disp('MODWT MRA')
X_MODWT_MRA = (modwtmra(X_MODWT', wname))';
disp(['Additive:' num2str(sum(abs(sum(X_MODWT_MRA,2)-X_tmp)))])
disp(['Variance:' num2str(sum(var(X_MODWT_MRA))-var(X_tmp))])

disp('AT')
X_AT = AT(X_tmp, wname, lev); 
disp(['Additive:' num2str(sum(abs(sum(X_AT,2)-X_tmp)))])
disp(['Variance:' num2str(sum(var(X_AT))-var(X_tmp))])

% plot wavelet decompositions of different wavelet transforms
figure
subplot(lev+2,1,1)
plot(X_tmp, 'k');
hold on
sgtitle(['Wavelet: ' num2str(wname)])
xlim([0 N]); ylim([min(X_tmp)*1.1,max(X_tmp)*1.1]); 
ylabel('X');
for is=1:lev+1
    subplot(lev+2,1,is+1)
 
    plot(X_DWT_MRA(:,is),'g');
    xlim([0 N]); ylim([min(X_tmp)*1.1,max(X_tmp)*1.1]); 
    hold on
    plot(X_MODWT(:,is),'r');
    hold on
%     plot(X_MODWT_MRA(:,is),'r');
%     hold on
    plot(X_AT(:,is),'b');
    hold off

    if is==lev+1
        ylabel(['A',int2str(is-1)])
    else
        ylabel(['D',int2str(is)])
    end
    
    legend('DWT MRA','MODWT','AT','NumColumns',1,'location','eastoutside')      
end
