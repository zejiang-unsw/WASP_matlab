clc
clear
close all

% Created by Ze Jiang on 15/09/2021
% Y: response = N x 1
% X: predictor= (N+N_fc) x n_var

N = 512;        % number of observation
N_fc=0;         % number of forecast (optimal)
n_var=4;        % number of variable
iseed = 101;    % seed number for random number generator

% initialize random number generator
rng(iseed,'twister')

% synthetic data generation
t = linspace(-pi,pi,N);
Y = (sin(t) + 1.0*randn(size(t)))'; % sine wave add noise
X = randn(N+N_fc,n_var); % random predictors

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

% wavelet transforms
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
