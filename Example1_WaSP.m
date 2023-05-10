clc
clear
close all

% Created by Ze Jiang on 15/09/2021 (ze.jiang@unsw.edu.au)
% Y: response = N x 1
% X: predictor= (N+N_fc) x n_var

N = 128;         % number of observation
N_fc=0;         % number of forecast (optional)
n_var=4;        % number of variable
iseed = 101;    % seed number for random number generator

%% initialize random number generator
rng(iseed,'twister')

%% synthetic data generation
t = linspace(-pi,pi,N);
Y = (sin(t) + 1.0*randn(size(t)))' %+ t'; % sine wave add noise

X = randn(N+N_fc,n_var) %- repmat(linspace(-pi,pi,N+N_fc)',1, n_var); % random predictors
% X = [randn(N,n_var); repmat(100,N_fc,n_var)]; % test on abnormal value from forecast

%% Daubechies wavelet with N vanishing moments 
% N is a positive integer from 1 to 45
n_vanish = 2; 
wname = ['db' num2str(n_vanish)] % db1 is equivalent to haar
flag_sign=0;


%% method of discrete wavelet transform 
switch 1
    case 1
        method = 'dwtmra' %dwtmra cannot be used in the forecast setting
    case 2
        method = 'modwt'
    case 3 
        method = 'at'
end

% maximum decomposition level: floor(log2(size(X,1)))
% or rule of thumb decomposition level: ceiling(log(n/(2*v-1))/log(2))-1 (Kaiser, 1994)
lev = floor(log2(size(X,1)))-1 ; 

% variance transformation using WaSP
[X_WaSP, C] = WaSP(Y, X, method, wname, lev, flag_sign); 

%% linear regression for each variable
RMSE=nan(1,n_var);
RMSE_WaSP=nan(1,n_var);
RMSE_opti=nan(1,n_var);
for i_var = 1:n_var
    % optimal RMSE - Eq. 12 in WRR2020 paper
    ratio=var(X(1:N,i_var))/var(X_WaSP(1:N,i_var));
    %disp(ratio)
    RMSE_opti(i_var) = sqrt((N-1)/N*(var(Y)-(norm(C(:,i_var))^2)*ratio)); 
    
    % Std model: linear regression using original predictor
    p_coeff = polyfit(X(1:N,i_var), Y, 1)  ;
    Y_fit = polyval(p_coeff, X(1:N, i_var)) ; 
    RMSE(i_var) = sqrt(mean((Y-Y_fit).^2));
    
    % VT model: linear regression using transformed predictor
    p_coeff = polyfit(X_WaSP(1:N,i_var), Y, 1)  ;
    Y_fit = polyval(p_coeff, X_WaSP(1:N, i_var)) ; 
    RMSE_WaSP(i_var) = sqrt(mean((Y-Y_fit).^2));

end

%% plot RMSE
figure
bar([RMSE',RMSE_WaSP',RMSE_opti']);
xlabel('No. of variable'); ylabel('RMSE');
legend('Std','VT','Optimal','NumColumns',1,'location','eastoutside')   
title(['Variance transformation based on ',num2str(method) ' using ' num2str(wname)])
saveas(gca,'RMSE.fig');

%% plot - Y, X and X_WaSP
figure
sgtitle(['Calibration: ',num2str(method) ' using ' num2str(wname)])
for i_var = 1:n_var
    subplot(n_var,1,i_var)
    
    plot(Y, 'k')     
    xlim([0,N+N_fc]); ylim([min(Y)*1.1,max(Y)*1.1]); 
    hold on 
    plot(X(:,i_var), 'b')
    hold on
    plot(X_WaSP(:,i_var), 'r')     
    hold off
    legend('Y','X','Xnew','NumColumns',1,'location','eastoutside')      
end 
saveas(gca,'comparision.fig');


%% validation
X_val=randn(N,n_var); % random predictors

% variance transformation using derived C from calibration period
X_WaSP_val = WaSP_val(X_val, C, method, wname,flag_sign); 

figure
sgtitle(['Validation: ',num2str(method) ' using ' num2str(wname)])
for i_var = 1:n_var
    subplot(n_var,1,i_var)
    
    plot(X_val(:,i_var), 'b')
    hold on
    plot(X_WaSP_val(:,i_var), 'r')     
    hold off
    legend('X','Xnew','NumColumns',1,'location','eastoutside')      
end 
saveas(gca,'comparision_val.fig');
