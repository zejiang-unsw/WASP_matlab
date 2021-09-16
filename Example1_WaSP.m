clc
clear
close all

% Created by Ze Jiang on 15/09/2021
% Y: response = N x 1
% X: predictor= (N+N_fc) x n_var

N = 400;        % number of observation
N_fc=0;         % number of forecast (optimal)
n_var=4;        % number of variable
iseed = 101;    % seed number for random number generator

% initialize random number generator
rng(iseed,'twister')

% synthetic data generation
t = linspace(-pi,pi,N);
Y = (sin(t) + 1.0*randn(size(t)))'; % sine wave add noise
X = randn(N+N_fc,n_var); % random predictors

% Daubechies wavelet with N vanishing moments 
% N is a positive integer from 1 to 45
n_vanish = 1; 
wname = ['db' num2str(n_vanish)] % db1 is equivalent to haar

% method of discrete wavelet transform 
switch 2
    case 1
        method = 'dwtmra'
    case 2
        method = 'modwt'
    case 3 
        method = 'at'
end

% maximum decomposition level: floor(log2(size(X,1)))
lev = floor(log2(size(X,1)))-1 ; 

% variance transformation using WaSP
[X_WaSP, C] = WaSP(Y, X, method, wname, lev); 

% linear regression for each variable
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

% plot RMSE
figure
bar([RMSE',RMSE_WaSP',RMSE_opti']);
xlabel('No. of variable'); ylabel('RMSE')
legend ('Std','VT','Optimal','NumColumns',1,'location','eastoutside')   

% plot - Y, X and X_WaSP
figure
sgtitle(['Variance transformation based on ',num2str(method) ' using ' num2str(wname)])
for i_var = 1:n_var
    subplot(n_var,1,i_var)
    
    plot(Y, 'k')     
    xlim([0,N+N_fc]); ylim([min(Y)*1.1,max(Y)*1.1]); 
    hold on 
    plot(X(:,i_var), 'b')
    hold on
    plot(X_WaSP(:,i_var), 'r')     
    hold off
    legend ('Y','X','Xnew','NumColumns',1,'location','eastoutside')      
end 

