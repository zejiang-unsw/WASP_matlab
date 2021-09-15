clc
clear
%close all

N = 40; 
N_fc=0;
n_var=4; 
iseed = 101; 
n_vanish = 1;

switch 2
    case 1
        method = 'dwtmra'
    case 2
        method = 'modwt'
    case 3 
        method = 'at'
end


switch 1
    case 1
        wname = ['db' num2str(n_vanish)]
    case 2
        wname = ['sym' num2str(n_vanish)]
    case 3 
        wname = ['coif' num2str(n_vanish)]
end

rand('seed',iseed);
Md1 = arima('Constant',0,'AR',{0.7 0.25},'Variance',0.5);
Y = simulate(Md1,N);

randn('seed',iseed);
X = randn(N+N_fc,n_var);

% maximum level floor(log2(size(X,1)))
lev = floor(log2(size(X,1))) ; 

% variance transformation using WaSP
[X_WaSP, C] = WaSP(Y, X, method, wname, lev); 

% cov([Y, X(1:N,:)])
% corrcoef([Y, X(1:N,:)])
% linear regression
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

RMSE_WaSP-RMSE_opti
% RMSE-RMSE_opti
% RMSE-RMSE_WaSP

% plot - Y, X and X_WaSP
figure
sgtitle([num2str(method) ': ' num2str(wname)])
for i_var = 1:n_var
    subplot(n_var,1,i_var)
    
    plot(Y, 'k')     
    ylim([-3,3])
    hold on 
    plot(X(:,i_var), 'b')
    hold on
    plot(X_WaSP(:,i_var), 'r')     
    hold off
    legend ('Y','X','Xnew','NumColumns',1,'location','eastoutside')      
end 

