
function [X_WaSP, lev] = WaSP(Y, X, method, wname, lev)

% Created by Ze Jiang on 13/09/2021
% Y: response = m x 1
% X: predictor= (m+l) x n
% OUTPUT:
% X_WaSP = (m+l) x n; m: no of obs, n: no of vars, l: no of forecasts (optional)
% USAGE:
% [X_WaSP] = WaSP(Y, X)

% method = 'modwt'
% wname = 'haar'
% Y = rain_YJ_train; 
% X = pred_YJ_train_left; 

[num_obs, num_var] = size(X) ;

% wavelet transform
if  isempty(lev) 
    lev = floor(log2(size(X,1))) ; 
end

X_WT = nan(num_obs,lev+1) ;
X_WaSP = nan(num_obs,num_var) ;

% disp(lev)

for var = 1 : num_var
    X_tmp = X(:,var)-mean(X(:,var)) ; 

    if isequal(method,'modwt')
        X_WT = (modwt(X_tmp, lev, wname))'; 
    else isequal(method,'dwt')
        X_WT = (dwt(X_tmp, lev, wname))'; 
    end
    
    corr = corrcoef([Y X_WT(1:length(Y),:)]) ;
    cov = normalize(corr(1,2:lev+2),'norm') ;
    
    X_WaSP(:,var) = normalize(X_WT)*(std(X(:,var)).*cov(:)) ; 
    
    X_WaSP(:,var) = X_WaSP(:,var) + mean(X(:,var)); 
end

% for var = 2
%     plot(pred_YJ_train_left(:,var), 'b')
%     hold on
%     plot(X_WaSP(:,var), 'r')     
%     hold on 
%     plot(rain_YJ_train-mean(rain_YJ_train), 'k')         
% end 
            
