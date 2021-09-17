
function [X_WaSP, C] = WaSP(Y, X, method, wname, lev)

% Created by Ze Jiang on 13/09/2021 (ze.jiang@unsw.edu.au)
% Note: This is a simplified version designed for short-term forecasts transformation. 

% Y: response = m x 1
% X: predictor= (m+l) x n;  m: no. of obs, n: no. of vars, l: no. of forecasts (optional)
% method: discrete wavelet transform, including dwtmar, modwt, modwtmar, and at
% wname: wavelet filter, Daubechies wavelets are widely used. 
% lev: decomposition levels

% OUTPUT:
% X_WaSP: variance transformed X = (m+l) x n
% C: covariance vector for each predictor variable C = (lev+1) x n

% USAGE:
% [X_WaSP, C] = WaSP(Y, X, method, wname, lev)

% REFERENCE:
% Jiang, Z., Sharma, A., & Johnson, F. (2020). Refining Predictor Spectral Representation Using Wavelet Theory for Improved Natural System Modeling. Water Resources Research, 56(3). https://doi.org/10.1029/2019WR026962
% Jiang, Z., Rashid, M. M., Johnson, F., & Sharma, A. (2020). A wavelet-based tool to modulate variance in predictors: An application to predicting drought anomalies. Environmental Modelling & Software, 135, 104907. https://doi.org/10.1016/j.envsoft.2020.104907
% Jiang, Z., Sharma, A., & Johnson, F. (2021). Variable transformations in the spectral domain - Implications for hydrologic forecasting. Journal of Hydrology, 603, 126816. https://doi.org/10.1016/j.jhydrol.2021.126816

% dimension of predictor
[num_obs, num_var] = size(X) ;
% dimension of response 
N = length(Y);

% output matrix
C = nan(lev+1, num_var); 
X_WaSP = nan(num_obs,num_var) ;

% variance transformation for each variable
for i_var = 1 : num_var
    % center
    X_tmp = X(:,i_var)-mean(X(:,i_var)) ; 

    % wavelet transform
    if isequal(method,'modwt')
        X_WT = (modwt(X_tmp, wname, lev))'; 
    elseif isequal(method,'dwtmra')
        X_WT = dwtmra(X_tmp, wname, lev); 
    elseif isequal(method,'modwtmra')
        X_WT_tmp = modwt(X_tmp, wname, lev); 
        X_WT = (modwtmra(X_WT_tmp, wname))'; 
    elseif isequal(method,'at')
        X_WT = AT(X_tmp, wname, lev); 
    else
        disp('This type of wavelet transform does not apply!')
    end
    
%     disp(['Additive:' num2str(sum(abs(sum(X_WT,2)-X_tmp)))])
%     disp(['Variance:' num2str(sum(var(X_WT))-var(X_tmp))])
 
    %standardization w.r.t. observed period
    X_WT_n = X_WT(1:N,:);
    %X_WT_norm=(X_WT-mean(X_WT_n))./std(X_WT_n);
    
    % alternative way to standardize data for older matlab versions
    %subtract mean of each column
    X_WT_c=X_WT-repmat(mean(X_WT_n),num_obs,1);
    %divide by the standard deviation of each column
    X_WT_norm=X_WT_c./repmat(std(X_WT_n,0,1),num_obs,1);
%     disp(sum(abs(X_WT_norm1-X_WT_norm)))   
    
%     disp(var(X_WT_norm(1:N,:)))  
%     disp(mean(X_WT_norm(1:N,:))) 
%     disp(sum(abs(normalize(X_WT_n)-X_WT_norm(1:N,:))))

    %corr1 = corrcoef([Y X_WT(1:length(Y),:)]) ; %correct
    %corr1 = cov([Y X_WT(1:length(Y),:)]) ; %wrong
    %C(:,i_var) = corr1(1,2:lev+2); 
    
    % covariance - Eq. 10 in WRR2020 paper
    corr = 1/(N-1)*Y'*X_WT_norm(1:N,:); 
    %corr = 1/(length(Y)-1)*Y'*normalize(X_WT(1:length(Y),:)); 
    %disp(normalize(corr1(1,2:lev+2),'norm') - normalize(corr,'norm'))
    C(:,i_var) = corr; 
    
    % normalization to get unit norm
    %C_norm = normalize(C(:,i_var),'norm') ;
    C_norm = C(:,i_var)./sqrt(sum(C(:,i_var).^2)) ;
%     disp(norm(C_norm))
%     disp(norm(C(:,i_var))-sqrt(sum(C(:,i_var).^2)))
    
    % variance transformation -  Eq. 9 in WRR2020 paper
    %X_WaSP(:,i_var) = normalize(X_WT)*(std(X(:,i_var)).*C_norm(:)) ; 
    X_WaSP(:,i_var) = X_WT_norm*(std(X(1:N,i_var)).*C_norm) ; 
    
    % add mean back
    X_WaSP(:,i_var) = X_WaSP(:,i_var) + mean(X(1:N,i_var)); 
   
end

