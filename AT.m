function  [X_AT] = AT(X, wname, lev)

% Created by Ze Jiang on 15/09/2021
% a trous wavelet transform

  s = nan(length(X),lev);
  for i = 1:lev
    tmp = (modwt(X, wname, i))' ;
    s(:,i) = tmp(:,i+1) ;
  end

  X_AT = nan(length(X),lev+1);
  X_AT(:,1) = X-s(:,1); 
  for i = 1:(lev-1) 
      X_AT(:,i+1) = s(:,i)-s(:,i+1) ;
  end
  X_AT(:,lev+1) =  s(:,lev); 

end
