clc
clear
close all

N = 512; 
n_var=4; 
iseed = 101; 
n_vanish = 10;

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
X = randn(N,n_var);

% maximum level floor(log2(size(X,1)))
lev = floor(log2(N)) ; 

% wavelet transform
X_tmp = X(:,1); 
X_DWT_MRA = dwtmra(X_tmp, wname, lev);
disp('DWT MRA')
disp(['Additive:' num2str(sum(abs(sum(X_DWT_MRA,2)-X_tmp)))])
disp(['Variance:' num2str(sum(var(X_DWT_MRA))-var(X_tmp))])

X_MODWT = (modwt(X_tmp, wname, lev))'; 
disp('MODWT')
disp(['Additive:' num2str(sum(abs(sum(X_MODWT,2)-X_tmp)))])
disp(['Variance:' num2str(sum(var(X_MODWT))-var(X_tmp))])

X_MODWT_MRA = (modwtmra(X_MODWT', wname))';
disp('MODWT MRA')
disp(['Additive:' num2str(sum(abs(sum(X_MODWT_MRA,2)-X_tmp)))])
disp(['Variance:' num2str(sum(var(X_MODWT_MRA))-var(X_tmp))])

X_AT = AT(X_tmp, wname, lev); 
disp('AT')
disp(['Additive:' num2str(sum(abs(sum(X_AT,2)-X_tmp)))])
disp(['Variance:' num2str(sum(var(X_AT))-var(X_tmp))])


% plot all wavelet decompositions
figure
subplot(lev+2,1,1)
plot(X_tmp, 'k');
hold on
sgtitle(['Wavelet: ' num2str(wname)])
xlim([0 N])
ylabel('X');
for is=1:lev+1
    subplot(lev+2,1,is+1)
 
    plot(X_DWT_MRA(:,is),'g');
    xlim([0 N])
    hold on
    plot(X_MODWT(:,is),'r');
    hold on
    plot(X_AT(:,is),'b');
    hold off

    if is==lev+1, ylabel(['A',int2str(is-1)]),
    else ylabel(['D',int2str(is)]), end;
    legend('DWT MRA','MODWT','AT','NumColumns',1,'location','eastoutside')      
end
