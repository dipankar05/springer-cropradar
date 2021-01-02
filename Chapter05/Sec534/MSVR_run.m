%% Multioutput Support Vector Rregression (MSVR)
%% References:
% Tuia, D. and Verrelst, J. and Alonso, L. and Perez-Cruz, F. and Camps-Valls, G. "Multioutput support vector regression for remote sensing biophysical parameter estimation"  IEEE Geoscience and Remote Sensing Letters 8 (4):804-808, 2011
% Ref: https://github.com/DSPKM/DSPKM/blob/master/ch08/msvr/msvr.m
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%% Load calibration and validation data from import
% Load following csv files as column vectors

% calibrationcrop.csv
% LUT_sigma0.csv
% validation_cropparam.csv
% validation_sigma0.csv

X = [C11c,C22c];
Y = [PAIc,WBc];

Xv = [C11v,C22v];
Yv = [PAIv,WBv];

%% Normalizing data sets
[X1,mu1,sig1] = scale(X);
[Y1,mu2,sig2] = scale(Y);
[X1v,mu3,sig3] = scale(Xv);
[Y1v,mu1v,sig1v] = scale(Yv);

Xtrain = X1;       % training set
Ytrain = Y1;       % observed training variables
Xtest  = X1v;   % validation set
Ytest  = Y1v;   % observed variables

%% Training with pairs Xtrain-Ytrain
% Parameters:
%   C: penalization (regularization) parameter
%   epsilon: tolerated error, defines the insensitivity zone
%   sigma: lengthscale of the RBF kernel
%   tol is a very small number (like 1e-15) for early stopping the optimization 

ker     = 'rbf';
epsilon = 0.05;
C       = 4;
sigma   = 1.0;
tol     = 1e-6;
[Beta,NSV,val] = msvr(Xtrain,Ytrain,ker,C,epsilon,sigma,tol);

% Prediction on new test data Xtest
% Build the test kernel matrix
Ktest = kernelmatrix('rbf',Xtest',Xtrain',sigma);
Ypred = Ktest*Beta;

%% Rescaling predicted data with min-max
PAIpred = (sig1v(1).*Ypred(:,1))+mu1v(1);
WBpred = (sig1v(2).*Ypred(:,2))+mu1v(2);

%% Error estimates on validation data
RMSE_PAI = sqrt(mean((PAIv - PAIpred).^2));
RMSE_WB = sqrt(mean((WBv - WBpred).^2));
MAE_PAI = mean(abs((PAIv - PAIpred).^2));
MAE_WB = mean(abs((WBv - WBpred).^2));
[rr, pp]     = corrcoef(PAIv,PAIpred);
R_PAI = rr(1,2);
[rr1, pp1]     = corrcoef(WBv,WBpred);
R_WB = rr1(1,2);

%% Plotting validation data
figure
scatter(PAIv,PAIpred,'ok','filled')
grid on;
axis([0 9 0 9]);
hold on
line([0 9],[0 9],'Color','blue','LineStyle','--')
axis equal;
xlim([0 9])
ylim([0 9])
xlabel('Observed PAI (m^{2} m^{-2})');
ylabel('Estimated PAI (m^{2} m^{-2})');
text(0.2,8.6,sprintf('R=%.2f',R_PAI))
text(0.2,8.2,sprintf('RMSE=%.3f',RMSE_PAI))
text(0.2,7.8,sprintf('MAE=%.3f',MAE_PAI))
box on
%%
figure
scatter(WBv,WBpred,'ok','filled')
grid on;
axis([0 6 0 6]);
axis equal;
hold on
line([0 6],[0 6],'Color','blue','LineStyle','--')
xlim([0 6])
ylim([0 6])
xlabel('Observed WB (kg m^{-2})');
ylabel('Estimated WB (kg m^{-2})');
text(0.2,5.6,sprintf('R=%.2f',R_WB))
text(0.2,5.3,sprintf('RMSE=%.3f',RMSE_WB))
text(0.2,5.0,sprintf('MAE=%.3f',MAE_WB))
box on

%%End
