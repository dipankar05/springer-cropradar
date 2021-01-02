%% Import calibration data from csv files
filename = 'calibrationHH.csv';
delimiter = '';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
HH = dataArray{:, 1};
clearvars filename delimiter formatSpec fileID dataArray ans;
%
filename = 'calibrationVV.csv';
delimiter = '';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
VV = dataArray{:, 1};
clearvars filename delimiter formatSpec fileID dataArray ans;
%
filename = 'calibrationcrop.csv';
delimiter = ',';
formatSpec = '%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
LAI = dataArray{:, 1};
SM = dataArray{:, 2};
clearvars filename delimiter formatSpec fileID dataArray ans;
%
%%
% Radar incident angle=35
theta = 35;
thr=(theta*pi)/180;
% Inputs, tbl=[LAI,SM
tbl = [LAI,SM];
%Defining Water Cloud Model (linear scale)
WCMfun = @(b,x) ((b(1).*(x(:,1).^b(5)).*cos(thr).*(1-exp((-2).*b(2).*x(:,1).*sec(thr)))) +... 
(b(4).*exp(b(3).*x(:,2)).*cos(thr).*exp((-2).*b(2).*x(:,1).*sec(thr))));

% Parameter initialization
% beta0 = [A B C D E]
beta0 = [0.14 1.31 5.76 0.03 0.35];
% Fitting function
mdlHH = fitnlm(tbl,HH,WCMfun,beta0)
mdlVV = fitnlm(tbl,VV,WCMfun,beta0)

%% Validation
yHH = predict(mdlHH,tbl);
rmseHH = sqrt( 1/length(yHH)* sum( (yHH-HH).^2  ));
yVV = predict(mdlVV,tbl);
rmseVV = sqrt( 1/length(yVV)* sum( (yVV-VV).^2  ));
%% Plotting
figure('DefaultAxesFontSize',14)
x0=10;
y0=10;
width=600;
height=600;
set(gcf,'position',[x0,y0,width,height])
scatter(HH,yHH,'sb', 'MarkerFaceColor','b')
% hold on
% scatter(VV,yVV)
xlim([0 0.35])
ylim([0 0.35])
xlabel('Observed $\sigma^\circ_{HH}$','interpreter','latex'); 
ylabel('Estimated $\sigma^\circ_{HH}$','interpreter','latex');
box on
hold on
plot([0,3.5],[0,3.5],'k--','LineWidth',0.5)
str1=sprintf('r = %.2f',corr(yHH,HH));
str2=sprintf('RMSE = %.3f',rmseHH);
annotation('textbox',...
    [0.144333333333334 0.823333333333333 0.194 0.0866666666666699],...
    'String',{str1, str2},...
    'FitBoxToText','on');