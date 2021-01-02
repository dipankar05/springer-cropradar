%% RADAR VEGETATION INDEX FROM FULL-POLARIMETRIC SAR DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[filename, path] = uigetfile('*.*', 'Path selection Time 1');
path
f0 = fopen([path 'config.txt']);
tmp = fgets(f0);
nrows = sscanf(fgets(f0),'%d');
tmp = fgets(f0);
tmp = fgets(f0);
ncols = sscanf(fgets(f0),'%d');

ep = 0;

f1 = fopen([path 'T11.bin'],'rb');
f2 = fopen([path 'T12_real.bin'],'rb');
f3 = fopen([path 'T12_imag.bin'],'rb');
f4 = fopen([path 'T13_real.bin'],'rb');
f5 = fopen([path 'T13_imag.bin'],'rb');
f6 = fopen([path 'T22.bin'],'rb');
f7 = fopen([path 'T23_real.bin'],'rb');
f8 = fopen([path 'T23_imag.bin'],'rb');
f9 = fopen([path 'T33.bin'],'rb');

t11_T1 = fread(f1,[ncols nrows],'float32') + ep;
t12_T1 = complex( fread(f2,[ncols nrows],'float32') , fread(f3,[ncols nrows],'float32')) + ep;
t21_T1 = conj(t12_T1);
t13_T1 = complex( fread(f4,[ncols nrows],'float32') , fread(f5,[ncols nrows],'float32')) + ep;
t31_T1 = conj(t13_T1);
t22_T1 = fread(f6,[ncols nrows],'float32') + ep;
t23_T1 = complex( fread(f7,[ncols nrows],'float32') , fread(f8,[ncols nrows],'float32')) + ep;
t32_T1 = conj(t23_T1);
t33_T1 = fread(f9,[ncols nrows],'float32') + ep;

fclose('all');
tic

%% Intitialization
span = zeros(ncols,nrows);
temp_rvi = zeros(ncols,nrows);
D = (1/sqrt(2)).*[1 0 1; 1 0 -1; 0 sqrt(2) 0];
fp22 = zeros(ncols,nrows);
GD_t1_t = zeros(ncols,nrows);
GD_t1_d = zeros(ncols,nrows);
GD_t1_rv = zeros(ncols,nrows);
GD_t1_nd = zeros(ncols,nrows);
GD_t1_c = zeros(ncols,nrows);
GD_t1_lh = zeros(ncols,nrows);
GD_t1_rh = zeros(ncols,nrows);
beta = zeros(ncols,nrows);
f = zeros(ncols,nrows);
a = zeros(ncols,nrows);
b = zeros(ncols,nrows);
temp_gamma = zeros(ncols,nrows);
t_d = zeros(46,1);
t_nd = zeros(46,1);
t_t = zeros(46,1);
t_c = zeros(46,1);
theta_map = zeros(ncols,nrows);

%% for window processing

wsi=input('Window Size: ');
wsj = wsi; % Number of columns in the window

inci=fix(wsi/2); % Up & down movement margin from the central row
incj=fix(wsj/2); % Left & right movement from the central column
% Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999

starti=fix(wsi/2)+1; % Starting row for window processing
startj=fix(wsj/2)+1; % Starting column for window processing

stopi= nrows-inci; % Stop row for window processing
stopj= ncols-incj; % Stop column for window processing

t = cputime;

for ii=startj:stopj
    for jj=starti:stopi
        
        t11s = mean2(t11_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t12s = mean2(t12_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t13s = mean2(t13_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        t21s = mean2(t21_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t22s = mean2(t22_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t23s = mean2(t23_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        t31s = mean2(t31_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t32s = mean2(t32_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t33s = mean2(t33_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        %Coherency matrix
        T_T1 = [t11s t12s t13s; t21s t22s t23s; t31s t32s t33s];
        
        %Coherency matrix
        C_T1 = (D')*T_T1*D;
        
        span(ii,jj) = t11s + t22s + t33s;
        temp_span = span(ii,jj);
        
        Temp_T1 = T_T1;
        
        t11 = Temp_T1(1,1); t12 = Temp_T1(1,2); t13 = Temp_T1(1,3);
        t21 = conj(t12); t22 = Temp_T1(2,2); t23 = Temp_T1(2,3);
        t31 = conj(t13); t32 = conj(t23); t33 = Temp_T1(3,3);
        
       
        %% Kennaugh Matrix
        
        m11 = t11+t22+t33; m12 = t12+t21; m13 = t13+t31; m14 = -1i*(t23 - t32);
        m21 = t12+t21; m22 = t11+t22-t33; m23 = t23+t32; m24 = -1i*(t13-t31);
        m31 = t13+t31; m32 = t23+t32; m33 = t11-t22+t33; m34 = 1i*(t12-t21);
        m41 = -1i*(t23-t32); m42 = -1i*(t13-t31); m43 = 1i*(t12-t21); m44 = -t11+t22+t33;
        
        M_T = 0.5.*[m11 m12 m13 m14; m21 m22 m23 m24; m31 m32 m33 m34; m41 m42 m43 m44];
        
        %% Elementary targets
        
        M_d = [1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 1];
        M_nd = [0.625 0.375 0 0; 0.375 0.625 0 0; 0 0 -0.5 0; 0 0 0 0.5];
        M_t = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1];
        M_c = [0.625 0.375 0 0; 0.375 0.625 0 0; 0 0 0.5 0; 0 0 0 -0.5];
        M_lh = [1 0 0 -1; 0 0 0 0; 0 0 0 0; -1 0 0 1];
        M_rh = [1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 1];
        
        M_T_theta = M_T;
        
        %% GVSM
        
        t011 = M_T_theta(1,1) + M_T_theta(2,2) + M_T_theta(3,3) - M_T_theta(4,4);
        t012 = M_T_theta(1,2) - 1i*M_T_theta(3,4);
        t013 = M_T_theta(1,3) + 1i*M_T_theta(2,4);
        t021 = conj(t012);
        t022 = M_T_theta(1,1) + M_T_theta(2,2) - M_T_theta(3,3) + M_T_theta(4,4);
        t023 = M_T_theta(2,3) +1i*M_T_theta(1,4);
        t031 = conj(t013);
        t032 = conj(t023);
        t033 = M_T_theta(1,1) - M_T_theta(2,2) + M_T_theta(3,3) + M_T_theta(4,4);
        
        %% T to C
        
        T0 = [t011./2 t012 t013; t021 t022./2 t023; t031 t032 t033./2];
        C0 = (D')*T0*D;
        
        %% Gamma/Rho
        
        gamma = C0(1,1)/C0(3,3); rho = 1/3;
        temp_gamma(ii,jj) = gamma; % variable to save
        
        %% Covariance matrix
        
        c11 = gamma; c12 = 0; c13 = rho*sqrt(gamma);
        c21 = 0; c22 = 0.5*(1 + gamma) - rho*sqrt(gamma); c23 = 0;
        c31 = conj(rho)*sqrt(gamma); c32 = 0; c33 = 1;
        
        R = (3/2)*(1 + gamma) - rho*sqrt(gamma);
        C1 = (1./R).*[c11 c12 c13; c21 c22 c23; c31 c32 c33];
        
        %% Coherency matrix
        
        T1 = D*C1*(D');
        
        t11 = T1(1,1); t12 = T1(1,2); t13 = T1(1,3);
        t21 = T1(2,1); t22 = T1(2,2); t23 = T1(2,3);
        t31 = T1(3,1); t32 = T1(3,2); t33 = T1(3,3);
        
        m11 = t11+t22+t33; m12 = t12+t21; m13 = t13+t31; m14 = -1i*(t23 - t32);
        m21 = t12+t21; m22 = t11+t22-t33; m23 = t23+t32; m24 = -1i*(t13-t31);
        m31 = t13+t31; m32 = t23+t32; m33 = t11-t22+t33; m34 = 1i*(t12-t21);
        m41 = -1i*(t23-t32); m42 = -1i*(t13-t31); m43 = 1i*(t12-t21); m44 = -t11+t22+t33;
        
        %% Generalized Random Volume (Antropov et al.)
        
        M_rv = real([m11 m12 m13 m14; m21 m22 m23 m24; m31 m32 ...
            m33 m34; m41 m42 m43 m44]);
        
        f(ii,jj) = 1;
        
        %% GD Volume
        
        num1 = ((M_T_theta)')*M_rv; % volume
        num = trace(num1);
        den1 = sqrt(abs(trace(((M_T_theta)')*M_T_theta)));
        den2 = sqrt(abs(trace(((M_rv)')*M_rv)));
        den = den1*den2;
        temp_aa = 2.*acosd(num./den);
        GD_t1_rv(ii,jj) = temp_aa./180;
        
        %% GD ALL
        
        num1 = ((M_T_theta)')*M_c; % cylinder
        num = trace(num1);
        den1 = sqrt(abs(trace(((M_T_theta)')*M_T_theta)));
        den2 = sqrt(abs(trace(((M_c)')*M_c)));
        den = den1*den2;
        temp_aa = 2.*acosd(num./den);
        GD_t1_c(ii,jj) = temp_aa./180;
        
        num1 = ((M_T_theta)')*M_t; % trihedral
        num = trace(num1);
        den1 = sqrt(abs(trace(((M_T_theta)')*M_T_theta)));
        den2 = sqrt(abs(trace(((M_t)')*M_t)));
        den = den1*den2;
        temp_aa = 2.*acosd(num./den);
        GD_t1_t(ii,jj) = temp_aa./180;
        
        num1 = ((M_T_theta)')*M_d; % dihedral
        num = trace(num1);
        den1 = sqrt(abs(trace(((M_T_theta)')*M_T_theta)));
        den2 = sqrt(abs(trace(((M_d)')*M_d)));
        den = den1*den2;
        temp_aa = 2.*acosd(num./den);
        GD_t1_d(ii,jj) = temp_aa./180;
        
        num1 = ((M_T_theta)')*M_nd; % n-dihedral
        num = trace(num1);
        den1 = sqrt(abs(trace(((M_T_theta)')*M_T_theta)));
        den2 = sqrt(abs(trace(((M_nd)')*M_nd)));
        den = den1*den2;
        temp_aa = 2.*acosd(num./den);
        GD_t1_nd(ii,jj) = temp_aa./180;
        
        %% VI
        
        t_t = GD_t1_t(ii,jj);
        t_d = GD_t1_d(ii,jj);
        t_c = GD_t1_c(ii,jj);
        t_nd = GD_t1_nd(ii,jj);
        
        a(ii,jj) = max([t_t, t_d, t_c, t_nd]);
        b(ii,jj) = min([t_t, t_d, t_c, t_nd]);
        beta(ii,jj) = (b(ii,jj)/a(ii,jj))^2;
        
        %% RVI
        
        e_v = sort(eig(T_T1),'descend');
        e_v1 = e_v(1); e_v2 = e_v(2); e_v3 = e_v(3);
        
        p1 = e_v1./(e_v1 + e_v2 + e_v3);
        p2 = e_v2./(e_v1 + e_v2 + e_v3);
        p3 = e_v3./(e_v1 + e_v2 + e_v3);
        
        p1(p1 < 0) = 0; p2(p2 < 0) = 0; p3(p3 < 0) = 0;
        
        p1(p1 > 1) = 1; p2(p2 > 1) = 1; p3(p3 > 1) = 1;
        
        temp_rvi(ii,jj) = (4.*p3)./(p1 + p2 + p3);
        
    end
    fprintf('Processing column: %d \n',ii);
end

%% GRVI

vi = power(beta, GD_t1_rv).*(1 - (1./f).*GD_t1_rv);

idx1 = (GD_t1_rv > f);
vi(idx1) = 0;
vi(~idx1) = vi(~idx1);

%% RVI scaled (0 - 1)

rvi = temp_rvi;
idx = (rvi > 1);
rvi(idx) = (3/4).*rvi(idx);
rvi(~idx) = rvi(~idx);


%% File Saving

f_name_100 = strcat(['GRVI','.bin']);
fileandpath_100=strcat([path, f_name_100]);
fid_100 = fopen(fileandpath_100,'wb');
fwrite(fid_100,vi, 'float32');

f_name_200 = strcat(['RVI','.bin']);
fileandpath_200=strcat([path, f_name_200]);
fid_200 = fopen(fileandpath_200,'wb');
fwrite(fid_200,rvi, 'float32');

%% Image display
f1 = figure('Name', 'GRVI');
set(gca,'FontSize',15)
imagesc(vi')
axis('image');
colormap('jet');
colorbar('FontSize', 15);
caxis([0 1]);

figname_fd = strcat([path,'GRVI.fig']);
savefig(figname_fd)
figname_fdpng = strcat([path,'GRVI.png']);
print(f1,figname_fdpng,'-dpng')
%%
f2 = figure('Name', 'RVI');
set(gca,'FontSize',15)
imagesc(rvi')
axis('image');
colormap('jet');
colorbar('FontSize', 15);
caxis([0 1]);

figname_fd1 = strcat([path,'RVI.fig']);
savefig(figname_fd1)
figname_fd1png = strcat([path,'RVI.png']);
print(f2,figname_fd1png,'-dpng')

fclose('all');
