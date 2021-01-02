%% RADAR VEGETATION INDEX FROM Compact-POLARIMETRIC SAR DATA

%% Select path of a C2 matrix
[filename, path] = uigetfile('*.*', 'Path selection Time 1');
path
f0 = fopen([path 'config.txt']);
tmp = fgets(f0);
nrows = sscanf(fgets(f0),'%d');
tmp = fgets(f0);
tmp = fgets(f0);
ncols = sscanf(fgets(f0),'%d');

ep = 0;

%% C2 matrix
f1 = fopen([path 'C11.bin'],'rb');
f2 = fopen([path 'C12_real.bin'],'rb');
f3 = fopen([path 'C12_imag.bin'],'rb');
f4 = fopen([path 'C22.bin'],'rb');

c11 = fread(f1,[ncols nrows],'float32') + ep;
c12 = complex( fread(f2,[ncols nrows],'float32') , fread(f3,[ncols nrows],'float32')) + ep;
c21 = conj(c12);
c22 = fread(f4,[ncols nrows],'float32') + ep;

fclose('all');
tic

%% Intitialization
fp22 = zeros(ncols,nrows);
lambda = zeros(ncols,nrows);
dop = zeros(ncols,nrows);

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

chi_in = input('Tau (Ellipticity: +: LC & -: RC): '); % RC = -theta and LC = +theta

for ii=startj:stopj
    for jj=starti:stopi
        
        %% C2 matrix
        
        c11c = mean2(c11(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c12c = mean2(c12(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c21c = mean2(c21(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c22c = mean2(c22(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        C0 = [c11c c12c; c21c c22c];
        
               
        %% LH-LV/RH-RV
        
        % Stokes Parameter
        
        s0 = c11c + c22c;
        s1 = c11c - c22c;
        s2 = c12c + c21c;
        
        if (chi_in >= 0)
            s3 = (1i.*(c12c - c21c)); %The sign is according to RC or LC sign !!
        end
        
        if (chi_in < 0)
            s3 = -(1i.*(c12c - c21c)); %The sign is according to RC or LC sign !!
        end
        
        k11 = s0; k12 = 0; k13 = s2; k14 = 0; 
		k21 = k12; k22 = 0; k23 = 0; k24 = s1;
		k31 = k13; k32 = k23; k33 = 0; k34 = 0;
		k41 = k14; k42 = k24; k43 = k34; k44 = s3;
        
		K_T = (0.5).*[k11 k12 k13 k14; k21 k22 k23 k24; ...
            k31 k32 k33 k34; k41 k42 k43 k44];
        
        %% Stokes vector child products
        
        SC = (s0 - s3)./2;
        OC = (s0 + s3)./2;
        
        dop(ii,jj) = sqrt((s1).^2 + (s2).^2 + (s3).^2)./(s0);
        
        min_sc_oc = min(SC,OC);
        max_sc_oc = max(SC,OC);
        
        %% Random target & Depolarizer
        
        K_depol = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

  
        %% GD_DEPOL
        
        num1 = ((K_T)')*K_depol;
        num = trace(num1);
        den1 = sqrt(abs(trace(((K_T)')*K_T)));
        den2 = sqrt(abs(trace(((K_depol)')*K_depol)));
        den = den1*den2;
        temp_aa = 2.*acosd(num./den);
        GD_t1_depol = temp_aa./180;
        
        %%
        lambda(ii,jj) = (3/2)*GD_t1_depol;
        
        %% GD_VI -- RH-RV/LH-LV
        fp22(ii,jj) = (min_sc_oc/max_sc_oc);
        
               
    end
    fprintf('Processing column: %d \n',ii);
end

%% CpRVI--RH-RV
vi_c = (1 - lambda).*power(fp22, 2.*lambda);

%% File Saving
f_name_100 = strcat(['CpRVI','.bin']);
fileandpath_100=strcat([path, f_name_100]);
fid_100 = fopen(fileandpath_100,'wb');
fwrite(fid_100,vi_c, 'float32');

%% Display VI image
f1 = figure('Name', 'New Vegetation Index');
set(gca,'FontSize',15)
imagesc(vi_c')
axis('image');
colormap('jet');
colorbar('FontSize', 15);
caxis([0 1]);

figname_fd = strcat([path,'CpRVI.fig']);
savefig(figname_fd)
figname_fdpng = strcat([path,'CpRVI.png']);
print(f1,figname_fdpng,'-dpng')



fclose('all');

%% End