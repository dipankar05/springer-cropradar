%% DUAL-POLARIMETRIC RADAR VEGETATION INDEX
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
dprvi = zeros(ncols,nrows);
dop_b = zeros(ncols,nrows);
fp22 = zeros(ncols,nrows);
% rvi = zeros(ncols,nrows);
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
        
        %% C2 matrix
        
        c11c = mean2(c11(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c12c = mean2(c12(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c21c = mean2(c21(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c22c = mean2(c22(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        C0 = [c11c c12c; c21c c22c];
        
        %% GD_VI -- VV-VH
        
        e_v = sort(eig(C0),'descend');
        e_v1 = e_v(1);
        e_v2 = e_v(2);
        x_1 = e_v1/(e_v1 + e_v2);
        fp22(ii,jj) = x_1;
        
        %% dop Barakat-DPRVI
        
        dop_b(ii,jj) = real(sqrt(1 - (4.*det(C0)./(trace(C0).^2))));
        
        dprvi(ii,jj) = 1 - dop_b(ii,jj).*(fp22(ii,jj)); % c22c for S1/c11c for others
        
        %rvi (ii, jj) = 4*c11c./span_c; %RVI dual-pol for comparison for S1
                
    end
    fprintf('Processing column: %d \n',ii);
end

xx = dprvi;

%% File Saving
f_name_100 = strcat(['DPRVI','.bin']);
fileandpath_100=strcat([path, f_name_100]);
fid_100 = fopen(fileandpath_100,'wb');
fwrite(fid_100,xx, 'float32');

%% Display image VI
f1 = figure('Name', 'DPRVI');
set(gca,'FontSize',15)
imagesc(dprvi')
axis('image');
colormap('jet');
colorbar('FontSize', 15);
caxis([0 1]);

figname_fd = strcat([path,'DPRVI.fig']);
savefig(figname_fd)
figname_fdpng = strcat([path,'DPRVI.png']);
print(f1,figname_fdpng,'-dpng')


fclose('all');

