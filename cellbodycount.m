function [] = cellbodycount(magnification)
%This function counts the number of cell bodies in a zstack
%   Detailed explanation goes here
clc
prompt='Please choose the first BEFORE image?';
[file,folder]=uigetfile({'*.tif';'*.jpg';'*.*'},prompt);
fb=[folder,file];
cd(folder)
[~,~,ext]=fileparts(fb);
fbefore=dir(['**/*',ext]);
prompt='Please choose the first AFTER image?';
[file,folder]=uigetfile({'*.tif';'*.jpg';'*.*'},prompt);
fa=[folder,file];
cd(folder)
[~,~,ext]=fileparts(fa);
fafter=dir(['**/*',ext]);
prompt='Please Choose the Stack Image';
[stack,fol]=uigetfile({'*.tif';'*.jpg';'*.*'},prompt);
stack=[fol,stack];
Zsta=imread(stack);
prompt='Would you like to draw an ROI?';
answer=questdlg(prompt);
switch answer
    case 'Yes'
        frame=imshow(Zsta);
        h=drawassisted(frame);
        ro=createMask(h);
        M=double(ro);
        prompt='Would you like to draw another ROI?';
        conti=questdlg(prompt);
         while ~isequal(conti,'No')
             frame=imshow(Zsta);
            h=drawassisted(frame);
            ro=createMask(h);
            M=M+double(ro);
            prompt='Would you like to draw another ROI?';
            conti=questdlg(prompt);
         end
             
    case 'No'
        [y,x,i]=size(Zsta);
        M=ones(y,x,i);
        clear y x i
end
clc
close all
z1=size(fbefore,1);
z2=size(fafter,2);
if z1>z2
    msg='There are more images in the before folder than in the after folder.  The program will use the minimum number of images.';
    warning(msg);
elseif z2>z1
    msg='There are more images in the after folder than in the before folder.  The program will use the minimum number of images.';
    warning(msg);
end
z=min([z1,z2]);
A=cell(1,z);
B=cell(1,z);
N=cell(1,z);
MN=cell(1,z);
marea=62.72509992*(magnification/40)^2;
sarea=25.70480614*(magnification/40)^2;
meccen=0.411596748;
seccen=0.258436806;
video = VideoWriter([fafter(1).folder, 'cellbody1.avi']);
video2= VideoWriter([fafter(1).folder,'cellbodyz.avi']);
fr=10;
video.FrameRate=fr;
open(video); %open the file for writing
parfor i=1:z
    %B{i}=imread([fbefore(i).folder,'\',fbefore(i).name]); %Reference image to extract background information
    A{i}=imread([fafter(i).folder,'\',fafter(i).name]);
    %MN{i}=A{i}-B{i};
    MN{i}=A{i}-imread([fbefore(i).folder,'\',fbefore(i).name]);
    N{i}=double(A{i}-B{i}).*M;
    %N{i}=double(A{i}).*M;
    N{i}=uint8(N{i});
end
% NABC=cat(1,N{:,1});
% [a,b,c]=size(NABC);
% NABC=reshape(NABC,a*b*c,1);
% NABC=double(NABC);
% mu=mean(NABC);
% sig=std(NABC);
threshold=55;
[row,col,~]=size(A{z});
N_gray=zeros(row,col,z);
MN_gray=zeros(row,col,z);
parfor i=1:z
    N_gray(:,:,i)=rgb2gray(N{i});
    MN_gray(:,:,i)=rgb2gray(MN{i});
    %imshowpair(N_gray(:,:,i),N{i},'montage')
end
%% thresholding section
imN=zeros(row,col,z);
parfor i=1:z
    %lvl=graythresh(N_gray(:,:,i));
    imN(:,:,i)=N_gray(:,:,i)>threshold;
   % imN(:,:,i)=edge(N_gray(:,:,i),'prewitt');
end
imN=bwareaopen(imN,8); %removes anything with an area less than 8 pixels
imN=imfill(imN,'holes'); %fills in missing pieces of cell bodies
%% Image processing
% maxarea=(marea+1.6*sarea)*(magnification/40);
% minarea=17*(magnification/40);
% ecenmax=.85; %maximum eccentricity value. The lower the eccentricity the more circular the object is
most=zeros(z,1);
supcenters=[];
suprad=[];
%imshow(all)
sigval=.05;
for f=1:z
    all=false(size(imN(:,:,f)));
    cc=bwconncomp(imN(:,:,f),8); %finds number of connected components
    n=cc.NumObjects;
    idx=cc.PixelIdxList;
    Area=zeros(n,1);
    majlength=zeros(n,1);
    minlength=zeros(n,1);
    eccen=zeros(n,1);
    most(f)=n;
    for g=1:n
        REG=regionprops(cc,'Area','MajorAxisLength','MinorAxisLength','Eccentricity');
        Area(g)=REG(g).Area;
        majlength(g)=REG(g).MajorAxisLength;
        minlength(g)=REG(g).MinorAxisLength;
        eccen(g)=REG(g).Eccentricity;
        if (normcdf((Area(g)-marea)/sarea)<sigval|| normcdf((Area(g)-marea)/sarea)>1-sigval) || (normcdf((eccen(g)-meccen)/seccen)<sigval|| normcdf((eccen(g)-meccen)/seccen)>(1-sigval))
            idx{g}=[];
        else
            all(idx{g})=true;
        end
    end
    these=regionprops(all,'Centroid','EquivDiameter');
    [n,~]=size(these);
    centers=[];
    radius=[];
    ma=find(most==max(most),1);
    for i=1:n
        cent=these(i).Centroid;
        centers=cat(1,centers,cent);
        rad=these(i).EquivDiameter/2;
        radius=cat(1,radius,rad);
    end
    supcenters=cat(1,supcenters,centers);
    suprad=cat(1,suprad,radius);
    imshowpair(MN{f},MN_gray(:,:,f),'montage')
    viscircles(supcenters,suprad);
    writeVideo(video,getframe(figure(1))); %write the image to file
end
[n,~]=size(supcenters);
cellnumber=numel(suprad);
car=zeros(1,n);
for i=1:n
    for zn=1:n
        if i~=zn
            if pdist2(supcenters(i,:),supcenters(zn,:),'euclidean')<=.25*max([suprad(i),suprad(zn)])
                car(1,i)=car(1,i)+1;
            end
        end
    end
end
cellnumber=cellnumber-sum(car)/2;
clc
close all
close(video)

imshowpair(MN{ma},MN_gray(:,:,ma),'montage')
viscircles(supcenters,suprad);
video2.FrameRate=fr;
open(video2)
for i=1:z
    imshow(MN{i})
    viscircles(supcenters,suprad);
    writeVideo(video2,getframe(figure(1)));    
end
imshow(Zsta)
viscircles(supcenters,suprad);
print( figure(1),'-bestfit',[fa,'zstacked'],'-dpdf')
close(video2)
disp(['MATLAb found ', num2str(cellnumber),' cell bodies'])
end