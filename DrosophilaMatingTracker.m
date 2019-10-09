clear
clc
close all
%warning('off','MATLAB:nearlySingularMatrix')
[file,path]=uigetfile({'*.avi';'*.mp4';'*.*'},'Choose The Video File');
filepath=[path,file];
obj.reader=vision.VideoFileReader(filepath,'VideoOutputDataType','uint8');
videoPlayer=vision.VideoPlayer;
%msav=[path,'Tracking Results'];
ap=strrep(file,'.avi','');
ap=strrep(ap,'.mp4','');
vidsav=[path,ap,'_Tracking_Results.avi'];
msav=[path,ap,'_Tracking Results'];
dummy=VideoReader(filepath);
numframes=dummy.NumberOfFrames; %#ok<VIDREAD>
fps=dummy.FrameRate;
dt=1/fps;
time=0:dt:dt*numframes;
cou=1;
wei=0;
Tmat=diag(ones(1,8));
Tmat(1,3)=dt;
Tmat(1,5)=dt^2/2;
Tmat(1,7)=dt^3/6;
Tmat(2,4)=dt;
Tmat(2,6)=dt^2/2;
Tmat(2,8)=dt^3/6;
Tmat(3,5)=dt;
Tmat(4,6)=dt;
Tmat(3,6)=dt^2/2;
Tmat(4,8)=dt^2/2;
Tmat(5,7)=dt;
Tmat(6,8)=dt;
mem=.25;
QQ=mem*diag(ones(1,8));
QQ(1,3)=mem/2;
QQ(3,1)=QQ(1,3);
QQ(2,4)=mem/2;
QQ(4,2)=QQ(2,4);
QQ(3,5)=mem/2;
QQ(5,3)=QQ(3,5);
QQ(4,6)=mem/2;
QQ(6,4)=QQ(4,6);
QQ(5,7)=mem/2;
QQ(7,5)=QQ(5,7);
QQ(6,8)=mem/2;
QQ(8,6)=QQ(6,8);
WindowSize=20;
while exist(vidsav,'file')
    vidsav=[msav,'_',num2str(cou),'.avi'];
    cou=cou+1;
end
name=ap;
clear cou
VID=VideoWriter(vidsav);
VID.FrameRate=10; %Frames per second
VID.Quality=100;
open(VID)
%Drawing ROI
prompt='Would you like to draw an ROI?';
answer=questdlg(prompt);
%lvl=.45;%Intensity Threshold
switch answer
    case  'Yes'
        frame=obj.reader.step();
        b=imshow(frame);
        h=drawassisted(b);
        ro=createMask(h);
        close Figure 1
        
    case 'No'
        frame=obj.reader();
        [bro,co,~]=size(frame);
        ro=ones(bro,co,1);
        
end
prompt={'Enter the number of animals being tracked:'};
dlgTitle='Input';
dims = [1 35];
definput ={'36'};
answer = inputdlg(prompt,dlgTitle,dims,definput);
nanimals=str2double(answer{1});
prompt='Would you like to cross validate the results? (Cross validation is highly recommended)';
valans=questdlg(prompt);
switch valans
    case 'Yes'
        valcheck=1;
    case 'No'
        valcheck=0;
end

count=1;
X=zeros(nanimals,numframes-1);
Y=zeros(nanimals,numframes-1);
AX=zeros(nanimals,numframes-1);
AY=zeros(nanimals,numframes-1);
VX=zeros(nanimals,numframes-1);
VY=zeros(nanimals,numframes-1);
JX=zeros(nanimals,numframes-1);
JY=zeros(nanimals,numframes-1);
fly=struct('ID',[],'position',[],'velocity',[],'acceleration',[],'jerk',[],'kalman',[],'cov',[],'error',[]);
lvl=.29;
while ~isDone(obj.reader)
    videoFrame = obj.reader();
    NFrame=rgb2gray(videoFrame);
    %Framenew=double(NFrame).*M;
    Framenew=bsxfun(@times,NFrame, cast(ro,class(NFrame)));
    Framenew(Framenew==0)=10000000000000000;
    BW=imbinarize(Framenew,lvl);
    mask=boundarymask(BW);
    mask=imfill(mask,'holes');
    %mask=mask==0;
    %BW=imbinarize(Framenew);
    %BW=BW==0;
    %BW=imfill(BW,'holes');
    %BW=bwareaopen(BW,20);
    Reg=regionprops(mask,'Area','Centroid','MajorAxisLength','MinorAxisLength');
    n=size(Reg,1);
    AS=cat(1,Reg.Area);
    CE=cat(1,Reg.Centroid);
    MajL=cat(1,Reg.MajorAxisLength);
    MinL=cat(1,Reg.MinorAxisLength);
    MArea=5000;
    MIrea=80;
    CE(AS>MArea,:)=[];
    MajL(AS>MArea)=[];
    MinL(AS>MArea)=[];
    AS(AS>MArea)=[];
    CE(AS<MIrea,:)=[];
    MajL(AS<MIrea)=[];
    MinL(AS<MIrea)=[];
    AS(AS<MIrea)=[];
    Diam=mean([MajL, MinL],2);
    rad=Diam/2;
    pos=[CE,rad];
    %frame=insertShape(videoFrame,'Circle',pos);
    numTracks=size(pos,1);
    prob=zeros(1,nanimals);
    if count==1
        IS=[];
        JcS=IS;
        if numTracks<nanimals
            clc
            disp('This frame does not contain all of the flies, searching for another reference frame....')
            %pause(.1)
            clc
            continue
        end
        %initialize fly tracks
        frame=videoFrame;
        for i=1:numTracks
            clc
            fly(i).ID=['Fly #', num2str(i)];
            x=CE(i,1);
            y=CE(i,2);
            vx=0;
            vy=0;
            ax=0;
            ay=0;
            jx=0;
            jy=0;
            fly(i).position=[x y];
            fly(i).velocity=[0 0];
            fly(i).acceleration=[0 0];
            fly(i).jerk=[0 0];
            X(i,count)=fly(i).position(1);
            Y(i,count)=fly(i).position(2);
            VX(i,count)=fly(i).velocity(1);
            VY(i,count)=fly(i).velocity(2);
            AX(i,count)=fly(i).velocity(1);
            AY(i,count)=fly(i).velocity(2);
            JX(i,count)=fly(i).jerk(1);
            JY(i,count)=fly(i).jerk(2);
            kalman=[x y vx vy ax ay jx jy]';
            observed=[x y vx vy ax ay jx jy]';
            fly(i).kalman=kalman;
            fly(i).observed=observed;
            %fly(i).cov=diag(ones(1,8));
            fly(i).cov=QQ;
            prob(i)=100;
            fly(i).error=[50,50,50,50,50,50];
        end
    elseif numTracks==nanimals
        IS=[];
        JcS=IS;
        prob=zeros(1,nanimals);
        for i=1:nanimals
            cost=zeros(1,numTracks);
            stx=mean(fly(i).error(:,1));
            sty=mean(fly(i).error(:,2));
            stvx=mean(fly(i).error(:,3));
            stvy=mean(fly(i).error(:,4));
            stax=mean(fly(i).error(:,5));
            stay=mean(fly(i).error(:,6));
            A=Tmat*fly(i).kalman;
            lg=size(CE,1);
            VfX=filter(ones(1,WindowSize)/WindowSize,1,VX(i,1:count));
            VfY=filter(ones(1,WindowSize)/WindowSize,1,VY(i,1:count));
            AfX=filter(ones(1,WindowSize)/WindowSize,1,AX(i,1:count));
            AfY=filter(ones(1,WindowSize)/WindowSize,1,AY(i,1:count));
            for j=1:lg
                vx=(CE(j,1)-A(1))/dt;
                vy=(CE(j,2)-A(2))/dt;
                ax=(vx-A(3))/dt;
                ay=(vy-A(4))/dt;
                ax=filter(ones(1,WindowSize)/WindowSize,1,[AfX, ax]);
                ax=ax(end);
                ay=filter(ones(1,WindowSize)/WindowSize,1,[AfY, ay]);
                ay=ay(end);
                vx=filter(ones(1,WindowSize)/WindowSize,1,[VfX, vx]);
                vx=vx(end);
                vy=filter(ones(1,WindowSize)/WindowSize,1,[VfY, vy]);
                vy=vy(end);
                cost(j)=1/((2*pi)^3*stx*sty*stvx*stvy*stax*stay)*exp(-1/2*(((CE(j,1)-A(1))/stx)^2+((CE(j,2)-A(2))/sty)^2+wei*((vx-A(3))/stvx)^2+wei*((vy-A(4))/stvy)^2+wei^2*((ax-A(5))/stax)^2+wei^2*((ay-A(6))/stay)^2));
            end
            b=find(cost==max(cost),1,'first');
            prob(i)=cost(b)/sum(cost)*100;
            x=CE(b,1);
            y=CE(b,2);
            CE(b,:)=[];
            vx=(fly(i).position(1)-x)/dt;
            vy=(fly(i).position(2)-y)/dt;
            ax=(vx-fly(i).velocity(1))/dt;
            ay=(vy-fly(i).velocity(2))/dt;
            jx=(ax-fly(i).acceleration(1))/dt;
            jy=(ay-fly(i).acceleration(2))/dt;
            AX(i,count)=ax;
            AY(i,count)=ay;
            VX(i,count)=vx;
            VY(i,count)=vy;
            X(i,count)=x;
            Y(i,count)=y;
            JX(i,count)=jx;
            JY(i,count)=jy;
            fly(i).position=[x y];
            fly(i).jerk=[jx jy];
            VfX=filter(ones(1,WindowSize)/WindowSize,1,VX(i,1:count));
            VfY=filter(ones(1,WindowSize)/WindowSize,1,VY(i,1:count));
            AfX=filter(ones(1,WindowSize)/WindowSize,1,AX(i,1:count));
            AfY=filter(ones(1,WindowSize)/WindowSize,1,AY(i,1:count));
            JfX=filter(ones(1,WindowSize)/WindowSize,1,JX(i,1:count));
            JfY=filter(ones(1,WindowSize)/WindowSize,1,JY(i,1:count));
            VX(i,count)=VfX(count);
            VY(i,count)=VfY(count);
            AX(i,count)=AfX(count);
            AY(i,count)=AfY(count);
            JX(i,count)=JfX(count);
            JY(i,count)=JfY(count);
            fly(i).velocity=[VX(i,count) VY(i,count)];
            fly(i).acceleration=[AX(i,count) AY(i,count)];
            fly(i).jerk=[JX(i,count) JY(i,count)];
            observed=[x y VX(i,count) VY(i,count) AX(i,count) AY(i,count) JX(i,count) JY(i,count)]';
            fly(i).observed=observed;
        end
    elseif numTracks<nanimals %one or more animals were not detected
        pot=nanimals-size(CE,1) ;
        set=zeros(1,pot);
        IS=1:nanimals;
        JS=IS;
        for i=1:numTracks
            PoS=CE(i,:);
            isE=numel(IS);
            cost=zeros(1,isE);
            for is=1:isE
                VfX=filter(ones(1,WindowSize)/WindowSize,1,VX(IS(is),1:count));
                VfY=filter(ones(1,WindowSize)/WindowSize,1,VY(IS(is),1:count));
                AfX=filter(ones(1,WindowSize)/WindowSize,1,AX(IS(is),1:count));
                AfY=filter(ones(1,WindowSize)/WindowSize,1,AY(IS(is),1:count));
                A=Tmat*fly(IS(is)).kalman;
                stx=mean(fly(IS(is)).error(:,1));
                sty=mean(fly(IS(is)).error(:,2));
                stvx=mean(fly(IS(is)).error(:,3));
                stvy=mean(fly(IS(is)).error(:,4));
                stax=mean(fly(IS(is)).error(:,5));
                stay=mean(fly(IS(is)).error(:,6));
                vx=(CE(i,1)-A(1))/dt;
                vy=(CE(i,2)-A(2))/dt;
                ax=(vx-A(3))/dt;
                ay=(vy-A(4))/dt;
                ax=filter(ones(1,WindowSize)/WindowSize,1,[AfX, ax]);
                ax=ax(end);
                ay=filter(ones(1,WindowSize)/WindowSize,1,[AfY, ay]);
                ay=ay(end);
                vx=filter(ones(1,WindowSize)/WindowSize,1,[VfX, vx]);
                vx=vx(end);
                vy=filter(ones(1,WindowSize)/WindowSize,1,[VfY, vy]);
                vy=vy(end);
                cost(is)=1/((2*pi)^3*stx*sty*stvx*stvy*stax*stay)*exp(-1/2*(((CE(i,1)-A(1))/stx)^2+((CE(i,2)-A(2))/sty)^2+wei*((vx-A(3))/stvx)^2+wei*((vy-A(4))/stvy)^2+wei^2*((ax-A(5))/stax)^2+wei^2*((ay-A(6))/stay)^2));
            end
            p=find(cost==max(cost),1,'first');
            P=IS(p);
            IS(p)=[];
            prob(P)=cost(p)/sum(cost)*100;
            x=CE(i,1);
            y=CE(i,2);
            vx=(fly(P).position(1)-x)/dt;
            vy=(fly(P).position(2)-y)/dt;
            ax=(vx-fly(P).velocity(1))/dt;
            ay=(vy-fly(P).velocity(2))/dt;
            jx=(ax-fly(P).acceleration(1))/dt;
            jy=(ay-fly(P).acceleration(2))/dt;
            fly(P).position=[x y];
            AX(P,count)=ax;
            AY(P,count)=ay;
            VX(P,count)=vx;
            VY(P,count)=vy;
            X(P,count)=x;
            Y(P,count)=y;
            JX(P,count)=jx;
            JY(P,count)=jy;
            VfX=filter(ones(1,WindowSize)/WindowSize,1,VX(P,1:count));
            VfY=filter(ones(1,WindowSize)/WindowSize,1,VY(P,1:count));
            AfX=filter(ones(1,WindowSize)/WindowSize,1,AX(P,1:count));
            AfY=filter(ones(1,WindowSize)/WindowSize,1,AY(P,1:count));
            JfX=filter(ones(1,WindowSize)/WindowSize,1,JX(P,1:count));
            JfY=filter(ones(1,WindowSize)/WindowSize,1,JY(P,1:count));
            VX(P,count)=VfX(count);
            VY(P,count)=VfY(count);
            AX(P,count)=AfX(count);
            AY(P,count)=AfY(count);
            JX(P,count)=JfX(count);
            JY(P,count)=JfY(count);
            fly(P).velocity=[VX(P,count) VY(P,count)];
            fly(P).acceleration=[AX(P,count) AY(P,count)];
            fly(P).jerk=[JX(P,count) JY(P,count)];
            observed=[x y VX(P,count) VY(P,count) AX(P,count) AY(P,count) JX(P,count) JY(P,count)]';
            fly(P).observed=observed;
            fly(P).position=[x y];
        end
        JcS=IS;
        mu=30;
        for i=1:pot
            xl=sqrt((X(IS(i),count-1)-X(JS(~(JS==IS(i))),count-1)).^2+(Y(IS(i),count-1)-Y(JS(~(JS==IS(i))),count-1)).^2);
            b=1/(mu*sqrt(2*pi))*exp((-1/2)*(xl/mu).^2);
            ps=1/(mu*sqrt(2*pi))*exp(-1/2);
            ix=find(b>=ps,1,'first');
            if ~isempty(ix)
                %disp('It works')
                set(i)=1;
            end
        end
        for a=1:pot
            i=IS(a);
            if set(a)==0
                fly(i).position=[fly(i).kalman(1) fly(i).kalman(2)];
                fly(i).velocity=[fly(i).kalman(3) fly(i).kalman(4)];
                X(i,count)=fly(i).kalman(1);
                Y(i,count)=fly(i).kalman(2);
                VX(i,count)=fly(i).kalman(3);
                VY(i,count)=fly(i).kalman(4);
                AX(i,count)=fly(i).kalman(5);
                AY(i,count)=fly(i).kalman(6);
                JX(i,count)=fly(i).kalman(7);
                JY(i,count)=fly(i).kalman(8);
                fly(i).acceleration=[fly(i).kalman(5) fly(i).kalman(6)];
                fly(i).jerk=[fly(i).kalman(7) fly(i).kalman(8)];
                prob(i)=85;
            else
                cost=zeros(1,numTracks);
                stx=mean(fly(i).error(:,1));
                sty=mean(fly(i).error(:,2));
                stvx=mean(fly(i).error(:,3));
                stvy=mean(fly(i).error(:,4));
                stax=mean(fly(i).error(:,5));
                stay=mean(fly(i).error(:,6));
                A=Tmat*fly(i).kalman;
                lg=numTracks;
                VfX=filter(ones(1,WindowSize)/WindowSize,1,VX(i,1:count));
                VfY=filter(ones(1,WindowSize)/WindowSize,1,VY(i,1:count));
                AfX=filter(ones(1,WindowSize)/WindowSize,1,AX(i,1:count));
                AfY=filter(ones(1,WindowSize)/WindowSize,1,AY(i,1:count));  
                for j=1:numTracks
                    vx=(CE(j,1)-A(1))/dt;
                    vy=(CE(j,2)-A(2))/dt;
                    ax=(vx-A(3))/dt;
                    ay=(vy-A(4))/dt;
                    ax=filter(ones(1,WindowSize)/WindowSize,1,[AfX, ax]);
                    ax=ax(end);
                    ay=filter(ones(1,WindowSize)/WindowSize,1,[AfY, ay]);
                    ay=ay(end);
                    vx=filter(ones(1,WindowSize)/WindowSize,1,[VfX, vx]);
                    vx=vx(end);
                    vy=filter(ones(1,WindowSize)/WindowSize,1,[VfY, vy]);
                    vy=vy(end);
                    cost(j)=1/((2*pi)^3*stx*sty*stvx*stvy*stax*stay)*exp(-1/2*(((CE(j,1)-A(1))/stx)^2+((CE(j,2)-A(2))/sty)^2+wei*((vx-A(3))/stvx)^2+wei*((vy-A(4))/stvy)^2+wei^2*((ax-A(5))/stax)^2+wei^2*((ay-A(6))/stay)^2));
                end
                b=find(cost==max(cost),1,'first');
                prob(i)=cost(b)/sum(cost)*100;
                x=CE(b,1);
                y=CE(b,2);
                vx=(fly(i).position(1)-x)/dt;
                vy=(fly(i).position(2)-y)/dt;
                ax=(vx-fly(i).velocity(1))/dt;
                ay=(vy-fly(i).velocity(2))/dt;
                jx=(ax-fly(i).acceleration(1))/dt;
                jy=(ay-fly(i).acceleration(2))/dt;
                AX(i,count)=ax;
                AY(i,count)=ay;
                VX(i,count)=vx;
                VY(i,count)=vy;
                X(i,count)=x;
                Y(i,count)=y;
                JX(i,count)=jx;
                JY(i,count)=jy;
                fly(i).position=[x y];
                fly(i).jerk=[jx jy];
                VfX=filter(ones(1,WindowSize)/WindowSize,1,VX(i,1:count));
                VfY=filter(ones(1,WindowSize)/WindowSize,1,VY(i,1:count));
                AfX=filter(ones(1,WindowSize)/WindowSize,1,AX(i,1:count));
                AfY=filter(ones(1,WindowSize)/WindowSize,1,AY(i,1:count));
                JfX=filter(ones(1,WindowSize)/WindowSize,1,JX(i,1:count));
                JfY=filter(ones(1,WindowSize)/WindowSize,1,JY(i,1:count));
                VX(i,count)=VfX(count);
                VY(i,count)=VfY(count);
                AX(i,count)=AfX(count);
                AY(i,count)=AfY(count);
                JX(i,count)=JfX(count);
                JY(i,count)=JfY(count);
                fly(i).velocity=[VX(i,count) VY(i,count)];
                fly(i).acceleration=[AX(i,count) AY(i,count)];
                fly(i).jerk=[JX(i,count) JY(i,count)];
                observed=[x y VX(i,count) VY(i,count) AX(i,count) AY(i,count) JX(i,count) JY(i,count)]';
                fly(i).observed=observed;
            end
        end
        IS(set==1)=[];
    elseif numTracks>nanimals %an artifact was detected or one fly was split in multiple pieces
        prob=zeros(1,nanimals);
        IS=[];
        JcS=IS;
        gai=numTracks-nanimals;
        T=1:numTracks;
        for i=1:nanimals
            cost=zeros(1,numTracks);
            VfX=filter(ones(1,WindowSize)/WindowSize,1,VX(i,1:count));
            VfY=filter(ones(1,WindowSize)/WindowSize,1,VY(i,1:count));
            AfX=filter(ones(1,WindowSize)/WindowSize,1,AX(i,1:count));
            AfY=filter(ones(1,WindowSize)/WindowSize,1,AY(i,1:count));
            for j=T %1:numTracks
                A=Tmat*fly(i).kalman;
                stx=mean(fly(i).error(:,1));
                sty=mean(fly(i).error(:,2));
                stvx=mean(fly(i).error(:,3));
                stvy=mean(fly(i).error(:,4));
                stax=mean(fly(i).error(:,5));
                stay=mean(fly(i).error(:,6));
                vx=(CE(j,1)-A(1))/dt;
                vy=(CE(j,2)-A(2))/dt;
                ax=(vx-A(3))/dt;
                ay=(vy-A(4))/dt;
                ax=filter(ones(1,WindowSize)/WindowSize,1,[AfX, ax]);
                ax=ax(end);
                ay=filter(ones(1,WindowSize)/WindowSize,1,[AfY, ay]);
                ay=ay(end);
                vx=filter(ones(1,WindowSize)/WindowSize,1,[VfX, vx]);
                vx=vx(end);
                vy=filter(ones(1,WindowSize)/WindowSize,1,[VfY, vy]);
                vy=vy(end);
                cost(j)=1/((2*pi)^3*stx*sty*stvx*stvy*stax*stay)*exp(-1/2*(((CE(j,1)-A(1))/stx)^2+((CE(j,2)-A(2))/sty)^2+wei*((vx-A(3))/stvx)^2+wei*((vy-A(4))/stvy)^2+wei^2*((ax-A(5))/stax)^2+wei^2*((ay-A(6))/stay)^2));
            end
            a=find(cost==max(cost),1,'first');
            T(T==a)=[];
            prob(i)=cost(a)/sum(cost)*100;
            x=CE(a,1);
            y=CE(a,2);
            vx=(fly(i).position(1)-x)/dt;
            vy=(fly(i).position(2)-y)/dt;
            ax=(vx-fly(i).velocity(1))/dt;
            ay=(vy-fly(i).velocity(2))/dt;
            jx=(ax-fly(i).acceleration(1))/dt;
            jy=(ay-fly(i).acceleration(2))/dt;
            fly(i).position=[x y];
            AX(i,count)=ax;
            AY(i,count)=ay;
            VX(i,count)=vx;
            VY(i,count)=vy;
            X(i,count)=x;
            Y(i,count)=y;
            JX(i,count)=jx;
            JY(i,count)=jy;
            VfX=filter(ones(1,WindowSize)/WindowSize,1,VX(i,1:count));
            VfY=filter(ones(1,WindowSize)/WindowSize,1,VY(i,1:count));
            AfX=filter(ones(1,WindowSize)/WindowSize,1,AX(i,1:count));
            AfY=filter(ones(1,WindowSize)/WindowSize,1,AY(i,1:count));
            JfX=filter(ones(1,WindowSize)/WindowSize,1,JX(i,1:count));
            JfY=filter(ones(1,WindowSize)/WindowSize,1,JY(i,1:count));
            VX(i,count)=VfX(count);
            VY(i,count)=VfY(count);
            AX(i,count)=AfX(count);
            AY(i,count)=AfY(count);
            JX(i,count)=JfX(count);
            JY(i,count)=JfY(count);
            fly(i).velocity=[VX(i,count) VY(i,count)];
            fly(i).acceleration=[AX(i,count) AY(i,count)];
            fly(i).jerk=[JX(i,count) JY(i,count)];
            observed=[x y VX(i,count) VY(i,count) AX(i,count) AY(i,count) JX(i,count) JY(i,count)]';
            fly(i).observed=observed;
        end
    end
    Brame=videoFrame;
    p=zeros(1,nanimals);
    CE=pos(:,1:2);
    
    Dros=zeros(nanimals,3);
    labels=cell(1,nanimals);
    ra=mean(rad);
    prob(isnan(prob))=100;
    for n=1:nanimals
        if ~ismember(n,IS)
            a=find(CE==fly(n).position,1);
            labels{n}=['Fly #', num2str(n),' ', num2str(prob(n)),'%'];
            Dros(n,:)=pos(a,:);
        else
            Dros(n,1:2)=fly(n).position;
            Dros(n,3)=ra;
            labels{n}=['Predicted location of Fly #', num2str(n),' ', num2str(prob(n)),'%'];
        end
        
    end
    p(p==0)=[];
    Brame=insertObjectAnnotation(Brame,'circle',Dros,labels,'TextBoxOpacity',0.45,'FontSize',12);
    if valcheck==1
        AID=videoFrame;
        v=vmagsign(VX(:,count),VY(:,count));
        apa=vmagsign(AX(:,count),AY(:,count));
        amax=1e9;
        vmax=1e8;
        vr=find(v>=vmax);
        ar=find(apa>=amax);
        %prr=find(prob<=75);
        if ~exist('br','var')
            br=[];
        end
        wr=union(ar,vr);
        wr=union(wr,br);
        %wr=union(wr,prr);
        if ~isempty(wr)
            Brame=videoFrame;
            uiwait(msgbox('Validation Required','User Validation','modal'));
            treee=figure('Position', get(0, 'Screensize'));
            copyfly=fly;
            for i=1:nanimals
                tframe=insertObjectAnnotation(AID,'circle',Dros(i,:),labels{i});
                Tframe=imresize(tframe,.5);
                AcTuP=imresize(ACTUP,.5);
                prompt={'What is the ID of this Fly?'};
                dlgTitle='ID Input';
                dims = [1 35];
                definput ={num2str(i)};
                APE=imshowpair(Tframe,AcTuP,'montage');
                answ = inputdlg(prompt,dlgTitle,dims,definput);
                opts.WindowStyle = 'normal';
                ID=str2double(answ{1});
                fly(ID).position=copyfly(i).position;
                fly(ID).observed=copyfly(i).observed;
                fly(ID).kalman=copyfly(i).kalman;
                prob(ID)=100;
                clc
                a=find(CE==fly(ID).position,1);
                labels{ID}=['Fly #', num2str(ID),' ', num2str(prob(ID)),'%'];
                Dros(ID,:)=pos(a,:);
            end
            close(treee)
            Brame=insertObjectAnnotation(Brame,'circle',Dros,labels,'TextBoxOpacity',0.45,'FontSize',12);
        end
        br=[];
    end
    writeVideo(VID,Brame);
    videoPlayer(Brame);
    
    ACTUP=Brame;
    %Update Kalman filters
    for i=1:nanimals
        bol=zeros(1,8);
        sol=[.5 .5 1 1 2.5 2.5 5 5];
        WW=mvnrnd(bol,sol);
        WW=transpose(WW);
        if isempty(QQ)
            QQ=0;
        end
        if ~ismember(i,JcS)
            X_t=Tmat*fly(i).kalman+WW;
            P_t=Tmat*fly(i).cov*transpose(Tmat)+QQ;
            K=P_t/(P_t+QQ);
            erx=sqrt((X_t(1)-fly(i).observed(1))^2);
            ery=sqrt((X_t(2)-fly(i).observed(2))^2);
            ervx=sqrt((X_t(3)-fly(i).observed(3))^2);
            ervy=sqrt((X_t(4)-fly(i).observed(4))^2);
            erax=sqrt((X_t(5)-fly(i).observed(5))^2);
            eray=sqrt((X_t(6)-fly(i).observed(6))^2);
            fly(i).kalman=X_t+K*(fly(i).observed-X_t);
            fly(i).cov=P_t-K*P_t;
            ABC=[erx ery ervx ervy erax eray];
            fly(i).error=cat(1,fly(i).error,ABC);
             X(i,count)=fly(i).kalman(1);
           Y(i,count)=fly(i).kalman(2);
        fly(i).position=[fly(i).kalman(1) fly(i).kalman(2)];
        fly(i).velocity=[fly(i).kalman(3) fly(i).kalman(4)];
        VX(i,count)=fly(i).kalman(3);
        VY(i,count)=fly(i).kalman(4);
        AX(i,count)=fly(i).kalman(5);
        AY(i,count)=fly(i).kalman(6);
        JX(i,count)=fly(i).kalman(7);
        JY(i,count)=fly(i).kalman(8);
        fly(i).acceleration=[fly(i).kalman(5) fly(i).kalman(6)];
        fly(i).jerk=[fly(i).kalman(7) fly(i).kalman(8)];
        
        end
       
    end
    count=count+1;
    if ~exist('amen','var')
        amen=prob;
    else
        amen=cat(1,prob,amen);
    end
end
mamen=mean(mean(amen));
release(videoPlayer);
close(VID)
%% Analysis Section
uiwait(msgbox(['The estimated accuracy of this tracking is ', num2str(mamen),'%']))
prompt={'Enter the number of regions being analyzed:'};
dlgTitle='Input';
dims = [1 35];
definput ={'18'};
answer = inputdlg(prompt,dlgTitle,dims,definput);
nregions=str2double(answer{1});
I=cell(1,nregions);
imshow(frame)
Regions=struct('Time',[],'In',[],'Xvar',[]);
for i=1:nanimals
    labels{i}=['Fly #', num2str(i)];
end
for i=1:nregions
    I{i}=drawpolygon;
end
st=find(Y(1,:)==0,1,'first');
if isempty(st)
    st=numframes;
end
VE=vmagsign(VX(:,1:st-1),VY(:,1:st-1));
for i=1:nregions
    R=inpolygon(X,Y,I{i}.Position(:,1),I{i}.Position(:,2));
    Regions(i).Time=sum(R,2);
    Regions(i).In= find(sum(R,2)>0);
    [r,l]=xcorr(VE(Regions(i).In(1),:),VE(Regions(i).In(2),:),WindowSize,'coeff');
    Regions(i).Xvar=max(r);
end
close all
traj=figure;
hold on
grid on
title('Trajectories')
for i=1:nanimals
    plot3(X(i,1:st-1),Y(i,1:st-1),time(1:st-1),'o-','LineWidth',2,'MarkerSize',10)
end
xlabel('X location (pixels)');
ylabel('Y location (pixels)');
zlabel('time (sec)');
legend(labels,'Location','best');
h=image(frame);
uistack(h,'bottom');
hold off
a=0; 
function v=vmagsign(V_x,V_y)
theta=atan(V_y./V_x);
v=V_y./cos(theta);
end