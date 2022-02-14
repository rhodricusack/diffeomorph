function diffeomorphic_movie_warpshift()
%===========================================================================
% Create diffeowarped movies.
%  Applies a gradually drifting warp field to each frame, achieved by the 
%  phase of the components taking a random walk.
%  Will operate on all .mov files in the directory defined by "moviepath"
%Please reference: Stojanoski, B., & Cusack, R (2013). Time to wave goodbye to phase scrambling � creating unrecognizable control stimuli using a diffeomorphic transform.  Abstract Vision Science Society

% =======================
% Rhodri Cusack and Bobby Stojanoski July 2013
% Rhodri Cusack 2020-10-12: adapted to allow gently drifting warp fields - phases take random walk

maxdistortion=60; % changes max amount of distortion
nsteps=20; % number of steps
ncomp=10;
%imsz= 1000; % size of output images (bigger or equal to 2 x input image)
phasedrift=pi/8; % phase drift by randn(0,phasedrift) each second

moviepath='/home/clionaodoherty/diffeomorph/towarp'; 
outvidpath='/home/clionaodoherty/diffeomorph/warped';

% Read in .mov
fns=dir(fullfile(moviepath,'*.mov'));
figure(10);


imsz=1280+4*maxdistortion;    

% Generate disortion field for all frames - try keeping it fixed to start
% Only need one not 3 distortion fields, as no longer need the
% continuous circle
tic
[cx, cy, a, b, ph]=getdiffeo(imsz,maxdistortion,nsteps,ncomp);
[cx, cy]=postprocess_diffeo(imsz,cx,cy);

% In one second well be here
% Alter phase
ph=ph+phasedrift*randn(ncomp,ncomp,4);

[nextcx, nextcy, a, b, ph]=getdiffeo(imsz,maxdistortion,nsteps,ncomp,a,b,ph);
[nextcx, nextcy]=postprocess_diffeo(imsz,nextcx,nextcy);

fprintf('Time to create distortion field: %f s\n',toc)




for i=1:length(fns) %This is the number of objects in the directory
    M=VideoReader(fullfile(moviepath,fns(i).name));
    z=0;
    v=VideoWriter((fullfile(outvidpath,sprintf('Vid_%02d',i)))); 
    outframenum=0;
    v.FrameRate=M.FrameRate; 
    open(v);
    
    oldsecs=0;
    while M.hasFrame() %to pass over
        z=z+1;
        outframenum=outframenum+1;
        P=M.readFrame();
        P=P(:,:,1:3);
        Psz=size(P);
        
        %moved below chunk from outside of the loop
        
        if (~exist('Xn','var'))
        end;
        
        Im=uint8(ones(imsz,imsz,3)*256); %128 for grey bckground and 256 for white
        
        % Pad image if necessary
        x1=round((imsz-Psz(1))/2);
        y1=round((imsz-Psz(2))/2);
        
        % Add fourth plane if necessary
        if (Psz(3)==4)
            Im(:,:,4)=0;
        end;
        Im((x1+1):(x1+Psz(1)),(y1+1):(y1+Psz(2)),:)=P;
        % Pad with mirrored extensions of image

        Im(1:x1,(y1+1):(y1+Psz(2)),:)=P(x1:-1:1,:,:);
        Im(x1+Psz(1)+1:end,(y1+1):(y1+Psz(2)),:)=P(end:-1:end-x1+1,:,:);

        Im(:,1:y1,:)=Im(:,(y1+1+y1):-1:(y1+2),:);
        Im(:,y1+Psz(2)+1:end,:)=Im(:,(y1+Psz(2):-1:Psz(2)+1),:);

        % Warp field is updated each second, by drifting the phase
        if phasedrift>0
            newsecs=floor(z/v.FrameRate);
            if newsecs~=oldsecs
                fprintf('Secs of movie processed: %d\n',newsecs)
                oldsecs=newsecs;
                % Alter phase
                ph=ph+phasedrift*randn(ncomp,ncomp,4);
                
                % Create new distortion field
                cx=nextcx;
                cy=nextcy;
                [nextcx, nextcy, a, b, ph]=getdiffeo(imsz,maxdistortion,nsteps,ncomp,a,b,ph);
                [nextcx, nextcy]=postprocess_diffeo(imsz,nextcx,nextcy);
            end
        end
        
        % Interpolate between now and next distortion field in 1 second
        f=z/v.FrameRate;
        f=f-floor(f);
        cyi=(1-f)*cy + f*nextcy;
        cxi=(1-f)*cx + f*nextcx;
        
        % Start off with undisorted image
        interpIm=Im;
        for j=1:nsteps %This is the number of steps - Total number of warps is nsteps * quadrant
            interpIm(:,:,1)=interp2(double(interpIm(:,:,1)),cyi,cxi);
            interpIm(:,:,2)=interp2(double(interpIm(:,:,2)),cyi,cxi);
            interpIm(:,:,3)=interp2(double(interpIm(:,:,3)),cyi,cxi);
        end;
        
        % Trim down again
        interpIm=interpIm((x1+1):(x1+Psz(1)),(y1+1):(y1+Psz(2)),:);
        movframe=uint8(interpIm);
        writeVideo(v,movframe);

    end
    
    close(v); %moved from below
  
end
end


% post process diffeo
function [cx,cy]=postprocess_diffeo(imsz,cx,cy)

[YI, XI]=meshgrid(1:imsz,1:imsz);
cy=YI+cy;
cx=XI+cx;
mask=(cx<1) | (cx>imsz) | (cy<1) | (cy>imsz) ;
cx(mask)=1;
cy(mask)=1;
end


% Create distortion field
% inputs
%  imsz - size of visual image
%  maxdistortion - max amount of displacement of any pixel after nsteps
%    applications of flow field 
%  nsteps - number of flowfield steps needed
%  ncomp - number of sin + cosine components along each axis
% outputs
%  XIn, YIn are imsz*imsz flow field arrays of displacement (in pixels) in X and Y direction
%

%% Other function at bottom of the script
function [XIn, YIn, a, b, ph]=getdiffeo(imsz,maxdistortion,nsteps,ncomp,a,b,ph)
if ~exist('ncomp','var')
    ncomp=6;
end;

[YI, XI]=meshgrid(1:imsz,1:imsz);

% make diffeomorphic warp field by adding random DCTs
if ~exist('ph','var')
	ph=rand(ncomp,ncomp,4)*2*pi;
end;
if ~exist('a','var')
	a=rand(ncomp,ncomp)*2*pi;
end;
if ~exist('b','var')
	b=rand(ncomp,ncomp)*2*pi; % different amplitudes for x and y DCT components
end;

Xn=zeros(imsz,imsz);
Yn=zeros(imsz,imsz);
for xc=1:ncomp
    for yc=1:ncomp
        Xn=Xn+a(xc,yc)*cos(xc*XI/imsz*2*pi+ph(xc,yc,1))*cos(yc*YI/imsz*2*pi+ph(xc,yc,2));
        Yn=Yn+b(xc,yc)*cos(xc*XI/imsz*2*pi+ph(xc,yc,3))*cos(yc*YI/imsz*2*pi+ph(xc,yc,4));
    end;
end;
% Normalise to RMS of warps in each direction
Xn=Xn/sqrt(mean(Xn(:).*Xn(:)));
Yn=Yn/sqrt(mean(Yn(:).*Yn(:)));

YIn=maxdistortion*Yn/nsteps;
XIn=maxdistortion*Xn/nsteps;

end
