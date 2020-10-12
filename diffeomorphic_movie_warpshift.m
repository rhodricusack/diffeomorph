function BS_RC_diffeomorphic_movie_v3()
%===========================================================================
%Create diffewarped images.
%Number of morphed steps (images) is set by 'nsteps' with the amout of morphing held constant between images.
%Amount of morphing is set by 'maxdistortion'
%Each morphed images is saved as a jpeg.
%In figure (11), each morphed images is also positioned along the outline of a circle
%Please reference: Stojanoski, B., & Cusack, R (213). Time to wave goodbye to phase scrambling – creating unrecognizable control stimuli using a diffeomorphic transform.  Abstract Vision Science Society
%Note: Mturk perceptual ratings of images are based on maxdistortion = 80; and nsteps = 20

% Rhodri Cusack and Bobby Stojanoski July 2013
%===========================================================================
% Rhodri Cusack 2020-10-12: adapted to allow gently drifting warp fields - phases take random walk

cd 'G:\Dropbox\MovieWarpingScript\'

maxdistortion=60; % changes max amount of distortion
nsteps=20; % number of steps
ncomp=10;
%imsz= 1000; % size of output images (bigger or equal to 2 x input image)
phasedrift=0.1; % phase drift by randn(0,phasedrift) each second

picpath='G:\Dropbox\MovieWarpingScript\'; %'D:\Bobby\Experiments\fMRI\VSTM_hierarchy\Images\AllImages';
outpicpath='G:\Dropbox\MovieWarpingScript\WarpedMovie'; %'D:\Bobby\Experiments\fMRI\VSTM_hierarchy\Images\WarpedImages1_80';
outvidpath='G:\Dropbox\MovieWarpingScript\WarpedMovie_Vid';
imgtype='png'; % stills for now, you'll need to use VideoWriter to write video  

% But read in .mov
fns=dir(fullfile(picpath,'*.mov'));
figure(10);


imsz=1280+4*maxdistortion;    

% Generate disortion field for all frames - try keeping it fixed to start
% Only need one not 3 distortion fields, as no longer need the
% continuous circle
[cx, cy, a, b, ph]=getdiffeo(imsz,maxdistortion,nsteps,ncomp);
[YI, XI]=meshgrid(1:imsz,1:imsz);
cy=YI+cy;
cx=XI+cx;
mask=(cx<1) | (cx>imsz) | (cy<1) | (cy>imsz) ;
cx(mask)=1;
cy(mask)=1;



for i=1:length(fns) %This is the number of objects in the directory
    M=VideoReader(fullfile(picpath,fns(i).name));
    z=0;
    v=VideoWriter((fullfile(outvidpath,sprintf('Vid_%02d',i)))); 
    outframenum=0;
    v.FrameRate=M.FrameRate; 
    open(v);
    
    lastsecs=0;
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

        
        % Start off with undisorted image
        interpIm=Im;
        for j=1:nsteps %This is the number of steps - Total number of warps is nsteps * quadrant
            interpIm(:,:,1)=interp2(double(interpIm(:,:,1)),cy,cx);
            interpIm(:,:,2)=interp2(double(interpIm(:,:,2)),cy,cx);
            interpIm(:,:,3)=interp2(double(interpIm(:,:,3)),cy,cx);
        end;
        
        % Trim down again
        interpIm=interpIm((x1+1):(x1+Psz(1)),(y1+1):(y1+Psz(2)),:);
        %imwrite(uint8(interpIm),fullfile(outpicpath,sprintf('Im_%02d_%02d.%s',i,z,imgtype)),imgtype);
        movframe=uint8(interpIm);
        writeVideo(v,movframe);

	% Drift phase each second
	if phasedrift>0
		newsecs=floor(z/V.FrameRate)
		if newsecs~=oldsecs
			oldsecs=newsecs
			ph=ph+phasedrift*randn(ncomp,ncomp,4);
			[cx, cy, a, b, ph]=getdiffeo(imsz,maxdistortion,nsteps,ncomp,a,b,ph);
		end;
	end;
    end;
    
    close(v); %moved from below
    
    
    
    %         end;
    %     end;
    
    %I think this is where we will want to write out the movie before we
    %move onto the next movie clip
    %M=VideoReader(fullfile(picpath,fns(i).name));
    %outframenum=0;
    %v=VideoWriter((fullfile(outvidpath,sprintf('Vid_%02d',i)))); 
    %open(v);
    %v.FrameRate=M.FrameRate; 
    %while M.hasFrame() %to load into the video
        %outframenum=outframenum+1
        %movframe=imread(fullfile(outpicpath,sprintf('Im_%02d_%02d.%s',i,outframenum,imgtype)));
        %writeVideo(v,movframe);
    %end;
    %close(v);
end
end



% Distortion field
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
