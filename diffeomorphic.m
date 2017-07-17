function diffeomorphic()
%===========================================================================
%Create diffewarped images.
%Number of morphed steps (images) is set by 'nsteps' with the amout of morphing held constant between images.
%Amount of morphing is set by 'maxdistortion'
%Each morphed images is saved as a jpeg.
%In figure (11), each morphed images is also positioned along the outline of a circle
%Please reference: Stojanoski, B., & Cusack, R. (2014). Time to wave good-bye to phase scrambling: Creating controlled scrambled images using
%diffeomorphic transformations. Journal of Vision, 14(12), 6. doi:10.1167/14.12.6

%Note: Mturk perceptual ratings of images are based on maxdistortion = 80; and nsteps = 20

%Rhodri Cusack and Bobby Stojanoski July 2013
%===========================================================================

maxdistortion=160; % changes max amount of distortion
nsteps=40; % number of steps
imsz= 1000; % size of output images (bigger or equal to 2 x input image)

picpath='directory with images'
outpicpath='output directory'
%create output directory
if ~exist(outpicpath, 'dir')
    mkdir(outpicpath);
end

imgtype='jpg'; %file type
fns=dir(fullfile(picpath,sprintf('*.%s',imgtype)));
figure(10);

[YI XI]=meshgrid(1:imsz,1:imsz);

if (~exist('Xn','var'))
end;


phaseoffset=floor(rand(1)*40);

for i=1:length(fns) %This is the number of objects in the directory
    
    Im=uint8(ones(imsz,imsz,3)*256); %128 for grey bckground; 256 for white
    P=imread(fullfile(picpath,fns(i).name));
    P=P(:,:,1:3);
    Psz=size(P);
    
    % Upsample by factor of 2 in two dimensions
    P2=zeros([2*Psz(1:2),Psz(3)]);
    P2(1:2:end,1:2:end,:)=P;
    P2(2:2:end,1:2:end,:)=P;
    P2(2:2:end,2:2:end,:)=P;
    P2(1:2:end,2:2:end,:)=P;
    P=P2;
    Psz=size(P);
    
    % Pad image if necessary
    x1=round((imsz-Psz(1))/2);
    y1=round((imsz-Psz(2))/2);
    
    % Add fourth plane if necessary
    if (Psz(3)==4)
        Im(:,:,4)=0;
    end;
    Im((x1+1):(x1+Psz(1)),(y1+1):(y1+Psz(2)),:)=P;
    
    [cxA cyA]=getdiffeo(imsz,maxdistortion,nsteps);
    [cxB cyB]=getdiffeo(imsz,maxdistortion,nsteps);
    [cxF cyF]=getdiffeo(imsz,maxdistortion,nsteps);
    
    interpIm=Im;
    figure(11);
    clf
    for quadrant=1:4
        switch (quadrant)
            case 1
                cx=cxA;
                cy=cyA;
                ind=1;
                indstep=1;
            case 2
                cx=cxF-cxA;
                cy=cyF-cyA;
            case 3
                ind=4*nsteps;
                indstep=-1;
                interpIm=Im;
                cx=cxB;
                cy=cyB;
            case 4
                cx=cxF-cxB;
                cy=cyF-cyB;
        end
        cy=YI+cy;
        cx=XI+cx;
        mask=(cx<1) | (cx>imsz) | (cy<1) | (cy>imsz) ;
        cx(mask)=1;
        cy(mask)=1;
        figure(10);
        subplot (4,2,quadrant*2-1)
        imagesc(cx)
        subplot (4,2,quadrant*2)
        imagesc(cy)
        w=0.1;
        for j=1:nsteps %This is the number of steps - Total number of warps is nsteps * quadrant
            centrex=0.5+(0.5-w/2)*cos((phaseoffset+ind)*2*pi/(4*nsteps));
            centrey=0.5+(0.5-w/2)*sin((phaseoffset+ind)*2*pi/(4*nsteps));
            figure(11);
            if (mod(ind,2)==0)
                axes('position',[centrex-w/2 centrey-w/2 w w]);
                imagesc(interpIm(:,:,1:3));
                axis off
            end;
            [pth randfn] = fileparts(tempname);
            randstr = randfn(27:end);
            imwrite(uint8(interpIm),fullfile(outpicpath,sprintf('Im_%02d_%02d.%s',i,ind,imgtype)),imgtype);
            randfn = sprintf('W_%02d_%02d_%s.%s',i,ind,randstr,imgtype); %This is the recoded name of the file
            interpIm(:,:,1)=interp2(double(interpIm(:,:,1)),cy,cx);
            interpIm(:,:,2)=interp2(double(interpIm(:,:,2)),cy,cx);
            interpIm(:,:,3)=interp2(double(interpIm(:,:,3)),cy,cx);
            filenames{i,ind} = randfn;
            origfn{i} = fns(i).name;
            save filename_mapping.mat filenames origfn
            ind=ind+indstep;
        end;
    end;
end
end


function [XIn YIn]=getdiffeo(imsz,maxdistortion,nsteps)
ncomp=6;

[YI XI]=meshgrid(1:imsz,1:imsz);

% make diffeomorphic warp field by adding random DCTs
ph=rand(ncomp,ncomp,4)*2*pi;
a=rand(ncomp,ncomp)*2*pi;
Xn=zeros(imsz,imsz);
Yn=zeros(imsz,imsz);
for xc=1:ncomp
    for yc=1:ncomp
        Xn=Xn+a(xc,yc)*cos(xc*XI/imsz*2*pi+ph(xc,yc,1))*cos(yc*YI/imsz*2*pi+ph(xc,yc,2));
        Yn=Yn+a(xc,yc)*cos(xc*XI/imsz*2*pi+ph(xc,yc,3))*cos(yc*YI/imsz*2*pi+ph(xc,yc,4));
    end;
end;
% Normalise to RMS of warps in each direction
Xn=Xn/sqrt(mean(Xn(:).*Xn(:)));
Yn=Yn/sqrt(mean(Yn(:).*Yn(:)));

YIn=maxdistortion*Yn/nsteps;
XIn=maxdistortion*Xn/nsteps;

end
