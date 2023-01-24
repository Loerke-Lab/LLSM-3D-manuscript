function [] = Figure2_T1examplesWithImages_LS(data,MovieNum,startInt)
%This function is used to generate T1 transition length plots and images of
%the cells

close all

% generate a blue-to-red colormap and then flip it
b2rcmap = b2r(-6,6);
b2rcmap = flip(b2rcmap,1);

% length scale factor (microns per pixel)
lsf = 0.1;

% font size
fs = 16;

% Move to directory with the data
cd(data(MovieNum).Source)
cd ..

% pull seconds per frame from data structure
spf = data(MovieNum).SecPerFrame;

% z-vector and z-step
zvec = 1:60;
zstep = 0.2635;
zvecMicrons = zvec*zstep;
Nlayers = numel(zvec);

% color map for each z-layer for plotting interface lengths
cmap = jet(Nlayers);


%load trackingMatrixZT_NodeFix that contains interface lengths
Tmatrix = load('trackingMatrixZT_NodeFix.mat').trackingMatrixZT_NodeFix;
Lengths3D = Tmatrix(:,1:8:end,zvec)*lsf;
[~,Nframes,Nlayers] = size(Lengths3D);

% load the indices of the T1 transitions
loadData = load('typeT1IntsV2.mat');
typeT1Ints = loadData.typeT1Ints;
typeT1Cells= loadData.typeT1Cells;
Nints = size(typeT1Ints,1);

    
% loop over all T1 transitions starting with 'startInt'
for int = startInt:Nints

        % form a Length matrix with both T1 and T3 interface lengths
        Length_T1 = squeeze(Lengths3D(typeT1Ints(int,1),:,:))';
        Length_T3 = squeeze(Lengths3D(typeT1Ints(int,2),:,:))';
        Length_T1andT3 = NaN(Nlayers,Nframes);
        finiteT1 = isfinite(Length_T1);
        finiteT3 = isfinite(Length_T3);
        Length_T1andT3(finiteT1) = Length_T1(finiteT1);
        Length_T1andT3(finiteT3) = -Length_T3(finiteT3);
        
        % find starting lengths and ending lengths and crop in time
        finite = isfinite(Length_T1andT3);
        startTime = find(any(finite,1),1,'first');
        endTime   = find(any(finite,1),1,'last');
        tvec = startTime:endTime;
        Length_T1andT3 = Length_T1andT3(:,tvec);

        
        

        % interpolate NaNs and if interpolated length is less than threshhold convert
        % to zeros
        Linterp =inpaint_nans(Length_T1andT3,2);
        Lnans = isnan(Length_T1andT3);
        Lthresh = 5*lsf; % threshhold is 5 pixels
        T2nans = Lnans & (Linterp < Lthresh);
        Linterp(T2nans) = 0;

        % smooth out the length matrix with Gaussian filter
        LinterpF = filterImage3DpaddedEdges(Linterp, 'Gauss', 2);
        maxL = ceil(max(abs(LinterpF(:))));
        
        Lmedian = median(LinterpF);
        
        fmedzero = find(Lmedian<=0,1,'first');
        if isempty(fmedzero)
            continue;
        end
        
        tvecMin = (tvec - tvec(fmedzero))*spf/60;
    
        
        figure(1)
        for j=1:Nlayers
            plot(tvecMin,Linterp(j,:),'color',cmap(j,:),'LineWidth',2)
            hold on
        end
        hold off
        grid on
        ylabel('Length (microns)','FontSize',16)
        xlabel('Time (min)','FontSize',16)
        title(num2str(int))
        set(gca,'FontSize',16)
        set(gca,'Ylim',[-10 10])
        
        figure(2)
        for j=1:Nlayers
            plot(tvec,Linterp(j,:),'color',cmap(j,:),'LineWidth',2)
            hold on
        end
        hold off
        grid on
        ylabel('Length (microns)','FontSize',16)
        xlabel('Time (F)','FontSize',16)
        title(num2str(int))
        set(gca,'FontSize',16)
        set(gca,'Ylim',[-10 10])
                
        choice = menu('Select Option','Select time point','Skip','Quit');

        switch choice
            case 1 % continue and keep data
                t = input('Pick a time frame: ');
                for z= [1 20 40 60]
                    clf
                    Iout=VisualizeCells(data,MovieNum,typeT1Cells(int,1:2),typeT1Cells(int,3:4),t,z,spf);
  
                    filename = strcat('~/Desktop/SampleVideo/Zlayer',num2str(z),'.tif');

                    imwrite(Iout,filename);
                end
                
            case 2
            case 3
                break

        end % switch

end % for
end % function

function [Icrop] = VisualizeCells(data,MovieNum,CLdis,CLapp,tvec,zlayer,spf)
% View cells for disappearing interface, and appearing interface

od = cd;


% generate color map for uniquely coloring each interface
colorlist = [0 0 1;1 0 0];
%h1 = figure;
for t=tvec
    
    % get image for current time frame
    cd(data(MovieNum).Source)
    image = tiffreadVolume(data(MovieNum).ImageFileListGap43{t});
    image = im2double(image(:,:,zlayer));
    image = imadjust(image);

    image3 = cat(3,image,image,image);

    [sx, sy] = size(image);

    % get tracking label matrix
    cd('SegmentationData/')
    trackmatrixpath = strcat(sprintf('frame%04d',t),'/ImageBWlabel.mat');
    loadtrackmatrix = load(trackmatrixpath);
    Tmatrix = loadtrackmatrix.ImageBWlabel;
    if ndims(Tmatrix)==3;Tmatrix = Tmatrix(:,:,zlayer);end

    RGBcells1 = zeros(sx,sy);
    RGBcells2 = zeros(sx,sy);
    RGBcells3 = zeros(sx,sy);
    % create color image of cells
    BWtrackcells = zeros(sx,sy);

    idx = Tmatrix == CLdis(1) | Tmatrix== CLdis(2);
    RGBcells1(idx) = colorlist(1,1);
    RGBcells2(idx) = colorlist(1,2);
    RGBcells3(idx) = colorlist(1,3);
    BWtrackcells(idx) = 1;

    idx = Tmatrix == CLapp(1) | Tmatrix == CLapp(2);
    RGBcells1(idx) = colorlist(2,1);
    RGBcells2(idx) = colorlist(2,2);
    RGBcells3(idx) = colorlist(2,3);
    BWtrackcells(idx) = 1;



    RGBcells = cat(3,RGBcells1,RGBcells2,RGBcells3);
    % transpose the color segmented image on the grayscale image
    imageRGB = 0.75*image3 + 0.25*RGBcells;


    % find centroid of the collection of the cells associated with the
    % interfaces. First, get BW of cells.
    stats = regionprops(BWtrackcells,'Centroid','BoundingBox');
    Centroid = stats.Centroid;
    % bounding box gives the smallest rectangle that holds all the cells. Since
    % we would like a square window we use the larger of the two dimensions and
    % enlarge it by 50%
    BoundingBox = stats.BoundingBox;

    % cropsize  
    %cropsize = round(2*max(BoundingBox(3:4)));
    
    % update cropsize  
    %if t == tvec(1)
        cropsize = 300;%round(2*max(BoundingBox(3:4)));
    %end

    % crop image
    rect = [round(Centroid(1) - cropsize/2) round(Centroid(2) - cropsize/2) cropsize cropsize];
    Icrop = imcrop(imageRGB,rect);
    Icrop = imrotate(Icrop,90);

    % show image
    imshow(Icrop);
    title(['T= ',num2str(t*spf/60),', Z=',num2str(zlayer*.5)],'FontSize',16)
    %figure(h1)
    pause(0.01)
end
% change back to original directory
cd(od)

end