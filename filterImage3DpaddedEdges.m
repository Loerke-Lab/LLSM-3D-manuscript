function [ImageMatrixFiltered] = filterImage3DpaddedEdges(image, filtername, filtersize)
% filter image in 3 dimension; pad edges to prevent boundary effects
% INPUT:    image       = 2D or 3D image
%           filtername  = filter method (see below) 
%           filtersize  = filter size (see below)
%          
%           OPTIONS:
%           filtername = 'Gauss' Gaussian filter
%           filtersize = sigma of Gauss, can be either 
%           [sigma] or [sigmax sigmay sigmaz];
%
%           filtername = 'ra' rolling average
%           filtersize = width of rolling average (+-ra) 
%           [width] or [widthx widthy widthz];
%
%           NOTE 1: If image is 2-dimensional, the filtering will be
%           performed two-dimensional, as well. If the image is
%           3-dimensional, then there are two options: If filtersize is of
%           length 1, then the filtering will be performed in a
%           two-dimensional way on successive layers of the image; if
%           filtersize has length 3 [sigx sigz sigz], the filtering will be
%           3-dimensional
%       
%           NOTE 2: The image matrix is padded before filtering to avoid 
%           edge effects; the padded pixel layer corresponds to filtersize 
%           in that dimension
%
%
% OUTPUT:   ImageMatrixFiltered: a matrix of the image with the specified
%           filter applied


%% determine filter size
if length(filtersize)==1
    fsx = filtersize; fsy = filtersize; fsz = 0;
    ftype = 2;
else
    fsx = filtersize(1); fsy = filtersize(2); fsz = filtersize(3);
    ftype = 3;
end

%% define filter

switch filtername
    
    % Rolling average
    case 'ra'
        fdimx = 2*fsx+1; 
        fdimy = 2*fsy+1; 
        fdimz = 2*fsz+1;
        
        cfilter     = ones(fdimx,fdimy,fdimz);
        cfilterNorm = cfilter/sum(cfilter(:));
       
    % Gaussian filter    
    case 'Gauss'
        sigmax = fsx; sigmay = fsy; sigmaz = fsz;
        if fsz==0, sigmaz = inf; end
        
        fsxfull = ceil(3*fsx); 
        fsyfull = ceil(3*fsy);  
        fszfull = ceil(3*fsz); 
        fdimx   = 2*fsxfull+1; 
        fdimy   = 2*fsyfull+1; 
        fdimz   = 2*fszfull+1;
        
        cfilter     = zeros(fdimx,fdimy,fdimz);
        for i = -fsxfull:fsxfull
            ex = exp( -(i^2) / (2*sigmax^2) );
            for j = -fsyfull:fsyfull
                ey = exp( -(j^2) / (2*sigmay^2) );
                for k= -fszfull:fszfull
                    ez = exp( -(k^2) / (2*sigmaz^2) );
                    cfilter(i+fsxfull+1, j+fsyfull+1, k+fszfull+1) = ex * ey * ez;
                end
            end
        end
        cfilterNorm     = cfilter/sum(cfilter(:)); 
        
        padx = (fdimx-1)/2; 
        pady = (fdimy-1)/2; 
        padz = (fdimz-1)/2; 
        
end


%% construct padded image

% size of image
[sx,sy,sz] = size(image);

% size of padded image
sxp = sx + 2*padx;
syp = sy + 2*pady;
szp = sz + 2*padz;

% initialize matrix for ImageMatrixPadded
ImageMatrixPadded = zeros(sxp,syp,szp);

% fill center
ImageMatrixPadded(padx+1:padx+sx,pady+1:pady+sy,padz+1:padz+sz) = image;

% pad x (rows) - pad top rows, then pad bottom rows
ImageMatrixPadded(1:padx,:,:) = repmat(ImageMatrixPadded(padx+1,:,:),[padx 1 1]);
ImageMatrixPadded(padx+sx+1:sxp,:,:) = repmat(ImageMatrixPadded(padx+sx,:,:),[padx 1 1]);

% pad y (columns) - pad left columns, then pad right columns
ImageMatrixPadded(:,1:pady,:) = repmat(ImageMatrixPadded(:,pady+1,:),[1 pady 1]);
ImageMatrixPadded(:,pady+sy+1:syp,:,:) = repmat(ImageMatrixPadded(:,pady+sy,:),[1 pady 1]);

% pad z (layers) - pad top layers, the pad bottom layers
ImageMatrixPadded(:,:,1:padz) = repmat(ImageMatrixPadded(:,:,padz+1),[1 1 padz]);
ImageMatrixPadded(:,:,padz+sz+1:szp) = repmat(ImageMatrixPadded(:,:,padz+sz),[1 1 padz]);


%% filter image

if (sz>1) & (ftype==2)
    for iz = 1:sz
        % filter layers individually
        cimage = ImageMatrixPadded(:,:,iz);
        cimage_filter = convn(cimage,cfilterNorm);
        ImageMatrixPaddedFiltered(:,:,iz) = cimage_filter;
    end
else
    % filtered image
    ImageMatrixPaddedFiltered = convn(ImageMatrixPadded,cfilterNorm); 
end

% crop out center
ImageMatrixFiltered = ImageMatrixPaddedFiltered(2*padx+1:2*padx+sx,2*pady+1:2*pady+sy,2*padz+1:2*padz+sz);

% figure; imshow(ImageMatrix(:,:,1),[]);
% figure; imshow(ImageMatrixFiltered(:,:,1),[]);

end

