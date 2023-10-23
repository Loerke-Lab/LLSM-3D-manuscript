function [autocorrAllcell] = AutocorrArea_LS(data,MovieNums,zvec)
% This function performs the area rate auto-correlation for each cell over
% all available movies
% INPUTS:
%           data: structure that contains Source, ImageFileList, and
%                   SecPerFrame fields for each movie
%           MovieNums: vector of movie indices (of data) to include in
%                      analysis
%           zvec: vector of z-layers to include. For the manuscript I used
%                   1:65.
%
%OUTPUTS:
%           autocorrAllCells: MxNxP matrix containing the 2D autocorrelations
%           of each cell. M is the number of Z-lags, N is the number of
%           Time-lags, and P is the number of cells.


close all;
fs = 16; % font size
XYsf = 0.1;% length scale factor for XY
Zsf  = 0.263; % length scale factor for Z
spf = data(MovieNums(1)).SecPerFrame; % seconds per frame
titleStr = strcat('Z=',num2str((zvec(1)-1)*Zsf),':',num2str((zvec(end)-1)*Zsf));

Nz = length(zvec); % number of z-layers

DeltaT = round(45/spf); % delta T for computing rate
maxlag_Z = round(Nz*0.75); % maximum lag value
maxlag_T = round(3*60/spf);

lagvecZ = -maxlag_Z:maxlag_Z; % lag vector in Z
lagzeroZ = find(lagvecZ == 0);
lagvecZMicrons = lagvecZ*Zsf;

lagvecT = -maxlag_T:maxlag_T; % lag vector in T
lagzeroT = find(lagvecT==0);
lagvecTMin = lagvecT*spf/60;


% initialize empty variable for storing autoCorr results
autocorrAllcell = [];

Nmovies = numel(MovieNums);
for mn=1:Nmovies

    %load data
    cd(data(MovieNums(mn)).Source)
    cd ..
    loaddata = load('trackingMatrixZT_NodeFix.mat');
    Tmatrix = loaddata.trackingMatrixZT_NodeFix(:,:,zvec);


    loadGdata = load('GeometryData.mat');
    Area3D = loadGdata.Area(:,:,zvec)*XYsf^2;
    CellNumbers = loadGdata.CellNumbers;
    Ncells = size(Area3D,1);


    for icell = 1:Ncells
        
        AreaZvT = squeeze(Area3D(icell,:,:))';
        
        
        [~, tvec] = FiniteSequenceFromNans(mean(AreaZvT));%,1);
        
        if length(tvec)*spf/60 < 3
            continue
        end
        
        AreaZvT = AreaZvT(:,tvec);
        
        % apply a little filtering to the area matrix
        AreaZvT = filterImage3DpaddedEdges(AreaZvT, 'Gauss', 2);

%         % remove all time points at the ends in which all layers are NaN.
%         frameStart = find(all(isfinite(AreaZvT),1),1,'first');
%         frameEnd = find(all(isfinite(AreaZvT),1),1,'last');
%         frameVec = frameStart:frameEnd;
%         if numel(frameVec)< range_T*2
%             continue;
%         end
%         AreaZvT = AreaZvT(:,frameStart:frameEnd);


        AreaZvT2 = AreaZvT(:,(DeltaT+1):end);
        AreaZvT1 = AreaZvT(:,1:(end-DeltaT));
        AreaRate = (AreaZvT2 - AreaZvT1)/DeltaT;

        % normalize rate
        AreaRateNorm= (AreaRate- nanmean(AreaRate(:)))/nanstd(AreaRate(:));
        
        
        %imagesc(AreaRateNorm)



        % perform auto-correlation and store result in in matrix
        [ autocorr1cell ] = autocorrelation2D(AreaRateNorm ,maxlag_Z, maxlag_T);
        autocorrAllcell = cat(3,autocorrAllcell,autocorr1cell);

        
        
%         subplot(1,2,1)  % Heatmap of Cell Area
%         fvec = 1:size(AreaZvT,1);
%         imagesc(fvec*spf/60,zvec*0.5,AreaZvT)
%         xlabel('Time(min)','FontSize',fs)
%         ylabel('Depth (microns)','FontSize',fs)
%         title('Area','FontSize',fs)
%         colorbar;
        
        
        %subplot(1,2,2) % Heatmap of Area Rate
%         fvec = 1:size(AreaDiff,2);
%         [AreaDiffFilt] = filterImage3DpaddedEdges(AreaDiff, 'Gauss', 2);
%         imagesc(fvec*spf/60,zvec*0.5,AreaDiffFilt)
%         xlabel('Time(min)','FontSize',fs)
%         ylabel('Depth (microns)','FontSize',fs)
%         title('Area Rate','FontSize',fs)
%         cmax = prctile(abs(AreaDiffFilt(:)),90);
%         b2rmap = b2r(-cmax,cmax);
%         colorbar;
%         colormap(b2rmap);
%         caxis([-cmax cmax])
%         hold on
%         
%         [FX,FY] = gradient(AreaDiffFilt);
%         [X,Y] = meshgrid(fvec,zvec);
%         quiver(X(1:3:end,1:3:end)*spf/60,Y(1:3:end,1:3:end)*0.5,FY(1:3:end,1:3:end),-FX(1:3:end,1:3:end))
        
%         subplot(2,3,3) % Heatmap of Autocorrelation
%         imagesc(colvecMin,rowvecMicrons,ac1cell)
%         xlabel('Time Lag (min)','FontSize',fs)
%         ylabel('Depth Lag (microns)','FontSize',fs)
%         title('2D Autocorrelation','FontSize',fs)
%         colorbar;
%        
%         
%         subplot(2,3,4)
%         
%         
%         imagesc(fvec*spf/60,zvec*0.5,-FY./FX)
%         xlabel('Time(min)','FontSize',fs)
%         ylabel('Depth (microns)','FontSize',fs)
%         title('-FY/FX','FontSize',fs)
%         
%         mat = -FY./FX;
%         cmax = prctile(abs(mat(:)),80);
%         b2rmap = b2r(-cmax,cmax);
%         colorbar;
%         colormap(b2rmap);
%         caxis([-cmax cmax])
       
        
%         subplot(2,3,5)
%         mag = sqrt(FX.^2+FY.^2);
%         imagesc(fvec*spf/60,zvec*0.5,mag)
%         xlabel('Time(min)','FontSize',fs)
%         ylabel('Depth (microns)','FontSize',fs)
%         title('Magnitude: sqrt(FX^2 + FY^2)','FontSize',fs)
%         colorbar;
%         colormap(jet)
%         
%         
%         subplot(2,3,6)
%         Mat = -FY./FX;
%         Thresh = prctile(mag(:),40);
%         mask = mag < Thresh;
%         Mat(mask) = NaN;
%         imagesc(fvec*spf/60,zvec*0.5,Mat)
%         xlabel('Time(min)','FontSize',fs)
%         ylabel('Depth (microns)','FontSize',fs)
%         title('-FY/Fx (with low Mag NaN''d)','FontSize',fs)
%         colorbar;
%         colormap(b2rmap);
%         caxis([-cmax cmax])
       % pause;clf;
        
        

         






    end % cell
    
    
end % movie

figure(1);
FullMatMean = nanmean(autocorrAllcell,3);
FullMatStd = nanstd(autocorrAllcell,[],3);
subplot(1,2,1)
imagesc(lagvecTMin,lagvecZMicrons,FullMatMean)
xlabel('Time Lag (min)','FontSize',fs)
ylabel('Depth Lag (microns)','FontSize',fs)
colorbar
cellcount = size(autocorrAllcell,3);
title([titleStr,' Auto-correlation, Ncells=',num2str(cellcount)],'FontSize',fs)

subplot(1,2,2)
plot(lagvecZMicrons,FullMatMean(:,lagzeroT),'LineWidth',2,'Color','b')
xlabel('Depth Lag (microns)','FontSize',fs,'Color','b')
set(gca,'Ylim',[-0.5 1])
hAx(1)=gca;
hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none');
hold(hAx(2),'on')
plot(hAx(2),lagvecTMin,FullMatMean(lagzeroZ,:),'LineWidth',2,'Color','g')
xlabel(hAx(2),'Time Lag (min)','FontSize',fs,'Color','g')
set(gca,'Ylim',[-0.5 1])
grid on


figure(2);
FullMatStd = nanstd(autocorrAllcell,[],3);
plot(lagvecZMicrons,FullMatMean(:,lagzeroT),'LineWidth',2,'Color','b')
xlabel('Depth Lag (microns)','FontSize',fs,'Color','b')
grid on
err = FullMatStd(:,lagzeroT);
y = FullMatMean(:,lagzeroT);
hold on
%errorenvelope(lagvecZMicrons,y,err,'b',0.2);
errorenvelope(lagvecZMicrons,y',err','b',0.2);
ylabel('Correlation Coeff.','FontSize',16)
%alpha(1)
print('1DautocorrAreaRate.eps','-depsc')

% subplot(1,3,3)
% 
% lagzeroZ = find(lagvecZ == 0);
% lagmaxZ  = find(lagvecZ == range_Z);
% lagminZ  = find(lagvecZ == -range_Z);
% 
% plot(lagvecTMin,FullMatMean(lagmaxZ,:),'LineWidth',2,'Color','b');hold on
% plot(lagvecTMin,FullMatMean(lagzeroZ,:),'LineWidth',2,'Color','r')
% plot(lagvecTMin,FullMatMean(lagminZ,:),'LineWidth',2,'Color','g');hold off
% xlabel('Time Lag (min)','FontSize',fs)
% set(gca,'Ylim',[-0.5 1])
% legend('max Z-lag','0 Z-lag','min Z-lag')
% grid on

%[psor,lsor] = findpeaks(PeakSig,x,'SortStr','descend');



cd('~/Desktop/')

end

