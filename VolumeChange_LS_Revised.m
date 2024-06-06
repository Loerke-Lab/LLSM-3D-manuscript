function [resultsPDfull,resultsPDapical,resultsPDbasal] = VolumeChange_LS_Revised(data,MovieNums)
% this function calculates the percent differences in cell volume over 1
% minute time change for full depth, apical region, and basal region

% INPUTS:
%           data: structure that contains Source, ImageFileList, and
%                   SecPerFrame fields for each movie
%           MovieNums: vector of movie indices (of data) to include in
%                      analysis

%OUTPUTS: 
%           resultsPDfull/apical/basal: vector of cell volume percent differences for
%           each cell and time change. Each cell will contribute multiple
%           percent difference values depending on how long it is present.


%% Movie parameters
XYsf = 0.1; % voxel size in X and Y
zsf  = 0.2635; % voxel size in Z

spf = data(MovieNums(1)).SecPerFrame; % seconds per time point

zvecFull = 1:65; % full depth layer vector
zvecApical = 1:19; % apical depth layer vector
zvecBasal = 47:65; % basal depth layer vector
Nlayers = numel(zvecFull); % total number of layers

DeltaT = 1; % time difference (in minutes) to measure percent change
DeltaF = round(DeltaT*60/spf); % time difference in Frames

% initialize empty variables to store Percent Difference values
resultsPDfull = [];
resultsPDapical = [];
resultsPDbasal = [];

% initialize count for number of cells + data points
numberCells = 0;
numberDataPoints = 0;

Nmovies = numel(MovieNums);
for mn=1:Nmovies

    % move to the movie's source to load Geometry data
    cd(data(MovieNums(mn)).Source)
    cd ..
    
    % load 3D matrix of cell areas (rows=cells, columns=time frames, layers=z-depth)
    loadGdata = load('GeometryData.mat');
    Area3D = loadGdata.Area(:,:,zvecFull)*XYsf^2*zsf; % convert to microns^3
    [Ncells,Nframes,~] = size(Area3D);
    
    % initialize a 2D volume matrix to store cell volume (sum over depth)
    VolumeMatFull = NaN(Ncells,Nframes);
    VolumeMatApical = NaN(Ncells,Nframes);
    VolumeMatBasal = NaN(Ncells,Nframes);
    
    % loop over each cell to get volumes
    for icell = 1:Ncells
        
        % take a single cell slice and squeeze to matrix with rows=z_depth
        % and columns=time_frames
        AreaZvT = squeeze(Area3D(icell,:,:))';
        

        numFinite = sum(isfinite(AreaZvT));
        % if number of Finite layers is not Nlayers, set frame to NaN
        numFinite(numFinite<Nlayers) = NaN;
        [~, tvec] = FiniteSequenceFromNans(numFinite);%,1);
        
        VolumeMatFull(icell,tvec) = sum(AreaZvT(:,tvec),1);
        VolumeMatApical(icell,tvec) = sum(AreaZvT(zvecApical,tvec),1);
        VolumeMatBasal(icell,tvec) = sum(AreaZvT(zvecBasal,tvec),1);
        
    end % cell
    
    
    % if the cell has no finite volume values remove its row
    idxallnan = all(isnan(VolumeMatFull),2);
    VolumeMatFull(idxallnan,:) = [];
    idxallnan = all(isnan(VolumeMatApical),2);
    VolumeMatApical(idxallnan,:) = [];
    idxallnan = all(isnan(VolumeMatBasal),2);
    VolumeMatBasal(idxallnan,:) = [];
    
    % compute the percent differences over DeltaF frames
    [moviePDfull] = PercentDiff(VolumeMatFull,DeltaF);
    [moviePDapical] = PercentDiff(VolumeMatApical,DeltaF);
    [moviePDbasal] = PercentDiff(VolumeMatBasal,DeltaF);


    % store percent differences in a vector
    resultsPDfull = [resultsPDfull;moviePDfull(isfinite(moviePDfull))];
    resultsPDapical = [resultsPDapical;moviePDapical(isfinite(moviePDapical))];
    resultsPDbasal = [resultsPDbasal;moviePDbasal(isfinite(moviePDbasal))];


    %%% these lines give the N values (# of cells & # of data points)
    % if the cell percent difference has no finite volume values remove its row
    idxallnan = all(isnan(moviePDfull),2);
    moviePDfull(idxallnan,:) = [];
    idxallnan = all(isnan(moviePDapical),2);
    moviePDapical(idxallnan,:) = [];
    idxallnan = all(isnan(moviePDbasal),2);
    moviePDbasal(idxallnan,:) = [];

    numberCells = numberCells + size(moviePDfull,1);
    numberDataPoints = numberDataPoints + sum(isfinite(moviePDfull),'all');


    
end % movie


fs = 16; % font size for figure




figure(1)

n1 = numel(resultsPDfull);
n2 = numel(resultsPDapical);
n3 = numel(resultsPDbasal);
g1 = repmat({'Full Depth'},n1,1);
g2 = repmat({'Apical'},n2,1);
g3 = repmat({'Basal'},n3,1);
g = [g1; g2; g3];
x = [resultsPDfull;resultsPDapical;resultsPDbasal];

boxplot(x,g,'Symbol','')
ylabel('Percent Diff in Volume over 1 min','FontSize',fs)
set(gca,'Ylim',[-100 100])
yticks([-100 -75 -50 -25 0 25 50 75 100])
grid on

% print the following p-values to command window
[~,p1] = kstest2(resultsPDfull,resultsPDapical)
[~,p2] = kstest2(resultsPDfull,resultsPDbasal)

set(findobj(gca,'type','line'),'linew',2)
set(gca,'FontSize',20)


% plot histograms for the data
PDedges = -75:5:75;

figure
% histogram(resultsPDfull,PDedges,'facealpha',.5,'edgecolor','none','Normalization','probability')
% hold on
% histogram(resultsPDapical,PDedges,'facealpha',.5,'edgecolor','none','Normalization','probability')
% histogram(resultsPDbasal,PDedges,'facealpha',.5,'edgecolor','none','Normalization','probability')
histogram(resultsPDfull,PDedges,'facealpha',.5,'Normalization','probability')
hold on
histogram(resultsPDapical,PDedges,'facealpha',.5,'Normalization','probability')
histogram(resultsPDbasal,PDedges,'facealpha',.5,'Normalization','probability')
legend('Full Depth','Apical','Basal')
ylabel('Probability')
xlabel('Percent Diff in Volume over 1 min')
set(gca,'Fontsize',20)
grid on


figure

[valuesFull, edgesFull] = histcounts(resultsPDfull,PDedges, 'Normalization', 'probability');
centersFull = (edgesFull(1:end-1)+edgesFull(2:end))/2;

plot(centersFull, valuesFull, 'LineWidth',2);hold on

[valuesApical, edgesApical] = histcounts(resultsPDapical,PDedges, 'Normalization', 'probability');
centersApical = (edgesApical(1:end-1)+edgesApical(2:end))/2;

plot(centersApical, valuesApical, 'LineWidth',2)

[valuesBasal, edgesBasal] = histcounts(resultsPDbasal,PDedges, 'Normalization', 'probability');
centersBasal = (edgesBasal(1:end-1)+edgesBasal(2:end))/2;

plot(centersBasal, valuesBasal, 'LineWidth',2)
legend('Full Depth','Apical','Basal')
ylabel('Probability')
xlabel('Percent Diff in Volume over 1 min')
set(gca,'Fontsize',20)
grid on
ylim([0 0.25])
xlim([-100 100])
xticks([-100 -80 -60 -40 -20 0 20 40 60 80 100])




end %function




function [PD] = PercentDiff(Mat,DeltaF)
% finds the percent differences between every two points separated by
% DeltaF in the vector M

Mat2 = Mat(:,(DeltaF+1):end);
Mat1 = Mat(:,1:(end-DeltaF));
%PD = 100*abs(Mat2 - Mat1)./((Mat1+Mat2)/2);
PD = 100*(Mat2 - Mat1)./((Mat1+Mat2)/2);

end
