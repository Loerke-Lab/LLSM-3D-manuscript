function [resultsPDfull,resultsPDapical,resultsPDbasal] = VolumeChange_LS_Rho_Wire(data,MovieNums)
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
XYsf = 0.163; % voxel size in X and Y
zsf  = 1; % voxel size in Z

spf = data(MovieNums(1)).SecPerFrame; % seconds per time point

zvecFull = 4:19; % full depth layer vector
zvecApical = 1:5; % apical depth layer vector
zvecBasal = 12:16; % basal depth layer vector
Nlayers = numel(zvecFull); % total number of layers

DeltaT = 1; % time difference (in minutes) to measure percent change
DeltaF = round(DeltaT*60/spf); % time difference in Frames

% initialize empty variables to store Percent Difference values
resultsPDfull = [];
resultsPDapical = [];
resultsPDbasal = [];

cellNumber = 0; % added 5/9/23 KB

Nmovies = numel(MovieNums);
for mn=1:Nmovies

    % move to the movie's source to load Geometry data
    cd(data(MovieNums(mn)).Source)
    %cd ..
    
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

        cellNumber = cellNumber+1; % added 5/9/23 KB
        
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
set(gca,'Ylim',[-5 25])
grid on

% print the following p-values to command window
[~,p1] = kstest2(resultsPDfull,resultsPDapical)
[~,p2] = kstest2(resultsPDfull,resultsPDbasal)


        
end %function

function [PD] = PercentDiff(Mat,DeltaF)
% finds the percent differences between every two points separated by
% DeltaF in the vector M

Mat2 = Mat(:,(DeltaF+1):end);
Mat1 = Mat(:,1:(end-DeltaF));
PD = 100*abs(Mat2 - Mat1)./((Mat1+Mat2)/2);


end
