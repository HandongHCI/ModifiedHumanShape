clear
clc
startup
expidx = 0;
p = expParams(expidx);

%% manipulate shape and pose
% fprintf('mean pose and shape');
% poseParams = zeros(1,31); % mean pose
% shapeParams = zeros(1,p.nPCA); % mean shape
% load([p.modelInDir '/evectors'],'evectors'); % shape space eigenvectors
% load([p.modelInDir '/evalues'],'evalues'); % shape space eigenvalues
% evectors = evectors(1:p.nPCA,:);
% evalues = evalues(1:p.nPCA);
% points = shapepose(poseParams,shapeParams,evectors,p.modelInDir);
% 
% % show model
% load(p.facesSM,'faces');
% clf;
% showmodel(points,faces,'r',[],0);
% axis equal; view(45,22.5); 
% fprintf(', press any key\n');
% pause;
% 
% fprintf('change pose and shape');
% % description of pose parameters in shapemodel/poseParamsDescript.m
% poseParams(24) = 45/180*pi; % rotate right shoulder by 45 degrees
% shapeParams(1) = 3*sqrt(evalues(1)); % change height to 3 st.d. from mean
% points = shapepose(poseParams,shapeParams,evectors,p.modelInDir);
% 
% % show model
% hold on;
% showmodel(points,faces,'g',[],0);
% axis equal; view(45,22.5); 
% fprintf(', press any key\n');
% pause;
% 
% fprintf('visualize eigenvector scaled by 3 st.d.');
% vpAngle = 90; % visualization viewpoint
% idxShape = 1; % first eigenvector
% sign = -1; % sign of st.d.
% bSave = false; % save figure
% clf;
% visModelVP(expidx,vpAngle,idxShape,sign,bSave);
% fprintf(', press any key\n');
% pause;

%% register human scan
fprintf('[[register human scan]]\n\n');

files = dir('data/*.ply');
path = 'data/';
starting = 16;

for i = starting:length(files)
% for i = starting
    fprintf('[%04d/%04d]\n', i, length(files));
    
    scanName = files(i).name(1:4);
    scanFilenames = {[path scanName '.ply']};
    landmarkFilenames = {[path scanName '.asc']};
    fitMesh(scanFilenames,landmarkFilenames,expidx);

%     % visualize registration
%         fprintf('results: red - fitted mesh, blue - overlaid scan\n');
% 
%         fprintf('result of pose fitting using landmarks');
%         clf; visFitDir([p.fitDir '/' scanName],0);
%         fprintf(', press any key\n');
%         % pause;
% 
%         fprintf('result of pose and shape fitting using all vertices');
%         clf; visFitDir([p.fitDir '/' scanName],1);
%         fprintf(', press any key\n');
%         % pause;
% 
%         fprintf('result of non-rigid deformation');
%         clf; visFitDir([p.fitDir '/' scanName],2);
%         fprintf(', press any key\n');
%         % pause;
    fprintf('\n--------------------------------------------------------------------------\n');
end


%% learn PCA model
learnPCA(expidx);
createModel(expidx);
fprintf('done\n');