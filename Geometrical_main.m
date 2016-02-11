function [error,iterations] = Geometry(a,b,c,d,t1,t2,tol,MAX_NUM_ITER,target_image)

% ***********************************************************************************************
% ***********************************************************************************************
% 
%       The set of given images should be in the following form, 
%           saved in an MS Excel spreadsheet.
% 
%  Row\Column      1          2         3      ..          n
%   1           x(1,1)     x(1,2)     x(1,3)   ...       x(1,n)       
%   2           y(1,1)     y(1,2)     y(1,3)   ...       y(1,n)       
%   :              :          :          :      :          :
% (2k-1)        x(k,1)     x(k,2)     x(k,3)   ...       x(k,n)
%   k           y(k,1)     y(k,2)     y(k,3)   ...       y(k,n)
% 
%   Here {x(i,j),y(i,j)} is the jth landmark point on the ith image.
% 
% ***********************************************************************************************
% ***********************************************************************************************

% DEFINING THE NOSE POINTS CONNECTIONS (in adjacent order -> for bisectors calculations)
connections=zeros(6,4); %6 points -- 4 possible connections
connections(1,1)=5;connections(1,2)=2;connections(1,3)=6;
connections(2,1)=1;connections(2,2)=3;
connections(3,1)=2;connections(3,2)=6;connections(3,3)=4;connections(3,4)=5;
connections(4,1)=5;connections(4,2)=3;connections(4,3)=6;
connections(5,1)=1;connections(5,2)=3;connections(5,3)=4;
connections(6,1)=1;connections(6,2)=3;connections(6,3)=4;

%-------------------------------------------------------------------------------
%   Reading in Training Data
%-------------------------------------------------------------------------------
% Reading in the original images from the MS Excel file
%     [fname,pname] = uigetfile(,'Reading in the original images from file ......\n\n');
%     filename = fullfile(pname,fname);
W = xlsread('Landmark Points for processing.xls');
sW = size(W);
k = sW(1)/2;  % The number of images in the set
n = sW(2);   % The number of landmark points in each image

index = 1;
for i = 1:k
    Xin1(1,:,i) = W(index,:); % FER - index could be 2*i -1
    index = index + 1;
    Xin1(2,:,i) = W(index,:); % FER - index could be 2*i 
    index = index + 1;
end

clear W;
NosePoints = GetStablePoints(Xin1); 

[NosePoints, Centers] = CentreThis(NosePoints);   

[NosePoints, scales,MeanRefDist] = ReScale(NosePoints);

clear Xin1;

szx = size(NosePoints);
m = szx(3);


%         % ------------ test - check that image '6' can be recostructed from T -----------------
%         X=getNewX(NosePoints(:,:,1),NosePoints(:,:,m),T);
%         X=UnScale(X,scales(6));
%         X=UnCentre(X,Centers(6,:));
%         % make line drawing of generated image
%         I=zeros(1,400,400); % FER - CHANGE - size of image?
%         I(1,:,:)=DrawLine(X(:,1),X(:,2),I(1,:,:));
%         I(1,:,:)=DrawLine(X(:,1),X(:,5),I(1,:,:));
%         I(1,:,:)=DrawLine(X(:,1),X(:,6),I(1,:,:));
%         I(1,:,:)=DrawLine(X(:,2),X(:,3),I(1,:,:));
%         I(1,:,:)=DrawLine(X(:,3),X(:,4),I(1,:,:));
%         I(1,:,:)=DrawLine(X(:,3),X(:,5),I(1,:,:));
%         I(1,:,:)=DrawLine(X(:,3),X(:,6),I(1,:,:));
%         I(1,:,:)=DrawLine(X(:,4),X(:,5),I(1,:,:));
%         I(1,:,:)=DrawLine(X(:,4),X(:,6),I(1,:,:));
%         % Showing result
%         figure(666);
%         colormap(gray);
%         imagesc(squeeze(I(1,:,:)));    

%---------------------------------- DRAWING - READ IMAGE AND COMPUTE CENTRES AND SCALE, just as a guess
I=imread(['drawing' num2str(target_image) '.bmp']);
[M,Points,centre,scale]=Image2Matrix(I,MeanRefDist); %we need MeanRefDist to scale
% % NOTE: Points in any order!!!!!!!!... connections inconsistent - Should not use Points - may use I or M
% %if we only had the image, we wouldnt be able to calculate the error

%this two following lines are to do the centred and scaled testing and the basins of attraction
centre=Centers(target_image,:); % test... remove this
scale=scales(target_image);     % test... remove this
%test for basis of attracton of AFFINE/CATT_translation in x ('y' in matlab)
centre(1)=centre(1)+t1;
centre(2)=centre(2)+t2;

target=NosePoints(:,:,target_image); %for error supervising purposes

%---------------------------------- FIND INITIAL POSE --------------------------

% get frontal view
%T=getNewT(NosePoints(:,:,1),NosePoints(:,:,m),NosePoints(:,:,6)); % using basis and frontal view (#6)
%FV=getNewX(NosePoints(:,:,1),NosePoints(:,:,m),T);
FV=NosePoints(:,:,6); 

%------------------------------------ PREVIOUS TEST ---- <> ---------------

%FV=FV(:,[1 4 5 6]); %QUAD
szfv=size(FV);
num_corners=szfv(2);
% defining the quad point connections
connections_quad=zeros(num_corners,2); %4 points, 2 connections
connections_quad(1,1)=3; connections_quad(1,2)=4;
connections_quad(2,1)=3; connections_quad(2,2)=4;
connections_quad(3,1)=1; connections_quad(3,2)=2;
connections_quad(4,1)=1; connections_quad(4,2)=2;
%        connections=connections_quad; %THIS LINE TO TEST FOR THE KITE (QUAD)
global NUM_MIDDLE_POINTS; NUM_MIDDLE_POINTS=2; %DEFINITION 
global SIZE_PROFILE; SIZE_PROFILE=4; %number of pixels in one wing of the profile (excluding the central one)

%warp starting point- SIMULATION
[FV,A]=AffineTrans(FV,a,b,c,d,t1,t2);

%testing with known inputs
%diamond within diamond
%         FV=[0 0 -20 20;20 -20 0 0];
%        target=[ 0 0 -10 10;10 -10 0 0];
%one shifted to the right
%       FV=[0 0 -20 20;20 -20 0 0];
%        target=[ 5 5 -15 25;10 -10 0 0];       

%target=AddNoise(target);

% SEARCH - looking for the "correctly affined" from each new Frontal View
FV_old=FV+50; %just to start the first time
iterations=1;

%centre to image space and then uncentre back with the guessed centre, just to show the initial error
Y=UnScale(FV,scales(target_image)); 
Y=UnCentre(Y,Centers(target_image,:));
Y=UnCentre(Y,-centre);
Y=UnScale(Y,1/scale);
errors(iterations)=norm(Y-target,'fro')/sqrt(num_corners);

% %image show
% S=I;
% Y=UnScale(FV,scale);
% Y=UnCentre(Y,centre); 
% S=DrawPoints(Y,S,40); 

% Option 1: Add middle points here (once, and then work with them... only taking care when calculating bisectors). Stopping rule may take middle points into account too
[FV,look_up]=AddMiddlePoints(FV,connections);
[basis1,look_up]=AddMiddlePoints(NosePoints(:,:,1),connections);
[basis2,look_up]=AddMiddlePoints(NosePoints(:,:,m),connections);

% Structure of FV [corner points    middle points between corner #1 and #2
while ((StoppingRule(FV_old(:,[1:num_corners]),FV(:,[1:num_corners]),tol) & (iterations < MAX_NUM_ITER)) | (iterations==1))
    FV_i=UnScale(FV,scale);   %taking FV to image space of the target image (guessed)
    FV_i=UnCentre(FV_i,centre); %taking FV to image space of the target image (guessed)
    points_found=FindPointsDrawing(FV_i,I,connections,look_up); %Drawing mode 
    %remove points not found from basis and 'points_found'
    [basis1_aux,basis2_aux,points_found,num_corners_found]=RemoveNotFound(basis1,basis2,points_found); %useless with greylevels
    if (num_corners_found<4) %for the AFFINE case we need 4 points IN TOTAL, for the CATT 3, we will generalise to 4
        fprintf('WARNING: less than 4 corner points found. EXITING\n');
        break; %
    elseif (num_corners_found<6) %some corner points not found
        %the centre and scale calculation will be very biased, so do not update them, use the old ones
        points_found=UnCentre(points_found,-centre); %centering
        points_found=UnScale(points_found,1/scale); %scaling        
    else
        centre=UpdateCentre(points_found(:,[1:num_corners]));%updating centre with the corner points only, !not centering
        points_found=UnCentre(points_found,-centre); %centering ALL the points
        if ( (points_found(:,1)==[0 0]') & (points_found(:,2)==[0 0]') ) %points found are [0 0]
            fprintf('ERROR: points found are [0 0], unknown reason. EXITING\n');
            break;
        end
        scale=UpdateScale(points_found(:,[1:num_corners]),MeanRefDist);%updating scale with the corner points only, !not scaling
        points_found=UnScale(points_found,1/scale); %scaling
    end
    %             % AFFINE - Find the transformation newA that is a correct affine trasformation of FV from point_found
    %             %make homogenous (to obtain a 3x3 newA - Homography)
    %             %points_found(3,:)=1; %allow for scaling
    %             %FV(3,:)=1; %allow translation
    %             newA=points_found*pinv(FV);
    %             FV_old=FV([1 2],:); %remove homogenity
    %             FV=newA*FV;
    %             FV=FV([1 2],:); %remove homogenity
    
    %Affine 2:1
    %     [a,b]=GetNewab(basis1_aux,basis2_aux,points_found);
    %     FV_old=FV;
    %     FV=GetNewXAffine(basis1,basis2,a,b);
    
    %CATT 
    T=GetNewT(basis1_aux,basis2_aux,points_found);
    FV_old=FV;
    FV=GetNewX(basis1,basis2,T);
    
    iterations=iterations+1;
    
    % % view FV_old FV - using centres and scales from '6', for example
    % if option 1, DOnt take into account middle points! (remove them temporarily)
    %     Y=UnScale(FV(:,[1:num_corners]),scale);
    %     %             Y=UnScale(FV,scales(6));
    %     % 
    %     Y=UnCentre(Y,centre);
    %     %              % I=DrawPoints(X,I,150);
    %     S=DrawPoints(Y,S,ceil((100+iterations)/((100+MAX_NUM_ITER)/255)));%the operations are just to let the pixels have fancy colours
    
    errors(iterations)=norm(FV(:,[1:num_corners])-target,'fro')/sqrt(num_corners);
end
error=errors(iterations);
