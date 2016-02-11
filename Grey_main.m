function [error,iterations] = Greylevels_fixed_without_middle()
tol=0.1; %tolerance of convergence (stopping rule)
MAX_NUM_ITER=500;

% DEFINING THE NOSE POINTS CONNECTIONS (in adjacent order -> for bisectors calculations)
connections=zeros(6,4); %6 points -- 4 possible connections
connections(1,1)=5;connections(1,2)=2;connections(1,3)=6;
connections(2,1)=1;connections(2,2)=3;
connections(3,1)=2;connections(3,2)=6;connections(3,3)=4;connections(3,4)=5;
connections(4,1)=5;connections(4,2)=3;connections(4,3)=6;
connections(5,1)=1;connections(5,2)=3;connections(5,3)=4;
connections(6,1)=1;connections(6,2)=3;connections(6,3)=4;

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
     NosePoints = GetStablePoints(Xin1); %FER - maybe we need to do this 
    [NosePoints, Centers] = CentreThis(NosePoints);   
    [NosePoints, scales,MeanRefDist] = ReScale(NosePoints);
    clear Xin1;
    szx = size(NosePoints);
    m = szx(3);
    
%---------------------------------- DRAWING - READ IMAGE AND COMPUTE CENTRES AND SCALE, just as a guess
target_image=6;
I=imread('images/0.bmp');
%guess centre and scale :   frontal view
centre=Centers(6,:); scale=scales(6);

target=NosePoints(:,:,target_image); %for error supervising purposes

FV=NosePoints(:,:,6); 

szfv=size(FV);
num_corners=szfv(2);

global NUM_MIDDLE_POINTS; NUM_MIDDLE_POINTS=2; %number of middle points between pair of corner points
global SIZE_PROFILE; SIZE_PROFILE=4; %number of pixels in one wing of the profile (excluding the central one)

FV_old=FV+50; %just to start the first time

%centre to image space and then uncentre back with the guessed centre, just to show the initial error correctly
iterations=1;
Y=UnScale(FV,scales(target_image)); 
Y=UnCentre(Y,Centers(target_image,:));
Y=UnCentre(Y,-centre);
Y=UnScale(Y,1/scale);
errors(iterations)=norm(Y-target,'fro')/sqrt(num_corners);

% Option 1: Add middle points here (once, and then work with them). Stopping rule may take middle points into account or not
% look_up is a table indicating the corner points involved in each CONNECTION (not actual middle points!)
% [FV,look_up]=AddMiddlePoints(FV,connections);
% [basis1,look_up]=AddMiddlePoints(NosePoints(:,:,1),connections);
% [basis2,look_up]=AddMiddlePoints(NosePoints(:,:,m),connections);
basis1=NosePoints(:,:,1);
basis2=NosePoints(:,:,m);
look_up=1; %just to ignore it
% Basis views profiles
szi=size(I);
% here the basis profiles are built. profiles1wing2, for example, means profiles of basis view 1, wing corresponding to the second segment (defined by segments1)
%the profilesXwingX are structured as follows: [point,bisector,profile], where the profile is of length SIZE_PROFILE + 1 (common pixel included on both wings, and always first)
%the segmentsX are structured as follows: [point, bisector, (v_wing1,w_wing1,v_wing2,w_wing2)], where v,w are vectors centred at the common pixel and ending where the profile should end
[profiles1wing1,profiles1wing2,segments1,profiles2wing1,profiles2wing2,segments2]=BuildBasisProfiles(basis1,basis2,scales,Centers,FV,look_up,connections,m,szi);
T=GetNewT(basis1,basis2,FV); %intitialize T (to join profile segments and calculate weights)
% Structure of FV [corner points, middle points between corner #1 and #2, middle points bt #1 #5, and so on] the exact order can be calculated using look_up
%%%%%%%%%%%%%%%%%%%% MAIN PROGRAM LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ((StoppingRule(FV_old(:,[1:num_corners]),FV(:,[1:num_corners]),tol) & (iterations < MAX_NUM_ITER)) | (iterations==1))
    FV_i=UnScale(FV,scale);   %taking FV to image space of the target image (guessed)
    FV_i=UnCentre(FV_i,centre); %taking FV to image space of the target image (guessed)
    points_found=FindPoints(FV_i,I,connections,look_up,profiles1wing1,profiles1wing2,profiles2wing1,profiles2wing2,segments1,segments2,T); %Greylevels mode
    centre=UpdateCentre(points_found(:,[1:num_corners]));%updating centre with the corner points only, !not centering
    points_found=UnCentre(points_found,-centre); %centering ALL the points
    if ( (points_found(:,1)==[0 0]') & (points_found(:,2)==[0 0]') ) %points found are [0 0]
        fprintf('ERROR: points found are [0 0], unknown reason. EXITING\n');
        break;
    end
    scale=UpdateScale(points_found(:,[1:num_corners]),MeanRefDist);%updating scale with the corner points only, !not scaling
    points_found=UnScale(points_found,1/scale); %scaling    
    % CATT 
    T=GetNewT(basis1,basis2,points_found);
    FV_old=FV;
    %FV=GetNewXAffine(basis1,basis2,a,b);
    FV=GetNewX(basis1,basis2,T);
    
    iterations=iterations+1;
    
    % % view FV_old FV - using centres and scales from '6', for example
    % if option 1, DOnt take into account middle points! (remove them temporarily)
    Y=UnScale(FV(:,[1:num_corners]),scale);
    Y=UnCentre(Y,centre);
    S=DrawPoints(Y,S,ceil((100+iterations)/((100+450)/255)));%substract 200 to have the image "zoomed" in the area of interest
    
    errors(iterations)=norm(FV(:,[1:num_corners])-target,'fro')/sqrt(num_corners);
end
