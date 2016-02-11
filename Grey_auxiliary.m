
%*******************************************************************************************************
% FindQuad - Bisectors (DoBisector), Perpendicular, Intersections (IntersectionsImage,CalculateIntersection,DoLine2Points,DoLine)
%*************************************************************************************************

%NOTE: if(errors(i)==999), do not take point into account (bisector didn't exist)
function best=GetBestMatch(point,fits,errors)
%OPTION 1 - The lowest error
[Y,I]=min(errors);
best=fits(:,I);
% OPTION 2 - closest - not implemented yet
% OPTION 3 - weighted interpolation  - not implemented yet

% this function does the image search for the greylevel case
function found=FindPoints(X,I,con,look_up,profiles1wing1,profiles1wing2,profiles2wing1,profiles2wing2,segments1,segments2,T); %Greylevels mode
global SIZE_PROFILE;
found=X; %initialization
[a,b]=CATT2Affine(T);
[w1,w2]=GetWeights(a,b);
segments=JointSegments(segments1,segments2,T); %calculate weighted profile segments 
wing1=JointProfiles(profiles1wing1,profiles2wing1,w1,w2); % calculate weighted profiles
wing2=JointProfiles(profiles1wing2,profiles2wing2,w1,w2); % calculate weighted profiles
szs=size(segments);
num_points=szs(1);
num_bis=szs(2);
for i=1:num_points
    point=X(:,i);
    errors=zeros(1,num_bis)+999; %initialisation - fit error for each bisector
    best_fits=zeros(2,num_bis); %initialization - best point that fit for each bisector
    for j=1:num_bis
        error1=999; error2=999; % initialization - fit error of each of the wings
        if ((segments(i,j,1)==0) & (segments(i,j,2)==0))%no more bisectors for this point
            break;
        end 
        %wing1
        for k=1:SIZE_PROFILE*3 %limit = 3xSIZE_PROFILE (one wing)
            v=squeeze(segments(i,j,[1:2])); %search vector
            %normalise... defining search step
            d=sqrt(v(1)^2+v(2)^2);
            vu=v*sqrt(2)/d; %search step is sqrt(2) (asure we advance one pixel)
            %match profiles(i,j,:) with image 'profile' (offseted)
            profile=GetProfile(I,point+vu*k,v);
            error=norm(double(profile)-squeeze(wing1(i,j,:))','fro')/sqrt(SIZE_PROFILE+1);
            if (error<error1)
                error1=error;
                best_fit1=point+vu*k;
            end
            %we must search in the other direction too
            profile=GetProfile(I,point-vu*k,v);
            error=norm(double(profile)-squeeze(wing1(i,j,:))','fro')/sqrt(SIZE_PROFILE+1);
            if (error<error1)
                error1=error;
                best_fit1=point-vu*k;
            end
        end %for search wing1
        %wing2
        for k=1:SIZE_PROFILE*3 %limit = 3xSIZE_PROFILE (one wing)
            v=squeeze(segments(i,j,[3:4])); %search vector
            %normalise... defining search step
            d=sqrt(v(1)^2+v(2)^2);
            vu=v*sqrt(2)/d; %search step is sqrt(2) (asure we advance one pixel)
            %match profiles(i,j,:) with image 'profile' (offseted)
            profile=GetProfile(I,point+vu*k,v);
            error=norm(double(profile)-squeeze(wing2(i,j,:))','fro')/sqrt(SIZE_PROFILE+1);
            if (error<error2)
                error2=error;
                best_fit2=point+vu*k;
            end
            %we must search in the other direction too
            profile=GetProfile(I,point-vu*k,v);
            error=norm(double(profile)-squeeze(wing2(i,j,:))','fro')/sqrt(SIZE_PROFILE+1);
            if (error<error2)
                error2=error;
                best_fit2=point-vu*k;
            end
        end %for search wing1   
        %Join both points found together using error1 error2
        best_fits(:,j)=(error1*best_fit2 + error2*best_fit1)/(error1+error2); %weighted average
        errors(j)=2*error1*error2/(error1+error2); %weighted average
    end %bisctors
    %chose the point found amongst found in each bisector
    found(:,i)=GetBestMatch(point,best_fits,errors);
end %points

function centre=UpdateCentre(Xin)
szx=size(Xin);
n=szx(2);
d = 0;
totX = 0;
totY = 0;
for j = 1:n
    d = d + 1;
    totX = totX + Xin(1,j);
    totY = totY + Xin(2,j);
end
centre(1) = totX/d; % FER -   d=n?
centre(2) = totY/d;        

function scale=UpdateScale(Xin,MeanRefDist)
szx = size(Xin);
n = szx(2);
r = 0;
dist = 0;
for i = 1:n
    dist = dist + sqrt((Xin(1,i)*Xin(1,i)) + (Xin(2,i)*Xin(2,i)));
end
if (dist==0)
    fprintf('\nERROR: distance 0!?\n');
end
scale = MeanRefDist*n/dist; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GREYLEVEL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute transformation for one triangle - return as one dimensional array
function T=ComputeTransformation(Origin,Target)
%put origin in homogeneous
Origin(3,:)=1;
Inv=pinv(Origin);
if (Inv==Inf) %was singular - I dont know if pinv check this itself
    Inv=pinv(Origin);
end
T=reshape(Target*Inv,1,6);

%compute affine transformations of all triangles of one view - return as a one dimensional array
function trans=ComputeTransformations(FV,basis,box,triangles);
szt=size(triangles);
num_triangles=szt(1);
trans=zeros(num_triangles,6);
for i=1:num_triangles
    tri=triangles(i,:); %get points of triangle (indices)
    a=tri(1);b=tri(2);c=tri(3);
    if (a<0) %point of the box
        A=box(:,abs(a));
        A_basis=box(:,abs(a));
    else %point of the nose
        A=FV(:,a);
        A_basis=basis(:,a);
    end
    if (b<0) %point of the box
        B=box(:,abs(b));
        B_basis=box(:,abs(b));
    else %point of the nose
        B=FV(:,b);
        B_basis=basis(:,b);
    end
    if (c<0) %point of the box
        C=box(:,abs(c));
        C_basis=box(:,abs(c));
    else %point of the nose
        C=FV(:,c);
        C_basis=basis(:,c);
    end  
    trans(i,:)=ComputeTransformation([A B C],[A_basis B_basis C_basis]);
end

function inside=PointInsideTriangle(point,A,B,C)
ab=cross([(B-A)' 0],[(point-A)' 0]);
z_ab=ab(3);
bc=cross([(C-B)' 0],[(point-B)' 0]);
z_bc=bc(3);
ca=cross([(A-C)' 0],[(point-C)' 0]);
z_ca=ca(3);
if (  ((z_ab<=0)&(z_bc<=0)&(z_ca<=0)) | ((z_ab>=0)&(z_bc>=0)&(z_ca>=0))  ) %checking if all positive just in case the points are not ordered clockwise
    inside=1;
else
    inside=0;
end

% checks for point being inside each of the triangels in 'possible' (only indices), triangles is the list of triangles to look for the points involved (only indices of points); FV and boox to look for the acutal coordinates
function t=WhichTriangle(point,possible,triangles,FV,box);
num_possible=length(possible);
t=0;
for i=1:num_possible
    tri=triangles(possible(i),:); %get points of triangle (indices)
    a=tri(1);b=tri(2);c=tri(3);
    if (a<0) %point of the box
        A=box(:,abs(a));
    else %point of the nose
        A=FV(:,a);
    end
    if (b<0) %point of the box
        B=box(:,abs(b));
    else %point of the nose
        B=FV(:,b);
    end
    if (c<0) %point of the box
        C=box(:,abs(c));
    else %point of the nose
        C=FV(:,c);
    end    
    if PointInsideTriangle(point,A,B,C) %found
        t=possible(i);
        break;
    end
end
if (t==0) %not inside any triangle - could be because end point gets out of triangle (in that case, use the closest)
    fprintf('WARNING: In WhichTriangle, point not inside expected triangles, getting closest instead\n');
    % first, find the triangle in wich the point really is
    all=[1 2 3 4 5 6 7 8 9 10 11 12]; %all the triangles
    rest=setxor(all,possible); %all except those already tried
    num_rest=length(rest);
    for i=1:num_rest
        tri=triangles(rest(i),:); %get points of triangle (indices)
        a=tri(1);b=tri(2);c=tri(3);
        if (a<0) %point of the box
            A=box(:,abs(a));
        else %point of the nose
            A=FV(:,a);
        end
        if (b<0) %point of the box
            B=box(:,abs(b));
        else %point of the nose
            B=FV(:,b);
        end
        if (c<0) %point of the box
            C=box(:,abs(c));
        else %point of the nose
            C=FV(:,c);
        end    
        if PointInsideTriangle(point,A,B,C) %found
            t_aux=rest(i);
            break;
        end    
    end %for rest
    if (t_aux==0)
        fprintf('ERROR: In WhichTriangle, point not inside ANY triangle!!!');
    end
    % point is inside 't_aux', but we want the closest among one of the possible. We suppose 't_aux' and the one we want are adjacent
    % if they are adjacent -> they share two vertex, no other triangle in 'possible' shares also two vertex (in our nose) 
    for i=1:num_possible
        if (intersect(triangles(t_aux,:),triangles(possible(i),:)) == 2) % two common vertex
            t=i;
            break;
        end
    end
    if (t==0) % 't_aux' wasn't adjacent to any triangle in 'possible' -> return 't_aux'
        t=t_aux;
    end
end %t==0

%returns vector from basis to transformed end_point
function segment_vector=BasisEndVector(end_point,basis,trans,possible,triangles,FV,box)
t=WhichTriangle(end_point,possible,triangles,FV,box);
T=reshape(trans(t,:),2,3); %buiding transformation matrix
end_point(3)=1; %homogenous
basis_end=T*end_point; %affine transformation
segment_vector=basis_end-basis;


%Ben's asymetry - returns the two vectors going from each central landmark point to the segment end 
function [segments1,segments2]=BasisSegments(FV,basis1,basis2,box,bis,trans1,trans2,look_up,look_up_triangles,triangles)
global SIZE_PROFILE;
% global NUM_MIDDLE_POINTS;
szb=size(bis);
num_points=szb(1);
num_bisectors=szb(2);
segments1=zeros(num_points,num_bisectors,4); %initialization - last component is 4 to contain both vector-segments
segments2=zeros(num_points,num_bisectors,4); %initialization - last component is 4 to contain both vector-segments
% szl=size(look_up);
% num_connections=szl(1);
% num_middle_points=NUM_MIDDLE_POINTS*num_connections;
% num_corners=num_points-num_middle_points;
for i=1:num_points
    %find possible triangles
%     if (i>num_corners) %middle point
%         con_num=ceil((i-num_corners)/NUM_MIDDLE_POINTS);
%         possible=intersect(look_up_triangles(look_up(con_num,1),:),look_up_triangles(look_up(con_num,2),:));
%     else
        possible=look_up_triangles(i,:);
%     end
    %remove zeros in 'possible'
    ind=find(possible==0);
    possible(ind)=[];
    %for each bisector get the two segment ends
    for j=1:num_bisectors
        if ((bis(i,j,1)==0) & (bis(i,j,2)==0))%no more bisectors for this point
            break;
        end 
        v(1)=-bis(i,j,2);
        v(2)=bis(i,j,1);
        %normalise... again?
        d=sqrt(v(1)^2+v(2)^2);
        v=v/d;
        %one end of the segment - for both views
        end_point=FV(:,i)+v'*sqrt(2)*SIZE_PROFILE;
        segments1(i,j,1:2)=BasisEndVector(end_point,basis1(:,i),trans1,possible,triangles,FV,box); 
        segments2(i,j,1:2)=BasisEndVector(end_point,basis2(:,i),trans2,possible,triangles,FV,box); 
        %the other end of the segment - for both views
        end_point=FV(:,i)-v'*sqrt(2)*SIZE_PROFILE;
        segments1(i,j,3:4)=BasisEndVector(end_point,basis1(:,i),trans1,possible,triangles,FV,box); 
        segments2(i,j,3:4)=BasisEndVector(end_point,basis2(:,i),trans2,possible,triangles,FV,box); 
    end
end %for points

%gets the greylevels of one wing of a particular point
function profile=GetProfile(I,point,v)
global SIZE_PROFILE;
profile(1)=I(round(point(1)),round(point(2)));
for k=1:SIZE_PROFILE
    profile(k+1)=I(  round( point(1)+v(1)*k/SIZE_PROFILE ),round( point(2)+v(2)*k/SIZE_PROFILE )  );  
end

%gets the greylevels of all the segments of all the points - both wings
function [wing1,wing2]=BasisProfiles(I,points,segments)
szs=size(segments);
num_points=szs(1);
num_seg=szs(2);
wing1=zeros(num_points,num_seg,5); %5 values counting the middle point
wing2=zeros(num_points,num_seg,5); %5 values counting the middle point
for i=1:num_points
    for j=1:num_seg
        if((segments(i,j,1)==0) & (segments(i,j,2)==0) & (segments(i,j,3)==0) & (segments(i,j,4)==0))  %no more segments for this point
            break;
        end
        wing1(i,j,:)=GetProfile(I,points(:,i),segments(i,j,[1:2]));
        wing2(i,j,:)=GetProfile(I,points(:,i),segments(i,j,[3:4]));
    end
end

function [profiles1wing1,profiles1wing2,segments1,profiles2wing1,profiles2wing2,segments2]=BuildBasisProfiles(basis1,basis2,scales,Centers,FV,look_up,connections,m,szi)
basis1_uncentred=UnScale(basis1,scales(1));
basis1_uncentred=UnCentre(basis1_uncentred,Centers(1,:));
FV_uncentred=UnScale(FV,scales(6));
FV_uncentred=UnCentre(FV_uncentred,Centers(6,:));
basis2_uncentred=UnScale(basis2,scales(m));
basis2_uncentred=UnCentre(basis2_uncentred,Centers(m,:));
%defining outer limits (image box)
box=[0 0 szi(1) szi(1);0 szi(2) szi(2) 0];
% points belonging to each triangle (clockwise) - for all views - negative means outer limits (image box)
triangles=[1 -1 -2;1 -2 6;1 6 3;1 3 5; 1 5 -1;3 6 4;3 4 5;4 6 -3;4 -3 -4;4 -4 5;5 -4 -1;6 -2 -3];
% compute affine transformations to use for each triangle (both basis views)
trans1=ComputeTransformations(FV_uncentred,basis1_uncentred,box,triangles);
trans2=ComputeTransformations(FV_uncentred,basis2_uncentred,box,triangles);
% look up table to know which triangles are adjacent to which points (used to reduce number of possible triangles)
look_up_triangles=[1 2 3 4 5;3 4 0 0 0;3 4 6 7 0;6 7 8 9 10;4 5 7 10 11;2 3 6 8 12];
%bisectors for FV
bisectors_FV=Bisectors(FV_uncentred,connections); %own lines -> Bisectors2
% bisectors_FV=Perpendiculars(FV_uncentred,connections,look_up,bisectors_FV); %maybe dont need con
%get segments for both views
[segments1,segments2]=BasisSegments(FV_uncentred,basis1_uncentred,basis2_uncentred,box,bisectors_FV,trans1,trans2,look_up,look_up_triangles,triangles);
%get profiles
BASIS1=imread('images/-25.bmp');
BASIS2=imread('images/25.bmp');
%each view has a profile wing for each segment-vector
[profiles1wing1,profiles1wing2]=BasisProfiles(BASIS1,basis1_uncentred,segments1);
[profiles2wing1,profiles2wing2]=BasisProfiles(BASIS2,basis2_uncentred,segments2);

function segments=JointSegments(segments1,segments2,T) %calculate weighted profile segments 
segments=segments1; % initialization
szs=size(segments1);
num_points=szs(1);
num_bis=szs(2);
for i=1:num_points
    for j=1:num_bis
        if ((segments(i,j,1)==0) & (segments(i,j,2)==0))%no more bisectors for this point
            break;
        end       
        segments_aux=GetNewX(reshape(segments1(i,j,:),2,2),reshape(segments2(i,j,:),2,2),T);
        segments(i,j,:)=reshape(segments_aux,1,4);
    end
end

function profiles=JointProfiles(profiles1,profiles2,w1,w2) % calculate weighted profiles
profiles=profiles1*w1+profiles2*w2;

function [w1,w2]=GetWeights(a,b)
d1=a(3)^2+a(4)^2+b(3)^2+b(4)^2;
d2=a(1)^2+a(2)^2+b(1)^2+b(2)^2;
w1=d2/(d1+d2);
w2=d1/(d1+d2);

%%%%%
% END OF GREYLEVEL FUNCTIONS
%%%%%%%%%%%%    

