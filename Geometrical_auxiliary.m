
%cuts out the points not found ([666 666] in 'points_found'), from points_found and the basis views
%returns number of corner points found (to update centre and scale with the corner points found only. 
%We suppose there are 6 corner points, could calculate it form num_points and NUM_MIDDLE_POINTS
function [basis1_aux,basis2_aux,points_found,num_corners_found]=RemoveNotFound(basis1,basis2,points_found)
indices=find(points_found(1,:)==666);
basis1_aux_1=basis1(1,:);
basis1_aux_1(indices)=[];
basis1_aux_2=basis1(2,:);
basis1_aux_2(indices)=[];
basis1_aux=[basis1_aux_1;basis1_aux_2];

basis2_aux_1=basis2(1,:);
basis2_aux_1(indices)=[];
basis2_aux_2=basis2(2,:);
basis2_aux_2(indices)=[];
basis2_aux=[basis2_aux_1;basis2_aux_2];

points_found_1=points_found(1,:);
points_found_1(indices)=[];
points_found_2=points_found(2,:);
points_found_2(indices)=[];
points_found=[points_found_1;points_found_2];

corners=find(indices<7);
num_corners_found=6-length(corners);

%OUT: look_up table of corner points involved in each connection number, used when calculating perpendicular lines to middle points
function [added,look_up]=AddMiddlePoints(P,con)
global NUM_MIDDLE_POINTS;
szc=size(con);
num_points=szc(1);
num_con=szc(2);
added=P;
count=1;
%calculate middle points
for i=1:num_points
    for j=1:num_con
        if(con(i,j)==0) %no more connections for this point
            break;
        end
        if(i<con(i,j)) %asure only once per connection
            step=(P(:,i)-P(:,con(i,j)))/(NUM_MIDDLE_POINTS+1);
            for k=1:NUM_MIDDLE_POINTS
                added(:,num_points+count)=P(:,con(i,j))+step*k;
                count=count+1;
            end
            look_up((count-1)/NUM_MIDDLE_POINTS,1)=i;
            look_up((count-1)/NUM_MIDDLE_POINTS,2)=con(i,j);
            %implicit else - connection already dealt with
        end
    end
    %cannot do intersections here because cannot use point info anymore, should check all lines
end


% Not taking middle points into account (could, if option 1)
function bool=StoppingRule(old,new,tol)
szn=size(new);
num_points=szn(2);
dist=norm((old-new),'fro')/sqrt(num_points); % RMS
%dist=sum(sqrt((dif'*dif)));
if (dist<tol)
    bool=0; %stop
else
    bool=1; %continue
end

%*******************************************************************************************************
% FindQuad - Bisectors (DoBisector), Perpendicular, Intersections (IntersectionsImage,CalculateIntersection,DoLine2Points,DoLine)
%*************************************************************************************************

function found=FindPointsDrawing(X,I,con,look_up)
bisectors=Bisectors(X,con); %own lines -> Bisectors2
bisectors=Perpendiculars(X,con,look_up,bisectors); %maybe dont need con
found=IntersectionsImage(bisectors,X,I); %X to know where to start 

function found=FindPoints(X,P,con,look_up) %pure geometry - check to accept middle points
bisectors=Bisectors(X,con); %own lines -> Bisectors2
bisectors=Perpendiculars(X,con,look_up,bisectors); %maybe dont need con
found=Intersections(bisectors,X,P,con); %X necessary for distances, con necessary to calculate line equ for P

function found=FindPoints2(X,P,con) %neighbourhood alternative - pure geometry
linesP=DoLineEqsFromPoints(P,con);
szx=size(X);
num_points=szx(2);
found=zeros(2,num_points);
bisectors_calculated=0; %if one point fail, only calculate "bisector alternative" once
for i=1:num_points
    aux=Neighbour(X(:,i),linesP,P)'; %P necessary to limit intersections
    if ((aux(1)==0) & (aux(2)==0) & bisectors_calculated==0) %neighbour not found, first time
        found=FindPoints(X,P,con);
        bisectors_calculated=1;
        %break; %witout 'break' rest of points calculated by neighbourhood if possible
    elseif ((aux(1)~=0) | (aux(2)~=0)) %neighbourhood found
        found(:,i)=aux;
        %implicit else - neighbourhood not found, second time-> do nothing
    end
end

function closest=NeighbourImage(start,I)
closest=start;
step=1;
limit=100/step;
iter=1;
done=0;
while ((iter<limit) & (done==0))
    belongs=CheckImage([start(1)+step*iter start(2)],I);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1)+step*iter start(2)];
        break;
    end
    belongs=CheckImage([start(1)-step*iter start(2)],I);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1)-step*iter start(2)];
        break;
    end
    belongs=CheckImage([start(1) start(2)+step*iter],I);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1) start(2)+step*iter];
        break;
    end
    belongs=CheckImage([start(1) start(2)-step*iter],I);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1) start(2)-step*iter];
        break;
    end
    belongs=CheckImage([start(1)+step*iter start(2)+step*iter],I);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1)+step*iter start(2)+step*iter];
        break;
    end
    belongs=CheckImage([start(1)-step*iter start(2)+step*iter],I);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1)-step*iter start(2)+step*iter];
        break;
    end
    belongs=CheckImage([start(1)-step*iter start(2)+step*iter],I);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1)-step*iter start(2)+step*iter];
        break;
    end
    belongs=CheckImage([start(1)+step*iter start(2)-step*iter],I);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1)+step*iter start(2)-step*iter];
        break;
    end
    iter=iter+1;
end
if (done==0)
    fprintf('WARNING: Neighbour not found...  switching to bisectors!');
    closest=[0 0];
end


function closest=Neighbour(start,linesP,P)
closest=start;
step=1;
limit=100/step;
iter=1;
done=0;
while ((iter<limit) & (done==0))
    belongs=Check([start(1)+step*iter start(2)],linesP,P);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1)+step*iter start(2)];
        break;
    end
    belongs=Check([start(1)-step*iter start(2)],linesP,P);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1)-step*iter start(2)];
        break;
    end
    belongs=Check([start(1) start(2)+step*iter],linesP,P);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1) start(2)+step*iter];
        break;
    end
    belongs=Check([start(1) start(2)-step*iter],linesP,P);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1) start(2)-step*iter];
        break;
    end
    belongs=Check([start(1)+step*iter start(2)+step*iter],linesP,P);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1)+step*iter start(2)+step*iter];
        break;
    end
    belongs=Check([start(1)-step*iter start(2)+step*iter],linesP,P);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1)-step*iter start(2)+step*iter];
        break;
    end
    belongs=Check([start(1)-step*iter start(2)+step*iter],linesP,P);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1)-step*iter start(2)+step*iter];
        break;
    end
    belongs=Check([start(1)+step*iter start(2)-step*iter],linesP,P);%P necessary to limit intersections
    if (belongs)
        done=1;
        closest=[start(1)+step*iter start(2)-step*iter];
        break;
    end
    iter=iter+1;
end
if (done==0)
    fprintf('WARNING: Neighbour not found...  switching to bisectors!');
    closest=[0 0];
end

function belongs=CheckImage(x,I)  %line drawing
belongs=I(round(x(1)),round(x(2)));

function belongs=Check(x,linesP,P) %pure geometry
belongs=0;
szl=size(linesP);
num_points=szl(1);
%loop through top part of line-matrix
for k=1:num_points-1
    for l=k+1:num_points
        if((linesP(k,l,1)~=0) & (belongs==0)) %there is a line between these points
            inside=LimitIntersection(x,P(:,k),P(:,l));
            if (inside)
                if (abs(linesP(k,l,1)*x(1)+linesP(k,l,2)*x(2)+linesP(k,l,3))<=abs(linesP(k,l,1)*0.5+linesP(k,l,2)*0.5)) % a*0.5+b*0.5 maximum error with step 1
                    belongs=1;
                    break;
                end
            end
        end
    end
end % for - lines

%%%%%%%%%%%%%%%%% Intersections (DoLine,DoLine2Poinnts,CalculateIntersection,DoLineEqusFromPoints,Intersections,LimitIntersection,IntersectionsImage)

%check the image along the line, starting from x, in a 4-connected way (not to miss a 8-connected drawing)
function int=CalculateIntersectionImage(x,line,I) %cheking in both directions (and with image limit control)
v(1)=-line(2);
v(2)=line(1);
%find the ratio (how many units we must go along one dimension for 1 unit of the other)     NOTE: maintaining the signs!
if (v(1)==0)
    ratio=[0 sign(v(2))]; %vertical line [0 1] [0 -1]
else
    ratio=[sign(v(1)) v(2)/abs(v(1))]; % [+/-1 +/-xxx] xxx can be '0' and 0.xxx
    %we want the smallest to be 1 (and the other more than 1) (except in the case [1 0]), so...
    if ((ratio(2)~=0) & (abs(ratio(2))<1))
        ratio=[sign(ratio(1))/abs(ratio(2)) sign(ratio(2))]; % [+/-xxx +/-1]
    end
end
%round the larger component and calculate remainder (probability of adding an extra pixel to the larger component)
if (abs(ratio(1))==1) %second component larger    NOTE: if [1 0], second component will be classified as larger, and remainder will be '0', this should be ok
    larger_component=2;
    int_part=floor(abs(ratio(2)));
    remainder=abs(ratio(2))-int_part; %probablility of adding an extra pixel to the largest component
    ratio(2)=int_part*sign(ratio(2));
else %first component larger
    larger_component=1;
    int_part=floor(abs(ratio(1)));
    remainder=abs(ratio(1))-int_part;
    ratio(1)=int_part*sign(ratio(1));   
end

szi=size(I);
found=0;
image_end=0;
offset=[0 0];
x=round(x');
acum=remainder; %will tell if we have to add an extra pixel (if >0.5)
if ( max(x>szi) | max(x<[0 0]) )  %limit not reached in one side - else exit
    fprintf('ERROR:out of the image!!!\n');
    int=[666 666];
    break;
end

if (I(x(1),x(2))>0) % check starting point
    found=1;
    int=x+offset;
end
while ((found==0) & (image_end==0))
    for i=1:abs(ratio(1))
        offset(1)=offset(1)+sign(ratio(1));
        if (   ( min(x+offset<=szi) & min(x+offset>[0 0]) ) | ( min(x-offset<=szi) & min(x-offset>[0 0]) )   ) %limit not reached in one side - else exit
            if ( min(x+offset<=szi) & min(x+offset>[0 0]) ) %limit not reached in this direction - else ignore this direction
                if (I(x(1)+offset(1),x(2)+offset(2)))
                    found=1;
                    int=x+offset;
                    break;
                end
            end
            if ( min(x-offset<=szi) & min(x-offset>[0 0]) ) %limit not reached in this direction - else ignore this direction
                if (I(x(1)-offset(1),x(2)-offset(2)))
                    found=1;
                    int=x-offset;
                    break;
                end
            end
        else
            image_end=1;
            break;
        end
    end %for
    if (found==0)
        for i=1:abs(ratio(2))
            offset(2)=offset(2)+sign(ratio(2));
            if (   ( min(x+offset<=szi) & min(x+offset>[0 0]) ) | ( min(x-offset<=szi) & min(x-offset>[0 0]) )   ) %limit not reached in one side - else exit
                if ( min(x+offset<=szi) & min(x+offset>[0 0]) ) %limit not reached in this direction - else ignore this direction
                    
                    if (I(x(1)+offset(1),x(2)+offset(2)))
                        found=1;
                        int=x+offset;
                        break;
                    end
                end
                if ( min(x-offset<=szi) & min(x-offset>[0 0]) ) %limit not reached in this direction - else ignore this direction
                    if (I(x(1)-offset(1),x(2)-offset(2)))
                        found=1;
                        int=x-offset;
                        break;
                    end
                end
            else
                image_end=1;
                break;
            end                
        end %for
    end %if found
    %check 'acum' to add an extra pixel-search
    if (found==0)
        if(acum>=0.5) %go one more
            acum=acum+remainder-1;
            offset(larger_component)=offset(larger_component)+sign(ratio(larger_component));
            if (   ( min(x+offset<=szi) & min(x+offset>[0 0]) ) | ( min(x-offset<=szi) & min(x-offset>[0 0]) )   ) %limit not reached in one side - else exit
                if ( min(x+offset<=szi) & min(x+offset>[0 0]) ) %limit not reached in this direction - else ignore this direction
                    
                    if (I(x(1)+offset(1),x(2)+offset(2)))
                        found=1;
                        int=x+offset;
                        break;
                    end
                end
                if ( min(x-offset<=szi) & min(x-offset>[0 0]) ) %limit not reached in this direction - else ignore this direction
                    if (I(x(1)-offset(1),x(2)-offset(2)))
                        found=1;
                        int=x-offset;
                        break;
                    end
                end
            else
                image_end=1;
                break;
            end       
        else %acum<1
            acum=acum+remainder;
        end %if acum
    end %if found
end %while
if (found==0)
    fprintf('WARNING: IntersectionImage not found!');
    int=[666 666]; %not found
end

%not found are returned as [666 666]
function found=IntersectionsImage(bis,X,I)
sz=size(bis);
num_points=sz(1);
num_bis=sz(2);
found=zeros(2,num_points);
szi=size(I);

%calculate intersections
for i=1:num_points
    closest=[666 666];
    closest_dist=99999;
    if ( max(X(:,i)>[szi(1) szi(2)]') | max(X(:,i)<[0 0]') )  %limit not reached in one side - else exit
        fprintf('ERROR:IntersectionsImage: starting point out of the image!!!');
    else
        for j=1:num_bis
            %calculate all intersections of current bisector with all lines, then get the minimum distance
            if ((bis(i,j,1)==0) & (bis(i,j,2)==0))%no more bisectors for this point
                break;
            end
            int=CalculateIntersectionImage(X(:,i),bis(i,j,:),I);
            % if point not found ... try with neighbour ... or do nothing and give less points to GetNewT        
            %         if (int==[666 666]) %not found 
            %             int=NeighbourImage(X(:,i),I);
            %         end
            if (int~=[666 666])
                %calculate distance between int and X(:,i) - necessary anyway because we might have several bisectors
                dist=sqrt((int(1)-X(1,i))^2 + (int(2)-X(2,i))^2);
                if (dist<closest_dist)
                    closest_dist=dist;
                    closest=int;
                end
            end
        end % for - bisectors
    end %if out of the image
    found(:,i)=closest';  
end % for - points

function eq=DoLine(v,p)
eq(1)=v(2);
eq(2)=-v(1);
eq(3)=v(1)*p(2)-v(2)*p(1);

function eq=DoLine2Points(p1,p2)
eq=DoLine(p1-p2,p1);

function p=CalculateIntersection(line1,line2)
if (line1(2)*line2(1)-line1(1)*line2(2)==0) % no intersection
    p=[999999 999999]; %just to make sure it is not used
else
    p(2)=(line2(3)*line1(1)-line1(3)*line2(1))/(line1(2)*line2(1)-line1(1)*line2(2));
    if (line1(1)==0) %avoid dividing by zero
        p(1)= - (line2(2)*p(2)+line2(3))/line2(1);
    else
        p(1)= - (line1(2)*p(2)+line1(3))/line1(1);
    end
end

function linesP=DoLineEqsFromPoints(P,con)
szc=size(con);
num_points=szc(1);
num_con=szc(2);
%calculate line equations of P
linesP=zeros(num_points,num_points,3+2+2); %simetric matrix with line_eq in its elements. If '0'-> no connection. Contains a,b,c,p1_x,p1_y,p2_x,p2_y
for i=1:num_points
    for j=1:num_con
        if(con(i,j)==0) %no more connections for this point
            break;
        end
        if(linesP(i,con(i,j),1)==0) %no line calculated yet (suppose 'a' in ax+by+c=0 is not '0', if it was, it will just be calculated again)
            linesP(i,con(i,j),:)=[DoLine2Points(P(:,i),P(:,con(i,j))) P(:,i)' P(:,con(i,j))'];
            linesP(con(i,j),i,:)=linesP(i,con(i,j),:); %to maintain symetry
            %implicit else - line already calculated
        end
    end
    %cannot do intersections here because cannot use point info anymore, should check all lines
end

function found=Intersections(bis,X,P,con)
sz=size(bis);
num_points=sz(1);
num_bis=sz(2);
szc=size(con);
num_corners=szc(1);
num_con=szc(2);
found=zeros(2,num_points);

linesP=DoLineEqsFromPoints(P,con);

%calculate intersections
for i=1:num_points
    closest=zeros(2);
    closest_dist=99999;
    for j=1:num_bis
        %calculate all intersections of current bisector with all lines, then get the minimum distance
        if ((bis(i,j,1)==0) & (bis(i,j,2)==0))%no more bisectors for this point
            break;
        end
        %loop through top part of line-matrix
        for k=1:num_corners-1
            for l=k+1:num_corners
                if(linesP(k,l,1)~=0) %there is a line between these points
                    int=CalculateIntersection(bis(i,j,:),linesP(k,l,:));
                    %check if it is inside the segment
                    inside=LimitIntersection(int,P(:,k),P(:,l));
                    if (inside)
                        %calculate distance between int and X(:,i)
                        dist=sqrt((int(1)-X(1,i))^2 + (int(2)-X(2,i))^2);
                        if (dist<closest_dist)
                            closest_dist=dist;
                            closest=int;
                        end
                    end
                end
            end
        end % for - lines
    end % for - bisectors
    found(1,i)=closest(1); found(2,i)=closest(2);
end % for - points

% checking in both coordinates (knowing that they must be colinear, checking for one coordinate is enough)
function inside = LimitIntersection(x,p1,p2)
inside=0;
if ((x(1)<=max(p1(1),p2(1))) & (x(1)>=min(p1(1),p2(1))) & (x(2)<=max(p1(2),p2(2))) & (x(2)>=min(p1(2),p2(2))) )
    inside=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bisectors (DoBisector,Perpendiculars)

function eq=DoBisector(P,X1,X2)
%vectors
v1=X1-P;
v2=X2-P;
%normalise
v1=v1/(sqrt(v1(1)*v1(1)+v1(2)*v1(2)));
v2=v2/(sqrt(v2(1)*v2(1)+v2(2)*v2(2)));
%bisector vector
v=(v1+v2)/2;
if (v==[0 0]') %opposite vectors
    eq=DoLine([-v1(2) v1(1)],P); %perpendicular vector
else
    eq=DoLine(v,P);
end

function added=Perpendiculars(X,con,look_up,bisectors) %maybe dont need con
global NUM_MIDDLE_POINTS;
added=bisectors;
%get number of corners ... from con or bisectors or just pass it as parameter
szb=size(bisectors);
num_corners=szb(1);
szx=size(X);
num_points=szx(2);
for i=num_corners+1:num_points
    con_num=ceil((i-num_corners)/NUM_MIDDLE_POINTS);
    added(i,:,:)=0; %initialisation
    vector=X(:,look_up(con_num,1)) - X(:,look_up(con_num,2)); %vector between the corner points (more robust than using the middle point)
    pvector(1)=-vector(2); %perpendicular vector
    pvector(2)=vector(1);
    added(i,1,:)=DoLine(pvector,X(:,i));
end


% option 2 (all possible combinations of adjacent connections to form bisectors)
% if '0' in con, no more connections for that point
% we want also the bisector of the last connection with the first (except when only 2 connections)
% % % % 'bisectors' structure:
% % % % bisectors(num_corners, num_bisectors (with the '0' at the end if no more), bisector equation (ax + by + c = 0)
% PRE: minimum 2 connections
function bisectors=Bisectors(X,con)
sz = size(con);
num_corners=sz(1);
num_con=sz(2);
if (num_con==2)
    bisectors=zeros(num_corners,1,3);
else
    bisectors=zeros(num_corners,num_con,3);
end

for i=1:num_corners
    do_last_one=1; %TRUE
    for j=2:num_con
        if (con(i,j)==0) %no more connections
            if(j==3) % only 2 connections
                do_last_one=0; %FALSE
            elseif (j<3) % j==2, something wrong
                fprintf('ERROR: only 1 connection !!???');
                do_last_one=0;
            else % 3 or more connections, do last bisector
                j=j-1; %to index correctly when doing the last one
                break;
            end
        else
            bisectors(i,j-1,:)=DoBisector(X(:,i),X(:,con(i,j-1)),X(:,con(i,j)));
        end % if - connections
    end %for - connections
    %do last one if necessary
    if ((do_last_one==1) & (num_con>2))
        bisectors(i,j,:)=DoBisector( X(:,i), X(:,con(i,j)), X(:,con(i,1)) );
    end
end %for - points

% option 2 (all possible combinations of adjacent connections to form bisectors) + OWN LINES
% if '0' in con, no more connections for that point
% we want also the bisector of the last connection with the first (except when only 2 connections)
% % % % 'bisectors' structure:
% % % % bisectors(num_points, num_bisectors (with the '0' at the end if no more), bisector equation (ax + by + c = 0)
% PRE: minimum 2 connections
function bisectors=Bisectors2(X,con)
sz = size(con);
num_points=sz(1);
num_con=sz(2);
if (num_con==2)
    bisectors=zeros(num_points,1 + num_con,3);
    last_bisector_index=1;
else
    bisectors=zeros(num_points,num_con*2,3);
    last_bisector_index=num_con;
end

for i=1:num_points
    do_last_one=1; %TRUE
    bisectors(i,last_bisector_index+1,:)=DoLine2Points(X(:,i),X(:,con(i,1))); %add first own line
    for j=2:num_con
        if (con(i,j)==0) %no more connections
            if(j==3) % only 2 connections
                do_last_one=0; %FALSE
            elseif (j<3) % j==2, something wrong
                fprintf('ERROR: only 1 connection !!???');
                do_last_one=0;
            else % 3 or more connections, do last bisector
                j=j-1; %to index correctly when doing the last one
                break;
            end
        else
            bisectors(i,j-1,:)=DoBisector(X(:,i),X(:,con(i,j-1)),X(:,con(i,j)));
            % add own lines
            bisectors(i,last_bisector_index+j,:)=DoLine2Points(X(:,i),X(:,con(i,j)));
        end % if - connections
    end %for - connections
    %do last one if necessary
    if ((do_last_one==1) & (num_con>2))
        bisectors(i,j,:)=DoBisector( X(:,i), X(:,con(i,j)), X(:,con(i,1)) );
    end
    %now the own lines
    
end %for - points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of FindQuad- Bisectors (DoBisector) , Intersections (IntersectionImage,CalculateIntersection,DoLine2Points,DoLine)
%*****************************************************************************************************

%------------------------------------------
%   Begining of "AffineTrans" and "AddNoise" function 
%-----------------------------------------

function [Xout,A]=AffineTrans(Xin,a,b,c,d,t1,t2)
if ((a==0) & (d==0)) 
    %define the transformation randomly
    %A=rand(3,3)*2 - 1; %random numbers from -1 .. 1
    %A=rand(3,3)*4 - 2; %random numbers from -2 .. 2
    A=rand(3,3) - 0.5; %random numbers from -0.5 .. 0.5
    A(1,3)=0; A(2,3)=0; %no translation
    A(3,1)=0; A(3,2)=0; %no scaling
    A(3,3)=1;
    A(1,1)=A(1,1)+1; A(2,2)=A(2,2)+1; %from 0.5 .. 1.5
else
    %use the one given
    A=zeros(3,3);
    A(1,1)=a;A(1,2)=b;A(2,1)=c;A(2,2)=d;A(1,3)=t1;A(2,3)=t2;A(3,3)=1;
end
szx = size(Xin);
m = szx(2);
Xout=zeros(2,m);
for i=1:m
    w= A(3,1)*Xin(1,i) + A(3,2)*Xin(2,i) + A(3,3);
    Xout(1,i)= (A(1,1)*Xin(1,i) + A(1,2)*Xin(2,i) + A(1,3))/w;
    Xout(2,i)= (A(2,1)*Xin(1,i) + A(2,2)*Xin(2,i) + A(2,3))/w;
end
%why not matrix form? Xout=A*Xin;

function X=AddNoise(X)
szx=size(X);
m=szx(2);
noise=rand(2,m)*10-5;
X=X+noise;

%------------------------------------------
%   End of "AffineTrans" and "AddNoise" function 
%-----------------------------------------

%------------------------------------------
%   Begining of "Image2Matrix" function - FER
%-----------------------------------------
% Calculate centre & scale, but results in image space
function [Xout,Points,centre,scale]=Image2Matrix(I,MeanRefDist)

szx = size(I);
m = szx(1);
n = szx(2);
totX=0; %to calculate centre
totY=0; %to calculate centre
[xp,yp]=find(I==150); % look for the 255 & 150 (sparse)
[xc,yc]=find(I==255); % look for the 255 & 150 (sparse)
x=[xc;xp];
y=[yc;yp];

num_pixels=length(x);
centre(1)=sum(x)/num_pixels;
centre(2)=sum(y)/num_pixels;

Xout=[x y]';
%centre to calculate scale
Xaux=[x-centre(1) y-centre(2)]'; %centre
Points=[xc yc]'; % Points=[xc-centre(1) yc-centre(2)]'; %centre

% calculate scale 
dist=sum(sqrt(diag(Xaux'*Xaux)));
scale = MeanRefDist*num_pixels/dist;

% % scale calculated, now scale
% Xout=Xout*scale;
% Points=Points*scale;

%------------------------------------------
%   End of "Image2Matrix" function - FER
%-----------------------------------------

%------------------------------------------
%   Begining of "DrawLine" function - FER
%-----------------------------------------
function Iout=DrawLine(p1,p2,Iin)

Iout=Iin;
eq=DoLine2Points(p1,p2);

%%this gives better results, but slower:
for i = min(p1(1),p2(1)):0.01:max(p1(1),p2(1))
    for j = min(p1(2),p2(2)):0.01:max(p1(2),p2(2))
        if ( round(eq(1)*i+eq(2)*j+eq(3)) == 0 )  
            
            % % this is faster
            % for i = min(p1(1),p2(1)):max(p1(1),p2(1))
            %     for j = min(p1(2),p2(2)):max(p1(2),p2(2))
            %         if ( abs(eq(1)*i+eq(2)*j+eq(3)) < 20 )     
            
            Iout(1,round(i),round(j))=150; %FER - CHANGE - size of image/2
        end
    end
end
% draw points - for viewing
Iout(1,round(p1(1)),round(p1(2)))=255;
Iout(1,round(p2(1)),round(p2(2)))=255;     

function Iout=DrawPoints(X,Iin,intensity)
Iout=Iin;
sz=size(X);
for i=1:sz(2)
    % draw points - for viewing
    Iout(round(X(1,i)),round(X(2,i)))=intensity;
    Iout(round(X(1,i)),round(X(2,i)))=intensity;   
end

%--------------------------------------
%   End of "DrawLine" function - FER
%--------------------------------------


%--------------------------------------------------
%   Begining of CENTRING AND SCALING FUNCTIONS
%-------------------------------------------------
function [Xout,centre]=Centre(Xin)
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

for j = 1:n
    Xout(:,j) = Xin(:,j) - centre';
end

function [Xout, Centers] = CentreThis(Xin)
szx = size(Xin);
m = szx(3);
for i = 1:m
    [Xout(:,:,i),Centers(i,:)]=Centre(Xin(:,:,i));
end

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

function [Xout,scale]=Scale(Xin,MeanRefDist)
szx = size(Xin);
n = szx(2);
r = 0;
dist = 0;
for i = 1:n
    dist = dist + sqrt((Xin(1,i)*Xin(1,i)) + (Xin(2,i)*Xin(2,i)));
end
scale = MeanRefDist*n/dist; 
for i = 1:n
    Xout(1,i) = scale*Xin(1,i);
    Xout(2,i) = scale*Xin(2,i);
end

function [Xout, scales,MeanRefDist] = ReScale(Xin)
Xout = Xin;
szx = size(Xin);
n = szx(2);
m = szx(3);
RefDist = 0;
for i = 1:n
    RefDist = RefDist + sqrt((Xin(1,i,1)*Xin(1,i,1)) + (Xin(2,i,1)*Xin(2,i,1)));
end
MeanRefDist=RefDist/n; 
scales(1) = 1;
for j = 2:m  
    [Xout(:,:,j),scales(j)]=Scale(Xin(:,:,j),MeanRefDist);
end

function Xout=UnCentre(Xin,centre)
Xout=Xin;
Xout(1,:) = Xin(1,:) + centre(1);
Xout(2,:) = Xin(2,:) + centre(2);

function Xout=UnScale(Xin, scale)
Xout=Xin;
Xout(1,:) = Xin(1,:)/scale;
Xout(2,:) = Xin(2,:)/scale;

%--------------------------------------------------
%   End of CENTRING AND SCALING FUNCTIONS
%-------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%% 2:1 AFFINE FUNCTIONS %%%%%%%%%%%%%%%% CATT2Affine, getNewXAffine, GetNewab
function [a,b]=CATT2Affine(T)
aux=T(1)+T(4);
bux=T(5)+T(8);
cux=T(2)+T(3);
dux=T(7)+T(6);
a(1)=T(12)*bux - T(11)*dux;
a(2)=T(11)*bux - T(12)*dux;
a(3)=T(10)*bux - T(9)*dux;
a(4)=T(9)*bux - T(10)*dux;
b(1)=T(11)*cux - T(12)*aux;
b(2)=T(12)*cux - T(11)*aux;
b(3)=T(9)*cux - T(10)*aux;
b(4)=T(10)*cux - T(9)*aux;

function X=GetNewXAffine(basis1,basis2,a,b)
% szb=size(basis1); delete
% num_points=szb(2); delete
x1=[a(1) a(2)]*basis1+[a(3) a(4)]*basis2;
x2=[b(1) b(2)]*basis1+[b(3) b(4)]*basis2;
X=[x1;x2];

function [a,b]=GetNewab(X1,X2,X)
Y=[X1;X2];
AB=X*pinv(Y);
a=AB(1,:);
b=AB(2,:);
