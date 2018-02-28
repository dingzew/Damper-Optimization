% This script belongs to University of British Columbia
% Smart Structural Team. Please don't use it for any object other than study

% Author: Dingze Wang
% Finished on May 31, 2017

%-------------------------------------------------------

%set parapmeter
%please make sure in this area
%The range of the thickness should be 10
%The range of the width is 100


Fslip=15;                  %friction force
length=250;                %distance between bolt to bolt
dbolt=22.225;              %diameter of bolt
dpad=100;                  %diameter of the friction pad
t_min=5;                   %minimum thickness of the strip
t_max=15;                  %maximum thickness of the strip
w_min=50;                  %minimum width of the strip
w_max=150;                 %maximum width of the strip
miu = 0.4;                 %The friction coefficient
alpha = pi/4;              %Half of the initial angle between two strip

% set the range
t1=linspace(t_min,t_max,100);
w1=linspace(w_min,w_max,100);

%create meshgrid 
[t,w]=meshgrid(t1,w1);

%calculate buckling and yielding force from equation, you may add safety
%factor after the equation
Pcr=pi^2 * 200 * 1000 /12 .* t .^3 .* w /(4*length^2)/1000 /sqrt(2) * 2 ;
Fy=sqrt(248.22^2 ./(( 6 * cos(alpha)* length  ./ t ./ w ./ w + sin(alpha) ./ w ./ t).^2 + 27 / 4 * ( cos(alpha) ./ w ./ t ).^2 ))/1000  ;
 
% initialize Fslip matrix
for i=1:100
    for j=1:100
        m(i,j)=Fslip;
    end
end
 
 %-------------------------------------------------------------------
 %This part draw buckling yielding and slip force in one 3D image
 
 subplot(3,2,1)
 
 %the last parameter is for color print, no actual meaning
 
 mesh(t,w,Pcr,Pcr);
 hold on;
 mesh(t,w,Fy,Pcr-100);
 hold on;
 mesh(t,w,m,Pcr+100);
 hold on;
 xlabel('Thickness(mm)');
 ylabel('width£¨mm)');
 zlabel('Vol(mm3)');
 grid on;
 
 %--------------------------------------------------------------------
 %here is for selecting feasible area (Flag(i,j)==1)
 
 subplot(3,2,2)
 
 %selecting feasible area
 %Flag the matrix
 %represent on 2D plot
 for i=1:100
     for j=1:100
         if (Pcr(i,j)>Fslip && Fy(i,j)>Fslip) %The buckling and yielding force should be above the Slip force
             Flag(i,j)=1; %Means it's feasible on this point
             scatter(t_min+10*j/100,w_min+100*i/100); %mark this point on the 2D plot
             hold on;
         else Flag(i,j)=0; %Means this point yielding or buckling or both
         end
         hold on;
     end
     hold on;
 end
 
 xlabel('Thickness(mm)');
 ylabel('width(mm)');
 axis([t_min,t_max,w_min,w_max]); %set the plot in range, delete the white area
 grid on;
 
 %---------------------------------------------------------------
 %Safety check, only for strip strenth,
 %If not strong enough, terminate the program
 %If that happens you may go back to change the range of width and
 %thickness, this depends on whether there's a Flag(i,j)==1?
 
 temp = 0;
 for i=1 :100
     for j=1:100
         if (Flag(i,j)==1 )
             temp = 1;
         end
     end
 end
 
 if(temp == 0) %in this case, no Flag is 1
     disp('the strip is not strong enough, you may change the range');
     return    %here the program is terminated
 end

%--------------------------------------------------------------
%this module is to compare different point in feasible area
%here we will get the minimum volume in feasible area
%And it will also return the dimension of this critical point

Volmin=10e8; %set a initial min value randomly, as it must be change in the future
for i=1:100
    for j=1:100
        if Flag(i,j)==1
            Vol(i,j)=(t_min+10*j/100) * (w_min+100*i/100) * length;
              if Vol(i,j)<Volmin
                  Volmin=Vol(i,j); imin=i; jmin=j;
              end
        else
            Vol(i,j)=0;
        end
    end
end

%this is for data check because imin and jmin will be overlapped
%----------------------------------------------------------------
%iminmemo=imin
%jminmemo=jmin
%output the data
%-----------------------------------------------------------------

thickness = t_min+10*jmin/100; %calculate the thickness using iterater
width = w_min+100*imin/100;    %calculate the width using iterater
volmin=Vol(imin,jmin);         %retract the minimum volume from the dataset

disp ('Thickness(original):');
disp(thickness);
disp('Width(original):');
disp(width);
disp('Volume(original):')
disp(volmin);

%------------------------------------------------------
%plot Volume

subplot(3,2,3)
c=[1,0,0]; %this color control is for the red point on the critical point
meshc(t,w,Vol); %draw the volume distribution in 3D plot
hold on;
scatter3(t_min+10*jmin/100,w_min+100*imin/100,Vol(imin,jmin),15,c); %emphsize the critical point on 3D plot
xlabel('Thickness(mm)');
ylabel('width£¨mm)');
zlabel('Vol(mm3)');
grid on;

%-----------------------------------------------------------------
%This method modified the thickness is not right, however it's a good template 
%The thickness is only available in several inch( eg 5/16'')

%set the value as the ceiling value
% remin = mod( thickness,0.5 );
% if(remin == 0)
%     thickness = thickness;
% else thickness = thickness + (0.5 - remin);
% end
% 
% remin=0;
% 
% remin = mod( width , 5 );
% if(remin == 0)
%     width = width;
% else width = width + (5 - remin);
% end

%------------------------------------------------
% Here we modified the thickness to availble value
% Please make sure the maximax thickness no more than 20
% because I didn't set function for more than 20mm thickness 

if(thickness>4.765 & thickness<6.35)
    thickness = 6.35;
    i=1;
    temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(6.35-t_min))) == 1) %To find the smallest feasible width with a given thickness
            temp=1; %once find, make the flag 1 to jump out the loop
        end
        i=i+1; %just iteritor
    end
    width = w_min + i; %once jump out the loop, calculate the result with iteritor
end

%-----------------------------------------------------
%-----------------------------------------------------
if(thickness>6.35 & thickness<7.9375)
    
    thickness1 = 6.35;
    i=1;
    temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(6.35-t_min))) == 1) %retrack the new width on the feasible area
            temp=1;
        end
        i=i+1;
    end
    width1 = w_min + i; %width1 is the width modified by the first time
    vol1=thickness1*width1*length; %vol1 is just used just for comparasion
    if(i==100) %This case will not get the data
        vol1=10e8;  %so set it as a big number in order to make it lose in the later comparision
    end;
    
    thickness2 = 7.93;
    i=1;
    temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(7.93-t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width2 = w_min + i; %width2 is the width modified by the first time
    vol2=thickness2*width2*length; %vol1 is just used just for comparasion
    
    if(vol1<vol2) %compare two colume and decide which dimension to choose
        thickness=thickness1;
        width=width1;
    else
        thickness=thickness2;
        width=width2;
    end
end

%------------------------------------------------------
%------------------------------------------------------
%the down part are simply repeat the first part

if(thickness>7.9375 & thickness<9.525)
    
    thickness1 = 7.9375;
    i=1;
    temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(7.93-t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width1 = w_min + i;
    vol1=thickness1*width1*length;
    if(i==100)
        vol1=10e8;
    end
    
    thickness2 = 9.52;
    i=1;
    temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(9.52-t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width2 = w_min + i;
    vol2=thickness2*width2*length;
    
    if(vol1<vol2)
        thickness=thickness1;
        width=width1;
    else
        thickness=thickness2;
        width=width2;
    end
end

%-----------------------------------------------
%-----------------------------------------------

if(thickness>9.525 & thickness<11.1125)
    
    thickness1 = 9.525;
    i=1;
    temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(9.52-t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width1 = w_min + i;
    vol1=thickness1*width1*length;
    if(i==100)
        vol1=10e8;
    end
    
    thickness2 = 11.1125;
    i=1;
    temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(11.11-t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width2 = w_min + i;
    vol2=thickness2*width2*length;
    
    if(vol1<vol2)
        thickness=thickness1;
        width=width1;
    else
        thickness=thickness2;
        width=width2;
    end
end

%-----------------------------------------------
%-----------------------------------------------

if(thickness>11.1125 & thickness<12.7)
    
    thickness1 = 11.1125;
    i=1;
    temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(11.11-t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width1 = w_min + i;
    vol1=thickness1*width1*length;
    if(i==100)
        vol1=10e8;
    end
    
    thickness2 = 12.7;
    i=1;temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(12.7-t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width2 = w_min + i;
    vol2=thickness2*width2*length;
    
    if(vol1<vol2)
        thickness=thickness1;
        width=width1;
    else
        thickness=thickness2;
        width=width2;
    end
end

%-----------------------------------------------
%-----------------------------------------------

if(thickness>12.7 & thickness<14.29)
    
    thickness1 = 12.7;
    i=1;
    temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(12.7-t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width1 = w_min + i;
    vol1=thickness1*width1*length;
    if(i==100)
        vol1=10e8;
    end
    
    thickness2 = 14.29;
    i=1;
    temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(14.29-t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width2 = w_min + i;
    vol2=thickness2*width2*length;
    
    if(vol1<vol2)
        thickness=thickness1;
        width=width1;
    else
        thickness=thickness2;
        width=width2;
    end
end

%-----------------------------------------------
%-----------------------------------------------

if(thickness>14.29 & thickness<15.875)
    
    thickness1 = 14.29;
    i=1;
    temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(14.29-t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width1 = w_min + i;
    vol1=thickness1*width1*length;
    if(i==100)
        vol1=10e8;
    end
    
    thickness2 = 15.875;
    i=1;temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(15.875-t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width2 = w_min + i;
    vol2=thickness2*width2*length;
    
    if(vol1<vol2)
        thickness=thickness1;
        width=width1;
    else
        thickness=thickness2;
        width=width2;
    end
end

%-----------------------------------------------
%-----------------------------------------------

if(thickness>15.875 & thickness<17.46)
    
    thickness1 = 15.875;
    i=1;temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(15.875-t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width1 = w_min + i;
    vol1=thickness1*width1*length;
    if(i==100)
        vol1=10e8;
    end
    
    thickness2 = 17.46;
    i=1;
    temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(17.46-t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width2 = w_min + i;
    vol2=thickness2*width2*length;
    
    if(vol1<vol2)
        thickness=thickness1;
        width=width1;
    else
        thickness=thickness2;
        width=width2;
    end
end

%-----------------------------------------------
%-----------------------------------------------

if(thickness>17.46 & thickness<19.05)
    
    thickness1 = 17.46;
    i=1;temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(17.46 - t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width1 = w_min + i;
    vol1=thickness1*width1*length;
    if(i==100) vol1=10e8; end;
    
    thickness2 = 19.05;
    i=1;temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(19.05 - t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width2 = w_min + i;
    vol2=thickness2*width2*length;
    
    if(vol1<vol2)
        thickness=thickness1;
        width=width1;
    else
        thickness=thickness2;
        width=width2;
    end
end

%-----------------------------------------------
%-----------------------------------------------

if(thickness>19.05 & thickness<20.64)
    
    thickness1 = 19.05;
    i=1;temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(19.05 - t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width1 = w_min + i;
    vol1=thickness1*width1*length;
    if(i==100)
        vol1=10e8;
    end
    
    thickness2 = 20.64;
    i=1;temp=0;
    while((i<100) & (temp == 0))
        if (Flag(i,ceil(10*(20.64 - t_min))) == 1)
            temp=1;
        end
        i=i+1;
    end
    width2 = w_min + i;
    vol2=thickness2*width2*length;
    
    if(vol1<vol2)
        thickness=thickness1;
        width=width1;
    else
        thickness=thickness2;
        width=width2;
    end
end

%-----------------------------------------------
%this is the end of comparision part
%-----------------------------------------------


%output again
%This is after modified the thickness to available size

disp ('Thickness(modified 1st):');
disp(thickness);
disp('Width(modified 1st):');
disp(width);
disp('Volume(modified 1st):')
volafter=thickness*width*length;
disp(volafter);

%----------------------------------------------
%here is to modified the width to a %5=0 number
remin=0;
remin = mod( width , 5 );
if(remin == 0)
    width = width;
else width = width + (5 - remin);
end

disp ('Thickness (modified 2nd):');
disp(thickness);
disp('Width (modified 2nd):');
disp(width);
disp('Volume (modified 2nd):')
volafter=thickness*width*length;
disp(volafter); %volafter is the volume after modified 2 times
%------------------------------------------------


disp('-------------------------------');
disp('Here we will examine the modified version (1)');
radius = width/10/2; %here assuming the diameter of friction pad equals to the width of the strip
%However, must attention that it may cause problem in the later cutting
moment = Fslip * length * cos(alpha);
P = moment/ (2/3) /miu / radius/ 10 / 2 ; %clamping force. May be used in abaqus model to set pressure
%Attention:there are two friction pads
radius1=dpad / 2 / 10;
P1 = moment/ (2/3) /miu / radius1/ 10 / 2 ;
Pressure = P1 * 1000 / ((dpad/10/2)^2*pi)
V=Fslip; %shear force
Vr= 0.6 * 0.8 * 2 * (dbolt/2)^2 * pi * 825; %shear resistance
Tr= 495 * (dbolt/2)^2 * pi /1000; %vertical force resistance

Pmax = pi * radius^2 * 3102 / 1000; %strength of the friction pad
Pbolt = 495 * (dbolt/2)^2 * pi /1000; %strengh of the bolt
 
if(Pmax<Pbolt)
    disp('it depends on friction pad strength')
    if (Pmax > P)
        disp(Pmax);disp('>'); disp(P)
        disp('The specimen passed the test' );
    end
end

if(Pmax>Pbolt)
    disp('it depends on bolt strength')
    Pmax=Pbolt;
    
    if ((P/(Tr))^2+(Fslip/(Vr/1000))^2 < 1)
        disp((P/(Tr))^2+(Fslip/(Vr/1000))^2);disp('<'); disp(1)
        disp('The specimen pass the test' );
    else
        disp('The bolt is not strong enough');
        dbolt_r = 2 * (((P*1000/495)^2+(Fslip*1000 /0.6 / 0.8 / 2 / 825)^2)/pi/pi)^(1/4); 
        %Assuming equation is 1, calculate the bolt size
        dbolt_inch = dbolt_r /25.4; %transfer to inch
        %modified to n/16 inch to make it available
        dbolt_modified = floor(dbolt_inch * 16) + 1; 
        dbolt_modified =  dbolt_modified / 16; 
        dbolt_modified =  dbolt_modified * 25.4; %back to mm
        disp('the bolt should no less than:')
        disp(dbolt_modified);
        min_b = dbolt_modified * 1.75 * 2; %the bolt edge requirement
        if ( width < min_b )
            disp ( 'the width is too small, Im going to make the width bigger');
            width =min_b;
        else
            disp('the width is big enough, Just make the bolt bigger');
        end
        %modified Vr and Tr
        Vr_modified = 0.6 * 0.8 * 2 * (dbolt_modified/2)^2 * pi * 825;
        Tr_modified = 495 * (dbolt_modified/2)^2 * pi /1000;
        disp('In this way the equation result will be');
        %check the number whether it's lower than 1
        Aim_1=(P/(Tr_modified))^2+(Fslip/(Vr_modified/1000))^2;
        disp('the modified check number is :');
        disp(Aim_1);
    end
end
 

 
 %------------------------------------------------------------------------
 %This is the end of rectangular optimization
 %start to cut
 %set parameters


b0 = width; %longer side
thetamax = atan (b0 /2 / length); %calculate maximum cutting angle
theta=linspace(0, thetamax, 100);
% This 100*100 matrix is for explore 2 parameters, since now only explore 1
% Maybe later there will be needed to explore 2 para problem in tappered section 
% ---------------------------------------------------------------------------
% [t,theta]=meshgrid(t1,theta1);
% for i=1:100
%     for j=1:100
%         t(i,j)=thickness;
%     end
% end
%-------------------------------------

b1 = b0 - 2 * tan (theta) * length ; %The upper side

% Calculate the yielding forcre for all the cross-section
for i=1:100 %i controls the theta to be cut
    temp=b1(i);
    w=linspace(temp+0.001, b0, 100);
    Fy1(i)=1e8;
    for k=1:100 %k controls the distance from the top to the cross-section
        side(k) =  ( w(k) - temp ) /2;
        dis(k) = side(k) / tan (theta(i));
        Fy2(k) = sqrt( 248.22 ^ 2 /(( 6 * cos(alpha) * dis(k)  / thickness / w(k) / w(k) + sin(alpha) / w(k) / thickness)^2 + 27 / 8 * ( cos(alpha) / w(k) / thickness )^2 ) )/1000 ;
        if(Fy2(k)<Fy1(i))
            Fy1(i)= Fy2(k); 
            FyminWidth(i)=w(k);
            %Pick up the most easily buckling cross-section
        end
    end
    % FyminWidth
end

I0= 1 / 12 * b0 * thickness^ 3;
bL= 2 * tan(theta) / b0 * length;
para = -0.3333.*(bL).*(bL).*(bL) + 0.15.*(bL).*(bL)-0.8117.*(bL)+2.473; %got in curve fitting tool box
%here just remember the original data, using 3 times polynomial
%bL= [0.1 0.3 0.5 0.7 0.9];
%para=[2.393 2.235 2.062 1.865 1.621];
Pcr1 = para * 200 * 1000 * I0 /(length^2) / 1000 /sqrt(2) *2;

% set the range
%---------------------draw the rectangular specimen
subplot(3,2,4)

%here is for visional coding
%This part is for rectangular visional coding
%------------------------------------------------------------------------
%set the track to draw an cuboid
x=[0 0 width width 0 0 0 0 0 0 width width 0 width width width width];
y=[0 0 0 0 0 length length 0 0 length length length length length 0 0 length];
z=[0 thickness thickness 0 0 0 thickness thickness 0 0 0 thickness thickness thickness thickness 0 0];

%---------------this part color the image
%color every face
c = [0.7 0.7 0.7]; %set color grey
plot3(x,y,z,'k');
%the four point of each rectangular
fill3([0 0 width width],[0 0 0 0],[0 thickness thickness 0],c);
hold on;
fill3([0 0 0 0],[length length 0 0],[0 thickness thickness 0],c);
hold on;
fill3([0 width width 0],[length length length length],[0 0 thickness thickness],c);
hold on;
fill3([width width width width],[length 0 0 length],[thickness thickness 0 0],c);
hold on;
fill3([0 0 width width],[0 length length 0],[thickness thickness thickness thickness],c);
hold on;
fill3([0 0 width width],[0 length length 0],[0 0 0 0],c);
hold on;
%-------------------------------------------------------------

axis([0,width,0,length,0,thickness]); %delete white area
axis equal;

%-------------------------------------------------
%-------------------------------------------------------------------
 
 subplot(3,2,5)
 
%-------------------------------------
%no longer use 2D matrix
%  mesh(t,theta,Pcr,Pcr);
%  hold on;
%  mesh(t,theta,Fy,Pcr-100);
%  hold on;
%  mesh(t,theta,m,Pcr+100);
%  hold on;
%--------------------------------------
for i=1:100
         m(i)=Fslip;
  end
%This is for 2D drawing
theta=theta/pi*180; %tansfer theta into degree units
plot(theta,Pcr1,'k');
hold on;
plot(theta,Fy1,'r'); %Yielding force
hold on;
plot(theta,m,'y'); %Buckling force
hold on;

legend('Buckling','Yielding','Fslip');
xlabel('theta£¨o)');
ylabel('Force(KN)');
grid on;

%------------------------------------------------
%------------------------------------------------
%initialize Flag1(i) which is different from Flag as a 100*100 matrix
for i=1:100
    Flag1(i)=1;
end

%use matrix to mark the buckling and yielding condition, if fail keep 0
%The 100*1 matrix cannot overlap 100*100 matrix so Flag1 is re-defined
%if not fail, set Flag as 1
for i=1:100
    if (Pcr1(i)>Fslip && Fy1(i)>Fslip)
        Flag1(i)=0;
    end
end

temp = 0;

%here retract the maximum theta that can be cut
i=2;
%because i=1 is a strange value and should not be considerated
while((i<100) & (temp == 0))
    if (Flag1(i) == 1)
        temp=1;
        %once find, set temp as 1 in order to jump out the loop
        imin=i;
    end
    i=i+1;
    %i++ the iterator
end

thetamin = thetamax*imin/100;

a1=b0-2*length*tan(thetamin);
V_56 = 0.6 * 0.8 * 2 * 1 *7.94 ^2 * pi * 825 /1000; %Vr of (5/8)'' diameter bolt in the short end
%--------------------------------------------------------------------
%check the short side bolt shear strength
if(V_56 < Fslip)
    disp('the bolt size needs to be changed');
    %use Vr as the Fslip to get the minimum bolt size in the short end
    dbolt1 = 2 * sqrt(Fslip * 1000 /(0.6 * 0.8 * 2 * 1 * pi * 825));
    disp('the new diameter of bolt should no less than');
    disp(dbolt1);
    %------------------------------
    %here will modified the diameter to available value
    dbolt1 = dbolt1 /25.4;
    dbolt1 = ceil(dbolt1 * 16);
    dbolt1 = dbolt1 / 16;
    dbolt1 = dbolt1 * 25.4;
    edge_min = 1.75 * 2 * dbolt1;
    %-----------------------------------------------
    %here will check the shear capacity of the strip
    Br=3 * 0.8 * 1 * thickness * dbolt1 * 825 /1000;
    if(Br<Fslip)
        disp('strip will fail(Br)');
    end
    %check end
    %------------------------------------------------
    if(a1<edge_min) %here we may cut too much
        disp('the upper edge depends on bolt size')
        Volmin=(edge_min+b0)*length/2*thickness;
        %here because of changing bolt, the edge should not be 56
        disp('theta(original):');
        disp(atan((b0-edge_min)/2/length)/pi*180);
        disp('Volume(original):');
        disp(Volmin);
        subplot(3,2,6)
        
        a0=width;
        a1=edge_min;
        %-----------------------------visualize coding
        s=(a0-edge_min)/2;
        x=[0 0 a0 a0 0 0 s s 0 s s+a1 s+a1 s s+a1 a0 a0 s+a1];
        y=[0 0 0 0 0 0 length length 0 length length length length length 0 0 length];
        z=[0 thickness thickness 0 0 thickness thickness 0 0 0 0 thickness thickness thickness thickness 0 0];
        plot3(x,y,z,'k');
        axis([0,width,0,length,0,thickness]);
        axis equal;
        %--------------this part is to color the image
        fill3([0 0 a0 a0],[0 0 0 0],[0 thickness thickness 0],c);
        hold on;
        fill3([0 s s 0],[0 length length 0],[thickness thickness 0 0],c);
        hold on;
        fill3([0 s s 0],[0 length length 0],[thickness thickness 0 0],c);
        hold on;
        fill3([s s+a1 s+a1 s],[length length length length],[0 0 thickness thickness],c);
        hold on;
        fill3([s+a1 a0 a0 s+a1],[length 0 0 length],[thickness thickness 0 0],c);
        hold on;
        fill3([0 s s+a1 a0],[0 length length 0],[thickness thickness thickness thickness],c);
        hold on;
        fill3([0 s s+a1 a0],[0 length length 0],[0 0 0 0],c);
        hold on;
        
        axis([0,width,0,length,0,thickness]);
        axis equal;
        
        disp('the up side is ');
        disp(a1);
    end
    %---------------------------------------
    
    if(a1>edge_min) %just use theta which we calculate to cut
        Volmin=(2*b0-2*length*tan(thetamin)) * length / 2 * thickness;
        
        disp('theta(original):');
        disp(thetamin/pi*180);
        disp('Volume(original):')
        disp(Volmin);
        
        subplot(3,2,6)
        a0=width;
        a1=b0-2*length*tan(thetamin);
        s=(a0-a1)/2;
        %visualize coding-------------------------------------------------
        x=[0 0 a0 a0 0 0 s s 0 s s+a1 s+a1 s s+a1 a0 a0 s+a1];
        y=[0 0 0 0 0 0 length length 0 length length length length length 0 0 length];
        z=[0 thickness thickness 0 0 thickness thickness 0 0 0 0 thickness thickness thickness thickness 0 0];
        plot3(x,y,z,'k');
        axis([0,width,0,length,0,thickness]);
        axis equal;
        
        %--------------this part is to color the image
        fill3([0 0 a0 a0],[0 0 0 0],[0 thickness thickness 0],c);
        hold on;
        fill3([0 s s 0],[0 length length 0],[thickness thickness 0 0],c);
        hold on;
        fill3([0 s s 0],[0 length length 0],[thickness thickness 0 0],c);
        hold on;
        fill3([s s+a1 s+a1 s],[length length length length],[0 0 thickness thickness],c);
        hold on;
        fill3([s+a1 a0 a0 s+a1],[length 0 0 length],[thickness thickness 0 0],c);
        hold on;
        fill3([0 s s+a1 a0],[0 length length 0],[thickness thickness thickness thickness],c);
        hold on;
        fill3([0 s s+a1 a0],[0 length length 0],[0 0 0 0],c);
        hold on;
        
        axis([0,width,0,length,0,thickness]);
        axis equal;
        
        disp('the up side is ');
        disp(a1);
    end
    %-------------------------------------------------
    
    
elseif(V_56 > Fslip && a1 < 56) %(5/8)'' diameter bolt is enough and we cut too much
    %not meet the steel structure code
   %These two data has no actual meaning, just for test
    width_memo = a1; 
    angle_memo = atan((b0-a1)/2/length)/pi*180;
    disp('the upper edge depends on bolt size')
    Volmin=(56+b0)*length/2*thickness;
    disp('theta(original):');
    disp(atan((b0-56)/2/length)/pi*180);
    disp('Volume(original):')
    disp(Volmin);
    subplot(3,2,6)
    %------------------------------------------------
    %always check the Br
    Br=3 * 0.8 * 1 * thickness * 15.875 * 825 /1000;
    if(Br<Fslip)
        disp('strip will fail(Br)');
    end
    %----------------------------------------------
    %initialize
    a0=width;
    a1=56;
    s=(a0-56)/2;
    %visualize coding
    x=[0 0 a0 a0 0 0 s s 0 s s+a1 s+a1 s s+a1 a0 a0 s+a1];
    y=[0 0 0 0 0 0 length length 0 length length length length length 0 0 length];
    z=[0 thickness thickness 0 0 thickness thickness 0 0 0 0 thickness thickness thickness thickness 0 0];
    plot3(x,y,z,'k');
    axis([0,width,0,length,0,thickness]);
    axis equal;
    
    %--------------this part is to color the image
    fill3([0 0 a0 a0],[0 0 0 0],[0 thickness thickness 0],c);
    hold on;
    fill3([0 s s 0],[0 length length 0],[thickness thickness 0 0],c);
    hold on;
    fill3([0 s s 0],[0 length length 0],[thickness thickness 0 0],c);
    hold on;
    fill3([s s+a1 s+a1 s],[length length length length],[0 0 thickness thickness],c);
    hold on;
    fill3([s+a1 a0 a0 s+a1],[length 0 0 length],[thickness thickness 0 0],c);
    hold on;
    fill3([0 s s+a1 a0],[0 length length 0],[thickness thickness thickness thickness],c);
    hold on;
    fill3([0 s s+a1 a0],[0 length length 0],[0 0 0 0],c);
    hold on;
    
    axis([0,width,0,length,0,thickness]);
    axis equal;
    
    disp('the up side is ');
    disp(a1);
    
    %thetamino=thetamin*180/pi;
    %-----------------------------------------------------------------
    
else(V_56 > Fslip && a1 > 56) %(5/8)'' bolt is enough and meet the code
    %just cut as it calculates
    Volmin=(2*b0-2*length*tan(thetamin)) * length / 2 * thickness;
    
    disp('theta(original):');
    disp(thetamin/pi*180);
    disp('Volume(original):')
    disp(Volmin);
    
    Br=3 * 0.8 * 1 * thickness * 15.875 * 825 /1000;
    if(Br<Fslip)
        disp('strip will fail(Br)');
    end
    subplot(3,2,6)
    a0=width;
    a1=b0-2*length*tan(thetamin);
    s=(a0-a1)/2;
    x=[0 0 a0 a0 0 0 s s 0 s s+a1 s+a1 s s+a1 a0 a0 s+a1];
    y=[0 0 0 0 0 0 length length 0 length length length length length 0 0 length];
    z=[0 thickness thickness 0 0 thickness thickness 0 0 0 0 thickness thickness thickness thickness 0 0];
    plot3(x,y,z,'k');
    axis([0,width,0,length,0,thickness]);
    axis equal;
    
    %--------------this part is to color the image
    fill3([0 0 a0 a0],[0 0 0 0],[0 thickness thickness 0],c);
    hold on;
    fill3([0 s s 0],[0 length length 0],[thickness thickness 0 0],c);
    hold on;
    fill3([0 s s 0],[0 length length 0],[thickness thickness 0 0],c);
    hold on;
    fill3([s s+a1 s+a1 s],[length length length length],[0 0 thickness thickness],c);
    hold on;
    fill3([s+a1 a0 a0 s+a1],[length 0 0 length],[thickness thickness 0 0],c);
    hold on;
    fill3([0 s s+a1 a0],[0 length length 0],[thickness thickness thickness thickness],c);
    hold on;
    fill3([0 s s+a1 a0],[0 length length 0],[0 0 0 0],c);
    hold on;
    
    axis([0,width,0,length,0,thickness]);
    axis equal;
    
    disp('the up side is ');
    disp(a1);
end
%-------------------------------------------------
%Last Brcheck
