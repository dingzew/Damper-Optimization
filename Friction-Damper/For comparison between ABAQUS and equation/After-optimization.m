%Notice

%1.Draw abaqus model result and numeric result.

%2.Before you use this program, make sure all the data and command has been
%cleared

%3.Open xls file or whatever you like to read the data

%4.Check is force named 'force' and disp named 'disp1'?if not copy the
%vector to force and disp1

%5.Change the number in two titles

%6.Never forget to change the Fmiu(friction coefficient) and Fnorm(clamping 
%force)
 %y1(i)=(2/3*Fmiu*Fnorm*((0.5*100)^3-(0.5*22.225)^3)/((0.5*100)^2-(0.5*22.225)^2)*2)*0.000001/0.25/cos(x1(i)/2/180*pi);
%--------------------------------------------------------------------------

Fmiu=0.4;
Fnorm=99428;

figure(1);
force=force/1000;
disp1=-disp1;
plot(disp1,force,'k'); %plot the abaqus 
grid on;
hold on;

disp1=-disp1;
x1=linspace(45,90,length(disp1));
for i=1:length(disp1)
    y1(i)=(2/3*Fmiu*Fnorm*((0.5*100)^3)/((0.5*100)^2)*2)*0.000001/0.25/cos(x1(i)/2/180*pi);
    %claculate the force of numeric model
end

for i=1:length(disp1)
    x1(i)=25*sqrt(2)-50*sin(x1(i)/2/180*pi); %change angle to displacement
 end
 
plot(x1,y1,'--k','linewidth',2); hold on; %plot numeric model


y1=-y1;  plot(x1,y1,'--k','linewidth',2);
y1=-y1;

 xlabel('Displacement(cm)');
 ylabel('Force(KN)');
axis([0,16,-50,50]); %adjust it after first figure comes out
%------------------------------------------
title('T3-15KN-Displacement-Force-Figure');
%do not forget to change it everytime
%------------------------------------------
%hleg1 = legend('ABAQUS Model','Analytical Equation','Orientation','horizontal');
%set(hleg1, 'Position', [.13,.94,.4,.05]);
 legend('ABAQUS Model','Analytical Equation','Location','NorthWest');
 %This is the end of the first figure
 
figure(2)

angle=asin(((disp1/2)+12.5*sqrt(2))/25);
angle_edit=angle/pi*180*2;
plot(angle_edit,force,'k');
hold on; grid on;
x1=linspace(45,90,length(disp1)); %reload angle
for i=1:length(disp1)
    y1(i)=(2/3*Fmiu*Fnorm*((0.5*100)^3/((0.5*100)^2)*2))*0.000001/0.25/cos(x1(i)/2/180*pi);
    %re-claculate the force of numeric model
end
plot(x1,y1,'--k','linewidth',2); hold on;
y1=-y1; plot(x1,y1,'--k','linewidth',2);  %draw negative image
y1=-y1; %make the data to it's orginal status
axis([45,90,-50,50]);
hold on; 
%y1=-y1; plot(angle_edit,y1,'r');
%y1=-y1;
xlabel('Angle(Degree)');
ylabel('Force(KN)');
%do not forget to change it everytime
%------------------------------------------
title('T3-15KN-Angle-Force-Figure');
%------------------------------------------
%hleg1 = legend('ABAQUS Model','Analytical Equation','Orientation','horizontal');
%set(hleg1, 'Position', [.13,.94,.4,.05]);

legend('ABAQUS Model','Analytical Equation','Location','NorthWest');

force=force*1000; %make the data to it's orginal status
