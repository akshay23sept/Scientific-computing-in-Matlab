%%%% This code plots the ensemble average velocity distribution at
%%%% different phases over one acoustic cycle
clear 
clc
close
load('Spatial Average dis_100cycles'); %%% This file is created by the code in Apendix C3
de=0.045; %% diameter of the duct (m);
nu=1.57*10^-5; %%kinematic viscosity of the working fluid (m^2/s);
fs=582.5; %% The sampling frequency Hz
f=23.3; %% The oscillating frequency Hziyyes qnd zill plot
dv=sqrt(2*nu/(2*pi*f)); %%% Viscous penetration depth
nPh=round(fs/f);%%% the number of phases in each cycle
% n=input('The number of files to be read=');
NOC=round(n/nPh); %%% The total number of cycles
video=VideoWriter('Ensemble_Average.avi');
video.FrameRate=1;
open(video);
for P=1:1:nPh;
    YEA(P,:)=Ya(P,:);
    uEA(P,:)=ua(P,:);
    vEA(P,:)=va(P,:);
    for g=nPh+P:nPh:n
        YEA(P,:)=YEA(P,:)+Ya(g,:);
        uEA(P,:)=uEA(P,:)+ua(g,:);
        vEA(P,:)=vEA(P,:)+va(g,:);
    end
     YEA(P,:)= YEA(P,:)/NOC;
     uEA(P,:)= uEA(P,:)/NOC;
     vEA(P,:)= vEA(P,:)/NOC;
end
save('Ensemble Average dis','YEA','uEA','vEA');
 % Plotting the measured data versus the theoretical values
 for Ph=1:1:nPh
   %%%% Calculating the theoretical distibution in a circular pipe (Reyt's%%%% work of  2013) Theoritical data
   tc=t(Ph)+(Phase/360)*(1/f)-0.25*(1/f);  %%% The last term is added because the phase is measured from the time at which the velocity is zero (sine wave),
   %%% Whereas the starting time of the theoretical equation is at
   %%% maximum velocity (Cos wave). 
   for D = 0:0.0001:1 %% D: dimensionless depth. This "For-loop" is to calculate the velocity at different depths.
   y = round(D*100+1);
   PsI(y)=D*de/2;
   wI(y)=Amp*exp(i*2*pi*f*tc)*(1-(besselj(0,PsI(y)*sqrt(-1*i*2*pi*f/nu))/besselj(0,(de/2)*sqrt(-1*i*2*pi*f/nu))));
   end
    %%%% Calculating the theoretical distibution in a square duct (Fan's work 1965)
   
for D = 0:0.01:1 %% D: dimensionless depth. This "For-loop" is to calculate the velocity at different depths.
    y = round(D*100+1);
    PsF(y)=D/2;
    k=1;
    wt=t(Ph)*f*2*pi+(Phase/360)*2*pi;
    W=0; %% dimensionless width;
    P=1; %% aspect ratio of the rectangular duct (width/depth);
    FP=de*de*(2*pi*f)/(4*nu); %% dimensionless frequency parameter = (width*depth*omega)/(4*Dynmaic viscosity);
 for m=0:1:100; %% these "For-loops" are to calculate the double summation in equation (17)
    for  n=0:1:100;
  t1=(-1)^(m+n)/((2*m+1)*(2*n+1)); %% first term of equ (17)
  t2=cos((2*m+1)*(pi*W/2)); %% second term of equ (17)
  t3=cos((2*n+1)*(pi*D/2)); %% third term of equ (17)
  t4n=(FP*pi^2*[((2*m+1)^2/P)+((2*n+1)^2*P)]*cos(wt)/4)+(FP^2*sin(wt)); %% numerator of the fourth term of equ (17)
  t4d=(pi^4*([((2*m+1)^2/P)+((2*n+1)^2)*P]^2)/16)+FP^2; %% denominator of the fourth term of equ (17)
  vv(k)=t1*t2*t3*(t4n/t4d); 
  k=k+1;
    end
 end
 wF(y)=Amp*16*sum(vv)/pi^2; %% dimensionless velocity (%% equation # 17);
end 
  Phs=Phase+(Ph-1)*360/nPh;
if Phs>360
    Phs=round(Phs-360);
else
    Phs=round(Phs);
end  
subplot(1,2,1,'position',[0.1 0.35 0.25 0.3]);
th=0:1:360;
ampl=sin(th*pi/180);
plot(th,ampl,'-k',Phs,sin(Phs*pi/180),'*r','markersize',17);
set(gca,'fontsize',16)
xlabel('Phase, deg');
ylabel('U/Amp');
ylim([-1.2 1.2])
grid on
   %%%%% Plotting the theoretical values of both Fan and Reyt versus the measured values
   Sh=-0.0021; %%% This value is used to shift the data to get the wall position by trial
   subplot(1,2,2,'position',[0.45 0.1 0.45 0.8]);
   plot(uEA(Ph,:)/Amp,(YEA(Ph,:)/(1000)+Sh)/dv,'*r',wI/Amp,(-1*PsI+(de/2))/dv,'-k',wF/Amp,(-1*PsF+0.5)*de/dv,'-g')
   set(gca,'fontsize',16)
   legend('Measured','Reyt,2013','Fan,1965','Location','northeast'); 
   xlim([-1.2*Amp 1.2*Amp]); 
   ylim([0.0 10]);
   name=num2str(Phs);
   xlabel('Acoustic Velocity (m/s)');
   ylabel('width/\delta_v');
   title(['Phase=',name,'deg']);
   saveas(gcf,['Phase#',name,'deg.jpg']);00
   saveas(gcf,['Phase#',name,'deg.fig']);
   img = imread(['Phase#',name,'deg.jpg']);
   writeVideo(video,img); 
 end
 close 
 clear
