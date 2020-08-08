clear
clc
close
load ('205.mat');
de=0.045;
nu= 1.57*10^-5;
fs= 582.5;
f=23.3;
dv
%Reyt's work
for Ph= 1:1:nph
tc=t(Ph)+(phase/360)*(1/f)-0.25*(1/f);
for D= 0:0.0001:1
y= round(D*100+1);
wI(y)=Amp*exp(i*2*pi*f*tc)*(1-(besselj(0,PsI(y)*sqrt(-1*i*2*pi*f/nu))/besselj(0,(de/2)*sqrt(-1*i*2*pi*f/nu))));
 end


Sh=-0.0021; %%% This value is used to shift the data to get the wall position by trial
   subplot(1,2,2,'position',[0.45 0.1 0.45 0.8]);
   plot(uEA(Ph,:),(YEA(Ph,:)/(1000)+Sh)/dv,'*r',wI,(-1*PsI+(de/2))/dv,'-k',wF,(-1*PsF+0.5)*de/dv,'-g')
   set(gca,'fontsize',14)
   legend('Measured','Reyt,2013','Fan,1965','Location','northeast'); 
   xlim([-1.2*Amp 1.2*Amp]); 
   ylim([0.0 10]);
   name=num2str(Phs);
   xlabel('Acoustic Velocity (m/s)');
   ylabel('width/\delta_v');
   title(['Phase=',name,'deg']);
   saveas(gcf,['Phase#',name,'deg.jpg']);
   saveas(gcf,['Phase#',name,'deg.fig']);
   img = imread(['Phase#',name,'deg.jpg']);
   
 end
 close 
 clear
