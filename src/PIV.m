clear
clc
close
load ('Re_418.mat');
nu= 1.57*10^-5;
nu0= 1.57*10^-5;
f=23.3; %% The oscillating frequency in Hz
de=0.045; %% the diameter of the duct
dv=sqrt(2*nu/(2*pi*f));
dv0=sqrt(2*nu0/(2*pi*f));
% %Reyt's work
for Ph= 1:1:25
    tc(Ph)=((Ph-1)/(25*f))-0.25*(1/f)+(Phase/360)*(1/f);
 %To calculate the velocity at different depths.
for D = 0:0.0001:1 %% D: dimensionless depth. This "For-loop" is to calculate the velocity at different depths.
   y = round(D*100+1);
   PsI(y)=D*de/2;
   wI(y)=1*exp(i*2*pi*f*tc(Ph))*(1-(besselj(0,PsI(y)*sqrt(-1*i*2*pi*f/nu))/besselj(0,(de/2)*sqrt(-1*i*2*pi*f/nu))));
end
    plot(U_Amp(Ph,:),y_dv(:),'*r',wI,-1*(PsI-(de/2))/dv0,'-k')
    xlim([-1.2 1.2]);
    ylim([0 7]);
    name=num2str(Ph);
   saveas(gcf,['Phase#',name,'.jpg']);
   saveas(gcf,['Phase#',name,'.fig']);
end