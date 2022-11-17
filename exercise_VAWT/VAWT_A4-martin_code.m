
clear all;
%close all;


airfoil=importdata('airfoil.txt');

aoatab=airfoil.data(:,1);
cltab=airfoil.data(:,2);
cdtab=airfoil.data(:,3);
% clinv=strucdata(:,4);
% clfs=strucdata(:,5);
% fstat=strucdata(:,6);

%plot(aoatab,cdtab)
[m,n]=size(aoatab);
for tab=1:m
    cntab(tab)=cltab(tab)*cosd(aoatab(tab))+cdtab(tab)*sind(aoatab(tab));
    cttab(tab)=cltab(tab)*sind(aoatab(tab))-cdtab(tab)*cosd(aoatab(tab));
end


dynflag=0       % 1 means a dynamic stall model used

rho=1.225;
%for k=1:30;
vo=8;
rad=3.0;
omega=14.0
%omega=1+(k-1)

tipspeedratio=omega*rad/vo

nb=3;
solidity=0.2
chord=solidity*rad/nb
pitch=0;

time_rot=2*pi/omega;
nrot=360;
dt=time_rot/nrot;
dt=0.0012;
ntime=nrot*20;
ntime=10000;

tau_airf=4*chord/(omega*rad);
tau=2*rad/vo;

            
a(1)=0.0;
theta(1,1)=0.0;
% for b=1:nb
%   fold(b)=0;
% end

for i=2:ntime;                 % time loop starts
  time(i)=(i-1)*dt;            % time
  theta(i,1)=theta(i-1,1)+omega*dt;    % azimuthal position blade #1

  % calculate the induced velocity in the mid plane from the
  % current value of axial induction factor
  wxcenter=a(i-1)*vo;
  
    for b=1:nb;        % Blade loop starts
   
    theta(i,b)=theta(i,1)+2*pi*(b-1)/nb; %azimuth position blade #b
    azimuth=mod(theta(i,b),2*pi);
    
    % Redistribute induced wind speed in mid plane to determine local
    % axial and tranverse induced wind at blade #b
    kk=0.4;
    wx(i,b)=wxcenter*(1.0-kk*sin(azimuth));
    wy(i,b)=kk*cos(theta(i,b))*wx(i,b);
                   
    % Calculate angle of attack and loads from blade #b    
    xpos(i,b)=-rad*sin(theta(i,b));
    ypos(i,b)=rad*cos(theta(i,b));
    vrelx(i,b)=omega*ypos(i,b)+vo-wx(i,b);
    vrely(i,b)=-omega*xpos(i,b)+wy(i,b);
    vtan(i,b)=omega*rad+(vo-wx(i,b))*cos(theta(i,b))+wy(i,b)*sin(theta(i,b));
    vnorm(i,b)=(vo-wx(i,b))*sin(theta(i,b))-wy(i,b)*cos(theta(i,b));
    vrel(i,b)=sqrt(vrelx(i,b)^2+vrely(i,b)^2);
    %vrel2(i,b)=sqrt(vnorm(i,b)^2+vtan(i,b)^2);
   
    aoa(i,b)=atan(vnorm(i,b)/vtan(i,b));
    aoadeg=(aoa(i,b)*180/pi)+pitch;
    
    % using static airfoil data only
    clstat(i,b)=interp1(aoatab,cltab,aoadeg);
%     clinvis(i,b)=interp1(aoatab,clinv,aoadeg);
%     clfullysep(i,b)=interp1(aoatab,clfs,aoadeg);
%     fstatic=interp1(aoatab,fstat,aoadeg);
%     f(b)=fstatic+(fold(b)-fstatic)*exp(-dt/tau_airf);
%     fold(b)=f(b);
%     
%     if(dynflag==1) 
%      cl(i,b)=f(b)*clinvis(i,b)+(1-f(b))*clfullysep(i,b);
%     else
%      cl(i,b)=clstat(i,b);
%     end
    cl(i,b)=clstat(i,b);
    cd(i,b)=interp1(aoatab,cdtab,aoadeg);
    
        
    lift(i,b)=0.5*rho*vrel(i,b)^2*chord*cl(i,b);
    drag(i,b)=0.5*rho*vrel(i,b)^2*chord*cd(i,b);
  
    cosbeta=vrely(i,b)/vrel(i,b);
    sinbeta=vrelx(i,b)/vrel(i,b);
    px(i,b)=lift(i,b)*cosbeta+drag(i,b)*sinbeta;    % load in flow direction
    py(i,b)=-lift(i,b)*sinbeta+drag(i,b)*cosbeta;  % load normal to flow direction

    ftan(i,b)=lift(i,b)*sin(aoa(i,b))-drag(i,b)*cos(aoa(i,b)); %tangential load
    fnorm(i,b)=lift(i,b)*cos(aoa(i,b))+drag(i,b)*sin(aoa(i,b));%normal to radius load
    cn(i,b)=cl(i,b)*cos(aoa(i,b))+cd(i,b)*sin(aoa(i,b));
    ct(i,b)=cl(i,b)*sin(aoa(i,b))-cd(i,b)*cos(aoa(i,b));
    
    cnplus(i,b)=cn(i,b)*(vrel(i,b)/vo)^2;
    
  
end;    % Blade loop ends

    % summing up the loads and torques from the all blades
    thrust=0;
    torque=0;
    for blade=1:nb
     thrust=thrust+px(i,blade);
     torque=torque+ftan(i,blade)*rad;
    end
    
    % compute the overall thrust coefficient
    cthrust(i)=thrust/(rho*vo^(2)*rad); % compute the overall thrust coefficient
    
    % compute power and power coefficient
    power=torque*omega;
    cp(i)=power/(0.5*rho*vo^3*2*rad);
    
%     fnplus(i)=cn(i,1)*(vrel(i,1)/vo)^(2);
%     ftplus(i)=ct(i,1)*(vrel(i,1)/vo)^(2);


    % Update the global axial induction factor     
    fg=1;
    if(a(i-1)>0.333)
        fg=0.25*(5-3*a(i-1));
    end
    astar=cthrust(i)/(4*(1-fg*a(i-1)));
    a(i)=astar+(a(i-1)-astar)*exp(-dt/tau);
    
    
end;                          % end time loop

for i=1:ntime
mom=0;
pxblade=0;
pyblade=0;
for b=1:nb
    mom=mom+rad*ftan(i,b);
    pxblade=pxblade+px(i,b);
    pyblade=pyblade+py(i,b);
end
moment(i)=mom;
pxtot(i)=pxblade;
pytot(i)=pyblade;
pt_time(i)=ftan(i,1);
pn_time(i)=fnorm(i,1);
end
   
cpmean=0.d0;
ctmean=0;
amean=0;
numav=0;
momentmean=0;

for i=1:ntime
  if( (theta(i,1)>5*2*pi)&&(theta(i,1)<6*2*pi) )
      cptime(i,1)=cp(i);
      
      cptime(i,2)=theta(i,1);
      cpmean=cpmean+cp(i);
      
      cttime(i)=cthrust(i);
      ctmean=ctmean+cttime(i);
      
      momentmean=momentmean+moment(i);
      
      amean=amean+a(i);
      numav=numav+1;
  end
end

numav;
cpmean=cpmean/numav
ctmean=ctmean/numav
amean=amean/numav;
moment=momentmean/numav;

% cplamda(k)=cpmean;
% ctlamda(k)=ctmean;
% lambda(k)=tipspeedratio;

%end

% dum=[lambda;cplamda;ctlamda]'

p=plot(time(2:ntime),ftan(2:ntime,1))%,'-',time(2:ntime),fnorm(2:ntime,1),'-')
p(1).LineWidth=2;
%p(2).LineWidth=2;
%legend('pt','pn')
xlabel('time [s]')
ylabel('tangential loads for blade#1 [N/m]')
title('B=3, R=3m, Vo=8m/s, \omega=14rad/s, S=0.2')
set(gca,'FontSize',16)
% ylim([-400 600])
% plot(time,pxtot,time,pytot)
% plot(time,cthrust)


% p=plot(time(2:ntime),cp(2:ntime),time(2:ntime),cthrust(2:ntime))
% p(1).LineWidth=3;
% p(2).LineWidth=3;
% legend('C_p','C_T')
% xlabel('time [s]')
% %ylabel('pn loads [N/m]')
% %title('Lateral loads, B=1')
% set(gca,'FontSize',16)
% xlim([0 4])

% p=plot(time,py(:,1),time,py(:,2),time,py(:,3))
% p(1).LineWidth=3;
% p(2).LineWidth=3;
% p(3).LineWidth=3;
% legend('pyB1','pyB2','pyB3')
% xlabel('time [s]')
% ylabel('py loads [N/m]')
% title('Lateral loads, B=3')
% set(gca,'FontSize',16)
% xlim([6 7])

% p=plot(time,py(:,1)+py(:,2)+py(:,3))
% % p(1).LineWidth=3;
% % p(2).LineWidth=3;
% % p(3).LineWidth=3;
% % legend('pyB1','pyB2','pyB3')
% xlabel('time [s]')
% ylabel('py total load [N/m]')
% title('Lateral loads, B=3')
% set(gca,'FontSize',16)
% xlim([6 7])

% fid=fopen('pnpt.txt','w');
% fprintf(fid,'time  pxtot   pytot   pxB1  pyB1\r\n')
% for i=1:ntime
%      fprintf(fid,'%12.5f %12.5f %12.5f %12.5f %12.5f\r\n',...
%         theta(i,1),pxtot(i),pytot(i),px(i,1),py(i,1));
% 
% if( (theta(i,1)>5*2*pi)&&(theta(i,1)<6*2*pi) )
%       fprintf(fid,'%12.5f %12.5f\r\n',...
%           theta(i,1)-5*2*pi,ftan(i,1)/(0.5*rho*chord*vo^2));
% end
%     fprintf(fid,'%12.5f %12.5f %12.5f\n',time(i),cp(i),cthrust(i));
%     fprintf(fid,'%12.5f %12.5f %12.5f\n',time(i),pn_time(i),pt_time(i));
% end
% fclose(fid);

