% PNP + LJ
%%orignal code is written by Alleb
%%revised by shixin xu
function argininesingletrap_xu
clear all;
global z1 z2 z3 N1 N2 N3 Dz1 Dz2 Dz3%% zi is the coordinate of node in 3 region; Ni is the node number; Dzi is the dirivative of basis
global DAz1 DAz2 DAz3 Az1 Az2 Az3 %% Azi area of cross section, DAi is 1./Azi Dzi used for matrix
global gij Vz %% gij is interaction Vz is the trap
global NaL NaR ClL ClR phiL phiR %% boundary conditions
global DNa Dq DCl %% diffusion constants
global zNa zCl zq Lap Gamma1 Gamma2 Gamma3 %% z* is the valence; Lap is the laplace matrix; Gammai is the dimensionless coefficient
global Dq1 Dq2 Dq3
global zq1 zq2 zq3
global iplot
global ton toff Vwidthval %% potent actives in ton and off in toff with smooth tanh function
format short e
iclose=input('Input 1 to close old figures, default is 1: ');
if isempty(iclose)
    iclose=1;
end
if iclose==1
    close all;
end
icon=input('Input 1 for new computing, 2 for consecutive computing, else for seeing result, default is 1: ');
if isempty(icon)
    icon=1;
end
if icon==1 | icon==2
    if icon==2
        dir argininesingletrapdefault*.mat
        dataini=input('Input the file to load: ','s');
        dataini=dataini(dataini~=' ');
        tit2=['load ' dataini];
        eval(tit2);
        u0=uall(end,:); u0=u0(:);
        u=u0;
        %Na1=u(1:N1); Cl1=u(N1+1:2*N1); q2=u(2*N1+1:2*N1+N2); Na3=u(2*N1+N2+1:2*N1+N2+N3); Cl3=u(2*N1+N2+N3+1:end);
        
        %%%%%%%%%%%%%%%%%%%%%%%% %% Revised by XU
        Na1=u(1:N1); Cl1=u(N1+1:2*N1);
        Aq1=u(2*N1+1:2*N1+N2); Aq2=u(2*N1+N2+1:2*N1+2*N2); Aq3=u(2*N1+2*N2+1:2*N1+3*N2);
        Na3=u(2*N1+3*N2+1:2*N1+3*N2+N3); Cl3=u(2*N1+3*N2+N3+1:end);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %         rhs1=-(zNa*Na1+zCl*Cl1)/Gamma1; rhs1(1)=phiL; rhs1(end)=0;
        %         rhs2=-(zq*q2)/Gamma2; rhs2([1 end])=0;
        %         rhs3=-(zNa*Na3+zCl*Cl3)/Gamma3; rhs3(1)=0; rhs3(end)=phiR;
        %         rhs=[rhs1;rhs2;rhs3]; warning on; phi=Lap\rhs; warning off; phi1=phi(1:N1); phi2=phi(N1+1:N1+N2); phi3=phi(N1+N2+1:end);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %% Revised by XU
        rhs1=-(zNa*Na1+zCl*Cl1)/Gamma1; rhs1(1)=phiL; rhs1(end)=0;
        
        rhs2=-(zq1*Aq1+zq2*Aq2+zq3*Aq3)/Gamma2; rhs2([1 end])=0;
        rhs3=-(zNa*Na3+zCl*Cl3)/Gamma3; rhs3(1)=0; rhs3(end)=phiR;
        rhs=[rhs1;rhs2;rhs3]; warning on; phi=Lap\rhs; warning off;
        phi1=phi(1:N1); phi2=phi(N1+1:N1+N2); phi3=phi(N1+N2+1:end);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(1); subplot(1,3,1);
        %         plot(z1,Na1,'b','LineWidth',2); hold on; plot(z3,Na3,'b','LineWidth',2);
        %         plot(z1,Cl1,'g','LineWidth',2); plot(z3,Cl3,'g','LineWidth',2);
        %         plot(z2,q2,'r','LineWidth',2);   hold off; grid on;
        %         xlabel('z','FontSize',20); ylabel('Na, Cl, q','FontSize',20); title('t=0','FontSize',20);
        
        plot(z1,Na1,'b','LineWidth',2); hold on; plot(z3,Na3,'b','LineWidth',2);
        plot(z1,Cl1,'g','LineWidth',2); plot(z3,Cl3,'g','LineWidth',2);
        plot(z2,Aq1,'r','LineWidth',2);   plot(z2,Aq2,'y','LineWidth',2);plot(z2,Aq3,'k','LineWidth',2);hold off; grid on;
        xlabel('z','FontSize',20); ylabel('Na, Cl, q','FontSize',20); title('t=0','FontSize',20);
        
        subplot(1,3,2);
        plot(z1,phi1,'k','LineWidth',2); hold on; plot(z2,phi2,'k','LineWidth',2); plot(z3,phi3,'k','LineWidth',2); hold off;
        grid on; xlabel('z','FontSize',20); ylabel('phi','FontSize',20);
        subplot(1,3,3); plot(z2,V,'Color','b','LineWidth',2); xlabel('z','FontSize',20); ylabel('V','FontSize',20); grid on;
        set(gcf,'Position',[37 198 1335 583]); drawnow;
        Gammas=[Gamma1 Gamma2 Gamma3]
        phiR=input('Input phiR, default is -60: ');
        if isempty(phiR)
            phiR=-60;
        end
        % ton=input('Input ton, default is 10: ');
        % if isempty(ton)
        ton=10;
        % end
        % toff=input('Input toff, default is 30: ');
        % if isempty(toff)
        toff=30;
        % end
        % Vwidthval=input('Input Vwidthval, default is 5: ');
        % if isempty(Vwidthval)
        Vwidthval=5;
        % end
        % Arginine trap:
        % double Gaussian
        % Vposr
        % Vposr=input('Input Vposr for right trap potential, default is 0.7: ');
        % if isempty(Vposr)
        Vposr=0.7;
        % end
        zrell=z2-LR; zrelr=z2-(LR+Vposr*L);
        % ampVl
        % ampVl=input('Input ampVl for left trap potential, default is 15: ');
        % if isempty(ampVl)
        ampVl=15;
        % end
        % ampVr
        % ampVr=input('Input ampVr for right trap potential, default is 0: ');
        % if isempty(ampVr)
        ampVr=0;
        % end
        % Vwidthl
        % Vwidthl=input('Input V width at left for trap potential, default is 2.5: ');
        % if isempty(Vwidthl)
        Vwidthl=2.5;
        % end
        % Vwidthr
        % Vwidthr=input('Input V width at right for trap potential, default is 2.5: ');
        % if isempty(Vwidthr)
        Vwidthr=2.5;
        % end
        % [Vposr ampVl ampVr Vwidthl Vwidthr]
        V=-ampVl*exp(-(Vwidthl*zrell).^2)-ampVr*exp(-(Vwidthr*zrelr).^2);
        Vz=2*ampVl*Vwidthl^2*zrell.*exp(-(Vwidthl*zrell).^2)+2*ampVr*Vwidthr^2*zrelr.*exp(-(Vwidthr*zrelr).^2);
        tspan=input('Input tspan, default is 0:0.01:50: ' );
        if isempty(tspan)
            tspan=(0:0.01:50)';
        end
    else
        % valence:
        % zNa=input('Input zNa, default is 1: ');
        % if isempty(zNa)
        zNa=1;
        % end
        % zCl=input('Input zCl, default is -1: ');
        % if isempty(zCl)
        zCl=-1;
        % end
        % grid:
        % Nz1=input('Input Nz1, default is 70: ');
        % if isempty(Nz1)
        Nz1=70;
        % end
        N1=Nz1+1; Nz3=Nz1; N3=N1;
        % Nz2=input('Input Nz2, default is 140: ');
        % if isempty(Nz2)
        Nz2=140;
        % end
        N2=Nz2+1;
        % [N1 N2 N3]
        % geometry:
        L=2;
        % LR=input('Input LR, default is 1: ');
        % if isempty(LR)
        LR=1;
        % end
        AR=input('Input AR, default is 1.5: ');
        if isempty(AR)
            AR=1.5;
        end
        Geometrydata=[L LR AR]
        % Diffusion coefficient:
        DNa=input('Input DNa, default is 1: ');
        if isempty(DNa)
            DNa=1;
        end
        DCl=input('Input DCl, default is 1: ');
        if isempty(DCl)
            DCl=1;
        end
        %         Dq=input('Input Dq, default is 1: ');
        %         if isempty(Dq)
        %             Dq=1;
        %         end
        %         Diffcoeff=[DNa DCl Dq]
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %% Revised by XU
        Dq1=input('Input Dq1, default is 1: ');
        if isempty(Dq1)
            Dq1=1;
        end
        Dq2=input('Input Dq2, default is 1: ');
        if isempty(Dq2)
            Dq2=1;
        end
        Dq3=input('Input Dq3, default is 1: ');
        if isempty(Dq3)
            Dq3=1;
        end
        Diffcoeff=[DNa DCl Dq1 Dq2 Dq3]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gamma:
        Gamma1=input('Input Gamma1, default is 1: ');
        if isempty(Gamma1)
            Gamma1=1;
        end
        Gamma3=Gamma1;
        % Gamma2=input('Input Gamma2, default is 1: ');
        % if isempty(Gamma2)
        Gamma2=1;
        % end
        Gammas=[Gamma1 Gamma2 Gamma3]
        % Boundary conditions:
        % NaR=input('Input NaR, default is 1: ');
        % if isempty(NaR)
        NaR=1;
        % end
        NaL=1; ClL=NaL; ClR=NaR;
        phiR=input('Input phiR, default is 0: ');
        if isempty(phiR)
            phiR=0;
        end
        phiL=0;
        BCdata=[NaL NaR ClL ClR phiL phiR]
        % LJ constants:
        % gij=input('Input gij, default is 0: ');
        % if isempty(gij)
        gij=0;
        % end
        % collocation matrix
        % Left reservoir length LR, radius AR
        [Dxi,xi]=cheb(Nz1); Dxi=fliplr(flipud(Dxi)); xi=flipud(xi); z1=(xi+1)*LR/2; Dz1=Dxi*2/LR; Dzz1=Dz1^2;
        
        eyez1=eye(N1); gz1=AR*ones(size(z1)); Az1=pi*gz1.^2; DAz1=diag(1./Az1)*Dz1;
        % port length L, radius 1
        [Dxi,xi]=cheb(Nz2); Dxi=fliplr(flipud(Dxi)); xi=flipud(xi); stp=LR; z2=(xi+1)*L/2+stp; Dz2=Dxi*2/L; Dzz2=Dz2^2;
        eyez2=eye(N2); gz2=ones(size(z2)); Az2=pi*gz2.^2; DAz2=diag(1./Az2)*Dz2;
        % right reservoir length LR, radius AR
        [Dxi,xi]=cheb(Nz3); Dxi=fliplr(flipud(Dxi)); xi=flipud(xi); stp=LR+L; z3=(xi+1)*LR/2+stp; Dz3=Dxi*2/LR; Dzz3=Dz3^2;
        eyez3=eye(N3); gz3=AR*ones(size(z3)); Az3=pi*gz3.^2; DAz3=diag(1./Az3)*Dz3;
        figure(1); plot(z1,gz1,'b','LineWidth',2); ylim([-8 8]); grid on; xlabel('z','FontSize',20); ylabel('r','FontSize',20); hold on;
        plot(z2,gz2,'b','LineWidth',2); plot(z3,gz3,'b','LineWidth',2);
        plot(z1,-gz1,'b','LineWidth',2); plot(z2,-gz2,'b','LineWidth',2); plot(z3,-gz3,'b','LineWidth',2); hold off;
        title('channel shape','FontSize',20); set(gcf,'Position',[626 289 746 492]); drawnow;
        % Laplace operator for electric potential:
        Lap1=Dzz1; Lap1(1,:)=eyez1(1,:); Lap1(end,:)=Gamma1*Az1(1)*Dz1(end,:);
        Lap2=Dzz2; Lap2(1,:)=eyez2(1,:); Lap2(end,:)=Gamma2*Az2(1)*Dz2(end,:);
        Lap3=Dzz3; Lap3(1,:)=eyez3(1,:); Lap3(end,:)=eyez3(end,:);
        Lap=blkdiag(Lap1,Lap2,Lap3);
        Lap(N1,N1+1:N1+N2)=-Gamma2*Az2(1)*Dz2(1,:); Lap(N1+1,1:N1)=-eyez1(end,:);
        Lap(N1+N2,N1+N2+1:end)=-Gamma3*Az3(1)*Dz3(1,:); Lap(N1+N2+1,N1+1:N1+N2)=-eyez2(end,:);
        figure(2); spy(Lap); title('Laplace','FontSize',20);
        % mass matrix: 3 species in pore
        mass1=ones(size(z1)); mass1([1 end])=0;
        mass2=ones(size(z2)); mass2([1 end])=0;
        mass3=ones(size(z3)); mass3([1 end])=0;
        % M=sparse(diag([mass1(:);mass1(:);mass2(:);mass3(:);mass3(:)]));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        M=sparse(diag([mass1(:);mass1(:);mass2(:);mass2(:);mass2(:);mass3(:);mass3(:)])); %% Revised by XU
        % initial condition:
        Na01=NaL*ones(size(z1)); Na03=NaR*ones(size(z3)); Cl01=Na01; Cl03=Na03;
        % Arginine initial distribution
        % Q=input('Input Q, default is 1: ');
        % if isempty(Q)
        Q=1;
        % end
        %         zq=input('Input zq, default is 1: ');
        %         if isempty(zq)
        %             zq=1;
        %         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Revised by XU
        zq1=input('Input zq1, default is 1: ');
        if isempty(zq1)
            zq1=1;
        end
        zq2=input('Input zq2, default is 1: ');
        if isempty(zq2)
            zq2=1;
        end
        zq3=input('Input zq3, default is 1: ');
        if isempty(zq3)
            zq3=1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        zrel=z2-(0.2+LR); sigma=0.03;
        
        % q02=Q*exp(-0.5*(zrel/sigma).^2)/(sigma*sqrt(2*pi));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %% Revised by XU
        Aq01=Q*exp(-0.5*(zrel/sigma).^2)/(sigma*sqrt(2*pi));
        Aq02=Aq01;
        Aq03=Aq01;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % warning off; tmp=Dz2\q02; warning on; Qcheck=tmp(end)-tmp(1);
        % [zq Q Qcheck]
        % Arginine trap:
        % double Gaussian
        % Vposr=input('Input Vposr for right trap potential, default is 0.7: ');
        % if isempty(Vposr)
        Vposr=0.7;
        % end
        zrell=z2-LR; zrelr=z2-(LR+Vposr*L);
        % ampVl=input('Input ampVl for left trap potential, default is 15: ');
        % if isempty(ampVl)
        ampVl=15;
        % end
        % ampVr=input('Input ampVr for right trap potential, default is 0: ');
        % if isempty(ampVr)
        ampVr=0;
        % end
        % Vwidthl=input('Input V width at left for trap potential, default is 2.5: ');
        % if isempty(Vwidthl)
        Vwidthl=2.5;
        % end
        % Vwidthr=input('Input V width at right for trap potential, default is 2.5: ');
        % if isempty(Vwidthr)
        Vwidthr=2.5;
        % end
        % [Vposr ampVl ampVr Vwidthl Vwidthr]
        V=-ampVl*exp(-(Vwidthl*zrell).^2)-ampVr*exp(-(Vwidthr*zrelr).^2);
        Vz=2*ampVl*Vwidthl^2*zrell.*exp(-(Vwidthl*zrell).^2)+2*ampVr*Vwidthr^2*zrelr.*exp(-(Vwidthr*zrelr).^2);
        % initial setup plotting
        figure(3); subplot(1,2,1); plot(z1,Na01,'Color','b','LineWidth',2); hold on; plot(z3,Na03,'Color','b','LineWidth',2);
        plot(z1,Cl01,'Color','g','LineWidth',2); plot(z3,Cl03,'Color','g','LineWidth',2);
        %plot(z2,q02,'Color','r','LineWidth',2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Revised by XU
        plot(z2,Aq01,'Color','r','LineWidth',2);plot(z2,Aq02,'Color','y','LineWidth',2);plot(z2,Aq03,'Color','k','LineWidth',2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        hold off; xlabel('z','FontSize',20); ylabel('initial Na, Cl, q','FontSize',20); grid on;
        subplot(1,2,2); plot(z2,V,'Color','b','LineWidth',2); xlabel('z','FontSize',20); ylabel('V','FontSize',20); grid on;
        set(gcf,'Position',[37 198 1129 583]); drawnow;
        
        %u0=[Na01(:);Cl01(:);q02(:);Na03(:);Cl03(:)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Revised by XU
        u0=[Na01(:);Cl01(:);Aq01(:);Aq02(:);Aq03(:);Na03(:);Cl03(:)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % ODE parameters:
        tspan=input('Input tspan, default is 0:0.01:10: ' );
        if isempty(tspan)
            tspan=(0:0.01:10)';
        end
    end
    % mxstep=input('Input max time step, default is 0.01: ' );
    % if isempty(mxstep)
    mxstep=0.01;
    % end
    iplot=input('Input 1 to plot during computation, default is 0: ' );
    if isempty(iplot)
        iplot=0;
    end
    options=odeset('RelTol',1e-3,'AbsTol',1e-5,'Mass',M,'MaxStep',mxstep);
    %tic;
    ton
    % pause;
    [t,uall]=ode15s(@pnp1d,tspan,u0,options);
    %et=toc;
    if ~isempty(ton)
        % tit=['argininesingletrapdefault-phiR=' num2str(phiR) '-Q=' num2str(Q) '-zq=' num2str(zq)  ...
        %     '-Dq=' num2str(Dq) '-ampVl=' num2str(ampVl) '-ampVr=' num2str(ampVr) '-Vposr=' num2str(Vposr) '-Vwidthl=' num2str(Vwidthl) '-Vwidthr=' num2str(Vwidthr) ...
        %     '-Nz1=' num2str(Nz1) '-Nz2=' num2str(Nz2) '-ton=' num2str(ton) '-toff=' num2str(toff) '-AR=' num2str(AR) ...
        %     '-t=' num2str(tspan(1))  '-' num2str(tspan(end))]
        tit=['argininesingletrapdefault-phiR=' num2str(phiR) '-DNa=' num2str(DNa) '-Dq1=' num2str(Dq1)  '-zq1=' num2str(zq1) '-Gamma1=' num2str(Gamma1) ...
            '-ton=' num2str(ton) '-toff=' num2str(toff) '-AR=' num2str(AR) '-t=' num2str(tspan(1))  '-' num2str(tspan(end))]
    else
        % tit=['argininesingletrapdefault-phiR=' num2str(phiR) '-Q=' num2str(Q) '-zq=' num2str(zq)  ...
        %     '-Dq=' num2str(Dq) '-ampVl=' num2str(ampVl) '-ampVr=' num2str(ampVr) '-Vposr=' num2str(Vposr) '-Vwidthl=' num2str(Vwidthl) '-Vwidthr=' num2str(Vwidthr) ...
        %     '-Nz1=' num2str(Nz1) '-Nz2=' num2str(Nz2) '-AR=' num2str(AR) '-t=' num2str(tspan(1))  '-' num2str(tspan(end))]
        tit=['argininesingletrapdefault-phiR=' num2str(phiR) '-DNa=' num2str(DNa) '-Dq1=' num2str(Dq1)  '-zq1=' num2str(zq1) '-Gamma1=' num2str(Gamma1) ...
            '-AR=' num2str(AR) '-t=' num2str(tspan(1)) '-' num2str(tspan(end))]
    end
    tit1=['save ' tit '.mat'];
    eval(tit1);
else
    close all;
    dir argininesingletrapdefault*.mat
    dataini=input('Input the file to load: ','s');
    dataini=dataini(dataini~=' ');
    tit2=['load ' dataini];
    eval(tit2);
end
%% mass conservation checking:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% u=uall(1,:); u=u.'; q2=u(2*N1+1:2*N1+N2); warning off; tmp=Dz2\q2; warning on; Q2initial=tmp(end)-tmp(1);
% u=uall(end,:); u=u.'; q2=u(2*N1+1:2*N1+N2); warning off; tmp=Dz2\q2; warning on; Q2final=tmp(end)-tmp(1);
% Q2initialfinalcheck=[Q2initial Q2final]


%%%%%%%%%%%%%%%%%%%%% Revised by XU%%%%%%%%%%%%%%%
u=uall(1,:); u=u.'; Aq1=u(2*N1+1:2*N1+N2); warning off; tmp=Dz2\Aq1; warning on; AQ1initial=tmp(end)-tmp(1);
u=uall(end,:); u=u.'; Aq1=u(2*N1+1:2*N1+N2); warning off; tmp=Dz2\Aq1; warning on; AQ1final=tmp(end)-tmp(1);
AQ1initialfinalcheck=[AQ1initial AQ1final]

u=uall(1,:); u=u.'; Aq2=u(2*N1+N2+1:2*N1+2*N2); warning off; tmp=Dz2\Aq2; warning on; AQ2initial=tmp(end)-tmp(1);
u=uall(end,:); u=u.'; Aq2=u(2*N1+N2+1:2*N1+2*N2); warning off; tmp=Dz2\Aq2; warning on; AQ2final=tmp(end)-tmp(1);
AQ2initialfinalcheck=[AQ2initial AQ2final]

u=uall(1,:); u=u.'; Aq3=u(2*N1+2*N2+1:2*N1+3*N2); warning off; tmp=Dz2\Aq3; warning on; AQ3initial=tmp(end)-tmp(1);
u=uall(end,:); u=u.'; Aq3=u(2*N1+2*N2+1:2*N1+3*N2); warning off; tmp=Dz2\Aq3; warning on; AQ3final=tmp(end)-tmp(1);
AQ3initialfinalcheck=[AQ1initial AQ1final]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Iqmax max current recorded in zone 2.
% AQi number of arginine left at trap.
% Q1, Q3 number of net charges in reservoirs (zone 1 and 3), will be used
% to calculate capacitance current IQ1 and IQ3.
% phiRv: added voltage
% I1=[]; I2=[]; I3=[]; Q1=[]; Q2=[]; Q3=[]; phiRv=[]; tv=[]; Na1all=[]; Cl1all=[]; q2all=[]; Na3all=[]; Cl3all=[];
% phi1all=[]; phi2all=[]; phi3all=[]; phiall=[];

I1=[]; I2=[]; I3=[]; %% current 
Q1=[]; AQ1=[]; AQ2=[]; AQ3=[];  Q3=[]; 
phiRv=[]; tv=[];
Na1all=[]; Cl1all=[];   Na3all=[]; Cl3all=[];
phi1all=[]; phi2all=[]; phi3all=[]; phiall=[]; 
Aq1all=[]; Aq2all=[]; Aq3all=[];  %% Aqiall the time series of Argin 

int=input('Input interval, default is 50: ' );
if isempty(int)
    int=50;
end
for i=1:int:length(tspan)
    tv=[tv;t(i)];
    u=uall(i,:); u=u.';
    
    %     Na1=u(1:N1); Cl1=u(N1+1:2*N1); q2=u(2*N1+1:2*N1+N2); Na3=u(2*N1+N2+1:2*N1+N2+N3); Cl3=u(2*N1+N2+N3+1:end);
    %     Na1all=[Na1all;Na1.']; Cl1all=[Cl1all;Cl1.']; q2all=[q2all;q2.']; Na3all=[Na3all;Na3.']; Cl3all=[Cl3all;Cl3.'];
    %     q2l=q2; q2l(z2>LR+0.5*L)=0; warning off; tmp=Dz2\q2l; warning on; Q2l=tmp(end)-tmp(1); Q2=[Q2;Q2l];
    
    %%%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%%%%%%%%
    Na1=u(1:N1); Cl1=u(N1+1:2*N1); 
    Aq1=u(2*N1+1:2*N1+N2); Aq2=u(2*N1+N2+1:2*N1+2*N2); Aq3=u(2*N1+2*N2+1:2*N1+3*N2);
    Na3=u(2*N1+3*N2+1:2*N1+3*N2+N3); Cl3=u(2*N1+3*N2+N3+1:end);
    Na1all=[Na1all;Na1.']; Cl1all=[Cl1all;Cl1.']; 
    Aq1all=[Aq1all;Aq1.']; Aq2all=[Aq2all;Aq2.']; Aq3all=[Aq3all;Aq3.'];
    Na3all=[Na3all;Na3.']; Cl3all=[Cl3all;Cl3.'];
    Aq1l=Aq1; Aq1l(z2>LR+0.5*L)=0; warning off; tmp=Dz2\Aq1l; warning on; AQ1l=tmp(end)-tmp(1); AQ1=[AQ1;AQ1l];
    Aq2l=Aq2; Aq2l(z2>LR+0.5*L)=0; warning off; tmp=Dz2\Aq2l; warning on; AQ2l=tmp(end)-tmp(1); AQ2=[AQ2;AQ2l];
    Aq3l=Aq3; Aq3l(z2>LR+0.5*L)=0; warning off; tmp=Dz2\Aq3l; warning on; AQ3l=tmp(end)-tmp(1); AQ3=[AQ3;AQ3l];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    rho1=zNa*Na1+zCl*Cl1; warning off; tmp=Dz1\(Az1.*rho1); warning on; Q1=[Q1;(tmp(end)-tmp(1))];
    rho3=zNa*Na3+zCl*Cl3; warning off; tmp=Dz3\(Az3.*rho3); warning on; Q3=[Q3;(tmp(end)-tmp(1))];
    Vwidth=Vwidthval;
    if isempty(ton)
        phiRactual=0;
    else
        phiRactual=(phiR*tanh(Vwidth*(t(i)-(ton+1)))-phiR*tanh(Vwidth*(t(i)-(toff-1))))/2;
    end
    phiRv=[phiRv;phiRactual];
    rhs1=-(zNa*Na1+zCl*Cl1)/Gamma1; rhs1(1)=phiL; rhs1(end)=0;
    
%     rhs2=-(zq*q2)/Gamma2; rhs2([1 end])=0;
    %%%%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%%%%%%%%
    rhs2=-(zq1*Aq1+zq2*Aq2+zq3*Aq3)/Gamma2; rhs2([1 end])=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rhs3=-(zNa*Na3+zCl*Cl3)/Gamma3; rhs3(1)=0; rhs3(end)=phiRactual;
    rhs=[rhs1;rhs2;rhs3]; warning off; phi=Lap\rhs; warning on; phi1=phi(1:N1); phi2=phi(N1+1:N1+N2); phi3=phi(N1+N2+1:end);
    phi1all=[phi1all;phi1.']; phi2all=[phi2all;phi2.']; phi3all=[phi3all;phi3.']; phiall=[phiall;phi.'];
    Naz1=Dz1*Na1; Naz3=Dz3*Na3;
    Clz1=Dz1*Cl1; Clz3=Dz3*Cl3;
    %qz2=Dz2*q2;
    
    %%%%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%%%%%%%%
    Aqz1=Dz2*Aq1; Aqz2=Dz2*Aq2; Aqz3=Dz2*Aq3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    phiz1=Dz1*phi1; phiz3=Dz3*phi3; phiz2=Dz2*phi2;
    JNaz1=-Az1.*DNa.*(Naz1+zNa*Na1.*phiz1);
    JClz1=-Az1.*DCl.*(Clz1+zCl*Cl1.*phiz1);
    JNaz3=-Az3.*DNa.*(Naz3+zNa*Na3.*phiz3);
    JClz3=-Az3.*DCl.*(Clz3+zCl*Cl3.*phiz3);
   % Jqz2=-Az2.*Dq.*(qz2+zq*q2.*phiz2+gij*q2.*qz2+q2.*Vz);
   
   %%%%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    JAqz1=-Az2.*Dq1.*(Aqz1+zq1*Aq1.*phiz2+gij*Aq1.*Aqz1+Aq1.*Vz);
    JAqz2=-Az2.*Dq2.*(Aqz2+zq2*Aq2.*phiz2+gij*Aq2.*Aqz2+Aq2.*Vz);
    JAqz3=-Az2.*Dq3.*(Aqz3+zq3*Aq3.*phiz2+gij*Aq3.*Aqz3+Aq3.*Vz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    J1=zNa*JNaz1+zCl*JClz1;
    J3=zNa*JNaz3+zCl*JClz3;
    %J2=zq*Jqz2;
    
    %%%%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J2=zq1*JAqz1+zq2*JAqz2+zq3*JAqz3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    I1=[I1;J1.']; I2=[I2;J2.']; I3=[I3;J3.'];
end
q1all=zNa*Na1all+zCl*Cl1all; q3all=zNa*Na3all+zCl*Cl3all;
% figure(4);
% subplot(3,3,1); mesh(z1,tv,I1); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('flux J1','FontSize',20);
% title(['Gamma1=' num2str(Gamma1) ', DNa=' num2str(DNa) ', zq=' num2str(zq) ', phiR=' num2str(phiR)],'FontSize',20);
% subplot(3,3,2); mesh(z2,tv,I2); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('flux J2','FontSize',20);
% subplot(3,3,3); mesh(z3,tv,I3); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('flux J3','FontSize',20);
% subplot(3,3,4); mesh(z2,tv,q2all); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('arginine','FontSize',20);
% subplot(3,3,5); mesh(z1,tv,Na1all); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('Na','FontSize',20);
% hold on; mesh(z3,tv,Na3all); hold off;
% subplot(3,3,6); mesh(z1,tv,Cl1all); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('Cl','FontSize',20);
% hold on; mesh(z3,tv,Cl3all); hold off;
% subplot(3,3,7); mesh(z1,tv,q1all); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('net charge q','FontSize',20);
% hold on; mesh(z3,tv,q3all); hold off;
% subplot(3,3,8); mesh(z2,tv,phi2all); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('phi','FontSize',20); hold on;
% mesh(z1,tv,phi1all); mesh(z3,tv,phi3all); hold off;
% subplot(3,3,9); plot(z1,gz1,'b','LineWidth',2); ylim([-8 8]); grid on; xlabel('z','FontSize',20); ylabel('r','FontSize',20); hold on;
% plot(z2,gz2,'b','LineWidth',2); plot(z3,gz3,'b','LineWidth',2);
% plot(z1,-gz1,'b','LineWidth',2); plot(z2,-gz2,'b','LineWidth',2); plot(z3,-gz3,'b','LineWidth',2); hold off;
% title('channel shape','FontSize',20); drawnow;

%%%%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
subplot(3,3,1); mesh(z1,tv,I1); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('flux J1','FontSize',20);
title(['Gamma1=' num2str(Gamma1) ', DNa=' num2str(DNa) ', zq=' num2str(zq) ', phiR=' num2str(phiR)],'FontSize',20);

subplot(3,3,2); mesh(z2,tv,I2); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('flux J2','FontSize',20);

subplot(3,3,3); mesh(z3,tv,I3); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('flux J3','FontSize',20);

subplot(3,3,4); mesh(z2,tv,Aq2all); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('arginine 2','FontSize',20);

subplot(3,3,5); mesh(z1,tv,Na1all); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('Na','FontSize',20);
hold on; mesh(z3,tv,Na3all); hold off;

subplot(3,3,6); mesh(z1,tv,Cl1all); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('Cl','FontSize',20);
hold on; mesh(z3,tv,Cl3all); hold off;

subplot(3,3,7); mesh(z1,tv,q1all); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('net charge q','FontSize',20);
hold on; mesh(z3,tv,q3all); hold off;

subplot(3,3,8); mesh(z2,tv,phi2all); xlabel('z','FontSize',20); ylabel('t','FontSize',20); zlabel('phi','FontSize',20); hold on;
mesh(z1,tv,phi1all); mesh(z3,tv,phi3all); hold off;

subplot(3,3,9); plot(z1,gz1,'b','LineWidth',2); ylim([-8 8]); grid on; xlabel('z','FontSize',20); ylabel('r','FontSize',20); hold on;
plot(z2,gz2,'b','LineWidth',2); plot(z3,gz3,'b','LineWidth',2);

plot(z1,-gz1,'b','LineWidth',2); plot(z2,-gz2,'b','LineWidth',2); plot(z3,-gz3,'b','LineWidth',2); hold off;
title('channel shape','FontSize',20); drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Iat0p5 Iat1p5 Iat2 Iat2p5 Iat3p5 currents measured at z=0.5 1.5 2 2.5 3.5
% z=0.5 at i1=36
% z=1.5065 at i2=48; z=2.0 at i2=71; z=2.4935 at i2=94
% z=3.5 at i3=36
Iat0p5=I1(:,36); Iat1p5=I2(:,48); Iat2=I2(:,71); Iat2p5=I2(:,94); Iat3p5=I3(:,36);
I2max=[];
for i=1:length(tv)
    [tmp1,tmp2]=max(abs(I2(i,:)));
    I2max=[I2max;I2(i,tmp2)];
end
IQ1=(Q1(3:end)-Q1(1:end-2))/0.02; IQ1=[IQ1(1);IQ1;IQ1(end)];
IQ3=(Q3(3:end)-Q3(1:end-2))/0.02; IQ3=[IQ3(1);IQ3;IQ3(end)];
figure(5);
subplot(3,2,1);
plot(tv,Iat1p5,'Color','b','LineWidth',2);
xlabel('t','FontSize',20); ylabel('I at z=1.5','FontSize',20); grid on;
title(['Gamma1=' num2str(Gamma1) ', DNa=' num2str(DNa) ', zq1=' num2str(zq1) ', zq2=' num2str(zq2) ',zq3=' num2str(zq3) ', phiR=' num2str(phiR)],'FontSize',20);
subplot(3,2,3);
plot(tv,Iat2,'Color','g','LineWidth',2);
xlabel('t','FontSize',20); ylabel('I at z=2','FontSize',20); grid on;
subplot(3,2,5);
plot(tv,Iat2p5,'Color','r','LineWidth',2);
xlabel('t','FontSize',20); ylabel('I at z=2.5','FontSize',20); grid on;
subplot(3,2,2);
plot(tv,I2max,'Color','k','LineWidth',2);
xlabel('t','FontSize',20); ylabel('max I_{arg}','FontSize',20); grid on;
subplot(3,2,4);
plot(tv,Iat0p5,'Color','b','LineWidth',2);
xlabel('t','FontSize',20); ylabel('I at z=0.5','FontSize',20); grid on;
subplot(3,2,6);
plot(tv,Iat3p5,'Color','g','LineWidth',2);
xlabel('t','FontSize',20); ylabel('I at z=3.5','FontSize',20); grid on;
figure(6);
% subplot(2,3,1); plot(tv,Q2,'Color','b','LineWidth',2); xlabel('t','FontSize',20); ylabel('arginine at trap','FontSize',20); grid on;

%%%%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,1); 
plot(tv,AQ1,'Color','b','LineWidth',2);
hold on; plot(tv,AQ2,'Color','y','LineWidth',2);
hold on; plot(tv,AQ3,'Color','k','LineWidth',2);
hold off
xlabel('t','FontSize',20); ylabel('arginine at trap','FontSize',20);
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

title(['Gamma1=' num2str(Gamma1) ', DNa=' num2str(DNa) ', zq=' num2str(zq) ', phiR=' num2str(phiR)],'FontSize',20);
subplot(2,3,2); plot(tv,Q1,'Color','g','LineWidth',2); xlabel('t','FontSize',20); ylabel('net charge Q1 of left capacitor','FontSize',20); grid on;
subplot(2,3,3); plot(tv,Q3,'Color','r','LineWidth',2); xlabel('t','FontSize',20); ylabel('net charge Q3 of right capacitor','FontSize',20); grid on;
subplot(2,3,4); plot(tv,-phiRv,'Color','b','LineWidth',2); xlabel('t','FontSize',20); ylabel('Voltage','FontSize',20); grid on;
subplot(2,3,5); plot(tv,IQ1,'Color','r','LineWidth',2);
xlabel('t','FontSize',20); grid on;
hold on; plot(tv,IQ3,'Color','k','LineWidth',2);
xlabel('t','FontSize',20); ylabel('I_{cap}=dQ1,3/dt','FontSize',20); grid on; hold off;
hh=legend('dQ1/dt','dQ3/dt'); set(hh,'FontSize',20);
% subplot(2,3,6); plot(z2,q2all(1,:).','r','LineWidth',2); hold on;

%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%
subplot(2,3,6); 
plot(z2,Aq1all(1,:).','r','LineWidth',2); hold on;
plot(z2,Aq2all(1,:).','y','LineWidth',2);
plot(z2,Aq3all(1,:).','k','LineWidth',2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(z1,Na1all(1,:).','b','LineWidth',2);
plot(z1,Cl1all(1,:).','g','LineWidth',2);
plot(z3,Na3all(1,:).','b','LineWidth',2);
plot(z3,Cl3all(1,:).','g','LineWidth',2);
hold off; grid on;
xlabel('z','FontSize',20); ylabel('Na, Cl, arg.','FontSize',20); title('t=0','FontSize',20);
hh=legend('Arg1','Arg2','Arg3','Na','Cl'); set(hh,'FontSize',20); drawnow;
iplot=input('Input 1 for animation, default is 1: ' );
if isempty(iplot)
    iplot=1;
end
if iplot==1
    int=input('Input interval, default is 1: ' );
    if isempty(int)
        int=1;
    end
    ipause=input('Input 1 to pause, default is 0: ' );
    if isempty(ipause)
        ipause=0;
    end
    %ymax=max(q2all(:)); ymax1=max([q2all(:);Na1all(:);Na3all(:);Cl1all(:);Cl3all(:)]);
    
    %%%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ymax=max([Aq1all(:);Aq2all(:);Aq3all(:)]);
    ymax1=max([Aq1all(:);Aq2all(:);Aq3all(:);Na1all(:);Na3all(:);Cl1all(:);Cl3all(:)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for i=1:int:length(tv)
%         figure(8); 
%         subplot(2,2,1); plot(z2,q2all(i,:).','r','LineWidth',2); grid on; xlabel('z','FontSize',20); ylabel('Arginine','FontSize',20); ylim([-1 ymax+5]);
%         title(['t=' num2str(tv(i)) ', ton=' num2str(ton) ', toff=' num2str(toff)],'FontSize',20);
%         subplot(2,2,2);
%         plot(z2,V+phi2all(i,:).','b','LineWidth',2); grid on; xlabel('z','FontSize',20); ylabel('V+phi','FontSize',20); title(['phiR on =' num2str(phiR)],'FontSize',20);
%         subplot(2,2,3); plot(z1,Na1all(i,:).','b','LineWidth',2); hold on;
%         plot(z1,Cl1all(i,:).','g','LineWidth',2);
%         plot(z2,q2all(i,:).','r','LineWidth',2);
%         plot(z3,Na3all(i,:).','b','LineWidth',2);
%         plot(z3,Cl3all(i,:).','g','LineWidth',2);
%         hold off; grid on; xlabel('z','FontSize',20); ylabel('Na, Cl, Arginine','FontSize',20); ylim([-1 ymax1]);
%         title(['Gamma1=' num2str(Gamma1) ', DNa=' num2str(DNa) ', zq=' num2str(zq) ],'FontSize',20);
%         hh=legend('Na','Cl','Arg'); set(hh,'FontSize',20);
%         subplot(2,2,4);
%         plot(z1,phi1all(i,:).','k','LineWidth',2); hold on;
%         plot(z2,phi2all(i,:).','k','LineWidth',2);
%         plot(z3,phi3all(i,:).','k','LineWidth',2); hold off;
%         grid on; xlabel('z','FontSize',20); ylabel('\phi','FontSize',20); drawnow;
%         if ipause==1
%             pause;
%         else
%             pause(0.1);
%         end
%     end


    for i=1:int:length(tv)
        figure(8); 
        subplot(2,3,1); 
        plot(z2,Aq1all(i,:).','r','LineWidth',2);hold on;
   
        grid on; xlabel('z','FontSize',20); ylabel('Arginine 1 ','FontSize',20); ylim([-1 ymax+5]);
        title(['t=' num2str(tv(i)) ', ton=' num2str(ton) ', toff=' num2str(toff)],'FontSize',20);
        
        subplot(2,3,2);
        plot(z2,Aq2all(i,:).','y','LineWidth',2);
        grid on; xlabel('z','FontSize',20); ylabel('Arginine 2','FontSize',20); ylim([-1 ymax+5]);
        title(['t=' num2str(tv(i)) ', ton=' num2str(ton) ', toff=' num2str(toff)],'FontSize',20);
        
        subplot(2,3,3);
        plot(z2,Aq3all(i,:).','k','LineWidth',2);
        grid on; xlabel('z','FontSize',20); ylabel('Arginine 3','FontSize',20); ylim([-1 ymax+5]);
        title(['t=' num2str(tv(i)) ', ton=' num2str(ton) ', toff=' num2str(toff)],'FontSize',20);
        
        
        subplot(2,3,4);
        plot(z2,V+phi2all(i,:).','b','LineWidth',2); grid on; xlabel('z','FontSize',20); ylabel('V+phi','FontSize',20); title(['phiR on =' num2str(phiR)],'FontSize',20);
        
        subplot(2,3,5); plot(z1,Na1all(i,:).','b','LineWidth',2); hold on;
        plot(z1,Cl1all(i,:).','g','LineWidth',2);
        plot(z2,Aq1all(i,:).','r','LineWidth',2); 
        plot(z2,Aq2all(i,:).','y','LineWidth',2);
        plot(z2,Aq3all(i,:).','k','LineWidth',2);
        plot(z3,Na3all(i,:).','b','LineWidth',2);
        plot(z3,Cl3all(i,:).','g','LineWidth',2);
        hold off; grid on; xlabel('z','FontSize',20); ylabel('Na, Cl, Arginine','FontSize',20); ylim([-1 ymax1]);
        title(['Gamma1=' num2str(Gamma1) ', DNa=' num2str(DNa) ', zq1=' num2str(zq1) ],'FontSize',20);
        hh=legend('Na','Cl','Arg'); set(hh,'FontSize',20);
        subplot(2,3,6);
        plot(z1,phi1all(i,:).','k','LineWidth',2); hold on;
        plot(z2,phi2all(i,:).','k','LineWidth',2);
        plot(z3,phi3all(i,:).','k','LineWidth',2); hold off;
        grid on; xlabel('z','FontSize',20); ylabel('\phi','FontSize',20); drawnow;
        if ipause==1
            pause;
        else
            pause(0.1);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dudt=pnp1d(t,u)
global z1 z2 z3 N1 N2 N3 Dz1 Dz2 Dz3
global DAz1 DAz2 DAz3 Az1 Az2 Az3
global gij Vz
global NaL NaR ClL ClR phiL phiR
global DNa Dq DCl
global zNa zCl zq Lap Gamma1 Gamma2 Gamma3
global iplot
global ton toff Vwidthval
global Dq1 Dq2 Dq3
global zq1 zq2 zq3
if isempty(ton)
    phiRactual=phiR;
else
    Vwidth=Vwidthval;
    phiRactual=(phiR*tanh(Vwidth*(t-(ton+1)))-phiR*tanh(Vwidth*(t-(toff-1))))/2;
end
% Na1=u(1:N1); Cl1=u(N1+1:2*N1); q2=u(2*N1+1:2*N1+N2); Na3=u(2*N1+N2+1:2*N1+N2+N3); Cl3=u(2*N1+N2+N3+1:end);
%Naz1=Dz1*Na1; Naz3=Dz3*Na3; Clz1=Dz1*Cl1; Clz3=Dz3*Cl3; qz2=Dz2*q2;

%%%%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%%%%%%%%
Na1=u(1:N1); Cl1=u(N1+1:2*N1); 
Aq1=u(2*N1+1:2*N1+N2); Aq2=u(2*N1+N2+1:2*N1+2*N2); Aq3=u(2*N1+2*N2+1:2*N1+3*N2);
Na3=u(2*N1+3*N2+1:2*N1+3*N2+N3); Cl3=u(2*N1+3*N2+N3+1:end);
Naz1=Dz1*Na1; Naz3=Dz3*Na3;
Clz1=Dz1*Cl1; Clz3=Dz3*Cl3;
Aqz1=Dz2*Aq1; Aqz2=Dz2*Aq2; Aqz3=Dz2*Aq3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% electric field:
rhs1=-(zNa*Na1+zCl*Cl1)/Gamma1; rhs1(1)=phiL; rhs1(end)=0;
%rhs2=-(zq*q2)/Gamma2; rhs2([1 end])=0;

%%%%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%%%%%%%%
rhs2=-(zq1*Aq1+zq2*Aq2+zq3*Aq3)/Gamma2; rhs2([1 end])=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhs3=-(zNa*Na3+zCl*Cl3)/Gamma3; rhs3(1)=0; rhs3(end)=phiRactual;
rhs=[rhs1;rhs2;rhs3]; warning on; phi=Lap\rhs; warning off;
phi1=phi(1:N1); phi2=phi(N1+1:N1+N2); phi3=phi(N1+N2+1:end);
phiz1=Dz1*phi1; phiz2=Dz2*phi2; phiz3=Dz3*phi3;
% flux:
% JNaz1=-Az1.*DNa.*(Naz1+zNa*Na1.*phiz1);
% JClz1=-Az1.*DCl.*(Clz1+zCl*Cl1.*phiz1);
% Jqz2=-Az2.*Dq.*(qz2+zq*q2.*phiz2+gij*q2.*qz2+q2.*Vz);
% JNaz3=-Az3.*DNa.*(Naz3+zNa*Na3.*phiz3);
% JClz3=-Az3.*DCl.*(Clz3+zCl*Cl3.*phiz3);

%%%%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
JNaz1=-Az1.*DNa.*(Naz1+zNa*Na1.*phiz1);
JClz1=-Az1.*DCl.*(Clz1+zCl*Cl1.*phiz1);
JNaz3=-Az3.*DNa.*(Naz3+zNa*Na3.*phiz3);
JClz3=-Az3.*DCl.*(Clz3+zCl*Cl3.*phiz3);
JAqz1=-Az2.*Dq1.*(Aqz1+zq1*Aq1.*phiz2+gij*Aq1.*Aqz1+Aq1.*Vz);
JAqz2=-Az2.*Dq2.*(Aqz2+zq2*Aq2.*phiz2+gij*Aq2.*Aqz2+Aq2.*Vz);
JAqz3=-Az2.*Dq3.*(Aqz3+zq3*Aq3.*phiz2+gij*Aq3.*Aqz3+Aq3.*Vz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% conservation law:
dNadt1=-DAz1*JNaz1; dNadt3=-DAz3*JNaz3; dCldt1=-DAz1*JClz1; dCldt3=-DAz3*JClz3;
%dqdt2=-DAz2*Jqz2;
%%%%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
dqdt1=-DAz2*JAqz1;
dqdt2=-DAz2*JAqz2;
dqdt3=-DAz2*JAqz3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BC:
dNadt1(1)=Na1(1)-NaL;
dNadt1(end)=JNaz1(end);
dNadt3(1)=JNaz3(1);
dNadt3(end)=Na3(end)-NaR;
dCldt1(1)=Cl1(1)-ClL;
dCldt1(end)=JClz1(end);
dCldt3(1)=JClz3(1);
dCldt3(end)=Cl3(end)-ClR;
%dqdt2(1)=Jqz2(1); dqdt2(end)=Jqz2(end);
% dudt=[dNadt1;dCldt1;dqdt2;dNadt3;dCldt3];

%%%%%%%%%%%%%%%%%%%%%%Revised by XU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
dqdt1(1)=JAqz1(1); dqdt1(end)=JAqz1(end);
dqdt2(1)=JAqz2(1); dqdt2(end)=JAqz2(end);
dqdt3(1)=JAqz3(1); dqdt3(end)=JAqz3(end);
dudt=[dNadt1;dCldt1;dqdt1;dqdt2;dqdt3;dNadt3;dCldt3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[t phiRactual max(abs(dudt))]
if iplot==1 & t>0
    figure(5); subplot(1,2,1);
    plot(z1,Na1,'b','LineWidth',2); hold on; plot(z3,Na3,'b','LineWidth',2); plot(z1,Cl1,'g','LineWidth',2); plot(z3,Cl3,'g','LineWidth',2);
%     plot(z2,q2,'r','LineWidth',2); 
    plot(z2,Aq1,'r','LineWidth',2); plot(z2,Aq2,'y','LineWidth',2); plot(z2,Aq3,'k','LineWidth',2); 
    hold off; grid on; xlabel('z'); ylabel('Na, Cl, q'); title(['t=' num2str(t)]);
    subplot(1,2,2);
    plot(z1,phi1,'k','LineWidth',2); hold on; plot(z2,phi2,'k','LineWidth',2); plot(z3,phi3,'k','LineWidth',2);
    grid on; xlabel('z'); ylabel('\phi'); hold off; drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEB  compute D = differentiation matrix, x = Chebyshev grid

function [D,x] = cheb(N)
if N==0, D=0; x=1; return, end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
D  = D - diag(sum(D'));                 % diagonal entries
