
%%  Research Topic : Modeling multilayer Building Envelopes Explicit Scheme, Temperature boundary                    
%%  Authors : Sajith Wijesuriy
%%  Institution : NREL 
%%  Research Group : 
%%  File Feature : Finite Difference Code while calling in supportive files 
%%  Current Author: Sajith Wijesuriya 
clear all
close all
clc
%%  Notes (by Sajith Wijesuriya)
% Solves the 1D heat equation with an explicit finite difference scheme
%%  2.PCM recognition and the material property mfile
% This script defines the layer material properties as user requires. It can,
% 1. Set any number of material layers
% 2. Manually input material data layer by layer or read from a pre defined
% excel sheet 
% 3. Request the user if there is a need to include a PCM layer and read
% data from a file. 
% 4. Differentiate between IDEAL and real PCMs based on the user
% requirement and relate to relavant enthalpy curves to calculate Cp 
%--------------------------------------------------------------------------
%% Run the data files 
%% Take basic inputs (queries and file reads)
ft2_m2=0.092903;
WWR=0.15; 
% dt=1;
% file 
[num,txt,raw] = xlsread('material_property_pcm.xlsx');
%% direct input
% numberoflayers = 'How many layers in the envelope? ';
% disp('------')
% numlay = input(numberoflayers);  
% disp(numlay)
numlay =3;
%% same or different materials? 
% samemat='If same material is in all layers 1, if not 0:';
% disp('------')
% samemat=input(samemat);
% if samemat==1
%     assemblytype='same_material';
% elseif samemat==0
%     assemblytype='diff_material';
% end
% disp('------') 
%% Temporal settings 
% dtime = 'Time step (s)? '; % define the time step (s)
% dt= input(dtime);
dt= 60;
% steps= 'Simulation steps (number)? '; % number of steps to simulate 
% sim_length= input(steps);
sim_length= 44640;
%% File inputs: 
filenamexpdata='Heat_transfer_balance_noPCMccnr.xlsx';
sim.extbc= xlsread(filenamexpdata,'P2:P44641')+273.15;
sim.intbc= xlsread(filenamexpdata,'K2:K44641')+273.15;
sim.intrf= xlsread(filenamexpdata,'I2:I44641')+273.15;
sim.intair1= xlsread(filenamexpdata,'N2:N44641')+273.15;
% Interior radiation 
sim.Qfsr= xlsread(filenamexpdata,'B2:B44641');% Output:Variable,*,Surface Inside Face Solar Radiation Heat Gain Rate per Area, hourly; !- Zone Average [W/m2]
sim.Qflr= xlsread(filenamexpdata,'C2:C44641');% Output:Variable,*,Surface Inside Face Lights Radiation Heat Gain Rate per Area, hourly; !- Zone Average [W/m2]
sim.Qfigr= xlsread(filenamexpdata,'D2:D44641');% Output:Variable,*,Surface Inside Face Internal Gains Radiation Heat Gain Rate per Area, hourly; !- Zone Average [W/m2]
sim.Qfstr= xlsread(filenamexpdata,'F2:F44641');% Output:Variable,*,Surface Inside Face Net Surface Thermal Radiation Heat Gain Rate per Area, hourly; !- Zone Average [W/m2]
sim.Qconvin= xlsread(filenamexpdata,'E2:E44641');% Output:Variable,*,Surface Inside Face Convection Heat Gain Rate per Area, hourly; !- Zone Average [W/m2]
sim.Qcondin= xlsread(filenamexpdata,'G2:G44641');% Output:Variable,*,Surface Inside Face Conduction Heat Transfer Rate per Area,hourly; !- Zone Average [W/m2]
% sim.Qfsref=sim.Qfsr/0.5;
% sim.Qflref=sim.Qflr/0.9;
% sim.Qfigref=sim.Qfigr/0.9;
% sim.Qfstref=sim.Qfstr/0.7;
% sim.Qfsref=sim.Qfsr*0.5;
% sim.Qflref=sim.Qflr*0.9;
% sim.Qfigref=sim.Qfigr*0.9;
% sim.Qfstref=sim.Qfstr*0.9;
% sim.Qradhf=sim.Qfsref+sim.Qflref+sim.Qfigref+sim.Qfstref;
sim.Qradhf=sim.Qfsr+sim.Qflr+sim.Qfigr+sim.Qfstr;
% sim.Qradhf=sim.Qfstr;
sim.Qallinhf=sim.Qfsr+sim.Qflr+sim.Qfigr+sim.Qfstr+sim.Qconvin;
sim.hconvin= xlsread(filenamexpdata,'S2:S44641');
% allocating material properties 
%% allocating material properties 
% defined
Cp_air=1.0035;% J/kg-K 
k_air=0.0248;% W/m-K
rho_air=1.225; % kg/m3
for i=1:numlay 
                                                                                                                                    lay(i).entpcr=num(8,i);  %reads existance
                                                                                                                                    lay(i).Tc_pcr=num(9,i);  %reads existance
                                                                                                                                    lay(i).dT_pcr=num(10,i); %reads existance
                                                                                                                %% Now to set up each layer properties other than PCMs 
                                                                                                                                    lay(i).winA=156*ft2_m2;
                                                                                                                                    lay(i).A=lay(i).winA/WWR; 
                                                                                                                                    lay(i).t=num(1,i);%m                                   
                                                                                                                                    lay(i).dens=num(2,i);                                     
                                                                                                                                    lay(i).Cp=num(3,i);                                  
                                                                                                                                    lay(i).k=num(4,i);
                                                                                                                                    lay(i).entpc=num(5,i);
                                                                                                                                    lay(i).Tc_pc=num(6,i); 
                                                                                                                                    lay(i).dT_pc=num(7,i);
                                                                                                                                    lay(i).nodes=10;
                                                                                                                                    lay(i).dx=lay(i).t/lay(i).nodes;      
                                                                                                                                    lay(i).emis=0.9;              %[] Emissivity of external wall
                                                                                                                                    lay(i).abs =0.3;              %[] Absorptivity of external wall  
                                                                                                                                    % lay(i).G = G.vert;            % Solar radiation (beam + diffuse + gnd reflect)
                                                                                                                                    lay(i).L=sqrt(lay(i).A);
end
% small PCM recognition segment 
for  i=1:numlay
if lay(i).entpcr>0
    pcmrec(i,1)=1;
    else
    pcmrec(i,1)=0;
end 
end
% set up layers
if numlay==1
layers=[lay(1)];
elseif numlay==2
layers=[lay(1),lay(2)];
elseif numlay==3
layers=[lay(1),lay(2),lay(3)];
elseif numlay==4
layers=[lay(1),lay(2),lay(3),lay(4)];    
elseif numlay==5
layers=[lay(1),lay(2),lay(3),lay(4),lay(5)];   
elseif numlay==6
layers=[lay(1),lay(2),lay(3),lay(4),lay(5),lay(6)];  
end
% calling discretization calculator (will set new number of nodes per layer and the new dx unless number of nodes are pre defined) 
layers=dxCALCex(layers,dt);
% disp(layers)
for i=1:numlay
    disp(lay(i).nodes)
    disp(lay(i).dx)
end
disp(lay(3))
%% Total Nodes  
% Aggrgate nodes represetns the total number of the nodes in the envelope assembly.
% This is calculated by adding nodes of each alyer and accounting for the interface nodes which should not be double counted
TotNodes=0;% initializing tyhe total nodes to zero here 
for i = 1:length(layers)    
    TotNodes=TotNodes+layers(i).nodes;
end
TotNodes=TotNodes-(length(layers)-1);
x=zeros(1,TotNodes);
node=1;
for i = 1:length(layers)
    if i==1
        for j=1:layers(i).nodes
            x(node)=layers(i).dx*(j-1);
            node=node+1;            
        end
    else
        for j = 2:layers(i).nodes
            x(node)=x(node-1)+layers(i).dx;
            node=node+1;
        end
    end
end
%% Setup initial temperature profile
Tinit=20+273.15;
T=zeros(TotNodes,sim_length+1);
% T=zeros(TotNodes,1); 
% Tnew=zeros(TotNodes,1);
T(:,:)=Tinit;
disp(T(:,1))
%% Boundary condition: Temperature profile 
% % method 1
% Text=(40+273.15);
% Tint=(30+273.15);
% % method 2
% Text=(40+273.15)*ones(1,sim_length+1);
% Tint=(30+273.15)*ones(1,sim_length+1);
% method 3 
Text=sim.extbc';
Tint=sim.intbc';
Tia=sim.intair1';
sigma=5.67*10^(-8);
Qradint=sim.Qradhf';
% Qradint=0*ones(1,sim_length);
Qconvint=sim.Qconvin';
timec=0;
% Check: Meeting the stability criteria 
disp(layers(1,1).dx)
disp(layers(1,2).dx)
disp(layers(1,3).dx)
disp(((dt*layers(1,1).k/(layers(1,1).Cp*layers(1,1).dens)))/layers(1,1).dx^2)
disp(((dt*layers(1,2).k/(layers(1,2).Cp*layers(1,2).dens)))/layers(1,2).dx^2)
disp(((dt*layers(1,3).k/(layers(1,3).Cp*layers(1,3).dens)))/layers(1,3).dx^2)
disp((((dt/(layers(1,3).Cp*layers(1,3).dens)))/layers(1,3).dx)*((layers(1,3).k/layers(1,3).dx)+6))
%% Heat transfer coefficients: Interior boundary 
% convection heat transfer 
h_convi=sim.hconvin';
% h_convi=5*ones(1,sim_length); % convection heat transfer coeeficient 
% longwave heat transfer 
Q_fac=0*ones(1,sim_length);
Qconvicalc=ones(1,sim_length);
%% Counters and timers 
% Counters and timers help to measure the performance of the code (spped,
% accuracy) 
count = 0;
%% Number of warmup loops 
% Warm up loops are important when the program needs to start after a
% training time when it runs a transient solution. Therfore, we need to
% allocate enough time steps for training loop.
% -------------------------------------------------------------------------
warmup=5;       % Number of warm up loops ( EnergyPlus runs 50 days of warm-up)
converge=25;    % Maximum number of convergence iterations in Gauss Siedel loop  
limit=1e-5;     % Residual Minimum for convergence loop 

%% PCM INCLUSIONS 
                                     for lay=1:length(layers)
                                                          if pcmrec(lay)==1
                                                                                    %% Setting the basic phase change material properties
                                                                                    disp('pcm inclusions')
                                                                                    % layers(1,lay).dE_pcm=layers(1,lay).entpc; % per unit mass storage is brought in 
                                                                                    k=9;% how many points to define the curve discretization 
                                                                                    Ecf=zeros(k,1);% curve fit                                                                                    
                                                                                    kcf=zeros(k,1);% curve fit 
                                                                                    Tcf=zeros(k,1);                                                                                 
                                                                                    %% Phase change material properties are used to calculate an enthalpy curve 
                                                                                    Tc_pc=layers(1,lay).Tc_pc; % the center value of the phase change 
                                                                                    dT_pc=layers(1,lay).dT_pc; % phase change temperature interval 
                                                                                    Th_pc=Tc_pc+(dT_pc/2); % where melting ends and solidification starts 
                                                                                    Tl_pc=Tc_pc-(dT_pc/2); % where solidication ends and melting starts
                                                                                    %% Separate calculation of the curves by discretization 
                                                                                    % until the 20000 base aggregate/sum enthalpy it is assumed that the cp=
                                                                                    % 2500 and then at the starting point of the temperature the cp converts to
                                                                                    % become cph if PCM is present. 
%                                                                                     for i= 1:k 
%                                                                                     Ecf(i)=20000+((i-1)*(layers(1,lay).entpc/(k-1))); %20000 s a base value set to start the curve slope 
%                                                                                     Tcf(i)=Tl_pc + ((Th_pc-Tl_pc)/(k-1))*(i-1); 
%                                                                                     end 
                                                                                   CPH=layers(1,lay).entpc/dT_pc;
                                                              elseif pcmrec(lay)==0
                                                                                       disp('no pcm inclusions')
                                                              end
                                     end
%% Initialize wait bars 
% waitbars indicate the temporal progress of the code and gives an
% indication that the code is in progress. 
wb1 = waitbar(0,'Model is progressing and iterating');             % Progress Update
wb2 = waitbar(0,'Warm-up status');                                 % Warm up status update Progress Update
for init=1:warmup  
    % init runs a waitbar  
    waitbar(init/warmup,wb2)    
    % why this condition? 
    % for the first warmup-1 loops the counter is set to 100 and runs
    % smaller number of time steps (100*(warmup-)) number of total steps
    % before the final run. 
    if init==1
        counter=100; 
        T(:,:) = Tinit;
    elseif init~=warmup
        counter=100;      
        T(:,1)=T(:,counter+1);
    else
        counter=sim_length;    
        T(:,1)=T(:,counter+1);
    end  
for n=1:sim_length % Timestep loop
    waitbar(n/sim_length,wb1); %Waitbar update    
    count=count+1;
% Convergence arrays initialized   
    Tconverge_P_e=zeros(TotNodes,converge);
    Tconverge_P_e(:,:)=27+273.15; %In C    
for j=1:converge
     i=1; % i is set to 1 here to assign an initial value 
    for lay = 1:length(layers)
%% Compute new temperature: method 3       
        for node = 1:layers(1,lay).nodes   
                dx=layers(1,lay).dx;
%                 (disp(dx))
                k=layers(1,lay).k;
                cp=layers(1,lay).Cp;                
                dens=layers(1,lay).dens;  
            if node == 1 && lay==1

               if  pcmrec(lay)==1
                if n>1  
                if T(i,n)>Tl_pc && T(i,n)<Th_pc
                        Cp_pc=CPH;
                        elseif T(i,n)<=Tl_pc
                        Cp_pc=layers(1,lay).Cp;
                        elseif T(i,n)>=Th_pc
                        Cp_pc=layers(1,lay).Cp;
                end
                end      
              elseif pcmrec(lay)==0
                Cp_pc=cp;
              end   
            T(i,n+1) = Text(1,n);
            i=i+1;
            elseif node > 1 && node < layers(1,lay).nodes  
                if  pcmrec(lay)==1
                if n>1  
                if T(i,n)>Tl_pc && T(i,n)<Th_pc
                        Cp_pc=CPH;
                        elseif T(i,n)<=Tl_pc
                        Cp_pc=layers(1,lay).Cp;
                        elseif T(i,n)>=Th_pc
                        Cp_pc=layers(1,lay).Cp;
                end
                end      
              elseif pcmrec(lay)==0
                Cp_pc=cp;
              end 
            T(i,n+1) = (((Cp_pc*layers(1,lay).dx*layers(1,lay).dens/dt)*T(i,n))+((layers(1,lay).k*T(i-1,n+1))+(layers(1,lay).k*T(i+1,n+1)))/layers(1,lay).dx)/(((layers(1,lay).k+layers(1,lay).k)/layers(1,lay).dx)+Cp_pc*layers(1,lay).dx*layers(1,lay).dens/dt);  
            i=i+1;
            elseif node == layers(1,lay).nodes && lay ~= length(layers)  
                if  pcmrec(lay)==1
                if n>1  
                if T(i,n)>Tl_pc && T(i,n)<Th_pc
                        Cp_pc=CPH;
                        elseif T(i,n)<=Tl_pc
                        Cp_pc=layers(1,lay).Cp;
                        elseif T(i,n)>=Th_pc
                        Cp_pc=layers(1,lay).Cp;
                end
                end      
              elseif pcmrec(lay)==0
                Cp_pc=cp;
                end 
                if  pcmrec(lay+1)==1
                if n>1  
                if T(i,n)>Tl_pc && T(i,n)<Th_pc
                        Cp_pc2=CPH;
                        elseif T(i,n)<=Tl_pc
                        Cp_pc2=layers(1,lay+1).Cp;
                        elseif T(i,n)>=Th_pc
                        Cp_pc2=layers(1,lay+1).Cp;
                end
                end      
              elseif pcmrec(lay)==0
                Cp_pc2=layers(1,lay+1).Cp;
                end                                                                                                                                                     
                    dx2=layers(1,lay+1).dx;
                    k2=layers(1,lay+1).k;               
                    dens2=layers(1,lay+1).dens;
                    cp1_fac=Cp_pc*dx^2*dx2*dens;
                    cp2_fac=Cp_pc2*dx2^2*dx*dens2;                                                                                                                             
                    cp_fac_P_e=cp1_fac+cp2_fac;                                                                                                                              
            T(i,n+1)=(2.0*(dt*dx2*layers(1,lay).k*T(i-1,n+1)+dt*layers(1,lay).dx*k2*T(i+1,n+1))+cp_fac_P_e*T(i,n))/(2.0*(dt*dx2*layers(1,lay).k+ dt*layers(1,lay).dx*k2)+cp_fac_P_e);       
            i=i+1;
            elseif node == layers(1,lay).nodes && lay == length(layers)
            T(i,n+1) = Tint(1,n);
%             T(i,n+1) =(((2*dt*layers(1,lay).dx)*(Qradint(1,n)+(h_convi(1,n)*Tia(1,n))))+((layers(1,lay).Cp*layers(1,lay).dx^2*layers(1,lay).dens)*T(i,n))+(2*dt*layers(1,lay).k*T(i-1,n)))/((2*dt*layers(1,lay).dx*h_convi(1,n))+(2*dt*layers(1,lay).k)+(layers(1,lay).Cp*layers(1,lay).dx^2*layers(1,lay).dens));
%             Qconvicalc(1,n)=(h_convi(1,n)*(Tia(1,n)-T(i,n)));
            end
        end 
%% Update temperature and time
% T = Tnew;
timec = timec+dt;

    end
            Tconverge_P_e(:,j+1)=T(:,n+1);
        if j>=5
            Tconverge_P_e(:,j+1)=Tconverge_P_e(:,j)-((Tconverge_P_e(:,j)-Tconverge_P_e(:,j+1))*0.5); % Under/Over relaxation 
        elseif j>=10 && j<15
            Tconverge_P_e(:,j+1)=Tconverge_P_e(:,j)-((Tconverge_P_e(:,j)-Tconverge_P_e(:,j+1))*0.25); % Under/Over relaxation 
        elseif j>=15
            Tconverge_P_e(:,j+1)=Tconverge_P_e(:,j)-((Tconverge_P_e(:,j)-Tconverge_P_e(:,j+1))*0.125); % Under/Over relaxation  
        end
        % setup the exit convergence limits
        sumdelT_pc=sum(Tconverge_P_e(:,j+1)-Tconverge_P_e(:,j));
        sumT_pc=sum(Tconverge_P_e(:,j));
        ConvergeDiff_P_e=abs(sumdelT_pc/sumT_pc);  
        if ConvergeDiff_P_e<=limit
            break
        end
        J(n)=j;
end
end
end

% hradint=hrad(TotNodes,:);
% hrad2int=hrad2(TotNodes,:);
thickenv=x(1,TotNodes); 
time=zeros(1,sim_length);
for n=1:sim_length
    time(n)=n;
end
%% Reading to plot 
Tint_sim=T((layers(1).nodes+layers(2).nodes-1),1:44640)-273.15;
Tint_simwrite=Tint_sim';

Tintb_sim=T(TotNodes,1:44640)-273.15;
Tintb_simwrite=Tintb_sim';
%% Writing to file
% simdataread='simdatareadT.xlsx';
% xlswrite(simdataread,(Tint_simwrite(1:44640,1)),'B1:B44641')
%% File inputs for verification
% filenamesimdata='output_comparison.xlsx';
% sim.epintr= xlsread(filenamesimdata,'C1:C44641');
% sim.epintr=sim.epintr';
%% Conversions
sim.extbcc=sim.extbc-273.15;
sim.intrfc=sim.intrf-273.15;
sim.intbcc=sim.intbc-273.15;
%% Accuracy 
RMSE_T_ibc = (sqrt(mean(((T(TotNodes,1:44640))' -  (sim.intbc(1:44640,1))).^2)));   % Root Mean Squared Er
RMSE_T = (sqrt(mean(((T(layers(1).nodes+layers(2).nodes-1,1:44640))' -  (sim.intrf(1:44640,1))).^2)));   % Root Mean Squared Er

disp('RMSE value at ib')
disp(RMSE_T_ibc)
disp('RMSE value')
disp(RMSE_T)
% CVRMSE_T_insdw = 100*(sqrt(mean(((T(layers(1).nodes,1:44640))' -  (exp.Tplypoup3psaiave(2641:4561,1))).^2)))/mean(exp.Tplypoup3psaiave(2641:4561,1));   % Root Mean Squared Error
% NMBE_T_insdw  = 100*(mean(((T(layers(1).nodes,1:44640))' -  (exp.Tplypoup3psaiave(2641:4561,1)))))/mean(exp.Tplypoup3psaiave(2641:4561,1));   % Root Mean Squared Error
%% color schemes for plotting 
lg1   = [0.85 0.85 0.85];
lg2   = [0.7 0.7 0.7];
dg1  = [0.4 0.4 0.4];
dg2  = [0.2 0.2 0.2];
b2=[0,0.4470, 0.7410];
g2=[0.4660, 0.6740, 0.1880];
pd=[0.75,0,0.75];
pdd=[0.4940, 0.1840, 0.5560];
bl=[0, 0.75, 0.75];
yd=[0.9290, 0.6940, 0.1250];
ydd=[0.8500, 0.3250, 0.0980]; 
grd=[0, 0.5, 0];
p2=[0.75, 0, 0.75];
pur=[0.4940, 0.1840, 0.5560];
mar=[0.6350, 0.0780, 0.1840];
%% Plots
figure
set(gcf,'color','w')
set(gca,'FontSize',20)
grid off
box on
% xlim([97 168]) % From
% xlim([10 110])
% ylim([0 60])    
ax = gca;
ax.XLim = [1 10080];
% ax.XTick = [1 1800 3000 4200 5400 44640];
% ax.XTickLabel = {'1','1200','2400','3600','4800','44640'};
% title('HF')
xlabel('Time (minutes)')
ylabel('Temperature(\circ{C})')
% legend={'No CLT Rockwool last presentation','No CLT rockwool','No CLT XPS,CLT Rockwool','CLT XPS','STUD with Rockwool'};
hold on 
% %% Setpoint 
% plot(time,sim.TRWSP,'Color','k','LineWidth',4)
%% BCs
% plot(time(1:30:end),sim.extbcc(1:30:end),':','Color','k','LineWidth',4,'MarkerSize',10)
%% Experimental 
plot(time(1:30:end),sim.intrfc(1:30:end),'-','Color',b2,'LineWidth',4,'MarkerSize',10)
%% MATLAB 
plot(time(1:30:end),Tint_sim(1:30:end),'-','Color',ydd,'LineWidth',4,'MarkerSize',10)
%% EP
% plot(time(1:30:end),sim.epintr(1:30:end),'-','Color',b2,'LineWidth',4,'MarkerSize',10)
%% BCs
% plot(time(1:30:end),sim.intbcc(1:30:end),'--','Color','k','LineWidth',4,'MarkerSize',10)
%% room 
% plot(sim.intair1-273.15,'b--o','Color','k','LineWidth',2,'MarkerSize',6)
hold off 
% lgd=legend( '1. Exterior surface', '2. Experimental (plywood-PCM)','3. CSMPCM (plywood-PCM)','4. Experimental (PCM-drywall)','5. CSMPCM (PCM-drywall)','6. Interior surface');
lgd=legend( '1. EnergyPlus', '2. MATLAB alg.');

figure
set(gcf,'color','w')
set(gca,'FontSize',20)
grid off
box on
% xlim([97 168]) % From
% xlim([10 110])
% ylim([0 60])    
ax = gca;
ax.XLim = [1 10080];
% ax.XTick = [1 1800 3000 4200 5400 44640];
% ax.XTickLabel = {'1','1200','2400','3600','4800','44640'};
% title('HF')
xlabel('Time (minutes)')
ylabel('Temperature(\circ{C})')
% legend={'No CLT Rockwool last presentation','No CLT rockwool','No CLT XPS,CLT Rockwool','CLT XPS','STUD with Rockwool'};
hold on 
% %% Setpoint 
% plot(time,sim.TRWSP,'Color','k','LineWidth',4)
%% BCs
% plot(time(1:30:end),sim.extbcc(1:30:end),':','Color','k','LineWidth',4,'MarkerSize',10)
%% Experimental 
plot(time(1:30:end),sim.intbcc(1:30:end),'-','Color',b2,'LineWidth',4,'MarkerSize',10)
%% MATLAB 
plot(time(1:30:end),Tintb_sim(1:30:end),'-','Color',ydd,'LineWidth',4,'MarkerSize',10)
%% EP
% plot(time(1:30:end),sim.epintr(1:30:end),'-','Color',b2,'LineWidth',4,'MarkerSize',10)
%% BCs
% plot(time(1:30:end),sim.intbcc(1:30:end),'--','Color','k','LineWidth',4,'MarkerSize',10)
%% room 
% plot(sim.intair1-273.15,'b--o','Color','k','LineWidth',2,'MarkerSize',6)
hold off 
% lgd=legend( '1. Exterior surface', '2. Experimental (plywood-PCM)','3. CSMPCM (plywood-PCM)','4. Experimental (PCM-drywall)','5. CSMPCM (PCM-drywall)','6. Interior surface');
% lgd=legend( '1. Exterior surface', '2. Experimental (plywood-PCM)','3. CSMPCM (plywood-PCM)','4. Experimental (PCM-drywall)','5. CSMPCM (PCM-drywall)','6. Interior surface');
lgd=legend( '1. EnergyPlus', '2. MATLAB alg.');

figure
set(gcf,'color','w')
set(gca,'FontSize',20)
grid off
box on
% xlim([97 168]) % From
% xlim([10 110])
% ylim([0 60])    
ax = gca;
ax.XLim = [1 10080];
% ax.XTick = [1 1800 3000 4200 5400 44640];
% ax.XTickLabel = {'1','1200','2400','3600','4800','44640'};
% title('HF')
xlabel('Time (minutes)')
ylabel('HF(W/m^2)')
% legend={'No CLT Rockwool last presentation','No CLT rockwool','No CLT XPS,CLT Rockwool','CLT XPS','STUD with Rockwool'};
hold on 
%% 
plot(time(1:30:end),sim.Qfsr(1:30:end),'-','Color',b2,'LineWidth',4,'MarkerSize',10)% Face System Radiation Heat Gain 
%% 
plot(time(1:30:end),sim.Qflr(1:30:end),'-','Color',ydd,'LineWidth',4,'MarkerSize',10)% Face Lights Radiation Heat Gain Rate per Area
%% 
plot(time(1:30:end),sim.Qfigr(1:30:end),'-','Color','r','LineWidth',4,'MarkerSize',10)% Face Internal Gains Radiation Heat Gain 
%% 
plot(time(1:30:end),sim.Qfstr(1:30:end),'-','Color','c','LineWidth',4,'MarkerSize',10) % Surface Thermal Radiation Heat Gain 
%% 
plot(time(1:30:end),sim.Qconvin(1:30:end),'-','Color','y','LineWidth',4,'MarkerSize',10) % convections
%% 
plot(time(1:30:end),sim.Qallinhf(1:30:end),'-','Color','g','LineWidth',4,'MarkerSize',10) % All above heat fluxes 
%% 
plot(time(1:30:end),-sim.Qcondin(1:30:end),'-','Color','k','LineWidth',4,'MarkerSize',10) % conductions
%% 
plot(time(1:30:end),Qconvicalc(1:30:end),'-','Color',pur,'LineWidth',4,'MarkerSize',10) % conductions
%% 
plot(time(1:30:end),sim.Qradhf(1:30:end),'--','Color','k','LineWidth',2,'MarkerSize',10) % conductions
%% room 
% plot(sim.intair1-273.15,'b--o','Color','k','LineWidth',2,'MarkerSize',6)
hold off 
% lgd=legend('Exterior surface', 'Experimental (ply-pou)','Experimental(pou-dry)','MATLAB (ply-pou)','MATLAB (pou-dry)','EnergyPlus Hysteresis(ply-pou)','EnergyPlus Hysteresis(pou-dry)','WUFI(ply-pou)','WUFI(pou-dry)','Interior surface');
% lgd=legend('Radiant wall supply setpoint', 'Exterior surface', 'Experimental (ply-pou)','Experimental (pou-dry)','MATLAB (ply-pou)','MATLAB (pou-dry)','Interior surface','Chamber temperature');
% lgd=legend( '1. Exterior surface', '2. Experimental (plywood-PCM)','3. CSMPCM (plywood-PCM)','4. Experimental (PCM-drywall)','5. CSMPCM (PCM-drywall)','6. Interior surface');
% lgd=legend( '1. Observed', '2. Calculated');




figure
set(gcf,'color','w')
set(gca,'FontSize',20)
grid off
box on
% xlim([97 168]) % From
% xlim([10 110])
% ylim([0 60])    
ax = gca;
ax.XLim = [1 10080];
% ax.XTick = [1 1800 3000 4200 5400 44640];
% ax.XTickLabel = {'1','1200','2400','3600','4800','44640'};
% title('HF')
xlabel('Time (minutes)')
ylabel('HF(W/m^2)')
% legend={'No CLT Rockwool last presentation','No CLT rockwool','No CLT XPS,CLT Rockwool','CLT XPS','STUD with Rockwool'};
hold on 
%% 
% plot(time(1:30:end),sim.Qfsr(1:30:end),'-','Color',b2,'LineWidth',4,'MarkerSize',10)% Face System Radiation Heat Gain 
%% 
% plot(time(1:30:end),sim.Qflr(1:30:end),'-','Color',ydd,'LineWidth',4,'MarkerSize',10)% Face Lights Radiation Heat Gain Rate per Area
%% 
% plot(time(1:30:end),sim.Qfigr(1:30:end),'-','Color','r','LineWidth',4,'MarkerSize',10)% Face Internal Gains Radiation Heat Gain 
%% 
% plot(time(1:30:end),sim.Qfstr(1:30:end),'-','Color','c','LineWidth',4,'MarkerSize',10) % Surface Thermal Radiation Heat Gain 
%% 
% plot(time(1:30:end),sim.Qconvin(1:30:end),'-','Color','y','LineWidth',4,'MarkerSize',10) % convections
%% 
plot(time(1:30:end),sim.Qallinhf(1:30:end),'-','Color','g','LineWidth',4,'MarkerSize',10) % All above heat fluxes 
%% 
plot(time(1:30:end),-sim.Qcondin(1:30:end),'-','Color','k','LineWidth',4,'MarkerSize',10) % conductions
%% 
% plot(time(1:30:end),Qconvicalc(1:30:end),'-','Color',pur,'LineWidth',4,'MarkerSize',10) % conductions
%% 
% plot(time(1:30:end),h_convi(1:30:end),'-','Color','k','LineWidth',4,'MarkerSize',10) % conductions
%% room 
% plot(sim.intair1-273.15,'b--o','Color','k','LineWidth',2,'MarkerSize',6)
hold off 
% lgd=legend('Exterior surface', 'Experimental (ply-pou)','Experimental(pou-dry)','MATLAB (ply-pou)','MATLAB (pou-dry)','EnergyPlus Hysteresis(ply-pou)','EnergyPlus Hysteresis(pou-dry)','WUFI(ply-pou)','WUFI(pou-dry)','Interior surface');
% lgd=legend('Radiant wall supply setpoint', 'Exterior surface', 'Experimental (ply-pou)','Experimental (pou-dry)','MATLAB (ply-pou)','MATLAB (pou-dry)','Interior surface','Chamber temperature');
% lgd=legend( '1. Exterior surface', '2. Experimental (plywood-PCM)','3. CSMPCM (plywood-PCM)','4. Experimental (PCM-drywall)','5. CSMPCM (PCM-drywall)','6. Interior surface');
lgd=legend( '1. Radiation and convection', '2. Conduction');

