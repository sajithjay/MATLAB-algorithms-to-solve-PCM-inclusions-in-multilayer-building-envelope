
%% Research Topic : Modeling multilayer Building Envelopes                     

%%  Author : Sajith Wijesuriya, Nick Kincaid, Paulo Tabares  
%%  Institution : Department of Mechanical Engineering, Colorado School of Mines 
%%  Research Group : Dr. Tabares Research Group
%%  File Feature : Discreetization intended for the Finite Difference Code 

function [ layers ] = dxCALCex(layers,dt)
%   Inputs layers calculates optimal number of nodes and dx for each layer
disp('------')
% c= 'add multiplier for discretization:';
% c=input(c);    %Coefficient of discretizaiton ( can be the multiplier to increase dx and therfoore, reduce #of nodes per layer) 
c=1;
Fo=1/c;        %Coefficient of discretizaiton relation to the fourier number 
disp('------')
nodeset='Manually set the nodes?';
nodeset=input(nodeset);
if nodeset==1
    nodesetup='manual_nodes';
elseif nodeset==0
    nodesetup='optimized_nodes';
end

if isequal(nodesetup,'optimized_nodes')==1
%     disp(length(layers))
for lay = 1:length(layers)
    %% Calculation of optimal dx
    k=layers(lay).k;
    cp=layers(lay).Cp;
    dens=layers(lay).dens;
    thick=layers(lay).t;    
    alpha=k/(cp*dens);
    dx=sqrt(alpha*dt)*c; %% note on multiplier
    %% Reassign number of nodes and dx to layers
    layers(lay).nodes=ceil(thick/dx);
    layers(lay).dx=thick/(layers(lay).nodes-1);  
    disp(layers(1).dx)
end 
elseif isequal(nodesetup,'manual_nodes')==1
for lay = 1:length(layers) 
    thick=layers(lay).t;  
    layers(lay).nodes= 'nodes this layer:';
    layers(lay).nodes=input(layers(lay).nodes);
    layers(lay).dx=thick/(layers(lay).nodes-1);
end    
end


end

