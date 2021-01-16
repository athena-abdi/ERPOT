function optimization_ERPOT_BUS
clc
clear all
close all
%##########     MPSOC Charatcteristics     ###############%
No_bins = 10;

% #### Objectives Definition #### % 

% Set the power Objectives Here!

% Set the Temperature Objectives Here!

% Set the Lambda (failure rate) Objectives Here!
 

final=zeros(No_bins,No_bins);  %final optimization results

%% Calling the List Scheduling Algoritm for different levels of constraints to build the Pareto front
for tet=1:No_bins
    for lamb=1:No_bins
        for pow=1:No_bins
           [final(tet, lamb,pow)] = ERPOT_NEW_BUS(teta_obj(1,tet),lambda_obj(1,lamb),power_obj(1,pow));
           cmax(tet, lamb,pow)=final;
           teta(tet, lamb,pow)=theta;
           lambda(tet, lamb,pow)=lambda*1e9;
           power_total(tet, lamb,pow)=power_final;
	end
    end
end
%% Applying the Pareto elimination on the result and then plot the final Pareto front
scatHand = scatter3(teta, lambda, power_total,'fill');
set(scatHand, 'CData', Cmax);
length(teta)
length(lambda)
length(power_total)
length(Cmax)
xlabel ('Temperature (K)');
ylabel ('GSFR');
zlabel ('Power Consumption');

% 3D figure based on two constarints and objective:
% figure,surf(teta,lambda,Cmax);