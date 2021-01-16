function [teta_final]=cooling_ERPOT_Bus (teta,P_dyn,teta_n1,teta_bus)
%% Initialization
C = 0.03;
alpha = 0.1;
G = 0.3;
teta_amb=293;
beta_cool =-40;  %NORMAL
t_slack=1;
Gn=0.1;
a=(G+(2*Gn)-alpha)/C;
%% Computation
% teta_cool = teta_ss_cool+((teta-teta_ss_cool)*(exp(-(a)*(((t_slack)*0.001)))));
% teta_final=teta_cool;
b=(P_dyn+beta_cool+(G*teta_amb)+(Gn*teta_n1)+(Gn*teta_bus))/C;
teta_ss = b/a; 
teta_cool=teta_ss+(teta-teta_ss)*(exp(-a*(t_slack*0.001)));
teta_final=teta_cool;
   