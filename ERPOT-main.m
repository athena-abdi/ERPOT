function [cmax_final] = ERPOT_NEW_BUS(teta_obj,lambda_obj,power_obj)

% ######### Preparing the task graph as an input file and load it ########
% file = '...\graph_40_4.xlsx';
%  [No_tasks,graph_matrix,local_deadlines,pred,succ,execution_matrix,is_leaf]= TGFF_new_4cores(file);
%  save('graph_40-4');
% load('graph_40-4');
%########## Definition of MPSOC Charatcteristics ###############%
teta_obj_tmp = teta_obj;
lambda_obj_tmp = lambda_obj;
power_obj_tmp = power_obj;
teta_obj = teta_obj_tmp;
lambda_obj = lambda_obj_tmp;
power_obj = power_obj_tmp;
No_cores=4;
mpsoc = zeros (No_tasks,10);   % Final result: [slot#,core#,V/F#,finish_slot#,start_slot#]
No_vflevels=3;
v = -1*ones(No_cores,No_vflevels);
f = -1*ones(No_cores,No_vflevels);
f_normal = -1*ones(No_cores,No_vflevels);

v(1,1)=1.06;
f(1,1)=300;
v(1,2)=1.1;
f(1,2)=600;
v(1,3)=1.2;
f(1,3)=900;

f_normal(1,1)=1/3;
f_normal(1,2)=2/3;
f_normal(1,3)=1;

%##########     List Scheduling Algorithm     ###############%
%% Initialization
fmin =0;
teta0=298;
teta_core (1,1)=teta0;
teta_core (1,2)=teta0;
teta_core (1,3)=teta0;
teta_core (1,4)=teta0;
teta = zeros (No_tasks,No_cores,No_vflevels);
power_total=zeros(No_tasks,No_cores,No_vflevels);
teta_replication = teta0*ones (No_tasks,No_cores,No_vflevels);
lambda = ones (No_tasks,No_cores,No_vflevels);
GSFR = ones (No_tasks,No_cores,No_vflevels);
ro_s=ones (1,No_vflevels);
ro_h=ones (No_tasks,No_cores,No_vflevels);
lambda_0=1e-5;
lambda_0_bus = 5e-4;
n=1;    % Index of current iteration
L_ready = zeros(No_tasks,No_tasks);
L_scheduled = zeros(1,No_tasks);   
schedule =zeros(No_tasks,2); % each task is scheduled or not (the first column:1/0 and the second:n)    
schedule_pressure = inf*ones(No_tasks,No_cores);
urgent_task = -1*ones(No_tasks,1);  %urgent task in each iteration
end_slot = zeros(1,No_tasks);
ETS = -1*ones (No_tasks,No_cores);
ETS_new = -1*ones (No_tasks,No_cores,No_vflevels);
LTE = zeros (1,No_tasks);
CPL = zeros (1,No_tasks);
wcet = zeros (1,No_tasks);
Alte=zeros(No_tasks,No_tasks); % temporal LTE
busy_cores =zeros(1e5,No_cores);           % 1 if a core is busy
no_pred=zeros(1,No_tasks);
tmp_CPL=zeros(No_tasks,No_tasks);   %temporal Critical path length
cmax=0;
t_final=zeros(No_tasks,No_cores,No_vflevels);
feasible_core_level=zeros(No_tasks,No_cores,No_vflevels);
t_communication=zeros(No_tasks,No_tasks,No_cores);
communication_cost = zeros(No_tasks,No_cores);
link_delay=1;
%local_deadlines=local_deadlines*600;    % the system deadline (if available)
temp_deadline=zeros(1,No_tasks);
main_core=zeros(No_tasks,No_cores,No_vflevels);
replicate_core_1=zeros(No_tasks,No_cores,No_vflevels);
replicate_core_2=zeros(No_tasks,No_cores,No_vflevels);
replicate_core_3=zeros(No_tasks,No_cores,No_vflevels);
exe_cost=zeros(No_tasks,No_cores,No_vflevels);
% Define the constant values required in power and temperature equations
% here ( such as beta, alpha, C, G , Ea, K , ....)
bus_active = zeros (1,No_tasks);
busy_bus = zeros(1,5000); 
teta_bus  = teta0*ones (1,No_tasks);
t1  = teta0*ones (1,No_tasks);
t2  = teta0*ones (1,No_tasks);
t3  = teta0*ones (1,No_tasks);
t4  = teta0*ones (1,No_tasks);

for k=1:No_vflevels
    P_dyn(1,k)=C_eff*v(1,k)*v(1,k)*f(1,k);
end
P_dyn_bus = 1e-2*1*1*100;

for k=1:No_vflevels
    ro_s(1,k)=(10^((3*(1-f_normal(1,k)))/(1-fmin)));
end

teta_ss = zeros(No_cores,No_vflevels);
neighbor(1,1)=2;
% neighbor(2,1)=3;
% neighbor(3,1)=4;

neighbor(1,2)=1;
% neighbor(2,2)=2;
% neighbor(3,2)=3;

% neighbor(1,3)=1;
% neighbor(2,3)=2;
neighbor(1,3)=4;

% neighbor(1,4)=1;
% neighbor(2,4)=2;
neighbor(1,4)=3;

a=(G+((2*Gn)-alpha))/C;
a_bus = (G_bus+(4*Gn_bus)-alpha_bus)/C_bus;
b=zeros(No_cores,No_vflevels);
b_over=zeros(No_cores,1);
tet_ss=zeros(1,No_vflevels);
teta_over_estimation=zeros(1,1);
VF_replication_1=zeros(No_tasks,No_cores,No_vflevels);
VF_replication_2=zeros(No_tasks,No_cores,No_vflevels);
VF_replication_3=zeros(No_tasks,No_cores,No_vflevels);
temperature_final = zeros(1,No_tasks);
power_final = zeros(1,No_tasks);
GSFR_final= zeros(1,No_tasks);
idle_core_1 = zeros (No_cores,No_tasks);
idle_time = zeros (No_cores,No_tasks);
teta_bus_last = teta0;

%% Constructing the Ready_List
No_tasks_ready=0;
for j=1:No_tasks
    all_pred_scheduled = 1;
    for k=1:No_tasks
        if(pred(j,k) ==1 && schedule(k,1)==0)
            all_pred_scheduled = 0;
            break;
        end
    end
    
    if (pred(j,:)==zeros(1,No_tasks))
        no_pred(1,j)=1;
    end
    
    if (no_pred(1,j)==1 || all_pred_scheduled == 1)   % if j has no pred or its pred is scheduled former   
          L_ready(n,j)=1;                         %1:ready tasks of each slot
          No_tasks_ready = No_tasks_ready+1;
    end
end


% checked
%WCET Computation
for k=1:No_tasks
    wcet(1,k) = execution_matrix(k,1,3);
end
% checked
for i=1:No_tasks
        if(local_deadlines(1,i)<inf)
            temp_deadline(1,i)=local_deadlines(1,i);            
        end
        %temp_deadline(1,i)=temp_deadline(1,i)*100;
end
[max_deadline_value,~]=max(temp_deadline(1,:));
% checked

%% List Scheduling Algoritm
while (sum(L_ready(n,:)) > 0)          % while ready list != 0
    %curr_step=n
    for i=1:No_tasks
        if (L_ready(n,i)==1)
                    %% LTE computation%  Longest Path From Here to the end %
                    if (is_leaf(1,i)==1)   %i is a leaf
                        Alte(i,i) =  wcet(1,i); 
                    else                        %i is not a leaf
                        for w=1:No_tasks
                            if ((succ(i,w) ==1) && is_leaf(1,w)==1) % w should be a leaf and each leaf has a local_deadline 
                                Alte(i,w)=wcet(1,w);
                                for z=1:No_tasks
                                    if(pred(w,z)==1)
                                        if (succ(i,z)==1)
                                            Alte(i,w)= Alte(i,w)+wcet(1,z);
                                        end
                                    end
                                end
                            end 
                            %if(A(i,w) ~= 0)
                            Alte(i,w)=Alte(i,w)+wcet(1,i);
                            %end
                        end
                     end
                     LTE (1,i) = max (Alte(i,:));
                     [~,urgent_task(n,1)] = max(LTE(1,:));
                     
        end
    end
                    %% Step 1 of while --Schedule Pressure Computation--%%
                    % ETS Construction
                    %%% considering end slot of preds in ETS computation
                    max_end_pred = -inf;
                    if (sum(pred(urgent_task(n,1),:))>0)
                        for k=1:No_tasks
                           if(pred(urgent_task(n,1),k)==1)
                                if(end_slot(1,k) > max_end_pred)
                                    max_end_pred = end_slot(1,k);
                                end
                           end
                        end
                        max_end_pred=max_end_pred+1;
                    else                             
                        max_end_pred = n;         
                    end
             
                    for j=1:No_cores
                       if (sum(busy_cores(max_end_pred,j))==0)
                            ETS(urgent_task(n,1),j)=max_end_pred;
                       else
                           ETS(urgent_task(n,1),j)=find(busy_cores(:,j),1,'last')+1;
                       end

                    end
    %checked
    
    
                    %% Step 2 of while   --best core for each ready task--

                    %%% Finding the cores and VF levels which meet our constraints %%%
                                        
                            for j=1:No_cores
                                 for k=1:No_vflevels
                                     if (sum(pred(urgent_task(n,1),:))>0)  %if urg has any pred
                                         for m=1:No_tasks
                                             if (graph_matrix(m,urgent_task(n,1))==1)
                                                 
                                                 if(mpsoc(m,2)==j)
                                                     t_communication (m,urgent_task(n,1),j)=0;
                                                 else
                                                     t_communication (m,urgent_task(n,1),j)=link_delay;
                                                 end
                                             end
                                         end
                                     end
                                     if (sum(t_communication(:,urgent_task(n,1),j))>0)
                                         communication_cost(urgent_task(n,1),j) = link_delay;
                                         bus_active (1,urgent_task(n,1))= 1;
                                     else 
                                         communication_cost(urgent_task(n,1),j) = 0;
                                     end
                                     
                                     if (bus_active (urgent_task(n,1)) == 1)
                                         if (sum (busy_bus(1,:))==0)
                                             free_bus = 1;
                                         else
                                             free_bus = find(busy_bus(1,:),1,'last');
                                         end
                                         
                                         for t=free_bus:free_bus+link_delay
                                             busy_bus (1,t) = 1;
                                         end
                                         %ETS_bus=find(busy_bus(1,:),1,'last')+1;
                                         
                                         ETS(urgent_task(n,1),j) = max(ETS(urgent_task(n,1),j) , free_bus);
                                     end
                                     %if(bus_active (1,urgent_task(n,1)) == 1)
                                         b_bus=(bus_active (1,urgent_task(n,1))*P_dyn_bus + beta_bus +(G_bus*(teta_amb))+(Gn_bus*teta_core(1,1))+(Gn_bus*teta_core(1,2))+(Gn_bus*teta_core(1,3))+(Gn_bus*teta_core(1,4)))/C_bus;
                                         teta_ss_bus = b_bus / a_bus; 
                                         teta_bus(1,urgent_task(n,1)) = teta_ss_bus + (teta_bus_last-teta_ss_bus)*(exp(-a_bus*link_delay*0.001));
                                         %teta(urgent_task(n,1),j,k) =  teta(urgent_task(n,1),j,k) + teta_bus(1,urgent_task(n,1));
                                     %end
                                     teta_bus_last = teta_bus(1,urgent_task(n,1));
                                     
                                     t_final(urgent_task(n,1),j,k)=ETS(urgent_task(n,1),j)+execution_matrix(urgent_task(n,1),j,k)+communication_cost(urgent_task(n,1),j);
                                     cmax_now_1=find(busy_cores(:,neighbor(1,j)),1,'last');
%                                      cmax_now_2=find(busy_cores(:,neighbor(2,j)),1,'last');
%                                      cmax_now_3=find(busy_cores(:,neighbor(3,j)),1,'last');
                                     %Overestimation
                                     if(t_final(urgent_task(n,1),j,k)>cmax_now_1)
                                         b_over(neighbor(1,j),3)=(P_dyn(1,3)+beta+(G*(teta_amb))+...
                                             (Gn*teta_core(1,neighbor(1,j)))+(Gn*teta_bus(1,urgent_task(n,1))))/C;
                                         tet_ss(neighbor(1,j),3)=b_over(neighbor(1,j),3)/a;
                                         teta_over_estimation(1,neighbor(1,j))=tet_ss(neighbor(1,j),3)+((teta_core(1,neighbor(1,j))-...
                                             tet_ss(neighbor(1,j),3))*(exp(-a*((t_final(urgent_task(n,1),j,k)-cmax_now_1)*0.001))));                                         
                                         teta_over_estimation(1,neighbor(1,j))=teta_over_estimation(1,neighbor(1,j));
                                     else
                                         teta_over_estimation(1,neighbor(1,j))=teta_core(1,neighbor(1,j));
                                     end
                                                                        
                                     
                                     
                                     b(j,k)=(P_dyn(1,k)+beta+(G*(teta_amb))+(Gn*teta_over_estimation(1,neighbor(1,j)))+(Gn*teta_bus(1,urgent_task(n,1))))/C;
                                     teta_ss (j,k)= b(j,k)/a; 
                                     teta(urgent_task(n,1),j,k)=teta_ss(j,k)+(teta_core(1,j)-teta_ss(j,k))*(exp(-a*((execution_matrix(urgent_task(n,1),j,k)+communication_cost(urgent_task(n,1),j))*0.001)));
                                     
                                     if (teta(urgent_task(n,1),j,k) > teta_obj)
                                         t_slack=0;
                                         while((teta(urgent_task(n,1),j,k) > teta_obj))%&&((t_final(i,j,k)+t_slack)<max_deadline_value))
                                            t_slack=t_slack+1;
                                            teta(urgent_task(n,1),j,k)=cooling_ERPOT_Bus(teta(urgent_task(n,1),j,k),P_dyn(1,k),teta_over_estimation(1,neighbor(1,j)),...
                                            teta_bus(1,urgent_task(n,1)));
                                         end                                     
                                     else
                                          t_slack=0;
                                     end
                                     
                                     t_final(urgent_task(n,1),j,k)=t_final(urgent_task(n,1),j,k)+t_slack;
                                     %checked
                                     power_total(urgent_task(n,1),j,k)=P_dyn(1,k)+((alpha*teta(urgent_task(n,1),j,k))+beta);
                                     
                                     %if(bus_active (1,urgent_task(n,1)) == 1)
                                         power_bus = (bus_active (1,urgent_task(n,1))*P_dyn_bus+((alpha_bus*teta_bus(1,urgent_task(n,1))+beta_bus)/10))/5;
                                         power_total(urgent_task(n,1),j,k) = power_total(urgent_task(n,1),j,k) + power_bus;
                                     %end
                                     
                                     %%  Lambda  %%
                                     ro_h(urgent_task(n,1),j,k)=(exp((-Ea_k).*((1/teta(urgent_task(n,1),j,k))-(1/teta0))));                                     
                                     lambda(urgent_task(n,1),j,k)=lambda_0*ro_s(1,k)*ro_h(urgent_task(n,1),j,k);
                                     if(bus_active (1,urgent_task(n,1)) == 1)
                                         lambda(urgent_task(n,1),j,k) =lambda(urgent_task(n,1),j,k)+lambda_0_bus;
                                     end
                                     GSFR(urgent_task(n,1),j,k)=lambda(urgent_task(n,1),j,k)/(execution_matrix(urgent_task(n,1),j,k)+communication_cost(urgent_task(n,1),j));
                                     exe_cost(urgent_task(n,1),j,k)=(execution_matrix(urgent_task(n,1),j,k)+communication_cost(urgent_task(n,1),j));
                                 end
                            end
                    
                            for j=1:No_cores
                                for k=1:No_vflevels
                                     if((teta(urgent_task(n,1),j,k)<teta_obj)&&(power_total(urgent_task(n,1),j,k)<power_obj))  % If teta and power constraints are met  
                                        if(GSFR(urgent_task(n,1),j,k)>lambda_obj)        %If lambda constraint is not met
                                            %&&(lambda_best<lambda_obj)
                                            % Adding Replications %
                                           [lamb_in,rc1,rc2,rc3,mc,tr1,tr2,tr3,vf1,vf2,vf3]=replication_new...
                                               (power_obj,lambda_obj,GSFR(urgent_task(n,1),j,k),exe_cost(urgent_task(n,1),j,k),urgent_task(n,1),j,teta_core,teta_obj,...
                                               No_cores,No_vflevels,ETS,execution_matrix,communication_cost);
                                           GSFR(urgent_task(n,1),j,k)=(lamb_in);
                                           replicate_core_1(urgent_task(n,1),j,k)=rc1;
                                           replicate_core_2(urgent_task(n,1),j,k)=rc2;
                                           replicate_core_3(urgent_task(n,1),j,k)=rc3;
                                           VF_replication_1(urgent_task(n,1),j,k)=vf1(1);
                                           VF_replication_2(urgent_task(n,1),j,k)=vf2(1);
                                           VF_replication_3(urgent_task(n,1),j,k)=vf3(1);
                                           main_core(urgent_task(n,1),j,k)=mc;
                                           teta_replication_1=tr1;
                                           teta_replication_2=tr2;
                                           teta_replication_3=tr3;
                                           if(replicate_core_1(urgent_task(n,1),j,k)>0)
                                             teta_replication(urgent_task(n,1),replicate_core_1(urgent_task(n,1),j,k),k)=teta_replication_1;
                                           end
                                           if(replicate_core_2(urgent_task(n,1),j,k)>0)
                                             teta_replication(urgent_task(n,1),replicate_core_2(urgent_task(n,1),j,k),k)=teta_replication_2;
                                           end
                                           if(replicate_core_3(urgent_task(n,1),j,k)>0)
                                             teta_replication(i,replicate_core_3(urgent_task(n,1),j,k),k)=teta_replication_3;
                                           end
                                        end                                         
                                     end
                                        
                                 end
                            end
                    %     end
                    %end
                
                    A_matrix=zeros(No_tasks,1);
                    %for i=1:No_tasks
                    %    if(L_ready(n,i)==1)
                            for j=1:No_cores
                                for k=1:No_vflevels
                                    ETS_new(urgent_task(n,1),j,k)=ETS(urgent_task(n,1),j);
                                    if((replicate_core_1(urgent_task(n,1),j,k)>0)&&(replicate_core_2(urgent_task(n,1),j,k)>0)&&(replicate_core_3(urgent_task(n,1),j,k)>0))    %3 replications
                                        A_matrix(urgent_task(n,1),1)=ETS(urgent_task(n,1),replicate_core_1(urgent_task(n,1),j,k));
                                        A_matrix(urgent_task(n,1),2)=ETS(urgent_task(n,1),replicate_core_2(urgent_task(n,1),j,k));
                                        A_matrix(urgent_task(n,1),3)=ETS(urgent_task(n,1),replicate_core_3(urgent_task(n,1),j,k));
                                        A_matrix(urgent_task(n,1),4)=ETS(urgent_task(n,1),main_core(urgent_task(n,1),j,k));
                                        A_max=max(A_matrix(urgent_task(n,1),:));
                                        ETS_new(urgent_task(n,1),replicate_core_1(urgent_task(n,1),j,k),k)=A_max;
                                        ETS_new(urgent_task(n,1),replicate_core_2(urgent_task(n,1),j,k),k)=A_max;
                                        ETS_new(urgent_task(n,1),replicate_core_3(urgent_task(n,1),j,k),k)=A_max;
                                        ETS_new(urgent_task(n,1),main_core(urgent_task(n,1),j,k))=A_max;                                        
                                    elseif((replicate_core_1(urgent_task(n,1),j,k)>0)&&(replicate_core_2(urgent_task(n,1),j,k)>0)&&(replicate_core_3(urgent_task(n,1),j,k)==0))   %2 replication
                                        A_matrix(urgent_task(n,1),1)=ETS(urgent_task(n,1),replicate_core_1(urgent_task(n,1),j,k));
                                        A_matrix(urgent_task(n,1),2)=ETS(urgent_task(n,1),replicate_core_2(urgent_task(n,1),j,k));
                                        A_matrix(urgent_task(n,1),4)=ETS(urgent_task(n,1),main_core(urgent_task(n,1),j,k));
                                        A_max=max(A_matrix(urgent_task(n,1),:));
                                        ETS_new(urgent_task(n,1),replicate_core_1(urgent_task(n,1),j,k),k)=A_max;
                                        ETS_new(urgent_task(n,1),replicate_core_2(urgent_task(n,1),j,k),k)=A_max;
                                        ETS_new(urgent_task(n,1),main_core(urgent_task(n,1),j,k),k)=A_max;
                                    elseif((replicate_core_1(urgent_task(n,1),j,k)>0)&&(replicate_core_2(urgent_task(n,1),j,k)==0)&&(replicate_core_3(urgent_task(n,1),j,k)==0))  %1 replications
                                        A_matrix(urgent_task(n,1),1)=ETS(urgent_task(n,1),replicate_core_1(urgent_task(n,1),j,k));
                                        A_matrix(urgent_task(n,1),4)=ETS(urgent_task(n,1),main_core(urgent_task(n,1),j,k));
                                        A_max=max(A_matrix(urgent_task(n,1),:));
                                        ETS_new(urgent_task(n,1),replicate_core_1(urgent_task(n,1),j,k),k)=A_max;
                                        ETS_new(urgent_task(n,1),main_core(urgent_task(n,1),j,k),k)=A_max;
                                    elseif((replicate_core_1(urgent_task(n,1),j,k)==0)&&(replicate_core_2(urgent_task(n,1),j,k)==0)&&(replicate_core_3(urgent_task(n,1),j,k)==0))    %No replicaton
                                        ETS_new(urgent_task(n,1),j,k)=ETS(urgent_task(n,1),j);
                                    end                                                                    
                                end
                            end
                    %    end
                    %end
                    
                            for j=1:No_cores
                                for k=1:No_vflevels
                                    delta_ETS=ETS_new(urgent_task(n,1),j,k)-ETS(urgent_task(n,1),j);
                                    ETS(urgent_task(n,1),j)= ETS_new(urgent_task(n,1),j,k);
                                end
                            end
                    %
                            for j=1:No_cores
                                for k=1:No_vflevels
                                    if((power_total(urgent_task(n,1),j,k)<power_obj)&&(teta(urgent_task(n,1),j,k)<teta_obj)&&(GSFR(urgent_task(n,1),j,k)<lambda_obj))
                                                                                                                     
                                            feasible_core_level(urgent_task(n,1),j,k)=1;           %feasible
                                    else                                            
                                            feasible_core_level(urgent_task(n,1),j,k)=0;           %infeasible
                                    end
                                end
                            end
                    
                            for j=1:No_cores
                                for k=1:No_vflevels
                                    %%%%%if(n>1)
                                        schedule_pressure(urgent_task(n,1),j,k)= t_final(urgent_task(n,1),j,k);
                                    
                                end
                            end
                            
                    
                    %%% Finding best core independent of VF level %%%
                    best_pressure = -inf*ones(1,No_tasks);
                    best_core_index= -inf*ones(1,No_tasks);
                    best_core = zeros(1,No_tasks);
                    best_vf= zeros(1,No_tasks);
                    ok=0;
                    
                    %for i=1:No_tasks
                    %    if (L_ready(n,i)==1)
                            best_pressure(1,urgent_task(n,1))=inf;
                            best_core_index(1,urgent_task(n,1)) =inf;
                            t_feasible = 0;
                            for j=1:No_cores
                                for k=1:No_vflevels
                                    if((feasible_core_level(urgent_task(n,1),j,k))==1)
                                        t_feasible = 1;
                                        if(schedule_pressure(urgent_task(n,1),j,k) <= best_pressure(1,urgent_task(n,1)))
                                            best_pressure(1,urgent_task(n,1)) = schedule_pressure(urgent_task(n,1),j,k);
                                            best_core(1,urgent_task(n,1))=j;

                                            
                                            best_vf (1,urgent_task(n,1))=k;
                                            ok=1;
                                        end
                                    end
                                end
                            end
                            if(t_feasible == 0)
                                best_pressure(1,urgent_task(n,1))=-inf;
                                ok=0;
                            end
                    %    end
                    %end

           

            %% Step 5 earlier            
            if (ok == 0)     %the objectives did not meet                
                cmax = -1;
                msg='infeasible'
                break
            elseif (ok ==1)
            %% Step 4 of while --schedule urgent operation--
            %[~,urgent_task(n,1)] = max(LTE(1,:));
            core_urgent = best_core(1,urgent_task(n,1));
            VF_urgent=best_vf(1,urgent_task(n,1));


            if ((replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent)>0)&&replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent)>0&&...
                (replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent)>0))    %if 3 replications
                for t=ETS_new(urgent_task(n,1),main_core(urgent_task(n,1),core_urgent,VF_urgent),VF_urgent)...                                             % Main core
                    :t_final(urgent_task(n,1),main_core(urgent_task(n,1),core_urgent,VF_urgent),VF_urgent)
                    busy_cores(t,main_core(urgent_task(n,1),core_urgent,VF_urgent))=1;
                end
                for t=ETS_new(urgent_task(n,1),replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_1(urgent_task(n,1),core_urgent,VF_urgent))...
                    :t_final(urgent_task(n,1),replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_1(urgent_task(n,1),core_urgent,VF_urgent))
                    busy_cores(t,replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent))=1;
                end
                for t=ETS_new(urgent_task(n,1),replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_2(urgent_task(n,1),core_urgent,VF_urgent))...
                    :t_final(urgent_task(n,1),replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_2(urgent_task(n,1),core_urgent,VF_urgent))
                    busy_cores(t,replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent))=1;
                end
                for t=ETS_new(urgent_task(n,1),replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_3(urgent_task(n,1),core_urgent,VF_urgent))...
                    :t_final(urgent_task(n,1),replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_3(urgent_task(n,1),core_urgent,VF_urgent))
                    busy_cores(t,replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent))=1;
                end
                
            elseif ((replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent)>0)&&replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent)>0&&...
                (replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent)==0))    %if 2 replications
                for t=ETS_new(urgent_task(n,1),main_core(urgent_task(n,1),core_urgent,VF_urgent),VF_urgent)...
                    :t_final(urgent_task(n,1),main_core(urgent_task(n,1),core_urgent,VF_urgent),VF_urgent)
                    busy_cores(t,main_core(urgent_task(n,1),core_urgent,VF_urgent))=1;
                end
                for t=ETS_new(urgent_task(n,1),replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_1(urgent_task(n,1),core_urgent,VF_urgent))...
                    :t_final(urgent_task(n,1),replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_1(urgent_task(n,1),core_urgent,VF_urgent))
                    busy_cores(t,replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent))=1;
                end
                for t=ETS_new(urgent_task(n,1),replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_2(urgent_task(n,1),core_urgent,VF_urgent))...
                    :t_final(urgent_task(n,1),replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_2(urgent_task(n,1),core_urgent,VF_urgent))
                    busy_cores(t,replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent))=1;
                end
            elseif((replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent)>0)&&replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent)==0&&...
                (replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent)==0))    %if 1 replication    
                for t=ETS_new(urgent_task(n,1),core_urgent,VF_urgent):+t_final(urgent_task(n,1),core_urgent,VF_urgent)
                    busy_cores(t,core_urgent) = 1;
                end
                for t=ETS_new(urgent_task(n,1),replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_1(urgent_task(n,1),core_urgent,VF_urgent))...
                    :t_final(urgent_task(n,1),replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_1(urgent_task(n,1),core_urgent,VF_urgent))
                    busy_cores(t,replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent))=1;
                end
            elseif((replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent)==0)&&replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent)==0&&...
                (replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent)==0))  %if No replication
                for t=ETS_new(urgent_task(n,1),core_urgent,VF_urgent):+t_final(urgent_task(n,1),core_urgent,VF_urgent)
                    busy_cores(t,core_urgent) = 1;
                end
            end
            
                                               
            if ((replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent)>0)&&(replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent)>0)&&...
                    (replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent))>0)  %if 3 replications
                teta_core(1,core_urgent)=teta(urgent_task(n,1),core_urgent,VF_urgent);
                teta_core(1,replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent))=teta_replication(urgent_task(n,1),replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent),VF_urgent);
                teta_core(1,replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent))=teta_replication(urgent_task(n,1),replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent),VF_urgent);
                teta_core(1,replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent))=teta_replication(urgent_task(n,1),replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent),VF_urgent);
                B_matrix(1,1)=t_final(urgent_task(n,1),core_urgent,VF_urgent);
                B_matrix(1,2)=t_final(urgent_task(n,1),replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_1(urgent_task(n,1),core_urgent,VF_urgent));
                B_matrix(1,3)=t_final(urgent_task(n,1),replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_2(urgent_task(n,1),core_urgent,VF_urgent));
                B_matrix(1,4)=t_final(urgent_task(n,1),replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_3(urgent_task(n,1),core_urgent,VF_urgent));
                Bmax=max(B_matrix(1,:));
                mpsoc(urgent_task(n,1),5)=Bmax;
                
            elseif ((replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent)>0)&&(replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent)>0)&&...
                    (replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent))==0) %if 2 replications
                teta_core(1,core_urgent)=teta(urgent_task(n,1),core_urgent,VF_urgent);
                teta_core(1,replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent))=teta_replication(urgent_task(n,1),replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent),VF_urgent);
                teta_core(1,replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent))=teta_replication(urgent_task(n,1),replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent),VF_urgent);
                B_matrix(1,1)=t_final(urgent_task(n,1),core_urgent,VF_urgent);
                B_matrix(1,2)=t_final(urgent_task(n,1),replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_1(urgent_task(n,1),core_urgent,VF_urgent));
                B_matrix(1,3)=t_final(urgent_task(n,1),replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_2(urgent_task(n,1),core_urgent,VF_urgent));
                Bmax=max(B_matrix(1,:));
                mpsoc(urgent_task(n,1),5)=Bmax;
                mm = 1;
                for idle_core =1:No_cores
                    if((core_urgent == idle_core) || (replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent) == idle_core) || ...
                            (replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent) == idle_core))
                        idle_core_1 (mm,n) = 0;
                    else
                        idle_core_1 (mm,n) = idle_core;
                        idle_time (idle_core,n) = execution_matrix(urgent_task(n,1),core_urgent,VF_urgent);
                    end
                end
                
            elseif ((replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent)>0)&&(replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent)==0)&&...
                    (replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent))==0)    %if 1 replication
                teta_core(1,core_urgent)=teta(urgent_task(n,1),core_urgent,VF_urgent);
                teta_core(1,replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent))=teta_replication(urgent_task(n,1),replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent),VF_urgent);
                B_matrix(1,1)=t_final(urgent_task(n,1),core_urgent,VF_urgent);
                B_matrix(1,2)=t_final(urgent_task(n,1),replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent),VF_replication_1(urgent_task(n,1),core_urgent,VF_urgent));
                Bmax=max(B_matrix(1,:));
                mpsoc(urgent_task(n,1),5)=Bmax;
                mm=1;
                for idle_core=1:No_cores
                    if((core_urgent == idle_core) || (replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent) == idle_core))
                    else
                        idle_core_1(mm,n)=idle_core;
                        idle_time (idle_core,n) = execution_matrix(urgent_task(n,1),core_urgent,VF_urgent);
                        mm = mm +1;
                    end
                end
                
            elseif ((replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent)==0)&&(replicate_core_2(urgent_task(n,1),core_urgent,VF_urgent)==0)&&...
                    (replicate_core_3(urgent_task(n,1),core_urgent,VF_urgent))==0)
                teta_core(1,core_urgent)=teta(urgent_task(n,1),core_urgent,VF_urgent);
                mpsoc(urgent_task(n,1),5)=t_final(urgent_task(n,1),core_urgent,VF_urgent);
                mm=1;
                for idle_core=1:No_cores
                    if((core_urgent == idle_core) || (replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent) == idle_core))
                    else
                        idle_core_1(mm,n)=idle_core;
                        idle_time (idle_core,n) = execution_matrix(urgent_task(n,1),core_urgent,VF_urgent);
                        mm = mm +1;
                    end
                end
            end
            
            time_cool=execution_matrix (urgent_task(n,1),core_urgent,VF_urgent);
            for co=1:No_cores
                for idle_core = 1:No_cores
                    if(idle_core_1(co,n) == idle_core)
                            teta_core (1,idle_core) =  cooling_idle_cores (teta_core (1,idle_core),P_dyn(1,3),time_cool);
                            if(teta_core (1,idle_core)<298)
                                teta_core (1,idle_core) = 298;
                            end
                    end 
                end
            end
            
            for tt=1:execution_matrix(urgent_task(n,1),core_urgent,VF_urgent)
                teta_bus(1,urgent_task(n,1)) = cooling_BUS(teta_bus(urgent_task(n,1)),P_dyn_bus,teta_core(1,1),teta_core(1,2),...
                    teta_core(1,1),teta_core(1,4));
            end
            if (teta_bus(1,urgent_task(n,1)) < teta0)
                teta_bus(1,urgent_task(n,1)) = teta0;
            end
           teta_bus_last = teta_bus(1,urgent_task(n,1));
           
           t_bus(1,urgent_task(n,1)) = teta_bus_last+5;
           t1(1,urgent_task(n,1)) = teta_core (1,1);
           t2(1,urgent_task(n,1))= teta_core (1,2);
           t3(1,urgent_task(n,1))= teta_core (1,3);
           t4(1,urgent_task(n,1))= teta_core (1,4);
            
            
            %mpsoc (urgent_task(n,1),1)=power_total(urgent_task(n,1),core_urgent,VF_urgent);             %slot of execution (start time)
            mpsoc (urgent_task(n,1),1)=n;
            mpsoc (urgent_task(n,1),2)=core_urgent; %core
            mpsoc (urgent_task(n,1),3)=VF_urgent; %v/f level
            mpsoc (urgent_task(n,1),4)=teta_bus_last+9;%replicate_core_1(urgent_task(n,1),core_urgent,VF_urgent);
            mpsoc(urgent_task(n,1),6)=teta_core (1,1);
            mpsoc(urgent_task(n,1),7)=teta_core (1,2);
            mpsoc(urgent_task(n,1),8)=teta_core (1,3);
            mpsoc(urgent_task(n,1),9)=teta_core (1,4);
            mpsoc(urgent_task(n,1),10)=(teta_core (1,1)+teta_core (1,2)+teta_core (1,3)+teta_core (1,4)+teta_bus_last+9)/5;
            temperature_final (1,n) = teta (urgent_task(n,1),core_urgent,VF_urgent);
            power_final (1,n) = power_total(urgent_task(n,1),core_urgent,VF_urgent);
            GSFR_final (1,n) = GSFR(urgent_task(n,1),core_urgent,VF_urgent);
            end_slot (1,urgent_task(n,1))=ceil(mpsoc (urgent_task(n,1),5));
            schedule (urgent_task(n,1),1)=1;    % 1: this task is scheduled
            schedule (urgent_task(n,1),2)=n;    % scheduled in the nth slot        
            
            %% Step 6 of while
            for x=1:No_tasks
                L_scheduled(1,urgent_task (n,1)) = 1;  % (L_scheduled) U {urgent task}
            end
            
            L_ready(n+1,:)=L_ready(n,:);
            temp_succ = succ (urgent_task(n,1),:);
            for x=1:No_tasks
                if(temp_succ(1,x) ==1)
                    ready_t = 1;
                    temp_pred=pred(x,:);
                    for y=1:No_tasks
                        if (temp_pred(1,y) == 1 && schedule(y,1)==0)
                            ready_t = 0;
                            break;
                        end
                    end
                    if(ready_t == 1)
                        L_ready (n+1,x)=1;
                    end
                end
            end
            L_ready (n+1,urgent_task (n,1))=0;

            %% Step 7 of while
            for i=1:No_tasks
                No_scheduled=sum(L_scheduled(1,:));
%                 if(L_scheduled(1,i)==1)
%                     tmp_CPL(i,n) = critical_path_delay(i,graph_matrix,wcet,No_scheduled);
%                 end
            end
            %CPL(1,n)=max (tmp_CPL(:,n));

            n=n+1;
            schedule_pressure=zeros(No_tasks,No_cores);
            LTE = zeros (1,No_tasks);
            end
end 
max_iterations=n-1;
%% Backtrack Part

if(cmax>-1)

    cmax_final=ceil(max(mpsoc(:,5)));
    
else
    cmax_final=inf;
end   


