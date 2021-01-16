function [lambda_final,replicate_core_1,replicate_core_2,replicate_core_3,main_core,teta_rep_1,teta_rep_2,teta_rep_3,vf1,vf2,vf3]=...
    replication_new(power_obj,lambda_obj,lambda,exe_cost,task,core,teta_core,teta_obj,...
    No_cores,No_vflevels,ETS,execution_matrix,communication_cost)

%% Initialization
available_cores_replication=zeros((No_cores-1),No_cores);
f_normal(1,1)=0;
f_normal(1,2)=0.5;
f_normal(1,3)=1;
v(1,1)=1.06;
f(1,1)=300;
v(1,2)=1.1;
f(1,2)=600;
v(1,3)=1.2;
f(1,3)=900;
%% Initialization
teta0=298;
fmin =0;
teta_amb=293;
ro_s_replication=inf (1,No_vflevels);
ro_h_replication=inf (No_cores,No_vflevels);
lambda_rep_1=inf (No_cores,No_vflevels);
lambda_replication=inf(No_cores,No_vflevels);
t_final_replication=zeros(No_cores,No_vflevels);
replicate_core_1=0;
replicate_core_2=0;
replicate_core_3=0;
main_core=0;
candidate_lambda=inf(No_cores-1,No_vflevels);
lambda_final=inf;
teta_replication=teta0*ones(No_cores-1,No_vflevels);
power_replication=zeros(No_cores-1,No_vflevels);
G=0.3;
C=0.03;
alpha=0.1;
beta=-20;
available_cores_replication(1,1)=2;
available_cores_replication(2,1)=3;
available_cores_replication(3,1)=4;
available_cores_replication(1,2)=1;
available_cores_replication(2,2)=3;
available_cores_replication(3,2)=4;
available_cores_replication(1,3)=1;
available_cores_replication(2,3)=2;
available_cores_replication(3,3)=4;
available_cores_replication(1,4)=1;
available_cores_replication(2,4)=2;
available_cores_replication(3,4)=3;

neighbor(1,1)=2;
neighbor(2,1)=3;
neighbor(3,1)=4;

neighbor(1,2)=1;
neighbor(2,2)=3;
neighbor(3,2)=4;

neighbor(1,3)=1;
neighbor(2,3)=2;
neighbor(3,3)=4;

neighbor(1,4)=1;
neighbor(2,4)=2;
neighbor(3,4)=3;

for j=1:No_cores-1
    for k=1:No_vflevels
        teta_replication(available_cores_replication(j,core),k)=teta_core(1,available_cores_replication(j,core));
    end
end
teta_rep_1=0;
teta_rep_2=0;
teta_rep_3=0;
vf1=0;
vf2=0;
vf3=0;
Gn=0.1;
P_dyn=zeros(1,No_vflevels);
C_eff=(1e-2);
for k=1:No_vflevels
    P_dyn(1,k)=C_eff*v(1,k)*v(1,k)*f(1,k);
end
a=(G+((3*Gn)-alpha))/C;
b_n=zeros(No_cores-1,No_vflevels);
teta_ss =zeros(No_cores-1,No_vflevels);

for m=1:No_cores-1
    for k=1:No_vflevels
        t_final_replication(available_cores_replication(m,core),k)=ETS(task,available_cores_replication(m,core))+execution_matrix(task,available_cores_replication(m,core),k)+communication_cost(task,available_cores_replication(m,core));
    end
end

for r=1:No_cores-1
    for k=1:No_vflevels
        b_n(available_cores_replication(r,core),k)=(P_dyn(1,k)+beta+(G*(teta_amb))+...
            (Gn*teta_core(1,neighbor(1,available_cores_replication(r,core))))+...
            (Gn*teta_core(1,neighbor(2,available_cores_replication(r,core))))+...
            (Gn*teta_core(1,neighbor(3,available_cores_replication(r,core)))))/C;    
        teta_ss (available_cores_replication(r,core),k)=b_n(available_cores_replication(r,core),k)/a; 
        teta_replication(available_cores_replication(r,core),k)=...
            teta_ss(available_cores_replication(r,core),k)+(teta_core(1,available_cores_replication(r,core))-...
            teta_ss(available_cores_replication(r,core),k))...
            *(exp(-a*((execution_matrix(task,available_cores_replication(r,core),k)+...
            communication_cost(task,available_cores_replication(r,core)))*0.001)));
    end
end
    
%% Calculating Lambda for 1 Replication and Main Core
for k=1:No_vflevels
    ro_s_replication(1,k)=(10^((3*(1-f_normal(1,k)))/(1-fmin)));
end
for j=1:No_cores-1
    for k=1:No_vflevels
        ro_h_replication(available_cores_replication(j,core),k)=...
            exp((-5570)*((1/teta_replication(available_cores_replication(j,core),k))-(1/teta0)));
    end
end
lambda_0_replication=1e-5;
        for m=1:No_cores-1
            for k=1:No_vflevels
                if(teta_replication(available_cores_replication(m,core),k)<teta_obj)
                    lambda_replication(available_cores_replication(m,core),k)=lambda_0_replication*ro_s_replication(1,k)*...
                        ro_h_replication(available_cores_replication(m,core),k);
                    %GSFR
                    lambda_replication(available_cores_replication(m,core),k)=lambda_replication(available_cores_replication(m,core),k)/...
                        ((execution_matrix(task,available_cores_replication(m,core),k))+...
                        communication_cost(task,available_cores_replication(m,core))); 
                    
                        %Total Lambda with 1 Replication
                        

% %                                 lambda_rep_1(available_cores_replication(m,core),k)=1/...
% %                                ((1/lambda)+(1/lambda_replication(available_cores_replication(m,core),k))-...
% %                                (1/(lambda+lambda_replication(available_cores_replication(m,core),k))));
                    lambda_rep_1(available_cores_replication(m,core),k)=...
                    (lambda*lambda_replication(available_cores_replication(m,core),k)* ...
                    ((execution_matrix(task,available_cores_replication(m,core),k))+...
                    communication_cost(task,available_cores_replication(m,core)))*exe_cost)/(execution_matrix...
                    (task,available_cores_replication(m,core),k)+...
                    communication_cost(task,available_cores_replication(m,core))+exe_cost);   
                else
                    % Cooling in Replications
                    %cc='cooling_replication';
                    t_slack_replication=0;
                    while((teta_replication(available_cores_replication(m,core),k)>teta_obj))
                        t_slack_replication=t_slack_replication+1;
                        teta_replication(available_cores_replication(m,core),k)=cooling_GSFR_Theta_Power_8cores...
                        (teta_replication(available_cores_replication(m,core)),P_dyn(1,k),...
                        teta_replication(neighbor(1,available_cores_replication(m,core))),...
                        teta_replication(neighbor(2,available_cores_replication(m,core))),...
                        teta_replication(neighbor(3,available_cores_replication(m,core))));                        
                    end
                    t_final_replication(available_cores_replication(m,core),k)=...
                        t_final_replication(available_cores_replication(m,core),k)+t_slack_replication;
                    power_replication(available_cores_replication(m,core),k)=...
                        (P_dyn(1,k)+(((alpha*teta_replication(available_cores_replication(m,core),k))+beta)/10))/5;
                    if((teta_replication(available_cores_replication(m,core),k)<teta_obj)&&...
                            (power_replication(available_cores_replication(m,core),k)<power_obj))
                        
                            ro_h_replication(available_cores_replication(m,core),k)=exp((-5570)*...
                                ((1/teta_replication(available_cores_replication(m,core),k))-(1/teta0)));
                            lambda_replication(available_cores_replication(m,core),k)=...
                                lambda_0_replication*ro_s_replication(1,k)*ro_h_replication(available_cores_replication(m,core),k);
                            %GSFR
                            lambda_replication(available_cores_replication(m,core),k)=lambda_replication(available_cores_replication(m,core),k)/...
                                ((execution_matrix(task,available_cores_replication(m,core),k))+...
                                communication_cost(task,available_cores_replication(m,core)));   
                            %Total Lambda with 1 Replication
                            lambda_rep_1(available_cores_replication(m,core),k)=...
                            (lambda*lambda_replication(available_cores_replication(m,core),k)* ...
                            ((execution_matrix(task,available_cores_replication(m,core),k))+...
                            communication_cost(task,available_cores_replication(m,core)))*exe_cost)/(execution_matrix...
                            (task,available_cores_replication(m,core),k)+...
                            communication_cost(task,available_cores_replication(m,core))+exe_cost);                         
                   end
                end
            end
        end
% End of Calculating Lambda for 1 Replication and the Main Core
    
%% Fisrt level of replication
        level=0;
        % Testing level one of Replication
    for m=1:No_cores-1
        for k=1:No_vflevels
            if((lambda_rep_1(available_cores_replication(m,core),k)<lambda_obj))
                candidate_lambda(available_cores_replication(m,core),k)=t_final_replication(available_cores_replication(m,core),k);
                level=1;
            end
        end
    end
        % Select the Best Replication Option Between all 
        [min_TF,min_TF_index]=min(candidate_lambda);  %core
        [~,Min_TF_index] = min(min_TF);   %VF level
        if(level==1)
            lambda_final=lambda_rep_1(min_TF_index(Min_TF_index),Min_TF_index);
            replicate_core_1=min_TF_index(Min_TF_index);
            main_core=core;
            teta_rep_1=teta_replication(min_TF_index(Min_TF_index),Min_TF_index);
            vf1=Min_TF_index;
        end
        
        %% Calculating Lambda for 2 Replications
if(level==0)    % 1 replication is not enough %
    t=1;
%     tt=1;
    lambda_rep_2=inf(No_cores-1,No_vflevels,No_cores-1,No_vflevels);
    lambda_2_rep_final=inf(No_cores-1,No_vflevels,No_cores-1,No_vflevels);
    candidate_2_lambda=zeros(1,1);
    candidate_2_rep1=zeros(1,1);
    candidate_2_rep2=zeros(1,1);
    candidate_2_rep1_vf=zeros(1,1);
    candidate_2_rep2_vf=zeros(1,1);
%     candidate_2_lambda_prime=zeros(1,1);
%     candidate_2_rep1_prime=zeros(1,1);
%     candidate_2_rep2_prime=zeros(1,1);
%     candidate_2_rep1_vf_prime=zeros(1,1);
%     candidate_2_rep2_vf_prime=zeros(1,1);
    for j=1:No_cores-1  %2nd Replica
        for k=1:No_vflevels  %2nd Replica vf
            for m=1:No_cores-1  %First Replica
                for n=1:No_vflevels     %Firts Replica vf
                    if(j==m)
                        tmp_56=0;
                    else
                        lambda_rep_2(available_cores_replication(j,core),k,available_cores_replication(m,core),n)=...
                        (lambda_replication(available_cores_replication(m,core),n)*lambda_replication(available_cores_replication(j,core),k)*...
                        ((execution_matrix(task,available_cores_replication(m,core),n))+communication_cost(task,available_cores_replication(m,core)))*...
                        ((execution_matrix(task,available_cores_replication(j,core),k))+communication_cost(task,available_cores_replication(j,core))))/...
                        (((execution_matrix(task,available_cores_replication(m,core),n))+communication_cost(task,available_cores_replication(m,core)))+...
                        ((execution_matrix(task,available_cores_replication(j,core),k))+communication_cost(task,available_cores_replication(j,core))));
                        
% % % %                         lambda_rep_2(available_cores_replication(j,core),k,available_cores_replication(m,core),n)=...
% % % %                         1/((1/lambda_replication(available_cores_replication(m,core),n))+...
% % % %                         (1/lambda_replication(available_cores_replication(j,core),k))...
% % % %                         -(1/((lambda_replication(available_cores_replication(j,core),k))+...
% % % %                         (lambda_replication(available_cores_replication(m,core),n)))));
                    
                        lambda_2_rep_final(available_cores_replication(j,core),k,available_cores_replication(m,core),n)=...
                        (lambda_rep_2(available_cores_replication(j,core),k,available_cores_replication(m,core),n)*lambda*...
                        exe_cost*(((execution_matrix(task,available_cores_replication(m,core),n))+communication_cost(task,available_cores_replication(m,core)))+...
                        ((execution_matrix(task,available_cores_replication(j,core),k))+communication_cost(task,available_cores_replication(j,core)))))/...
                        (exe_cost+(((execution_matrix(task,available_cores_replication(m,core),n))+communication_cost(task,available_cores_replication(m,core)))+...
                        ((execution_matrix(task,available_cores_replication(j,core),k))+communication_cost(task,available_cores_replication(j,core)))));
                    
                    
                    
                        if(lambda_2_rep_final(available_cores_replication(j,core),k,available_cores_replication(m,core),n)<lambda_obj)
                            candidate_2_lambda(t,1)=lambda_2_rep_final(available_cores_replication(j,core),k,available_cores_replication(m,core),n);
                            candidate_2_rep2(t,1)=available_cores_replication(j,core);
                            candidate_2_rep1(t,1)=available_cores_replication(m,core);
                            candidate_2_rep1_vf(t,1)=n;
                            candidate_2_rep2_vf(t,1)=k;
                            t=t+1;
                            level=2;
%                         else
%                             candidate_2_lambda_prime(tt,1)=lambda_2_rep_final(available_cores_replication(j,core),k,available_cores_replication(m,core),n);
%                             candidate_2_rep1_prime(tt,1)=available_cores_replication(j,core);
%                             candidate_2_rep2_prime(tt,1)=available_cores_replication(m,core);
%                             candidate_2_rep1_vf_prime(tt,1)=n;
%                             candidate_2_rep2_vf_prime(tt,1)=k;
%                             tt=tt+1;
                        end
                    end
                end
            end
        end
    end
end
%%%%Second level of replication%%%%%
if(level==2)
    [~,selected_index]=min(candidate_2_lambda);
    selected_index = selected_index(1);
    % Select the Best Replication Option Between all 
    lambda_final=candidate_2_lambda(selected_index,1);
    replicate_core_2=candidate_2_rep2(selected_index,1);
    replicate_core_1=candidate_2_rep1(selected_index,1);
    main_core=core;
    vf1=candidate_2_rep1_vf(selected_index,1);
    vf2=candidate_2_rep2_vf(selected_index,1);
    teta_rep_1=teta_replication(replicate_core_1,vf1);
    teta_rep_2=teta_replication(replicate_core_2,vf2);        
end

%%      Calculating Lambda for 3 Replications
if(level==0)    % 2 replication is not enough %
    t_1=1;    
    lambda_rep_3=inf(No_cores-1,No_vflevels,No_cores-1,No_vflevels,No_cores-1,No_vflevels);
    lambda_rep_3_2=inf(No_cores-1,No_vflevels,No_cores-1,No_vflevels,No_cores-1,No_vflevels);
    lambda_3_rep_final=inf(No_cores-1,No_vflevels,No_cores-1,No_vflevels,No_cores-1,No_vflevels);
    candidate_3_lambda=zeros(1,1);
    candidate_3_rep1=zeros(1,1);
    candidate_3_rep2=zeros(1,1);
    candidate_3_rep3=zeros(1,1);
    candidate_3_rep1_vf=zeros(1,1);
    candidate_3_rep2_vf=zeros(1,1);
    candidate_3_rep3_vf=zeros(1,1);
    for c3=1:No_cores-1 %3rd Replica
        for v3=1:No_vflevels    %3rd Replica vf
            for j2=1:No_cores-1  %2nd Replica
                for k2=1:No_vflevels  %2nd Replica vf
                    for m1=1:No_cores-1  %First Replica
                        for n1=1:No_vflevels     %Firts Replica vf
                            if((c3==j2)||(j2==m1)||(c3==m1))
                                tmp_56=0;
                            else
                                lambda_rep_3(available_cores_replication(c3,core),v3,available_cores_replication(j2,core),k2,...
                                    available_cores_replication(m1,core),n1)=...
                                    (lambda_replication(available_cores_replication(m1,core),n1)*...
                                    lambda_replication(available_cores_replication(j2,core),k2)*...
                                    (execution_matrix(task,available_cores_replication(m1,core),n1)+...
                                    communication_cost(task,available_cores_replication(m1,core)))*...
                                    (execution_matrix(task,available_cores_replication(j2,core),k2)+...
                                    communication_cost(task,available_cores_replication(j2,core))))/...
                                    ((execution_matrix(task,available_cores_replication(m1,core),n1)+...
                                    communication_cost(task,available_cores_replication(m1,core)))+...
                                    (execution_matrix(task,available_cores_replication(j2,core),k2)+...
                                    communication_cost(task,available_cores_replication(j2,core))));
                                time_exe=((execution_matrix(task,available_cores_replication(m1,core),n1)+...
                                    communication_cost(task,available_cores_replication(m1,core)))+...
                                    (execution_matrix(task,available_cores_replication(j2,core),k2)+...
                                    communication_cost(task,available_cores_replication(j2,core))));
                                
                                lambda_rep_3_2(available_cores_replication(c3,core),v3,available_cores_replication(j2,core),k2,...
                                    available_cores_replication(m1,core),n1)=...
                                    (lambda_rep_3(available_cores_replication(c3,core),v3,available_cores_replication(j2,core),k2,...
                                    available_cores_replication(m1,core),n1)*lambda_replication(available_cores_replication(c3,core),v3)*...
                                    time_exe*(execution_matrix(task,available_cores_replication(c3,core),v3)+communication_cost...
                                    (task,available_cores_replication(c3,core))))/(time_exe+...
                                    (execution_matrix(task,available_cores_replication(c3,core),v3)+communication_cost...
                                    (task,available_cores_replication(c3,core))));
                                                                
                                lambda_3_rep_final(available_cores_replication(c3,core),v3,available_cores_replication(j2,core),k2,...
                                    available_cores_replication(m1,core),n1)=...
                                (lambda_rep_3_2(available_cores_replication(c3,core),v3,available_cores_replication(j2,core),k2,...
                                    available_cores_replication(m1,core),n1)*lambda*exe_cost*(time_exe+...
                                    (execution_matrix(task,available_cores_replication(c3,core),v3)+communication_cost...
                                    (task,available_cores_replication(c3,core)))))/((time_exe+...
                                    (execution_matrix(task,available_cores_replication(c3,core),v3)+communication_cost...
                                    (task,available_cores_replication(c3,core))))+exe_cost);
                                
                                if(lambda_3_rep_final(available_cores_replication(c3,core),v3,available_cores_replication(j2,core),k2,...
                                    available_cores_replication(m1,core),n1)<lambda_obj)
                                    candidate_3_lambda(t_1,1)=lambda_3_rep_final(available_cores_replication(c3,core),v3,available_cores_replication(j2,core),k2,...
                                    available_cores_replication(m1,core),n1);
                                    candidate_3_rep1(t_1,1)=available_cores_replication(m1,core);
                                    candidate_3_rep2(t_1,1)=available_cores_replication(j2,core);
                                    candidate_3_rep3(t_1,1)=available_cores_replication(c3,core);
                                    candidate_3_rep1_vf(t_1,1)=n1;
                                    candidate_3_rep2_vf(t_1,1)=k2;
                                    candidate_3_rep3_vf(t_1,1)=v3;
                                    t_1=t_1+1;
                                    level=3;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if(level==3)
    [~,selected_index_3]=min(candidate_3_lambda);
    lambda_final=candidate_3_lambda(selected_index_3,1);
    selected_index_3=selected_index_3(1);
    replicate_core_1=candidate_3_rep1(selected_index_3,1);
    replicate_core_2=candidate_3_rep2(selected_index_3,1);
    replicate_core_3=candidate_3_rep3(selected_index_3,1);
    main_core=core;
    vf1=candidate_3_rep1_vf(selected_index_3,1);
    vf2=candidate_3_rep2_vf(selected_index_3,1);
    vf3=candidate_3_rep3_vf(selected_index_3,1);
    teta_rep_1=teta_replication(replicate_core_1,vf1);
    teta_rep_2=teta_replication(replicate_core_2,vf2); 
    teta_rep_3=teta_replication(replicate_core_3,vf3);
end

        