function  genExpPathLossResults()
clc;close all; clear;
%format long e;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.P = 1;
params.simulation_area_side = [-100 100];    % simulation area side to side
params.space_realizations = 100;
params.time_slots = 10;
        params.beta = 1;
        params.alpha = 1.037;
%sim = 'CoverageProbability_TAU';
%sim = 'CoverageProbability_CD';
%sim = 'CoverageProbability_CD_Ref';
%sim = 'NetworkThroughput_CD';
%sim = 'NetworkThroughput_CD_zoomed';
sim = 'ASE_CD';

switch(sim)
    case 'ASE_CD'
        %params.la_s =  [1e-3 3e-3 5e-3 8e-3 1e-2 3e-2 5e-2 8e-2 1e-1] ;
        params.la_s = [1e-6 10e-6 100e-6 1000e-6 10000e-6]
        params.la_u =   500e-6 ;               % users density (users/m2)
%         params.alpha = 5e-3;
%         params.beta = 2;
        tau_dB = 10;
        params.tau = 10^(tau_dB/10);
        
        %ASE_math  = com_ASE_with_smallcell_density(params);
        %ASE_simul = gen_ASE_with_smallcell_density(params);
        ASE_ref = com_ASE_with_smallcell_density_Ref(params);
                
        f1 = semilogx(params.la_s ,ASE_simul,'ko' ,params.la_s,ASE_math,'k-',params.la_s,ASE_ref,'r--');
        %title('$$\alpha = 3$$','Interpreter','LaTex')
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 2]),{'Simulation' , 'Analysis'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Area spectral efficiency (nats/s/Hz/$km^2$)','Interpreter','LaTex');
        
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
    case 'NetworkThroughput_CD'
        %params.la_s =  [1e-3 3e-3 5e-3 8e-3 1e-2 3e-2 5e-2 8e-2 1e-1] ;
        params.la_s =  [1e-6:1e-6:1e-5 2e-5:1e-5:1e-4 2e-4:1e-4:1e-3 2e-3:1e-3:1e-2 2e-2:1e-2:1e-1] ;
        params.la_u =   200e-6 ;               % users density (users/m2)
        %params.alpha = 5e-2;
        
        tau_dB = 10;
        params.tau = 10^(tau_dB/10);
        
        Nt_math_200  = com_network_throughput_with_smallcell_density(params);
        Nt_simul_200 = gen_network_throughput_with_smallcell_density(params);
        
       
 
        params.la_u =   600e-6 ;               % users density (users/m2)
        Nt_math_600  = com_network_throughput_with_smallcell_density(params);
        Nt_simul_600 = gen_network_throughput_with_smallcell_density(params);
        
        params.la_u =   7000e-6 ;
        Nt_math_7000  = com_network_throughput_with_smallcell_density(params);
        
        params.la_u =   20000e-6 ;
        Nt_math_20000  = com_network_throughput_with_smallcell_density(params);
        params.la_u =   50000e-6 ;
        Nt_math_50000  = com_network_throughput_with_smallcell_density(params);
        
        Nt_ref = com_network_throughput_with_smallcell_density_Ref(params);
        
        f1 = semilogx(params.la_s ,Nt_simul_200,'ko', ...
                      params.la_s,Nt_math_200,'k-' ,...
                      params.la_s ,Nt_simul_600,'b^' ,...
                      params.la_s,Nt_math_600,'b:' ,...
                      params.la_s,Nt_math_7000,'r:>' ,...
                      params.la_s,Nt_math_20000,'r:<' ,...
                      params.la_s,Nt_math_50000,'r:s' ,...
                      params.la_s,Nt_ref,'r--');
        %title('$$\alpha = 3$$','Interpreter','LaTex')
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend({'Simulation \lambda_u = 200 users/km^2' , 'Analysis \lambda_u = 200 users/km^2' , ...
                'Simulation \lambda_u = 600 users/km^2' , 'Analysis \lambda_u = 600 users/km^2',...
                'Analysis \lambda_u = 7000 users/km^2',...
                'Analysis \lambda_u = 20000 users/km^2',...
                'Analysis \lambda_u = 500000 users/km^2',...
                'Ref [11]',...
                'Interpreter','LaTex'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Network throughput (nats/s/Hz/$km^2$)','Interpreter','LaTex');
        
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold'); 
        
        case 'NetworkThroughput_CD_zoomed'
        %params.la_s =  [1e-3 3e-3 5e-3 8e-3 1e-2 3e-2 5e-2 8e-2 1e-1] ;
        params.la_s =  [1e-4:1e-4:1e-3 1e-3:1e-3:1e-2 2e-2:1e-2:1e-1] ;
        params.la_u =   100e-6 ;               % users density (users/m2)
        %params.alpha = 5e-2;
        
        tau_dB = 10;
        params.tau = 10^(tau_dB/10);
        
        Nt_math_100  = com_network_throughput_with_smallcell_density(params);
        Nt_simul_100 = gen_network_throughput_with_smallcell_density(params);
        
       
 
        params.la_u =   300e-6 ;               % users density (users/m2)
        Nt_math_300  = com_network_throughput_with_smallcell_density(params);
        Nt_simul_300 = gen_network_throughput_with_smallcell_density(params);
        
        params.la_u =   600e-6 ;               % users density (users/m2)
        Nt_math_600  = com_network_throughput_with_smallcell_density(params);
        Nt_simul_600 = gen_network_throughput_with_smallcell_density(params);
        
        f1 = semilogx(params.la_s ,Nt_simul_100,'ko', ...
                      params.la_s,Nt_math_100,'k-' ,...
                      params.la_s ,Nt_simul_300,'bs' ,...
                      params.la_s,Nt_math_300,'b:' ,...
                      params.la_s ,Nt_simul_600,'r^' ,...
                      params.la_s,Nt_math_600,'r--' );
        %title('$$\alpha = 3$$','Interpreter','LaTex')
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend({'Simulation \lambda_u = 100 users/km^2' , 'Analysis \lambda_u = 100 users/km^2' , ...
                'Simulation \lambda_u = 300 users/km^2' , 'Analysis \lambda_u = 300 users/km^2',...
                'Simulation \lambda_u = 600 users/km^2' , 'Analysis \lambda_u = 600 users/km^2',...
                'Interpreter','LaTex'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Network throughput (nats/s/Hz/$km^2$)','Interpreter','LaTex');
        
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
    case 'CoverageProbability_CD_Ref'
        params.la_s =  [1e-3 3e-3 5e-3 8e-3 1e-2 3e-2 5e-2 8e-2 1e-1] ;
        params.la_u =   500e-6 ;               % users density (users/m2)
%         params.alpha = 5e-2;
%         params.beta = 2;
        tau_dB = 10;
        params.tau = 10^(tau_dB/10);
        
        Pc_math  = com_coverage_probability_with_smallcell_density(params);
        Pc_ref   = com_coverage_probability_with_smallcell_density_Ref(params);
        
        f1 = semilogx(params.la_s ,Pc_ref,'r--' ,params.la_s,Pc_math,'k-');
        %title('$$\alpha = 3$$','Interpreter','LaTex')
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 2]),{'Ref[]' , 'Analysis'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Coverage probability $\mathcal{P}_c (\tau)$','Interpreter','LaTex');
        
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');           
        
    case 'CoverageProbability_CD'
       % params.la_s =  [1e-3 3e-3 5e-3 8e-3 1e-2 3e-2 5e-2 8e-2 1e-1] ;
        params.la_s =  [1e-3:1e-3:1e-2 2e-2:1e-2:1e-1] ;
        params.la_u =   300e-6 ;               % users density (users/m2)
%         params.alpha = 5e-2;
%         params.beta = 2;
        tau_dB = 10;
        params.tau = 10^(tau_dB/10);
        
        Pc_math  = com_coverage_probability_with_smallcell_density(params);
        Pc_simul = gen_coverage_probability_with_smallcell_density(params);
        Pc_ref = com_coverage_probability_with_smallcell_density_Ref(params);
        
        params.la_u =   600e-6 ;               % users density (users/m2)
        
        Pc_math_  = com_coverage_probability_with_smallcell_density(params);
        Pc_simul_ = gen_coverage_probability_with_smallcell_density(params);
        Pc_ref_ = com_coverage_probability_with_smallcell_density_Ref(params);
        
        f1 = semilogx(params.la_s ,Pc_simul,'ko' ,params.la_s,Pc_math,'k-',params.la_s,Pc_ref,'r--', params.la_s ,Pc_simul_,'b^' ,params.la_s,Pc_math_,'b-',params.la_s,Pc_ref_,'r-.');
        %title('$$\alpha = 3$$','Interpreter','LaTex')
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 2 4 5 3]),{'Simulation \lambda_u = 300 users/km^2' , 'Analysis \lambda_u = 300 users/km^2' , 'Simulation \lambda_u = 600 users/km^2' , 'Analysis \lambda_u = 600 users/km^2','Ref []','Interpreter','LaTex'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('Small cell density (cells/$m^2$) ','Interpreter','LaTex');
        ylabel('Coverage probability $\mathcal{P}_c (\tau)$','Interpreter','LaTex');
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');   
        case 'CoverageProbability_TAU'
        params.la_s =  1e-2;
        params.la_u =   300e-6 ;               % users density (users/m2)
%         params.alpha = 5e-2;
%         params.beta = 2;
        tau_dB = -10:5:30;
        params.tau = 10.^(tau_dB/10);
        
        Pc_math  = com_coverage_probability_with_tau(params);
        Pc_ref  = com_coverage_probability_with_tau_Ref(params);
        Pc_simul = gen_coverage_probability_with_tau(params);
        
        params.la_u =  600e-6;
        
        Pc_math_  = com_coverage_probability_with_tau(params);
        Pc_ref_  = com_coverage_probability_with_tau_Ref(params);
        Pc_simul_ = gen_coverage_probability_with_tau(params);
        
        f1 = plot(tau_dB ,Pc_simul,'ko' ,tau_dB,Pc_math,'k-',tau_dB,Pc_ref,'r--', tau_dB ,Pc_simul_,'b^' ,tau_dB,Pc_math_,'b-',tau_dB,Pc_ref_,'r-.');
        %title('$$\alpha = 3$$','Interpreter','LaTex')
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 2 4 5 3]),{'Simulation \lambda_u = 300 users/km^2' , 'Analysis \lambda_u = 300 users/km^2' , 'Simulation \lambda_u = 600 users/km^2' , 'Analysis \lambda_u = 600 users/km^2','Ref []','Interpreter','LaTex'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('$\tau$(dB)  ','Interpreter','LaTex');
        ylabel('Coverage probability $\mathcal{P}_c (\tau)$','Interpreter','LaTex');
        
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');     
end
end




%% Coverage Probability with Small cells Density
function [Pc] = com_coverage_probability_with_smallcell_density(params)

points = numel(params.la_s);

Pc = zeros(points,1);


for p = 1:points
    k = params.la_s(p) / params.la_u;   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    pa = 1 - po;
    switch params.beta
        case 2
            Pc(p) = exp((- pi * pa * params.la_s(p) / params.alpha) * log(1 + params.tau)  );
    
    end
end

end
function [Pc] = com_coverage_probability_with_smallcell_density_Ref(params)

points = numel(params.la_s);

Pc = zeros(points,1);


for p = 1:points

    switch params.beta
        case 2
            Pc(p) = exp((- pi * params.la_s(p) / params.alpha) * log(1 + params.tau)  );
    
    end
end

end
function [Pc] = gen_coverage_probability_with_smallcell_density(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_s);
Pc_U = zeros(points, params.space_realizations , params.time_slots);
for p = 1:points
    fprintf('\n')
    disp(['Users Density: ' , num2str(params.la_u)]);
    disp(['Smallcell Density: ' , num2str(params.la_s(p))]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s(p) * simulation_area;
        mu_u = params.la_u * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        
        Ac = zeros(N_cells,1);  % association matrix
        S_u = zeros(N_users,1);
        I_u = zeros(N_users,1);
        SIR_u = zeros(N_users,1);
        Server_u = zeros(N_users,1);
        r_u = zeros(N_users,1);

        [r_u , Server_u] = min(r_su);
        Ac(Server_u) = 1;  
        
        for t = 1:params.time_slots
            Hui = exprnd(1,N_cells,N_users) ;
            R_u =  params.P *  Hui.* exp(- params.alpha.*r_su.^params.beta);

            
            for i = 1:N_users;

                    S_u(i) = R_u(Server_u(i),i);
                    I_u(i) = sum(Ac .* R_u(:,i)) - S_u(i);
            end
            
            SIR_u = S_u ./ I_u; 

             
            Pc_U(p,m,t) = sum(SIR_u > params.tau) / N_users;
            
        end
    end
    
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;

Pc = sum(sum(Pc_U,3),2) / normfact;
end

%% Coverage Probability with tau
function [Pc] = com_coverage_probability_with_tau(params)

points = numel(params.tau);

Pc = zeros(points,1);


for p = 1:points
    k = params.la_s / params.la_u;   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    pa = 1 - po;
    switch params.beta
        case 2
            Pc(p) = exp((- pi * pa * params.la_s / params.alpha) * log(1 + params.tau(p))  );
    
    end
end

end
function [Pc] = com_coverage_probability_with_tau_Ref(params)

points = numel(params.tau);

Pc = zeros(points,1);


for p = 1:points

    switch params.beta
        case 2
            Pc(p) = exp((- pi * params.la_s / params.alpha) * log(1 + params.tau(p))  );
    
    end
end

end
function [Pc] = gen_coverage_probability_with_tau(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.tau);
Pc_U = zeros(points, params.space_realizations , params.time_slots);
for p = 1:points
    fprintf('\n')
    disp(['Users Density: ' , num2str(params.la_u)]);
    disp(['Smallcell Density: ' , num2str(params.la_s)]);
    disp(['Tau: ' , num2str(params.tau(p))]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s * simulation_area;
        mu_u = params.la_u * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        
        Ac = zeros(N_cells,1);  % association matrix
        S_u = zeros(N_users,1);
        I_u = zeros(N_users,1);
        SIR_u = zeros(N_users,1);
        Server_u = zeros(N_users,1);
        r_u = zeros(N_users,1);

        [r_u , Server_u] = min(r_su);
        Ac(Server_u) = 1;  
        
        for t = 1:params.time_slots
            Hui = exprnd(1,N_cells,N_users) ;
            R_u =  params.P *  Hui.* exp(- params.alpha.*r_su.^params.beta);

            
            for i = 1:N_users;

                    S_u(i) = R_u(Server_u(i),i);
                    I_u(i) = sum(Ac .* R_u(:,i)) - S_u(i);
            end
            
            SIR_u = S_u ./ I_u; 

             
            Pc_U(p,m,t) = sum(SIR_u > params.tau(p)) / N_users;
            
        end
    end
    
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;

Pc = sum(sum(Pc_U,3),2) / normfact;
end
%% Network Throughput with Small cells Density
function [Nt] = com_network_throughput_with_smallcell_density(params)

points = numel(params.la_s);

Nt = zeros(points,1);
 simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;

for p = 1:points
    k = params.la_s(p) / params.la_u;   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    pa = 1 - po;
    switch params.beta
        case 2
            Nt(p) = pa .* params.la_s(p) * 1e6.* log(1 + params.tau) .* exp((- pi * pa * params.la_s(p) / params.alpha) * log(1 + params.tau)  );
    
    end
end

end
function [Nt] = com_network_throughput_with_smallcell_density_Ref(params)

points = numel(params.la_s);

Nt = zeros(points,1);
 simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;

for p = 1:points
    k = params.la_s(p) / params.la_u;   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    pa = 1 - po;
    switch params.beta
        case 2
            Nt(p) =  params.la_s(p) *1e6 .* log(1 + params.tau) .* exp((- pi  * params.la_s(p) / params.alpha) * log(1 + params.tau)  );
    
    end
end

end
function [Nt] = gen_network_throughput_with_smallcell_density(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_s);
Nt_U = zeros(points, params.space_realizations , params.time_slots);
for p = 1:points
    fprintf('\n')
    disp(['Users Density: ' , num2str(params.la_u)]);
    disp(['Smallcell Density: ' , num2str(params.la_s(p))]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s(p) * simulation_area;
        mu_u = params.la_u * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        A = zeros(N_cells,N_users);  % association matrix
        Ac = zeros(N_cells,1);  % association matrix
        S_u = zeros(N_users,1);
        I_u = zeros(N_users,1);
        SIR_u = zeros(N_users,1);
        Server_u = zeros(N_users,1);
        r_u = zeros(N_users,1);
        N_connected = zeros(N_users,1);

        [r_u , Server_u] = min(r_su);
        
        for i = 1:N_users
            A(Server_u(i),i) = 1;
        end
        
        Aconn = sum(A,2);  
        Ac = Aconn >=1;
        
        for t = 1:params.time_slots
            Hui = exprnd(1,N_cells,N_users) ;
            R_u =  params.P *  Hui.* exp(- params.alpha.*r_su.^params.beta);
                                            
            for i = 1:N_users;
                    S_u(i) = R_u(Server_u(i),i);
                    I_u(i) = sum(Ac .* R_u(:,i)) - S_u(i);
                    N_connected(i) = Aconn(Server_u(i));
            end
            
            SIR_u = S_u ./ I_u; 
            
            It = SIR_u >= params.tau;  % indicator for outage  0 = outage
            Nt = log(1 + params.tau) * It ./ N_connected ; 
            Nt_U(p,m,t) = sum(Nt) ;
            
        end
    end
    
end
normfact = params.space_realizations * params.time_slots * simulation_area*1e-6 ;% * simulation_area;
Nt = sum(sum(Nt_U,3),2) / normfact;
end
%% ASE with Small cells Density
function [ASE] = com_ASE_with_smallcell_density(params)

points = numel(params.la_s);

ASE = zeros(points,1);

for p = 1:points
    k = params.la_s(p) / params.la_u;   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    pa = 1 - po;
    switch params.beta
        case 2
            ASE(p) = params.alpha / pi  *  1e6 ;
    
    end
end

end
function [ASE] = com_ASE_with_smallcell_density_Ref(params)

points = numel(params.la_s);

ASE = zeros(points,1);

for p = 1:points

    switch params.beta
        case 2
            ASE(p) = params.alpha / pi *  1e6 ;
        case 1
            F = @(w,a,la) (exp((2 .* pi .* la ./ a^2) .* polylog(2,-w)) .* ...
                (1 - ((2.*pi.*sqrt(la) ./ a) .* log(1 + w) .* exp(pi .* la .* (log(1+w)).^2 ./ a.^2) .* qfunc((sqrt(2.*pi.*la) ./ a) .* log(1+w))))) .* (1 + w).^-1  ;
            ASE(p) = params.la_s(p) .* integral(@(w)F(w,params.alpha,params.la_s(p)),0,inf) .* 1e6;
    
    end
end

end
function [ASE] = gen_ASE_with_smallcell_density(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_s);
Rate_U = zeros(points, params.space_realizations , params.time_slots);
for p = 1:points
    fprintf('\n')
    disp(['Users Density: ' , num2str(params.la_u)]);
    disp(['Smallcell Density: ' , num2str(params.la_s(p))]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s(p) * simulation_area;
        mu_u = params.la_u * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        Ac = zeros(N_cells,1);  % association matrix
        S_u = zeros(N_users,1);
        I_u = zeros(N_users,1);
        SIR_u = zeros(N_users,1);
        Server_u = zeros(N_users,1);
        r_u = zeros(N_users,1);
        N_connected = zeros(N_users,1);
                
        [r_u , Server_u] = min(r_su);
        for i = 1:N_users
            A(Server_u(i),i) = 1;
        end
        Aconn = sum(A,2);  
        Ac = Aconn >=1;
        
        for t = 1:params.time_slots
            Hui = exprnd(1,N_cells,N_users) ;
            R_u =  params.P *  Hui.* exp(- params.alpha.*r_su.^params.beta);

            
            for i = 1:N_users;

                    S_u(i) = R_u(Server_u(i),i);
                    I_u(i) = sum(Ac .* R_u(:,i)) - S_u(i);
                    N_connected(i) = Aconn(Server_u(i));
            end
            
            SIR_u = S_u ./ I_u; 
            
            Rate_u = log(1 + SIR_u) ./ N_connected ;
            Rate_u(isinf(Rate_u)) = 0;
            Rate_u(Rate_u == 0) = max(Rate_u(:));
            Rate_U(p,m,t) = sum(Rate_u);
            
        end
    end
    
end
normfact = params.space_realizations * params.time_slots * simulation_area*1e-6 ;% * simulation_area;
Rate_U(isinf(Rate_U)) = 0;
Rate_U(Rate_U == 0) = max(Rate_U(:));
ASE = sum(sum(Rate_U,3),2) / normfact;
end
%%
function [servers,r_servers] = getServers(user,distances, N)
servers = zeros(N,1);
r_servers = zeros(N,1);
for i = 1:N
    [r_servers(i) , servers(i)] = min(distances(:,user));
    distances(servers(i),user) = inf;
end
end
