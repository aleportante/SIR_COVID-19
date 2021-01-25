%{
SIR Model
JM Mahoney (c) 2020
%}

clear; close all; clc;

samps = 100;    % repetitions
dayz = 180;     % number of days
N = 2e2;        % number of agents

% Parameters
initially_infected_percentage = 0.05;
recoverytimemean = 12; 
recoverytimestd = 1.5; 
infectradius = 0.1;
maxinfectchance = 0.85;
f = @(x)(maxinfectchance/ infectradius)*x;

percentage_social_distancing=0.5;
social_distancing =true;
social_distance = 0.05;

quarantine = false;
central_hub = true;
isodays = 2;
percentagesymptoms = 0.5;
probability_central_hub = 0.02;

different_speeds = true;
percentage_fast=0.3;
percentage_slow=0.7;

timesteps_per_day = 4;

% initialize 
TI = zeros(samps,dayz);
TR = TI;
TS = TI;
R0s = TI;
T_max = 0;
for JJ = 1:samps    % repeat Monte Carlo
    
    for k=1:N
    agents(k).pos = rand(1,2)*4; % 1 by 2 array with uniformly distributed numbers between 0 and 1 %#ok<*SAGROW>
    agents(k).infect = 0; % By default the agent is not infected
    agents(k).infectday = 1;
    agents(k).symptoms = true;
    end
    
    i0 = randi(N,N*initially_infected_percentage,1);  % choose 5 initially-infected agents. Random int with max N

    for k = 1:length(i0)
        agents(i0(k)).infect = 1; % The Random agents are labeled as infected
    end

    i1 = randi(N, N-N*percentagesymptoms, 1);

    for g = 1:length(i1)
        agents(i1(g)).symptoms = false;
    end

    for T = 1:dayz 
        daily_infections=0;
        index_of_infected = find([agents.infect]==1); % Returns indice of infected agents
        if isempty(index_of_infected)      % stop when no more infected agents
            break  
        end
        
        for T2 = 1:timesteps_per_day % halfdays - time span
            
            for k=1:N %Iterates over each agent
                current_no_infected = find([agents.infect]==1);
                
                if central_hub == true && 3>= agents(k).pos(1) >= 2.5 && 3>= agents(k).pos(2) >= 2.5
                    agents(k).pos = rand(1,2)*4;
                end
                % random walk agents
                th = 2*pi*rand; % angle

                if different_speeds ==true && rand(1,1) < percentage_fast
                    r = 0.8*rand; 
                    agents(k).pos = [agents(k).pos] + r*[cos(th) sin(th)]; 
                end

                if different_speeds ==true && rand(1,1) < percentage_slow
                    r = 0.075*rand;
                    agents(k).pos = [agents(k).pos] + r*[cos(th) sin(th)];
                end

                if different_speeds ==false
                    r = 0.075*rand;
                    agents(k).pos = [agents(k).pos] + r*[cos(th) sin(th)];
                end
            % keep in box

                if agents(k).pos(1) > 4
                    agents(k).pos(1) = 5 - agents(k).pos(1); % If an agent has a position > 1 (outside box) the position is changed to 2 - the position
                elseif agents(k).pos(1) < 0
                    agents(k).pos(1) = abs(agents(k).pos(1)); % if pos<0 absolute value is taken
                end

            % Same for other coordinate
                if agents(k).pos(2) > 4
                    agents(k).pos(2) = 5 - agents(k).pos(2);
                elseif agents(k).pos(2) < 0
                    agents(k).pos(2) = abs(agents(k).pos(2));
                end



                if central_hub == true && rand(1,1) < probability_central_hub
                     agents(k).pos = 2.5 + rand(1,2)*0.5;
                end


                % social distancing
                if social_distancing ==true && rand(1,1) < percentage_social_distancing
                    for j = 1:N           % iterating over all the 'i+1' particle which is interacting with the other
                    d = distance(agents(k), agents(j));  % calculating the diameter of 2 particles and then calculating the distance of separation between them
                        if d < social_distance 
                            Mx = (agents(k).pos(1) + agents(j).pos(1))/2;
                            My = (agents(k).pos(2) + agents(j).pos(2))/2;
                            M = [Mx, My];
                            th = atan((agents(k).pos(2)-agents(j).pos(2)/(agents(k).pos(1)-agents(j).pos(1))));
                            agents(k).pos = [M] + social_distance*[cos(th) sin(th)];
                        end
                    end
                end


                if agents(k).infect == 1    % infected agents

                    if quarantine == true && T >= agents(k).infectday + isodays && agents(k).symptoms == true
                        agents(k).pos = 4.1+ rand(1,2)*0.9;
                    end
                % recover
                    if T - agents(k).infectday > ... %T is current day infect day, is the day they got infected
                        recoverytimemean + recoverytimestd*randn
                        agents(k).infect = 2;
                        agents(k).pos = rand(1,2)*4;
                    end
                end

                if agents(k).infect == 0    % susceptible agents
                    % infected by neighbors
                    for j = 1:length(current_no_infected) %ii = array with indices of infected agents
                    
                        if distance(agents(k), agents(current_no_infected(j))) < infectradius ... %Norm returns the distance between the two agents
                            && rand < f(distance(agents(k), agents(current_no_infected(j))))  % random value is below the infection chance
                            agents(k).infect = 1;
                            daily_infections = daily_infections + 1;
                            agents(k).infectday = T; % Setting infect day to the current day
                        end
                    end
                end
            end
            
           if T2 == 4
             R_0 = (daily_infections/ recoverytimemean); 
           end
           
        end
        
        totalsus(T) = sum([agents.infect]==0);
        totalinfect(T) = sum([agents.infect]==1);
        totalrecover(T) = sum([agents.infect]==2);
        R0(T) = (daily_infections/ recoverytimemean);
    end
    
    if T > T_max
       T_max = T;
    end
    
    TI(JJ,1:size(totalinfect, 2)) = totalinfect;
    TR(JJ,1:size(totalrecover, 2)) = totalrecover;
    TS(JJ,1:size(totalsus, 2)) = totalsus;
    R0s(JJ,1:size(R0, 2)) = R0;
    disp(JJ)
end

%% Plots
%Cutting arrays
TI = TI(:, 1:T_max);
TR = TR(:, 1:T_max);
TS = TS(:, 1:T_max);
R0s = R0s(:, 1:T_max);
disp(size(TI))
disp(size(TR))
disp(size(TS))
disp(T_max)

%Infected
figure(1); clf;

subplot(2,1,1);

plot(1:T_max,quantile(TI,0.5),'k-',...
    1:T_max,quantile(TI,0.975),'r-',...
    1:T_max,quantile(TI,0.025),'g-');
grid on;

xlabel(' ','fontsize',1);
ylabel('Population Infected','fontsize',14);
legend('50%','97.5%','2.5%','fontsize',10);
set(gca,'fontsize',11);
axis([1 T_max 0 N]);
% Plot Recovered
subplot(2,1,2);

plot(1:T_max,quantile(TR,0.5),'k-',...
    1:T_max,quantile(TR,0.975),'r-',...
    1:T_max,quantile(TR,0.025),'g-');
grid on;

xlabel('Days Since Patient 0','fontsize',14);
ylabel('Population Removed','fontsize',14);
% legend('50%','97.5%','2.5%','fontsize',14);
set(gca,'fontsize',11);
axis([1 T_max 0 N]);

% Convergence
% figure(2);
% clf;
% 
% semilogx(cumsum(TI(:,10))./(1:samps)')
% grid on;
% axis tight;
% xlabel('N','fontsize',14);
% ylabel('Mean Percent Infected on Day 20','fontsize',14);
% set(gca,'fontsize',14);

figure(2);
plot(1:T_max,quantile(R0s,0.5),'k-',...
    1:T_max,quantile(R0s,0.975),'r-',...
    1:T_max,quantile(R0s,0.025),'g-');
grid on;
xlabel('Days Since Patient 0','fontsize',14);
ylabel('R0','fontsize',14);
legend('50%','97.5%','2.5%','fontsize',14);
set(gca,'fontsize',11);
axis ([1 T_max 0 3])
axis tight;

return