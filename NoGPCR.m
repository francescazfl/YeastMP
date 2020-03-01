%% this block fetches the framework: do not edit
urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
    'part_composition_setup.m');
rehash;
part_composition_setup('v5');

%% the system is constructed & simulated in the block below: feel free to edit
% dCas9 Cooperativity = 1, 1.055 and 1.6
% define a BioSystem: this object will hold the parts and compositors
 % plot
figure();

dCas9Corp = [1,1.055,1.6];
spiket = [1*10^4,5*10^5,1.75*10^3];
stopt = [2.5*10^4,10*10^5,3.5*10^3];
for i = 1:3

    sys = BioSystem();

    % define constants. assume that concentration units are in nM and time in s

    sys.AddConstant(Const('k_B', 5)); % 50
    sys.AddConstant(Const('k_A', 5)); % 50
    sys.AddConstant(Const('k_B_deg', 0.01)); % 0.1
    sys.AddConstant(Const('k_A_deg', 0.01)); % 0.1 
    sys.AddConstant(Const('corpB', dCas9Corp(i))); % 
    sys.AddConstant(Const('corpA', dCas9Corp(i))); % 2
    sys.AddConstant(Const('K', 1)); % 2

    % add compositors: these are for the concentrations of species in the system.
    % we use the short-hand version to add the compositors right away,
    % the second parameter is the initial concentration

    dAdt = sys.AddCompositor('A', 0.1);
    dBdt = sys.AddCompositor('B', 0);

    % Define parts and add them to the system
    P1 = Part('AB', [dAdt dBdt], ...
              [Rate('k_A/(1+(B/K)^corpB) - k_A_deg*A'), ... % dAdt 
               Rate('k_B/(1+(A/K)^corpA) - k_B_deg*B')]);  % dBdt
    % add parts
    sys.AddPart(P1);
    
    [T, Y] = sys.run_pulses([Pulse(0, 'B', 0), ...     % initial conditions
                             Pulse(spiket(i), 'A', 0), ... % spike in 1@ t = 
                             Pulse(stopt(i), '', 0), ...     % stop the simulation at time 
                             ])
  
    subplot(2,2,i)
    plot(T, Y(:, sys.CompositorIndex('A')),...
         T, Y(:, sys.CompositorIndex('B')));
    legend('GFP','mKate','Location', 'best');
    xlabel('Time (s)');
    ylabel('Concentration (nM)');
    title(['dCas9 Cooperativity = ',num2str(dCas9Corp(i))]);
end

sgtitle('dCas9 Needs Cooperativity (>1), starting {\geq} 1.055 could work');