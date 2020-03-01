%% this block fetches the framework: do not edit
urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
    'part_composition_setup.m');
rehash;
part_composition_setup('v5');
%% the system is constructed & simulated in the block below: feel free to edit

% define a BioSystem: this object will hold the parts and compositors
sys = BioSystem();

% define constants. assume that concentration units are in nM and time in s

sys.AddConstant(Const('k_PETase', 1000)); % rand
sys.AddConstant(Const('k_MHETase', 4.25*10^3));
sys.AddConstant(Const('k_TPA_meta', 1)); % rand
sys.AddConstant(Const('k_on_MHET', 10^7/(10^9*60))); % 10^7 M-1 min-1 = 10^7/(10^9*60)
sys.AddConstant(Const('k_off_MHET', 0.05/60)); % 0.05 min-1
sys.AddConstant(Const('k_MHET_tln', 1000)); % rand
sys.AddConstant(Const('k_tln', 100)); % rand
sys.AddConstant(Const('k_txn', 100)); % rand
sys.AddConstant(Const('k_port_PETase', 1000)); % rand
sys.AddConstant(Const('k_port_MHETase', 1000)); % rand
sys.AddConstant(Const('k_GFP_deg', 1)); % rand
sys.AddConstant(Const('k_mKate_deg', 1000)); % rand
sys.AddConstant(Const('k_PETase_outdeg', 0)); % rand
sys.AddConstant(Const('k_MHETase_outdeg', 0)); % rand
sys.AddConstant(Const('k_Cas9byBonA_deg', 10000)); % rand
sys.AddConstant(Const('k_Cas9byAonB_deg', 1)); % rand
sys.AddConstant(Const('corpCas9byBonA', 1)); % rand
sys.AddConstant(Const('corpCas9byAonB', 1)); % rand

% add compositors: these are for the concentrations of species in the system.
% we use the short-hand version to add the compositors right away,
% the second parameter is the initial concentration

dPETndt = sys.AddCompositor('PETn', 0);
dPETnmin1dt = sys.AddCompositor('PETn_1', 0);
dPETaseoutdt = sys.AddCompositor('PETaseout', 1);
dMHETaseoutdt = sys.AddCompositor('MHETaseout', 0);
dPETaseindt = sys.AddCompositor('PETasein', 0);
dMHETaseindt = sys.AddCompositor('MHETasein', 0);
dMHETdt = sys.AddCompositor('MHET', 0);
dTPAdt = sys.AddCompositor('TPA', 0);
dEGdt = sys.AddCompositor('EG', 0);
dGPCRdt = sys.AddCompositor('GPCR', 1000);
dMHETGPCRdt = sys.AddCompositor('MHETGPCR', 0);
dmRNAPETaseMHETaseGFPCas9byAonBdt = sys.AddCompositor('mRNAPETaseMHETaseGFPCas9byAonB', 0);
dmKatedt = sys.AddCompositor('mKate', 0);
dGFPdt = sys.AddCompositor('GFP', 0);
dCas9byAonBdt = sys.AddCompositor('Cas9byAonB', 0);
dCas9byBonAdt = sys.AddCompositor('Cas9byBonA', 0);

% TODO: define parts and add them to the system
% Remember to use the same parts as in the Groundwork section of this problem set
P1 = Part('All the process', [dPETndt dPETnmin1dt dPETaseoutdt dMHETaseoutdt dPETaseindt dMHETaseindt dMHETdt dTPAdt dEGdt dGPCRdt dMHETGPCRdt dmRNAPETaseMHETaseGFPCas9byAonBdt dmKatedt dGFPdt dCas9byAonBdt dCas9byBonAdt], ...
          [Rate('-k_PETase*PETn*PETaseout'), ... % dPETndt
           Rate('k_PETase*PETn*PETaseout'), ... % dPETnmin1dt 
           Rate('k_port_PETase*PETasein - k_PETase_outdeg*PETaseout'), ... % dPETaseoutdt 
           Rate('k_port_MHETase*MHETasein - k_MHETase_outdeg*MHETaseout'), ... % dMHETaseoutdt 
           Rate('k_txn*mRNAPETaseMHETaseGFPCas9byAonB - k_port_PETase*PETasein'), ... % dPETaseindt 
           Rate('k_txn*mRNAPETaseMHETaseGFPCas9byAonB - k_port_MHETase*MHETasein'), ... % dMHETaseindt 
           Rate('k_PETase*PETn*PETaseout - k_MHETase*MHET*MHETaseout - k_on_MHET*MHET*GPCR + k_off_MHET*MHETGPCR'), ... % dMHETdt 
           Rate('k_MHETase*MHET*MHETaseout - k_TPA_meta*TPA'), ... % dTPAdt 
           Rate('k_MHETase*MHET*MHETaseout'), ... % dEGdt 
           Rate('- k_on_MHET*MHET*GPCR + k_off_MHET*MHETGPCR'), ... % dGPCRdt 
           Rate('k_on_MHET*MHET*GPCR - k_off_MHET*MHETGPCR - k_MHET_tln*MHETGPCR'), ... % dMHETGPCRdt 
           Rate('k_MHET_tln*MHETGPCR - k_txn*mRNAPETaseMHETaseGFPCas9byAonB'), ... % dmRNAPETaseMHETaseGFPCas9byAonBdt 
           Rate('k_tln*k_txn/(1+Cas9byAonB^corpCas9byAonB) - k_mKate_deg*mKate'), ... % dmKatedt 
           Rate('k_tln*k_txn*mRNAPETaseMHETaseGFPCas9byAonB /(1+Cas9byBonA^corpCas9byBonA) - k_GFP_deg*GFP'), ... % dGFPdt 
           Rate('k_tln*k_txn*mRNAPETaseMHETaseGFPCas9byAonB /(1+Cas9byBonA^corpCas9byBonA) - k_Cas9byAonB_deg*Cas9byAonB'), ... % dCas9byAonBdt 
           Rate('k_tln*k_txn/(1+Cas9byAonB^corpCas9byAonB) - k_Cas9byBonA_deg*Cas9byBonA')]);  % dCas9byBonAdt
% add parts
sys.AddPart(P1);


% [T , Y] = sys.run([0 0.01])

[T, Y] = sys.run_pulses([...
    Pulse(0, 'PETn', 0), ...     % initial conditions
    Pulse(0.000000001, 'PETn', 1000), ... % spike in 10 nM of inducer @ t = 1000
    Pulse(1, '', 0), ...     % stop the simulation at time 2000
]);


% plot
figure();
plot(T, Y(:, sys.CompositorIndex('PETn')), ...
     T, Y(:, sys.CompositorIndex('MHET')), ...
     T, Y(:, sys.CompositorIndex('GFP')), ...
     T, Y(:, sys.CompositorIndex('mKate')))
     % T, Y(:, sys.CompositorIndex('PETn_1')), ...
     % T, Y(:, sys.CompositorIndex('MHETaseout')), ...
     % T, Y(:, sys.CompositorIndex('PETaseout')), ...
     % T, Y(:, sys.CompositorIndex('Cas9byBonA')), ...
     % T, Y(:, sys.CompositorIndex('Cas9byAonB')), ...
% ylim([0 120]);
% legend('PETn','PETn_1', 'MHET', 'MHETaseout', 'PETaseout', 'Cas9byBonA', 'Cas9byAonB', 'GFP', 'mKate', 'Location', 'best')
legend('PETn', 'MHET', 'GFP', 'mKate', 'Location', 'best')
xlabel('Time (s)');
ylabel('Concentration (nM)');



% plot all compositors, useful for debugging:
% figure();
% num_compositors = size(Y, 2);
% for i = 1:num_compositors
%     subplot(num_compositors, 1, i);
%     plot(T, Y(:, i));
%     xlabel('Time');
%     ylabel(sprintf('[%s]', sys.compositors(i).name), 'Rotation', 0);
% end

% 
% DNA2_values = [0.0625 0.125 0.25 0.5 1 2 4];
% fractional_decreases = [];
% for DNA2 = DNA2_values
%     sys.ChangeInitialValue('DNA2', DNA2);
%     % solve/simulate the system. change the amount of inducer at 3 times
%     [T, Y] = sys.run_pulses([...
%         Pulse(0, 'Ind', 0), ...     % initial conditions
%         Pulse(1000, 'Ind', 10), ... % spike in 10 nM of inducer @ t = 1000
%         Pulse(2000, '', 0), ...     % stop the simulation at time 2000
%     ]);
%     highest_P2 = max(Y(:, sys.CompositorIndex('P2')));
%     last_P2 = Y(end, sys.CompositorIndex('P2'));
%     fractional_decreases(end + 1) = (highest_P2 - last_P2) / highest_P2;
% end
% figure();
% plot(DNA2_values, fractional_decreases, 'o-');
% title('Effect of changing DNA2')
% xlabel('Amount of DNA2 (nM)')
% ylabel('Fractional decrease in P2 after induction')