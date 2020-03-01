%% this block fetches the framework: do not edit
urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
    'part_composition_setup.m');
rehash;
part_composition_setup('v5');
%% the system is constructed & simulated in the block below: feel free to edit
figure();
PETnConc = [0,0.01,0.1,1,10,100]; %in um
casenumb = length(PETnConc);

for i = 1:casenumb

% define a BioSystem: this object will hold the parts and compositors
sys = BioSystem();
% define constants. assume that concentration units are in nM and time in s

sys.AddConstant(Const('k_PETase', 10^-4)); % rand
sys.AddConstant(Const('k_MHETase', 4.25*10^-3)); % kcat/kM  = (31 ± 0.8 s-1) / (7.3 ± 0.6 ?M) = 4.25 ?M-1 ? s-1
sys.AddConstant(Const('k_TPA_meta', 10^-3)); % rand
sys.AddConstant(Const('k_on_MHET', 10^7/(10^9*60))); % 10^7 M-1 min-1 = 10^7/(10^9*60)
sys.AddConstant(Const('k_off_MHET', 0.05/60)); % 0.05 min-1
sys.AddConstant(Const('k_MHET_tln',2.24)); % assume same as normal
sys.AddConstant(Const('k_tln', 2.24)); % 0.1 nM/h http://homepages.ulb.ac.be/~dgonze/BIONUMBERS/bionumbers.html?fbclid=IwAR0js_SN9QkwW5e0r4ZeJnF96TTh_LqhtKdFVsmaN2l5hMuDOUVHqQM4W0I#transcriptionrate
sys.AddConstant(Const('k_txn', 0.224)); % 0.01 nM/h
sys.AddConstant(Const('k_port_PETase', 3.4/60)); % 3.4 min-1 https://onlinelibrary.wiley.com/doi/epdf/10.1002/bit.260330305
sys.AddConstant(Const('k_port_MHETase', 3.4/60)); % 3.4 min-1 https://onlinelibrary.wiley.com/doi/epdf/10.1002/bit.260330305
sys.AddConstant(Const('k_GFP_deg', 10^-3)); % 5 hr half life, k = ln2/5hr = ln2/(5*3600), 0.1x
sys.AddConstant(Const('k_mKate_deg', 0.01)); % same as above
sys.AddConstant(Const('k_PETase_outdeg', 0.001)); % 0.1x
sys.AddConstant(Const('k_MHETase_outdeg', 0.001)); %
sys.AddConstant(Const('k_Cas9byBonA_deg', 0.1)); % 10x faster 
sys.AddConstant(Const('k_Cas9byAonB_deg', 0.001)); % 0.1x 
sys.AddConstant(Const('corpCas9byBonA', 1.6)); % rand
sys.AddConstant(Const('corpCas9byAonB', 1.6)); % rand
% sys.AddConstant(Const('cells', 3*10^7*10^3)); % Say yeast OD600 1.0  ~3 x 10^7 cells/ml.
% Say the device 1L
% So everything below is for pet cell, now we have 3*10^7*10^3 cells

% add compositors: these are for the concentrations of species in the system.
% we use the short-hand version to add the compositors right away,
% the second parameter is the initial concentration

    dPETndt = sys.AddCompositor('PETn', PETnConc(i)*10^3); % need in nm
    dPETnmin1dt = sys.AddCompositor('PETn_1', 0);
    dPETaseoutdt = sys.AddCompositor('PETaseout', 10);
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

    % Define parts and add them to the system
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
               Rate('-k_on_MHET*MHET*GPCR + k_off_MHET*MHETGPCR'), ... % dGPCRdt 
               Rate('k_on_MHET*MHET*GPCR - k_off_MHET*MHETGPCR - k_MHET_tln*MHETGPCR'), ... % dMHETGPCRdt 
               Rate('k_MHET_tln*MHETGPCR - k_txn*mRNAPETaseMHETaseGFPCas9byAonB'), ... % dmRNAPETaseMHETaseGFPCas9byAonBdt 
               Rate('k_tln*k_txn/(1+(Cas9byAonB)^corpCas9byAonB) - k_mKate_deg*mKate'), ... % dmKatedt 
               Rate('k_tln*k_txn*mRNAPETaseMHETaseGFPCas9byAonB /(1+(Cas9byBonA)^corpCas9byBonA) - k_GFP_deg*GFP'), ... % dGFPdt 
               Rate('k_tln*k_txn*mRNAPETaseMHETaseGFPCas9byAonB /(1+(Cas9byBonA)^corpCas9byBonA) - k_Cas9byAonB_deg*Cas9byAonB'), ... % dCas9byAonBdt 
               Rate('k_tln*k_txn/(1+(Cas9byAonB)^corpCas9byAonB) - k_Cas9byBonA_deg*Cas9byBonA')]);  % dCas9byBonAdt
    % add parts
    sys.AddPart(P1);

    [T, Y] = sys.run_pulses([...
        Pulse(0, 'PETn', PETnConc(i)*10^3), ... % initial conditions, need in nm
        % Pulse(0.000000001, 'PETn', 100), ... % spike in 1@ t = 
        Pulse(10000, '', 0), ...     % stop the simulation at time 
    ]);

    subplot(2,casenumb,i)
    plot(T, Y(:, sys.CompositorIndex('PETn')))
    legend('PETn','Location', 'best')
    xlabel('Time (s)');
    ylabel('PETn Concentration (uM)');
    title(['With',num2str(PETnConc(i)),'uM PETn']);

    subplot(2,casenumb,casenumb+i)
    plot(T, Y(:, sys.CompositorIndex('GFP')), ...
         T, Y(:, sys.CompositorIndex('mKate')))
         % T, Y(:, sys.CompositorIndex('Cas9byAonB')),...
         % T, Y(:, sys.CompositorIndex('Cas9byBonA')),...
         % T, Y(:, sys.CompositorIndex('MHET'))) ...
         % T, Y(:, sys.CompositorIndex('PETn_1')), ...
         % T, Y(:, sys.CompositorIndex('MHETaseout')), ...
         % T, Y(:, sys.CompositorIndex('PETaseout')), ...
         % , ...
         % T, Y(:, sys.CompositorIndex('Cas9byAonB')), ...
    % ylim([0 120]);
    % legend('PETn','PETn_1', 'MHET', 'MHETaseout', 'PETaseout', 'Cas9byBonA', 'Cas9byAonB', 'GFP', 'mKate', 'Location', 'best')
    % legend('PETn', 'MHET', 'GFP', 'mKate', 'Location', 'best')
    % 'Cas9-|mKate','Cas9-|GFP',
    %legend('Cas9byAonB','Cas9byBonA','GFP', 'mKate', 'Location', 'best')
    legend('GFP','mKate','Location', 'best')
    xlabel('Time (s)');
    ylabel('Protein Concentration (nM)');
    title(['With',num2str(PETnConc(i)),'uM PETn']);
end

sgtitle('Examples of PETn (1st Row) and Protein Concentrations (2nd Row) Consumption Under Different PETn Starting Concentrations');