%% this block fetches the framework: do not edit
urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
    'part_composition_setup.m');
rehash;
part_composition_setup('v5');

%% 
PETnConc=[0,0.001,0.01,0.1,1,10,100]; %uM
tstop = 15000;

subplot(2,2,1) % no toggel different PETn conc.
for i = 1:length(PETnConc)
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
    sys.AddConstant(Const('k_PETase_outdeg', 0.001)); % 0.1x
    sys.AddConstant(Const('k_MHETase_outdeg', 0.001)); %

    % sys.AddConstant(Const('cells', 3*10^7*10^3)); % Say yeast OD600 1.0  ~3 x 10^7 cells/ml.
    % Say the device 1L
    % So everything below is for pet cell, now we have 3*10^7*10^3 cells

    % add compositors: these are for the concentrations of species in the system.
    % we use the short-hand version to add the compositors right away,
    % the second parameter is the initial concentration

    dPETndt = sys.AddCompositor('PETn', PETnConc(i)*10^3); %um
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
    dGFPdt = sys.AddCompositor('GFP', 0);

    % Define parts and add them to the system
    P1 = Part('All the process', [dPETndt dPETnmin1dt dPETaseoutdt dMHETaseoutdt dPETaseindt dMHETaseindt dMHETdt dTPAdt dEGdt dGPCRdt dMHETGPCRdt dmRNAPETaseMHETaseGFPCas9byAonBdt dGFPdt], ...
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
               Rate('k_tln*k_txn*mRNAPETaseMHETaseGFPCas9byAonB - k_GFP_deg*GFP')]);  % dGFP
    % add parts
    sys.AddPart(P1);

    [T, Y] = sys.run_pulses([...
        Pulse(0, 'PETn', PETnConc(i)*10^3), ...     % initial conditions, in nm
        % Pulse(0.000000001, 'PETn', 100), ... % spike in 1@ t = 
        Pulse(tstop, '', 0), ...     % stop the simulation at time 
    ]);

    T = T./60;

    % plot
    hold on
    % plot(T, Y(:, sys.CompositorIndex('PETn')))
    plot(T, Y(:, sys.CompositorIndex('GFP')))
end

legend('0 nM','1 nM','10 nM','100 nM','1 uM','10 uM','100 uM', 'Location', 'best')
xlim([0, 150]);
xlabel('Time (mins)');
ylabel('GFP Concentration (nM)');
title('No Toggle with Different PETn Concentrations No Distrubance');

%%
subplot(2,2,2) % with toggel different PETn conc.
storelen = zeros(1,length(PETnConc));
storeT = zeros(tstop,length(PETnConc));
storeGFPval = zeros(tstop,length(PETnConc));

for j = 1:length(PETnConc)
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

    dPETndt = sys.AddCompositor('PETn', PETnConc(j)*10^3); %in nm
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
        Pulse(0, 'PETn', PETnConc(j)*10^3), ...     % initial conditions, in nm
        % Pulse(0.000000001, 'PETn', 100), ... % spike in 1@ t = 
        Pulse(tstop, '', 0), ...     % stop the simulation at time 
    ]);

    T = T./60;

    hold on
    % plot(T, Y(:, sys.CompositorIndex('PETn')))
    plot(T, Y(:, sys.CompositorIndex('GFP')))
end

legend('0 nM','1 nM','10 nM','100 nM','1 uM','10 uM','100 uM', 'Location', 'best');
xlim([0,150]);
xlabel('Time (mins)');
ylabel('GFP Concentration (nM)');
title('Toggle System with Different PETn Concentrations No Distrubance');

%% 
subplot(2,2,3) % 'GFP with Toggle','mKate  with Toggle','GFP at 0.8 overall rate with Toggle','mKate at 0.8 overall rate with Toggle', 'GFP without Toggle', 'GFP at 0.8 overall rate without Toggle'
ratex=[1,0.8,0.6];
my_output = zeros(length(ratex)*2,8);
for k = 1:length(ratex)
    % define a BioSystem: this object will hold the parts and compositors
    sys = BioSystem();

    % define constants. assume that concentration units are in nM and time in s
    sys.AddConstant(Const('k_PETase', ratex(k)*10^-4)); % rand
    sys.AddConstant(Const('k_MHETase', ratex(k)*4.25*10^-3)); % kcat/kM  = (31 ± 0.8 s-1) / (7.3 ± 0.6 ?M) = 4.25 ?M-1 ? s-1
    sys.AddConstant(Const('k_TPA_meta', ratex(k)*10^-3)); % rand
    sys.AddConstant(Const('k_on_MHET', ratex(k)*10^7/(10^9*60))); % 10^7 M-1 min-1 = 10^7/(10^9*60)
    sys.AddConstant(Const('k_off_MHET', ratex(k)*0.05/60)); % 0.05 min-1
    sys.AddConstant(Const('k_MHET_tln',ratex(k)*2.24)); % assume same as normal
    sys.AddConstant(Const('k_tln', ratex(k)*2.24)); % 0.1 nM/h http://homepages.ulb.ac.be/~dgonze/BIONUMBERS/bionumbers.html?fbclid=IwAR0js_SN9QkwW5e0r4ZeJnF96TTh_LqhtKdFVsmaN2l5hMuDOUVHqQM4W0I#transcriptionrate
    sys.AddConstant(Const('k_txn', ratex(k)*0.224)); % 0.01 nM/h
    sys.AddConstant(Const('k_port_PETase', ratex(k)*3.4/60)); % 3.4 min-1 https://onlinelibrary.wiley.com/doi/epdf/10.1002/bit.260330305
    sys.AddConstant(Const('k_port_MHETase', ratex(k)*3.4/60)); % 3.4 min-1 https://onlinelibrary.wiley.com/doi/epdf/10.1002/bit.260330305
    sys.AddConstant(Const('k_GFP_deg', ratex(k)*10^-3)); % 5 hr half life, k = ln2/5hr = ln2/(5*3600), 0.1x
    sys.AddConstant(Const('k_mKate_deg', ratex(k)*0.01)); % same as above
    sys.AddConstant(Const('k_PETase_outdeg', ratex(k)*0.001)); % 0.1x
    sys.AddConstant(Const('k_MHETase_outdeg', ratex(k)*0.001)); %
    sys.AddConstant(Const('k_Cas9byBonA_deg', ratex(k)*0.1)); % 10x faster 
    sys.AddConstant(Const('k_Cas9byAonB_deg',ratex(k)* 0.001)); % 0.1x 
    sys.AddConstant(Const('corpCas9byBonA', 1.6)); % rand
    sys.AddConstant(Const('corpCas9byAonB', 1.6)); % rand
    % sys.AddConstant(Const('cells', 3*10^7*10^3)); % Say yeast OD600 1.0  ~3 x 10^7 cells/ml.
    % Say the device 1L
    % So everything below is for pet cell, now we have 3*10^7*10^3 cells

    % add compositors: these are for the concentrations of species in the system.
    % we use the short-hand version to add the compositors right away,
    % the second parameter is the initial concentration

    dPETndt = sys.AddCompositor('PETn', 100);
    dPETnmin1dt = sys.AddCompositor('PETn_1', 0);
    dPETaseoutdt = sys.AddCompositor('PETaseout', ratex(k)*10);
    dMHETaseoutdt = sys.AddCompositor('MHETaseout', 0);
    dPETaseindt = sys.AddCompositor('PETasein', 0);
    dMHETaseindt = sys.AddCompositor('MHETasein', 0);
    dMHETdt = sys.AddCompositor('MHET', 0);
    dTPAdt = sys.AddCompositor('TPA', 0);
    dEGdt = sys.AddCompositor('EG', 0);
    dGPCRdt = sys.AddCompositor('GPCR', ratex(k)*1000);
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
        Pulse(0, 'PETn', 100), ...     % initial conditions
        % Pulse(0.000000001, 'PETn', 100), ... % spike in 1@ t = 
        Pulse(20000, '', 0), ...     % stop the simulation at time 
    ]);

    T = T./3600;

    % plot
    hold on
    % plot(T, Y(:, sys.CompositorIndex('PETn')))
    plot(T, Y(:, sys.CompositorIndex('GFP')), ...
         T, Y(:, sys.CompositorIndex('mKate')))

    [mKate_max_value,mKate_max_index] = max(Y(:, sys.CompositorIndex('mKate')));
    mKate_max_t = T(mKate_max_index);

    mKate_conc = Y(:, sys.CompositorIndex('mKate'));
    [mKate_min_value,mKate_min_index] = min(mKate_conc);
    mKate_min_t = T(mKate_min_index);

    [GFP_max_value,GFP_max_index] = max(Y(:, sys.CompositorIndex('GFP')));
    GFP_max_t = T(GFP_max_index);

    GFP_conc = Y(:, sys.CompositorIndex('GFP'));
    [GFP_min_value,GFP_min_index] = min(GFP_conc(GFP_max_index:end,1));
    GFP_min_t = T(GFP_min_index)+GFP_max_index;

    my_output(k,1) = mKate_min_value;
    my_output(k,2) = mKate_min_t;
    my_output(k,3) = mKate_max_value;
    my_output(k,4) = mKate_max_t;
    my_output(k,5) = GFP_min_value;
    my_output(k,6) = GFP_min_t;
    my_output(k,7) = GFP_max_value;
    my_output(k,8) = GFP_max_t; 

end

for l = 1:length(ratex)
    sys = BioSystem();

    % define constants. assume that concentration units are in nM and time in s

    sys.AddConstant(Const('k_PETase', ratex(l)*10^-4)); % rand
    sys.AddConstant(Const('k_MHETase', ratex(l)*4.25*10^-3)); % kcat/kM  = (31 ± 0.8 s-1) / (7.3 ± 0.6 ?M) = 4.25 ?M-1 ? s-1
    sys.AddConstant(Const('k_TPA_meta', ratex(l)*10^-3)); % rand
    sys.AddConstant(Const('k_on_MHET', ratex(l)*10^7/(10^9*60))); % 10^7 M-1 min-1 = 10^7/(10^9*60)
    sys.AddConstant(Const('k_off_MHET', ratex(l)*0.05/60)); % 0.05 min-1
    sys.AddConstant(Const('k_MHET_tln',ratex(l)*2.24)); % assume same as normal
    sys.AddConstant(Const('k_tln', ratex(l)*2.24)); % 0.1 nM/h http://homepages.ulb.ac.be/~dgonze/BIONUMBERS/bionumbers.html?fbclid=IwAR0js_SN9QkwW5e0r4ZeJnF96TTh_LqhtKdFVsmaN2l5hMuDOUVHqQM4W0I#transcriptionrate
    sys.AddConstant(Const('k_txn', ratex(l)*0.224)); % 0.01 nM/h
    sys.AddConstant(Const('k_port_PETase', ratex(l)*3.4/60)); % 3.4 min-1 https://onlinelibrary.wiley.com/doi/epdf/10.1002/bit.260330305
    sys.AddConstant(Const('k_port_MHETase', ratex(l)*3.4/60)); % 3.4 min-1 https://onlinelibrary.wiley.com/doi/epdf/10.1002/bit.260330305
    sys.AddConstant(Const('k_GFP_deg', ratex(l)*10^-3)); % 5 hr half life, k = ln2/5hr = ln2/(5*3600), 0.1x
    sys.AddConstant(Const('k_PETase_outdeg', ratex(l)*0.001)); % 0.1x
    sys.AddConstant(Const('k_MHETase_outdeg', ratex(l)*0.001)); %

    % sys.AddConstant(Const('cells', 3*10^7*10^3)); % Say yeast OD600 1.0  ~3 x 10^7 cells/ml.
    % Say the device 1L
    % So everything below is for pet cell, now we have 3*10^7*10^3 cells

    % add compositors: these are for the concentrations of species in the system.
    % we use the short-hand version to add the compositors right away,
    % the second parameter is the initial concentration

    dPETndt = sys.AddCompositor('PETn', 100);
    dPETnmin1dt = sys.AddCompositor('PETn_1', 0);
    dPETaseoutdt = sys.AddCompositor('PETaseout', ratex(l)*10);
    dMHETaseoutdt = sys.AddCompositor('MHETaseout', 0);
    dPETaseindt = sys.AddCompositor('PETasein', 0);
    dMHETaseindt = sys.AddCompositor('MHETasein', 0);
    dMHETdt = sys.AddCompositor('MHET', 0);
    dTPAdt = sys.AddCompositor('TPA', 0);
    dEGdt = sys.AddCompositor('EG', 0);
    dGPCRdt = sys.AddCompositor('GPCR', ratex(l)*1000);
    dMHETGPCRdt = sys.AddCompositor('MHETGPCR', 0);
    dmRNAPETaseMHETaseGFPCas9byAonBdt = sys.AddCompositor('mRNAPETaseMHETaseGFPCas9byAonB', 0);
    dGFPdt = sys.AddCompositor('GFP', 0);

    % TODO: define parts and add them to the system
    % Remember to use the same parts as in the Groundwork section of this problem set
    P1 = Part('All the process', [dPETndt dPETnmin1dt dPETaseoutdt dMHETaseoutdt dPETaseindt dMHETaseindt dMHETdt dTPAdt dEGdt dGPCRdt dMHETGPCRdt dmRNAPETaseMHETaseGFPCas9byAonBdt dGFPdt], ...
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
               Rate('k_tln*k_txn*mRNAPETaseMHETaseGFPCas9byAonB - k_GFP_deg*GFP')]);  % dGFP
    % add parts
    sys.AddPart(P1);

    [T, Y] = sys.run_pulses([...
        Pulse(0, 'PETn', 100), ...     % initial conditions
        % Pulse(0.000000001, 'PETn', 100), ... % spike in 1@ t = 
        Pulse(20000, '', 0), ...     % stop the simulation at time 
    ]);

    T = T./3600;

    % plot

    hold on
    % plot(T, Y(:, sys.CompositorIndex('PETn')))
    plot(T, Y(:, sys.CompositorIndex('GFP')))

    [GFP_max_value,GFP_max_index] = max(Y(:, sys.CompositorIndex('GFP')));
    GFP_max_t = T(GFP_max_index);

    GFP_conc = Y(:, sys.CompositorIndex('GFP'));
    [GFP_min_value,GFP_min_index] = min(GFP_conc(GFP_max_index:end,1));
    GFP_min_t = T(GFP_min_index)+GFP_max_index;

    my_output(l+length(ratex),5) = GFP_min_value;
    my_output(l+length(ratex),6) = GFP_min_t;
    my_output(l+length(ratex),7) = GFP_max_value;
    my_output(l+length(ratex),8) = GFP_max_t;
end

legend('GFP + toggle','mKate  + toggle','GFP at 0.8 rate + toggle','mKate at 0.8 rate + toggle','GFP at 0.6 rate + toggle','mKate at 0.6 rate + toggle', 'GFP no toggle', 'GFP  at 0.8 rate no toggle','GFP at 0.6 rate no toggle' ,'Location', 'best')
xlim([0, 3.5]);
xlabel('Time (hrs)');
ylabel('Concentration (nM)');
title('+/- Toggle for GFP Sensing 100 nM PETn +/- Overall Rate Distrubance');
%% 
subplot(2,2,4)
wtoggleop = my_output((1:length(ratex)),:);
wotoggleop = my_output(((length(ratex)+1):end),:);
GFP_to_mKate = wtoggleop(:,7)./wtoggleop(:,3);
GFP_to_PET = wotoggleop(:,7)./100; % PETn = 100 nM
err = zeros(2,2);
err(1,1)=[GFP_to_mKate(1)-GFP_to_mKate(2)]/GFP_to_mKate(1)*100;
err(2,1)=[GFP_to_mKate(1)-GFP_to_mKate(3)]/GFP_to_mKate(1)*100;
err(1,2)=[GFP_to_PET(1)-GFP_to_PET(2)]/GFP_to_PET(1)*100;
err(2,2)=[GFP_to_PET(1)-GFP_to_PET(3)]/GFP_to_PET(1)*100;

X = categorical({'0.8 Times','0.6 Times'});
X = reordercats(X,{'0.8 Times','0.6 Times'});
bar(X,err)
legend('With Toggle','Without toggle' ,'Location', 'best')
xlabel('Overall Rate Distrubance')
ylabel('Percent Error (%)')
title('Overall Rate Distrubance Percent Error +/- Toggle')

sgtitle('Protein Concentrations With and Without mKate Toggle System With or Without Overall Rate Distrubance');