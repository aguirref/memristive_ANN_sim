%% Cross-Point-Array (CPA) settings
    dualSide_connect_h='yes';
    dualSide_connect_v='no';    
    pos_neg={'pos' 'neg'};
    partitions_N=[8 1];
    R_shunt=10;
    R_cs=1;
    series_resistance=[1e-3 1e-2 0.1];
    read_voltage=0.2;%[0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75];
    write_voltage=[1.3];
    dev_polarity='positive';
    wrFrequency_vec=2.5e4;
    Dcycle_vec=[0.25];
    Dcycle_v=0.25;
    trise=1e-6;
    tdelay_gral=10e-6;
    
    digits_order=[1 1;
                  2 2;
                  3 3;
                  4 4;
                  5 5;
                  6 6;
                  7 7;
                  8 8;
                  9 9;
                  10 10]; 
              
     digits_order=[[1:1:outputs_n];[1:1:outputs_n]]'; 

%% Memristor model

    % model_ver: This option defines the kind of model being used.
    % supported options are 'QMM'/'DMM'/'QMM_SBSF' / 'yakopic' / 'laiho_biolek' /
    % 'univ_michigan' / 'lineal'
    model_ver='DMM'; 
    
    % memdiode_model_array: List of memdiode models used for the simulation
    memdiode_model_array={'memdiode_HSPICE_DMM_oxidation_0p1','memdiode_HSPICE_DMM_IGZO-sputtered'};
    
    % timeSTEP_sim_memdiode: used to define the maximum timestep allowed
    % for the simulation of the quasi-static IV loop of the device.
    timeSTEP_sim_memdiode=1e-3;

    %mode: mapping mode. It sets the mapping mode and normalization mode.
    %Available modes are: memdiode<1|2>, log<1|2> and lin<1|2>
    mode='memdiode2'; %other options are 'log' o 'memdiode' 'lin'
    numSigma=4;
    H_max=1;
    H_min=1e-2;

%% Neuron models    

    neurons='analytic';
    senseMode='TIA';

    %TIA_gain: This is the "amplification" of the output from the crossbar.
    %suggested values are 1e6 (for rather normal I-V loops). For the
    %samples from KAUST this values also works in the SLP. For the MLP
    %higher values should be considered. For example with TIA=1e9. When
    %using he Rsense mode lower values can be used (1e4 for instance). Also for the high current model from KAUST, 1e4 should be ok.even for TIA MODE
    TIA_gain=1e8; 

    %Rsense: sense the current and translates it into a voltage. the
    %recomended value is 1e3. For cases of very low currents (KAUST
    %devices), higher values are recomended
    Rsense=1e3;

    %This is the input impedance of the TIA. with high current devices,
    %this could be set to 1. In fact, the behavioural logsig neuron 
    R_in_TIA=0.001;
    meas_mode='voltage';
    sel_type='ideal';

    %This sets the threshold at which changes in the circuit signals are
    %ignored. This is very importante for devices with very low currents.
    print_I_tolerance = 1e-12;
    print_V_tolerance = 1e-4;

%% latency calculations

    Csweep='one2one';
    use_ftSize='no';
    C_memdiode=45e-18;
    feature_size=[38e-9];
    C_interline=[11.48e-15];

%% Variability simulations
    
    % sim_variability: This option enables the simulation of a D2D
    % variability
    sim_variability='no';    
    
    variability_mode='mode1';
    variability_MC=5;  
    variabilities=[0.05:0.15:0.3];
    stdImax_array=0.3;
    stdImin_array=0.3;

%% Reliability simulations

    % sim_degradation: This option enables the simulation of a degraded CPA
    sim_degradation='no';    

    STCKFratio=[0.01 0.15 0.3];
    STCKF_total_ratio=0.15;
    flr_type='STCKF_ON';
    MC_total=1;
