%% Simulation Settings
    % run_sim: This enables (1) / disables (0) the simulation procedures.
    % If set to 0, the scripts only runs the simulation if not simulation
    % results are available for the selected CPA/Sim settings. Otherwise
    % (run_sim=1) it overrides the existing simulation data.
    run_sim=1;

%% Plot options

    % gen_plot: This enables (1) / disables (0) the plot generation during
	% runtime
    gen_plot=0;


%% Saving options
    
% close_after_save: This enables closing the confusion matrix plots one
    % after having save the image to avoid lagging the system.
    close_after_save=0;
    
    save_write_signals='yes';
    save_h='yes';
    save_vpn='no';
    save_inner_neuron='yes';
    save_imemd='no';
    save_all=1;
    debug_mode='no';

    % remove_RAW_after_sim: This options enables/disables the removal of
    % the .fsdb RAW data generated during the simulation. It is extremely
    % usefull during MC runs, to avoid generating huges amounts of binary
    % data that saturates the HD capacity.
    remove_RAW_after_sim='yes';    

    test_write='diagonal';
    
%% Correction algorithms

    % correct_rs: This enables ('yes') / disables ('no') the series resistance
    % compensation of the ex-situ calculated synaptic Weights
    correct_rs='no';  
    correct_mode='mode6';
    init_cal_criterion=[0.1];

%% Remapping algorithms
    
    %remap_W:
    remap_W='no';
    
    remap_mode='mode7';

%% Analysis     

    % calcl_pwr: Writes the meas statements to log the power consumption in
    % each CPA device (parasitic and memdiodes)
    calc_pwr='no';

    % calc_WR_margin: calculates the write margins for the NN
    calc_WR_margin='no';

    % calc_RE_margin: calculates the read margins for the NN
    calc_RE_margin='no';

    % calc_RE_latency: calculates the read latency the NN
    calc_RE_latency='no';    
    
    % calc_inference_accuracy
    calc_inference_accuracy='yes';

    %Noise Analysis
    noise_transient='no';
    noisefmax=1e5; %This should be 10e6
    noisefmin=1e2;
    noisescale=1; 
    
    % scale: This variable (ranging from 0 to 1) defines the ratio between
    % currents in the corrected CPA with Rs with regard to the CPA with NO Rs
    % effects. Only relevant if correct_rs='yes'
    scale=1;
    
    %timing_stats: allows repeating the inference simulation for the
    %statistical analysis of the runing time.
    timing_stats='no';
    
    % write: This option defines whether the scripts writes the CPA (1) or
    % performes a read (inference) with the already loadad synaptic weights
    % (0)
    write=0;

    % use_write_data: This options defines wheter the inference operation
    % is done with the ideal calculated synaptic weights (0), or with the
    % synaptic weights effectivelly written to the memdiodes (if they are not
    % yet written, the scripts performs a write operation (write=1), obtains
    % the G matrix, loads in the CPA, and then runs the inference). If the 
    % CPA is large, this might be very time consuming.
    use_write_data=0;

    % simEngine: This option defines the SPICE simulator being used.
    % SPICE engines supported are HSPICE, FineSim and LTSpice
    simEngine='FineSim';%'HSPICE';%
    
   
    %N_Proc= Number of CPU processors used for running the simulation
    N_Proc=4;
       
    simulation_read_time=(number_of_images(1)*(1./input_vector_freq(1)+in_vect_r_f_time))*1.05+read_delay;
    simulation_WR_margin_time=65e-3;  
    simulation_write_time=300e-3; 
    
   
