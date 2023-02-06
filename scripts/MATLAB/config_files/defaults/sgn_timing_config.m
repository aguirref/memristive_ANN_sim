    read_delay=10e-6;
    sys_clk_WR_margin=[30e-6 60e-6];
    sys_clk=[100e-6 200e-6];
    time_offset=1e-6;

    % input_vector_freq: Frequency (Hz) at which MNIST images are supplied
    % to the CPA being simulated
    input_vector_freq=[100e3];
    
    % in_vect_r_f_time: Rise/falling edges of the voltage signals
    % comprising the MNIST image (n² signals) being supplied to the CPA
    in_vect_r_f_time=1e-7;
    
    error_band=0.01;
    latency_tr=10e-9;