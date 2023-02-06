function I_neuron=avg_noisy_signal(time_vector, I_neuron, sim_settings, prb_points)
    time_image_acc=sim_settings.time_offset;
    in_vect_r_f_time=sim_settings.in_vect_r_f_time;
    input_vector_freq=sim_settings.input_vector_freq;
    number_of_images=sim_settings.number_of_images;

    for image_cnt_i=1:number_of_images
        t_start_image=time_image_acc+in_vect_r_f_time;
        t_end_image=time_image_acc+in_vect_r_f_time+1/input_vector_freq;
        time_image_acc=t_end_image;
        
        time_idx=find(time_vector>t_start_image & time_vector<t_end_image);

        for neuron_i=1:size(I_neuron,1)
            I_neuron(neuron_i,time_idx)=ones(size(I_neuron(neuron_i,time_idx)))*mean(I_neuron(neuron_i,time_idx(round(linspace(1,length(time_idx),prb_points)))));
        end
    end
end