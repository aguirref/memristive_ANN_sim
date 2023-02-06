set 	fileID 			[sx_open_sim_file_read "/home/users/aguirref/uab_rs_neural_networks/results/NW_784by10_closed_loop/read/FineSim/NW_784by10_closed_loop.fsdb"]
set	i_neuron_file		"/home/users/aguirref/uab_rs_neural_networks/results/NW_784by10_closed_loop/read/FineSim/i_neurons.txt"
set 	i_neuron 		[open $i_neuron_file w+]
set 	signals 		[ sx_signal *neuronb* *tim* ] 
sx_export_mfile "/home/users/aguirref/uab_rs_neural_networks/results/NW_784by10_closed_loop/read/FineSim/NW_784by10_closed_loop.m" $signals
foreach signal $signals {
	set	 neuron_name		[sx_signal_attribute $signal -name]
	set	 neuron_xunit		[sx_signal_attribute $signal -xunit]
	puts	-nonewline		$i_neuron		"$neuron_name "
	puts	-nonewline		"$neuron_name "	
	puts	-nonewline		"$neuron_xunit "		
	}
#puts 	$signals


close $i_neuron

#puts $signal
#sx_export_mfile "/home/users/aguirref/uab_rs_neural_networks/results/NW_19by19_closed_loop/FineSim/NW_19by19_closed_loop.m" $signals
#set avg [ sx_equation “mean($signal)” ]
#puts “average value is $avg”



