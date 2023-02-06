%% Configuration of the Dataset

    % database: MNIST / F-MNIST / K-MNIST / CIFAR-10 / SVHN / YALE-FACE /
    % YALE-FACE-B / ORL
    database='MNIST';
    
    % limit_output
    limit_output=1;
    
    %input polarity: Input signals are bipoar or unipolar
    in_polarity='unipolar';
    
    % number_of_images: Number of MNIST images being supplied during the
    % simulation procedure. MNIST-like: 10000 max, ORL 80 max, Yale 42 max,
    % Yale B 387 max
    number_of_images=[1000];

%% Configuration of the training procedure

    %Learning algorithm; defined the learning algorithm to be used for
    %training the Network. So far we are using the Scaled Conjugate
    %Gradient (trainscg)
    learn_algorithm='trainscg';
    
    %all_logsig:
    all_logsig=0;
    
    train_tool='MATLAB';
    
    learning_rate=0.1;
    quantization=64;
    num_Epochs=10;
    
    train_MC_vector=[1];

%% Structure of the network

    ANN_type='DNN';
    neural_layers={'A' 'B'};

    %MNIST-like / CIFAR-10 / SVHN: 10 clases, Yale: 15 clases, YaleB: 38
    %clases, ORL:40 classes
    outputs_n=10;

    %Hidden layers
    hidden_layer_1=54;
    hidden_layer_2=50;
    hidden_layer_3=25;
    hidden_layer_4=14;
    hidden_layer_5=14;
    hidden_layer_6=14;

    %Resolution of the input layer
    image_size_mat=[
%                     4 4;...
% 		            %5 5;...
%                     6 6;...
% 		            %7 7;...
                     8 8;...
%                     9 9;...
%                     10 10;...
% 		            %11 11;...
                     12 12;...
% 		            %13 13;...
%                     14 14;...
% 		            %15 15;...
                     16 16;...
% 		            %17 17;...
%                     18 18;...
% 		            %19 19;...
                     20 20;...
% 		            %21 21;...
%                     22 22;...
% 		            %23 23;...
                     24 24;...
% 		            %25 25;...
%                     26 26;...
% 		            %27 27;...
                     28 28
                   ];

    
    
