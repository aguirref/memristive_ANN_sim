# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 11:45:42 2021

@author: ferag
"""
import sys
import tensorflow as tf
from tensorflow import keras
import numpy as np
from tensorflow.keras.constraints import max_norm
from os import path
import os
import scipy.io as sio
import larq

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return rray[idx]


# clip model weights to a given hypercube
class ClipConstraint(tf.keras.constraints.Constraint):
    # set clip value when initialized
    def __init__(self, clip_value):
        self.clip_value = clip_value
        
    # clip model weights to hypercube
    def __call__(self, weights):
        return tf.keras.backend.clip(weights, -self.clip_value, self.clip_value)
    # get the config
    def get_config(self):
        return {'clip_value': self.clip_value}


class WeightClip(tf.keras.constraints.Constraint):
    '''Clips the weights incident to each hidden unit to be inside a range
    https://stackoverflow.com/questions/42264567/keras-ml-library-how-to-do-weight-clipping-after-gradient-updates-tensorflow-b/42264773
    '''
    def __init__(self, c=2):
        self.c = c

    def __call__(self, p):
        return tf.keras.backend.clip(p, -self.c, self.c)

    def get_config(self):
        return {'name': self.__class__.__name__,
                'c': self.c}

class quantizeWeights(tf.keras.constraints.Constraint):
    '''Clips the weights incident to each hidden unit to be inside a range
    https://stackoverflow.com/questions/42264567/keras-ml-library-how-to-do-weight-clipping-after-gradient-updates-tensorflow-b/42264773
    '''
    def __init__(self, c=1, levels=16):
        self.c = c
        self.levels = int(0 if levels<2 else levels)
        states_pos = np.linspace(0,self.c,self.levels)
        states_neg = np.linspace(-self.c,0,levels)
        states = np.concatenate((states_neg, states_pos), axis=0)
        states = np.unique(states).T
        states = np.reshape(states,(states.size,1))
        self.states = tf.keras.backend.constant(states)
        #states = np.reshape(states,(7,1))
        
    def __call__(self, weights):
        clipped_weights = tf.keras.backend.clip(weights, -self.c, self.c)
        
        if self.levels>2:
            num_rows, num_cols = tf.keras.backend.int_shape(clipped_weights)
            
            aux_weights = tf.keras.backend.reshape(weights, shape=(1,tf.keras.backend.count_params(weights)))
    
            difference_tensor = tf.keras.backend.abs(tf.math.subtract(aux_weights,self.states))
            min_positions = tf.keras.backend.argmin(difference_tensor, axis=0)
            aux_weights_nearest = tf.transpose(tf.keras.backend.gather(self.states, min_positions))
    
            q_weights =  tf.keras.backend.reshape(aux_weights_nearest, shape=tf.keras.backend.int_shape(weights))
    
            return q_weights
        else:
            return clipped_weights

    def get_config(self):
        return {'name': self.__class__.__name__,
                'c': self.c}
    
def train_ANN(mat_file_DB, mat_file_weights, layers, ANN_type, flatten_layer_required, numEpochs, learning_rate, quantization):

    #data_dir = pjoin(dirname(sio.__file__), 'matlab', 'tests', 'data')
    #mat_fname = pjoin(data_dir, mat_file_DB)
    mat_contents = sio.loadmat(mat_file_DB)

    x_train=mat_contents['images_'+mat_file_DB.split('_')[-1].split('.')[0]].T
    y_train=mat_contents['labels_'+mat_file_DB.split('_')[-1].split('.')[0]].T
    x_test=mat_contents['images_t10k_'+mat_file_DB.split('_')[-1].split('.')[0]].T
    y_test=mat_contents['labels_t10k_'+mat_file_DB.split('_')[-1].split('.')[0]].T
    
    tf.compat.v1.enable_eager_execution()
    model = tf.keras.models.Sequential()
    
    if flatten_layer_required==1:
        model.add(tf.keras.layers.Flatten(input_shape=(28, 28)))
    
    for idx,layer_i in enumerate(layers, start=1):    
        if idx==len(layers):
            activation_fcn='softmax'
        else:
            activation_fcn='sigmoid'
        
        if ANN_type=='BNN':
            model.add(larq.layers.QuantDense(layer_i,
                                  #input_quantizer="ste_sign",
                                  kernel_quantizer="ste_sign",
                                  kernel_constraint="weight_clip",
                                  use_bias='False',
                                  activation=activation_fcn))    
            #model.add(tf.keras.layers.BatchNormalization(scale=False))
            #model.add(tf.keras.layers.Dropout(0.1))
            #model.add(tf.keras.layers.Dense(10, activation='softmax'))
        else:
            model.add(tf.keras.layers.Dense(layer_i,
                                            use_bias='False',
                                            #kernel_constraint = WeightClip(1),
                                            kernel_constraint = quantizeWeights(1,quantization),
                                            activation=activation_fcn))    
            #model.add(tf.keras.layers.BatchNormalization(scale=False))
            #model.add(tf.keras.layers.Dropout(0.1))
            #model.add(tf.keras.layers.Dense(10, activation='softmax'))
    
    opt = keras.optimizers.Adam(lr=learning_rate, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.0)
    model.compile(loss='sparse_categorical_crossentropy',
                  optimizer=opt,
                  metrics=['accuracy'])
    
    model.fit(x_train, y_train, epochs=numEpochs)
    
    model.evaluate(x_test, y_test, verbose=5)
    
    for idx,layer_i in enumerate(layers, start=1):   
        
        if idx==1:
            if ANN_type=='BNN':
                q_weights=[larq.math.sign(model.layers[idx-1].weights[0]).numpy()]
            else:
                q_weights=[model.layers[idx-1].weights[0].numpy()]

            positive_weights=[np.where(q_weights[idx-1]<0, 0, q_weights[idx-1])]
            negative_weights=[np.where(q_weights[idx-1]>0, 0, q_weights[idx-1])*-1]
        else:
            if ANN_type=='BNN':    
                q_weights.append(larq.math.sign(model.layers[idx-1].weights[0]).numpy())
            else:
                q_weights.append(model.layers[idx-1].weights[0].numpy())
            
            positive_weights.append(np.where(q_weights[idx-1]<0, 0, q_weights[idx-1]))
            negative_weights.append(np.where(q_weights[idx-1]>0, 0, q_weights[idx-1])*-1)
    
    for idx,layer_i in enumerate(q_weights):
        if idx==0:
            if flatten_layer_required==1:
                prediction=tf.math.sigmoid(np.matmul(np.reshape(x_test,(10000,784)),positive_weights[idx])-np.matmul(np.reshape(x_test,(10000,784)),negative_weights[idx]))
            else:
                prediction=tf.math.sigmoid(np.matmul(x_test,positive_weights[idx])-np.matmul(x_test,negative_weights[idx]))                
        elif idx<len(q_weights)-1:
            prediction=tf.math.sigmoid(np.matmul(prediction,positive_weights[idx])-np.matmul(prediction,negative_weights[idx]))
        else:
            prediction=np.matmul(prediction,positive_weights[idx])-np.matmul(prediction,negative_weights[idx])
    
    prediction=np.argmax(prediction,1)
    accuracy=np.sum(prediction==y_test.T)/y_test.size
    
    if path.exists(mat_file_weights)==False:
        os.makedirs(mat_file_weights)
    sio.savemat(''+mat_file_weights+'/G_values.mat', dict(G_real=q_weights, G_pos=positive_weights, G_neg=negative_weights, Accuracy_mat=accuracy*100, labels_test=y_test.T, images_test=x_test.T))
