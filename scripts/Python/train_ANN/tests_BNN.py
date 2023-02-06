# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 11:45:42 2021

@author: ferag
"""
import train_ANN
import tensorflow as tf

train_ANN=train_ANN.train_ANN

layers=[10]
flatten_layer_required=0
ANN_type='DNN'
learning_rate=0.2
quantization=4
num_Epochs=100
 
mnist = tf.keras.datasets.mnist

(x_train, y_train), (x_test, y_test) = mnist.load_data()
x_train, x_test = x_train / 255, x_test / 255

train_ANN('C:\\Users\\ferag\\Downloads\\train-images_16x16.mat', './prueba', layers, ANN_type, flatten_layer_required, num_Epochs, learning_rate, quantization)
