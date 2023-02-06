#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 11:14:51 2021

@author: aguirref
"""

import numpy as np
import tensorflow as tf

c=1
levels=4

weights = tf.keras.backend.random_uniform((4,4), minval=-1.0, maxval=1.0)

states_pos = np.linspace(0,c,levels)
states_neg = np.linspace(-c,0,levels)
states = np.concatenate((states_neg, states_pos), axis=0)
states = np.unique(states).T
states = np.reshape(states,(7,1))
states = tf.keras.backend.constant(states)

aux_weights = tf.keras.backend.reshape(weights, shape=(1,tf.keras.backend.count_params(weights)))

difference_tensor = tf.keras.backend.abs(tf.math.subtract(aux_weights,states))
min_positions = tf.keras.backend.argmin(difference_tensor, axis=0)
                            

aux_weights_nearest = tf.transpose(tf.keras.backend.gather(states, min_positions))
q_weights =  tf.keras.backend.reshape(aux_weights_nearest, shape=tf.keras.backend.int_shape(weights))

