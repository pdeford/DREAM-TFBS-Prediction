#!/usr/bin/env python
"""
Bayesian neural network using mean-field variational inference.
(see, e.g., Blundell et al. (2015); Kucukelbir et al. (2016))
Inspired by autograd's Bayesian neural network example.

Probability model:
  Bayesian neural network
  Prior: Normal
  Likelihood: Normal with mean parameterized by fully connected NN
Variational model
  Likelihood: Mean-field Normal
"""
#from __future__ import absolute_import
#from __future__ import division
#from __future__ import print_function

import edward as ed
import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf
import tables

from edward.models import Normal
from edward.stats import norm
from sklearn.linear_model import LogisticRegression as logit
from sklearn.metrics import roc_auc_score, average_precision_score
#plt.style.use('ggplot')


def neural_network(x, W_0, W_1, b_0, b_1):
    h = tf.nn.tanh(tf.matmul(x, W_0) + b_0)
    h = tf.matmul(h, W_1) + b_1
    return tf.reshape(h, [-1])

# DATA
# name = "/Users/pdeford/Desktop/DREAM/DREAM-TFBS-Prediction/ATF2_GM12878_train.h5"
name = "/Users/pdeford/Desktop/DREAM/DREAM-TFBS-Prediction/ATF2_down.h5"

f = tables.open_file(name, mode='r')
y = f.root.Y[:]
x = f.root.X[:]

n1 = 1000

all_indices = np.asarray(range(y.shape[0]))
np.random.shuffle(all_indices)

training_indices = list(all_indices[:n1])
full_train = list(all_indices[n1:])

x_train, y_train = x[full_train,:], y[full_train]
x_test, y_test = x[training_indices,:], y[training_indices]

N = x_train.shape[0]   # num data ponts
D = x_train.shape[1]   # num features

# MODEL
W_0 = Normal(mu=tf.zeros([D, 2]), sigma=tf.ones([D, 2]))
W_1 = Normal(mu=tf.zeros([2, 1]), sigma=tf.ones([2, 1]))
b_0 = Normal(mu=tf.zeros(2), sigma=tf.ones(2))
b_1 = Normal(mu=tf.zeros(1), sigma=tf.ones(1))

x = tf.convert_to_tensor(x_train, dtype=tf.float32)
y = Normal(mu=neural_network(x, W_0, W_1, b_0, b_1),
           sigma=0.1 * tf.ones(N))

# INFERENCE
qW_0 = Normal(mu=tf.Variable(tf.random_normal([D, 2])),
              sigma=tf.nn.softplus(tf.Variable(tf.random_normal([D, 2]))))
qW_1 = Normal(mu=tf.Variable(tf.random_normal([2, 1])),
              sigma=tf.nn.softplus(tf.Variable(tf.random_normal([2, 1]))))
qb_0 = Normal(mu=tf.Variable(tf.random_normal([2])),
              sigma=tf.nn.softplus(tf.Variable(tf.random_normal([2]))))
qb_1 = Normal(mu=tf.Variable(tf.random_normal([1])),
              sigma=tf.nn.softplus(tf.Variable(tf.random_normal([1]))))

data = {y: y_train}
inference = ed.MFVI({W_0: qW_0, b_0: qb_0,
                     W_1: qW_1, b_1: qb_1}, data)


sess = ed.get_session()
init = tf.initialize_all_variables()
init.run()

# RUN MEAN-FIELD VARIATIONAL INFERENCE
inference.run(n_iter=500, n_samples=5, n_print=100)


# GET FITS, AND LEARN LOGISTIC REGRESSION MODEL ON OUTPUT
mus = neural_network(x_train, qW_0.sample(), qW_1.sample(),
                     qb_0.sample(), qb_1.sample())
outputs = mus.eval()

clf = logit()
clf.fit(outputs.reshape(-1, 1), y_train)

# SCORE THE PERFORMANCE OF THE FULL MODEL
mus = neural_network(x_test, qW_0.sample(), qW_1.sample(),
                     qb_0.sample(), qb_1.sample())
outputs = mus.eval()

Y2 = clf.predict_proba(outputs.reshape(-1, 1))[:,1]
score1 = average_precision_score(y_test, Y2)
score2 = roc_auc_score(y_test, Y2)

print "auPRC:", score1
print "auROC:", score2