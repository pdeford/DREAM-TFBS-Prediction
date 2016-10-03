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

#/home/pdeford1/anaconda2/envs/DREAM-TFBS/lib/python2.7/site-packages/edward/inferences.py
#Modified `edward` to allow `feed_dict` to be passed to Inference.run()
import edward as ed
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
print "LOADING DATA"
name = "/cache/DREAM-tfbs-pred/out/ATF2_GM12878_train.h5"
n1 = 100000
n2 = 1000000

f = tables.open_file(name, mode='r')
y = f.root.data.Y[:]
x = f.root.data.X[:]
all_indices = np.asarray(range(y.shape[0]))
pos_sample = all_indices[y==1]
neg_sample = all_indices[y==0]

np.random.shuffle(neg_sample)
np.random.shuffle(all_indices)
training_indices = (list(all_indices[:n1]))
pos_sample = list(pos_sample[~np.in1d(pos_sample, training_indices)])
neg_sample = list(neg_sample[~np.in1d(neg_sample, training_indices)])
left = n2 - len(pos_sample)
full_train = sorted(pos_sample + neg_sample[:left])
x_test = x[training_indices,:]
y_test = y[training_indices]
x_train = x[full_train,:]
y_train = y[full_train]
f.close()

x_place = tf.placeholder(tf.float32, shape=x_train.shape)
y_place = tf.placeholder(tf.int8, shape=y_train.shape)

N = x_train.shape[0]   # num data ponts
D = x_train.shape[1]   # num features

# MODEL
print "DEFINING MODEL"
W_0 = Normal(mu=tf.zeros([D, 2]), sigma=tf.ones([D, 2]))
W_1 = Normal(mu=tf.zeros([2, 1]), sigma=tf.ones([2, 1]))
b_0 = Normal(mu=tf.zeros(2), sigma=tf.ones(2))
b_1 = Normal(mu=tf.zeros(1), sigma=tf.ones(1))

x = tf.convert_to_tensor(tf.Variable(x_place), dtype=tf.float32)
y = Normal(mu=neural_network(x, W_0, W_1, b_0, b_1),
           sigma=0.1 * tf.ones(N))

# INFERENCE
print "PREPARING TO CARRY OUT INFERENCE"
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
init.run(feed_dict={x_place: x_train})

# RUN MEAN-FIELD VARIATIONAL INFERENCE
print "INFERENCE"
inference.run(n_iter=5000, n_samples=5, n_print=5001, feed_dict={x_place: x_train})

# GET FITS, AND LEARN LOGISTIC REGRESSION MODEL ON OUTPUT
print "TRAIN LOGIT"
mus = neural_network(x_train, qW_0.sample(), qW_1.sample(),
                     qb_0.sample(), qb_1.sample())
outputs = mus.eval()

clf = logit()
clf.fit(outputs.reshape(-1, 1), y_train)

# SCORE THE PERFORMANCE OF THE FULL MODEL
print "EVAL PERFORMANCE"
mus = neural_network(x_test, qW_0.sample(), qW_1.sample(),
                     qb_0.sample(), qb_1.sample())
outputs = mus.eval()

Y2 = clf.predict_proba(outputs.reshape(-1, 1))[:,1]
score1 = average_precision_score(y_test, Y2)
score2 = roc_auc_score(y_test, Y2)

print "auPRC:", score1
print "auROC:", score2
