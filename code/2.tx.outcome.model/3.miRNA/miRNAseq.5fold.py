#### Step2 Model Training and Predictions

### five-fold cross validation without PCA 
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import pandas as pd
import sklearn
from sklearn import tree, preprocessing, ensemble
from sklearn.cross_validation import train_test_split 
from sklearn.metrics import accuracy_score
import tensorflow as tf
from tensorflow.contrib import learn
import numpy as np
from sklearn import svm
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn import feature_selection


import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline
from sklearn import datasets, linear_model
from sklearn import cross_validation

def main():
	for k in range(20):
		for n in range(1,6):
			## deep learning
			def dnnclassifier():
				tf.logging.set_verbosity(tf.logging.INFO)
			 	traindata = pd.read_csv("~/Project/ov_prognosis_drugprecision/result/4.prediction/3.miRNAseq/classifier.trainset.5fold" + str(n) +".csv")
			 	y_train = traindata['txoutcome']
			 	X_train = traindata[list(range(2,len(traindata.columns)))]
			 	testdata = pd.read_csv("~/Project/ov_prognosis_drugprecision/result/4.prediction/3.miRNAseq/classifier.testset.5fold" + str(n) +".csv")
			 	y_test = testdata['txoutcome']	
			 	X_test = testdata[list(range(2,len(traindata.columns)))]
			 	feature_columns=learn.infer_real_valued_columns_from_input(X_train)
			 	dnn_classifier = learn.DNNClassifier(hidden_units=[20, 40, 20], n_classes=5,feature_columns=feature_columns)
			 	dnn_classifier.fit(X_train, y_train, steps = 2000)
			 	dnn_prediction = dnn_classifier.predict(X_test)
			 	print('DNN Prediction Score: {0}'.format( accuracy_score(dnn_prediction, y_test)))
			 	print(len(dnn_prediction))
			 	print(len(y_test))
			 	print(dnn_prediction[4])
			 	print(y_test[4])
				# save the predicted value for the next step of C-index calculation by R
			 	fout = open("/Users/Zhangzhe/Project/ov_prognosis_drugprecision/result/4.prediction/3.miRNAseq/dnn_classifier/dnn_classifier.txoutcome.5fold"+str(5*k+n)+".txt","w")
			 	for j in range(len(dnn_prediction)):
			 		fout.write(str(y_test[j])+'\t' + str(dnn_prediction[j])+'\n')
			dnnclassifier()
			print(k)


			## Decision Tree
			def treeclassifier():
				tf.logging.set_verbosity(tf.logging.INFO)
				traindata = pd.read_csv("~/Project/ov_prognosis_drugprecision/result/4.prediction/3.miRNAseq/classifier.trainset.5fold" + str(n) +".csv")
				y_train = traindata['txoutcome']
				X_train = traindata[list(range(2,len(traindata.columns)))]
				testdata = pd.read_csv("~/Project/ov_prognosis_drugprecision/result/4.prediction/3.miRNAseq/classifier.testset.5fold" + str(n) +".csv")
				y_test = testdata['txoutcome']	
				X_test = testdata[list(range(2,len(traindata.columns)))]
				dt_classifier = tree.DecisionTreeClassifier(max_depth=3) 
				dt_classifier.fit(X_train, y_train)
				dt_prediction = dt_classifier.predict(X_test)
				print (dt_prediction)
				print(y_test)
				# save the predicted value for the next step of C-index calculation by R
				fout = open("/Users/Zhangzhe/Project/ov_prognosis_drugprecision/result/4.prediction/3.miRNAseq/dt_classifier/dt_classifier.txoutcome.5fold"+str(5*k+n)+".txt","w")
				for j in range(len(dt_prediction)):
					fout.write(str(y_test[j]) + '\t' + str(dt_prediction[j]) +'\n')
			treeclassifier()	

			## Random Forest
			def rfclassifier():
				tf.logging.set_verbosity(tf.logging.INFO)
				traindata = pd.read_csv("~/Project/ov_prognosis_drugprecision/result/4.prediction/3.miRNAseq/classifier.trainset.5fold" + str(n) +".csv")
				y_train = traindata['txoutcome']
				X_train = traindata[list(range(2,len(traindata.columns)))]
				testdata = pd.read_csv("~/Project/ov_prognosis_drugprecision/result/4.prediction/3.miRNAseq/classifier.testset.5fold" + str(n) +".csv")
				y_test = testdata['txoutcome']	
				X_test = testdata[list(range(2,len(traindata.columns)))]
				rf_classifier = ensemble.RandomForestClassifier(n_estimators=40)
				rf_classifier.fit(X_train, y_train)
				rf_prediction = rf_classifier.predict(X_test)
				# save the predicted value for the next step of C-index calculation by R
				fout = open("/Users/Zhangzhe/Project/ov_prognosis_drugprecision/result/4.prediction/3.miRNAseq/rf_classifier/rf_classifier.txoutcome.5fold"+str(5*k+n)+".txt","w")
				for j in range(len(rf_prediction)):
					fout.write(str(y_test[j]) + '\t'+str(rf_prediction[j]) +'\n')
			rfclassifier()



main()