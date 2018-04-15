import pandas
from sklearn import model_selection
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFE
import numpy as np
import pickle
import statsmodels.api as sm
from sklearn.metrics import precision_score, recall_score, confusion_matrix, classification_report, accuracy_score, f1_score
from sklearn.metrics import precision_recall_fscore_support

def calculate_recall_precision(label, prediction):
    true_positives = 0
    false_positives = 0
    true_negatives = 0
    false_negatives = 0
 
    for i in range(0, len(label)):
        if prediction[i] == 1:
            if prediction[i] == label[i]:
                true_positives += 1
            else:
                false_positives += 1
        else:
            if prediction[i] == label[i]:
                true_negatives += 1
            else:
                false_negatives += 1
    print("TP=")
    print(true_positives)
    print("FP")
    print(false_positives)
    print("TN")
    print(true_negatives) 
    print("FN")
    print(false_negatives)
    # a ratio of correctly predicted observation to the total observations
    accuracy = (true_positives + true_negatives) \
               / (true_positives + true_negatives + false_positives + false_negatives)
 
    # precision is "how useful the search results are"
    precision = true_positives / (true_positives + false_positives)
    # recall is "how complete the results are"
    recall = true_positives / (true_positives + false_negatives)
 
    f1_score = 2 / ((1 / precision) + (1 / recall))
 
    return accuracy, precision, recall, f1_score

df = pandas.read_csv("real1_training.csv", header = 0)
headers = list(df.columns.values)
array = df.values
x_train=array[:,2:22]
y_train=array[:,22]
train_data = np.asarray(x_train,dtype=np.float32)
train_labels = np.asarray(y_train, dtype=np.int32)

print(np.any(np.isnan(train_data)))
dftest = pandas.read_csv('real1_test.csv', header = 0)
marray = dftest.values
x_test=marray[:,2:22]
y_test=marray[:,22]
test_data = np.asarray(x_test,dtype=np.float32)
test_label = np.asarray(y_test, dtype=np.int32)
print(np.any(np.isnan(test_data)))


logistic = LogisticRegression()
logistic.fit(train_data, train_labels)
rfe = RFE(logistic, 18)
rfe = rfe.fit(train_data, train_labels )
print(rfe.support_)
print(rfe.ranking_)

logit_model=sm.Logit(train_labels,train_data)
#result=logit_model.fit()
result = logit_model.fit(method='bfgs')

print(result.summary())

results = logistic.predict(test_data)
print('Accuracy of logistic regression classifier on test set: {:.4f}'.format(logistic.score(test_data, test_label)))

accuracy, precision, recall, f1_score = calculate_recall_precision(y_test, results)
 
print("Accuracy: ", accuracy)
print("Precision: ", precision)
print("Recall: ", recall)
print("F1 score: ", f1_score)


filename = 'myresultreal1.txt'
with open(filename, 'w') as f:
    for row in results:
        f.write("%s\n" % str(row))

