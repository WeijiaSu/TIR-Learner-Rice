import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, VotingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
from sklearn.neural_network import MLPClassifier
from sklearn import preprocessing as pp
import warnings
import time
import os


# remove warnings
from sklearn.tree import DecisionTreeClassifier


def warn(*args, **kwargs):
    pass


warnings.warn = warn
warnings.filterwarnings("ignore")

train = pd.read_csv('Rice_all.csv',header=0)
test = pd.read_csv('Rice_TestSetFinal.csv',header=0)

#train = pd.read_csv('/Users/weijiasu/Dropbox/Research/BioinformaticsProject/TE_ML/SorghumData/SorghumTrain_0712.csv',header=0)
#test = pd.read_csv('/Users/weijiasu/Dropbox/Research/BioinformaticsProject/TE_ML/SorghumData/SorghumTest_0712.csv', header=0)

all=train.append(test,ignore_index=True)

trainx = train.drop(['ID', 'Target'], axis=1)
trainy = train['Target']
#
testx = test.drop(['ID','Target'], axis=1)
testy = test['Target']

all_x=all.drop(['ID', 'Target'], axis=1)
all_y=all["Target"]

###################################################### Data_To Predict ######################################################
# topre = pd.read_csv('C:/Users/Zijian Zhao/Downloads/chr1_NonHomopre.csv', header=0)
# toprex = topre.drop(['ID'], axis=1)



def knn(nb, w, a, s):
    trx = trainx
    tex = testx
    if s == 'True':
        trx = pd.DataFrame(pp.MinMaxScaler().fit_transform(trainx))
        tex = pd.DataFrame(pp.MinMaxScaler().fit_transform(testx))
    model = KNeighborsClassifier(n_neighbors=nb, weights=w, algorithm=a)
    model.fit(trx, trainy)
    result = model.score(tex, testy)
    results = cross_val_score(model, trx, trainy, cv=10)
    pre = model.predict(tex)
    print('=================================================')
    print('Model Name: K Nearest Neighbors\nCV Accuracy:', results.mean())
    print('Test Accuracy:', result)
    print('Best Model Parameter(s):\n- n_neighbors:', nb, '\n- weights:', w, '\n-algorithm:', a)
    print('Confusion Matrix: ')
    print(pd.crosstab(pre, testy))

    return results.mean()



nb = [3, 5, 10, 30, 50, 70]
w = ['uniform', 'distance']
al = ['ball_tree', 'kd_tree', 'brute']
sc = ['True', 'False']


def knntune(nb, w, a, b):
    max = 0
    for i in nb:
        for j in w:
            for k in a:
                for l in b:
                    temp = knn(i, j, k, l)
                    if temp > max:
                        result = []
                        max = temp
                        result.append([max, i, j, k, l])
    print(result)



seed = 7


def en5(file,vo, w, s, num):
    topre = pd.read_csv(file, header=0)
    toprex = topre.drop(['ID'], axis=1)
    trx = all_x
    ty = all_y
    if s == 'True':
        trx = pd.DataFrame(pp.MinMaxScaler().fit_transform(trainx))
    m1 = MLPClassifier(hidden_layer_sizes=(30, 100), activation='tanh', learning_rate='adaptive',
                       learning_rate_init=0.001, momentum=0.9, random_state=seed)
    m2 = KNeighborsClassifier(n_neighbors=3, weights='distance', algorithm='ball_tree')
    m3 = RandomForestClassifier(n_estimators=200, max_features='auto', max_depth=None,
                                criterion='gini', bootstrap=True, random_state=seed)
    m4 = AdaBoostClassifier(DecisionTreeClassifier(criterion='gini', max_depth=50, max_features='auto'),
                            n_estimators=20,
                            learning_rate=0.0001, random_state=seed)
    ensemblecf = VotingClassifier(
        estimators=[('nn', m1), ('knn', m2), ('rf', m3), ('ada', m4)], voting=vo,
        weights=w)
    ensemblecf = ensemblecf.fit(trx, ty)
    pre = ensemblecf.predict(toprex)
    p = pd.DataFrame(pre,columns=['prediction'])
    r = pd.DataFrame(topre['ID']).join(p)
    r.to_csv('prediction_chr_'+str(num)+'.csv',index=None)


for i in range(1,13):
    print("predicting chr %s"%(i))
    t = time.time()
    en5("chr%s_NonHomoToPre.csv"%(i),'soft', [1, 1, 2, 1], 'False', i)
    print("chr %s finished, time:"%(i)+str(time.time() - t))



vo = ['hard', 'soft']
w5 = [[1, 1, 1, 1], [1, 1, 5, 5], [1, 1, 2, 1], [1, 1, 10, 1], [1, 1, 1.5, 2], [1, 10, 1, 1],
      [1, 1, 1, 5], [5, 1, 1, 1], [10, 4, 1, 4]]
s = ['True', 'False']


def eutune(vo, w, sc):
    max = 0
    for i in vo:
        for j in w:
            for k in sc:
                temp = en5(i, j, k)
                if temp > max:
                    result = []
                    max = temp
                    result.append([max, i, j, k])
    print(result)



def nn(parameter):
    trx = trainx
    tex = testx
    if parameter[5] == 'True':
        trx = pd.DataFrame(pp.MinMaxScaler().fit_transform(trainx))
        tex = pd.DataFrame(pp.MinMaxScaler().fit_transform(testx))
    model = MLPClassifier(hidden_layer_sizes=parameter[0], activation=parameter[1], learning_rate=parameter[2],
                          learning_rate_init=parameter[3],
                          momentum=parameter[4],
                          validation_fraction=0.1, random_state=seed)
    results = cross_val_score(model, trx, trainy, cv=10)
    var = np.var(results)
    print('Variance:', var)
    model.fit(trx, trainy)
    result = model.score(tex, testy)
    pre = model.predict(tex)
    print('=================================================')
    print('Model Name: Neural Network\nCV Accuracy:', results.mean(), '\nTest Accuracy:', result)
    print('Best Model Parameter(s):\n- hidden_layer_size:', parameter[0], '\n- activation :', parameter[1],
          '\n- learning rate:', parameter[2], '\n- learning_rate_init:', parameter[3], '\n- momentum:', parameter[4],'\n- Scale:',parameter[5])
    print('Confusion Matrix: ')
    print(pd.crosstab(pre, testy))
    return results.mean()


size = [(10, 10), (100, 100), (30, 100)]
act = ['tanh', 'relu']
lr = ['constant', 'adaptive']
lri = [0.001, 0.01]
mo = [0.1, 0.9]
s = ['True', 'False']


def nntune(size, act, lr, lri, mo, scale):
    max = 0
    for i in size:
        for j in act:
            for k in lr:
                for a in lri:
                    for b in mo:
                        for c in scale:
                            temp = nn([i, j, k, a, b, c])
                            if temp > max:
                                result = []
                                max = temp
                                result.append([max, i, j, k, a, b, c])
    print(result)


