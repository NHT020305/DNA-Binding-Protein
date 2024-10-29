import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random as random
import seaborn as sns
import statsmodels.api as sm


from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier, plot_tree
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier


from amino_acids import plot_amino_acids_distribution
from hydropathy import hydropathy
from net_charge import net_charge
from k_mer import Kmer_featurization
from k_mer import plot_kmers_distribution
from feature_distribution import plot_feature_distribution
from performance_metrics import sketch_confusion_matrix, get_performance_metrics
from performance_metrics import plot_roc_curve
from feature_scaling import normalize, standardize, robust_scaling






####################  READ AND SET UP DATASET  ####################


amino_acids_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 
                    'G', 'H', 'I', 'L', 'K', 'M', 'F', 
                    'P', 'S', 'T', 'W', 'Y', 'V']


characters_to_remove = ['X', 'U', "O", "B"]


def remove_character(input_string, chars_to_remove):
    for char in chars_to_remove:
        input_string = input_string.replace(char, '')
    return input_string


def set_up_fasta_data(file_path, name):
    
    data = []
    records = SeqIO.parse(file_path, "fasta")

    for record in records:
        data.append([record.id, record.seq])
    
    data = pd.DataFrame(data, columns=["Label", "Sequence"])
    data['Label'] = data['Label'].apply(lambda x: x[-1])
    data['Length'] = data['Sequence'].apply(len)

    for amino_acid in amino_acids_list:
        data['Frequency_' + amino_acid] = data['Sequence'].apply(lambda x: 
                                            x.count(amino_acid) / len(x))
        
    data['Sequence'] = data['Sequence'].apply(lambda x: 
                                remove_character(x, characters_to_remove))
    
    # Create new column: Hydropathy
    data['Hydropathy'] = data['Sequence'].apply(lambda x: 
                        hydropathy(x, amino_acids_list))

    # Create new column: Net_charge
    data['Net_charge'] = data['Sequence'].apply(lambda x: net_charge(x))

    # Define empty lists to store the computed values for each column
    molecular_weights = []
    instability_indexes = []
    aromaticities = []
    flexibilities = []
    isoelectric_points = []
    sec_structure_fractions = []

    # Iterate through each protein sequence in the 'Sequence' column
    for seq in data['Sequence']:
        X = ProteinAnalysis(seq)
    
        # Append the computed values for each property to the respective list
        molecular_weights.append(X.molecular_weight())
        instability_indexes.append(X.instability_index())
        aromaticities.append(X.aromaticity())
        flexibilities.append(X.flexibility())  # Flexibility is a list
        isoelectric_points.append(X.isoelectric_point())
        sec_structure_fractions.append(X.secondary_structure_fraction())

    # Add these lists as new columns to your DataFrame
    data['Molecular Weight'] = molecular_weights
    data['Instability Index'] = instability_indexes
    data['Aromaticity'] = aromaticities
    data['Flexibility(Mean)'] = [sum(flex)/len(flex) if len(flex) > 0 
                                 else 0 for flex in flexibilities]
    data['Isoelectric Point'] = isoelectric_points
    data['Secondary Structure Fraction'] = [sum(fraction)/len(fraction) 
                                for fraction in sec_structure_fractions]
    
    data.to_csv(name, index=False)

    return data



# Read dataset(.fasta file) and convert to csv
#dna_train = set_up_fasta_data("Data/Train.fasta", "Data/DNA_Train.csv")
#dna_test = set_up_fasta_data('Data/Test.fasta', "Data/DNA_Test.csv")
dna_train = pd.read_csv('Data/DNA_Train.csv')
dna_test = pd.read_csv('Data/DNA_Test.csv')


# Convert 'Label' to a categorical variable
dna_train['Label'] = pd.Categorical(dna_train['Label'])
dna_test['Label'] = pd.Categorical(dna_test['Label'])


# Remove missing values from the training data
dna_train = dna_train.dropna()


# Ensure the labels are treated as binary
dna_train['Label'] = dna_train['Label'].cat.codes  # 0 for "0" and 1 for "1"
dna_test['Label'] = dna_test['Label'].cat.codes






####################  FEATURE DISTRIBUTION  #################### 

#plot_feature_distribution(dna_train, "Length", 0, 5000)

#plot_amino_acids_distribution(dna_train, 1, amino_acids_list)

#plot_amino_acids_distribution(dna_train, 0, amino_acids_list)

#plot_kmers_distribution(dna_train, 5000, 1000, amino_acids_list, 2)

""""
k_mer = Kmer_featurization(3, amino_acids_list)
data = dna_train["Sequence"][dna_train["Label"] == "1"]
array = np.mean(k_mer.obtain_kmer_feature_for_a_list_of_sequences(data), axis=0)
plt.hist(pd.DataFrame(array), bins=40, edgecolor='black')
plt.show()

"""
for aa in amino_acids_list:
    plot_feature_distribution(dna_train, "Frequency_" + aa, -0.01, 0.2)

#plot_feature_distribution(dna_train, 'Hydropathy', -3, 3)

#plot_feature_distribution(dna_train, 'Net_charge', -500, 500)

#plot_feature_distribution(dna_train, 'Molecular Weight')

#plot_feature_distribution(dna_train, 'Instability Index', -50, 180)

#plot_feature_distribution(dna_train, 'Aromaticity', -0.05, 0.4)

#plot_feature_distribution(dna_train, 'Flexibility(Mean)', 0.93, 1.06)

#plot_feature_distribution(dna_train, 'Isoelectric Point')

#plot_feature_distribution(dna_train, 'Secondary Structure Fraction', 0.1, 0.45)






####################  SET UP FEATURES FOR ALL MODELS  #################### 


# List of column names

columns = [ # 0
            'Length',
            
            # 1 - 20
            'Frequency_A', 'Frequency_R', 'Frequency_N', 'Frequency_D', 
            'Frequency_C', 'Frequency_E', 'Frequency_Q', 'Frequency_G', 
            'Frequency_H', 'Frequency_I', 'Frequency_L', 'Frequency_K', 
            'Frequency_M', 'Frequency_F', 'Frequency_P', 'Frequency_S', 
            'Frequency_T', 'Frequency_W', 'Frequency_Y', 'Frequency_V',
            
            # 21 - 22
            'Hydropathy', 
            'Net_charge',

            # 23 - 24
            'Molecular Weight',
            'Instability Index',
            
            # 25 - 26
            'Aromaticity',
            'Flexibility(Mean)',
            
            # 27 - 28
            'Isoelectric Point',
            'Secondary Structure Fraction' ]


# Choose features

columns_to_find = columns[21:]

indices = [dna_train.columns.get_loc(col) for col 
           in columns_to_find if col in dna_train.columns]

features = indices


# Set up dna_train and dna_test

train_x = dna_train.iloc[:, features]
test_x = dna_test.iloc[:, features]
train_y = dna_train.iloc[:, 0]
test_y = dna_test.iloc[:, 0]






####################  LOGISTIC REGRESSION  #################### 


# Add a constant (intercept) to the model
train_x = sm.add_constant(train_x)
test_x = sm.add_constant(test_x)


# Fit logistic regression model
logit_model = sm.Logit(train_y, train_x)
result = logit_model.fit()


# Print summary (which includes p-values)
print(result.summary())
#print(result.pvalues)


# Predict probabilities on the test set
predicted_probs = result.predict(test_x)


# Convert probabilities to class labels (0 or 1) using 0.5 threshold
predicted_classes = np.where(predicted_probs > 0.5, 1, 0)


# Confusion matrix
#conf_matrix = sketch_confusion_matrix(test_y, predicted_classes, 
                        #"Logistic Regression - Physiochemical Properties")
 

# Performance metrics
print()
print("PERFORMANCE METRICS OF LOGISTIC REGRESSION:")
print()
print(get_performance_metrics(test_y, predicted_classes))
print()
print(classification_report(test_y, predicted_classes))
print() 
plot_roc_curve(test_y, predicted_probs, "LOGISTIC REGRESSION")






####################  K NEAREST NEIGHBOURS  ####################


# Feature scaling
train_x_scaled = standardize(train_x)
test_x_scaled = standardize(test_x)


# Create the k-NN classifier with value k
k = 99
knn = KNeighborsClassifier(n_neighbors=k)


# Fit the model on the training data
knn_pred = knn.fit(train_x_scaled, train_y)


# Make predictions on the test data
knn_pred = knn.predict(test_x_scaled)
knn_scores = knn.predict_proba(test_x_scaled)[:,1]


# Confusion matrix
conf_matrix = sketch_confusion_matrix(test_y, knn_pred, "K-NN Model with k = " 
                                + str(k) + " - Physiochemical Properties")


# Performance metrics
print()
print("PERFORMANCE METRICS OF K-NN MODEL WITH K = " + str(k) + ":" )
print()
print(get_performance_metrics(test_y, knn_pred))
print()
print(classification_report(test_y, knn_pred)) 
print() 
plot_roc_curve(test_y, knn_scores, "K-NN with k = " + str(k))






####################  DECISION TREE  ####################


# Initialize and fit Decision Tree
dt_model = DecisionTreeClassifier(criterion='entropy', min_samples_split=2)
dt_model.fit(train_x, train_y)


# Predict
dt_pred = dt_model.predict(test_x)
dt_scores = dt_model.predict_proba(test_x)[:,1]


'''''
# Visualize the Decision Tree
plt.figure(figsize=(20,10)) 
plot_tree(dt_model, filled=True, feature_names=train_x.columns, 
          class_names=["0", "1"], rounded=True, fontsize=14)
plt.title("Decision Tree")
plt.show()
'''''


# Confusion matrix
#conf_matrix = sketch_confusion_matrix(test_y, dt_pred, 
      #"Decision Tree - Physiochemical Properties")


# Performance metrics
print()
print("PERFORMANCE METRICS OF DECISION TREE:")
print()
print(get_performance_metrics(test_y, dt_pred))
print()
print(classification_report(test_y, dt_pred)) 
print()
#plot_roc_curve(test_y, dt_scores, "DECISION TREE")





####################  NAIVE BAYES  ####################


# Initialize and fit Naive Bayes
nb_model = GaussianNB()
nb_model.fit(train_x, train_y)


# Predict
nb_pred = nb_model.predict(test_x)
nb_scores = nb_model.predict_proba(test_x)[:,1]


# Confusion matrix
#conf_matrix = sketch_confusion_matrix(test_y, nb_pred, 
      #"Naive Bayes - Physiochemical Properties")


# Performance metrics
print()
print("PERFORMANCE METRICS OF NAIVE BAYES:")
print()
print(get_performance_metrics(test_y, nb_pred))
print()
print(classification_report(test_y, nb_pred)) 
print()
#plot_roc_curve(test_y, nb_scores, "NAIVE BAYES")






####################  RANDOM FOREST  ####################


# Initialize and fit Random Forest
rf_model = RandomForestClassifier()
rf_model.fit(train_x, train_y)


# Performance metrics
rf_pred = rf_model.predict(test_x)
rf_scores = rf_model.predict_proba(test_x)[:,1]


# Confusion matrix
#conf_matrix = sketch_confusion_matrix(test_y, rf_pred, 
      #"Random Forest - Physiochemical Properties")


# Performance metrics
print()
print("PERFORMANCE METRICS OF RANDOM FOREST:")
print()
print(get_performance_metrics(test_y, rf_pred))
print()
print(classification_report(test_y, rf_pred)) 
print()
plot_roc_curve(test_y, rf_scores, "RANDOM FOREST")






####################  SUPPORT VECTOR MACHINE  ####################


# Feature scaling
train_x_scaled = robust_scaling(train_x)
test_x_scaled = robust_scaling(test_x)


# Initialize and fit SVM
svm_model = SVC()
svm_model.fit(train_x_scaled, train_y)


# Predict
svm_pred = svm_model.predict(test_x_scaled)


# Confusion matrix
conf_matrix = sketch_confusion_matrix(test_y, svm_pred, 
      "Support Vector Machine - Physiochemical Properties")


# Performance metrics
print()
print("PERFORMANCE METRICS OF SUPPORT VECTOR MACHINE:")
print()
print(get_performance_metrics(test_y, svm_pred))
print()
print(classification_report(test_y, svm_pred)) 
print()