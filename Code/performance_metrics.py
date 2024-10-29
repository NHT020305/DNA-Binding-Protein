import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, classification_report
from sklearn.metrics import roc_curve, auc, matthews_corrcoef






####################  CONFUSION MATRIX SKETCHING FUNCTION  ####################

def sketch_confusion_matrix(test_set, predicted_set, name):
    
    conf_matrix = confusion_matrix(test_set, predicted_set)
    
    plt.figure(figsize=(8, 6))
    
    sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', 
            xticklabels=['Non DNA-Binding', 'DNA-Binding'], 
            yticklabels=['Non DNA-Binding', 'DNA-Binding'])
    
    plt.title('Confusion Matrix of ' + name)
    
    plt.xlabel('Predicted Label')
    
    plt.ylabel('True Label')
    
    plt.show()
    
    return conf_matrix






####################  PERFORMANCE METRICS FUNCTION  #################### 

def get_performance_metrics(test, pred):
    
    # Calculate metrics
    accuracy = accuracy_score(test, pred)
    precision = precision_score(test, pred)
    recall = recall_score(test, pred)
    f1 = f1_score(test, pred)
    mcc = matthews_corrcoef(test, pred)
    
    tn, fp, fn, tp = confusion_matrix(test, pred).ravel()

    # Calculate specificity
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    
    return {
        "Accuracy": round(accuracy, 2),
        "Precision": round(precision, 2),
        "Recall(Sensitivity)": round(recall, 2),
        "F1 Score": round(f1, 2),
        "Specificity": round(specificity, 2),
        "Matthews Correlation Coefficient(MCC)": round(mcc, 2)
    }






####################  AUC - ROC  ####################

def plot_roc_curve(y_test, y_probs, name):
    
    # Step 5: Compute ROC curve and AUC
    fpr, tpr, thresholds = roc_curve(y_test, y_probs)
    roc_auc = auc(fpr, tpr)

    # Step 6: Plot ROC curve
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--') # Diagonal line (random guess)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('RECEIVER OPERATING CHARACTERISTIC of '+ name)
    plt.legend(loc="lower right")
    plt.show()