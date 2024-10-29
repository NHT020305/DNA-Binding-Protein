import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



class Kmer_featurization:


  def __init__(self, k, letters):
    
    """
    seqs: a list of DNA sequences
    k: the "k" in k-mer
    """
    
    self.k = k
    self.letters = letters
    self.multiplyBy = len(letters) ** np.arange(k-1, -1, -1) 
    # the multiplying number for each digit position in the k-number system
    self.n = len(letters) ** k # number of possible k-mers


  def obtain_kmer_feature_for_a_list_of_sequences(self, seqs, 
                          write_number_of_occurrences=False):
    
    """
    Given a list of m DNA sequences, 
    return a 2-d array with shape (m, 4**k) for the 1-hot representation of the kmer features.

    Args:
      write_number_of_occurrences:
        a boolean. If False, then in the 1-hot representation, 
        the percentage of the occurrence of a kmer will be recorded; 
        otherwise the number of occurrences will be recorded. Default False.    
    """
    
    kmer_features = []
    for seq in seqs:
      this_kmer_feature = self.obtain_kmer_feature_for_one_sequence(seq.upper(), 
                          write_number_of_occurrences=write_number_of_occurrences)
      kmer_features.append(this_kmer_feature)

    kmer_features = np.array(kmer_features)

    return kmer_features


  def obtain_kmer_feature_for_one_sequence(self, seq, 
                            write_number_of_occurrences=False):
    
    """
    Given a DNA sequence, return the 1-hot representation of its kmer feature.

    Args:
      seq: 
        a string, a DNA sequence
      write_number_of_occurrences:
        a boolean. If False, then in the 1-hot representation, 
        the percentage of the occurrence of a kmer will be recorded; 
        otherwise the number of occurrences will be recorded. Default False.
    """
    
    number_of_kmers = len(seq) - self.k + 1

    kmer_feature = np.zeros(self.n)

    for i in range(number_of_kmers):
      this_kmer = seq[i:(i+self.k)]
      this_numbering = self.kmer_numbering_for_one_kmer(this_kmer)
      kmer_feature[this_numbering] += 1

    if not write_number_of_occurrences:
      kmer_feature = kmer_feature / number_of_kmers

    return kmer_feature


  def kmer_numbering_for_one_kmer(self, kmer):
    
    """
    Given a k-mer, return its numbering (the 0-based position in 1-hot representation)
    """
    
    digits = []
    for letter in kmer:
      digits.append(self.letters.index(letter))

    digits = np.array(digits)

    numbering = (digits * self.multiplyBy).sum()

    return numbering



def plot_kmers_distribution(data, start, space, amino_acids_list, k):
    
  k_mer = Kmer_featurization(k, amino_acids_list)

  column_list = {"1": "DNA-binding Sequences", "0": "non DNA-binding Sequences"}

  fig, axes = plt.subplots(nrows=1, ncols=2) 
  axes = axes.flatten()

  for i, column in enumerate(column_list):
    used_data = k_mer.obtain_kmer_feature_for_a_list_of_sequences(
            data["Sequence"][data["Label"] == column][start : start + space, ])
    used_data = pd.DataFrame(np.transpose(used_data))
    #used_data = np.mean(used_data, axis=0)
    axes[i].hist(used_data) 
    axes[i].set_xlabel('Normalized frequency of the k-mers')
    axes[i].set_title("K-mer frequency distribution " + column_list[column])
    axes[i].set_xlim(0, 0.03)

  plt.show()