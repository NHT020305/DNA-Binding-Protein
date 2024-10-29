import matplotlib.pyplot as plt



####################  AMINO ACIDS DISTRIBUTION  ####################


#max_index = data['Length'].idxmax()
#amino_acid_list = list(set(data["Sequence"][max_index]))


def plot_amino_acids_distribution(data, binding_state, amino_acids_list):

    column_list = amino_acids_list

    fig, axes = plt.subplots(nrows=5, ncols=4) 
    axes = axes.flatten()

    title = "Distributions of Amino Acids"

    for i, column in enumerate(column_list):
        used_data = data[data["Label"] == binding_state]["Frequency_" + column]
        counts, bins, _  = axes[i].hist(used_data, bins=20, color='blue', edgecolor='black') 
        axes[i].set_xlabel('Frequency')
        axes[i].set_ylabel('Count')
        axes[i].set_title(column)
        #max_y_value = counts.max()
        #axes[i].set_ylim(0, max_y_value + 2000)
        #axes[i].set_xlim(0, 0.2)

    if binding_state:
        plt.suptitle(title + " in DNA-Binding Sequences")
    else:
        plt.suptitle(title + " in non DNA-Binding Sequences")

    plt.subplots_adjust(wspace=0.5, hspace=2.5)
    plt.show()