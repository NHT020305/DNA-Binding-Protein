import matplotlib.pyplot as plt


def plot_feature_distribution(data, feature, start=None, end=None):
    
    column_list = {1: "DNA-binding Sequences", 0: "non DNA-binding Sequences"}

    fig, axes = plt.subplots(nrows=1, ncols=2) 
    axes  = axes.flatten()

    for i, column in enumerate(column_list):
        used_data = data[data["Label"] == column][feature]
        axes[i].boxplot(used_data, patch_artist=True, 
                        boxprops=dict(facecolor='lightblue', 
                                      edgecolor='black'))
        axes[i].set_xlabel(feature)
        axes[i].set_title(feature + " of " + column_list[column])
        axes[i].set_ylim(start, end)

    plt.show()  