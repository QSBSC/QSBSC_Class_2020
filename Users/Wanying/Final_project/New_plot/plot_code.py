import pandas as pd
import matplotlib.pyplot as plt

# ----------- Start of helper function -----------
# This fucntion plot cluster and overly another set of values (such as libraseq score) on the cluster
def overlay(cluster_df,
            x_col_name,
            y_col_name,
            overlay_df,
            overlay_val_col_name,
            cluster_ID_col_name ='Unnamed: 0',
            overlay_ID_col_name='BARCODE'):
    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    x_vals = cluster_df[x_col_name]
    y_vals = cluster_df[y_col_name]
    ax.scatter(x_vals, y_vals, marker=".", color="0.9")

    overlay_IDs = overlay_df[overlay_ID_col_name]
    cluster_mask = (cluster_df[cluster_ID_col_name].isin(overlay_IDs))
    overlay_in_cluster_df = cluster_df[cluster_mask]
    
    overlay_mask = (overlay_df[overlay_ID_col_name].isin(overlay_in_cluster_df[cluster_ID_col_name]))
    vals_to_plot_df = overlay_df[overlay_mask]

    vals = list(vals_to_plot_df[overlay_val_col_name])
    norm = mpl.colors.Normalize(vmin=sorted(vals)[0], vmax=sorted(vals)[-1])
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Blues)
    cmap.set_array([])

    cluster_dim1 = list(overlay_in_cluster_df[x_col_name])
    cluster_dim2 = list(overlay_in_cluster_df[y_col_name])

    for i in range(len(cluster_dim1)):
        ax.scatter(cluster_dim1[i], cluster_dim2[i], marker=".", color=cmap.to_rgba(vals[i]), alpha=0.6)

    ax.set_title(overlay_val_col_name + " libraseq score on RNA seq clusters")
    ax.set_xlabel(x_col_name)
    ax.set_ylabel(y_col_name)
    fig.colorbar(cmap)
    return(fig, ax)
    
# ----------- End of helper function -----------


tsne_df = pd.read_csv("reduction_tsne.csv")

# contruct libraseq score dataframe for ech antigen
libraScore_file = "/home/libraseq/wennie/3602/AG_libraseq_scores.tab"
libraScore_df = pd.read_csv(libraScore_file, sep="\t")

libraScore_dic["HIV_BG505_N332T_df"] = libraScore_df.loc[:,("BARCODE","HIV_BG505_N332T")]

HIV_BG505_N332T_fig, HIV_BG505_N332T_ax = overlay(cluster_df = tsne_df,
                                  x_col_name = "tSNE_1",
                                  y_col_name = "tSNE_2",
                                  overlay_df = libraScore_dic["HIV_BG505_N332T_df"],
                                  overlay_val_col_name = "HIV_BG505_N332T")

HIV_BG505_N332T_fig.savefig(fname = "HIV_BG505_N332T libraseq score on RNA seq clusters.jpg")