import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as net
import lammps_logfile
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
from itertools import count
from matplotlib.pyplot import cm

# Plot functions
from scipy import stats
from scipy.stats import gaussian_kde

def set_size(width, fraction=1):
    """ Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float
            Document textwidth or columnwidth in pts
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio

    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim


def freedman_diaconis(data, returnas="width"):
    """
    Use Freedman Diaconis rule to compute optimal histogram bin width. 
    ``returnas`` can be one of "width" or "bins", indicating whether
    the bin width or number of bins should be returned respectively. 


    Parameters
    ----------
    data: np.ndarray
        One-dimensional array.

    returnas: {"width", "bins"}
        If "width", return the estimated width for each histogram bin. 
        If "bins", return the number of bins suggested by rule.
    """
    
    data = np.asarray(data, dtype=np.float_)
    IQR = stats.iqr(data, rng=(25, 75), scale=1.0, nan_policy="omit")
    if(IQR == 0 or np.isinf(data.max())):
        # int(freedman_diaconis(data[data!=0], returnas="bins")/(np.sum(data!= 0)/len(data)))
        NBRBINS = -1
        bw = -1
    else:
        N = data.size
        bw = (2 * IQR) / np.power(N, 1/3)
        datmin, datmax = data.min(), data.max()
        datrng = datmax - datmin
        NBRBINS = int((datrng / bw) + 1)
    return(NBRBINS, bw)

def get_hist(array):
    figsize=set_size(500.484)
    fig, axs = plt.subplots(1, 1, figsize=(8,8))
    NBR_BINS, bw = freedman_diaconis(data=array, returnas="bins")
    Prob_dist, bin_edges, _ = plt.hist(array, bins=NBR_BINS, 
                                       histtype='bar',density=True, 
                                       fill=False, cumulative=False, 
                                       log=False)
    # plt.clf()
    print(np.dot(Prob_dist,np.diff(bin_edges)))
    assert abs(np.dot(Prob_dist,np.diff(bin_edges))-1) < 1e-5
    
    bin_avgs = (bin_edges[:-1]+bin_edges[1:])/2
    y = Prob_dist[Prob_dist!=0]
    x = bin_avgs[Prob_dist!=0]
    
    return x, y 

def scale01(M):

    New_M = np.zeros((M.shape))
    MIN = np.min(M)
    MAX = np.max(M)
    for i in range(M.shape[0]):
        M_ = M[i]
        if (MAX == MIN):
            New_M[i] = 0.0 * M_
        else:
            New_M[i] = (M_ - MIN) / (MAX - MIN)

    return New_M


def Plot_force_chains(Grain, Criteria, OUTPUT_DIR):
    par_type = Grain[:, 1]
    x = Grain[:, 3]
    y = Grain[:, 4]
    R = Grain[:, 2]
    #h = 1*(Criteria > 1)
    h = np.array([i if i >= 1 else 0 for i in Criteria])
    s = scale01(h)
    colors = [matplotlib.cm.gray_r(color) for color in s]

    fig, ax1 = plt.subplots(1, 1, figsize=(15, 15))
    for xi, yi, r, color, t, lh in zip(x, y, R, colors, par_type, h):
        if(lh):
            circle = plt.Circle((xi, yi), r, color=color)
            ax1.add_patch(circle)

    ax1.set_aspect('equal')
    ax1.set_xlim(min(x)-.5, max(x)+.5)
    ax1.set_ylim(min(y)-.5, max(y)+.5)
    plt.axis('off')
    #sc = ax1.scatter(x, y,c=h, s=20);
    sc = ax1.scatter(x,
                     y,
                     s=0,
                     c=h,
                     cmap='gray_r',
                     vmin=min(h),
                     vmax=max(h),
                     facecolors='none')

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(sc, cax=cax)
    cbar.set_label('C', rotation=270, labelpad=10)
    cbar.ax.tick_params(labelsize=25)
    plt.savefig(OUTPUT_DIR + "Force_chains" + '.png',
                format='png',
                bbox_inches='tight')
    plt.close(fig)

def plot_force_chains_network(G, OUTPUT_DIR):

    fig, ax1 = plt.subplots(1, 1, figsize=(15, 12))
    # edges = G.edges()
    # weights = [G[u][v]['weight'] for u, v in edges]
    pos = net.get_node_attributes(G, 'pos')
    new_cmap = plt.get_cmap('jet')
    # net.draw(Force_chains,pos,edge_color='black',node_size=1,node_color='blue',edge_cmap=new_cmap)
    # net.draw(G, pos, edge_color=weights,
    #         width=2,
    #          node_size=2, node_color='blue', edge_cmap=new_cmap)

    weights = set(net.get_node_attributes(G, 'weight').values())
    mapping = dict(zip(sorted(weights), count()))
    nodes = G.nodes
    colors = [mapping[G.nodes[n]['weight']] for n in nodes]
    net.draw(G, pos,
             edge_color='black', node_size=20,
             node_color=colors, cmap=new_cmap)

    vmin = min(weights)
    vmax = max(weights)

    sm = plt.cm.ScalarMappable(
        cmap=new_cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_label('C', rotation=270, labelpad=10)
    cbar.ax.tick_params(labelsize=25)
    # plt.colorbar(sm, shrink=0.9)
    plt.savefig(OUTPUT_DIR + "Network" + '.png',
                format='png',
                bbox_inches='tight')
    plt.close(fig)


def plot_subgraph(sG_L, label, OUTPUT_DIR):

    # Plot the node community
    n = len(sG_L)
    color = cm.rainbow(np.linspace(0, 1, n))
    fig, ax = plt.subplots(1, 1, figsize=set_size(500.484))
    for index, sg in enumerate(sG_L):
        pos = net.get_node_attributes(sg, 'pos')
        net.draw(sg,
                 pos=pos,
                 node_size=0,
                 edge_color=color[index]
                 )
    plt.savefig(OUTPUT_DIR + label + '.png',
                format='png',
                bbox_inches='tight')
    plt.close(fig)


def plot_node_Community(Graph, community_label, OUTPUT_DIR):
    # Plot the node community
    fig, ax = plt.subplots(1, 1, figsize=set_size(500.484))
    edges = Graph.edges()
    #weights = [Graph[u][v]['weight'] for u, v in edges]
    pos = net.get_node_attributes(Graph, 'pos')
    groups = set(net.get_node_attributes(Graph, community_label).values())
    mapping = dict(zip(sorted(groups), count()))
    nodes = Graph.nodes
    colors = [mapping[Graph.nodes[n][community_label]] for n in nodes]
    net.draw(Graph, pos, edge_color='black', node_size=20,
             node_color=colors, cmap=plt.cm.hot)
    plt.savefig(OUTPUT_DIR + "community_label" + '.pngs',
                format='png',
                bbox_inches='tight')
    plt.close(fig)

# def plot_force_chains_network(G, OUTPUT_DIR):

#     fig, ax1 = plt.subplots(1, 1, figsize=(15, 15))
#     edges = G.edges()
#     weights = [G[u][v]['weight'] for u, v in edges]
#     pos = net.get_node_attributes(G, 'pos')
#     new_cmap = ax1.get_cmap('jet')
#     # net.draw(Force_chains,pos,edge_color='black',node_size=1,node_color='blue',edge_cmap=new_cmap)
#     net.draw(G, pos, edge_color=weights,
#             node_size=1, node_color='blue', edge_cmap=new_cmap)
#     vmin = min(weights)
#     vmax = max(weights)
#     sm = ax1.cm.ScalarMappable(
#         cmap=new_cmap, norm=ax1.Normalize(vmin=vmin, vmax=vmax))
#     sm._A = []
#     ax1.colorbar(sm, shrink=0.9)
#     plt.savefig(OUTPUT_DIR + "Netowrk" + '.eps',
#                 format='eps',
#                 dpi=1000,
#                 bbox_inches='tight')
#     plt.close(fig)
