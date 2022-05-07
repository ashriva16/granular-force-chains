import networkx.algorithms.community as nx_comm
from networkx.algorithms.community import greedy_modularity_communities
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as net
import lammps_logfile
import plots as PLT

def read_pairwise(dir_,STEP):

    file1 = open(dir_+'neigh.deform', 'r')
    lines = [line.rstrip() for line in file1]

    for num, line in enumerate(lines, 0):
        if 'ITEM' in line:
            if 'ITEM: TIMESTEP' in line:
                curr_timestep = int(lines[num+1])
            if 'ITEM: NUMBER OF ENTRIES' in line:
                num_entries = int(lines[num+1])
                table = []
            if 'ITEM: ENTRIES' in line:
                if(curr_timestep != STEP):
                    continue
                x = lines[num+1:num+1+num_entries]
                for i in range(len(x)):
                    table.append(list(map(float, x[i].split(" "))))

                if(curr_timestep == STEP):
                    break

    print("read_pairwise:", curr_timestep, "STEP:", STEP, flush=True)
    pairwise = np.array(table.copy())

    return pairwise


def read_grain_data(dir_,STEP):
    file2 = open(dir_+'visualize.deform', 'r')
    lines = [line.rstrip() for line in file2]

    for num, line in enumerate(lines, 0):
        if 'ITEM' in line:
            if 'ITEM: TIMESTEP' in line:
                curr_timestep = int(lines[num+1])
            if 'ITEM: NUMBER OF ATOMS' in line:
                num_entries = int(lines[num+1])
                table = []
            if 'ITEM: ATOMS' in line:
                if(curr_timestep != STEP):
                    continue
                x = lines[num+1:num+1+num_entries]
                for i in range(len(x)):
                    table.append(list(map(float, x[i].split(" "))))

                if(curr_timestep == STEP):
                    break

    print("read_grain_data:", curr_timestep, "STEP:", STEP, flush=True)

    grains = np.array(table.copy())
    grains = grains[grains[:, 1] != 3]
    grains = grains[grains[:, 0].argsort()]

    return grains


def create_chain_from_stress(Grain, Pairwise, Criteria):
    # Calculate force between particles
    fx = Pairwise[:, 5]
    fy = Pairwise[:, 6]
    fn = np.sqrt(fx**2 + fy**2)
    fN1 = fn / np.mean(fn)

    Force_chains = net.Graph()

    nx, ny = Grain.shape
    for node, x, y, h in zip(Grain[:, 0], Grain[:, 3], Grain[:, 4], Criteria):
        if(h >= 1):
            Force_chains.add_node(node, pos=(x, y), weight=h)

    nx, ny = Pairwise.shape
    for i in range(nx):
        id1 = int(Pairwise[i, 1])-1
        id2 = int(Pairwise[i, 2])-1
        dist = np.sqrt((Grain[id1, 3] - Grain[id2, 3])**2 +
                    (Grain[id1, 4] - Grain[id2, 4])**2)
        if (dist <= 5.0 and Force_chains.has_node(id1+1) and Force_chains.has_node(id2+1) and fN1[i] >= 1):
            Force_chains.add_edge(id1+1, id2+1, weight=fN1[i]/10)

    return Force_chains


def get_subgraph_from_community(G, community):
    
    sG_L = []
    for index, comm in enumerate(community):
        sG = net.Graph.copy(G)
        for node in G:
            if(sG.nodes[node]['community'] != index):
                sG.remove_node(node)
        sG_L.append(net.Graph.copy(sG))

    return sG_L


def get_independent_connected_components(G):
    # Get the number of independent components
    # https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.components.number_connected_components.html
    sG = [G.subgraph(c) for c in net.connected_components(G)]
    n = len(sG)

    return sG

def Community_detection(G1, best_res):

    G3 = net.Graph.copy(G1)
    m = G3.number_of_edges()
    for (u, v, d) in G3.edges(data=True):
        G3[u][v]['weight'] = G3[u][v]['weight'] + (G3.degree[u]*G3.degree[v])/(2*m) - best_res*1

    L_C = greedy_modularity_communities(G3,
                                        weight='weight')

    mod = nx_comm.modularity(G3, L_C, weight='weight')

    # print(best_res, mod, len(L_C))
    net.set_node_attributes(G3, 1, 'community')
    for index, comm in enumerate(L_C):
        for node in comm:
            G3.nodes[node]['community'] = index

    return get_subgraph_from_community(G3, L_C)

def posprocess_force_chains(G):
    G1 = net.Graph.copy(G)
    sG = [G1.subgraph(c) for c in net.connected_components(G1)]
    n = len(sG)

    for index, sg in enumerate(sG):
        e = list(sg.nodes)
        if(len(e) <= 25):
            G1.remove_nodes_from(e)

    return G1

def height_graph(G):
    max_distance = 0
    for u in G.nodes:
        for v in G.nodes:
            dist = abs(G.nodes[u]['pos'][1]-G.nodes[v]['pos'][1])
            if(dist > max_distance):
                max_distance = dist

    return max_distance

def dfs_traverse(G):
    chain_list = list(net.dfs_edges(G))

    Graph = net.Graph()
    for index, edge in enumerate(chain_list):
        u = edge[0]
        v = edge[1]
        Graph.add_edge(u, v, weight=G[u][v]['weight'])

    net.set_node_attributes(Graph, net.get_node_attributes(G, 'pos'), 'pos')
    net.set_node_attributes(
        Graph, net.get_node_attributes(G, 'weight'), 'weight')

    return Graph


def bfs_traverse(G):
    chains = net.algorithms.chains.chain_decomposition(G)
    chain_list = list(net.bfs_edges(G))

    Graph = net.Graph()
    for index, edge in enumerate(chain_list):
        u = edge[0]
        v = edge[1]
        Graph.add_edge(u, v, weight=G[u][v]['weight'])

    net.set_node_attributes(Graph, net.get_node_attributes(G, 'pos'), 'pos')
    net.set_node_attributes(
        Graph, net.get_node_attributes(G, 'weight'), 'weight')

    return Graph


def force_chain_analysis(dir_):
    log2 = lammps_logfile.File(dir_+"log.deform")
    STEP = int(log2.data_dict['Step'][-1])

    grain0 = read_grain_data(dir_, 0)
    grains = read_grain_data(dir_, STEP)
    pairwise0 = read_pairwise(dir_, 0)
    pairwise = read_pairwise(dir_, STEP)

    index = 15
    Criteria = -grains[:, index]/np.percentile(-grains[:, index], 60)
    PLT.Plot_force_chains(grains, Criteria, dir_)

    Force_chains = create_chain_from_stress(grains, pairwise, Criteria)
    G = posprocess_force_chains(Force_chains)

    number_of_nodes = G.number_of_nodes()
    number_of_edges = G.number_of_edges()
    density = net.density(G)
    # diameter = net.diameter(G)

    G_deg = net.degree_histogram(G)
    G_deg_sum = [a * b for a, b in zip(G_deg, range(0, len(G_deg)))]
    avg_degree = sum(G_deg_sum) / G.number_of_nodes()

    PLT.plot_force_chains_network(G, dir_)

    Independent_Subgraph = get_independent_connected_components(G)
    PLT.plot_subgraph(Independent_Subgraph, "Independent_Subgraph", dir_)

    number_of_Independent_components = len(Independent_Subgraph)
    Independent_component_dia = np.mean(
        [net.diameter(sg) for sg in Independent_Subgraph])
    Independent_component_ht = np.mean(
        [height_graph(sg) for sg in Independent_Subgraph])

    # ------------------------------------------------------------
    Community_Subgraph = Community_detection(G, 0.2)
    PLT.plot_subgraph(Community_Subgraph, "Community_Subgraph", dir_)

    Community_Subgraph_dia = np.mean(
        [net.diameter(sg) for sg in Community_Subgraph])

    Community_Subgraph_ht = np.mean(
        [height_graph(sg) for sg in Community_Subgraph])

    # ------------------------------------------------------------
    dfs_G = dfs_traverse(G)
    PLT.plot_force_chains_network(dfs_G, dir_+"dfs_")
    density_dfs = net.density(dfs_G)
    # ------------------------------------------------------------
    # Note working
    # bfs_G = bfs_traverse(G)
    # PLT.plot_force_chains_network(dfs_G, dir_+"bfs_")

    return [number_of_nodes, number_of_edges, density, avg_degree, 
            number_of_Independent_components, 
            Independent_component_dia,
            Independent_component_ht, 
            Community_Subgraph_dia, 
            Community_Subgraph_ht,
            density_dfs]

if __name__ == '__main__':
    dir_ = "../LAMMPS/uniaxial/results/Hooke_H-wall_V-periodic_Particle-indent/Hooke_H-wall_V-periodic_Particle-indent/"
    # dir_ = dir_ + "1_0.01_0.0_1.0/"
    # dir_ = dir_ + "2_0.1_0.0_1.0/"
    dir_ = dir_ + "1_0.0001_0_1.0/"
    
    n1, n2, n3, n4, n5, n6, n7, n8, n9 ,n10, n11 = force_chain_analysis(dir_)
