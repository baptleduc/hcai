import sys
import networkx as nx
import random
import matplotlib.pyplot as plt
from copy import deepcopy

# Simulation parameters
K_VALUES = [0, 5, 10, 20, 30, 40, 50, 100] 
BETA = 0.04   # infection rate
R = 0.2       # recovery rate
NUM_RUNS = 500  # number of simulation runs per K (to average results)
MAX_TIME_STEPS = 10000  # maximum time steps for simulation


def vaccinate_nodes(nodes_to_vaccinate, graph):
    """Remove specified nodes from the graph (simulate vaccination)."""
    graph_copy = deepcopy(graph)
    graph_copy.remove_nodes_from(nodes_to_vaccinate)
    return graph_copy


def high_degree_selection(graph, K, initial_infected):
    """Select the K nodes with highest out-degree centrality, excluding initial_infected."""
    if K <= 0:
        return []
    degree_centrality = nx.out_degree_centrality(graph)
    candidates = [n for n in graph.nodes() if n != initial_infected]
    sorted_nodes = sorted(candidates, key=lambda n: degree_centrality[n], reverse=True)
    return sorted_nodes[:K]



def random_selection(graph, K, initial_infected):
    """Select K random nodes excluding the initial_infected."""
    candidates = [n for n in graph.nodes() if n != initial_infected]
    return random.sample(candidates, k=K)


def launch_simulation(graph, initial_infected):
    """Run SIR simulation on the graph and return epidemic duration and final infected count."""
    inf = {initial_infected: 1}
    num_r = 0

    for t in range(1, MAX_TIME_STEPS):
        inf_list = list(inf.keys())
        random.shuffle(inf_list)

        for v in inf_list:
            for nbr in list(graph[v]):
                if nbr not in inf and random.random() < BETA:
                    inf[nbr] = 1

            if random.random() < R:
                del inf[v]
                graph.remove_node(v)
                num_r += 1
        num_inf = num_r + len(inf)

        if len(inf) == 0:
            return t, num_inf  # duration, total infected (I+R)

    return MAX_TIME_STEPS, num_inf


def compare_strategies(graph, K_values):
    hd_durations, hd_totals = [], []
    rand_durations, rand_totals = [], []

    for K in K_values:
        hd_runs_dur, hd_runs_total = [], []
        rand_runs_dur, rand_runs_total = [], []

        for _ in range(NUM_RUNS):
            g_base = deepcopy(graph)

            # common initial infected node
            initial_infected = random.choice(list(g_base.nodes()))

            # --- High-degree ---
            g_hd = vaccinate_nodes(high_degree_selection(g_base, K, initial_infected), g_base)
            dur_hd, tot_hd = launch_simulation(deepcopy(g_hd), initial_infected)
            hd_runs_dur.append(dur_hd)
            hd_runs_total.append(tot_hd)

            # --- Random ---
            g_base2 = deepcopy(graph)
            g_rand = vaccinate_nodes(random_selection(g_base2, K, initial_infected), g_base2)
            dur_rand, tot_rand = launch_simulation(deepcopy(g_rand), initial_infected)
            rand_runs_dur.append(dur_rand)
            rand_runs_total.append(tot_rand)

        # Average results over runs
        hd_durations.append(sum(hd_runs_dur) / NUM_RUNS)
        hd_totals.append(sum(hd_runs_total) / NUM_RUNS)
        rand_durations.append(sum(rand_runs_dur) / NUM_RUNS)
        rand_totals.append(sum(rand_runs_total) / NUM_RUNS)

    return hd_durations, hd_totals, rand_durations, rand_totals



def plot_results(K_values, hd_values, rand_values, ylabel, title, save_path=None):
    plt.figure(figsize=(9, 6))
    plt.plot(K_values, hd_values, marker="o", label="High-degree vaccination")
    plt.plot(K_values, rand_values, marker="s", label="Random vaccination")
    plt.xlabel("K (number of vaccinated nodes)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid()
    if save_path:
        plt.savefig(save_path)
        plt.close()
    else: 
        plt.show()


if __name__ == "__main__":
    network = sys.argv[1]
    fig_output_path = sys.argv[2]
    graph = nx.read_edgelist(network, create_using=nx.DiGraph)

    
    

    print("Evaluating both strategies fairly...")
    hd_durations, hd_totals, rand_durations, rand_totals = compare_strategies(graph, K_VALUES)

    # Plot epidemic duration
    plot_results(K_VALUES, hd_durations, rand_durations,
                 ylabel="Epidemic duration (timesteps)",
                 title="Epidemic duration vs K",
                 save_path=fig_output_path + "_duration.png")

    # Plot total infected
    plot_results(K_VALUES, hd_totals, rand_totals,
                 ylabel="Total infected (I+R)",
                 title="Total infected vs K",
                 save_path=fig_output_path + "_total_infected.png")
