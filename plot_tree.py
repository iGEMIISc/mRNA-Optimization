import networkx as nx
import numpy as np

# imports are only for networkx plotting in real-time (not important in sherlock)
try:
    from networkx.drawing.nx_agraph import graphviz_layout
    import matplotlib.pyplot as plt
except ImportError:
    pass


def plot_add_children(G, parent, parent_label):
    children = parent.get_children()
    for i, child in enumerate(children):
        #child_label = parent_label + '.' + str(i+1)
        child_label = child.get_id()
        if child.check_solution_found():
            type_color = 'blue'
        elif child.get_dead():
            type_color = 'green'
        else:
            type_color = 'red'
        score = child.get_MC_result()[-1]
        method = score[-1].split(':')
        method = method[0][0] + ':' + method[1].replace('->','-')
        G.add_node(child_label, color=type_color, method=method, score=score, UTC1='{:.1f}'.format(child.get_UCT1()))
        G.add_edge(parent_label, child_label)
        plot_add_children(G, child, child_label)


def plot_pickle_save(root, filename):
    G = nx.DiGraph()
    #parent_label = 'R'
    parent_label = root.get_id()
    score = root.get_MC_result()[-1]
    G.add_node(parent_label, color='red', method='Root', score=score, UTC1=np.NaN)
    plot_add_children(G, root, parent_label)
    nx.write_gpickle(G, filename)


def plot_pickle_load(filename, plot = True):
    G = nx.read_gpickle(filename)
    # for node, attr in G.nodes(data=True):
    #     if attr['color'] == 'green':
    #         G.remove_node(node)

    pos=graphviz_layout(G, prog='dot')
    color_list = []
    for node, attr in G.nodes(data=True):
        temp_color = attr['color']
        color_list.append(temp_color)
        if temp_color == 'blue':
            print(node)

    if plot:
        method_d = {}
        for node, attr in G.nodes(data=True):
            method_d[node] = attr['method'] 
        nx.draw(G, pos, labels=method_d, node_color=color_list, with_labels=True, arrows=False)
        #plt.savefig('nx_test.png')
        plt.show()
    return len(color_list), 'blue' in color_list


def plot_MC(root):
    G = nx.DiGraph()
    parent_label = 'R'
    score = root.get_MC_result()[-1]
    G.add_node(parent_label, color='red', method='Root', score=score, UTC1=np.NaN)
    plot_add_children(G, root, parent_label)
    pos=graphviz_layout(G, prog='dot')
    color_list = [attr['color'] for node, attr in G.nodes(data=True)]
    method_d = {}
    for node, attr in G.nodes(data=True):
        method_d[node] = attr['method'] 
    nx.draw(G, pos, labels=method_d, node_color=color_list, with_labels=True, arrows=False)
    plt.show()
