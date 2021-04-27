import networkx as nx
import matplotlib.pyplot as plt

def dfs(graph, k: int, path: list, vis: list):
    """因果关系链深搜"""
    flag = False

    for i in range(len(graph)):
        if (graph[i][k] == 1) and (vis[i] != True) :
            flag =True
            vis[i] = True
            path.append(i)
            dfs(graph, i, path, vis)
            path.pop()
            vis[i] = False

    if flag == False:
        print(path)

def plot(graph, labels: list, path: str):
    """可视化学习出的贝叶斯网络"""
    G = nx.DiGraph()  # 创建空有向图

    for i in range(len(graph)):
        G.add_node(labels[i])
        for j in range(len(graph[i])):
            if graph[i][j] == 1:
                G.add_edges_from([(labels[i], labels[j])])

    nx.draw(G, with_labels=True)
    plt.savefig(path)
    plt.show()
