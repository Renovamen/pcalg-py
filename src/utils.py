import networkx as nx
import matplotlib.pyplot as plt
from typing import List, Optional

chains = []

def dfs(graph, k: int, chain: List[int], visit: List[bool]):
    """因果关系链深搜"""
    flag = False

    for i in range(len(graph)):
        if graph[i][k] == 1 and not visit[i] :
            flag =True
            visit[i] = True
            chain.append(i)

            dfs(graph, i, chain, visit)

            chain.pop()
            visit[i] = False

    if not flag:
        chains.append(chain.copy())

def get_causal_chains(graph, start: int, labels: List[str]):
    global chains
    chains = []

    visit = [False for _ in range(len(labels))]
    visit[start] = True

    chain = [start]

    dfs(graph, start, chain, visit)

    return "\n".join([" <- ".join(list(map(lambda x: labels[x] + f" ({x})", c))) for c in chains])

def plot(graph, labels: List[str], path: str = Optional[None]):
    """可视化学习出的贝叶斯网络"""
    G = nx.DiGraph()  # 创建空有向图

    for i in range(len(graph)):
        G.add_node(labels[i])
        for j in range(len(graph[i])):
            if graph[i][j] == 1:
                G.add_edge(labels[i], labels[j])

    nx.draw(G, with_labels=True)

    if path:
        plt.savefig(path)

    plt.show()
