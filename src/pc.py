from itertools import combinations
from scipy.stats import norm
import numpy as np
import math
from typing import List

def get_neighbors(G, x: int, y: int):
    return [i for i in range(len(G)) if G[x][i] == True and i != y]

def gauss_ci_test(suff_stat, x: int, y: int, K: List[int], cut_at: float = 0.9999999):
    """条件独立性检验"""
    C = suff_stat["C"]
    n = suff_stat["n"]

    # ------ 偏相关系数 ------
    if len(K) == 0:  # K 为空
        r = C[x, y]

    elif len(K) == 1:  # K 中只有一个点，即一阶偏相关系数
        k = K[0]
        r = (C[x, y] - C[x, k] * C[y, k]) / math.sqrt((1 - C[y, k] ** 2) * (1 - C[x, k] ** 2))

    else:  # 其实我没太明白这里是怎么求的，但 R 语言的 pcalg 包就是这样写的
        m = C[np.ix_([x] + [y] + K, [x] + [y] + K)]
        p = np.linalg.pinv(m)
        r = -p[0, 1] / math.sqrt(abs(p[0, 0] * p[1, 1]))

    r = min(cut_at, max(-cut_at, r))

    # Fisher's z-transform
    z = 0.5 * math.log1p((2 * r) / (1 - r))
    z_standard = z * math.sqrt(n - len(K) - 3)

    # Φ^{-1}(1-α/2)
    p_value = 2 * (1 - norm.cdf(abs(z_standard)))

    return p_value

def skeleton(suff_stat, alpha: float):
    n_nodes = suff_stat["C"].shape[0]

    # 分离集
    O = [[[] for _ in range(n_nodes)] for _ in range(n_nodes)]

    # 完全无向图
    G = [[i != j for i in range(n_nodes)] for j in range(n_nodes)]

    # 点对（不包括 i -- i）
    pairs = [(i, (n_nodes - j - 1)) for i in range(n_nodes) for j in range(n_nodes - i - 1)]

    done = False
    l = 0  # 节点数为 l 的子集

    while done != True and any(G):
        done = True

        # 遍历每个相邻点对
        for x, y in pairs:
            if G[x][y] == True:
                neighbors = get_neighbors(G, x, y)  # adj(C,x) \ {y}

                if len(neighbors) >= l:  # |adj(C, x) \ {y}| > l
                    done = False

                    # |adj(C, x) \ {y}| = l
                    for K in set(combinations(neighbors, l)):
                        # 节点 x, y 是否被节点数为 l 的子集 K d-seperation
                        # 条件独立性检验，返回 p-value
                        p_value = gauss_ci_test(suff_stat, x, y, list(K))

                        # 条件独立
                        if p_value >= alpha:
                            G[x][y] = G[y][x] = False  # 去掉边 x -- y
                            O[x][y] = O[y][x] = list(K)  # 把 K 加入分离集 O
                            break

        l += 1

    return np.asarray(G, dtype=int), O

def extend_cpdag(G, O):
    n_nodes = G.shape[0]

    def rule1(g):
        """Rule 1: 如果存在链 i -> j - k ，且 i, k 不相邻，则变为 i -> j -> k"""
        pairs = [(i, j) for i in range(n_nodes) for j in range(n_nodes) if g[i][j] == 1 and g[j][i] == 0]  # 所有 i - j 点对

        for i, j in pairs:
            all_k = [k for k in range(n_nodes) if (g[j][k] == 1 and g[k][j] == 1) and (g[i][k] == 0 and g[k][i] == 0)]

            if len(all_k) > 0:
                for k in all_k:
                    g[j][k] = 1
                    g[k][j] = 0

        return g

    def rule2(g):
        """Rule 2: 如果存在链 i -> k -> j ，则把 i - j 变为 i -> j"""
        pairs = [(i, j) for i in range(n_nodes) for j in range(n_nodes) if g[i][j] == 1 and g[j][i] == 1]  # 所有 i - j 点对

        for i, j in pairs:
            all_k = [k for k in range(n_nodes) if (g[i][k] == 1 and g[k][i] == 0) and (g[k][j] == 1 and g[j][k] == 0)]

            if len(all_k) > 0:
                g[i][j] = 1
                g[j][1] = 0

        return g

    def rule3(g):
        """Rule 3: 如果存在 i - k1 -> j 和 i - k2 -> j ，且 k1, k2 不相邻，则把 i - j 变为 i -> j"""
        pairs = [(i, j) for i in range(n_nodes) for j in range(n_nodes) if g[i][j] == 1 and g[j][i] == 1]  # 所有 i - j 点对

        for i, j in pairs:
            all_k = [k for k in range(n_nodes) if (g[i][k] == 1 and g[k][i] == 1) and (g[k][j] == 1 and g[j][k] == 0)]

            if len(all_k) >= 2:
                for k1, k2 in combinations(all_k, 2):
                    if g[k1][k2] == 0 and g[k2][k1] == 0:
                        g[i][j] = 1
                        g[j][i] = 0
                        break

        return g

    # Rule 4: 如果存在链 i - k -> l 和 k -> l -> j，且 k 和 l 不相邻，把 i - j 改为 i -> j
    # 显然，这种情况不可能存在，所以不需要考虑 rule4

    # 相邻点对
    pairs = [(i, j) for i in range(n_nodes) for j in range(n_nodes) if G[i][j] == 1]

    # 把 x - y - z 变为 x -> y <- z
    for x, y in sorted(pairs, key=lambda x:(x[1], x[0])):
        all_z = [z for z in range(n_nodes) if G[y][z] == 1 and z != x]

        for z in all_z:
            if G[x][z] == 0 and y not in O[x][z]:
                G[x][y] = G[z][y] = 1
                G[y][x] = G[y][z] = 0

    # Orientation rule 1 - rule 3
    old_G = np.zeros((n_nodes, n_nodes))

    while not np.array_equal(old_G, G):
        old_G = G.copy()

        G = rule1(G)
        G = rule2(G)
        G = rule3(G)

    return np.array(G)

def pc(suff_stat, alpha: float = 0.05, verbose: bool = False):
    G, O = skeleton(suff_stat, alpha)  # 骨架
    cpdag = extend_cpdag(G, O)  # 扩展为 CPDAG

    if verbose:
        print(cpdag)

    return cpdag
