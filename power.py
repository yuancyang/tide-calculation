# 导入复数
import cmath
import math
from copy import copy

import numpy as np
from b import Y

# from b import calc_Y_matrix

class Node:
    """
        节点类
    """

    def __init__(self, id, P, Q, U, phase, n_type):
        self.id = id
        self.n_type = n_type
        self.P = P
        self.Q = Q
        self.U = U
        self.phase = phase
        self.delta_U = 0
        self.delta_phase = 0
        self.delta_P = 0
        self.delta_Q = 0


class Line:
    """
        线路类
    """

    def __init__(self, from_node, to_node, G, B):
        self.id = str(from_node.id) + "-" + str(to_node.id)
        self.from_node = from_node
        self.to_node = to_node
        self.G = G
        self.B = B


class PowerFlow:
    """
        潮流计算类
    """

    def __init__(self, Y: list[list[complex]], nodes: list[Node]):
        """
            潮流计算类初始化
            param：
                Y -> 节点导纳矩阵
                nodes -> 节点列表
        """
        self.Y = Y
        self.nodes = nodes
        # 所有PQ节点组成的列表
        self.PQ_nodes = []
        # 所有PQ节点组成的列表
        self.PV_nodes = []
        # n是除平衡节点外所有节点的数量
        self.n = len(nodes) - 1
        # m是PQ节点的数量
        self.m = 0

        # 后面会用到 y的列向量中的P分量和Q分量
        self.delta_P = []
        self.delta_Q = []
        self.pure_nodes = []
        # 初始化node列表
        self.init_node_list()



    def init_node_list(self):
        """
            初始化node列表
        """
        for node in self.nodes:
            if node.n_type == "PQ":
                self.PQ_nodes.append(node)
                self.m += 1
                #print(self.m)
            elif node.n_type == "PV":
                self.PV_nodes.append(node)
                

        self.pure_nodes = copy(self.nodes)
        for node in self.pure_nodes:
            if node.n_type == "平衡":
                self.pure_nodes.remove(node)



    def start_calculation(self, accuracy):
        """
            开始计算
            param：
                accuracy -> 计算的精度选择
        """
        # 开始计算
        n = 0
        while True:
            # 计算y的不平衡量
            self.calc_unbalance()
            mo_P = []
            mo_Q = []
            for i in self.delta_P:
                mo_P.append(i.delta_P)
            for i in self.delta_Q:
                mo_Q.append(i.delta_Q)
            # 判断精度
            if max(mo_P) < accuracy and max(mo_Q) < accuracy:
                break
            # print(max(mo_P))
            J = self.calc_Jacobian()
            j = J.A
            # print(J)
            self.calc_delta_U(J)
            self.revise()

            #print("------------")
            #for node in self.nodes:
            #         print(node.U)
                #     print(node.delta_U)
                #     print(node.P)
                #     print(node.Q)
                #print(node.delta_P)
            #     print(node.delta_Q)
            #     print(node.phase)
            #     print(node.n_type)
            #     print(node.id)
            #     print("-----------")
            n += 1
            if n > 200:
                #print("不收敛")
                self.nodes[4].U = 0
                break
            # print(f"\n第{n}次迭代")
            # print(f"{'节点ID':<10}{'相位':<10}{'电压':<10}{'ΔP':<10}{'ΔQ':<10}")
            # for node in self.nodes:
            #     print(f"{node.id:<10}{node.phase:<10.4f}{node.U:<10.4f}{node.delta_P:<10.4f}{node.delta_Q:<10.4f}")
    # 计算不平衡量
    def calc_unbalance(self):
        """
            计算delta_y
        """
        # 计算不平衡量
        for node in self.nodes:
            # 确保不是平衡节点
            if node.n_type == "平衡":
                continue
            a = 0
            b = 0
            c = 0
            i = node.id
            for no in self.nodes:
                j = no.id
                a = a + self.get_node_by_id(j).U * (
                            Y[i-1][j-1].real * math.cos(self.get_node_by_id(i).phase - self.get_node_by_id(j).phase) + Y[i-1][j-1].imag * math.sin(
                        self.get_node_by_id(i).phase - self.get_node_by_id(j).phase))
                b = b + self.get_node_by_id(j).U * (
                            Y[i-1][j-1].real * math.sin(self.get_node_by_id(i).phase - self.get_node_by_id(j).phase) - Y[i-1][j-1].imag * math.cos(
                        self.get_node_by_id(i).phase - self.get_node_by_id(j).phase))

                node.delta_P = node.P - node.U * a
                node.delta_Q = node.Q - node.U * b
        if len(self.delta_P) > 0:
            return
        for node in self.nodes:
            if node.n_type == "PQ":
                self.delta_P.append(node)
                self.delta_Q.append(node)
            elif node.n_type == "PV":
                self.delta_P.append(node)

    def get_node_by_id(self,id):
        for node in self.nodes:
            if node.id == id:
                return node

    def get_attr(self,i, j):
        Gij = Y[i-1][j-1].real
        Bij = Y[i-1][j-1].imag
        Sinij = math.sin(self.get_node_by_id(i).phase - self.get_node_by_id(j).phase)
        Cosij = math.cos(self.get_node_by_id(i).phase - self.get_node_by_id(j).phase)
        Ui = self.get_node_by_id(i).U
        Uj = self.get_node_by_id(j).U
        return Gij, Bij, Sinij, Cosij, Ui, Uj

    # 求雅克比矩阵
    def calc_Jacobian(self):

        n = len(self.pure_nodes)  # 使用实际的纯节点数量
        m = len(self.PQ_nodes)    # PQ节点数量
       #print(n,m)
        H = [[0 for i in range(n)] for j in range(n)]
        N = [[0 for i in range(m)] for j in range(n)]
        M = [[0 for i in range(n)] for j in range(m)]
        L = [[0 for i in range(m)] for j in range(m)]

        # 计算矩阵 H
        for from_node in self.delta_P:
            i = from_node.id
            x = self.pure_nodes.index(from_node)
            for to_node in self.delta_P:
                j = to_node.id
                y = self.pure_nodes.index(to_node)
                Gij, Bij, Sinij, Cosij, Ui, Uj = self.get_attr(i, j)
                if from_node == to_node:
                    H[x][y] = Ui * Ui * Bij + from_node.Q - from_node.delta_Q
                else:
                    H[x][y] = -Ui * Uj * (Gij * Sinij - Bij * Cosij)
        for from_node in self.delta_P:  # 计算矩阵N
            i = from_node.id
            x = self.pure_nodes.index(from_node)
            for to_node in self.delta_Q:
                j = to_node.id
                y = self.pure_nodes.index(to_node)
                Gij, Bij, Sinij, Cosij, Ui, Uj = self.get_attr(i, j)
                if from_node == to_node:
                    N[x][y] = -Ui * Ui * Gij - from_node.P + from_node.delta_P
                else:
                    N[x][y] = -Ui * Uj * (Gij * Cosij + Bij * Sinij)
        for from_node in self.delta_Q:  # 计算矩阵M
            i = from_node.id
            x = self.pure_nodes.index(from_node)
            for to_node in self.delta_P:
                j = to_node.id
                y = self.pure_nodes.index(to_node)
                Gij, Bij, Sinij, Cosij, Ui, Uj = self.get_attr(i, j)
                if from_node == to_node:
                    M[x][y] = Ui * Ui * Gij - from_node.P + from_node.delta_P
                else:
                    M[x][y] = Ui * Uj * (Gij * Cosij + Bij * Sinij)
        for from_node in self.delta_Q:  # 计算矩阵L
            i = from_node.id
            x = self.pure_nodes.index(from_node)
            for to_node in self.delta_Q:
                j = to_node.id
                y = self.pure_nodes.index(to_node)
                Gij, Bij, Sinij, Cosij, Ui, Uj = self.get_attr(i, j)
                if from_node == to_node:
                    L[x][y] = Ui * Ui * Bij - from_node.Q + from_node.delta_Q
                else:
                    L[x][y] = -Ui * Uj * (Gij * Sinij - Bij * Cosij)
        # 分块矩阵合成J

        J = np.bmat([[H, N], [M, L]])

        return J


    # 计算修正方程
    def calc_delta_U(self,J):
        j_inv = np.linalg.inv(J)

        delta_num = []
        for node in self.delta_P:
            delta_num.append(node.delta_P)
        for node in self.delta_Q:
            delta_num.append(node.delta_Q)
        vector = np.array(delta_num)
        # print(vector)
        vector = vector.reshape(self.m + self.n, 1)
        # print(vector)
        delta_x = - j_inv * vector
        for i in range(len(self.delta_P)):
            self.delta_P[i].delta_phase = delta_x[i, 0]
        for j in range(len(self.delta_Q)):
            self.delta_P[j].delta_U = delta_x[len(self.delta_P) + j, 0]


    # 修正
    def revise(self):
        # for i in range(len(self.delta_P)):
        #     self.nodes[i].phase = self.nodes[i].phase + self.nodes[i].delta_phase
        # for j in range(len(self.delta_Q)):
        #     self.nodes[j].U = self.nodes[j].U + self.nodes[j].delta_U * self.nodes[j].U
        #
        for node in self.delta_P:
            node.phase = node.phase + node.delta_phase
        for node in self.delta_Q:
            node.U = node.U + node.delta_U * node.U
def power_flow(Y,nodes,ACCURACY):
    powerFlow = PowerFlow(Y,nodes)
    powerFlow.start_calculation(ACCURACY)


if __name__ == "__main__":
    # Y = [
    #     [complex(2.904949, -11.503060), complex(0, 5.318187),complex(-1.244978, 2.371387),complex(-1.659971, 3.161849)],
    #     [complex(0, 5.318187), complex(0, -4.663848), complex(0, 0), complex(0, 0)],
    #     [complex(-1.244978, 2.371387),complex(0, 0),complex(2.074963, -3.909204) , complex(-0.829985, 1.580924)],
    #     [complex(-1.659971, 3.161849), complex(0, 0), complex(-0.829985, 1.580924), complex(2.489956, -4.703977)]
    # ]
    Y = [[complex(0,-10),complex(0,9.5238),complex(0,0),complex(0,0),complex(0,0)],
        [complex(0,9.5238),complex(10.3422,-36.5123),complex(-4.1096,10.9589),complex(-6.2326,16.5161),complex(0,0)],
        [complex(0,0),complex(-4.1096,10.9589),complex(10.9589,-34.9130),complex(-6.8493,18.2648),complex(0,5.7143)],
        [complex(0,0),complex(-6.2326,16.5161),complex(-6.8493,18.2648),complex(13.0819,-34.7059),complex(0,0)],
        [complex(0,0),complex(0,0),complex(0,5.7143),complex(0,0),complex(0,-5.7143)]]
    nodes = [
        Node(1, 0, 0, 1.05, 0, "平衡"),
        Node(2, -0.16, -0.1, 1.05, 0, "PQ"),
        Node(3, 0, 0, 1.05, 0, "PQ"),
        Node(4, -0.55, -0.25, 1.05, 0, "PQ"),
        Node(5, 0.2, 0, 1.05, 0, "PV")

    ]

    # nodes =[
    # Node(1, 0, 0, 1, 0, "PQ"),
    # Node(2, 0, 0, 1, 0,"PQ"),
    # Node(3, -322, -2.4, 1, 0, "PQ"),
    # Node(4, -500, -184, 1, 0, "PQ"),
    # Node(5, 0, 0, 1, 0, "PQ"),
    # Node(6, 0, 0, 1, 0, "PQ"),
    # Node(7, -233.8, -84, 1, 0, "PQ"),
    # Node(8, -522, -176, 1, 0, "PQ"),
    # Node(9, 0, 0, 1, 0, "PQ"),
    # Node(10, 0, 0, 1, 0, "PQ"),
    # Node(11, 0, 0, 1, 0, "PQ"),
    # Node(12, -8.5, -88, 1, 0, "PQ"),
    # Node(13, 0, 0, 1, 0, "PQ"),
    # Node(14, 0, 0, 1, 0, "PQ"),
    # Node(15, -320, -153, 1, 0, "PQ"),
    # Node(16, -329.4, -32.3, 1, 0, "PQ"),
    # Node(17, 0, 0, 1, 0, "PQ"),
    # Node(18, -158, -30, 1, 0, "PQ"),
    # Node(19, 0, 0, 1, 0, "PQ"),
    # Node(20, -680, -103, 1, 0, "PQ"),
    # Node(21, -274, -115, 1, 0, "PQ"),
    # Node(22, 0, 0, 1, 0, "PQ"),
    # Node(23, -247.5, -84.6, 1, 0, "PQ"),
    # Node(24, -308.6, -92.2, 1, 0, "PQ"),
    # Node(25, -224, -47.2, 1, 0, "PQ"),
    # Node(26, -139, -17, 1, 0, "PQ"),
    # Node(27, -281, -75, 1, 0, "PQ"),
    # Node(28, -206, -27.6, 1, 0, "PQ"),
    # Node(29, -283.3, -26.9, 1, 0, "PQ"),
    # Node(30, 250, 143.7, 1, 0, "PQ"),
    # Node(31, 563.9, 202.5, 1, 0, "平衡"),
    # Node(32, 650, 205.4, 1, 0, "PQ"),
    # Node(33, 632, 108.6, 1, 0, "PQ"),
    # Node(34, 508, 166.8, 1, 0, "PQ"),
    # Node(35, 650, 211.6, 1, 0, "PQ"),
    # Node(36, 560, 100.4, 1, 0, "PQ"),
    # Node(37, 540, 2.4, 1, 0, "PQ"),
    # Node(38, 830, 22.5, 1, 0, "PQ"),
    # Node(39, -104, -160.7, 1, 0, "PQ")
    # ]

    #ACCURACY = 0.00001

    power_flow(Y,nodes,0.00001)




