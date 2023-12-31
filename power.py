# 导入复数
import cmath
import math
import numpy as np

class Node:
    """
        节点类
    """
    def __init__(self,id,P,Q,U,phase,n_type):
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
    def __init__(self,from_node,to_node,G,B):
        self.id = str(from_node.id)+"-"+str(to_node.id)
        self.from_node = from_node
        self.to_node = to_node
        self.G = G
        self.B = B







# 计算不平衡量
def calc_unbalance():
    # 计算不平衡量
    a = 0
    b = 0
    for node in nodes:
        for j in range(node.id-1):
            for i in range(len(nodes)-1):
                if nodes[i].id == j :
                    a = a + nodes[i].U * (Y[node.id-1][j].real * math.cos(nodes[node.id-1].phase - nodes[j].phase))
                    b = b + nodes[i].U * (Y[node.id-1][j].imag * math.sin(nodes[node.id-1].phase - nodes[j].phase))


        if node.n_type == "PQ":
            node.delta_P = node.P - node.U * (a+b)
            node.delta_Q = node.Q - node.U * (a-b)
        elif node.n_type == "PV":
            node.delta_P = node.P - node.U * (a+b)

    # 不平衡量列向量
    delta_P = []
    delta_Q = []
    for node in nodes:
        if node.n_type == "PQ":
            delta_P.append(node)
            delta_Q.append(node)
        elif node.n_type == "PV":
            delta_P.append(node)
    return delta_P, delta_Q
def get_node_by_id(id):
    for node in nodes:
        if node.id == id:
            return node
# 求雅克比矩阵
def calc_Jacobian(delta_P, delta_Q,nodes,PV_nodes,PQ_nodes):
    H = [[0 for i in range(n)] for j in range(n)]
    N = [[0 for i in range(m)] for j in range(n)]
    M = [[0 for i in range(n)] for j in range(m)]
    L = [[0 for i in range(m)] for j in range(m)]
    for from_node in delta_P:
        for to_node in delta_P:
            if from_node == to_node:
                # print(nodes[delta_P.index(from_node)].U )
                H[delta_P.index(from_node)][delta_P.index(to_node)] = nodes[delta_P.index(from_node)].U * nodes[delta_P.index(from_node)].U * Y[delta_P.index(from_node)][delta_P.index(to_node)].imag + nodes[delta_P.index(from_node)].delta_Q - nodes[delta_P.index(from_node)].delta_Q
            else:
                H[delta_P.index(from_node)][delta_P.index(to_node)] = 0-nodes[delta_P.index(from_node)].U * nodes[delta_P.index(to_node)].U * (Y[delta_P.index(from_node)][delta_P.index(to_node)].real * math.sin(nodes[delta_P.index(from_node)].phase - delta_P[nodes.index(to_node)].phase) - Y[delta_P.index(from_node)][delta_P.index(to_node)].imag * math.cos(nodes[delta_P.index(from_node)].phase - nodes[delta_P.index(to_node)].phase))
    for from_node in delta_P:
        for to_node in delta_Q:
            if from_node == to_node:
                N[delta_P.index(from_node)][delta_P.index(to_node)] = 0 -nodes[delta_P.index(from_node)].U * nodes[delta_P.index(from_node)].U * Y[delta_P.index(from_node)][delta_P.index(to_node)].real - nodes[delta_P.index(from_node)].delta_P + nodes[delta_P.index(from_node)].delta_P
            else:
                N[delta_P.index(from_node)][delta_P.index(to_node)] = nodes[delta_P.index(from_node)].U * nodes[delta_P.index(to_node)].U * (Y[delta_P.index(from_node)][delta_P.index(to_node)].real * math.cos(nodes[delta_P.index(from_node)].phase - nodes[delta_P.index(to_node)].phase) + Y[delta_P.index(from_node)][delta_P.index(to_node)].imag * math.sin(nodes[delta_P.index(from_node)].phase - nodes[delta_P.index(to_node)].phase))

    for from_node in delta_Q:
        for to_node in delta_P:
            if from_node == to_node:
                M[delta_Q.index(from_node)][delta_P.index(to_node)] = nodes[delta_Q.index(from_node)].U * nodes[delta_Q.index(from_node)].U * Y[delta_Q.index(from_node)][delta_P.index(to_node)].real - nodes[delta_Q.index(from_node)].delta_P + nodes[delta_Q.index(from_node)].delta_P
            else:
                M[delta_Q.index(from_node)][delta_P.index(to_node)] = nodes[delta_Q.index(from_node)].U * nodes[delta_P.index(to_node)].U * (Y[delta_Q.index(from_node)][delta_P.index(to_node)].real * math.cos(nodes[delta_Q.index(from_node)].phase - nodes[delta_P.index(to_node)].phase) + Y[delta_Q.index(from_node)][delta_P.index(to_node)].imag * math.sin(nodes[delta_Q.index(from_node)].phase - nodes[delta_P.index(to_node)].phase))
    for from_node in delta_Q:
        for to_node in delta_Q:
            if from_node == to_node:
                L[delta_Q.index(from_node)][delta_Q.index(to_node)] = nodes[delta_Q.index(from_node)].U * nodes[delta_Q.index(from_node)].U * Y[delta_Q.index(from_node)][delta_Q.index(to_node)].imag - nodes[delta_Q.index(from_node)].delta_Q + nodes[delta_P.index(from_node)].delta_Q

            else:
                L[delta_Q.index(from_node)][delta_Q.index(to_node)] = 0-nodes[delta_Q.index(from_node)].U * nodes[delta_Q.index(to_node)].U * (Y[delta_Q.index(from_node)][delta_Q.index(to_node)].real * math.sin(nodes[delta_Q.index(from_node)].phase - nodes[delta_Q.index(to_node)].phase) - Y[delta_Q.index(from_node)][delta_Q.index(to_node)].imag * math.cos(nodes[delta_Q.index(from_node)].phase - nodes[delta_Q.index(to_node)].phase))

    # 分块矩阵合成J

    J = np.bmat([[H,N],[M,L]])

    return J

# 计算修正方程
def calc_delta_U(J,delta_P,delta_Q,nodes,PV_nodes,PQ_nodes):
    # print(J)
    j_inv = np.linalg.inv(J)

    delta_num = []
    for node in delta_P:
        delta_num.append(node.delta_P)
    for node in delta_Q:
        delta_num.append(node.delta_Q)
    vector = np.array(delta_num)
    # print(vector)
    vector = vector.reshape(m+n,1)
    # print(vector)
    delta_x = - j_inv * vector
    for i in range(len(delta_P)):
        nodes[i].delta_phase = delta_x[i,0]
    for j in range(len(delta_Q)):
        nodes[j].delta_U = delta_x[len(delta_P)+j,0]


# 修正
def revise(nodes,delta_U,PV_nodes,PQ_nodes):
    for i in range(len(delta_P)):
        nodes[i].phase = nodes[i].phase + nodes[i].delta_phase
    for j in range(len(delta_Q)):
        nodes[j].U = nodes[j].U + nodes[j].delta_U * nodes[j].U




if __name__ == "__main__":
    Y = [
        [complex(2.904949,-11.503060),complex(0,5.318187),complex(-1.659971,3.161849),complex(-1.659971,3.161849)],
        [complex(2.904949,-11.503060),complex(0,5.318187),complex(-1.659971,3.161849),complex(-1.659971,3.161849)],
        [complex(2.904949,-11.503060),complex(0,5.318187),complex(-1.659971,3.161849),complex(-1.659971,3.161849)],
        [complex(2.904949,-11.503060),complex(0,5.318187),complex(-1.659971,3.161849),complex(-1.659971,3.161849)]
        ]

    n = 4
    m = 2
    nodes = [
        Node(1,5,1,1,4, "PQ"),
        Node(2,2,3,2,3,"PQ"),
        Node(3,7,1,3,1,"PV"),
        Node(4,8,1,4,2,"PV")
        ]
    PQ_nodes = []
    PV_nodes = []
    for node in nodes:
        if node.n_type == "PQ":
            PQ_nodes.append(node)
        elif node.n_type == "PV":
            PV_nodes.append(node)
    
    ACCURACY = 0.0001

    while True:
        delta_P,delta_Q = calc_unbalance()
        mo_P = []
        mo_Q = []
        for i in delta_P:
            mo_P.append(i.delta_P)
        for i in delta_Q:
            mo_Q.append(i.delta_Q)

        if max(mo_P) < ACCURACY and max(mo_Q) < ACCURACY:
            break
        J = calc_Jacobian(delta_P,delta_Q,nodes,PV_nodes,PQ_nodes)
        calc_delta_U(J,delta_P,delta_Q,nodes,PV_nodes,PQ_nodes)
        revise(nodes,None,PV_nodes,PQ_nodes)

        for node in nodes:
            print(node.U)
            print(node.delta_U)
            print(node.P)
            print(node.Q)
            print(node.delta_P)
            print(node.delta_Q)
            print(node.phase)
            print(node.n_type)
            print(node.id)
            print("-----------")