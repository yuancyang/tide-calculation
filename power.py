# 导入复数
import cmath
import math
from copy import copy
import cmath
import numpy as np
import pandas as pd
class Line:
    '''
    线路类
    '''
    def __init__(self,id,from_node_id,to_node_id,R,X,B1_2):
        Z = complex(R,X)
        Y = 1/Z
        self.id = id
        self.from_node_id = from_node_id - 1
        self.to_node_id = to_node_id - 1
        self.R = complex(R,0)
        self.X = complex(0,X)
        self.Z = complex(R,X)
        self.Y0 = complex(0,B1_2)
        self.Y = Y
        self.pass_I_2 = 0
        self.pass_I_16 = 0

    def Z_id(self):
        return self.Z
    def from_node_id_line(self):
        return self.from_node_id

    def to_node_id_line(self):
        return self.to_node_id

class Transformer:
    '''
    变压器类
    '''
    def __init__(self,id, from_node,to_node, R,X,k):
        self.id = id
        self.from_node = from_node - 1
        self.to_node = to_node - 1
        self.R = R
        self.X = X
        self.k = k
        z = complex(R, X)
        self.z = z
        y = 1/z
        Y = y/k
        Yt = y * (1-k)/(k*k)
        Yf = y * (k-1)/k
        self.Y = Y
        self.Yf = Yf
        self.Yt = Yt
        self.pass_I_2 = 0
        self.pass_I_16= 0

    # def z_id(self):
    #     return self.z
    def from_node_id(self):
        return self.from_node

    def to_node_id(self):
        return self.to_node

    # def i_byq(self):
    #     a =get_node_by_id(self.from_node).U/get_node_by_id(self.to_node).U
    #     i = get_node_by_id(self.from_node).U - get_node_by_id(self.to_node).U
    #     return i

class fadianji:
    def __init__(self,id, node, xd):
        self.id = id
        self.node = node
        self.xd = xd
        j =complex(0,xd)
        yi = 1/j
        self.yi = yi

class fuhe:
    def __init__(self,id, node, Pldi,Qldi,U):
        self.id = id
        self.node = node
        self.Pldi = Pldi
        self.Qldi = Qldi
        self.U = U
        j = complex(0,1)
        yi = (Pldi - j*Qldi)/(U*U)
        self.yi = yi


def polar_to_rectangular(amplitude, phase):
    # 将极坐标形式转换为复数形式
    complex_number = amplitude * cmath.exp(1j * phase)  # 使用cmath库中的exp函数表示复指数形式

    return complex_number


def angle_of_complex_number(complex_num):
    angle = cmath.phase(complex_num)  # 使用 cmath.phase() 函数计算复数的相角（幅角）
    return angle

#计算线路的导纳矩阵
def daona_line(numbers1,lines):
    daona1 = [[0 for i in range(numbers1)] for j in range(numbers1)]
    for line in lines:
        i = line.from_node_id
        j = line.to_node_id
        daona1[i][i] += line.Y + line.Y0
        daona1[j][j] += line.Y + line.Y0
        daona1[i][j] += -line.Y
        daona1[j][i] += -line.Y
    return daona1

#计算变压器的导纳矩阵
def daona_transformer(numbers,transformers):
    daona2 = [[0 for i in range(numbers)] for j in range(numbers)]
    for transformer in transformers:
        i = transformer.from_node
        j = transformer.to_node
        daona2[i][i] += transformer.Y + transformer.Yf
        daona2[j][j] += transformer.Y + transformer.Yt
        daona2[i][j] += -transformer.Y
        daona2[j][i] += -transformer.Y
        # print(daona2[i][i])
    return daona2


def daona_fadianji(numbers2,fadianjis):
    daona3 = [[0 for i in range(numbers2)] for j in range(numbers2)]
    for fadianji in fadianjis:
        i = fadianji.node-1
        daona3[i][i] += fadianji.yi
    return daona3

def daona_fuhe(numbers3,fuhes):
    daona4 = [[0 for i in range(numbers3)] for j in range(numbers3)]
    for fuhe in fuhes:
        i = fuhe.node-1
        daona4[i][i] += fuhe.yi
    return daona4



#计算节点的导纳矩阵
def daona_node_fadianji_fuhe(numbers,daona1,daona2,daona3,daona4):
    daona = [[0 for i in range(numbers)] for j in range(numbers)]
    for i in range(numbers):
        for j in range(numbers):
            daona[i][j] = daona1[i][j] + daona2[i][j] + daona3[i][j] + daona4[i][j]
            # print(daona[i][j])
    return daona


def daona_node(numbers,daona1,daona2):
    daona = [[0 for i in range(numbers)] for j in range(numbers)]
    for i in range(numbers):
        for j in range(numbers):
            daona[i][j] = daona1[i][j] + daona2[i][j]
            # print(daona[i][j])
    return daona
def start1():
    daona1 = daona_line(39,lines)
    daona2 = daona_transformer(39,transformers)
    daona3 = daona_fadianji(39,fadianjis)
    daona4 = daona_fuhe(39,fuhes)
    daona5 = daona_node(39,daona1,daona2)
    daona6 = daona_node_fadianji_fuhe(39,daona1,daona2,daona3,daona4)
    J = np.mat(daona5)

    inverse_matrix = np.linalg.inv(daona6)
    # nijuzhen

    print(inverse_matrix)
    global second_column
    global sixteen_column
    second_column = [row[1] for row in inverse_matrix]

    sixteen_column = [row[15] for row in inverse_matrix]

    element = daona6[1][1]
    element16 = daona6[15][15]
    # dierhangdierlie
    global i_duanlu_2_1,i_duanlu_2_2,i_duanlu_16_1,i_duanlu_16_2
    i_duanlu_2_1 = polar_to_rectangular(1.0410625158471356,-0.12206768127774982)/element
    i_duanlu_2_1_modulus = abs(i_duanlu_2_1)
    i_duanlu_2_1_jiaodu = angle_of_complex_number(i_duanlu_2_1)
    print('11')
    print(i_duanlu_2_1,i_duanlu_2_1_modulus,i_duanlu_2_1_jiaodu)
    # duanludianliu
    i_duanlu_2_2 = 1 / element
    i_duanlu_2_2_modulus = abs(i_duanlu_2_2)
    i_duanlu_2_2_jiaodu = angle_of_complex_number(i_duanlu_2_2)
    print('12')
    print(i_duanlu_2_2,i_duanlu_2_2_modulus,i_duanlu_2_2_jiaodu)

    i_duanlu_16_1 = polar_to_rectangular(1.016590799899003,-0.13506286701523007)/element16
    i_duanlu_16_1_modulus = abs(i_duanlu_16_1)
    i_duanlu_16_1_jiaodu = angle_of_complex_number(i_duanlu_16_1)
    print('21')
    print(i_duanlu_16_1,i_duanlu_16_1_modulus,i_duanlu_16_1_jiaodu)
    i_duanlu_16_2 = 1 / element16
    i_duanlu_16_2_modulus = abs(i_duanlu_16_2)
    i_duanlu_16_2_jiaodu = angle_of_complex_number(i_duanlu_16_2)
    print('22')
    print(i_duanlu_16_2,i_duanlu_16_2_modulus,i_duanlu_16_2_jiaodu)
    print('mowucha')
    print((i_duanlu_2_1_modulus-i_duanlu_2_2_modulus)/i_duanlu_2_1_modulus)
    print((i_duanlu_16_1_modulus-i_duanlu_16_2_modulus)/i_duanlu_16_1_modulus)
    print('jiaoduwucha')
    print((i_duanlu_2_1_jiaodu-i_duanlu_2_2_jiaodu)/i_duanlu_2_1_jiaodu)
    print((i_duanlu_16_1_jiaodu-i_duanlu_16_2_jiaodu)/i_duanlu_16_1_jiaodu)
    print('111111111')
    return daona5,element16,element

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
        self.phasor = 0
        self.precision_culc_U_2 = 0
        self.precision_culc_phasor_U_2 = 0
        self.approximate_culc_U_2 = 0
        self.error_2 = 0
        self.precision_culc_U_16 = 0
        self.precision_culc_phasor_U_16 = 0
        self.approximate_culc_U_16 = 0
        self.error_16 = 0

    def culc_phasor(self):
        a = self.U * cmath.cos(self.phase)
        b = self.U * cmath.sin(self.phase)
        self.phasor = complex(a,b)

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
        while True:
            # 计算y的不平衡量
            self.calc_unbalance()
            # todo 看不懂,改！
            mo_P = []
            mo_Q = []
            for i in self.delta_P:
                mo_P.append(i.delta_P)
            for i in self.delta_Q:
                mo_Q.append(i.delta_Q)
            # 判断精度
            if max(mo_P) < accuracy and max(mo_Q) < accuracy:
                break
            print(max(mo_P))
            J = self.calc_Jacobian()
            j = J.A
            # print(J)
            self.calc_delta_U(J)
            self.revise()

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
        H = [[0 for i in range(self.n)] for j in range(self.n)]
        N = [[0 for i in range(self.m)] for j in range(self.n)]
        M = [[0 for i in range(self.n)] for j in range(self.m)]
        L = [[0 for i in range(self.m)] for j in range(self.m)]
        for from_node in self.delta_P:  # 计算矩阵H
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

        for node in self.delta_P:
            node.phase = node.phase + node.delta_phase
        for node in self.delta_Q:
            node.U = node.U + node.delta_U * node.U

if __name__ == "__main__":

    lines = [
        Line(1, 1, 2, 0.0035, 0.0411, 0.3494),
        Line(2, 1, 39, 0.001, 0.025, 0.375),
        Line(3, 2, 3, 0.0013, 0.0151, 0.1286),
        Line(4, 2, 25, 0.007, 0.0086, 0.073),
        Line(5, 3, 4, 0.0013, 0.0213, 0.1107),
        Line(6, 3, 18, 0.0011, 0.0133, 0.1069),
        Line(7, 4, 5, 0.0008, 0.0128, 0.0671),
        Line(8, 4, 14, 0.0008, 0.0129, 0.0691),
        Line(9, 5, 6, 0.0002, 0.0026, 0.0217),
        Line(10, 5, 8, 0.0008, 0.0112, 0.0738),
        Line(11, 6, 7, 0.0006, 0.0092, 0.0565),
        Line(12, 6, 11, 0.0007, 0.0082, 0.0695),
        Line(13, 7, 8, 0.0004, 0.0046, 0.039),
        Line(14, 8, 9, 0.0023, 0.0363, 0.1902),
        Line(15, 9, 39, 0.001, 0.025, 0.6),
        Line(16, 10, 11, 0.0004, 0.0043, 0.0365),
        Line(17, 10, 13, 0.0004, 0.0043, 0.0365),
        Line(18, 13, 14, 0.0009, 0.0101, 0.0862),
        Line(19, 14, 15, 0.0018, 0.0217, 0.183),
        Line(20, 15, 16, 0.0009, 0.0094, 0.0855),
        Line(21, 16, 17, 0.0007, 0.0089, 0.0671),
        Line(22, 16, 19, 0.0016, 0.0195, 0.152),
        Line(23, 16, 21, 0.0008, 0.0135, 0.1274),
        Line(24, 16, 24, 0.0003, 0.0059, 0.034),
        Line(25, 17, 18, 0.0007, 0.0082, 0.066),
        Line(26, 17, 27, 0.0013, 0.0173, 0.1608),
        Line(27, 21, 22, 0.0008, 0.014, 0.1283),
        Line(28, 22, 23, 0.0006, 0.0096, 0.0923),
        Line(29, 23, 24, 0.0022, 0.035, 0.1805),
        Line(30, 25, 26, 0.0032, 0.0323, 0.2565),
        Line(31, 26, 27, 0.0014, 0.0147, 0.1198),
        Line(32, 26, 28, 0.0043, 0.0474, 0.3901),
        Line(33, 26, 29, 0.0057, 0.0625, 0.5145),
        Line(34, 28, 29, 0.0014, 0.0151, 0.1245)
    ]
    transformers = [
        Transformer(1, 11, 12, 0.0016, 0.0435, 1.006),
        Transformer(2, 13, 12, 0.0016, 0.0435, 1.006),
        Transformer(3, 31, 6, 0, 0.025, 1.070),
        Transformer(4, 32, 10, 0, 0.02, 1.070),
        Transformer(5, 33, 19, 0.0007, 0.0142, 1.070),
        Transformer(6, 34, 20, 0.0009, 0.018, 1.009),
        Transformer(7, 35, 22, 0, 0.0143, 1.025),
        Transformer(8, 36, 23, 0.0005, 0.0272, 1.000),
        Transformer(9, 37, 25, 0.0006, 0.0232, 1.025),
        Transformer(10, 30, 2, 0, 0.0181, 1.025),
        Transformer(11, 38, 29, 0.0008, 0.0156, 1.025),
        Transformer(12, 20, 19, 0.0007, 0.0138, 1.060)
    ]
    fadianjis = [
        fadianji(1, 30, 0.0060),
        fadianji(2, 31, 0.0647),
        fadianji(3, 32, 0.0531),
        fadianji(4, 33, 0.0436),
        fadianji(5, 34, 0.1320),
        fadianji(6, 35, 0.0500),
        fadianji(7, 36, 0.0490),
        fadianji(8, 37, 0.0570),
        fadianji(9, 38, 0.0570),
        fadianji(10, 39, 0.0310)
    ]
    fuhes = [
        fuhe(1, 1, 0, 0, 1.0551855338427487),
        fuhe(2, 2, 0, 0, 1.0410625158471356),
        fuhe(3, 3, -3.22, -0.024, 1.0218271006632291),
        fuhe(4, 4, -5.0, -1.84, 0.9984512533536497),
        fuhe(5, 5, 0, 0, 1.0028235163736314),
        fuhe(6, 6, 0, 0, 1.005385054105983),
        fuhe(7, 7, -2.338, -0.84, 0.9958038020953689),
        fuhe(8, 8, -5.22, -1.76, 0.9953828368870602),
        fuhe(9, 9, 0, 0, 1.0381400432581023),
        fuhe(10, 10, 0, 0, 1.0138566306933867),
        fuhe(11, 11, 0, 0, 1.0097219660055474),
        fuhe(12, 12, -0.085, -0.88, 0.9964594481593316),
        fuhe(13, 13, 0, 0, 1.0100967079621959),
        fuhe(14, 14, 0, 0, 1.0053277873195388),
        fuhe(15, 15, -3.2, -1.53, 1.0027326711399587),
        fuhe(16, 16, -3.294, -0.323, 1.016590799899003),
        fuhe(17, 17, 0, 0, 1.021320243739838),
        fuhe(18, 18, -1.58, -0.3, 1.0201664920733124),
        fuhe(19, 19, 0, 0, 1.0440874052625837),
        fuhe(20, 20, -6.8, -1.03, 0.9875928506808153),
        fuhe(21, 21, -2.74, -1.15, 1.0203348923408901),
        fuhe(22, 22, 0, 0, 1.0426269885806785),
        fuhe(23, 23, -2.475, -0.846, 1.0362355337590596),
        fuhe(24, 24, -3.086, -0.922, 1.013928737022894),
        fuhe(25, 25, -2.24, -0.472, 1.0515690439476792),
        fuhe(26, 26, -1.39, -0.17, 1.0448174365254768),
        fuhe(27, 27, -2.81, -0.75, 1.028119493146151),
        fuhe(28, 28, -2.06, -0.276, 1.046092983568797),
        fuhe(29, 29, -2.833, -0.269, 1.0470131349982312)
    ]

    nodes:list[Node] = [
        Node(1, 0, 0, 1, 0, "PQ"),
        Node(2, 0, 0, 1, 0, "PQ"),
        Node(3, -3.22, -0.024, 1, 0, "PQ"),
        Node(4, -5.00, -1.84, 1, 0, "PQ"),
        Node(5, 0, 0, 1, 0, "PQ"),
        Node(6, 0, 0, 1, 0, "PQ"),
        Node(7, -2.338, -0.84, 1, 0, "PQ"),
        Node(8, -5.22, -1.76, 1, 0, "PQ"),
        Node(9, 0, 0, 1, 0, "PQ"),
        Node(10, 0, 0, 1, 0, "PQ"),
        Node(11, 0, 0, 1, 0, "PQ"),
        Node(12, -0.085, -0.88, 1, 0, "PQ"),
        Node(13, 0, 0, 1, 0, "PQ"),
        Node(14, 0, 0, 1, 0, "PQ"),
        Node(15, -3.20, -1.53, 1, 0, "PQ"),
        Node(16, -3.294, -0.323, 1, 0, "PQ"),
        Node(17, 0, 0, 1, 0, "PQ"),
        Node(18, -1.58, -0.30, 1, 0, "PQ"),
        Node(19, 0, 0, 1, 0, "PQ"),
        Node(20, -6.80, -1.03, 1, 0, "PQ"),
        Node(21, -2.74, -1.15, 1, 0, "PQ"),
        Node(22, 0, 0, 1, 0, "PQ"),
        Node(23, -2.475, -0.846, 1, 0, "PQ"),
        Node(24, -3.086, -0.922, 1, 0, "PQ"),
        Node(25, -2.24, -0.472, 1, 0, "PQ"),
        Node(26, -1.39, -0.17, 1, 0, "PQ"),
        Node(27, -2.81, -0.75, 1, 0, "PQ"),
        Node(28, -2.06, -0.276, 1, 0, "PQ"),
        Node(29, -2.833, -0.269, 1, 0, "PQ"),
        Node(30, 2.50, 1.437, 1.03, 0, "PV"),
        Node(31, 5.639, 2.025, 0.982, 0, "平衡"),
        Node(32, 6.50, 2.054, 0.983, 0, "PV"),
        Node(33, 6.32, 1.086, 0.997, 0, "PV"),
        Node(34, 5.08, 1.668, 1.012, 0, "PV"),
        Node(35, 6.50, 2.116, 1.049, 0, "PV"),
        Node(36, 5.60, 1.004, 1.063, 0, "PV"),
        Node(37, 5.40, 0.024, 1.028, 0, "PV"),
        Node(38, 8.30, 0.225, 1.026, 0, "PV"),
        Node(39, -1.04, -1.607, 1.047, 0, "PV")
    ]
    ACCURACY = 0.00001
    Y,element16,element = start1()
    powerFlow = PowerFlow(Y,nodes)
    powerFlow.start_calculation(ACCURACY)
    for node in nodes:
        node.culc_phasor()

    for node in nodes:
        # 2短路
        node.precision_culc_phasor_U_2 = node.phasor - i_duanlu_2_1*second_column[node.id-1]
        node.precision_culc_U_2 = np.abs(node.phasor - i_duanlu_2_1*second_column[node.id-1])
        node.approximate_culc_U_2=np.abs( 1- second_column[node.id-1]/element)
        node.error_2 = abs((node.precision_culc_U_2 - node.approximate_culc_U_2)/node.precision_culc_U_2)

        #16duanlu
        node.precision_culc_phasor_U_16 = node.phasor - i_duanlu_16_1 * sixteen_column[node.id - 1]
        node.precision_culc_U_16 = np.abs(node.phasor - i_duanlu_16_1*sixteen_column[node.id-1])
        node.approximate_culc_U_16=np.abs( 1- sixteen_column[node.id-1]/element16)
        node.error_16 = abs((node.precision_culc_U_16 - node.approximate_culc_U_16)/node.precision_culc_U_16)

    for transformer in transformers:
        #2duanlu
        transformer.pass_I_2 = (powerFlow.get_node_by_id(transformer.from_node).precision_culc_phasor_U_2 - powerFlow.get_node_by_id(transformer.to_node).precision_culc_phasor_U_2/ transformer.k)/transformer.z
        #16duanlu
        transformer.pass_I_16 = (powerFlow.get_node_by_id(transformer.from_node).precision_culc_phasor_U_16 - powerFlow.get_node_by_id(
            transformer.to_node).precision_culc_phasor_U_16 / transformer.k) / transformer.z

    for line in lines:
        #2duanlu
        line.pass_I_2 = (powerFlow.get_node_by_id(line.from_node_id+1).precision_culc_phasor_U_2 - powerFlow.get_node_by_id(line.to_node_id+1).precision_culc_phasor_U_2)/line.Z
        #16duanlu
        line.pass_I_16 = (powerFlow.get_node_by_id(line.from_node_id+1).precision_culc_phasor_U_16 - powerFlow.get_node_by_id(line.to_node_id+1).precision_culc_phasor_U_16)/line.Z