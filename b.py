import cmath


import math
import numpy as np
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
        y = 1/z
        Y = y/k
        Yt = y * (1-k)/(k*k)
        Yf = y * (k-1)/k
        self.Y = Y
        self.Yf = Yf
        self.Yt = Yt

class Line:
    """
        线路类
    """

    def __init__(self,id,from_node_id,to_node_id,R,X,B1_2):
        Z = complex(R,X)
        Y = 1/Z
        self.id = id
        self.from_node_id = from_node_id - 1
        self.to_node_id = to_node_id - 1
        # self.R = complex(R,0)
        # self.X = complex(0,X)
        # self.Z = complex(R,X)
        self.Y0 = complex(0,B1_2)
        self.Y = Y

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

#计算节点的导纳矩阵
def daona_node(numbers,daona1,daona2):
    daona = [[0 for i in range(numbers)] for j in range(numbers)]
    for i in range(numbers):
        for j in range(numbers):
            daona[i][j] = daona1[i][j] + daona2[i][j]
            #print(daona[i][j])
    return daona

def start():
    lines = [
        Line(1, 1, 2, 0.0035, 0.0411, 0.3494),
        Line(2, 2, 3, 0.001, 0.025, 0.375),
        Line(3, 2, 4, 0.0013, 0.0151, 0.1286),
        Line(4, 4, 3, 0.007, 0.0086, 0.073),
        Line(5, 3, 5, 0.0013, 0.0213, 0.1107),
        # Line(6, 3, 18, 0.0011, 0.0133, 0.1069),
        # Line(7, 4, 5, 0.0008, 0.0128, 0.0671),
        # Line(8, 4, 14, 0.0008, 0.0129, 0.0691),
        # Line(9, 5, 6, 0.0002, 0.0026, 0.0217),
        # Line(10, 5, 8, 0.0008, 0.0112, 0.0738),
        # Line(11, 6, 7, 0.0006, 0.0092, 0.0565),
        # Line(12, 6, 11, 0.0007, 0.0082, 0.0695),
        # Line(13, 7, 8, 0.0004, 0.0046, 0.039),
        # Line(14, 8, 9, 0.0023, 0.0363, 0.1902),
        # Line(15, 9, 31, 0.001, 0.025, 0.6),
        # Line(16, 10, 11, 0.0004, 0.0043, 0.0365),
        # Line(17, 10, 13, 0.0004, 0.0043, 0.0365),
        # Line(18, 13, 14, 0.0009, 0.0101, 0.0862),
        # Line(19, 14, 15, 0.0018, 0.0217, 0.183),
        # Line(20, 15, 16, 0.0009, 0.0094, 0.0855),
        # Line(21, 16, 17, 0.0007, 0.0089, 0.0671),
        # Line(22, 16, 19, 0.0016, 0.0195, 0.152),
        # Line(23, 16, 21, 0.0008, 0.0135, 0.1274),
        # Line(24, 16, 24, 0.0003, 0.0059, 0.034),
        # Line(25, 17, 18, 0.0007, 0.0082, 0.066),
        # Line(26, 17, 27, 0.0013, 0.0173, 0.1608),
        # Line(27, 21, 22, 0.0008, 0.014, 0.1283),
        # Line(28, 22, 23, 0.0006, 0.0096, 0.0923),
        # Line(29, 23, 24, 0.0022, 0.035, 0.1805),
        # Line(30, 25, 26, 0.0032, 0.0323, 0.2565),
        # Line(31, 26, 27, 0.0014, 0.0147, 0.1198),
        # Line(32, 26, 28, 0.0043, 0.0474, 0.3901),
        # Line(33, 26, 29, 0.0057, 0.0625, 0.5145),
        # Line(34, 28, 29, 0.0014, 0.0151, 0.1245)
    ]
    # transformers = [
    #     Transformer(1, 11, 12, 0.0016, 0.0435, 1.006),
    #     Transformer(2, 13, 12, 0.0016, 0.0435, 1.006),
    #     Transformer(3, 39, 6, 0, 0.025, 1.070),
    #     Transformer(4, 32, 10, 0, 0.02, 1.070),
    #     Transformer(5, 33, 19, 0.0007, 0.0142, 1.070),
    #     Transformer(6, 34, 20, 0.0009, 0.018, 1.009),
    #     Transformer(7, 35, 22, 0, 0.0143, 1.025),
    #     Transformer(8, 36, 23, 0.0005, 0.0272, 1.000),
    #     Transformer(9, 37, 25, 0.0006, 0.0232, 1.025),
    #     Transformer(10, 30, 2, 0, 0.0181, 1.025),
    #     Transformer(11, 38, 29, 0.0008, 0.0156, 1.025),
    #     Transformer(12, 20, 19, 0.0007, 0.0138, 1.060)
    # ]
    # daona1 = daona_line(39,lines)
    # daona2 = daona_transformer(39,transformers)
    # daona = daona_node(39,daona1,daona2)
    # 把daona转换成矩阵
    # J = np.mat(daona)

    # j = J.A
    # j = np.linalg.inv(j)
    # #print(J)
    # return daona
    # print(daona1)
    # print(daona2)

#Y = start()
Y = [[complex(0,-10),complex(0,9.5238),complex(0,0),complex(0,0),complex(0,0)],
     [complex(0,9.5238),complex(10.3422,-36.5123),complex(-4.1096,10.9589),complex(-6.2326,16.5161),complex(0,0)],
     [complex(0,0),complex(-4.1096,10.9589),complex(10.9589,-34.9130),complex(-6.8493,18.2648),complex(0,5.7143)],
     [complex(0,0),complex(-6.2326,16.5161),complex(-6.8493,18.2648),complex(13.0819,-34.7059),complex(0,0)],
     [complex(0,0),complex(0,0),complex(0,5.7143),complex(0,0),complex(0,-5.7143)]]
#print(Y)







