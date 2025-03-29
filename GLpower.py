import cmath
import math
from copy import copy
from power import Node,power_flow
import numpy as np
from b import Y
num_unsuitable=0

def check_voltage():
    V = []
    global nodes
    for i in nodes:
        V.append(i.U)
    #print(V)
    if max(V) > 1.05 or min(V) < 0.95:
        global num_unsuitable
        num_unsuitable+=1

def probabilistic_power_flow(PESS):
    N = 200
    global num_unsuitable
    num_unsuitable = 0  # 每次调用时重置计数器
    
    for i in range(N):
        mean = 2
        std_dev = 1  # 均值为2，方差为 1，标准差为 sqrt(1) = 1
        random_number = np.random.normal(mean, std_dev)
        PDG = round(random_number, 2)# 保留两位小数
        #print(PDG)
        global nodes
        nodes = [
            Node(1, 5, 0, 1.05, 0, "平衡"),#平衡节点
            Node(2, -0.16, -0.1, 1.05, 0, "PQ"),#负荷
            Node(3, -0.55, -0.25, 1.05, 0, "PQ"),#负荷
            Node(4, PESS, 0, 1.05, 0, "PV"),#储能
            Node(5, PDG, 0, 1.05, 0, "PV")#光伏
        ]

        power_flow(Y,nodes,0.00001)
        check_voltage()
        
    result = num_unsuitable/N
    return result

if __name__ == "__main__":
    #print(num_unsuitable)
    print(probabilistic_power_flow(PESS=1))
    #power_flow(Y,nodes,0.00001)







    



