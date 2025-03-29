class Node:
    """
        节点类
    """

    def __init__(i,id):
        i.id = id
    def A(i):
        i.id = i.id + 1
        print(i.id)
    def B(i):
        i.A()
A = Node(1)
A.A()
A.B()
