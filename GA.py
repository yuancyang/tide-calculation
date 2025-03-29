import numpy as np
import matplotlib.pyplot as plt
from GLpower import probabilistic_power_flow

# 定义遗传算法类
class GeneticAlgorithm:
    def __init__(self):
        self.popsize = 10      # 群体大小
        self.chromlength = 9  # 串的长度（个体长度）
        self.pc = 0.6         # 交叉概率
        self.pm = 0.1         # 变异概率
        self.xlim = [-5, 5]   # x范围
        self.G = 10          # 迭代次数

    # 定义目标函数
    def objective_function(self, x):
        # 目标函数
        k = round(probabilistic_power_flow(PESS=x),4)
        #k=x**2
        print(f"k = {k}")
        m = 1-k
        print(f"m = {m}\n")
        
        return m
        #return abs(x**2)
    # 二进制转十进制函数
    def bin_to_dec(self, pop):
        # 获取种群大小
        popsize = len(pop)
        # 创建存储十进制值的数组
        dec = np.zeros(popsize)
        # 遍历每个个体
        for i in range(popsize):
            # 计算二进制串的十进制值
            # 假设 pop[i] 是一个二进制数组
            bin_str = ''.join(map(str, pop[i].astype(int)))  # 将二进制数组转换为字符串
            # 检查第一位以确定符号
            sign = -1 if bin_str[0] == '1' else 1  # 如果第一位是1，sign为-1；否则为1
            # 提取剩余的二进制位并转换为十进制数
            remaining_bin_str = bin_str[1:]  # 获取从第二位开始的剩余位
            dec_value = int(remaining_bin_str, 2)  # 将剩余的二进制字符串转换为十进制数
            # 根据符号调整十进制值
            dec[i] = sign * dec_value  # 如果sign为-1，则dec[i]为负值
        #     bin_str = ''.join(map(str, pop[i].astype(int)))
        #     # 将二进制字符串转换为十进制数
        #     dec[i] = int(bin_str, 2)
        # # 将十进制数映射到实际范围[xlim[0], xlim[1]]
        dec = dec / 100
        print(f"十进制个体集合{dec}")
        # 返回映射后的十进制数组
        return dec

    # 计算适应度函数
    def calculate_fitness(self, dec_pop):
        # 对每个个体计算目标函数值并返回数组
        print(f"适应度{np.array([self.objective_function(x) for x in dec_pop])}")
        return np.array([self.objective_function(x) for x in dec_pop])

    # 复制（轮盘赌选择）
    def selection(self, pop, fitness):
        p = fitness / np.sum(fitness)
        # 计算累积概率
        cumsum_p = np.cumsum(p)
        # 创建新种群数组
        newpop = np.zeros_like(pop)
        # 生成随机数数组并排序
        r = np.sort(np.random.rand(self.popsize))
        i, j = 0, 0
        # 轮盘赌选择过程
        while j < self.popsize:
            if r[j] < cumsum_p[i]:
                newpop[j] = pop[i]
                j += 1
            else:
                i += 1
        return newpop

    # 交叉
    def crossover(self, pop):
        # 复制原始种群
        newpop = pop.copy()
        # 每次选择两个个体进行交叉
        for i in range(0, self.popsize-1, 2):
            # 根据交叉概率决定是否进行交叉
            if np.random.rand() < self.pc:
                # 随机选择两个交叉点
                points = np.sort(np.random.choice(self.chromlength, 2, replace=False))
                # 保存第一个个体的部分基因
                temp = newpop[i, points[0]:points[1]].copy()
                # 将第二个个体的对应部分基因复制给第一个个体
                newpop[i, points[0]:points[1]] = newpop[i+1, points[0]:points[1]]
                # 将保存的基因复制给第二个个体
                newpop[i+1, points[0]:points[1]] = temp
        return newpop

    # 变异操作
    def mutation(self, pop):
        # 复制原始种群
        newpop = pop.copy()
        # 遍历每个个体
        for i in range(self.popsize):
            # 根据变异概率决定是否进行变异
            if np.random.rand() < self.pm:
                # 随机选择一个变异位置
                point = np.random.randint(0, self.chromlength)
                # 将选中位置的基因取反（0变1，1变0）
                newpop[i, point] = 1 - newpop[i, point]
        return newpop

    # 绘制结果函数
    def plot_results(self, dec_pop, fx, generation):
        # 设置子图1
        plt.subplot(1, 2, 1)
        # 只绘制当前种群个体的位置
        plt.plot(dec_pop, fx, 'o')
        plt.title(f'第{generation}次迭代进化')
        plt.grid(True)
        plt.xlim(self.xlim)
        plt.ylim([0, 1])
        
        plt.pause(0.2)

    # 主运行函数
    def run(self):
        # 随机初始化二进制种群
        pop = np.random.randint(0, 2, (self.popsize, self.chromlength))
        # 创建存储最优x值的历史记录数组
        x_history = np.zeros(self.G)
        # 创建存储最优适应度值的历史记录数组
        y_history = np.zeros(self.G)
        
        # 创建图形窗口
        plt.figure(figsize=(12, 5))
        
        # 计算初始种群的十进制值
        dec_pop = self.bin_to_dec(pop)
        # 计算初始种群的适应度
        fx = self.calculate_fitness(dec_pop)
        # 绘制初始种群的分布
        self.plot_results(dec_pop, fx, 1)
        
        # 找到最优个体的索引
        best_idx = np.argmax(fx)
        # 记录初始最优适应度
        y_history[0] = fx[best_idx]
        # 记录初始最优x值
        x_history[0] = dec_pop[best_idx]
        
        # 开始迭代进化
        for i in range(1, self.G):
            # 选择操作
            newpop = self.selection(pop, fx)
            # 交叉操作
            newpop = self.crossover(newpop)
            # 变异操作
            newpop = self.mutation(newpop)
            
            # 计算新种群的十进制值
            new_dec_pop = self.bin_to_dec(newpop)
            # 计算新种群的适应度
            new_fx = self.calculate_fitness(new_dec_pop)
            
            # 找到新种群中优于原种群的个体
            better_indices = new_fx > fx
            # 用优秀个体替换原种群中的对应个体
            pop[better_indices] = newpop[better_indices]
            
            # 更新当前种群的十进制值
            dec_pop = self.bin_to_dec(pop)
            # 更新当前种群的适应度
            fx = self.calculate_fitness(dec_pop)
            
            # 绘制当前代的结果
            self.plot_results(dec_pop, fx, i+1)
            
            # 找到当前代最优个体
            best_idx = np.argmax(fx)
            # 记录当前代最优适应度
            y_history[i] = fx[best_idx]
            # 记录当前代最优x值
            x_history[i] = dec_pop[best_idx]
            print(x_history)
            print(y_history)
            # 绘制进化曲线
            plt.subplot(1, 2, 2)
            plt.plot(range(1, i+2), y_history[:i+1])
            plt.title('适应度进化曲线')
            plt.grid(True)
            
        # 找到所有代数中的最优解
        best_gen = np.argmax(y_history)
        print(f'找到的最优解位置为：{x_history[best_gen]:.4f}')
        print(f'对应最优解为：{y_history[best_gen]:.4f}')
        
        plt.show()

if __name__ == "__main__":
    ga = GeneticAlgorithm()
    ga.run()