# 导入随机数生成模块
import random
# 导入数值计算模块
import numpy as np
# 导入绘图模块
import matplotlib.pyplot as plt
# 从GLpower模块导入概率潮流计算函数
from GLpower import probabilistic_power_flow

# 定义遗传算法类
class GeneticAlgorithm:
    # 初始化类的属性
    def __init__(self):
        # 增加种群大小
        self.pop_size = 100  # 从50增加到100
        # 增加迭代次数
        self.n_generations = 50  # 从10增加到50
        # 增加变异概率
        self.mutation_rate = 0.3  # 从0.2增加到0.3
        # 减小交叉概率
        self.crossover_rate = 0.8  # 从0.7增加到0.8
    
    # 初始化种群的方法
    def initialize_population(self):
        # 确保初始种群更加多样化，且保留两位小数
        population = []
        for _ in range(self.pop_size):
            # 在[-5,5]范围内均匀生成两位小数的随机数
            value = round(random.uniform(-5, 5), 2)
            population.append(value)
        return population
    
    # 计算适应度的方法
    def fitness(self, individual):
        # 调用概率潮流计算函数计算个体的适应度
        return probabilistic_power_flow(individual)
    
    # 选择父代个体的方法
    def select_parent(self, population, fitness_values):
        # 减小锦标赛规模，增加选择压力
        tournament_size = 5  # 从10减少到5
        tournament = random.sample(list(enumerate(population)), tournament_size)
        tournament_fitness = [(i, fitness_values[i]) for i, _ in tournament]
        winner_idx = min(tournament_fitness, key=lambda x: x[1])[0]
        return population[winner_idx]
    
    # 交叉操作的方法
    def crossover(self, parent1, parent2):
        if random.random() < self.crossover_rate:
            # 使用更多样的交叉方式
            alpha = random.uniform(0.2, 0.8)  # 限制权重范围，避免子代过于接近父代
            child1 = round(alpha * parent1 + (1 - alpha) * parent2, 2)
            child2 = round(alpha * parent2 + (1 - alpha) * parent1, 2)
            # 确保在[-5,5]范围内
            child1 = round(max(min(child1, 5), -5), 2)
            child2 = round(max(min(child2, 5), -5), 2)
            return child1, child2
        return parent1, parent2
    
    # 变异操作的方法
    def mutate(self, individual):
        if random.random() < self.mutation_rate:
            # 增加变异幅度
            mutation = random.gauss(0, 1.0)  # 增加标准差
            # 确保变异后的值保留两位小数
            individual = round(individual + mutation, 2)
            individual = round(max(min(individual, 5), -5), 2)
        return individual
    
    # 运行遗传算法的主方法
    def run(self):
        # 初始化种群
        population = self.initialize_population()
        # 初始化最优适应度历史记录列表
        best_fitness_history = []
        # 初始化最优解为None
        best_solution = None
        # 初始化最优适应度为正无穷
        best_fitness = float('inf')
        
        # 开始迭代，迭代n_generations次
        for generation in range(self.n_generations):
            # 在每代开始时重新计算适应度
            fitness_values = [self.fitness(ind) for ind in population]
            
            # 添加种群多样性检查
            unique_individuals = len(set(population))
            if unique_individuals < self.pop_size * 0.5:  # 如果独特个体少于50%
                # 注入新的随机个体
                num_new = int(self.pop_size * 0.2)  # 注入20%新个体
                population.extend([round(random.uniform(-5, 5), 2) for _ in range(num_new)])
                population = population[:self.pop_size]  # 保持种群大小不变
            
            # 找到当前种群中适应度最小的个体的索引
            min_fitness_idx = fitness_values.index(min(fitness_values))
            # 如果找到更好的解，更新最优解和最优适应度
            if fitness_values[min_fitness_idx] < best_fitness:
                best_fitness = fitness_values[min_fitness_idx]
                best_solution = population[min_fitness_idx]
            
            # 将当前代的最优适应度添加到历史记录中
            best_fitness_history.append(best_fitness)
            
            # 打印当前代的最优解和适应度
            print(f"第 {generation} 代: 最优PESS = {best_solution:.4f}, 适应度 = {best_fitness:.4f}")
            
            # 创建新的种群列表
            new_population = []
            
            # 设置精英个体数量为2
            elite_size = 2
            # 获取适应度最好的两个个体的索引
            elite_indices = sorted(range(len(fitness_values)), 
                                key=lambda i: fitness_values[i])[:elite_size]
            # 将精英个体添加到新种群中
            new_population.extend([population[i] for i in elite_indices])
            
            # 生成新的子代，直到种群大小达到要求
            while len(new_population) < self.pop_size:
                # 选择第一个父代
                parent1 = self.select_parent(population, fitness_values)
                # 选择第二个父代
                parent2 = self.select_parent(population, fitness_values)
                
                # 对两个父代进行交叉操作
                child1, child2 = self.crossover(parent1, parent2)
                
                # 对第一个子代进行变异操作
                child1 = self.mutate(child1)
                # 对第二个子代进行变异操作
                child2 = self.mutate(child2)
                
                # 将两个子代添加到新种群中
                new_population.extend([child1, child2])
            
            # 如果新种群大小超过要求，截取前pop_size个个体
            population = new_population[:self.pop_size]
        
        # 打印最终结果
        print("\n最终结果:")
        print(f"最优PESS值 = {best_solution:.4f}")
        print(f"最小适应度 = {best_fitness:.4f}")
        
        # 创建一个10x6大小的图形
        plt.figure(figsize=(10, 6))
        # 绘制最优适应度历史曲线
        plt.plot(best_fitness_history)
        # 设置图标题
        plt.title('进化历史')
        # 设置x轴标签
        plt.xlabel('代数')
        # 设置y轴标签
        plt.ylabel('最优适应度值')
        # 添加网格
        plt.grid(True)
        # 显示图形
        plt.show()
        
        # 返回最优解和最优适应度
        return best_solution, best_fitness

# 如果直接运行此文件
if __name__ == "__main__":
    # 设置随机数种子，确保结果可重复
    random.seed(43)
    # 设置numpy的随机数种子
    np.random.seed(42)
    
    # 创建遗传算法实例
    ga = GeneticAlgorithm()
    # 运行算法并获取结果
    best_solution, best_fitness = ga.run()