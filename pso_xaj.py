# Program Name: NSGA-II\PSO\NSGA II
# Description: This is a python implementation of Prof. Kalyanmoy Deb's popular NSGA-II algorithm
# Author: Haris Ali Khan 
# Supervisor: Prof. Manoj Kumar Tiwari

#Importing required modules
import math
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
class Xaj:
    WUM = 20  #上层张力水容量5-20，0-50（通常值，GLUE法取值范围）
    WLM = 60  #下层张力水容量60-90,0-150
    WDM = 40  #深层张力水容量,0-100
    WM = WUM+WLM+WDM #流域平均张力水容量
    C = 0.18 #深层蒸散发系数 0.1-0.2,0-0.5
    EX = 1.5 #表层土自由水蓄水容量曲线的方次1.5，0.5-2.5
    SM = 20 #表层土自由水蓄水容量 5-45,0-200
    KG = 0.13 #表层土自由水蓄水量对地下水的出流系数0.7,0<KI+KG<1
    KI = 0.7-KG #表层土自由水蓄水量对壤中流的出流系数
    CS = 0.09 #河网蓄水消退系数,0-1
    CI = 0.7 #壤中流消退系数0.95-0.995,0.5-1
    CG = 0.98 #地下水库消退系数0.95-0.995，0.5-1
    L = 1
    K = 24
    X = 0.35
    Kc = 0.80 #蒸散发折算系数1,0-2
    IM = 0.01 #不透水面积占全流域的比值0.01-0.05,0-0.1
    B = 0.3 #张力水蓄水容量曲线的方次0.2-0.4,0-1
    F = 290.1 #
    WMM = WM*(1+B)
    # 三水源划分
    def water_division(self, S1, FR0, PE, R):
        if R == 0: #透水面积产流R
            RS = 0 #地面径流
            RG = S1 * self.KG * FR0 #地下径流
            RI = S1 * self.KI * FR0 #壤中流
            S2 = S1 * 0.3 #
        else:
            FR = R / PE #产流面积FR
            Smm = (1 + self.EX) * self.SM
            AU = Smm * (1 - (1 - S1 * FR0 / FR / self.SM) ** (1 / (1 + self.EX)))
            if PE+AU < Smm:
                RS = FR*(PE+S1*FR0/FR-self.SM+self.SM*(1-(PE+AU)/Smm)**(self.EX+1))
            else:
                RS = FR*(PE+S1*FR0/FR-self.SM)
            S = S1*FR0/FR+(R-RS)/FR
            RI = self.KI*S*FR
            RG = self.KG*S*FR
            S2 = S*0.3
        if RS < 0:
            RS = 0
        if RI < 0:
            RI = 0
        if RG < 0:
            RG = 0
        if S2 < 0:
            S2 = 0
        if S2 > self.SM:
            S2 = self.SM
        return RS, RI, RG, S2

    # 坡面汇流
    # 地面径流汇流
    def qs_calculation(self,FR,QS0,RS):
        QS = self.CS*QS0+(1-self.CS)*RS*self.F*FR/(24*3.6)
        return QS
    # 壤中流汇流
    def qi_calculation(self,FR,QI0,RI):
        QI = self.CI*QI0+(1-self.CI)*RI*self.F*FR/(24*3.6)
        return QI
    # 地下径流汇流
    def qg_calculation(self,FR,QG0,RG):
        QG = self.CG*QG0+(1-self.CG)*RG*self.F*FR/(24*3.6)
        return QG

    # 河网汇流
    def q_calculation(self,Q0,QT):
        Q = self.CS*Q0+(1-self.CS)*QT
        return Q

    # 三层蒸散发计算
    def e_calculation(self,P,WU,EP,WL):
        EU = WU + P
        if (WU + P) >= EP:
            EU = EP
            EL = 0
            ED = 0
        else:
            if WL >= self.C * self.WLM:
                EU = WU + P
                EL = (EP - EU) * WL / self.WLM
                ED = 0
            else:
                if WL >= self.C*(EP - EU):
                    EU = WU + P
                    EL = self.C * (EP - EU)
                    ED = 0
                else:
                    EU = WU + P
                    EL = WL
                    ED = self.C * (EP - EU) - EL
        E = EU + EL + ED
        return EU, EL, ED, E
    # 土壤含水量计算
    def w_calculation(self,WU,WL,WD,EU,EL,ED,P,R):
        WU = WU + P - EU - R
        WL = WL - EL
        WD = WD - ED
        if WU > self.WUM:
            WL = WU - self.WUM + WL
            WU = self.WUM
        if WL > self.WLM:
            WD = WD + WL - self.WLM
            WL = self.WLM
        if WD > self.WDM:
            WD = self.WDM
        W = WU + WL +WD
        return WU, WL, WD, W

    # 产流量计算
    def r_calculation(self,W,P,E):
        x = 1 - W / self.WM
        y = 1 / 1 + self.B
        a = self.WMM*(1-x**y)
        if P == 0:
            R = 0
        else:
            if (a + P - E) <= self.WMM:
                m = 1 - (P - E + a) / self.WMM
                n = 1+self.B
                R = P - E + W - self.WM + self.WM * m**n
            else:
                R = P - E + W - self.WM
        RR = (1 - self.IM) * R + self.IM * (P - E)
        if RR < 0:
            RR = 0
        return RR

def objFunction(X):
    # 数据读取
    df = pd.read_excel(r"D:\\广西大学资料\\博一资料整理\\现代水文模拟预报\\111.xlsx", header=None)
    arr = df.values
    P = arr[0:2922, 0] #降雨量
    E0 = arr[0:2922, 1] #蒸发量
    Q_measurement = arr[0:2922, 2] #流域出流量实测值
    l = len(P)

    # 产流计算
    # 初值设置
    EU = np.zeros(l)
    EL = np.zeros(l)
    ED = np.zeros(l)
    E = EU + EL + ED
    WU = np.zeros(l)
    WL = np.zeros(l)
    WD = np.zeros(l)
    W = np.zeros(l)
    R = np.zeros(l)
    PE = np.zeros(l)
    xaj = Xaj()
    #paramete,参数
    xaj.Kc=X[0]
    xaj.SM=X[1]
    xaj.KG=X[2]
    xaj.CS=X[3]
    xaj.L=X[4]
    EP = xaj.Kc * E0
    WU[0] = 0
    WL[0] = 30
    WD[0] = 40
    W[0] = WU[0] + WL[0] + WD[0]
    # 计算
    for i in range(1, l):
        WU[i], WL[i], WD[i], W[i] = xaj.w_calculation(WU[i-1], WL[i - 1], WD[i - 1], EU[i - 1], EL[i - 1], ED[i - 1],
                                       P[i - 1], R[i - 1])
        EU[i], EL[i], ED[i], E[i]= xaj.e_calculation(P[i], WU[i], EP[i], WL[i])
        R[i] = xaj.r_calculation(W[i - 1], P[i], E[i])
        PE[i] = P[i] - E[i]
    # 汇流计算
    # 初值设置
    RS = np.zeros(l)
    RI = np.zeros(l)
    RG = np.zeros(l)
    S = np.zeros(l)
    FR = np.zeros(l)
    QS = np.zeros(l)
    QI = np.zeros(l)
    QG = np.zeros(l)
    QT = np.zeros(l)
    Q = np.zeros(l)
    Q[0] = 10
    QS[0] = 10 #地面总入流
    QI[0] = 10 #壤中总入流
    QG[0] = 10 #地下总入流
    QT[0] = QS[0] + QI[0] + QG[0]
    FR[0] = 0.1
    # 计算
    for i in range(1, l):
        if PE[i] <= 0:
            FR[i] = 0
        else:
            FR[i] = R[i]/PE[i]
        RS[i], RI[i], RG[i], S[i] = xaj.water_division(S[i-1], FR[i-1], PE[i], R[i])
        QS[i] = xaj.qs_calculation(FR[i], QS[i-1], RS[i])
        QI[i] = xaj.qi_calculation(FR[i], QI[i-1], RI[i])
        QG[i] = xaj.qg_calculation(FR[i], QG[i-1], RG[i])
        QT[i] = QS[i] + QI[i] + QG[i]
        Q[i] = xaj.q_calculation(Q[i-1], QT[i-1])
    # 模型评价
    # 计算相对误差
    # 洪峰相对误差
    peak_err = (np.max(Q) - np.max(Q_measurement)) / np.max(Q_measurement)
    print(peak_err)
    # 洪量相对误差
    quantity_err = (np.sum(Q) - np.sum(Q_measurement)) / np.sum(Q_measurement)
    print(quantity_err)
    # 计算确定性系数
    up = np.zeros(l)
    down = np.zeros(l)
    measurement_average = np.average(Q_measurement)
    for i in range(l):
        up[i] = (Q[i] - Q_measurement[i])**2
        down[i] = (Q_measurement[i] - measurement_average)**2
    DC = 1 - np.sum(up)/np.sum(down)
    return abs(1-DC)  #限定性元素

# if __name__ == "__main__":
#     # 数据读取
#     df = pd.read_excel(r"D:\广西大学资料\博一资料整理\现代水文模拟预报\111.xlsx", header=None)
#     arr = df.values
#     P = arr[0:2922, 0] #降雨量
#     E0 = arr[0:2922, 1] #蒸发量
#     Q_measurement = arr[0:2922, 2] #流域出流量实测值
#     l = len(P)

#     # 产流计算
#     # 初值设置
#     EU = np.zeros(l)
#     EL = np.zeros(l)
#     ED = np.zeros(l)
#     E = EU + EL + ED
#     WU = np.zeros(l)
#     WL = np.zeros(l)
#     WD = np.zeros(l)
#     W = np.zeros(l)
#     R = np.zeros(l)
#     PE = np.zeros(l)
#     xaj = Xaj()
#     EP = xaj.Kc * E0
#     WU[0] = 0
#     WL[0] = 30
#     WD[0] = 40
#     W[0] = WU[0] + WL[0] + WD[0]
#     # 计算
#     for i in range(1, l):
#         WU[i], WL[i], WD[i], W[i] = xaj.w_calculation(WU[i-1], WL[i - 1], WD[i - 1], EU[i - 1], EL[i - 1], ED[i - 1],
#                                        P[i - 1], R[i - 1])
#         EU[i], EL[i], ED[i], E[i]= xaj.e_calculation(P[i], WU[i], EP[i], WL[i])
#         R[i] = xaj.r_calculation(W[i - 1], P[i], E[i])
#         PE[i] = P[i] - E[i]
#     # 汇流计算
#     # 初值设置
#     RS = np.zeros(l)
#     RI = np.zeros(l)
#     RG = np.zeros(l)
#     S = np.zeros(l)
#     FR = np.zeros(l)
#     QS = np.zeros(l)
#     QI = np.zeros(l)
#     QG = np.zeros(l)
#     QT = np.zeros(l)
#     Q = np.zeros(l)
#     Q[0] = 10
#     QS[0] = 10 #地面总入流
#     QI[0] = 10 #壤中总入流
#     QG[0] = 10 #地下总入流
#     QT[0] = QS[0] + QI[0] + QG[0]
#     FR[0] = 0.1
#     # 计算
#     for i in range(1, l):
#         if PE[i] <= 0:
#             FR[i] = 0
#         else:
#             FR[i] = R[i]/PE[i]
#         RS[i], RI[i], RG[i], S[i] = xaj.water_division(S[i-1], FR[i-1], PE[i], R[i])
#         QS[i] = xaj.qs_calculation(FR[i], QS[i-1], RS[i])
#         QI[i] = xaj.qi_calculation(FR[i], QI[i-1], RI[i])
#         QG[i] = xaj.qg_calculation(FR[i], QG[i-1], RG[i])
#         QT[i] = QS[i] + QI[i] + QG[i]
#         Q[i] = xaj.q_calculation(Q[i-1], QT[i-1])
#     # 模型评价
#     # 计算相对误差
#     # 洪峰相对误差
#     peak_err = (np.max(Q) - np.max(Q_measurement)) / np.max(Q_measurement)
#     print(peak_err)
#     # 洪量相对误差
#     quantity_err = (np.sum(Q) - np.sum(Q_measurement)) / np.sum(Q_measurement)
#     print(quantity_err)
#     # 计算确定性系数
#     up = np.zeros(l)
#     down = np.zeros(l)
#     measurement_average = np.average(Q_measurement)
#     for i in range(l):
#         up[i] = (Q[i] - Q_measurement[i])**2
#         down[i] = (Q_measurement[i] - measurement_average)**2
#     DC = 1 - np.sum(up)/np.sum(down)
#     print(DC)

class PSO(object):
    def __init__(self,particle_num,particle_dim,iter_num,c1,c2,w,max_value,min_value):
        '''参数初始化
        particle_num(int):粒子群的粒子数量
        particle_dim(int):粒子维度，对应待寻优参数的个数
        iter_num(int):最大迭代次数
        c1(float):局部学习因子，表示粒子移动到该粒子历史最优位置(pbest)的加速项的权重
        c2(float):全局学习因子，表示粒子移动到所有粒子最优位置(gbest)的加速项的权重
        w(float):惯性因子，表示粒子之前运动方向在本次方向上的惯性
        max_value(float):参数的最大值
        min_value(float):参数的最小值
        '''
        self.particle_num = particle_num
        self.particle_dim = particle_dim
        self.iter_num = iter_num
        self.c1 = c1  ##通常设为2.0
        self.c2 = c2  ##通常设为2.0
        self.w = w    
        self.max_value = max_value
        self.min_value = min_value
        
        
### 2.1 粒子群初始化
    def swarm_origin(self):
        '''粒子群初始化
        input:self(object):PSO类
        output:particle_loc(list):粒子群位置列表
               particle_dir(list):粒子群方向列表
        '''
        particle_loc = []
        particle_dir = []
        for i in range(self.particle_num):
            tmp1 = []
            tmp2 = []
            for j in range(self.particle_dim):
                a = random.random()
                b = random.random()
                tmp1.append(a * (self.max_value - self.min_value) + self.min_value)
                tmp2.append(b)
            particle_loc.append(tmp1)
            particle_dir.append(tmp2)
        
        return particle_loc,particle_dir

## 2.2 计算适应度函数数值列表;初始化pbest_parameters和gbest_parameter   
    def fitness(self,particle_loc):
        '''计算适应度函数值
        input:self(object):PSO类
              particle_loc(list):粒子群位置列表
        output:fitness_value(list):适应度函数值列表
        '''
        fitness_value = []
        ### 1.适应度函数为RBF_SVM的3_fold交叉校验平均值
        for i in range(self.particle_num):
            temp=objFunction(particle_loc[i][:])
            fitness_value.append(temp)


        ### 2. 当前粒子群最优适应度函数值和对应的参数
        current_fitness = 0.0
        current_parameter = []
        for i in range(self.particle_num):
            if current_fitness < fitness_value[i]:
                current_fitness = fitness_value[i]
                current_parameter = particle_loc[i]

        return fitness_value,current_fitness,current_parameter 
        

## 2.3  粒子位置更新 
    def updata(self,particle_loc,particle_dir,gbest_parameter,pbest_parameters):
        '''粒子群位置更新
        input:self(object):PSO类
              particle_loc(list):粒子群位置列表
              particle_dir(list):粒子群方向列表
              gbest_parameter(list):全局最优参数
              pbest_parameters(list):每个粒子的历史最优值
        output:particle_loc(list):新的粒子群位置列表
               particle_dir(list):新的粒子群方向列表
        '''
        ## 1.计算新的量子群方向和粒子群位置
        for i in range(self.particle_num): 
            a1 = [x * self.w for x in particle_dir[i]]
            a2 = [y * self.c1 * random.random() for y in list(np.array(pbest_parameters[i]) - np.array(particle_loc[i]))]
            a3 = [z * self.c2 * random.random() for z in list(np.array(gbest_parameter) - np.array(particle_dir[i]))]
            particle_dir[i] = list(np.array(a1) + np.array(a2) + np.array(a3))
#            particle_dir[i] = self.w * particle_dir[i] + self.c1 * random.random() * (pbest_parameters[i] - particle_loc[i]) + self.c2 * random.random() * (gbest_parameter - particle_dir[i])
            particle_loc[i] = list(np.array(particle_loc[i]) + np.array(particle_dir[i]))
            
        ## 2.将更新后的量子位置参数固定在[min_value,max_value]内 
        ### 2.1 每个参数的取值列表
        parameter_list = []
        for i in range(self.particle_dim):
            tmp1 = []
            for j in range(self.particle_num):
                tmp1.append(particle_loc[j][i])
            parameter_list.append(tmp1)
        ### 2.2 每个参数取值的最大值、最小值、平均值   
def function_name(parameter_list):
    value = []
    for i in range(self.particle_dim):
        tmp2 = []
        tmp2.append(max(parameter_list[i]))
        tmp2.append(min(parameter_list[i]))
        value.append(tmp2)
    
    # Add the limits for parameters X1, X2, X3, X4, X5
    value.append([0.1, 1.5])  # Parameter x1
    value.append([5, 50])  # Parameter x2
    value.append([0.01, 0.7])  # Parameter x3
    value.append([0.01, 0.9])  # Parameter x4
    value.append([0,10])  # Parameter x5
    
    for i in range(self.particle_num):
        for j in range(self.particle_dim):
            particle_loc[i][j] = (particle_loc[i][j] - value[j][1])/(value[j][0] - value[j][1]) * (self.max_value - self.min_value) + self.min_value
            
    return particle_loc, particle_dir


## 2.4 画出适应度函数值变化图
    def plot(self,results):
        '''画图
        '''
        X = []
        Y = []
        for i in range(self.iter_num):
            X.append(i + 1)
            Y.append(results[i])
        plt.plot(X,Y)
        plt.xlabel('Number of iteration',size = 15)
        plt.ylabel('Value of CV',size = 15)
        plt.title('PSO_RBF_SVM parameter optimization')
        plt.show() 
        
## 2.5 主函数        
    def main(self):
        '''主函数
        '''
        results = []
        best_fitness = 0.0 
        ## 1、粒子群初始化
        particle_loc,particle_dir = self.swarm_origin()
        ## 2、初始化gbest_parameter、pbest_parameters、fitness_value列表
        ### 2.1 gbest_parameter
        gbest_parameter = []
        for i in range(self.particle_dim):
            gbest_parameter.append(0.0)
        ### 2.2 pbest_parameters
        pbest_parameters = []
        for i in range(self.particle_num):
            tmp1 = []
            for j in range(self.particle_dim):
                tmp1.append(0.0)
            pbest_parameters.append(tmp1)
        ### 2.3 fitness_value
        fitness_value = []
        for i in range(self.particle_num):
            fitness_value.append(0.0)
    
        ## 3.迭代
        for i in range(self.iter_num):
            ### 3.1 计算当前适应度函数值列表
            current_fitness_value,current_best_fitness,current_best_parameter = self.fitness(particle_loc)
            ### 3.2 求当前的gbest_parameter、pbest_parameters和best_fitness
            for j in range(self.particle_num):
                if current_fitness_value[j] > fitness_value[j]:
                    pbest_parameters[j] = particle_loc[j]
            if current_best_fitness > best_fitness:
                best_fitness = current_best_fitness
                gbest_parameter = current_best_parameter
            
            print('iteration is :',i+1,';Best parameters:',gbest_parameter,';Best fitness',best_fitness)
            results.append(best_fitness)
            ### 3.3 更新fitness_value
            fitness_value = current_fitness_value
            ### 3.4 更新粒子群
            particle_loc,particle_dir = self.updata(particle_loc,particle_dir,gbest_parameter,pbest_parameters)
        ## 4.结果展示
        results.sort()
        self.plot(results)
        print('Final parameters are :',gbest_parameter)
            

if __name__ == '__main__':
    print('----------------1.Load Data-------------------')
    print('----------------2.Parameter Seting------------')
    particle_num = 20
    particle_dim = 6
    iter_num = 10
    c1 = 2
    c2 = 2
    w = 0.8
    max_value = 50
    min_value = 0
    print('----------------3.PSO_RBF_SVM-----------------')
    pso = PSO(particle_num,particle_dim,iter_num,c1,c2,w,max_value,min_value)
    pso.main()
