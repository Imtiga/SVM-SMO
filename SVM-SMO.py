import numpy as np
from sympy import *

#5个alpha，初始化
Alpha = [1,0,0,0,1]
y = [1,1,1,-1,-1]
X = np.array([[1,2],[2,3],[3,3],[2,1],[3,2]])
keep = 0

# 计算拉格朗日对偶函数的值的方程
def W(alpha,X,y):
    result = 0
    for k in range(5):
        result = result + alpha[k]
    
    for i in range(5):
        for j in range(5):
            result = result - 0.5 * (alpha[i] * alpha[j] * y[i] * y[j] * (X[i].T.dot(X[j]) ))

    return result


while(keep < 2):
    for i in range(5):
        for j in range(5):
            if i != j:
                Alpha_temp = [0,0,0,0,0] 
                Alpha_copy = [0,0,0,0,0]
                for q in range(5):
                    Alpha_temp[q] = Alpha[q]
                    Alpha_copy[q] = Alpha[q]

                c = 0
                W_max = W(Alpha,X,y)

                for k in range(5):
                    if (k != i) and (k != j):
                        c = c*10 - (Alpha[k]*y[k])*10 #由于python精度有限，因此适使用*10倍 再除以10的形式计算小数
                        c = c/10
            
                # 每次筛选alpha_i和alpha_j进行优化，利用alpha_i*y_i的累积为0这个性质，使得函数变为关于alpha_i的单变量的方程
                a = Symbol("a")
                Alpha_temp[i] = a
                Alpha_temp[j] = (c - (a*y[i]))/y[j]
                D = diff(W(Alpha_temp,X,y),a)
                a0 = float(solve(D,a)[0])

                if a0 < 0:
                    a0 = 0
            
                # print(i,j,c,a0)

                Alpha[i] = a0

                Alpha[j] = (c - (a0*y[i]))
                # print(Alpha[j])
                if Alpha[j] != 0:
                    Alpha[j] = Alpha[j]*y[j]  #因为除法精度问题，所以这里需要先针对分母为零的情况，使alpha_j直接为0，而不用之后的步骤

                # print(Alpha[i],Alpha[j])
            
                W_temp = W(Alpha,X,y)

                if W_temp <= W_max:
                    Alpha = Alpha_copy #如果此次迭代结果计算出来的W函数值没有上一次的好，此次迭代忽略
            
                if (Alpha[j] < 0.0):
                    Alpha = Alpha_copy # 保证alpha全部大于0

                print(Alpha)
    
    keep = keep + 1 # 一次循环不足以解出最优解，实践证明需要循环两次以上
            

#计算omega的值
omega = 0
for q in range(5):
    omega = omega + Alpha[q]*y[q]*X[q]

print("omega*为:",omega)
print(Alpha)

#计算b*的值
b = 0
max_B = 0
min_B = 0.00001
for p in range(5):
    if y[p] == 1:
        max_B = max(max_B,omega.T.dot(X[p]))
    else:
        min_B = min(min_B,omega.T.dot(X[p]))

print(max_B,min_B)
b = -1 * (max_B + min_B)/2


print("b*为:",b)
print("斜率:",-1*float(omega[0]/omega[1]))




# # # # x = Symbol("x")
# # # # print(solve('x*k+b',x)[0]) #可以计算

# # # # a = [2,1]
# # # # a[0] = x
# # # # print(type(a[0]))


# # # for i in range(5):
# # #     for j in range(i,5):
# # #         if i != j:
# # #             print(i,j)


# # print((-0.9+1.2))
# # print((-9+12)/10)

# print((1.0*10-0.8*10)/10)




#copy from Kark.Li's SMO method
# for i in range(5):
#     for j in range(5):
#         if i != j:
#             Alpha_temp = [0,0,0,0,0]
#             for q in range(5):
#                 Alpha_temp[q] = Alpha[q]

#             if (y[i]*y[j] > 0):
#                 t = Symbol("t")
#                 Alpha_temp[i] = Alpha_temp[i] + t
#                 Alpha_temp[j] = Alpha_temp[j] - t

#                 D = diff(W(Alpha_temp,X,y),t)
#                 t_best = float(solve(D,t)[0])

#                 if ((Alpha[i] + t_best) < 0) or ((Alpha[j] - t_best) < 0):
#                     if t_best > 0:
#                         t_best = Alpha[j]
#                     else:
#                         t_best = Alpha[i]

#                 Alpha[i] = Alpha[i] + t_best
#                 Alpha[j] = Alpha[j] - t_best
            
#             else:
#                 t = Symbol("t")
#                 Alpha_temp[i] = Alpha_temp[i] + t
#                 Alpha_temp[j] = Alpha_temp[j] + t
            
#                 D = diff(W(Alpha_temp,X,y),t)
#                 t_best = float(solve(D,t)[0])

#                 if ((Alpha[i] + t_best) < 0) or ((Alpha[j] + t_best) < 0):
#                     t_best = max(-1*Alpha[i],-1*Alpha[j])

#                 Alpha[i] = Alpha[i] + t_best
#                 Alpha[j] = Alpha[j] + t_best
            
#             print(Alpha)