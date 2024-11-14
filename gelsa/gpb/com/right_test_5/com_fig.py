import matplotlib.pyplot as plt

LS_elsa = []
p_value_elsa = []
xs_elsa = []
ys_elsa = []
Len_elsa = []
delay_elsa = []

with open("result_elsa.txt", "r") as file_:
    next(file_)
    for line in file_:
        nums = list(map(float,line.split(",")))
        LS_elsa.append(nums[0])
        p_value_elsa.append(nums[1])
        xs_elsa.append(nums[2])
        ys_elsa.append(nums[3])
        Len_elsa.append(nums[4])
        delay_elsa.append(nums[5])

# 读取Gelsa数据
LS_Gelsa = []
p_value_Gelsa = []
xs_Gelsa = []
ys_Gelsa = []
Len_Gelsa = []
delay_Gelsa = []

with open("result_Gelsa.txt", "r") as file_:
    next(file_)
    for line in file_:
        nums = list(map(float,line.split(",")))
        LS_Gelsa.append(nums[0])
        p_value_Gelsa.append(nums[1])
        xs_Gelsa.append(nums[2])
        ys_Gelsa.append(nums[3])
        Len_Gelsa.append(nums[4])
        delay_Gelsa.append(nums[5])


fig, axes = plt.subplots(3, 2, figsize=(8, 10))


axes[0, 0].scatter(LS_elsa, LS_Gelsa, color='red', label='LS')
axes[0, 0].set_xlabel('LS_elsa')
axes[0, 0].set_ylabel('LS_Gelsa')
axes[0, 0].legend()
# axes[0, 0].grid(True)

# 绘制p_value图
axes[2, 0].scatter(p_value_elsa, p_value_Gelsa, color='blue', label='p_value')
axes[2, 0].set_xlabel('p_value_elsa')
axes[2, 0].set_ylabel('p_value_Gelsa')
axes[2, 0].legend()
# axes[2, 0].grid(True)

# 绘制xs图
axes[0, 1].scatter(xs_elsa, xs_Gelsa, color='green', label='xs')
axes[0, 1].set_xlabel('xs_elsa')
axes[0, 1].set_ylabel('xs_Gelsa')
axes[0, 1].legend()
# axes[0, 1].grid(True)

# 绘制ys图
axes[1, 0].scatter(ys_elsa, ys_Gelsa, color='orange', label='ys')
axes[1, 0].set_xlabel('ys_elsa')
axes[1, 0].set_ylabel('ys_Gelsa')
axes[1, 0].legend()
# axes[1, 0].grid(True)

# 绘制Len图
axes[1, 1].scatter(Len_elsa, Len_Gelsa, color='purple', label='Len')
axes[1, 1].set_xlabel('Len_elsa')
axes[1, 1].set_ylabel('Len_Gelsa')
axes[1, 1].legend()
# axes[1, 1].grid(True)

# 绘制delay图
axes[2, 1].scatter(delay_elsa, delay_Gelsa, color='cyan', label='delay')
axes[2, 1].set_xlabel('delay_elsa')
axes[2, 1].set_ylabel('delay_Gelsa')
axes[2, 1].legend()
# axes[2, 1].grid(True)

# 调整子图之间的间距
plt.tight_layout()
plt.show()
