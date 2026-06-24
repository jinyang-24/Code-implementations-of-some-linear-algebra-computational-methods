def validate_results(filename):
    with open(filename, 'r') as f:
        times = [float(line.strip()) for line in f]
    
    print(f"总行数: {len(times)}")
    print(f"最小值: {min(times)}")
    print(f"最大值: {max(times)}")
    print(f"平均值: {sum(times)/len(times)}")
    
    # 检查是否单调递增（大致趋势）
    for i in range(1, len(times)):
        if times[i] < times[i-1] * 0.5:  # 允许波动，但不能下降太多
            print(f"警告: 行{i} ({times[i]}) 小于前一行的一半 ({times[i-1]})")
    
    # 检查n³增长趋势
    for n in [1, 10, 100, 500, 1000]:
        idx = n-1
        print(f"n={n}: {times[idx]:.10f} 秒")

validate_results('results1.num')