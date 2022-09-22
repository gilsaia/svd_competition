# 矩阵奇异值分解竞赛
## 准备环境需要做什么
- 安装docker
- 安装make（linux自带，mac应该可以通过homebrew安装）
- 创建data文件夹 将每个task数据分别放到data下（即格式为data/task1/data.mat,data/task2/data.mat...,不需要改名)
- 本地的python环境安装wandb包（线上看板包，用于后面传输测速结果）
## 自定义环境变量
项目提供了`.env.example` 提供一些个人环境的配置项 需要使用的话将其复制为`.env`即可 Makefile在运行时会自动读取
- `local` 为`on`时表示测试不在docker环境下运行 会直接运行对应的python文件/octave/matlab脚本
- `matlab` 为`on`时表示测试命令不使用octave而是替换为matlab **目前的测试命令是我查询推测的 可能需要进一步修改**
## 测试流程
第一次测试时由于有依赖目标时间可能会较长 之后重复文件不更新时依赖目标不会运行
```shell
make build-docker # 构建容器 正常无变动的情况下只需要运行一次
make one-matrix-task1 # 简单测试第一个任务 只运行一个矩阵 不查询正确性 不运行python脚本 检查是否有编程错误
make simple-check # 简单测试 每个任务只使用最小的矩阵测试是否满足要求
make complete-check # 尚未实现！！ 完整测试 每个任务测试所有矩阵是否满足要求
make measure # 尚未实现！！ 测速测试 预期会每个任务每个矩阵重复跑约200-1000次 计算平均时间
```
## 脚本对应关系
### task1
测试中会调用`utils/run_svd_task_with_r.m` 在该脚本中调用对应的求解函数`code/my_svd_1.m`
### task2
测试中会调用`utils/run_svd_task.m` 在该脚本中调用对应的求解函数`code/my_svd_2.m`
### task3
测试中会调用`utils/run_svd_task_dense.m` 在该脚本中调用对应的求解函数`code/my_svd_3.m`
### task4
测试中会调用`utils/run_inv_task.m` 在该脚本中调用对应的求解函数`code/my_inverse.m`

## Makefile目标解释
下列目标均可配合`make`命令使用
### build-docker
根据目录下的Dockerfile构建容器 目录下的`.dockerignore`是减小构建时上下文体积 加速构建使用的

当前环境使用ubuntu20.04 安装octave和python对应包 
### test-docker
进入容器 方便进行各项简单测试 输入`exit`退出容器 退出容器自动删除
### transpose
由于官方提供数据维度非常离谱（batch在最后） 通过简单命令手动调整了一下维度

实际运行的过程是在容器中运行了utils/gen_truth.py的代码

生成的数据会放在`data/trans/task[1234]/dataname.mat`中 会自动寻找`data`文件夹下所有存在的mat文件
### svd
**实际并不需要 没有必须的场景没必要运行** 

通过python库生成正确的结果并存储下来 后来发现测试其实不需要这部分结果
### simple-check
会分别运行每个任务的simple-check
### simple-check-task1/2/3/4
每个任务会查看transpose目标 即当不存在正确的数据源时会自动生成

之后会根据任务类别运行对应的脚本获取测试结果 

实际运行过程是在容器中运行了utils/check.py的代码
### one-matrix-task1/2/3/4
用于检测脚本是否能正确运行 会直接通过octave/matlab调用utils/one_matrix_test.m 生成一个符合要求的随机矩阵后调用指定函数 仅推荐用于检测脚本是否会产生运行错误 在调试时也可以直接通过命令行运行对应脚本