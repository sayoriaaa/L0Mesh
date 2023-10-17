# L0 复现

这是论文Mesh Denoising via L0 Minimization的C++复现，使用了`libigl`和`Eigen`作为依赖，便于编译。

## 编译

- 安装[libigl](https://libigl.github.io/)和[Eigen]()的源代码并将cmake下它们的路径替换为你安装的路径
- 双击`build.bat`
- 自动生成'build'文件夹，其下包含`L0min.exe`

## 运行
- 双击`run.bat`
- 在`examples`下会生成含噪声的cube网格以及L0降噪后的cube网格

你可以参考`run.bat`编写你自己需要的命令，更多帮助请查看
```
L0min.exe -h
noise.exe -h
```
目前`noise`支持均匀噪声和高斯噪声