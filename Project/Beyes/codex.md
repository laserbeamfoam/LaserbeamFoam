# 贝叶斯优化算法自动迭代参数
## 核心需求
智能化地运行laserbeamFoam的计算程序和后处理程序，获取不同参数下的熔池深度和宽度，并把这些数据和实验数据的差值作为优化对象输入到贝叶斯算法中，智能获取最符合实验结果的参数可能，然后智能地修改这些参数，并再次运行。
## 基本步骤
以下步骤均是在命令行中运行
1. 运行laserbeamFoam
cp -r initial 0
blockMesh
setSolidFraction
decomposePar
mpirun -np 12 laserbeamFoam -parallel
reconstructPar
2. 运行后处理程序
conda activate meltpool-postproc
python $FOAM_USER_APPBIN/postProcessing/characterise_meltpool.py
3. 计算差值
.... 这一部分以及后续地部分之后再撰写