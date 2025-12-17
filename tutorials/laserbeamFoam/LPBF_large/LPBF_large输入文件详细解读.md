# LPBF_large 案例输入文件详细解读文档

## 目录
1. [案例概述](#案例概述)
2. [目录结构](#目录结构)
3. [system 目录详解](#system-目录详解)
4. [constant 目录详解](#constant-目录详解)
5. [初始条件和边界条件](#初始条件和边界条件)
6. [粉末床生成(DEM模拟)](#粉末床生成)
7. [运行流程](#运行流程)
8. [关键参数总结](#关键参数总结)

---

## 案例概述

**案例名称**: LPBF_large (激光粉末床熔融 - 大规模案例)

**物理模型**:
- 多相流VOF方法追踪金属-气体界面
- 激光热源模型(高斯热源)
- 相变模拟(固-液-气相变)
- Marangoni对流效应
- 反冲压力和蒸发
- 粉末床结构(通过DEM预生成)

**求解器**: `Flint_multiphaseEulerFoamD` (定制版多相流求解器)

**计算域尺寸**:
- X方向: 0 ~ 200 μm
- Y方向: 0 ~ 800 μm
- Z方向: 0 ~ 200 μm

**网格规模**: 80 × 320 × 80 = 2,048,000 个网格单元

**并行计算**: 12个处理器 (3×2×2 分解)

---

## 目录结构

```
LPBF_large/
├── 0/                      # 初始时刻场文件(从initial复制)
├── initial/                # 初始场模板文件
├── constant/               # 物理常数和材料属性
│   ├── polyMesh/          # 网格文件
│   ├── LaserProperties    # 激光参数
│   ├── transportProperties # 输运属性
│   ├── timeVsLaserPosition # 激光路径
│   ├── timeVsLaserPower   # 激光功率时间表
│   ├── location           # 粉末颗粒位置数据
│   └── g                  # 重力加速度
├── system/                 # 求解控制文件
│   ├── controlDict        # 时间控制
│   ├── fvSchemes          # 离散格式
│   ├── fvSolution         # 求解器设置
│   ├── blockMeshDict      # 网格生成
│   ├── decomposeParDict   # 并行分解
│   ├── setFieldsDict      # 场初始化
│   └── bedPlateDict       # 基板定义
├── DEM_large/              # DEM粉末床生成
│   └── input.liggghts     # LIGGGHTS输入文件
├── processor0-11/          # 并行计算域
├── Allrun                  # 运行脚本
└── Allclean               # 清理脚本
```

---

## system 目录详解

### 1. controlDict - 时间控制和求解器设置

**文件路径**: `system/controlDict`

#### 主要参数解析:

```cpp
application     Flint_multiphaseEulerFoamD;
```
- **求解器**: 定制的多相Euler求解器，专门用于激光熔融模拟

```cpp
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         1e-3;              // 总模拟时间: 1毫秒
```
- **时间控制**: 从最新时刻继续运行，模拟1毫秒的激光熔融过程
- 这对应于激光扫描的一个典型时间尺度

```cpp
deltaT          1.0e-8;            // 初始时间步长: 10纳秒
```
- **初始时间步长**: 极小的时间步长确保捕捉快速的相变和流动现象

```cpp
writeControl    adjustableRunTime;
writeInterval   5e-6;              // 每5微秒输出一次
```
- **输出控制**: 使用可调整运行时间控制，每5微秒保存一次结果
- 在1毫秒内会生成约200个时间步的输出

```cpp
adjustTimeStep  yes;
maxCo           0.2;               // 最大Courant数
maxAlphaCo      0.2;               // 最大界面Courant数
maxDeltaT       1.0e-6;            // 最大时间步长: 1微秒
```
- **自适应时间步长**:
  - 启用自适应时间步调整
  - Courant数限制为0.2(保守设置，确保稳定性)
  - 界面追踪的Courant数也限制为0.2
  - 最大时间步限制为1微秒(防止步长过大)

**物理意义**:
- 激光加热导致极快的温度变化(~10^6 K/s)
- 熔池流动速度可达1-10 m/s
- 需要小时间步长捕捉这些快速现象

---

### 2. fvSchemes - 数值离散格式

**文件路径**: `system/fvSchemes`

#### 离散格式详解:

```cpp
ddtSchemes
{
    default         Euler;         // 一阶隐式Euler时间格式
}
```
- **时间离散**: 一阶Euler格式，稳定但精度相对较低
- 适合处理激光熔融中的强非线性问题

```cpp
gradSchemes
{
    default         Gauss linear;  // 高斯线性梯度格式
}
```
- **梯度计算**: 二阶精度的Gauss-线性格式
- 在非结构网格上表现良好

```cpp
divSchemes
{
    div(rhoPhi,U)               Gauss linearUpwind grad(U);
    div(phi,alpha)              Gauss interfaceCompression vanLeer 1;
    div(phirb,alpha)            Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
    div(rhoPhi,epsilon1)        Gauss linear;
    div(rhophicp,T)             Gauss upwind;
}
```
- **对流项离散**:
  - **动量方程**: linearUpwind格式，具有二阶精度和良好的稳定性
  - **界面追踪**: interfaceCompression + vanLeer限制器
    - `compression factor = 1`: 适度的界面压缩，保持界面锐利
  - **温度方程**: upwind格式，确保稳定性(激光加热导致强对流)

**关键技术 - 界面压缩**:
```cpp
div(phi,alpha) Gauss interfaceCompression vanLeer 1;
```
- 这是VOF方法的核心，压缩界面防止数值扩散
- vanLeer限制器确保有界性(alpha ∈ [0,1])
- 压缩因子1表示标准压缩强度

```cpp
laplacianSchemes
{
    default         Gauss linear corrected;
}
```
- **扩散项**: Gauss线性格式 + 非正交修正
- `corrected`: 处理网格非正交性的修正项

```cpp
interpolationSchemes
{
    default         linear;        // 线性插值
}

snGradSchemes
{
    default         corrected;     // 修正的表面法向梯度
}
```

**数值稳定性考虑**:
- 温度场使用upwind: 处理激光热源引起的强对流
- 界面使用压缩格式: 保持金属-气体界面清晰
- 动量使用linearUpwind: 平衡精度和稳定性

---

### 3. fvSolution - 求解器和算法设置

**文件路径**: `system/fvSolution`

#### 线性求解器设置:

```cpp
solvers
{
    "alpha.metal.*"
    {
        nAlphaCorr      2;         // 相分数修正次数
        nAlphaSubCycles 1;         // 子循环次数
        cAlpha          1;         // 界面压缩系数
        MULESCorr       yes;       // MULES修正
        nLimiterIter    5;         // 限制器迭代次数

        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0;
    }
```
- **相分数求解器(alpha.metal)**:
  - `nAlphaCorr=2`: 每个时间步内修正2次，提高界面精度
  - `MULESCorr=yes`: 启用MULES(Multidimensional Universal Limiter)
    - 保证相分数有界性: 0 ≤ α ≤ 1
    - 保守性: 质量守恒
  - `nLimiterIter=5`: 限制器迭代5次确保有界性
  - 绝对收敛容差: 10^-8 (非常严格，确保相界面精确)

```cpp
    "pcorr.*"
    {
        solver          PCG;       // 预条件共轭梯度法
        preconditioner  DIC;       // 对角不完全Cholesky预条件
        tolerance       1e-5;
        relTol          0;
    }
```
- **压力修正求解器**: 用于PIMPLE算法的初始压力修正
- PCG + DIC: 对称正定矩阵的高效求解

```cpp
    p_rgh                          // 修正压力(p - ρgh)
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-07;
        relTol          0.05;      // 5%相对容差
    }

    p_rghFinal
    {
        $p_rgh;
        relTol          0;         // 最终求解无相对容差
    }
```
- **压力求解器**:
  - 普通迭代: 相对容差5% (加快收敛)
  - 最终迭代: 无相对容差，绝对容差10^-7
  - 使用`p_rgh`而非`p`减少浮点误差

```cpp
    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-06;
        relTol          0;
    }
```
- **速度求解器**: Gauss-Seidel光滑求解器
- 容差10^-6，对流场精度要求

```cpp
    T
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-9;      // 极高精度
        relTol          0.0;
    }
```
- **温度求解器**:
  - 容差10^-9: 温度场需要极高精度
  - 原因: 温度决定相变、物性、Marangoni对流
  - 注释显示曾考虑PBiCG求解器(非对称)

#### MELTING算法 - 相变模型:

```cpp
MELTING
{
    minTempCorrector 1;            // 最小温度修正次数
    maxTempCorrector 20;           // 最大温度修正次数
    epsilonTolerance 1e-6;         // 液相分数容差
    epsilonRelaxation 0.9;         // 松弛因子
    damperSwitch true;             // 启用阻尼器
}
```
- **相变迭代控制**:
  - 每个时间步最多迭代20次求解相变
  - `epsilonTolerance`: 液相分数(ε)收敛判据
  - `epsilonRelaxation=0.9`: 欠松弛稳定非线性迭代
  - `damperSwitch`: 阻尼器处理固-液转变的数值振荡

**物理背景**:
- 液相分数 ε 定义:
  - ε = 0: 完全固态
  - 0 < ε < 1: 糊状区
  - ε = 1: 完全液态
- 固液转变区间: Tsolidus = 1658K, Tliquidus = 1723K

#### PIMPLE算法:

```cpp
PIMPLE
{
    momentumPredictor   no;        // 不使用动量预测
    nOuterCorrectors    1;         // 外循环1次 (PISO模式)
    nCorrectors         3;         // 压力修正3次
    nNonOrthogonalCorrectors 0;    // 无非正交修正
}
```
- **压力-速度耦合**:
  - `nOuterCorrectors=1`: 实际上是PISO算法而非PIMPLE
  - `nCorrectors=3`: 每个时间步3次压力修正
  - 无非正交修正: 网格质量良好(结构网格)

```cpp
relaxationFactors
{
    equations
    {
        ".*" 1;                    // 所有方程无松弛
    }
}
```
- **松弛因子**: 所有为1，表示不使用欠松弛
- 原因: 小时间步长已经保证稳定性

---

### 4. blockMeshDict - 网格生成

**文件路径**: `system/blockMeshDict`

#### 计算域定义:

```cpp
convertToMeters 1.0;               // 单位: 米

vertices
(
    (0      0      0)              //0 - 原点
    (200e-6 0      0)              //1 - X方向 200μm
    (200e-6 800e-6 0)              //2 - XY平面对角
    (0      800e-6 0)              //3 - Y方向 800μm
    (0      0      200e-6)         //4 - Z方向 200μm
    (200e-6 0      200e-6)         //5
    (200e-6 800e-6 200e-6)         //6 - 最远角点
    (0      800e-6 200e-6)         //7
);
```
- **计算域**: 200μm × 800μm × 200μm
- **物理意义**:
  - X(宽度): 200μm - 横跨激光束直径的几倍
  - Y(长度): 800μm - 激光扫描方向
  - Z(高度): 200μm - 包含粉末层和部分基板

```cpp
blocks
(
    hex (0 1 2 3 4 5 6 7) (80 320 80) simpleGrading (1 1 1)
);
```
- **网格划分**:
  - X方向: 80个网格 → 网格尺寸 2.5μm
  - Y方向: 320个网格 → 网格尺寸 2.5μm
  - Z方向: 80个网格 → 网格尺寸 2.5μm
  - 均匀网格: `simpleGrading (1 1 1)`

- **网格分辨率考虑**:
  - 激光束半径: 25μm
  - 每个激光束直径: 约20个网格
  - 粉末颗粒直径: 9-81μm
  - 可以分辨单个粉末颗粒

#### 边界定义:

```cpp
boundary
(
    back            { type wall;  }    // X=0平面
    front           { type wall;  }    // X=200μm平面
    leftWall        { type wall;  }    // Y=0平面
    rightWall       { type patch; }    // Y=800μm平面 (开放边界)
    topWall         { type patch; }    // Z=200μm平面 (激光入射面)
    bottomWall      { type patch; }    // Z=0平面 (基板)
)
```
- **边界类型**:
  - `wall`: 物理壁面(无流动穿过)
  - `patch`: 一般边界(可设置各种BC)

- **物理边界**:
  - `topWall`: 自由表面，激光从此入射
  - `bottomWall`: 基板
  - `rightWall`: 出口/开放边界
  - 其他: 对称或壁面边界

---

### 5. decomposeParDict - 并行分解

**文件路径**: `system/decomposeParDict`

```cpp
numberOfSubdomains 12;             // 12个计算域

method          simple;            // 简单几何分解

simpleCoeffs
{
    n           ( 3 2 2 );         // X×Y×Z方向分解
    delta       0.001;
}
```
- **并行策略**:
  - 总核心数: 12
  - 分解方式: 3×2×2
    - X方向: 3个子域(每个~67μm)
    - Y方向: 2个子域(每个400μm)
    - Z方向: 2个子域(每个100μm)

- **负载均衡**:
  - 每个子域: 2,048,000 / 12 ≈ 170,667 个网格
  - 简单分解保证负载均衡

---

### 6. setFieldsDict - 场初始化

**文件路径**: `system/setFieldsDict`

```cpp
defaultFieldValues
(
    volScalarFieldValue alpha.metal 0    // 默认为气相
);

regions
(
    boxToCell
    {
        box (-0.85 0.001 -0.5) (0.85 0.004 0.5);
        fieldValues
        (
            volScalarFieldValue alpha.metal 1    // 设为金属
        );
    }
```
- **初始化策略**:
  - 整个计算域默认为气相(α=0)
  - 特定区域设为金属(α=1)
  - 注意: 实际粉末床通过读取`location`文件设置

- **注释部分**: 显示了其他可能的初始化方法
  - `cylinderToCell`: 圆柱形区域
  - 可用于设置轨道扫描路径的预熔金属

---

### 7. bedPlateDict - 基板定义

**文件路径**: `system/bedPlateDict`

```cpp
Bed true;                          // 启用基板

xmin 0.0;
xmax 0.0002;                       // 200μm

ymin 0.0;
ymax 0.0008;                       // 800μm

zmin 0.0;
zmax 0.0001;                       // 100μm (基板高度)
```
- **基板范围**:
  - 占据计算域底部100μm
  - 与网格Z方向一半重合

- **物理意义**:
  - 基板提供热沉效应
  - 影响熔池冷却速率
  - 粉末层在基板之上

---

## constant 目录详解

### 1. LaserProperties - 激光参数

**文件路径**: `constant/LaserProperties`

```cpp
radialPolarHeatSource no;          // 不使用径向极化热源
```
- 使用标准Gaussian热源而非复杂的极化模型

#### 激光路径和功率:

```cpp
timeVsLaserPosition
{
    file    "$FOAM_CASE/constant/timeVsLaserPosition";
    outOfBounds clamp;             // 超出范围时固定在边界值
}

timeVsLaserPower
{
    file    "$FOAM_CASE/constant/timeVsLaserPower";
    outOfBounds clamp;
}
```
- **查表方法**: 从文件读取激光位置和功率随时间变化
- `clamp`: 插值时钳位在范围内

#### 激光光学参数:

```cpp
V_incident (0 1 0);                // 入射方向: Y轴正方向(垂直向下)
laserRadius 25e-6;                 // 激光束半径: 25μm
N_sub_divisions 1;                 // 光束细分数
```
- **束斑尺寸**:
  - 半径25μm → 直径50μm
  - 典型LPBF工艺参数

```cpp
wavelength 1.064e-6;               // 波长: 1.064μm (Nd:YAG/光纤激光器)
e_num_density 5.83e29;             // 电子数密度(m^-3)
Radius_Flavour 2.0;                // 味半径因子
```
- **波长**: 1064nm是常用的固体激光器波长
- **电子密度**: 用于Fresnel吸收率计算

```cpp
PowderSim true;                    // 启用粉末模拟模式
```
- 考虑粉末颗粒的多重反射和散射

---

### 2. transportProperties - 输运和相性质

**文件路径**: `constant/transportProperties`

#### 界面追踪方法:

```cpp
interfaceTrackingScheme isoAdvector;    // 几何VOF方法
phases (metal gas);                     // 两相: 金属和气体
```
- **isoAdvector**:
  - 几何VOF方法，保持界面锐利
  - 比代数MULES方法更精确
  - 质量守恒性好

#### 金属相属性:

```cpp
metal
{
    transportModel      Newtonian;
    nu                  5e-07;         // 运动粘度: 0.5 mm²/s
    rho                 8000;          // 密度: 8000 kg/m³
```
- **流体模型**: 牛顿流体
- **物性参数**:
  - 密度8000 kg/m³: 典型铁基合金
  - 粘度0.5 mm²/s: 液态金属的典型值

**热物理性质**:

```cpp
    Tsolidus            1658;          // 固相线温度: 1658 K
    Tliquidus           1723;          // 液相线温度: 1723 K
    LatentHeat          2.7e5;         // 潜热: 270 kJ/kg
```
- **相变温度**:
  - 固相线1658K (1385°C): 开始熔化
  - 液相线1723K (1450°C): 完全熔化
  - ΔT = 65K: 糊状区温度范围
- **潜热**: 270 kJ/kg (钢的典型值~247 kJ/kg)

```cpp
    beta                5.0e-6;        // 热膨胀系数: 5×10⁻⁶ K⁻¹
```
- 密度随温度变化: ρ = ρ₀[1 - β(T - T₀)]
- 影响浮力驱动流动

**温度相关物性(多项式)**:

```cpp
    poly_kappa  (10 0.015 0 0 0 0 0 0);    // 热导率多项式系数
    poly_cp     (520 0.075 0 0 0 0 0 0);   // 比热容多项式系数
```
- **热导率**: k(T) = 10 + 0.015×T (W/m·K)
  - 室温(300K): k ≈ 14.5 W/m·K
  - 熔点(1700K): k ≈ 35.5 W/m·K

- **比热容**: cp(T) = 520 + 0.075×T (J/kg·K)
  - 室温: cp ≈ 542 J/kg·K
  - 熔点: cp ≈ 647 J/kg·K

```cpp
    elec_resistivity    1.0e-6;        // 电阻率: 1 μΩ·m
}
```
- 用于计算激光吸收率(Fresnel公式)

#### 气相属性:

```cpp
gas
{
    transportModel      Newtonian;
    nu                  1.48e-05;      // 运动粘度: 14.8 mm²/s (空气)
    rho                 1;             // 密度: 1 kg/m³
    Tsolidus            1.0;
    Tliquidus           10.0;
    LatentHeat          1;
    beta                4.0e-5;        // 热膨胀系数
    poly_kappa          (0.04 0.0 0 0 0 0 0 0);     // 热导率~0.04 W/m·K
    poly_cp             (520 0.0 0 0 0 0 0 0);      // 比热~520 J/kg·K
}
```
- 气相性质设为空气的近似值
- 密度1 kg/m³: 简化处理，实际应为~1.2 kg/m³

#### 界面和相变性质:

```cpp
sigma                   0.07;          // 表面张力: 0.07 N/m
Marangoni_Constant      -0.5e-4;       // Marangoni系数: -5×10⁻⁵ N/m·K
```
- **表面张力**: 1.8 N/m是液态钢的典型值
- **Marangoni效应**: dσ/dT = -5×10⁻⁵ N/m·K
  - 负值表示温度升高时表面张力降低
  - 驱动表面流动从热区到冷区
  - 这是熔池对流的主要驱动力之一

**蒸发参数**:

```cpp
p0                      100000.0;      // 参考压力: 1 atm
Tvap                    3068.0;        // 蒸发温度: 3068 K
Mm                      5.58e-2;       // 摩尔质量: 55.8 g/mol (Fe)
LatentHeatVap           7.45e6;        // 蒸发潜热: 7.45 MJ/kg
```
- **蒸发温度**: 3068K (2795°C) - 铁的沸点
- **蒸发潜热**: 7.45 MJ/kg
- **反冲压力**: P_recoil = 0.54×P_sat
  - 蒸发产生的压力推动熔池形成凹陷(keyhole)

---

### 3. turbulenceProperties - 湍流模型

**文件路径**: `constant/turbulenceProperties`

```cpp
simulationType  laminar;               // 层流模拟
```
- **不使用湍流模型**的原因:
  1. 熔池尺寸小(~100μm)，雷诺数不高
  2. 主要流动驱动力是Marangoni效应(层流)
  3. 直接数值模拟捕捉所有流动特征
  4. 减少计算复杂度

---

### 4. g - 重力

**文件路径**: `constant/g`

```cpp
dimensions  [0 1 -2 0 0 0 0];          // 加速度量纲
value       (0 9.81 0);                // Y方向: +9.81 m/s²
```
- **重力方向**: Y轴正方向(垂直向上)
- **坐标系**:
  - X: 横向
  - Y: 垂直(重力方向)
  - Z: 激光扫描方向

---

### 5. timeVsLaserPosition - 激光路径

**文件路径**: `constant/timeVsLaserPosition`

```cpp
(
    (0          (100e-6 20e-6 100e-6))     // t=0: 起始位置
    (600e-6     (100e-6 20e-6 700e-6))     // t=600μs: 扫描至终点
    (1000e-6    (100e-6 20e-6 700e-6))     // t=1000μs: 保持
)
```
- **激光轨迹**:
  - t = 0: 位置(100, 20, 100)μm
  - t = 600μs: 位置(100, 20, 700)μm
  - 沿Z方向扫描600μm
  - 扫描速度: 600μm / 600μs = 1 m/s

- **空间位置**:
  - X固定在100μm (域中心)
  - Y固定在20μm (靠近表面)
  - Z从100μm扫描到700μm

---

### 6. timeVsLaserPower - 激光功率

**文件路径**: `constant/timeVsLaserPower`

```cpp
(
    (0              10)                // t=0: 10W (预热)
    (1e-8           150)               // t=10ns: 150W (开启)
    (600e-6         150)               // t=600μs: 保持150W
    (600.001e-6     0)                 // t=600.001μs: 关闭
    (1000e-6        0)                 // t=1000μs: 关闭
)
```
- **功率控制**:
  - 初始10W: 极短的预热
  - 主加工: 150W持续600μs
  - 快速关闭(1ns切换)

- **工艺参数**:
  - 功率: 150W
  - 速度: 1 m/s
  - 线能量密度: E = P/v = 150 J/m
  - 体能量密度: E_v = P/(v×h×d) (需要层厚和间距)

---

### 7. location - 粉末颗粒位置

**文件路径**: `constant/location`

这个文件由DEM模拟生成，包含:
- 每个粉末颗粒的位置(x, y, z)
- 颗粒半径
- 用于初始化alpha.metal场

**粉末床特征**(从DEM参数推断):
- 颗粒尺寸分布: 9-81μm (18种尺寸)
- 颗粒密度: 4.43 g/cm³
- 粉末床孔隙率: ~40-50% (典型值)
- 粉末层厚度: ~50-100μm

---

## 初始条件和边界条件

### 初始场文件(0目录)

#### 1. T - 温度场

**文件路径**: `0/T` (从`initial/T`复制)

```cpp
dimensions      [0 0 0 1 0 0 0];       // 温度量纲(K)
internalField   uniform 300.0;         // 初始温度: 300K (室温)
```

**边界条件**:
```cpp
boundaryField
{
    back/front/leftWall/rightWall/bottomWall
    {
        type    zeroGradient;          // 绝热边界
    }

    topWall
    {
        type    zeroGradient;          // 自由表面绝热
    }
}
```
- **所有边界**: zeroGradient(绝热)
  - 假设计算域足够大，边界上热损失可忽略
  - 实际物理: 对流和辐射散热通过源项处理

**物理考虑**:
- 室温预热: 300K
- 激光加热将局部温度提升至3000K以上
- 温度梯度: ΔT/Δx ~ 10⁷ K/m

---

#### 2. U - 速度场

**文件路径**: `0/U`

```cpp
dimensions      [0 1 -1 0 0 0 0];      // 速度量纲(m/s)
internalField   uniform (0 0 0);       // 初始静止
```

**边界条件**:
```cpp
boundaryField
{
    leftWall/rightWall/bottomWall/front/back
    {
        type    fixedValue;
        value   uniform (0 0 0);       // 无滑移壁面
    }

    topWall
    {
        type    pressureInletOutletVelocity;
        phi     phi;
        rho     rho;
        value   uniform (0 0 0);       // 压力出口速度
    }
}
```
- **壁面**: 无滑移边界(v=0)
- **顶部**: 压力出口
  - 允许蒸发的金属蒸汽流出
  - 根据压力梯度自动调整速度

**流动特征**:
- Marangoni流: 表面张力梯度驱动
- 浮力流: 密度梯度驱动
- 反冲压力: 蒸发驱动
- 典型速度: 0.1-10 m/s

---

#### 3. p_rgh - 修正压力

**文件路径**: `0/p_rgh`

```cpp
dimensions      [1 -1 -2 0 0 0 0];     // 压力量纲(Pa)
internalField   uniform 0;              // 相对压力为0
```

**为什么使用p_rgh**:
- `p_rgh = p - ρgh`: 静水压力修正
- 减少数值误差，特别是密度差异大的多相流
- 改善求解器收敛性

**边界条件**:
```cpp
boundaryField
{
    back/front/leftWall/rightWall/bottomWall
    {
        type    fixedFluxPressure;     // 固定压力梯度
        value   uniform 0;
    }

    topWall
    {
        type    totalPressure;         // 总压边界
        p0      uniform 0;             // 大气压(相对)
        U       U;
        phi     phi;
    }
}
```
- **壁面**: fixedFluxPressure - 保证压力梯度与速度一致
- **顶部**: 大气压边界

---

#### 4. alpha.metal - 金属相分数

**文件路径**: `0/alpha.metal`

这个文件实际从DEM生成的粉末床数据读取，包含:
- 粉末颗粒区域: α = 1(考虑孔隙率)
- 气体区域: α = 0
- 约2,048,000个网格点的非均匀分布

**边界条件**:
```cpp
boundaryField
{
    所有边界
    {
        type    zeroGradient;          // 相分数梯度为0
    }
}
```

---

#### 5. TRHS - 温度方程右端项

**文件路径**: `0/TRHS`

```cpp
dimensions      [1 -1 -3 0 0 0 0];     // 功率密度量纲(W/m³)
internalField   uniform 0;              // 初始无热源
```
- 在求解过程中由激光模型更新
- 包含: 激光吸收、相变潜热、表面散热等

---

#### 6. Laser_boundary - 激光边界标识

**文件路径**: `0/Laser_boundary`

```cpp
dimensions      [0 0 0 0 0 0 0];       // 无量纲
internalField   uniform 0;

boundaryField
{
    topWall
    {
        type    fixedValue;
        value   uniform 1.0;           // 激光入射面标记为1
    }

    其他边界
    {
        type    fixedValue;
        value   uniform -1.0;          // 其他面标记为-1
    }
}
```
- **标识符作用**:
  - 1.0: 激光可照射的表面
  - -1.0: 激光不可照射
- 用于激光吸收计算和射线追踪

---

## 粉末床生成

### DEM模拟(LIGGGHTS)

**文件路径**: `DEM_large/input.liggghts`

#### 模拟设置:

```liggghts
atom_style      granular               // 颗粒模型
units           cgs                     // 厘米-克-秒单位
newton          off                     // 简化牛顿第三定律

region  domain block 0.0 0.02 0.0 0.08 0.0 0.05 units box
```
- **计算域**: 200μm × 800μm × 500μm (比CFD域大)
- **单位**: cgs系统(需转换到SI)

#### 材料属性:

```liggghts
fix m1 all property/global youngsModulus peratomtype 5e7        // 杨氏模量: 50 MPa
fix m2 all property/global poissonsRatio peratomtype 0.45       // 泊松比: 0.45
fix m3 all property/global coefficientRestitution peratomtypepair 1 0.1    // 恢复系数
fix m4 all property/global coefficientFriction peratomtypepair 1 0.065     // 摩擦系数
```
- **接触模型**: Hertz无粘接模型
- **摩擦系数**: 0.065 (金属颗粒间)

#### 时间步长:

```liggghts
timestep    0.00000005                 // 50 ns
```
- DEM时间步长远小于CFD(需要分辨颗粒碰撞)

#### 颗粒生成:

定义了18种颗粒模板(pts1-pts18):
```liggghts
fix pts1 all particletemplate/sphere ... radius constant 0.0004625    # ~9.25μm
...
fix pts18 all particletemplate/sphere ... radius constant 0.0040505   # ~81μm
```

**粒径分布** (概率分布):
```liggghts
fix pdd all particledistribution/discrete 78593 18
    pts1  1.97883E-06      # <0.001%
    pts2  0.000158863      # 0.016%
    pts3  0.002123547      # 0.21%
    ...
    pts10 0.138350101      # 13.8% (峰值)
    ...
    pts18 5.1597E-06       # <0.001%
```
- **D50**: 约29μm (中位粒径)
- **分布**: 接近正态分布
- **粉末规格**: 典型LPBF粉末15-53μm范围

#### 颗粒插入:

```liggghts
region factory block 0.0 0.02 0.0 0.08 0.02 0.05 units box
fix ins all insert/rate/region seed 51869
    distributiontemplate pdd
    nparticles 2000                # 总共2000个颗粒
    particlerate 200000            # 插入速率
    insert_every 5                 # 每5步插入
    overlapcheck yes               # 检查重叠
    vel constant 0. 0. 0.0         # 初速度为0
    region factory                 # 在factory区域插入
    ntry_mc 10000                  # 蒙特卡洛尝试次数
```

#### 边界:

```liggghts
fix box all mesh/surface file meshes/domain.stl type 1 scale 0.1
fix plate all mesh/surface file meshes/plate.stl type 1 scale 0.1
fix wall all wall/gran model hertz tangential history mesh n_meshes 2 meshes box plate
```
- 从STL文件读取边界几何
- `scale 0.1`: 从mm转换到cm(cgs单位)

#### 重力:

```liggghts
fix grav all gravity 981 vector 0.0 0.0 -1.0    // 981 cm/s² = 9.81 m/s²
```
- Z方向负向: 颗粒在重力下落

#### 运行序列:

```liggghts
run 1                              # 初始化
run 500000 upto                    # 运行到t=50万步 = 0.025s
                                   # 颗粒沉降形成粉末床

# 删除过高的颗粒
variable co atom "z+c_1 > 0.015"   # z+半径 > 150μm
group layer variable co
delete_atoms group layer           # 删除这些颗粒

run 750000 upto                    # 继续运行到0.0375s
                                   # 粉末床稳定
```

**输出**:
```liggghts
dump mydump all custom 750000 post/location v_x1 v_y1 v_z1 v_rad1
```
- 输出颗粒位置和半径到`post/location`
- 这个文件被CFD读取用于初始化alpha.metal场

#### 粉末床特征:

通过DEM生成的粉末床具有:
1. **真实的颗粒堆积**: 随机密堆积
2. **孔隙率**: 约40-50% (通过颗粒重叠和空隙)
3. **粒径分布**: 符合实际粉末规格
4. **表面粗糙度**: 颗粒级别的粗糙度
5. **层厚控制**: 通过删除过高颗粒实现

---

## 运行流程

### Allrun脚本解析

**文件路径**: `Allrun`

```bash
#!/bin/sh

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions
```
- 加载OpenFOAM运行函数库

```bash
echo "Copying 'initial' to 0"
cp -r initial 0
```
- **步骤1**: 复制初始场模板到0目录
- 原因: 保持initial作为干净模板

```bash
runApplication blockMesh
```
- **步骤2**: 生成网格
  - 读取`system/blockMeshDict`
  - 生成80×320×80结构网格
  - 输出到`constant/polyMesh/`

```bash
runApplication setSolidFraction
```
- **步骤3**: 设置粉末床
  - 读取`constant/location`(DEM生成)
  - 将颗粒位置映射到alpha.metal场
  - 考虑颗粒半径和孔隙率
  - 生成真实的粉末床结构

```bash
runApplication transformPoints -rotate '((0 1 0) (0 0 1))'
```
- **步骤4**: 坐标转换
  - 旋转网格: Y轴 → Z轴
  - 确保重力方向正确
  - 激光扫描方向对齐

```bash
runApplication laserbeamFoam
```
- **步骤5**: 运行CFD求解器
  - 串行运行 或
  - 如果已分解使用并行运行

```bash
# run DEM manually before the simulation if needed
```
- **注释**: DEM需要单独预先运行

### 完整运行流程:

1. **预处理(DEM)**:
   ```bash
   cd DEM_large
   liggghts < input.liggghts
   # 生成post/location文件
   cp post/location ../constant/
   ```

2. **预处理(CFD)**:
   ```bash
   ./Allrun                       # 或手动执行各步骤
   ```

3. **并行计算**:
   ```bash
   decomposePar                   # 分解域
   mpirun -np 12 laserbeamFoam -parallel
   ```

4. **后处理**:
   ```bash
   reconstructPar                 # 重建结果
   paraFoam                       # 可视化
   ```

---

## 关键参数总结

### 几何和网格

| 参数 | 数值 | 说明 |
|------|------|------|
| 计算域尺寸 | 200×800×200 μm³ | X×Y×Z |
| 网格数量 | 80×320×80 | 2,048,000 网格 |
| 网格尺寸 | 2.5 μm | 均匀网格 |
| 并行核心数 | 12 | 3×2×2分解 |

### 时间参数

| 参数 | 数值 | 说明 |
|------|------|------|
| 总模拟时间 | 1 ms | 激光扫描时间 |
| 初始时间步 | 10 ns | Δt₀ |
| 最大时间步 | 1 μs | Δt_max |
| 最大Courant数 | 0.2 | CFL条件 |
| 输出间隔 | 5 μs | 200个输出文件 |

### 激光参数

| 参数 | 数值 | 说明 |
|------|------|------|
| 激光功率 | 150 W | 加工功率 |
| 束斑半径 | 25 μm | 高斯半径 |
| 扫描速度 | 1 m/s | 沿Z方向 |
| 扫描长度 | 600 μm | 轨道长度 |
| 线能量密度 | 150 J/m | P/v |
| 波长 | 1.064 μm | Nd:YAG |

### 材料属性(金属)

| 参数 | 数值 | 说明 |
|------|------|------|
| 密度 | 8000 kg/m³ | 铁基合金 |
| 固相线 | 1658 K | 开始熔化 |
| 液相线 | 1723 K | 完全熔化 |
| 熔化潜热 | 270 kJ/kg | 固-液相变 |
| 蒸发温度 | 3068 K | 沸点 |
| 蒸发潜热 | 7.45 MJ/kg | 液-气相变 |
| 表面张力 | 1.8 N/m | 液态金属 |
| Marangoni系数 | -5×10⁻⁵ N/m·K | dσ/dT |
| 热导率(室温) | ~14.5 W/m·K | 温度相关 |
| 比热容(室温) | ~542 J/kg·K | 温度相关 |

### 粉末特性

| 参数 | 数值 | 说明 |
|------|------|------|
| 颗粒数量 | 2000 | DEM模拟 |
| 粒径范围 | 9-81 μm | 18种尺寸 |
| 中位粒径D50 | ~29 μm | 正态分布 |
| 颗粒密度 | 4.43 g/cm³ | 金属粉末 |
| 粉末层厚 | ~50-100 μm | 单层 |
| 孔隙率 | ~40-50% | 堆积密度 |

### 数值方法

| 项目 | 方法 | 说明 |
|------|------|------|
| 时间格式 | Euler隐式 | 一阶 |
| 界面追踪 | isoAdvector | 几何VOF |
| 梯度格式 | Gauss linear | 二阶 |
| 动量对流 | linearUpwind | 二阶迎风 |
| 界面对流 | interfaceCompression | 界面锐化 |
| 温度对流 | upwind | 一阶迎风(稳定) |
| 压力求解器 | PCG + DIC | 共轭梯度 |
| 速度求解器 | smoothSolver | Gauss-Seidel |
| 压力-速度耦合 | PISO | 3次修正 |

### 收敛容差

| 场变量 | 容差 | 说明 |
|--------|------|------|
| alpha.metal | 10⁻⁸ | 相分数 |
| p_rgh | 10⁻⁷ | 压力 |
| U | 10⁻⁶ | 速度 |
| T | 10⁻⁹ | 温度(最严格) |
| 液相分数ε | 10⁻⁶ | 相变迭代 |

---

## 物理现象和模型

### 1. 多相流(VOF方法)

**控制方程**:
```
∂α/∂t + ∇·(αU) = 0
```
- α: 金属相体积分数
- 追踪金属-气体界面

**界面压缩**:
```
∂α/∂t + ∇·(αU) + ∇·(α(1-α)U_c) = 0
```
- U_c: 人工压缩速度
- 保持界面锐利

### 2. 激光热源

**高斯热源模型**:
```
Q(r) = (2AP)/(πR²) exp(-2r²/R²)
```
- A: 吸收率(Fresnel公式计算)
- P: 激光功率(150W)
- R: 束斑半径(25μm)
- r: 到光束中心距离

**Fresnel吸收率**:
```
A = f(θ, λ, ρ_e, T)
```
- θ: 入射角
- λ: 波长(1.064μm)
- ρ_e: 电阻率
- T: 温度

### 3. 相变模型

**固-液相变(焓法)**:
```
ε = {
    0,                          T < T_solidus
    (T-T_s)/(T_l-T_s),         T_solidus ≤ T ≤ T_liquidus
    1,                          T > T_liquidus
}
```
- ε: 液相分数
- 潜热通过焓包含

**液-气相变(蒸发)**:
```
ṁ_evap = C √(M/(2πRT)) P_sat(T)
```
- C: 常数
- M: 摩尔质量(55.8 g/mol)
- P_sat: 饱和蒸汽压(Clausius-Clapeyron)

### 4. Marangoni效应

**表面张力随温度变化**:
```
σ(T) = σ₀ + (dσ/dT)(T - T₀)
```
- dσ/dT = -5×10⁻⁵ N/m·K

**Marangoni应力**:
```
τ_s = ∇_s σ = (dσ/dT)∇_s T
```
- 驱动从热区到冷区的表面流动

### 5. 反冲压力

**蒸发引起的压力**:
```
P_recoil = 0.54 P_sat(T)
```
- 推动熔池向下凹陷
- keyhole形成的主要机制

### 6. 动量方程

**考虑相变的NS方程**:
```
∂(ρU)/∂t + ∇·(ρUU) = -∇p + ∇·τ + ρg + F_σ + F_Marangoni + F_recoil + F_Darcy
```
其中:
- F_σ: 表面张力
- F_Marangoni: Marangoni力
- F_recoil: 反冲压力
- F_Darcy: 糊状区阻力(固-液)

### 7. 能量方程

**带相变的能量方程**:
```
∂(ρc_pT)/∂t + ∇·(ρc_pUT) = ∇·(k∇T) + Q_laser - Q_evap - Q_conv - Q_rad
```
其中:
- Q_laser: 激光吸收
- Q_evap: 蒸发热损失
- Q_conv: 对流散热
- Q_rad: 辐射散热(Stefan-Boltzmann)

---

## 典型输出和后处理

### 输出文件

模拟完成后生成约200个时间目录:
```
0/
2.080279358e-05/
2.5e-05/
...
0.0001/
0.000105/
...
0.001/
```

每个时间目录包含:
- `alpha.metal`: 相分数场
- `T`: 温度场
- `U`: 速度场
- `p_rgh`: 压力场
- 其他计算场

### 关键观测量

1. **熔池尺寸**:
   - 深度: 通过等温线(T=Tliquidus)确定
   - 宽度: 横向熔化范围
   - 长度: 纵向熔化范围

2. **温度分布**:
   - 峰值温度: 可能超过3000K(蒸发温度)
   - 温度梯度: dT/dx ~ 10⁷ K/m
   - 冷却速率: dT/dt ~ 10⁶ K/s

3. **流场特征**:
   - Marangoni涡流: 表面从热到冷
   - 回流: 熔池内部循环
   - 蒸发流: 垂直向上的蒸汽流

4. **界面演化**:
   - 自由表面变形
   - Keyhole形成和闭合
   - 表面波动

### ParaView可视化建议

1. **温度场**:
   ```
   - 使用"温度-彩虹"色标
   - 显示等温面(Tliquidus, Tvap)
   - 提取熔池轮廓
   ```

2. **流场**:
   ```
   - 速度矢量(箭头或流线)
   - 在液相区域(α>0.5, T>Tliquidus)显示
   - 涡量等值面
   ```

3. **相分数**:
   ```
   - 等值面α=0.5显示界面
   - 体渲染显示粉末结构
   ```

4. **时间序列**:
   ```
   - 动画显示熔池演化
   - 跟踪激光移动
   - 观察keyhole动态
   ```

---

## 常见问题和调试

### 1. 数值不稳定

**现象**: 计算发散，Courant数爆炸

**可能原因**:
- 时间步长过大
- 激光功率过高导致极端温度梯度
- 相变区域数值振荡

**解决方案**:
```cpp
// controlDict
maxDeltaT       5e-7;      // 减小最大时间步
maxCo           0.1;       // 降低Courant数限制

// fvSolution/MELTING
epsilonRelaxation 0.8;     // 增强欠松弛
damperSwitch true;         // 启用阻尼
```

### 2. 界面模糊

**现象**: 金属-气体界面扩散，不清晰

**解决方案**:
```cpp
// fvSchemes
div(phi,alpha) Gauss interfaceCompression vanLeer 2;  // 增大压缩因子

// fvSolution
nAlphaCorr      3;         // 增加修正次数
cAlpha          1.5;       // 增大压缩系数
```

### 3. 计算过慢

**现象**: 每个时间步耗时过长

**优化策略**:
```cpp
// fvSolution
// 降低求解器容差(在精度可接受范围内)
p_rgh { tolerance 1e-6; }
U     { tolerance 1e-5; }

// 减少相变迭代
maxTempCorrector 15;

// 增加并行核心数
// decomposeParDict
numberOfSubdomains 24;     // 例如改为24核
```

### 4. 内存不足

**现象**: 内存溢出，计算中断

**解决方案**:
```cpp
// controlDict
writeControl    adjustableRunTime;
writeInterval   1e-5;      // 增大输出间隔
purgeWrite      3;         // 只保留最近3个时间步

// 或减小网格数量(降低精度)
blocks (60 240 60);        // 从80×320×80降低
```

---

## 进一步修改和扩展

### 1. 多道扫描

修改`timeVsLaserPosition`实现蛇形扫描:
```cpp
(
    (0          (100e-6 20e-6 100e-6))
    (600e-6     (100e-6 20e-6 700e-6))  // 第一道
    (700e-6     (150e-6 20e-6 700e-6))  // 移动到下一道
    (1300e-6    (150e-6 20e-6 100e-6))  // 反向扫描
    ...
)
```

### 2. 不同材料

修改`transportProperties`:
```cpp
metal
{
    // 例如: 钛合金Ti-6Al-4V
    rho         4430;
    Tsolidus    1878;
    Tliquidus   1928;
    LatentHeat  2.86e5;
    ...
}
```

### 3. 变功率扫描

修改`timeVsLaserPower`:
```cpp
(
    (0          100)         // 低功率预热
    (100e-6     200)         // 增加到200W
    (400e-6     200)         // 保持
    (500e-6     100)         // 降低功率
    (600e-6     0)           // 关闭
)
```

### 4. 细化网格

局部加密熔池区域:
```cpp
// blockMeshDict
blocks
(
    hex (0 1 2 3 4 5 6 7) (120 480 120)  // 更细网格
    simpleGrading
    (
        1                    // X方向均匀
        (                    // Y方向分段加密
            (0.2 0.3 4)      // 前20%用30%网格,扩展比4
            (0.6 0.4 1)      // 中60%用40%网格,均匀
            (0.2 0.3 0.25)   // 后20%用30%网格,收缩
        )
        1                    // Z方向均匀
    )
);
```

### 5. 添加保护气流

在`0/U`中设置初始速度:
```cpp
internalField   uniform (0 -0.5 0);  // Y方向-0.5m/s的保护气流
```

---

## 参考文献和资源

### 物理模型参考

1. **VOF方法**:
   - Hirt, C.W. & Nichols, B.D. (1981). Volume of Fluid (VOF) Method

2. **激光熔融建模**:
   - Khairallah, S.A. et al. (2016). Laser powder-bed fusion additive manufacturing

3. **Marangoni对流**:
   - Mills, K.C. et al. (2006). Marangoni effects in welding

### OpenFOAM资源

- **官方文档**: https://www.openfoam.com/documentation/
- **用户手册**: https://www.openfoam.org/docs/user/
- **CFD Online**: https://www.cfd-online.com/Forums/openfoam/

### LaserbeamFoam特定

- 求解器源代码: 检查`$FOAM_APP/solvers/`或定制路径
- 案例库: `$FOAM_TUTORIALS/multiphase/`

---

## 附录: 完整输入文件清单

### System目录
- [x] controlDict - 时间和IO控制
- [x] fvSchemes - 离散格式
- [x] fvSolution - 求解器设置
- [x] blockMeshDict - 网格生成
- [x] decomposeParDict - 并行分解
- [x] setFieldsDict - 场初始化
- [x] bedPlateDict - 基板定义

### Constant目录
- [x] LaserProperties - 激光参数
- [x] transportProperties - 材料和输运性质
- [x] turbulenceProperties - 湍流模型(层流)
- [x] g - 重力
- [x] timeVsLaserPosition - 激光路径
- [x] timeVsLaserPower - 激光功率
- [x] location - 粉末颗粒位置(DEM生成)

### 初始场(0或initial目录)
- [x] T - 温度
- [x] U - 速度
- [x] p_rgh - 修正压力
- [x] alpha.metal - 金属相分数
- [x] TRHS - 温度源项
- [x] Laser_boundary - 激光边界标识

### DEM相关
- [x] DEM_large/input.liggghts - LIGGGHTS输入文件
- [x] DEM_large/meshes/*.stl - 边界STL文件

### 脚本
- [x] Allrun - 运行脚本
- [x] Allclean - 清理脚本

---

**文档版本**: 1.0
**生成日期**: 2025-11-04
**作者**: Claude (AI Assistant)
**案例**: LPBF_large - LaserbeamFoam

---

## 总结

本文档详细解读了LPBF_large案例的所有输入文件，涵盖:

1. **物理模型**: 多相流、激光加热、相变、Marangoni效应、反冲压力
2. **数值方法**: VOF界面追踪、PISO压力速度耦合、相变迭代
3. **计算设置**: 网格、时间步长、边界条件、求解器参数
4. **工艺参数**: 激光功率、扫描速度、粉末床结构
5. **运行流程**: 从DEM粉末生成到CFD模拟的完整流程

希望这份文档能帮助您深入理解LPBF激光熔融模拟的方方面面!
