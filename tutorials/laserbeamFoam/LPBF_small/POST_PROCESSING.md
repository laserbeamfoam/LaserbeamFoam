# LPBF案例后处理指南 (Post-processing Guide)

## 在ParaView中打开案例

### 方法1: 使用.foam文件（推荐）
```bash
# 在案例目录中，打开现有的.foam文件
paraview PowderBed.foam
```

### 方法2: 创建.foam文件
```bash
# 如果没有.foam文件，创建一个空文件
touch case.foam
paraview case.foam
```

## 重要的流体状态场变量

根据你的模拟结果，以下是关键的场变量及其含义：

### 1. **相态识别变量**

#### `alpha.metal` (金属相体积分数)
- **范围**: 0-1
- **含义**:
  - `1` = 纯金属相
  - `0` = 纯气相
  - `0-1` = 界面区域
- **用途**: 识别金属-气体界面，粉末床形态

#### `condition` (相态条件)
- **含义**: 指示材料的状态
- **典型值**:
  - 固态区域
  - 熔融区域
  - 气化区域

### 2. **温度相关变量**

#### `T` (温度场)
- **单位**: Kelvin (K)
- **关键温度点**:
  - `Tsolidus = 1658 K` (固相线温度)
  - `Tliquidus = 1723 K` (液相线温度)
  - `Tvap = 3068 K` (气化温度)

#### `TLiquidus` 和 `TSolidus`
- 这些是参考温度场，用于判断相变

#### `meltHistory` (熔化历史)
- 记录材料是否经历过熔化

### 3. **热源相关变量**

#### `Qv` (体积热源)
- **单位**: W/m³
- **含义**: 激光能量沉积率
- **用途**: 显示激光加热区域

#### `Laser_boundary`
- 激光边界条件

#### `rayQ`
- 光线追踪的热量

### 4. **流体动力学变量**

#### `U` (速度场)
- **单位**: m/s
- **含义**: 流体速度矢量
- **用途**:
  - 显示熔池内的对流
  - Marangoni流动
  - 反冲压力引起的流动

#### `p_rgh` 和 `p` (压力)
- `p_rgh`: 修正压力 (用于浮力计算)
- `p`: 总压力

#### `phi` (通量)
- 面通量场

### 5. **材料属性变量**

#### `kappa` (热导率)
- **单位**: W/(m·K)
- 温度相关的热导率

#### `cp` (比热容)
- **单位**: J/(kg·K)
- 温度相关的比热容

## ParaView后处理步骤

### 1. 基本可视化设置

#### 显示熔池形态：
1. 加载 `alpha.metal` 场
2. 使用 **Contour** 滤波器：
   - Contour By: `alpha.metal`
   - Value: `0.5` (界面位置)
3. 或使用 **Threshold** 滤波器：
   - Scalar: `alpha.metal`
   - Range: `0.5-1.0` (只显示金属相)

#### 显示温度分布：
1. 选择 `T` 场作为颜色映射
2. 设置颜色范围：
   - Minimum: `300 K` (环境温度)
   - Maximum: `3500 K` (超过气化温度)
3. 推荐色标: **Rainbow** 或 **Black-Body Radiation**

#### 识别不同相态：
1. 使用 **Calculator** 滤波器创建相态指示器：
   ```
   # 固态 (T < Tsolidus)
   Result Array Name: "Solid"
   Formula: T < 1658

   # 液态 (Tsolidus < T < Tliquidus)
   Result Array Name: "Liquid"
   Formula: (T >= 1658) * (T <= 1723)

   # 过热液态 (T > Tliquidus)
   Result Array Name: "Superheated"
   Formula: (T > 1723) * (T < 3068)

   # 气化 (T > Tvap)
   Result Array Name: "Vapor"
   Formula: T >= 3068
   ```

### 2. 显示激光光线

激光光线已经导出到VTK文件：
1. 打开: `File -> Open -> VTKs/rays_laser0.vtk.series`
2. 这会自动加载所有时间步的光线，时间同步正确
3. 显示设置：
   - 颜色: 黑色或红色
   - **增加线宽**: Properties -> Line Width = 3
   - 或使用 **Tube** 滤波器使光线更明显：
     - Radius: `1e-5` 到 `1e-4` (根据网格大小调整)

### 3. 创建切片观察内部

1. 使用 **Slice** 滤波器：
   - Origin: 选择合适的位置
   - Normal:
     - `[1, 0, 0]` - YZ平面
     - `[0, 1, 0]` - XZ平面
     - `[0, 0, 1]` - XY平面

2. 在切片上显示：
   - 温度分布 (`T`)
   - 速度场 (`U` with **Glyph** 或 streamlines)
   - 熔池形态 (`alpha.metal` contour)

### 4. 速度场可视化

#### 流线 (Streamlines):
1. 选择数据集
2. Filters -> **Stream Tracer**
3. Vectors: `U`
4. Seed Type: Point Source 或 Line Source
5. 调整 Maximum Streamline Length

#### 矢量箭头 (Glyphs):
1. Filters -> **Glyph**
2. Glyph Type: Arrow
3. Vectors: `U`
4. Scale Mode: vector
5. Scale Factor: 调整箭头大小

### 5. 等值面 (Isosurfaces)

#### 熔池边界：
- **Contour**: `T = 1723 K` (液相线)
- **Contour**: `T = 1658 K` (固相线)

#### 金属-气体界面：
- **Contour**: `alpha.metal = 0.5`

#### 高温区域：
- **Threshold**: `T > 3068 K` (气化区)

## 常用分析流程

### 分析熔池尺寸和形态：
1. 提取 `T = Tliquidus` 等值面
2. 使用 **Integrate Variables** 滤波器计算熔池体积
3. 使用 **Plot Over Line** 提取熔池深度和宽度

### 分析温度历史：
1. 使用 **Plot Selection Over Time**
2. 选择关键点位置
3. 导出温度-时间曲线

### 分析流场：
1. 在熔池中心创建切片
2. 显示速度矢量或流线
3. 观察Marangoni对流模式

## 常见问题

### Q: 看不到金属相？
**A**: 检查 `alpha.metal` 的范围。使用 Threshold 滤波器，范围设为 `0.3-1.0`。

### Q: 只显示粉末，基板消失了？⚠️ **重要**
**A**: 这是因为基板和粉末的 `alpha.metal` **都是1**，无法用 `alpha.metal` 区分！

**解决方案1 - 使用温度场区分（推荐）：**
1. 不使用 Threshold，直接显示整个域
2. 用**温度场 `T`** 作为颜色：
   - 基板：初始温度 ~300K（冷）
   - 粉末：可能稍高或相同
   - 熔池：>1723K（热）
3. 这样可以看到完整的几何形态

**解决方案2 - 使用 `Deposition` 或 `meltHistory` 场：**
1. `Deposition` 场可能标记了粉末区域
2. `meltHistory` 显示经历过熔化的区域
3. 使用 Threshold: `Deposition > 0.5` 只显示粉末
4. 或 `meltHistory < 0.5` 显示未熔化区域

**解决方案3 - 使用切片观察：**
1. 创建 **Slice** (切片) 滤波器
2. 选择 Z方向切片，在不同高度观察：
   - z < 100μm: 基板区域
   - z ≈ 100-120μm: 粉末床区域
3. 在切片上显示温度和速度场

**解决方案4 - 结合几何信息：**
根据 `blockMeshDict`，几何为：
- X: 0-100μm, Y: 0-150μm, Z: 0-200μm
- 粉末床大约在 Z ≈ 100-120μm
- 使用 **Clip** 滤波器裁剪显示区域：
  - Clip Type: Plane
  - Origin: [0, 0, 100e-6]
  - Normal: [0, 0, -1]
  - 这样只显示 z > 100μm 的部分（粉末层）

### Q: 温度场显示异常？
**A**:
- 检查色标范围是否合适
- 确认单位 (应该是Kelvin)
- 尝试对数色标如果温度范围很大

### Q: 光线不显示？
**A**:
- 确认打开了 `rays_laser0.vtk.series` 文件
- 增加线宽或使用Tube滤波器
- 检查光线是否在当前时间步存在 (从7e-05后光线消失)

### Q: 如何区分固态和液态？
**A**: 使用Calculator滤波器和温度阈值：
- 固态: `T < 1658`
- 液态: `1658 < T < 1723`
- 过热液态: `T > 1723`

## 导出数据

### 导出图像：
- `File -> Save Screenshot`
- 设置分辨率 (推荐 1920x1080 或更高)

### 导出动画：
- `File -> Save Animation`
- 选择格式 (AVI, MP4等)
- 设置帧率和分辨率

### 导出数据：
- `File -> Save Data`
- 支持CSV, VTK等格式

## 推荐可视化组合

### 方案1: 完整显示基板+粉末+熔池（最推荐）✨
1. **不使用任何Threshold/Contour滤波器**，直接显示原始数据
2. **颜色映射**选择 `T` (温度场)
3. **色标设置**:
   - 范围: 300K - 3500K
   - 色标: Rainbow 或 Black-Body Radiation
4. **添加不透明度**（可选）:
   - Edit -> Color Map Editor
   - 在低温区域(300-1000K)设置半透明，可以看到内部熔池
5. **添加激光光线**:
   - 打开 `VTKs/rays_laser0.vtk.series`
   - 红色，线宽3-5
6. **效果**: 基板（冷，蓝色）+ 粉末（温热）+ 熔池（热，红色）全部可见

### 方案2: 使用 Deposition 场区分基板和粉末
1. 显示 `Deposition` 场（标记粉末区域）
2. **Threshold** 滤波器:
   - Scalar: `Deposition`
   - Range: `0.5 - 1.0` (只显示粉末)
3. 颜色用 `T` 温度场
4. 或创建两个显示对象：
   - 对象1: Threshold `Deposition > 0.5` (粉末) - 显示温度
   - 对象2: Threshold `Deposition < 0.5` (基板) - 灰色

### 方案3: 切片组合视图
1. **Slice 1** (纵向切片，XZ面):
   - Origin: [50e-6, 75e-6, 0]
   - Normal: [0, 1, 0]
   - 颜色: 温度 `T`
2. **Slice 2** (横向切片，XY面):
   - Origin: [0, 0, 110e-6] (粉末床高度)
   - Normal: [0, 0, 1]
   - 颜色: 温度 `T`
3. 叠加速度矢量/流线
4. 添加激光光线

### 方案4: 标准熔池可视化（仅关注熔池）
1. **Contour**: `T = 1723 K` (液相线) - 显示熔池边界，红色
2. **Contour**: `alpha.metal = 0.5` - 显示金属-气体界面，灰色
3. **Clip** 裁剪显示区域（只看粉末层）:
   - Plane origin: [0, 0, 100e-6]
   - Normal: [0, 0, -1]
4. 温度场 (volume rendering 或 slice) - Rainbow色标
5. 激光光线 (从VTK加载) - 红色，粗线
6. 速度流线 - 白色或黑色

### 方案5: 相态分析
1. 使用Calculator创建相态场
2. 用不同颜色表示固态、液态、气态
3. 叠加温度等值线

### 方案6: 动态演化
1. 播放动画显示熔池演化
2. 同时显示光线追踪
3. 叠加温度分布

## 快速操作步骤：同时显示基板和粉末

**最简单的方法** (30秒完成):

1. 打开 ParaView，加载 `PowderBed.foam`
2. 点击 "Apply"
3. 在左侧 Properties 面板:
   - **Coloring**: 选择 `T` (温度)
   - **Rescale to Data Range** (点击温度色标旁的按钮)
4. 调整色标:
   - 点击色标编辑器
   - 设置范围: Min=300, Max=3500
5. **完成！** 现在能看到：
   - 基板（底部，蓝色/绿色，低温）
   - 粉末床（中间，颜色渐变）
   - 熔池（顶部，红色/黄色，高温）
6. 可选：添加激光光线
   - File → Open → `VTKs/rays_laser0.vtk.series`
   - 改成红色，线宽=5
