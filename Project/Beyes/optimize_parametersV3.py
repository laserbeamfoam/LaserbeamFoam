import numpy as np
import pandas as pd
import os
from skopt import Optimizer
from skopt.space import Real
from skopt.sampler import Lhs
import matplotlib.pyplot as plt
import shutil
import json

# 添加中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# --- 1. 定义优化问题：参数空间 ---
# 根据您的需求，定义需要优化的参数及其合理的物理范围
# 格式: Real(下界, 上界, name='参数名')
PARAMETER_SPACE = [
    Real(293.15, 973.15, name='substrate_temp'),
    Real(0.3, 1, name='absorptivity'),
    Real(900, 1800, name='tension_A'),
    Real(0.1, 1, name='tension_B'),
    Real(8000, 20000, name='vapor_A'),
    Real(6, 12, name='vapor_B'),
    Real(0.0010, 0.0050, name='laser_spot_radius'),
    Real(0.1, 1, name='laser_focus_radius')
]

# 定义需要对齐的激光功率 (W)
POWER_SETTINGS = [140, 170, 200, 230, 260]

# --- 2. 文件路径定义 ---
EXPERIMENTAL_DATA_FILE = 'experimental_data.csv'
SIMULATION_RESULTS_FILE = 'simulation_results.csv'
OPTIMAL_PARAMS_FILE = 'optimal_parameters.json'

def generate_initial_data_if_not_exist(num_initial_points=15):
    """
    如果数据文件不存在，则创建包含初始拉丁超立方抽样点的模板文件。
    【已更新】现在会包含面积(area)列。
    """
    if not os.path.exists(EXPERIMENTAL_DATA_FILE):
        print(f"警告: 未找到实验数据文件 '{EXPERIMENTAL_DATA_FILE}'。")
        print("正在创建一个模板文件，请填入您的真实实验数据。")
        exp_data = {'power_W': POWER_SETTINGS}
        for p in POWER_SETTINGS:
            # 生成示例数据
            depth = p + 50
            width = p + 150
            exp_data.setdefault('depth_um', []).append(depth)
            exp_data.setdefault('width_um', []).append(width)
            # 使用半椭圆面积作为示例估算值，您需要替换为真实测量值
            exp_data.setdefault('area_um2', []).append(np.pi * depth * width / 2) 
        pd.DataFrame(exp_data).to_csv(EXPERIMENTAL_DATA_FILE, index=False)
        print(f"模板 '{EXPERIMENTAL_DATA_FILE}' 已创建。")

    if not os.path.exists(SIMULATION_RESULTS_FILE):
        print(f"警告: 未找到模拟结果文件 '{SIMULATION_RESULTS_FILE}'。")
        print(f"正在使用拉丁超立方抽样 (LHS) 生成 {num_initial_points} 个初始点用于模拟。")
        
        param_names = [p.name for p in PARAMETER_SPACE]
        lhs_sampler = Lhs(lhs_type="classic", criterion=None)
        initial_samples = lhs_sampler.generate(PARAMETER_SPACE, num_initial_points)
        
        # 【已更新】增加面积(area)的列
        columns = param_names + \
                  [f'depth_{p}W' for p in POWER_SETTINGS] + \
                  [f'width_{p}W' for p in POWER_SETTINGS] + \
                  [f'area_um2_{p}W' for p in POWER_SETTINGS]
        
        sim_df = pd.DataFrame(columns=columns)
        for i, sample in enumerate(initial_samples):
            row = {name: val for name, val in zip(param_names, sample)}
            # 用 "NaN" 填充结果列，提示用户需要运行模拟
            for col in columns[len(param_names):]:
                row[col] = np.nan
            sim_df.loc[i] = row
            
        sim_df.to_csv(SIMULATION_RESULTS_FILE, index=False)
        print("-" * 60)
        print(f"重要提示: 模板文件 '{SIMULATION_RESULTS_FILE}' 已创建。")
        print("请您根据此文件中的前15行参数，运行Flow3D模拟。")
        print("完成后，请将对应的熔池深度(um)、宽度(um)和面积(um^2)结果")
        print("填入该文件，然后重新运行此脚本。")
        print("-" * 60)
        exit() # 第一次生成文件后退出，等待用户填充数据

def calculate_objective_function(sim_row, exp_df):
    """
    计算单个模拟结果行与实验数据之间的加权均方根误差 (WRMSE)。
    【已更新】增加了深宽比(aspect ratio)和面积(area)的误差项。
    
    Args:
        sim_row (pd.Series): 模拟结果文件中的一行数据。
        exp_df (pd.DataFrame): 包含实验数据的DataFrame。
        
    Returns:
        float: 计算出的总误差。
    """
    errors_sq = []
    # --- 权重定义 ---
    # 您可以根据对各项指标的重视程度调整权重
    w_d = 1.0    # 深度权重
    w_w = 1.0    # 宽度权重
    w_ar = 0.5   # 深宽比权重 (aspect ratio)
    w_area = 0.8 # 面积权重

    for power in POWER_SETTINGS:
        # --- 1. 获取模拟和实验数据 ---
        sim_depth = sim_row[f'depth_{power}W']
        sim_width = sim_row[f'width_{power}W']
        sim_area = sim_row[f'area_um2_{power}W']
        
        exp_row = exp_df[exp_df['power_W'] == power].iloc[0]
        exp_depth = exp_row['depth_um']
        exp_width = exp_row['width_um']
        exp_area = exp_row['area_um2']
        
        # 检查数据是否有效 (防止用户未填写)
        if pd.isna(sim_depth) or pd.isna(sim_width) or pd.isna(sim_area):
            return np.inf # 返回无穷大误差，表示此行数据无效

        # --- 2. 计算深宽比 ---
        # 避免除以零的错误
        sim_ar = sim_width / sim_depth if sim_depth > 0 else 0
        exp_ar = exp_width / exp_depth if exp_depth > 0 else 0

        # --- 3. 计算各项相对误差的平方 ---
        err_depth_sq = ((sim_depth - exp_depth) / exp_depth) ** 2
        err_width_sq = ((sim_width - exp_width) / exp_width) ** 2
        err_area_sq = ((sim_area - exp_area) / exp_area) ** 2
        # 对于深宽比，如果实验值为0，则使用绝对误差
        err_ar_sq = ((sim_ar - exp_ar) / exp_ar) ** 2 if exp_ar > 0 else (sim_ar - exp_ar)**2

        # --- 4. 累加加权误差 ---
        total_err_sq_for_power = (w_d * err_depth_sq + 
                                  w_w * err_width_sq + 
                                  w_ar * err_ar_sq + 
                                  w_area * err_area_sq)
        
        errors_sq.append(total_err_sq_for_power)
        
    # 返回所有工况下总误差的均方根
    return np.sqrt(np.mean(errors_sq))


def add_params_to_csv(params_dict_1, params_dict_2, current_sim_count):
    """
    将新的参数组合添加到simulation_results.csv文件中
    
    Args:
        params_dict_1 (dict): 第一组参数
        params_dict_2 (dict): 第二组参数
        current_sim_count (int): 当前模拟数量
    """
    import pandas as pd
    
    # 读取现有的CSV文件
    sim_df = pd.read_csv(SIMULATION_RESULTS_FILE)
    
    # 准备新的行数据
    new_rows = []
    
    # 第一组参数
    row1 = {
        'substrate_temp': params_dict_1['substrate_temp'],
        'absorptivity': params_dict_1['absorptivity'],
        'tension_A': params_dict_1['tension_A'],
        'tension_B': params_dict_1['tension_B'],
        'vapor_A': params_dict_1['vapor_A'],
        'vapor_B': params_dict_1['vapor_B'],
        'laser_spot_radius': params_dict_1['laser_spot_radius'],
        'laser_focus_radius': params_dict_1['laser_focus_radius']
    }
    # 添加空的结果列（等待模拟完成后填写）
    for power in [140, 170, 200, 230, 260]:
        row1[f'depth_{power}W'] = None
        row1[f'width_{power}W'] = None
        row1[f'area_um2_{power}W'] = None
    
    new_rows.append(row1)
    
    # 第二组参数
    row2 = {
        'substrate_temp': params_dict_2['substrate_temp'],
        'absorptivity': params_dict_2['absorptivity'],
        'tension_A': params_dict_2['tension_A'],
        'tension_B': params_dict_2['tension_B'],
        'vapor_A': params_dict_2['vapor_A'],
        'vapor_B': params_dict_2['vapor_B'],
        'laser_spot_radius': params_dict_2['laser_spot_radius'],
        'laser_focus_radius': params_dict_2['laser_focus_radius']
    }
    # 添加空的结果列（等待模拟完成后填写）
    for power in [140, 170, 200, 230, 260]:
        row2[f'depth_{power}W'] = None
        row2[f'width_{power}W'] = None
        row2[f'area_um2_{power}W'] = None
    
    new_rows.append(row2)
    
    # 将新行添加到DataFrame
    new_df = pd.DataFrame(new_rows)
    updated_df = pd.concat([sim_df, new_df], ignore_index=True)
    
    # 保存更新后的CSV文件
    updated_df.to_csv(SIMULATION_RESULTS_FILE, index=False)
    
    print(f"\n📝 已将两组新参数添加到 {SIMULATION_RESULTS_FILE}")
    print(f"   第 {current_sim_count + 1} 行和第 {current_sim_count + 2} 行已添加（参数部分）")
    print("   模拟完成后，请填写对应的深度、宽度、面积结果")


def generate_flow3d_files(next_params_dict_1, next_params_dict_2, current_sim_count):
    """
    根据两组优化参数生成Flow3D模拟文件
    
    Args:
        next_params_dict_1 (dict): 第一组参数字典
        next_params_dict_2 (dict): 第二组参数字典
        current_sim_count (int): 当前模拟数据行数（不包括表头）
    """
    base_folder = "Base"
    if not os.path.exists(base_folder):
        print(f"错误: 未找到基准文件夹 '{base_folder}'")
        return
    
    # 生成两组文件
    param_sets = [next_params_dict_1, next_params_dict_2]
    sim_numbers = [current_sim_count + 1, current_sim_count + 2]
    
    for i, (params, sim_num) in enumerate(zip(param_sets, sim_numbers)):
        print(f"\n正在生成第 {sim_num} 次模拟文件...")
        
        # 为每个功率等级创建文件夹和文件
        for power in POWER_SETTINGS:
            folder_name = f"{sim_num}-{power}W"
            
            # 复制Base文件夹
            if os.path.exists(folder_name):
                shutil.rmtree(folder_name)
            shutil.copytree(base_folder, folder_name)
            
            # 重命名文件
            old_files = {
                f"{sim_num}-{power}W.Flow3dProj": "19-140W.Flow3dProj",
                f"cfgplot_f3d.{sim_num}-{power}W": "cfgplot_f3d.19-140W",
                f"prepin.{sim_num}-{power}W": "prepin.19-140W",
                f"versionOpts.{sim_num}-{power}W": "versionOpts.19-140W"
            }
            
            for new_name, old_name in old_files.items():
                old_path = os.path.join(folder_name, old_name)
                new_path = os.path.join(folder_name, new_name)
                if os.path.exists(old_path):
                    os.rename(old_path, new_path)
            
            # 修改prepin文件
            prepin_file = os.path.join(folder_name, f"prepin.{sim_num}-{power}W")
            modify_prepin_file(prepin_file, params, power, f"{sim_num}-{power}W")
            
            print(f"  已生成: {folder_name}/")
    
    print(f"\n✅ 成功生成了 {len(POWER_SETTINGS) * 2} 个Flow3D模拟文件夹")
    print("每组参数对应5个功率等级 (140W, 170W, 200W, 230W, 260W)")



def modify_prepin_file(file_path, params, power, project_name):
    """
    修改prepin文件中的参数
    
    Args:
        file_path (str): prepin文件路径
        params (dict): 参数字典
        power (int): 功率值
        project_name (str): 项目名称
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 1. 修改project名称
    content = content.replace("project='19-140W'", f"project='{project_name}'")
    
    # 2. 修改substrate_temp (treg(1)和treg(2))
    substrate_temp = params['substrate_temp']
    content = content.replace(
        "treg(1)=701.417894, remark='Fluid Tempearture',",
        f"treg(1)={substrate_temp:.6f}, remark='Fluid Tempearture',"
    )
    content = content.replace(
        "treg(2)=701.417894, remark='Fluid Tempearture',",
        f"treg(2)={substrate_temp:.6f}, remark='Fluid Tempearture',"
    )
    
    # 3. 修改absorptivity (cbmengy)
    absorptivity = params['absorptivity']
    content = content.replace(
        "cbmengy=0.3951423, remark='temperature dependent fluid-absorption rate',",
        f"cbmengy={absorptivity:.6f}, remark='temperature dependent fluid-absorption rate',"
    )
    
    # 4. 修改vapor_A和vapor_B
    vapor_A = params['vapor_A']
    vapor_B = params['vapor_B']
    content = content.replace(
        "avprslbm=15672.305017, remark='evaporation pressure pscoefficient A',",
        f"avprslbm={vapor_A:.6f}, remark='evaporation pressure pscoefficient A',"
    )
    content = content.replace(
        "bvprslbm=11.622112, remark='evaporation pressure pscoefficient B',",
        f"bvprslbm={vapor_B:.6f}, remark='evaporation pressure pscoefficient B',"
    )
    
    # 5. 修改laser_spot_radius (r0lbm(1)和rflbm(1))
    laser_spot_radius = params['laser_spot_radius']
    content = content.replace(
        "r0lbm(1)=0.00181, remark='lens radius',",
        f"r0lbm(1)={laser_spot_radius:.6f}, remark='lens radius',"
    )
    content = content.replace(
        "rflbm(1)=0.00181, remark='spot radius',",
        f"rflbm(1)={laser_spot_radius:.6f}, remark='spot radius',"
    )
    
    # 6. 修改rblbm(1) = laser_spot_radius * 0.707107 * laser_focus_radius
    laser_focus_radius = params['laser_focus_radius']
    rblbm_value = laser_spot_radius * 0.707107 * laser_focus_radius
    content = content.replace(
        "rblbm(1)=0.001200382, remark='Gauss dist. (rb)',",
        f"rblbm(1)={rblbm_value:.6f}, remark='Gauss dist. (rb)',"
    )
    
    # 7. 修改功率 (powlbm)
    power_value = power * 1e7  # 转换为科学计数法格式
    content = content.replace(
        "powlbm(1,1)=2.6e9,",
        f"powlbm(1,1)={power_value:.1e},"
    )
    content = content.replace(
        "powlbm(2,1)=2.6e9,",
        f"powlbm(2,1)={power_value:.1e},"
    )
    
    # 8. 修改sigma表格中的tension_A和tension_B
    tension_A = params['tension_A']
    tension_B = params['tension_B']
    
    # 计算第二行和第三行的值: max(0, tension_A - (3090-1708) * tension_B)
    sigma_2_3_value = max(0, tension_A - (3090 - 1708) * tension_B)
    
    # 替换sigma表格
    sigma_old = "#sigma\n1708\t1278.72\n3090\t1278.72\n5000\t1278.72\n#end sigma"
    sigma_new = f"#sigma\n1708\t{tension_A:.6f}\n3090\t{sigma_2_3_value:.6f}\n5000\t{sigma_2_3_value:.6f}\n#end sigma"
    content = content.replace(sigma_old, sigma_new)
    
    # 保存修改后的文件
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(content)


def write_optimal_parameters(params_dict, error_value):
    """
    将最优参数写入JSON文件
    
    Args:
        params_dict (dict): 最优参数字典
        error_value (float): 对应的误差值
    """
    output_data = {
        'optimal_parameters': params_dict,
        'wrmse_error': error_value,
        'timestamp': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'description': '贝叶斯优化得到的最优参数组合'
    }
    
    with open(OPTIMAL_PARAMS_FILE, 'w', encoding='utf-8') as f:
        json.dump(output_data, f, indent=4, ensure_ascii=False)
    
    print(f"\n📝 最优参数已写入文件: {OPTIMAL_PARAMS_FILE}")

def add_single_param_to_csv(params_dict, current_sim_count):
    """
    将单个参数组合添加到simulation_results.csv文件中
    
    Args:
        params_dict (dict): 参数字典
        current_sim_count (int): 当前模拟数量
    """
    import pandas as pd
    
    # 读取现有的CSV文件
    sim_df = pd.read_csv(SIMULATION_RESULTS_FILE)
    
    # 准备新的行数据
    new_row = {
        'substrate_temp': params_dict['substrate_temp'],
        'absorptivity': params_dict['absorptivity'],
        'tension_A': params_dict['tension_A'],
        'tension_B': params_dict['tension_B'],
        'vapor_A': params_dict['vapor_A'],
        'vapor_B': params_dict['vapor_B'],
        'laser_spot_radius': params_dict['laser_spot_radius'],
        'laser_focus_radius': params_dict['laser_focus_radius']
    }
    # 添加空的结果列（等待模拟完成后填写）
    for power in [140, 170, 200, 230, 260]:
        new_row[f'depth_{power}W'] = None
        new_row[f'width_{power}W'] = None
        new_row[f'area_um2_{power}W'] = None
    
    # 将新行添加到DataFrame
    new_df = pd.DataFrame([new_row])
    updated_df = pd.concat([sim_df, new_df], ignore_index=True)
    
    # 保存更新后的CSV文件
    updated_df.to_csv(SIMULATION_RESULTS_FILE, index=False)
    
    print(f"\n📝 已将新参数添加到 {SIMULATION_RESULTS_FILE}")
    print(f"   第 {current_sim_count + 1} 行已添加（参数部分）")
    print("   模拟完成后，请填写对应的深度、宽度、面积结果")


def generate_single_flow3d_files(next_params_dict, current_sim_count):
    """
    根据单组优化参数生成Flow3D模拟文件
    
    Args:
        next_params_dict (dict): 参数字典
        current_sim_count (int): 当前模拟数据行数（不包括表头）
    """
    base_folder = "Base"
    if not os.path.exists(base_folder):
        print(f"错误: 未找到基准文件夹 '{base_folder}'")
        return
    
    sim_num = current_sim_count + 1
    print(f"\n正在生成第 {sim_num} 次模拟文件...")
    
    # 为每个功率等级创建文件夹和文件
    for power in POWER_SETTINGS:
        folder_name = f"{sim_num}-{power}W"
        
        # 复制Base文件夹
        if os.path.exists(folder_name):
            shutil.rmtree(folder_name)
        shutil.copytree(base_folder, folder_name)
        
        # 重命名文件
        old_files = {
            f"{sim_num}-{power}W.Flow3dProj": "19-140W.Flow3dProj",
            f"cfgplot_f3d.{sim_num}-{power}W": "cfgplot_f3d.19-140W",
            f"prepin.{sim_num}-{power}W": "prepin.19-140W",
            f"versionOpts.{sim_num}-{power}W": "versionOpts.19-140W"
        }
        
        for new_name, old_name in old_files.items():
            old_path = os.path.join(folder_name, old_name)
            new_path = os.path.join(folder_name, new_name)
            if os.path.exists(old_path):
                os.rename(old_path, new_path)
        
        # 修改prepin文件
        prepin_file = os.path.join(folder_name, f"prepin.{sim_num}-{power}W")
        modify_prepin_file(prepin_file, next_params_dict, power, f"{sim_num}-{power}W")
        
        print(f"  已生成: {folder_name}/")
    
    print(f"\n✅ 成功生成了 {len(POWER_SETTINGS)} 个Flow3D模拟文件夹")
    print("对应5个功率等级 (140W, 170W, 200W, 230W, 260W)")


def run_bayesian_optimization():
    """
    主优化流程。
    """
    # 1. 加载数据
    try:
        exp_df = pd.read_csv(EXPERIMENTAL_DATA_FILE)
        sim_df = pd.read_csv(SIMULATION_RESULTS_FILE)
    except FileNotFoundError as e:
        print(f"错误: 缺少数据文件 {e.filename}。请先运行脚本生成模板。")
        exit()

    # 2. 准备已观测的数据点 (X_observed, y_observed)
    param_names = [p.name for p in PARAMETER_SPACE]
    
    # 筛选出已经完成模拟并有结果的行
    completed_sim_df = sim_df.dropna()
    
    if len(completed_sim_df) == 0:
        print("错误: 模拟结果文件不包含任何已完成的模拟数据。")
        print(f"请先运行 '{SIMULATION_RESULTS_FILE}' 中建议的模拟并填写结果。")
        exit()

    X_observed = completed_sim_df[param_names].values.tolist()
    y_observed = [calculate_objective_function(row, exp_df) for _, row in completed_sim_df.iterrows()]
    
    print(f"成功加载 {len(X_observed)} 组已完成的模拟数据。")

    # 找到当前最优参数
    best_index = np.argmin(y_observed)
    best_params_so_far = X_observed[best_index]
    best_error_so_far = y_observed[best_index]
    
    # 将最优参数转换为字典格式
    best_params_dict = {p.name: val for p, val in zip(PARAMETER_SPACE, best_params_so_far)}
    
    # 写入最优参数到文件
    write_optimal_parameters(best_params_dict, best_error_so_far)

    # 检查是否有未完成的模拟
    if len(sim_df) > len(completed_sim_df):
        print("检测到有未完成的模拟参数组。请先填写所有模拟结果，然后重新运行程序。")
        print("将直接显示当前结果的可视化。")
    else:
        # 3. 初始化贝叶斯优化器
        optimizer = Optimizer(
            dimensions=PARAMETER_SPACE,
            base_estimator="GP",  # 使用高斯过程作为代理模型
            acq_func="EI",       # 使用预期提升(Expected Improvement)作为采集函数
            random_state=42
        )

        # 4. "告知"优化器我们已经观测到的所有点
        optimizer.tell(X_observed, y_observed)
        
        print("优化器已根据现有数据完成代理模型训练。")
        
        # 5. "询问"优化器下一个最有价值去模拟的点
        next_point = optimizer.ask()
        
        # 6. 显示结果
        print("\n" + "="*60)
        print("           >>> 贝叶斯优化结果 <<<")
        print("="*60)

        print("\n--- 目前最优参数组合 ---")
        for param, value in best_params_dict.items():
            print(f"  {param:<20}: {value:.6f}")
        print(f"  > 对应的最低误差 (WRMSE): {best_error_so_far:.4f}")

        print("\n--- ★★★ 下一步行动建议 ★★★ ---")
        print("建议您使用以下参数运行Flow3D模拟：")
        
        next_params_dict = {p.name: val for p, val in zip(PARAMETER_SPACE, next_point)}
        for name, value in next_params_dict.items():
            print(f"  {name:<20}: {value:.6f}")
        
        print("\n模拟完成后，请在CSV文件中填写对应的熔池深度、宽度、面积结果，")
        print(f"参数已自动添加到 '{SIMULATION_RESULTS_FILE}' 文件中，然后再次运行本程序。")
        print("="*60)
        
        # 6.5. 自动生成Flow3D模拟文件
        try:
            current_sim_count = len(sim_df)  # 当前CSV文件中的数据行数（不包括表头）
            print(f"\n📊 检测到当前有 {current_sim_count} 行模拟数据")
            print(f"将生成第 {current_sim_count + 1} 次模拟文件")
            generate_single_flow3d_files(next_params_dict, current_sim_count)
            add_single_param_to_csv(next_params_dict, current_sim_count)
        except Exception as e:
            print(f"\n⚠️  生成Flow3D文件时出现错误: {e}")
            print("请手动检查Base文件夹是否存在，或联系开发者。")

    # 7. 按参数组顺序可视化WRMSE结果（保留的唯一可视化功能）
    if len(y_observed) > 1:
        plt.figure(figsize=(12, 6))
        x_indices = list(range(1, len(y_observed) + 1))  # 横坐标为参数组序号
        plt.bar(x_indices, y_observed, color='skyblue')
        plt.xlabel('参数组序号')
        plt.ylabel('WRMSE')
        plt.title('各参数组的WRMSE值')
        plt.xticks(x_indices)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        # 在每个柱子上标注WRMSE值
        for i, v in enumerate(y_observed):
            plt.text(i + 1, v + 0.02, f'{v:.4f}', ha='center', va='bottom', fontsize=9)
            
        plt.tight_layout()
        plt.show()
    
    print(f"\n📄 最优参数详情已保存到: {OPTIMAL_PARAMS_FILE}")

if __name__ == '__main__':
    # 检查并生成模板文件（如果需要）
    generate_initial_data_if_not_exist()
    
    # 运行主程序
    run_bayesian_optimization()