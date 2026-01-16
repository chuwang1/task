#!/usr/bin/env python3
"""
读取 EQGS1D 导出的 CSV 文件
从 eq 程序导出的 23 个等离子体平衡参数剖面

使用方法:
    python read_eqgs1d_csv.py

或在 Python 中:
    from read_eqgs1d_csv import load_all_profiles, load_profile
    data = load_all_profiles()
    pps = load_profile(1)  # 读取第1个图的数据 (PPS)
"""

import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt

# 23 个图形的描述信息
PLOT_INFO = {
    1:  {'file': 'eqgs1d_01_PPS.csv',      'name': 'PPS',     'desc': '压强剖面',                  'unit': 'MPa'},
    2:  {'file': 'eqgs1d_02_TTS.csv',      'name': 'TTS',     'desc': '环向磁通剖面',              'unit': 'Wb'},
    3:  {'file': 'eqgs1d_03_RLEN.csv',     'name': 'RLEN',    'desc': '磁力线长度',                'unit': 'm'},
    4:  {'file': 'eqgs1d_04_QPS.csv',      'name': 'QPS',     'desc': '安全因子 q',                'unit': '-'},
    5:  {'file': 'eqgs1d_05_VPS.csv',      'name': 'VPS',     'desc': '体积剖面',                  'unit': 'm³'},
    6:  {'file': 'eqgs1d_06_SPS.csv',      'name': 'SPS',     'desc': '磁面面积剖面',              'unit': 'm²'},
    7:  {'file': 'eqgs1d_07_MINMAX.csv',   'name': 'MINMAX',  'desc': 'R,Z,B 最小/最大值',         'unit': 'm,m,T'},
    8:  {'file': 'eqgs1d_08_BB.csv',       'name': 'BB',      'desc': '磁场强度最小/最大值',        'unit': 'T'},
    9:  {'file': 'eqgs1d_09_AVERR2.csv',   'name': 'AVERR2',  'desc': '⟨R²⟩ 磁面平均',            'unit': 'm²'},
    10: {'file': 'eqgs1d_10_AVEIR2.csv',   'name': 'AVEIR2',  'desc': '⟨1/R²⟩ 磁面平均',          'unit': '1/m²'},
    11: {'file': 'eqgs1d_11_AVEBB2.csv',   'name': 'AVEBB2',  'desc': '⟨B²⟩ 磁面平均',            'unit': 'T²'},
    12: {'file': 'eqgs1d_12_AVEIB2.csv',   'name': 'AVEIB2',  'desc': '⟨1/B²⟩ 磁面平均',          'unit': '1/T²'},
    13: {'file': 'eqgs1d_13_AVEGV2.csv',   'name': 'AVEGV2',  'desc': '⟨|∇V|²⟩ 磁面平均',         'unit': 'm⁴'},
    14: {'file': 'eqgs1d_14_AVEGVR2.csv',  'name': 'AVEGVR2', 'desc': '⟨|∇V|²/R²⟩ 磁面平均',      'unit': 'm²'},
    15: {'file': 'eqgs1d_15_AVEGP2.csv',   'name': 'AVEGP2',  'desc': '⟨|∇P|²⟩ 磁面平均',         'unit': 'Pa²/m²'},
    16: {'file': 'eqgs1d_16_AVEJ.csv',     'name': 'AVEJ',    'desc': '平均极向/环向电流密度',      'unit': 'MA/m²'},
    17: {'file': 'eqgs1d_17_RRPSI.csv',    'name': 'RRPSI',   'desc': '主半径剖面',                'unit': 'm'},
    18: {'file': 'eqgs1d_18_RSPSI.csv',    'name': 'RSPSI',   'desc': '次半径剖面',                'unit': 'm'},
    19: {'file': 'eqgs1d_19_ELIPPSI.csv',  'name': 'ELIPPSI', 'desc': '椭圆度剖面',                'unit': '-'},
    20: {'file': 'eqgs1d_20_TRIGPSI.csv',  'name': 'TRIGPSI', 'desc': '三角形度剖面',              'unit': '-'},
    21: {'file': 'eqgs1d_21_DVDPSIP.csv',  'name': 'DVDPSIP', 'desc': 'dV/dψp',                   'unit': 'm³/Wb'},
    22: {'file': 'eqgs1d_22_DVDPSIT.csv',  'name': 'DVDPSIT', 'desc': 'dV/dψt',                   'unit': 'm³/Wb'},
    23: {'file': 'eqgs1d_23_AVEGV.csv',    'name': 'AVEGV',   'desc': '⟨|∇V|⟩ 磁面平均',          'unit': 'm²'},
}


def load_profile(plot_num, directory='.'):
    """
    读取单个剖面数据

    参数:
        plot_num: 图号 (1-23)
        directory: CSV 文件所在目录

    返回:
        pandas.DataFrame: 包含数据的 DataFrame
    """
    if plot_num not in PLOT_INFO:
        raise ValueError(f"Invalid plot number: {plot_num}. Must be 1-23.")

    info = PLOT_INFO[plot_num]
    filepath = os.path.join(directory, info['file'])

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    # 读取 CSV 文件
    df = pd.read_csv(filepath)

    # 添加元数据
    df.attrs['plot_num'] = plot_num
    df.attrs['name'] = info['name']
    df.attrs['description'] = info['desc']
    df.attrs['unit'] = info['unit']

    return df


def load_all_profiles(directory='.'):
    """
    读取所有 23 个剖面数据

    参数:
        directory: CSV 文件所在目录

    返回:
        dict: 键为图号 (1-23)，值为 DataFrame
    """
    data = {}

    for plot_num in range(1, 24):
        try:
            df = load_profile(plot_num, directory)
            data[plot_num] = df
            print(f"✓ [{plot_num:2d}] {PLOT_INFO[plot_num]['name']:10s} - {PLOT_INFO[plot_num]['desc']}")
        except FileNotFoundError:
            print(f"✗ [{plot_num:2d}] {PLOT_INFO[plot_num]['name']:10s} - 文件不存在")

    return data


def plot_profile(df_or_num, directory='.', ax=None, **kwargs):
    """
    绘制单个剖面

    参数:
        df_or_num: DataFrame 或图号
        directory: CSV 文件所在目录（如果 df_or_num 是图号）
        ax: matplotlib axes 对象（可选）
        **kwargs: 传递给 plot 的额外参数
    """
    # 如果输入是图号，先加载数据
    if isinstance(df_or_num, int):
        df = load_profile(df_or_num, directory)
    else:
        df = df_or_num

    # 创建 axes（如果没有提供）
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    # 获取列名（跳过第一列 PSIP）
    x_col = df.columns[0]  # PSIP
    y_cols = df.columns[1:]

    # 绘制所有 Y 列
    for y_col in y_cols:
        ax.plot(df[x_col], df[y_col], label=y_col, **kwargs)

    # 设置标签
    ax.set_xlabel('归一化极向磁通 (PSIP)')
    ax.set_ylabel(f"{df.attrs.get('description', '')} ({df.attrs.get('unit', '')})")
    ax.set_title(f"{df.attrs.get('name', '')} - {df.attrs.get('description', '')}")
    ax.grid(True, alpha=0.3)

    if len(y_cols) > 1:
        ax.legend()

    return ax


def plot_all_profiles(directory='.', save_fig=False):
    """
    绘制所有 23 个剖面（6x4 网格）

    参数:
        directory: CSV 文件所在目录
        save_fig: 是否保存图像
    """
    data = load_all_profiles(directory)

    # 创建 6x4 子图
    fig, axes = plt.subplots(6, 4, figsize=(20, 24))
    axes = axes.flatten()

    for i, plot_num in enumerate(range(1, 24)):
        if plot_num in data:
            plot_profile(data[plot_num], ax=axes[i])
        else:
            axes[i].text(0.5, 0.5, f'数据缺失\n图 {plot_num}',
                        ha='center', va='center', transform=axes[i].transAxes)
            axes[i].set_xticks([])
            axes[i].set_yticks([])

    # 隐藏最后一个空白子图
    axes[23].axis('off')

    plt.tight_layout()

    if save_fig:
        plt.savefig('eqgs1d_all_profiles.png', dpi=150, bbox_inches='tight')
        print(f"\n图像已保存: eqgs1d_all_profiles.png")

    return fig, axes


def export_to_excel(directory='.', output_file='eqgs1d_all_data.xlsx'):
    """
    将所有数据导出到 Excel 文件（每个剖面一个 sheet）

    参数:
        directory: CSV 文件所在目录
        output_file: 输出的 Excel 文件名
    """
    data = load_all_profiles(directory)

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for plot_num, df in data.items():
            sheet_name = f"{plot_num:02d}_{PLOT_INFO[plot_num]['name']}"
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"\n数据已导出到: {output_file}")


def summary_statistics(directory='.'):
    """
    打印所有剖面的统计摘要

    参数:
        directory: CSV 文件所在目录
    """
    data = load_all_profiles(directory)

    print("\n" + "="*80)
    print("EQGS1D 剖面数据统计摘要")
    print("="*80)

    for plot_num, df in data.items():
        info = PLOT_INFO[plot_num]
        print(f"\n[{plot_num:2d}] {info['name']:10s} - {info['desc']}")
        print(f"     单位: {info['unit']}")

        # 统计信息（跳过 PSIP 列）
        for col in df.columns[1:]:
            data_col = df[col].dropna()
            if len(data_col) > 0:
                print(f"     {col:15s}: min={data_col.min():12.4e}, "
                      f"max={data_col.max():12.4e}, "
                      f"mean={data_col.mean():12.4e}")


# 示例用法
if __name__ == '__main__':
    print("="*80)
    print("EQGS1D CSV 数据读取工具")
    print("="*80)

    # 检查 CSV 文件是否存在
    csv_files = glob.glob('eqgs1d_*.csv')
    if not csv_files:
        print("\n错误: 未找到 CSV 文件!")
        print("请先运行: ./eq <in/eq.ITER01_export.in")
        print("或在 eq 程序中执行命令 'se' (export)")
        exit(1)

    print(f"\n找到 {len(csv_files)} 个 CSV 文件\n")

    # 读取所有数据
    print("正在读取数据...")
    data = load_all_profiles()

    # 显示统计信息
    summary_statistics()

    # 示例: 读取并绘制压强剖面
    print("\n" + "="*80)
    print("示例 1: 绘制压强剖面 (PPS)")
    print("="*80)
    pps = load_profile(1)
    print(f"\n数据维度: {pps.shape}")
    print(f"列名: {list(pps.columns)}")
    print(f"\n前 5 行:")
    print(pps.head())

    # 绘制单个图
    plt.figure(figsize=(10, 6))
    plot_profile(1)
    plt.savefig('example_PPS.png', dpi=150, bbox_inches='tight')
    print("\n示例图已保存: example_PPS.png")

    # 导出到 Excel
    print("\n" + "="*80)
    print("导出所有数据到 Excel...")
    print("="*80)
    export_to_excel()

    # 绘制所有图
    print("\n" + "="*80)
    print("绘制所有 23 个剖面...")
    print("="*80)
    plot_all_profiles(save_fig=True)

    print("\n" + "="*80)
    print("完成!")
    print("="*80)
    print("\n可用函数:")
    print("  - load_profile(plot_num)        : 读取单个剖面")
    print("  - load_all_profiles()            : 读取所有剖面")
    print("  - plot_profile(plot_num)         : 绘制单个剖面")
    print("  - plot_all_profiles()            : 绘制所有剖面")
    print("  - export_to_excel()              : 导出到 Excel")
    print("  - summary_statistics()           : 显示统计摘要")
