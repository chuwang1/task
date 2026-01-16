#!/usr/bin/env python3
"""
读取 EQGS2D 导出的 CSV 文件
从 eq 程序导出的 2D 等离子体平衡数据

使用方法:
    python read_eqgs2d_csv.py

或在 Python 中:
    from read_eqgs2d_csv import load_psirz_grid, plot_psirz_contour
    psirz = load_psirz_grid()
    plot_psirz_contour(psirz)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os

def load_psirz_grid(filename='eqgs2d_01_PSIRZ_grid.csv'):
    """
    读取 PSIRZ(R,Z) 网格数据

    返回:
        tuple: (R_grid, Z_grid, PSIRZ_data)
    """
    df = pd.read_csv(filename, index_col=0)

    R_grid = df.columns.astype(float).values
    Z_grid = df.index.astype(float).values
    PSIRZ = df.values

    return R_grid, Z_grid, PSIRZ


def load_derivative_grid(filename):
    """读取 DPSIDR 或 DPSIDZ 网格数据"""
    df = pd.read_csv(filename, index_col=0)

    R_grid = df.columns.astype(float).values
    Z_grid = df.index.astype(float).values
    data = df.values

    return R_grid, Z_grid, data


def load_flux_surfaces(r_filename='eqgs2d_04_RPS_flux_surfaces.csv',
                      z_filename='eqgs2d_05_ZPS_flux_surfaces.csv'):
    """
    读取磁面数据 RPS(theta,psi) 和 ZPS(theta,psi)

    返回:
        tuple: (theta, R_surfaces, Z_surfaces)
    """
    df_r = pd.read_csv(r_filename, index_col=0)
    df_z = pd.read_csv(z_filename, index_col=0)

    theta = df_r.index.astype(float).values
    R_surfaces = df_r.values
    Z_surfaces = df_z.values

    return theta, R_surfaces, Z_surfaces


def load_separatrix(filename='eqgs2d_07_separatrix.csv'):
    """读取分界面数据"""
    if not os.path.exists(filename):
        return None

    df = pd.read_csv(filename)
    return df


def load_parameters(filename='eqgs2d_08_parameters.csv'):
    """读取关键参数"""
    df = pd.read_csv(filename, skipinitialspace=True)
    params = {}
    for _, row in df.iterrows():
        param_name = row['parameter'].strip() if isinstance(row['parameter'], str) else row['parameter']
        params[param_name] = {
            'value': float(row['value']) if isinstance(row['value'], str) else row['value'],
            'unit': row['unit'].strip() if isinstance(row['unit'], str) else row['unit'],
            'description': row['description'].strip() if isinstance(row['description'], str) else row['description']
        }
    return params


def plot_psirz_contour(R_grid=None, Z_grid=None, PSIRZ=None,
                       levels=20, figsize=(12, 10), save_fig=False):
    """
    绘制 PSIRZ 等高线图

    参数:
        R_grid, Z_grid, PSIRZ: 如果为 None，则自动加载
        levels: 等高线数量
        figsize: 图形大小
        save_fig: 是否保存图像
    """
    if R_grid is None or Z_grid is None or PSIRZ is None:
        R_grid, Z_grid, PSIRZ = load_psirz_grid()

    # 创建 2D 网格
    R, Z = np.meshgrid(R_grid, Z_grid)

    fig, ax = plt.subplots(figsize=figsize)

    # 绘制等高线
    contour = ax.contour(R, Z, PSIRZ, levels=levels, colors='blue', linewidths=0.5)
    ax.clabel(contour, inline=True, fontsize=8)

    # 绘制填充等高线
    contourf = ax.contourf(R, Z, PSIRZ, levels=levels, cmap='viridis', alpha=0.6)
    fig.colorbar(contourf, ax=ax, label='PSIRZ (Wb)')

    # 尝试加载并绘制分界面
    sep = load_separatrix()
    if sep is not None:
        ax.plot(sep['RSU_m'], sep['ZSU_m'], 'r-', linewidth=2, label='Separatrix (upper)')
        ax.plot(sep['RSW_m'], sep['ZSW_m'], 'r--', linewidth=2, label='Separatrix (lower)')

    # 加载参数并绘制磁轴
    params = load_parameters()
    if 'RAXIS' in params and 'ZAXIS' in params:
        raxis = params['RAXIS']['value']
        zaxis = params['ZAXIS']['value']
        ax.plot(raxis, zaxis, 'r*', markersize=15, label='Magnetic axis')

    ax.set_xlabel('R (m)', fontsize=12)
    ax.set_ylabel('Z (m)', fontsize=12)
    ax.set_title('Poloidal Flux Contours (PSIRZ)', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_aspect('equal')

    if save_fig:
        plt.savefig('eqgs2d_psirz_contour.png', dpi=150, bbox_inches='tight')
        print("图像已保存: eqgs2d_psirz_contour.png")

    return fig, ax


def plot_flux_surfaces(theta=None, R_surfaces=None, Z_surfaces=None,
                      plot_every=5, figsize=(12, 10), save_fig=False):
    """
    绘制磁面

    参数:
        theta, R_surfaces, Z_surfaces: 如果为 None，则自动加载
        plot_every: 每隔几条磁面绘制一条
        figsize: 图形大小
        save_fig: 是否保存图像
    """
    if theta is None or R_surfaces is None or Z_surfaces is None:
        theta, R_surfaces, Z_surfaces = load_flux_surfaces()

    fig, ax = plt.subplots(figsize=figsize)

    n_surfaces = R_surfaces.shape[1]

    # 绘制磁面（每隔 plot_every 条）
    for i in range(0, n_surfaces, plot_every):
        r = R_surfaces[:, i]
        z = Z_surfaces[:, i]
        ax.plot(r, z, '-', linewidth=1.5, alpha=0.7, label=f'Surface {i+1}' if i < 10 else None)

    # 绘制分界面
    sep = load_separatrix()
    if sep is not None:
        ax.plot(sep['RSU_m'], sep['ZSU_m'], 'r-', linewidth=3, label='Separatrix')
        ax.plot(sep['RSW_m'], sep['ZSW_m'], 'r-', linewidth=3)

    # 绘制磁轴
    params = load_parameters()
    if 'RAXIS' in params and 'ZAXIS' in params:
        raxis = params['RAXIS']['value']
        zaxis = params['ZAXIS']['value']
        ax.plot(raxis, zaxis, 'r*', markersize=15, label='Magnetic axis')

    ax.set_xlabel('R (m)', fontsize=12)
    ax.set_ylabel('Z (m)', fontsize=12)
    ax.set_title(f'Flux Surfaces (every {plot_every} surfaces)', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=8, ncol=2)
    ax.set_aspect('equal')

    if save_fig:
        plt.savefig('eqgs2d_flux_surfaces.png', dpi=150, bbox_inches='tight')
        print("图像已保存: eqgs2d_flux_surfaces.png")

    return fig, ax


def plot_all_2d_data(save_fig=False):
    """绘制所有 2D 数据的组合图"""

    # 加载数据
    R_grid, Z_grid, PSIRZ = load_psirz_grid()
    _, _, DPSIDR = load_derivative_grid('eqgs2d_02_DPSIDR_grid.csv')
    _, _, DPSIDZ = load_derivative_grid('eqgs2d_03_DPSIDZ_grid.csv')
    theta, R_surfaces, Z_surfaces = load_flux_surfaces()
    params = load_parameters()

    # 创建 2x2 子图
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    R, Z = np.meshgrid(R_grid, Z_grid)

    # 子图 1: PSIRZ 等高线
    ax = axes[0, 0]
    contour = ax.contour(R, Z, PSIRZ, levels=20, colors='blue', linewidths=0.5)
    contourf = ax.contourf(R, Z, PSIRZ, levels=20, cmap='viridis', alpha=0.6)
    fig.colorbar(contourf, ax=ax, label='PSIRZ (Wb)')
    ax.set_title('Poloidal Flux (PSIRZ)')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    ax.set_aspect('equal')

    # 子图 2: DPSIDR
    ax = axes[0, 1]
    contourf = ax.contourf(R, Z, DPSIDR, levels=20, cmap='RdBu_r')
    fig.colorbar(contourf, ax=ax, label='∂ψ/∂R (Wb/m)')
    ax.set_title('Radial Flux Gradient (∂ψ/∂R)')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    ax.set_aspect('equal')

    # 子图 3: DPSIDZ
    ax = axes[1, 0]
    contourf = ax.contourf(R, Z, DPSIDZ, levels=20, cmap='RdBu_r')
    fig.colorbar(contourf, ax=ax, label='∂ψ/∂Z (Wb/m)')
    ax.set_title('Vertical Flux Gradient (∂ψ/∂Z)')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    ax.set_aspect('equal')

    # 子图 4: 磁面
    ax = axes[1, 1]
    n_surfaces = R_surfaces.shape[1]
    for i in range(0, n_surfaces, 5):
        ax.plot(R_surfaces[:, i], Z_surfaces[:, i], '-', linewidth=1, alpha=0.7)

    # 绘制分界面
    sep = load_separatrix()
    if sep is not None:
        ax.plot(sep['RSU_m'], sep['ZSU_m'], 'r-', linewidth=2, label='Separatrix')
        ax.plot(sep['RSW_m'], sep['ZSW_m'], 'r-', linewidth=2)

    # 绘制磁轴
    if 'RAXIS' in params and 'ZAXIS' in params:
        raxis = params['RAXIS']['value']
        zaxis = params['ZAXIS']['value']
        ax.plot(raxis, zaxis, 'r*', markersize=15, label='Magnetic axis')

    ax.set_title('Flux Surfaces')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_fig:
        plt.savefig('eqgs2d_all_plots.png', dpi=150, bbox_inches='tight')
        print("图像已保存: eqgs2d_all_plots.png")

    return fig, axes


def print_parameters():
    """打印关键参数"""
    params = load_parameters()

    print("\n" + "="*70)
    print("EQGS2D 关键参数")
    print("="*70)

    for key, val in params.items():
        print(f"{key:15s} = {val['value']:12.4e} {val['unit']:5s} | {val['description']}")

    print("="*70)


# 示例用法
if __name__ == '__main__':
    import glob

    print("="*70)
    print("EQGS2D CSV 数据读取和可视化工具")
    print("="*70)

    # 检查文件是否存在
    csv_files = glob.glob('eqgs2d_*.csv')
    if not csv_files:
        print("\n错误: 未找到 2D CSV 文件!")
        print("请先运行: ./eq <in/eq.ITER01_export2d.in")
        print("或在 eq 程序中执行命令 'ce' (计算后导出2D数据)")
        exit(1)

    print(f"\n找到 {len(csv_files)} 个 CSV 文件\n")

    # 打印参数
    print_parameters()

    # 绘制单独的 PSIRZ 等高线图
    print("\n绘制 PSIRZ 等高线图...")
    plot_psirz_contour(save_fig=True)

    # 绘制磁面
    print("绘制磁面...")
    plot_flux_surfaces(save_fig=True)

    # 绘制组合图
    print("绘制所有 2D 数据组合图...")
    plot_all_2d_data(save_fig=True)

    plt.show()

    print("\n" + "="*70)
    print("完成!")
    print("="*70)
    print("\n可用函数:")
    print("  - load_psirz_grid()          : 读取 PSIRZ(R,Z) 网格")
    print("  - load_flux_surfaces()        : 读取磁面数据")
    print("  - load_separatrix()           : 读取分界面")
    print("  - load_parameters()           : 读取关键参数")
    print("  - plot_psirz_contour()        : 绘制 PSIRZ 等高线")
    print("  - plot_flux_surfaces()        : 绘制磁面")
    print("  - plot_all_2d_data()          : 绘制所有 2D 数据")
