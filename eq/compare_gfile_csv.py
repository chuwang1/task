#!/usr/bin/env python3
"""
对比 gfile_efit 原始数据和导出的 CSV 数据

验证数据导出的准确性
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

# 导入用户提供的 geqdsk 读取器
sys.path.insert(0, '.')

def read_gfile_efit(filename):
    """使用 geqdsk 类读取 gfile_efit"""
    import re

    # 简化版 geqdsk 读取器（只读取我们需要的数据）
    class SimpleGEQDSK:
        def __init__(self, file):
            f = open(file, 'r', encoding='utf-8')

            # 第1行：header, nw, nh
            header1 = re.compile('^\s*(.*)\s+(\d+)\s+(\d+)\s*$', re.I)
            header2 = re.compile('^\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)')

            line = f.readline()
            pat = header1.match(line)
            self.header = pat.group(1)
            self.nw = int(pat.group(2))
            self.nh = int(pat.group(3))

            # 第2行：rdim, zdim, rcentr, rleft, zmid
            line = f.readline()
            line = re.sub('([^Ee])-', '\\1 -', line)
            pat = header2.match(line)
            self.rdim = float(pat.group(1))
            self.zdim = float(pat.group(2))
            self.rcentr = float(pat.group(3))
            self.rleft = float(pat.group(4))
            self.zmid = float(pat.group(5))

            # 第3行：rmaxis, zmaxis, psimag, psibdy, bcentr
            line = f.readline()
            line = re.sub('([^Ee])-', '\\1 -', line)
            pat = header2.match(line)
            self.rmaxis = float(pat.group(1))
            self.zmaxis = float(pat.group(2))
            self.psimag = float(pat.group(3))
            self.psibdy = float(pat.group(4))
            self.bcentr = float(pat.group(5))

            # 第4行：current
            line = f.readline()
            line = re.sub('([^Ee])-', '\\1 -', line)
            pat = header2.match(line)
            self.currentA = float(pat.group(1))

            # 跳过第5行
            f.readline()

            # 读取剩余数据
            lines = f.read()
            lines = re.sub('\n', ' ', lines)
            lines = re.sub('(\d)\-', '\\1 -', lines)
            fields = re.split('\s+', lines)
            numbers = []
            for n in fields:
                try:
                    numbers.append(float(n))
                except:
                    pass

            nw, nh = self.nw, self.nh

            # 读取 1D 剖面
            self.fpol = numbers[0:nw]
            self.pres = numbers[nw:2*nw]
            self.ffprime = numbers[2*nw:3*nw]
            self.pprime = numbers[3*nw:4*nw]

            # 读取 2D PSIRZ
            ia = 4*nw
            self.psirz = np.zeros([nw, nh])
            for ir in range(nw):
                self.psirz[ir, :] = numbers[ia:ia+nh]
                ia = ia + nh

            # 读取 q 剖面
            self.qpsi = numbers[ia:ia+nw]
            ia = ia + nw

            # 读取边界点
            self.nbbbs = int(numbers[ia])
            ia = ia + 1
            self.limitr = int(numbers[ia])
            ia = ia + 1

            self.rbbbs = []
            self.zbbbs = []
            for i in range(self.nbbbs):
                self.rbbbs.append(numbers[ia])
                self.zbbbs.append(numbers[ia+1])
                ia = ia + 2

            f.close()

    return SimpleGEQDSK(filename)

def load_csv_data():
    """从 CSV 文件加载数据"""
    from read_eqgs2d_csv import load_psirz_grid, load_parameters

    R_grid, Z_grid, PSIRZ = load_psirz_grid()
    params = load_parameters()

    return R_grid, Z_grid, PSIRZ, params

def compare_data():
    """对比原始 gfile 和 CSV 数据"""

    print("="*70)
    print("数据对比：gfile_efit vs CSV")
    print("="*70)
    print()

    # 读取原始 gfile
    print("1. 读取原始 gfile_efit...")
    gfile = read_gfile_efit('in/gfile_efit')
    print(f"   网格尺寸: {gfile.nw} × {gfile.nh}")
    print(f"   磁轴位置: R = {gfile.rmaxis:.4f} m, Z = {gfile.zmaxis:.6e} m")
    print(f"   磁通范围: {gfile.psimag:.4f} ~ {gfile.psibdy:.4f} Wb")
    print()

    # 读取 CSV 数据
    print("2. 读取 CSV 数据...")
    R_csv, Z_csv, PSIRZ_csv, params = load_csv_data()
    print(f"   网格尺寸: {len(R_csv)} × {len(Z_csv)}")
    print(f"   磁轴位置: R = {params['RAXIS']['value']:.4f} m, Z = {params['ZAXIS']['value']:.6e} m")
    print()

    # 构建 gfile 的网格坐标
    R_gfile = np.linspace(gfile.rleft, gfile.rleft + gfile.rdim, gfile.nw)
    Z_gfile = np.linspace(gfile.zmid - 0.5*gfile.zdim, gfile.zmid + 0.5*gfile.zdim, gfile.nh)

    # 对比网格
    print("3. 对比网格坐标...")
    r_diff = np.abs(R_gfile - R_csv).max()
    z_diff = np.abs(Z_gfile - Z_csv).max()
    print(f"   R 网格最大差异: {r_diff:.2e} m")
    print(f"   Z 网格最大差异: {z_diff:.2e} m")

    if r_diff < 1e-6 and z_diff < 1e-6:
        print("   ✓ 网格坐标完全一致")
    else:
        print("   ✗ 网格坐标有差异")
    print()

    # 对比 PSIRZ 数据
    print("4. 对比 PSIRZ 数据...")

    # gfile 的 psirz 是 [nw, nh]，需要转置为 [nh, nw] 才能匹配 (Z, R) 布局
    # 但是我们的 CSV 是 [nw, nh] 格式 (R, Z)
    # 需要检查索引顺序

    # CSV: PSIRZ_csv[j, i] 对应 (R_csv[i], Z_csv[j])
    # gfile: gfile.psirz[i, j] 对应 (R_gfile[i], Z_gfile[j])

    # 直接对比
    psirz_diff = np.abs(gfile.psirz - PSIRZ_csv.T)
    max_diff = psirz_diff.max()
    mean_diff = psirz_diff.mean()
    rel_diff = max_diff / (np.abs(gfile.psirz).max())

    print(f"   PSIRZ 最大绝对差异: {max_diff:.2e} Wb")
    print(f"   PSIRZ 平均绝对差异: {mean_diff:.2e} Wb")
    print(f"   PSIRZ 最大相对差异: {rel_diff:.2e}")

    if rel_diff < 1e-6:
        print("   ✓ PSIRZ 数据完全一致（相对误差 < 1e-6）")
    elif rel_diff < 1e-4:
        print("   ✓ PSIRZ 数据基本一致（相对误差 < 1e-4）")
    else:
        print("   ✗ PSIRZ 数据有明显差异")
    print()

    # 统计对比
    print("5. 统计对比...")
    print(f"   gfile PSIRZ:  min={gfile.psirz.min():.4f}, max={gfile.psirz.max():.4f}, mean={gfile.psirz.mean():.4f}")
    print(f"   CSV PSIRZ:    min={PSIRZ_csv.min():.4f}, max={PSIRZ_csv.max():.4f}, mean={PSIRZ_csv.mean():.4f}")
    print()

    # 可视化对比
    print("6. 生成对比图...")
    create_comparison_plot(R_gfile, Z_gfile, gfile.psirz,
                          R_csv, Z_csv, PSIRZ_csv.T,
                          gfile, params)

    print()
    print("="*70)
    print("对比完成！")
    print("="*70)

    return rel_diff < 1e-4

def create_comparison_plot(R_gfile, Z_gfile, PSIRZ_gfile,
                          R_csv, Z_csv, PSIRZ_csv,
                          gfile, params):
    """创建对比图"""

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    R2D_gfile, Z2D_gfile = np.meshgrid(R_gfile, Z_gfile)
    R2D_csv, Z2D_csv = np.meshgrid(R_csv, Z_csv)

    # 左上：gfile 原始数据
    ax = axes[0, 0]
    cf = ax.contourf(R2D_gfile, Z2D_gfile, PSIRZ_gfile.T, levels=30, cmap='viridis')
    plt.colorbar(cf, ax=ax, label='PSIRZ (Wb)')
    ax.plot([gfile.rmaxis], [gfile.zmaxis], 'r*', markersize=15, label='Magnetic axis')
    ax.plot(gfile.rbbbs, gfile.zbbbs, 'r-', linewidth=2, label='Boundary')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    ax.set_title('gfile_efit 原始数据')
    ax.set_aspect('equal')
    ax.legend()

    # 中上：CSV 数据
    ax = axes[0, 1]
    cf = ax.contourf(R2D_csv, Z2D_csv, PSIRZ_csv, levels=30, cmap='viridis')
    plt.colorbar(cf, ax=ax, label='PSIRZ (Wb)')
    ax.plot([params['RAXIS']['value']], [params['ZAXIS']['value']], 'r*', markersize=15, label='Magnetic axis')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    ax.set_title('CSV 导出数据')
    ax.set_aspect('equal')
    ax.legend()

    # 右上：差异图
    ax = axes[0, 2]
    diff = PSIRZ_gfile.T - PSIRZ_csv
    cf = ax.contourf(R2D_csv, Z2D_csv, diff, levels=30, cmap='RdBu_r')
    plt.colorbar(cf, ax=ax, label='Difference (Wb)')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    ax.set_title(f'差异 (max={np.abs(diff).max():.2e} Wb)')
    ax.set_aspect('equal')

    # 左下：沿 R 方向的剖面对比（Z=0 附近）
    ax = axes[1, 0]
    iz_mid = len(Z_csv) // 2
    ax.plot(R_gfile, PSIRZ_gfile[:, iz_mid], 'b-', linewidth=2, label='gfile')
    ax.plot(R_csv, PSIRZ_csv[iz_mid, :], 'r--', linewidth=2, label='CSV')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('PSIRZ (Wb)')
    ax.set_title(f'沿 R 方向剖面 (Z≈{Z_csv[iz_mid]:.2f} m)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 中下：沿 Z 方向的剖面对比（R=磁轴位置附近）
    ax = axes[1, 1]
    ir_axis = np.argmin(np.abs(R_csv - params['RAXIS']['value']))
    ax.plot(Z_gfile, PSIRZ_gfile[ir_axis, :], 'b-', linewidth=2, label='gfile')
    ax.plot(Z_csv, PSIRZ_csv[:, ir_axis], 'r--', linewidth=2, label='CSV')
    ax.set_xlabel('Z (m)')
    ax.set_ylabel('PSIRZ (Wb)')
    ax.set_title(f'沿 Z 方向剖面 (R≈{R_csv[ir_axis]:.2f} m)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 右下：差异统计直方图
    ax = axes[1, 2]
    diff_flat = diff.flatten()
    ax.hist(diff_flat, bins=50, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Difference (Wb)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'差异分布 (mean={diff_flat.mean():.2e} Wb)')
    ax.axvline(0, color='r', linestyle='--', linewidth=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('gfile_csv_comparison.png', dpi=150, bbox_inches='tight')
    print("   对比图已保存: gfile_csv_comparison.png")

    # 也绘制等高线叠加对比
    fig2, ax = plt.subplots(figsize=(12, 10))
    ax.contour(R2D_gfile, Z2D_gfile, PSIRZ_gfile.T, levels=20, colors='blue',
              linewidths=1.5, linestyles='-', alpha=0.7, label='gfile')
    ax.contour(R2D_csv, Z2D_csv, PSIRZ_csv, levels=20, colors='red',
              linewidths=1.5, linestyles='--', alpha=0.7, label='CSV')
    ax.plot([gfile.rmaxis], [gfile.zmaxis], 'b*', markersize=15, label='gfile axis')
    ax.plot([params['RAXIS']['value']], [params['ZAXIS']['value']], 'r*', markersize=15, label='CSV axis')
    ax.plot(gfile.rbbbs, gfile.zbbbs, 'b-', linewidth=2, alpha=0.5, label='gfile boundary')
    ax.set_xlabel('R (m)', fontsize=12)
    ax.set_ylabel('Z (m)', fontsize=12)
    ax.set_title('等高线叠加对比：蓝色实线=gfile，红色虚线=CSV', fontsize=14)
    ax.set_aspect('equal')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.savefig('gfile_csv_contour_overlay.png', dpi=150, bbox_inches='tight')
    print("   等高线叠加图已保存: gfile_csv_contour_overlay.png")

if __name__ == '__main__':
    success = compare_data()

    if success:
        print()
        print("✓ 验证通过：CSV 导出数据与原始 gfile 数据一致！")
        sys.exit(0)
    else:
        print()
        print("✗ 验证失败：CSV 数据与原始 gfile 有差异")
        sys.exit(1)
