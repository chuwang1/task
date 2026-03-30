#!/usr/bin/env python3
"""
从 TASK/EQ 二进制数据文件中导出 2D CSV 数据

这个脚本读取由 `EQSAVE` 保存的二进制文件，并导出 2D 数据到 CSV 格式。
用于从 EQDSK 文件加载后无法使用 'ce' 命令的情况。

使用方法:
    python export_2d_from_eqdata.py <eqdata_file>

例如:
    python export_2d_from_eqdata.py in/gfile_efit_output
"""

import struct
import numpy as np
import pandas as pd
import sys
import os

def read_fortran_real8(f):
    """读取 Fortran REAL*8 (双精度浮点数)"""
    data = f.read(8)
    if len(data) < 8:
        return None
    return struct.unpack('d', data)[0]

def read_fortran_int(f):
    """读取 Fortran INTEGER"""
    data = f.read(4)
    if len(data) < 4:
        return None
    return struct.unpack('i', data)[0]

def read_fortran_record_header(f):
    """读取 Fortran 记录头（4字节长度标记）"""
    data = f.read(4)
    if len(data) < 4:
        return None
    return struct.unpack('i', data)[0]

def read_fortran_record(f, dtype='d'):
    """读取一个完整的 Fortran 记录"""
    # 读取记录头
    rec_len = read_fortran_record_header(f)
    if rec_len is None:
        return None

    # 计算元素数量
    if dtype == 'd':  # REAL*8
        elem_size = 8
    elif dtype == 'i':  # INTEGER
        elem_size = 4
    else:
        raise ValueError(f"Unsupported dtype: {dtype}")

    n_elem = rec_len // elem_size

    # 读取数据
    data = []
    for _ in range(n_elem):
        if dtype == 'd':
            val = read_fortran_real8(f)
        elif dtype == 'i':
            val = read_fortran_int(f)
        if val is None:
            return None
        data.append(val)

    # 读取记录尾
    rec_len_end = read_fortran_record_header(f)
    if rec_len != rec_len_end:
        print(f"Warning: Record length mismatch: {rec_len} != {rec_len_end}")

    return data

def export_2d_from_eqdata(filename):
    """
    从 TASK/EQ 二进制文件中导出 2D 数据

    根据 eqfile.f 中的 EQSAVE 子程序结构读取：
      WRITE(21) RR,BB,RIP
      WRITE(21) NRGMAX,NZGMAX
      WRITE(21) (RG(NRG),NRG=1,NRGMAX)
      WRITE(21) (ZG(NZG),NZG=1,NZGMAX)
      WRITE(21) ((PSIRZ(NRG,NZG),NRG=1,NRGMAX),NZG=1,NZGMAX)
      ...
    """

    if not os.path.exists(filename):
        print(f"Error: File not found: {filename}")
        return False

    print(f"Reading TASK/EQ data from: {filename}")

    try:
        with open(filename, 'rb') as f:
            # Record 1: RR, BB, RIP
            rec = read_fortran_record(f, 'd')
            if rec is None or len(rec) < 3:
                print("Error: Failed to read record 1 (RR, BB, RIP)")
                return False
            RR, BB, RIP = rec[0], rec[1], rec[2]
            print(f"RR = {RR:.4f} m, BB = {BB:.4f} T, RIP = {RIP:.4f} MA")

            # Record 2: NRGMAX, NZGMAX
            rec = read_fortran_record(f, 'i')
            if rec is None or len(rec) < 2:
                print("Error: Failed to read record 2 (NRGMAX, NZGMAX)")
                return False
            NRGMAX, NZGMAX = rec[0], rec[1]
            print(f"Grid size: NRGMAX = {NRGMAX}, NZGMAX = {NZGMAX}")

            # Record 3: RG(NRG)
            RG = read_fortran_record(f, 'd')
            if RG is None or len(RG) != NRGMAX:
                print(f"Error: Failed to read RG array (expected {NRGMAX} elements)")
                return False

            # Record 4: ZG(NZG)
            ZG = read_fortran_record(f, 'd')
            if ZG is None or len(ZG) != NZGMAX:
                print(f"Error: Failed to read ZG array (expected {NZGMAX} elements)")
                return False

            # Record 5: PSIRZ(NRG,NZG) - 2D array
            psirz_data = read_fortran_record(f, 'd')
            if psirz_data is None or len(psirz_data) != NRGMAX * NZGMAX:
                print(f"Error: Failed to read PSIRZ array (expected {NRGMAX*NZGMAX} elements)")
                return False

            # 重塑为 2D 数组 (Fortran 列优先)
            PSIRZ = np.array(psirz_data).reshape((NRGMAX, NZGMAX), order='F')

            print(f"Successfully read PSIRZ grid: {PSIRZ.shape}")

            # 导出到 CSV
            export_psirz_grid(RG, ZG, PSIRZ)

            return True

    except Exception as e:
        print(f"Error reading file: {e}")
        import traceback
        traceback.print_exc()
        return False

def export_psirz_grid(R_grid, Z_grid, PSIRZ):
    """导出 PSIRZ 网格到 CSV 文件"""

    filename = 'eqgs2d_01_PSIRZ_grid.csv'

    # 创建 DataFrame (Z 为行，R 为列)
    df = pd.DataFrame(
        PSIRZ.T,  # 转置使 Z 为行，R 为列
        index=Z_grid,
        columns=R_grid
    )

    # 保存到 CSV
    with open(filename, 'w') as f:
        f.write('Z_m / R_m')
        for r in R_grid:
            f.write(f',{r:.8E}')
        f.write('\n')

        for i, z in enumerate(Z_grid):
            f.write(f'{z:.8E}')
            for j in range(len(R_grid)):
                f.write(f',{PSIRZ[j,i]:.8E}')
            f.write('\n')

    print(f"✓ Exported: {filename} ({df.shape[0]}×{df.shape[1]} grid)")

    # 也创建参数文件
    export_parameters(R_grid, Z_grid, PSIRZ)

def export_parameters(R_grid, Z_grid, PSIRZ):
    """导出基本参数到 CSV 文件"""

    filename = 'eqgs2d_08_parameters_basic.csv'

    # 找到最大值和最小值位置
    psi_min = PSIRZ.min()
    psi_max = PSIRZ.max()

    params = {
        'NRGMAX': len(R_grid),
        'NZGMAX': len(Z_grid),
        'R_min': R_grid[0],
        'R_max': R_grid[-1],
        'Z_min': Z_grid[0],
        'Z_max': Z_grid[-1],
        'PSI_min': psi_min,
        'PSI_max': psi_max,
    }

    with open(filename, 'w') as f:
        f.write('parameter,value,unit,description\n')
        f.write(f"NRGMAX,{params['NRGMAX']:.8E},-,R pts\n")
        f.write(f"NZGMAX,{params['NZGMAX']:.8E},-,Z pts\n")
        f.write(f"R_min,{params['R_min']:.8E},m,min R\n")
        f.write(f"R_max,{params['R_max']:.8E},m,max R\n")
        f.write(f"Z_min,{params['Z_min']:.8E},m,min Z\n")
        f.write(f"Z_max,{params['Z_max']:.8E},m,max Z\n")
        f.write(f"PSI_min,{params['PSI_min']:.8E},Wb,min flux\n")
        f.write(f"PSI_max,{params['PSI_max']:.8E},Wb,max flux\n")

    print(f"✓ Exported: {filename}")

    # 也导出网格向量
    export_grids(R_grid, Z_grid)

def export_grids(R_grid, Z_grid):
    """导出 R 和 Z 网格向量"""

    filename = 'eqgs2d_09_RZ_grids.csv'

    max_len = max(len(R_grid), len(Z_grid))

    with open(filename, 'w') as f:
        f.write('index,R_m,Z_m\n')
        for i in range(max_len):
            r_val = f'{R_grid[i]:.8E}' if i < len(R_grid) else ''
            z_val = f'{Z_grid[i]:.8E}' if i < len(Z_grid) else ''
            f.write(f'{i+1},{r_val},{z_val}\n')

    print(f"✓ Exported: {filename}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python export_2d_from_eqdata.py <eqdata_file>")
        print("")
        print("Example:")
        print("  python export_2d_from_eqdata.py in/gfile_efit_output")
        sys.exit(1)

    filename = sys.argv[1]

    print("="*70)
    print("TASK/EQ 2D Data Exporter")
    print("="*70)
    print("")

    success = export_2d_from_eqdata(filename)

    if success:
        print("")
        print("="*70)
        print("Export completed successfully!")
        print("="*70)
        print("")
        print("Exported files:")
        print("  - eqgs2d_01_PSIRZ_grid.csv")
        print("  - eqgs2d_08_parameters_basic.csv")
        print("  - eqgs2d_09_RZ_grids.csv")
        print("")
        print("You can now use read_eqgs2d_csv.py to visualize the data.")
    else:
        print("")
        print("="*70)
        print("Export failed!")
        print("="*70)
        sys.exit(1)
