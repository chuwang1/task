# G-File (EQDSK) 使用说明

## 问题描述

当使用 `MODELG=5` 从 EQDSK 格式文件（如 `gfile_efit`）读取平衡数据时，保存操作（`s` 命令）会覆盖原始输入文件，因为程序默认使用 `KNAMEQ` 参数指定的文件作为输入和输出。

## 解决方案

在保存前修改 `KNAMEQ` 参数，将输出重定向到新文件。

---

## 方法 1: 不保存数据（推荐用于快速查看）

**输入文件**: `in/eq.CFEDR_nosave.in`

```
0
f
eq.CFEDR.gs
c
p
 &eq
   MODELG=5
   KNAMEQ='/Users/dengxiaoya/TASK/latest/task/eq/in/gfile_efit'
 &end
l
g
c
x
q
```

**运行方式**:
```bash
./eq <in/eq.CFEDR_nosave.in
```

**说明**:
- 加载 `gfile_efit` 文件
- 显示图形
- 不执行保存操作，直接退出
- 原始文件保持不变

---

## 方法 2: 保存到新文件（推荐用于数据处理）

**输入文件**: `in/eq.CFEDR_save.in`

```
0
f
eq.CFEDR.gs
c
p
 &eq
   MODELG=5
   KNAMEQ='/Users/dengxiaoya/TASK/latest/task/eq/in/gfile_efit'
 &end
l
knameq='/Users/dengxiaoya/TASK/latest/task/eq/in/gfile_efit_output'
g
c
x
s
q
```

**运行方式**:
```bash
./eq <in/eq.CFEDR_save.in
```

**说明**:
- 加载 `gfile_efit` 文件（输入）
- 在主菜单中修改 `KNAMEQ` 参数指向新文件名
- **导出 2D 数据**到 `eqgs2d_*.csv` 文件（10 个文件）
- **导出 1D 剖面**到 `eqgs1d_*.csv` 文件（23 个文件）
- 保存图形到 `eq.CFEDR.gs`
- 保存数据到 `gfile_efit_output` 文件（输出）
- 原始 `gfile_efit` 文件保持不变

---

## 方法 3: 使用环境变量（适合脚本自动化）

```bash
export GSGDP="0f"
export EQOUT="in/gfile_efit_output"
./eq <<EOF
c
p
 &eq
   MODELG=5
   KNAMEQ='in/gfile_efit'
 &end
l
knameq='${EQOUT}'
g
c
x
s
q
EOF
```

---

## 方法 4: 交互式使用（手动控制）

```bash
./eq
```

然后按以下顺序输入命令：

```
0                # 图形模式选择（0 = 无图形输出）
f                # 保存图形到文件
eq.CFEDR.gs      # 图形文件名
c                # 继续
p                # 参数设置
 &eq
   MODELG=5
   KNAMEQ='in/gfile_efit'
 &end
l                # 加载数据
```

此时数据已加载。**重要**: 在执行保存前修改输出文件名：

```
knameq='in/gfile_efit_output'    # 修改保存文件名
g                                 # 图形菜单
c                                 # 计算
x                                 # 退出图形菜单
s                                 # 保存（现在会保存到新文件）
q                                 # 退出
```

---

## 参数说明

### MODELG 参数

根据 `eq.pdf` 文档（第9页）：

```
MODELG: Control plasma geometry model
  0: Slab geometry
  1: Cylindrical geometry
  2: Toroidal geometry
  3: TASK/EQ output geometry
  4: VMEC output geometry
  5: EQDSK output geometry       ← 用于读取 g-file
  6: Boozer output geometry
```

### KNAMEQ 参数

根据 `eq.pdf` 文档（第10页）：

```
KNAMEQ: Filename of equilibrium data
```

- **加载时**: 指定要读取的输入文件
- **保存时**: 指定要写入的输出文件
- **默认值**: `'eqdata'`

---

## 文件格式

### 输入格式 (EQDSK/G-file)

EQDSK 格式是托卡马克平衡数据的标准格式，由 `eq-eqdsk.f` 中的 `EQDSKR` 子程序读取。

包含的数据：
- 磁通网格 PSIRZ(R,Z)
- 压强剖面 P(ψ)
- FF' 和 P' 剖面
- 磁轴位置
- 边界点坐标
- 限制器坐标

### 输出格式 (TASK/EQ)

由 `eqfile.f` 中的 `EQSAVE` 子程序保存，使用 Fortran 二进制格式（非 EQDSK 格式）。

包含的数据：
- 基本参数：RR, BB, RIP
- 网格数据：PSIRZ, RG, ZG
- 剖面数据：PSIPS, PPPS, TTPS, TEPS, OMPS
- 计算结果：PSI, DELPSI, HJT, HJTRZ
- 几何参数：RAXIS, ZAXIS, PSITA, PSIPA, PSI0

---

## 常见使用场景

### 场景 1: 从 EFIT 读取数据并可视化

```bash
./eq <in/eq.CFEDR_nosave.in
gsview eq.CFEDR.gs
```

### 场景 2: 从 EFIT 读取，导出 CSV 数据

```bash
cat > in/eq.CFEDR_export.in <<'EOF'
0
f
eq.CFEDR.gs
c
p
 &eq
   MODELG=5
   KNAMEQ='in/gfile_efit'
 &end
l
g
c
ce
se
x
s
q
EOF

./eq <in/eq.CFEDR_export.in
```

此命令会：
- 加载 `gfile_efit`
- 导出 2D 数据（`ce` 命令）
- 导出 1D 剖面（`se` 命令）
- 保存图形文件

### 场景 3: 读取并修改参数后保存

如果需要调整参数并保存为 TASK/EQ 格式：

```bash
./eq
# 加载 g-file
l
# 读取时使用原始文件
knameq='in/gfile_efit'
# ... 进行计算或修改 ...
# 保存前修改输出文件名
knameq='in/gfile_modified'
s
```

---

## 注意事项

1. **文件路径**:
   - 可以使用绝对路径或相对路径
   - 相对路径是相对于 `eq` 可执行文件的工作目录

2. **文件格式不同**:
   - 输入: EQDSK 格式（ASCII，标准格式）
   - 输出（`s` 命令）: TASK/EQ 格式（二进制，专用格式）
   - 如果需要 EQDSK 格式输出，需要使用其他工具或修改代码

3. **备份原始文件**:
   ```bash
   cp in/gfile_efit in/gfile_efit_backup
   ```

4. **检查文件是否被覆盖**:
   ```bash
   ls -lh in/gfile_efit*
   stat in/gfile_efit  # 查看文件修改时间
   ```

---

## 快速参考

| 操作 | 命令序列 | 输入文件会被覆盖？ |
|------|----------|-------------------|
| 只读取和查看 | `l → g → c → x → q` | 否 |
| 读取后直接保存 | `l → s → q` | **是** |
| 读取后保存到新文件 | `l → knameq='新文件' → s → q` | 否 |
| 读取并导出 CSV | `l → g → c → ce → se → x → q` | 否 |

---

## 相关文件

- **输入文件**: `in/eq.CFEDR.in` （原始，会覆盖）
- **安全输入文件**: `in/eq.CFEDR_nosave.in` （不保存）
- **带保存的输入文件**: `in/eq.CFEDR_save.in` （保存到新文件）
- **源代码**:
  - `eqfile.f` - 保存/加载子程序
  - `eq-eqdsk.f` - EQDSK 格式读写
  - `eqinit.f` - 参数初始化
- **文档**: `doc/eq.pdf`

---

## 总结

**推荐做法**:

✅ **读取后不保存**: 使用 `in/eq.CFEDR_nosave.in`

✅ **读取后保存到新文件**: 使用 `in/eq.CFEDR_save.in` 或手动修改 `KNAMEQ`

❌ **不推荐**: 使用原始 `in/eq.CFEDR.in` 并执行保存命令（会覆盖输入文件）

---

**创建日期**: 2026-01-14
**版本**: 1.0
