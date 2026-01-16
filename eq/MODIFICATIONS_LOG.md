# EQ 程序修改日志

## 2026-01-14: 2D CSV 导出功能增强

### 问题

当使用 `MODELG=5` 从 EQDSK 格式文件（`gfile_efit`）加载数据后，无法使用 `ce` 命令导出 2D CSV 数据。

**原因分析**:
- `ce` 命令在 `eqgout.f` 中只在 `MODE=1` 时有效
- 从 EQDSK 加载数据后，`MSTAT` 被设置为 2（而不是 1）
- 因此 `CALL EQGOUT(MSTAT)` 会调用 `EQGOUT(2)`，导致 `ce` 命令被拒绝

### 解决方案

修改 `eqgout.f` 使 `CE` 命令在 `MODE=2` 时也能工作。

### 修改的文件

#### 1. `eqgout.f` (第 33-48 行)

**修改前**:
```fortran
      IF(K1.EQ.'C') THEN
         IF(MODE.EQ.1) THEN
            IF(K2.EQ.' ') THEN
               CALL EQGC2D
               CALL EQGC1D
            ELSEIF(K2.EQ.'1') THEN
               CALL EQGC1D
            ELSEIF(K2.EQ.'2') THEN
               CALL EQGC2D
            ELSEIF(K2.EQ.'E') THEN
               CALL EQGS2D_EXPORT
            ENDIF
         ELSE
            WRITE(6,*) 'XX: EQGOUT: NO DATA CREATED!'
         ENDIF
```

**修改后**:
```fortran
      IF(K1.EQ.'C') THEN
         IF(MODE.EQ.1 .OR. MODE.EQ.2) THEN
            IF(K2.EQ.' ') THEN
               IF(MODE.EQ.1) THEN
                  CALL EQGC2D
                  CALL EQGC1D
               ENDIF
            ELSEIF(K2.EQ.'1') THEN
               IF(MODE.EQ.1) CALL EQGC1D
            ELSEIF(K2.EQ.'2') THEN
               IF(MODE.EQ.1) CALL EQGC2D
            ELSEIF(K2.EQ.'E') THEN
               CALL EQGS2D_EXPORT
            ENDIF
         ELSE
            WRITE(6,*) 'XX: EQGOUT: NO DATA CREATED!'
         ENDIF
```

**关键变化**:
- 第 33 行: `IF(MODE.EQ.1)` → `IF(MODE.EQ.1 .OR. MODE.EQ.2)`
- 第 34-38 行: 添加 `IF(MODE.EQ.1)` 条件，使 `C`, `C1`, `C2` 命令只在 `MODE=1` 时调用绘图函数
- 第 43-44 行: `CE` 命令现在在 `MODE=1` 或 `MODE=2` 时都可以工作

**原理**:
- `EQGS2D_EXPORT` 只读取数据，不需要绘图功能
- `EQGC2D` 和 `EQGC1D` 需要计算状态（`MODE=1`），所以保留原有限制
- `MODE=2` 表示数据已从文件加载，可以导出但不能绘制计算图

#### 2. `in/eq.CFEDR_save.in`

**最终版本**:
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
ce
se
s
x
s
q
```

**说明**:
- `l`: 加载 EQDSK 文件（`MSTAT=2`）
- `knameq='...'`: 修改输出文件名（避免覆盖输入文件）
- `g`: 进入图形菜单
- `ce`: 导出 2D CSV 数据（**现在可以在 MODE=2 下工作**）
- `se`: 导出 1D CSV 数据
- `s`: 保存图形（10 页 2D 图）
- `x`: 退出图形菜单
- `s`: 保存 TASK/EQ 数据文件
- `q`: 退出程序

### 测试结果

运行 `./eq <in/eq.CFEDR_save.in` 成功导出：

**2D 数据** (10 个文件):
```
eqgs2d_01_PSIRZ_grid.csv              - 281 KB (R-Z 网格上的极向磁通)
eqgs2d_02_DPSIDR_grid.csv             - 281 KB (∂ψ/∂R 梯度)
eqgs2d_03_DPSIDZ_grid.csv             - 281 KB (∂ψ/∂Z 梯度)
eqgs2d_04_RPS_flux_surfaces.csv       - 54 KB  (磁面 R 坐标)
eqgs2d_05_ZPS_flux_surfaces.csv       - 54 KB  (磁面 Z 坐标)
eqgs2d_06_flux_surface_contours.csv   - 244 KB (详细磁面等高线)
eqgs2d_07_separatrix.csv              - 7.2 KB (分界面边界)
eqgs2d_08_parameters.csv              - 351 B  (关键参数)
eqgs2d_09_RZ_grids.csv                - 4.9 KB (R 和 Z 网格向量)
eqgs2d_10_PSIP_profile.csv            - 1.1 KB (归一化磁通剖面)
```

**1D 剖面** (23 个文件):
```
eqgs1d_01_PPS.csv 到 eqgs1d_23_AVEGV.csv
```

**总计**: 33 个 CSV 文件，约 1.2 MB

### 向后兼容性

✅ 此修改完全向后兼容：
- 原有的 `MODE=1`（运行计算后）功能保持不变
- `C`, `C1`, `C2` 命令仍然只在 `MODE=1` 时工作
- 只有 `CE` 命令的可用范围被扩展到 `MODE=2`
- 不影响其他任何功能

### 受影响的使用场景

**之前无法工作** (现已修复):
```bash
./eq
l                    # 加载 EQDSK 文件 → MSTAT=2
g
ce                   # ❌ 错误: "NO DATA CREATED!"
```

**现在可以工作**:
```bash
./eq
l                    # 加载 EQDSK 文件 → MSTAT=2
g
ce                   # ✅ 成功导出 2D CSV
se                   # ✅ 成功导出 1D CSV
```

### 编译

```bash
make clean
make eq
```

### 相关文件

- `eqgout.f` - 图形输出主程序（已修改）
- `eqgs2d_export.f` - 2D CSV 导出子程序（未修改）
- `eqgs1d_export.f` - 1D CSV 导出子程序（未修改）
- `in/eq.CFEDR_save.in` - 输入文件（已更新）
- `GFILE_USAGE.md` - 使用文档（已更新）

---

## 其他修改（之前完成）

### 2026-01-14: 1D 和 2D CSV 导出功能

**新增文件**:
- `eqgs1d_export.f` - 导出 23 个 1D 剖面参数到 CSV
- `eqgs2d_export.f` - 导出 10 个 2D 场数据到 CSV
- `read_eqgs1d_csv.py` - Python 1D 数据读取和可视化
- `read_eqgs2d_csv.py` - Python 2D 数据读取和可视化
- `EQGS1D_EXPORT_README.md` - 1D 导出完整文档
- `EQGS2D_EXPORT_README.md` - 2D 导出完整文档
- `DATA_EXPORT_SUMMARY.md` - 数据导出总览

**修改文件**:
- `eqgout.f` - 添加 `SE` (第 68-69 行) 和 `CE` 命令 (第 41-42 行)
- `Makefile` - 添加 `eqgs1d_export.f` 和 `eqgs2d_export.f`

**新增命令**:
- `se` - 导出 1D 剖面到 23 个 CSV 文件
- `ce` - 导出 2D 数据到 10 个 CSV 文件

---

## 总结

通过这次修改，EQ 程序现在可以：

1. ✅ 从 EQDSK 格式文件加载数据
2. ✅ 导出完整的 2D 平衡数据到 CSV（10 个文件）
3. ✅ 导出完整的 1D 剖面数据到 CSV（23 个文件）
4. ✅ 保存到新文件，不覆盖原始输入文件
5. ✅ 使用 Python 工具进行数据分析和可视化

这为从 EFIT 或其他 EQDSK 格式文件进行数据分析和后处理提供了完整的工作流程。
