# 3_1 Knot Arm Type 功能实现计划

## 摘要 (Summary)
本项目旨在为 `rust_knot` 增加一个额外的功能：判断 3_1 纽结（Trefoil knot）是 "two-arm" 还是 "three-arm" 状态。该功能将在获取 3_1 纽结核心（Knot core，即一条开链）后运行。算法通过 KMT 简化获取最关键的 3 个交叉点，按在链上的首次出现顺序编号为 1、2、3。每个交叉点在链上对应 2 个位置，通过比较 1 号和 3 号交叉点之间的 4 种连接路径（按原始链上的原子/点数量计算物理距离），选取最短的 2 条路径，如果 2 号交叉点都落在这两条路径内，则为 "two-arm"，否则为 "three-arm"。

## 当前状态分析 (Current State Analysis)
- 目前 `rust_knot` 能够通过 `get_knottype` 识别各种纽结类型，并通过 `find_knot_core` 定位纽结核心（一段开链）。
- 现有的 `kmt_open_chain` 会在简化过程中直接删除多余的点，这会导致丢失点在原始链上的索引（Index），从而无法计算准确的物理距离（以原始点数为度量）。
- `knottype.rs` 中的交叉点检测基于 2D 投影（`cal_intersection`），能够找到所有的交叉点，但并没有专门针对 3_1 纽结提取 3 个关键交叉点并进行拓扑距离分析的功能。

## 拟议更改 (Proposed Changes)

### 1. `src/kmt.rs`
**修改内容与原因：**
- 新增 `pub fn kmt_open_chain_with_indices(points: &mut Vec<(Point3, usize)>)`。
- **原因：** 需要在 KMT 简化过程中保留原始链的 index，以便后续准确计算路径的物理长度（即跨越的原始点数）。
- **实现方式：** 逻辑与现有的 `kmt_open_chain` 完全一致，只是操作的数据结构从 `Point3` 变为 `(Point3, usize)`。

### 2. `src/arm_type.rs` (新建模块) 或 `src/knottype.rs`
**修改内容与原因：**
- 新增 `pub fn get_31knot_arm_type(points: &[Point3], table: &AlexanderTable, config: &KnotConfig) -> Result<String>`。
- **原因：** 封装新的 3_1 纽结形态判断逻辑。
- **实现步骤：**
  1. 调用 `get_knottype` 验证输入片段确实是 `3_1` 纽结（如果不是则返回 Error）。
  2. 将输入的 `points` 映射为带索引的 `(Point3, usize)` 数组。
  3. 调用 `kmt_open_chain_with_indices` 进行简化。
  4. （可选）使用 `hull_ends` 对两端进行延伸，防止端点干扰交叉检测（延伸点沿用两端的 index）。
  5. 提取坐标，调用 `find_max_span` 确定最佳 2D 投影平面。
  6. 遍历简化后的线段，使用 `cal_intersection` 计算所有有效的交叉点（期望恰好得到 3 个交叉对）。若不是 3 个则返回 Error。
  7. 收集每个交叉点对应的两条线段的原始 index，形成 3 个交叉点、共 6 个位置（Location）。
  8. 对 6 个位置排序，按首次出现顺序将交叉点命名为 1、2、3。
  9. 计算 1 号和 3 号交叉点之间的 4 种连接路径（绝对差值 $|A - B|$），并选出最短的 2 条。
  10. 检查 2 号交叉点（具有 2 个位置）是否至少有一个位置落在这 2 条最短路径的区间内。
  11. 如果 2 号交叉点在两条最短路径内都出现过，则返回 `"two-arm"`，否则返回 `"three-arm"`。

### 3. `src/lib.rs`
**修改内容：**
- 导出新编写的 `get_31knot_arm_type` 函数，使其作为公共 API 可用。

### 4. `tests/integration.rs`
**修改内容：**
- 增加一个新的集成测试：读取 `L300_knot3_1_ring.xyz`，获取其 `3_1` 纽结核心，然后将核心传入 `get_31knot_arm_type` 进行测试，验证其能正确返回 "two-arm" 或 "three-arm" 而不报错。

## 假设与决策 (Assumptions & Decisions)
1. **输入限制**：该功能假设用户传入的是一段 **Open chain（纽结核心）**。即使原始数据是环，用户也应当先调用 `find_knot_core` 获取核心的连续点序列。
2. **恰好 3 个交叉点**：依据用户的确认，经过 KMT 简化后的 3_1 纽结 2D 投影应恰好包含 3 个交叉点。如果由于特殊构象导致简化后仍有多余的平凡交叉点，程序将安全地返回 Error，而不会产生误判。
3. **最短路径度量**：依据用户的确认，最短路径以**原始链上的原子数量**（即原始 Index 的绝对差值）为度量标准，这具有明确的物理意义。

## 验证步骤 (Verification Steps)
1. 编译项目 (`cargo build`) 确保无语法错误。
2. 运行现有的所有测试 (`cargo test`) 确保没有破坏已有功能。
3. 运行新增的集成测试，验证 `get_31knot_arm_type` 能够针对给定的 3_1 核心正确返回状态结果。