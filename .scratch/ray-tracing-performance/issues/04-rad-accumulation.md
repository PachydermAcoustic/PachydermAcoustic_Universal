# 04: RAD 并发累加修复（D5／审计 #8）

**What to build:** 修复单频带（BandCount=1）OctaveRay 命中时镜面、散射、散射后镜面三类 Reflection_Analysis_Data 的无锁读改写：采用与现有宽带/多分量路径一致的每频带累加同步；累加前校验目标样本在有效范围内，窗口外样本不写入，保留现有有效样本索引规则。直方图累加与 RAD 累加分别完成，不嵌套持有两者的锁。保持主接收器直接累计，不引入线程接收器副本或末尾合并方案。

**Blocked by:** 01（断言工程承载 T5 断言）

**Status:** done

- [x] 单频带三类 RAD 累加与宽带/多分量路径采用一致的每频带同步
- [x] 累加前校验有效样本范围；窗口外样本不写入；保留现有有效样本索引规则
- [x] 不嵌套持有直方图锁与 RAD 锁；检查相关写入及读取边界，无锁顺序反转
- [x] T5 断言通过：可重复命中同时累加三类别，多工作者结果与单线程参考一致（受控精确输入严格断言，不用浮点容差掩盖丢更新）；覆盖单频带与第三倍频程多分量路径、样本 0、末样本、窗口外输入
- [x] 固定射线数与最小收敛模式行为不回归；详细收敛模式因纠正丢更新导致的终止时机/最终 IR 变化已明确说明，不回退修复迁就旧竞争结果

## Comments

**实现（`Classes_Receivers.cs` `Spherical_Receiver.CheckRay` 单频带分支）：** 三类写入合并为与多分量分支及 `Pach_Map_Receiver` 完全一致的形态——`kindOctave = R.IsSpecular ? (R.WasScattered ? 2 : 0) : 1`，先校验 `sample >= 0 && sample < RAD[kindOctave][R.Octave].Length`，再 `lock (RAD[kindOctave][R.Octave])` 内累加。原代码只检查 `sample < Length`：负样本索引会直接 `IndexOutOfRangeException`（已由 T5 "before window" 场景在修复前复现）。`Recs.Add(...)` 在取 RAD 锁之前完成返回，直方图锁与 RAD 锁从不嵌套持有（与既有两条路径一致，代码形态自明，未加注释）。未引入线程接收器副本或末尾合并。

**锁顺序核查：** 全仓库 RAD 写入点（`CheckBroadbandRay`、`CheckRay` 两分支、另两个接收器类）均为单层 `lock` 且先完成 `Recs.Add`；`Receiver_Bank` 合并写入（第 478 行）位于任务完成后的合并流程，不在追踪工作者内。无锁顺序反转路径。

**T5（`tests/AnalysisCancellation/RadAccumulationSuite.cs`）：** 单频带（8 band 全频带 × 三类别）与第三倍频程（BandStart=21，Fractions={0.5,0.25,0.25} × 三类别）各覆盖样本 0、中部样本、末样本（999）、窗口外（t_sum 超窗与负时间）。夹具使每次命中贡献为精确可表示的二进有理数（零衰减谱 + 恰好为 1 的 tsphere，Intensity=π ⇒ 单频带贡献恰为 1.0，第三倍频程每 band 恰为 0.5/0.25/0.25），8 工作者 × 5000 次命中后按无浮点容差的严格相等断言：多工作者 == 单线程参考 == 闭式总数（如 40000×1.0、20000×0.5）。修复前实测丢失更新（39982 vs 40000）即被捕获；窗口外输入断言不写入且不抛异常。net8.0 与 net48 均通过，重复运行稳定。

**模式行为说明（D5）：** 固定射线数与最小收敛（负值）模式行为不回归——单线程算术表达式与求值顺序未变，`RayTraceDeterminismSuite` 修复前后均通过；本修复未触碰直方图、滤波与收敛循环代码。详细收敛模式（Raycount = 0）的 `SpecMin_Convergence_Check` 读取三类 RAD：修复前单频带命中存在丢更新，修复后累计值完整正确，其终止轮次/最低射线数与最终 IR 可能与修复前有竞争的旧结果不同。这是纠正丢更新的合法结果变化，不回退修复迁就旧竞争输出（验收以串行正确参考为准）。
