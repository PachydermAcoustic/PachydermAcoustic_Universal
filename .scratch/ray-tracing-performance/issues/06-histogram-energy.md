# 06: 方向直方图单一 Energy 分配（D3／审计 #2）

**What to build:** 两个方向直方图构造函数删除子类中重复的 Energy 分配，只保留基类已完成的初始化——构造后只有一份最终 Energy 数据。基类初始化、Pressure、正负三轴方向数据及三类 Reflection_Analysis_Data 的结构和零值语义保持不变。兼容 8/24 band、样本数及毫秒时长构造；结果读写格式不变。

**Blocked by:** 01（断言工程承载 T3 断言）

**Status:** done

- [x] 两种构造方式均只保留一份最终 Energy 数据（基类初始化那份），子类不再替换
- [x] Pressure、正负三轴方向、三类 RAD 结构与零值语义不变；8/24 band、样本数、毫秒时长构造兼容；结果读写格式不变
- [x] T3 断言通过：两种构造方式 × 两种频带数下 Energy/Pressure/方向结构一致；零值及末样本正确
- [x] 分配对照（池已暖、输入规模固定）证明不再付出第二份 Energy 成本；不写死精确字节承诺、不声称稳定存活内存等量下降

## Comments

Done on branch ticket/06. Classes_Receivers.cs Directional_Histogram: both constructors no longer re-allocate Energy; the base Histogram allocation is the single surviving set (length equivalence verified: int ctor uses the same SampleCT, ms ctor base SampleCT == subclass (int)(SampleRate*CO_Time) by FP commutativity; no discrepancy found). All other initialisation (SampleRate, CO_Time, SampleCT, Dir_Rec_Pos/Neg) untouched; RAD untouched. T3 suite: HistogramEnergySuite (one Register line in Program.cs) covers both ctors x 8/24 bands via public surface (Spherical_Receiver double-time ctor incl. bank path and fractional-cutoff truncation 333Hz/1500ms->499, PachMapReceiver.Map_Receiver int-sample ctor, Receiver_Bank.Write_Data/Read_Data round-trip), asserts live-accessor identity, base-semantics dimensions, all-zero content, writable last sample, RAD zero semantics. Allocation proof: differential GC.GetAllocatedBytesForCurrentThread measurement on warm ctors, fixed N/bands, S=2000 vs 4000 -> exactly 11 band-sample arrays for Spherical_Receiver (1E+1P+6dir+3RAD; was 11.999 before fix) and 8 for Map_Receiver (no RAD); no byte counts hardcoded. Harness csproj: net48 flavor now builds x64 (bundled libfftw3-3.dll is x64-only; suite round-trip runs the FFTW-backed filter creation). Red (11.999/FAIL) -> green verified on net8.0 (exit 0) and net48 (exit 0); library builds 0 errors net7.0/net48.
