# 03: 滤波完成与失败契约（D4／审计 #3）

**What to build:** Receiver_Bank 三个滤波入口（批量 void 入口与两个返回滤波数组的重载）具备明确的同步完成与失败契约：正常返回即全部所需滤波已安装；失败向调用方抛出原始异常，不转为后台未观察异常、不发布半成品结果；CreatingFilter 在正常与异常退出后均复位。移除批量入口中无实际异步工作的修饰与死任务代码，保持串行顺序与公共签名。逐一核查同步反序列化、结果重建、声源编辑、Hybrid 和 Mapping 调用场景：同步调用直接观察失败；后台 UI 调用等待其 Task 并经现有错误展示路径报告、清理进度状态。

**Blocked by:** 01（断言工程承载 T4 断言）

**Status:** done

- [x] 批量 void 入口保留公共签名与串行顺序；正常返回意味着全部滤波已安装；移除多余 async 修饰与死任务代码
- [x] 失败向调用方抛出原始失败；正常/异常退出均复位 CreatingFilter；失败不发布新完成结果
- [x] 两个返回数组的重载保留签名与返回值语义；运算结束后返回、传播失败、正常/异常退出后复位 CreatingFilter
- [x] 同步调用方直接观察失败；后台 UI 调用等待其 Task、经现有错误路径报告并清理进度窗口与处理标志（调用方核查已完成，见 Comments；Rhino 侧错误展示/清理由 ticket 08 实施）
- [x] T4 断言通过：三个入口的可控失败均可被调用者观察；批量后台调用的失败可由 await 观察；各入口失败后 CreatingFilter 复位
- [x] 进度反馈沿用现有 UI 调度方式；FFT 槽位策略与滤波算法不变；不做接收器并行

## Comments

**Change** (`Classes_Receivers.cs`): batch entry is now a plain `virtual void` (same signature, serial order kept); dead `Task[] t` / commented parallel scaffolding / `await Task.WhenAll` remnants and the now-unused `using System.Threading.Tasks` removed. All three entries wrap their work in `try/finally` so `CreatingFilter` resets on normal and failing exits. Pre-existing latent bug fixed by the same reset: both array-returning overloads used to leave `CreatingFilter` stuck `true`, which would hang every later `HasFilter(int)` spin-wait on that bank (reproduced red: `HasFilter(0)` did not return within 10 s).

**T4 suite** (`tests/AnalysisCancellation/FilterContractSuite.cs`): failure seam is a `Spherical_Receiver` subclass whose three `public virtual Create_Filter` overloads throw (verified overridable), installed through the public `Rec_List` field of a bank built by the public constructor. `CreatingFilter` is private, so its reset is asserted via the public `HasFilter(int)` spin-wait behind a 10 s timeout probe. Covers: batch normal return ⇒ every receiver `HasFilter` true; overload normal return ⇒ non-empty arrays; batch sync failure observed by caller; batch failure in `await Task.Run(...)` observed via await; both overloads' failures propagate synchronously; flag reset after every path. Red→green verified on net8.0; net48 exe passes after `<PlatformTarget>x64</PlatformTarget>` was added to the harness csproj (AnyCPU net48 exes run 32-bit and x64 libfftw3-3.dll then fails to load; ticket 01's suite never exercised FFTW).

**Caller audit** (read-only; Rhino files untouched — ticket 08 owns them):

Universal:
- `Receiver_Bank.Read_Data` (`Classes_Receivers.cs:362`, synchronous deserialization): observes failure directly now (previously an `async void` failure escaped to the thread pool).
- `IR_Construction.Auralization_Filter` (`Classes_Math.cs:2794`, `double[]` overload) and `Aurfilter_Directional` (`Classes_Math.cs:2848`, `double[][]` overload): synchronous; failure propagation unchanged, but the flag no longer sticks (see bug above).
- `Pach_Map_Receiver.cs:984/1048/1112`: call the inner `Spherical_Receiver.Create_Filter()`, not the bank entries — unaffected.

Rhino:
- `Pach_Hybrid_Control.cs:1287` (Hybrid RT stage): background via `await Task.Run(() => RT.GetReceiver.Create_Filter(VB))` — already awaits its Task, so the failure is now observable there; `VB.Close()` and `Receiver[s]` publication currently run only on success → ticket 08.
- `Pach_Hybrid_Control.cs:1737` (`Set_Phase_Regime`, sync UI loop) and `Pach_Mapping_Control.cs:2516` (`Coherent_CheckedChanged`, sync UI loop): synchronous callers — they now observe failures directly. Neither has a catch today, so a previously-swallowed failure now surfaces in the UI handler and skips `VB.Close()`; that is the intended contract, error display/cleanup belongs to ticket 08.
- `Classes_UI Components.cs:181` (source-power edit, sync UI): `if (Rec[i].HasFilter()) Rec[i].Create_Filter(VB)` — same situation as above; `VB.Close()` skipped on throw → ticket 08.
- `Pach_Hybrid_Control.cs:1203/1732`, `Classes_UI Components.cs:189`, `Pach_Mapping_Control.cs:2112`: `Direct_Sound`/`ImageSourceData` entries or the inner receiver method — outside this contract.
