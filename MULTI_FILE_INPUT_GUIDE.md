# 多文件输入配置指南 (Multi-File Input Guide)

## 文件输入配置 (7个文件)

### 配置文件格式

在 `postppp_25065_gamg.conf` 中按如下格式配置：

```ini
# 观测文件 (2天) - Observation files (2 days)
infile0 = /path/to/obs/gamg0650.25h  # 第1天观测 (DOY 065)
infile1 = /path/to/obs/gamg0660.25h  # 第2天观测 (DOY 066)

# 导航文件 (3天) - Navigation files (3 days)
infile2 = /path/to/nav/brd40640.25p  # 前一天导航 (DOY 064) - 可选但推荐
infile3 = /path/to/nav/brd40650.25p  # 第1天导航 (DOY 065)
infile4 = /path/to/nav/brd40660.25p  # 第2天导航 (DOY 066)

# B2b改正文件 (2天) - B2b correction files (2 days)
infile5 = /path/to/B2b/2025065.B2b   # 第1天B2b (DOY 065)
infile6 = /path/to/B2b/2025066.B2b   # 第2天B2b (DOY 066)
```

## 工作原理 (How It Works)

### 1. 文件处理顺序
- postppp **顺序读取**所有观测文件 (infile0 → infile1)
- 所有历元按**时间顺序**处理，即使来自不同文件
- 导航和B2b文件根据历元时间自动匹配

### 2. 跨天差分模糊度解算流程

```
时间线:
│
├─ DOY 065 ────────────────┐
│  (第1天)                  │
│  - 读取 gamg0650.25h     │
│  - 计算浮点解            │
│  - 保存模糊度到内存      │ check_day_change()
│                          │ 检测到 DOY 变化!
├─ DOY 066 ────────────────┤
│  (第2天)                  │
│  - 读取 gamg0660.25h     │
│  - 计算浮点解            │
│  - 形成跨天差分方程:     │
│    ΔN = N(day2) - N(day1)│
│  - 固定整数模糊度        │
│  - 输出固定解            │
└──────────────────────────┘
```

### 3. 代码中的关键机制

#### `check_day_change()` 函数
```c
static int check_day_change(gtime_t time) {
    int doy = (int)time2doy(time);  // 获取年积日 (day of year)

    if (doy != current_doy) {
        current_doy = doy;
        return 1;  // 天数改变
    }
    return 0;  // 同一天
}
```

- 自动检测 DOY 变化
- **不依赖**文件名或文件边界
- 只依赖观测历元的时间戳

#### 处理流程
```c
pppamb() {
    if (check_day_change(obs[0].time)) {
        // 第1天结束 → 保存当前模糊度作为"前一天参考"
        save_amb_for_next_day(rtk, obs, n);
        return 0;  // 第1天无法固定
    }

    // 第2天及以后 → 形成跨天差分
    ndd = form_dd_equations(rtk, obs, n, ...);

    // 固定整数模糊度
    if (resolve_dd_ambiguity(...)) {
        apply_fixed_ambiguity(...);
        update_clock_fixed(...);
        return 1;  // 固定解成功
    }
}
```

## 为什么当前代码无需修改?

### ✅ 已支持的特性

1. **多文件自动合并**: postppp内部会顺序读取所有infile0-6
2. **时间连续性**: 历元按时间顺序处理，与文件边界无关
3. **DOY检测**: `check_day_change()` 基于时间戳，不是文件索引
4. **状态持久化**: `prev_day_amb[]` 静态数组在文件切换时保持不变

### 📋 数据流示例

```
文件读取:
  infile0 (065.25h) → epoch t1, t2, ... t_end_of_day1
                      ↓ (DOY = 065)
                      save_amb_for_next_day() 保存模糊度

  infile1 (066.25h) → epoch t_start_of_day2, ...
                      ↓ (DOY = 066, DOY改变!)
                      check_day_change() = 1
                      form_dd_equations() 使用保存的模糊度
                      ↓
                      固定解输出
```

## 注意事项

### ⚠️ 重要提醒

1. **文件时间连续性**: 确保 infile1 的起始时间紧接 infile0 的结束时间
2. **观测数据完整性**: 跨天边界附近需要有共同可见卫星
3. **导航文件覆盖**: 建议导航文件多覆盖1天（前一天），确保历元0点附近有有效星历
4. **B2b文件同步**: B2b改正文件必须与观测文件时间段匹配

### 文件命名示例

```
观测文件 (RINEX 3.x):
- gamg0650.25h  → 2025年第65天 (3月6日)
- gamg0660.25h  → 2025年第66天 (3月7日)

导航文件:
- brd40640.25p  → 2025年第64天广播星历
- brd40650.25p  → 2025年第65天广播星历
- brd40660.25p  → 2025年第66天广播星历

B2b文件:
- 2025065.B2b   → 2025年第65天B2b改正
- 2025066.B2b   → 2025年第66天B2b改正
```

## 验证配置

### 运行测试
```bash
cd /home/user/RTKLIB-B2b/example/postppp
../../bin/postppp -k conf/postppp_25065_gamg.conf
```

### 检查输出
```bash
# 查看定位结果
cat out/2025065_GAMG.pos

# 检查是否有固定解 (Q=1)
grep "^2025" out/2025065_GAMG.pos | awk '{print $7}' | sort | uniq -c
#  Q=5: 浮点解
#  Q=1: 固定解 (期望在第2天出现)
```

## 常见问题

**Q: 为什么第1天全是浮点解?**
A: 跨天差分AR需要2天数据作差分，第1天只能计算浮点解并保存模糊度。

**Q: 如果某些时刻两天没有共同卫星怎么办?**
A: 代码会跳过这些卫星，只要有≥4颗共同卫星就能进行AR。

**Q: 是否可以只输入2天导航文件?**
A: 可以，但推荐3天确保午夜0点附近星历可用。

**Q: infile顺序能否调整?**
A: 建议保持 obs→nav→B2b 的顺序，postppp会正确识别文件类型。

## 总结

✅ **当前 ppp_ar.c 无需修改**
✅ **只需正确配置 7 个 infile**
✅ **代码自动处理多文件和跨天检测**

代码已针对多文件输入优化，核心设计确保无论观测数据来自单文件还是多文件，AR算法都能正确工作。
