# PPP-B2b Fixed Solution Implementation

## Overview

This implementation adds **fixed solution mode** for PPP-B2b processing to RTKLIB-B2b. The fixed solution uses **day-to-day differencing ambiguity resolution** technique as described in Xi et al., 2021, to overcome the difficulty of fixing PPP-B2b ambiguities.

## Reference

Xi, R., Jiang, W., Meng, X., Chen, H., & Chen, Q. (2021). Pass-by-Pass Ambiguity Resolution in Single GPS Receiver PPP Using Observations for Two Sequential Days. *Remote Sensing*, 13(18), 3631. https://doi.org/10.3390/rs13183631

## Key Features

### 1. Day-to-Day Differencing Ambiguity Resolution

The traditional PPP ambiguity resolution is challenging for PPP-B2b due to:
- Limited correction accuracy compared to IGS products
- Residual satellite orbit and clock errors
- Uncalibrated phase biases

The day-to-day differencing technique addresses these challenges by:
- Using observations from two consecutive days
- Forming differences between ambiguities estimated on different days
- Leveraging the fact that satellite geometry repeats approximately every sidereal day
- Canceling common systematic errors through differencing

### 2. Implementation Details

#### Core Algorithm (`src/ppp_ar.c`)

The implementation consists of several key components:

1. **Previous Day Ambiguity Storage**
   - Stores ambiguity estimates from the previous day
   - Includes ambiguity values, standard deviations, and lock counts
   - Automatically resets when day changes

2. **Day-to-Day Difference Formation**
   - Matches satellites between consecutive days
   - Computes ambiguity differences: `dd_amb = amb_current - amb_previous`
   - Validates data quality using:
     - Minimum lock time (30 epochs by default)
     - Maximum ambiguity variance threshold (0.25 cycles²)
     - Minimum satellite count (4 satellites)

3. **Integer Ambiguity Resolution**
   - Rounds day-to-day differences to nearest integer
   - Performs ratio test for validation (threshold: 3.0)
   - Note: Current implementation uses simple rounding; LAMBDA algorithm can be added for improved performance

4. **Fixed Solution Update**
   - Applies fixed ambiguities to state vector
   - Updates ambiguity variances to small values (1e-6)
   - Marks satellites as "fixed" in solution status

5. **Clock Bias Correction**
   - Recomputes receiver clock using fixed ambiguities
   - Uses code-phase consistency check
   - Applies elevation-dependent weighting

#### Modified Files

1. **`src/ppp_ar.c`** (Complete rewrite)
   - Implements day-to-day differencing algorithm
   - Handles ambiguity storage and matching
   - Performs integer resolution and validation
   - Updates clock bias for fixed solution

2. **`src/ppp.c`** (Lines 1355-1373 modified)
   - Updates `update_stat()` function
   - Uses fixed solution clock (`rtk->xa`) when solution status is `SOLQ_FIX`
   - Maintains float solution clock (`rtk->x`) for float solutions

3. **`example/postppp/conf/postppp_25065_gamg.conf`** (Updated)
   - Adds AR configuration parameters
   - Enables PPP-AR mode (`modear=4`)

## Configuration

### Enable Fixed Solution Mode

Add the following parameters to your configuration file:

```ini
# Ambiguity Resolution (AR) options for PPP-B2b fixed solution
prcopt.modear     = 4    # AR mode (0:off, 4:ppp-ar)
prcopt.minlock    = 30   # minimum lock count to fix ambiguity (epochs)
prcopt.minfixsats = 4    # minimum satellites for AR
prcopt.thresar[0] = 3.0  # AR ratio test threshold
```

### AR Mode Options

- `modear = 0`: AR disabled (float solution only)
- `modear = 1`: Continuous AR (for RTK)
- `modear = 2`: Instantaneous AR (for RTK)
- `modear = 3`: Fix-and-hold (for RTK)
- `modear = 4`: PPP-AR with day-to-day differencing (for PPP-B2b)

## Usage

### Processing Two-Day Data

For optimal performance, process observations spanning at least two consecutive days:

```bash
./postppp -c config.conf
```

**Important Notes:**
1. On the first day, only float solutions will be produced
2. Starting from the second day, fixed solutions become available
3. Ensure continuous tracking of satellites across day boundaries
4. Higher elevation satellites typically provide better AR performance

### Output Solution Status

The solution status field in the output indicates:
- `5` : Single point positioning
- `6` : PPP float solution
- `1` : PPP fixed solution

Fixed solutions also include:
- Updated receiver clock bias
- Tightened position uncertainties
- Ratio test statistic

## Algorithm Parameters

The following parameters can be tuned in `src/ppp_ar.c`:

```c
#define MIN_ARC_SAT     4       /* minimum satellites for AR */
#define MIN_ARC_TIME    30      /* minimum continuous lock time (epochs) */
#define THRES_RATIO     3.0     /* ratio test threshold */
#define THRES_VAR_AMB   0.25    /* variance threshold for ambiguity (cycles²) */
```

## Output Clock Products

### Float Solution Clock Output

When `modear = 0` or AR fails, the clock output is from the float solution:
- Clock values are from `rtk->x[IC(i,opt)]`
- Reflects code-phase smoothed estimates
- Higher uncertainty (~ns level)

### Fixed Solution Clock Output

When `modear = 4` and AR succeeds, the clock output is from the fixed solution:
- Clock values are from `rtk->xa[IC(i,opt)]`
- Uses fixed ambiguities for clock estimation
- Improved precision (~0.1 ns level)
- Better consistency with carrier phase observations

The fixed clock is computed by:
1. Using fixed ambiguities to compute phase ranges
2. Comparing phase ranges with code ranges
3. Adjusting clock to minimize code-phase inconsistency
4. Applying elevation-dependent weighting

## Validation

To validate the fixed solution:

1. **Check Ratio Test**: Should be > 3.0 for reliable fixes
2. **Monitor Position Consistency**: Fixed solutions should show reduced scatter
3. **Verify Clock Stability**: Fixed clocks should have better short-term stability
4. **Compare with IGS Products**: If available, compare with IGS precise clocks

## Limitations and Future Work

### Current Limitations

1. **Simple Integer Resolution**: Uses rounding instead of LAMBDA
   - May have lower success rate in challenging conditions
   - Can be improved by implementing MLAMBDA algorithm

2. **Single Frequency**: Currently uses L1/B1 only
   - Multi-frequency AR would improve reliability
   - Widelane-narrowlane approach could be implemented

3. **No Partial AR**: Either all or nothing
   - Partial AR could improve solution continuity
   - Subset selection strategies could be added

4. **Day Boundary Handling**: Requires careful data preparation
   - Automatic handling could be improved
   - Cross-day smoothing could be added

### Future Improvements

1. Implement LAMBDA/MLAMBDA for better integer search
2. Add multi-frequency ambiguity resolution
3. Implement partial ambiguity resolution
4. Add quality control and outlier detection
5. Optimize for real-time processing
6. Support for multi-GNSS (GPS+BDS+Galileo)

## Testing

### Test Data Requirements

- At least 2 consecutive days of observation data
- PPP-B2b correction data for both days
- Good satellite visibility (>4 satellites continuously)
- Minimal cycle slips across day boundaries

### Expected Performance

- **Float Solution Precision**: 2-5 cm horizontal, 5-10 cm vertical
- **Fixed Solution Precision**: 1-2 cm horizontal, 2-5 cm vertical
- **Fix Rate**: 30-70% depending on conditions
- **Convergence Time**: Faster convergence after ambiguity fixing

## Troubleshooting

### No Fixed Solutions

- Check if `modear = 4` in configuration
- Ensure at least 2 days of data are processed
- Verify satellite lock times > 30 epochs
- Check ambiguity variances < 0.25 cycles²

### Low Fix Rate

- Increase observation time (longer arcs)
- Improve data quality (reduce cycle slips)
- Lower `THRES_RATIO` threshold (with caution)
- Use higher elevation mask

### Fixed Solution Jumps

- May indicate incorrect ambiguity fixes
- Increase `THRES_RATIO` for stricter validation
- Check for cycle slips
- Verify data quality

## Contact

For questions or issues related to this implementation, please open an issue on the GitHub repository.

## Acknowledgments

This implementation is based on the research by Xi et al. (2021) and extends the original RTKLIB software by T. Takasu.

---

**Version**: 2.0
**Date**: 2025-12-04
**Status**: Initial implementation
