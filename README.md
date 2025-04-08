# RTKLIB-B2b Application
An open-source software toolkit for decoding and positioning using BeiDou-3 (BDS-3) PPP-B2b service signals, developed based on RTKLIB. Built with C/C++ and CMake for high efficiency and cross-platform compatibility.

## 🚀 Core Features

### PPP-B2b Decoding
- ​**Hardware Support**
  - ComNav(SinoGNSS) receivers (e.g. K803W)
  - Unicore receivers (e.g. UM980)
  
- ​**Data Formats**
  - Proprietary binary PPP-B2b correction streams
  - GPS LNAV/BDS CNAV1 broadcast ephemerides
  
- ​**Operation Modes**
  - Real-time decoding
  - Post-processing decoding
  
- ​**Output Format**
  - Extended ASCII text compatible with BNC (BKG Ntrip Client)

### PPP Positioning
- ​**Positioning Modes**
  - Real-time positioning
  - Post-processing positioning
  
- ​**Satellite Systems**
  - BDS-only mode
  - GPS+BDS combined mode

## 🌍 Cross-Platform Support
| Platform       | Real-time | Embedded Platforms  | Recommended Usage            |
|----------------|-----------|---------------------|------------------------------|
| ​**Linux**      | ✔️         | Rk3566, Imx6ull     | Realtime and Post-processing |
| ​**Windows**    | ✔️         | -                   | Post-processing only         |

> 💡 ​**Recommendation**: For best real-time performance, use Unix-like systems (Linux) especially for console-based operations.

## 📊 GNSS Datasets for RTKLIB-B2b

### Complementary Datasets for PPP Research & Validation
This collection supports two key application scenarios: ​**real-time playback simulation** and ​**post-processing analysis** of PPP-B2b corrections. Both datasets are archived on Zenodo:

#### a. Real-time Playback Dataset
**Purpose**: Simulate real-time streams using authentic receiver outputs  
**Equipment**:  
- Receivers: SinoGNSS K803W & Unicore UM980  
- Antenna: CNT AT360  
**Location**: APM Building 1 Rooftop (Wuhan, China)  
**Coverage**:  
- Time: DOY 080-082, 2025 (Mar 21-23, 2025)  
- Data: Raw binary logs (SINO/UNICORE formats)  

🔗 [Download on Zenodo](https://doi.org/10.5281/zenodo.15105929)

#### b. Post-processing Analysis Dataset
**Components**:  
- 1 Hz MGEX observations (gamg, ulab, wuh2)  
- DLR (German Aerospace Center) broadcast ephemeris  
- Unicore PPP-B2b corrections  
**Coverage**:  
- Observations: DOY 065-067, 2025 (Mar 6-8)  
- Ephemeris/Corrections: DOY 065-068, 2025 (Mar 6-9)  

🔗 [Download on Zenodo](https://doi.org/10.5281/zenodo.15110885)

**Usage Recommendations**:  
1. Dataset 1 for testing real-time module via `--playback` mode  
2. Dataset 2 for post-processing validation with precise reference  

## ❓ Support & Feedback
🚧 ​**Beta Testing Notice**: This software is currently in active development phase.
ℹ️ ​**Detailed are description in the [User Manual](manual/RTKLIB-B2b_UserManual.pdf).**
ℹ️ ​**Developers in China:https://gitee.com/liu-chunbobo/RTKLIB-B2b.git**

## 📜 License

[![GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


RTKLIB-B2b - GNSS PPP-B2b decoding and positioning toolkit
Copyright (C) 2024 [Chunbo Liu/APM,CAS]

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.