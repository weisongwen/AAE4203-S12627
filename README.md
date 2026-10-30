# AAE4203 - Guidance and Navigation

## Course Repository
**The Hong Kong Polytechnic University**  
**Department of Aeronautical and Aviation Engineering**  
**Academic Year 2026-2027 Semester 1**

---

## 📚 Course Information

- **Course Code**: AAE4203
- **Course Title**: Guidance and Navigation
- **Credit Units**: 3
- **Lecture Time**: Friday 3:30pm - 6:30pm at CD302
- **Instructor**: Prof. Weisong WEN ([welson.wen@polyu.edu.hk](mailto:welson.wen@polyu.edu.hk))
- **Teaching Assistants**: 
  - Mr. Jiahao Hu
  - Mr. Yingming Chen
  - Mr. Xupeng Zhai

---

## 🎯 Subject Intended Learning Outcomes (SILOs)

Upon completion of this subject, students will be able to:

1. Understand and explain the working principles of navigation and guidance systems for unmanned autonomous systems (UAS)
2. Competently apply the fundamental mathematical concepts of UAS navigation
3. Critically evaluate the characteristics, purposes, and design procedures of UAS navigation and guidance systems
4. Identify the technological and design trends of future UAS navigation

---

## 📊 Assessment Scheme

| Assessment Component | Weighting |
|---------------------|-----------|
| Assignment 1        | 7.5%      |
| Assignment 2        | 7.5%      |
| Mid-term Test       | 15%       |
| Group Project       | 20%       |
| Final Examination   | 50%       |
| **Total**           | **100%**  |

---

## 📁 Repository Structure

```
AAE4203-S12627/
├── Lab1/                                        # Lab 1 (Week 9): GNSS Positioning and Accuracy Challenge
│   ├── AAE4203_26_27_S1_Lab_Tutorial.pdf        # Lab tutorial 2026/27 (theory, workflow, challenge rules)
│   ├── AAE4203_Week9_GNSS_Lab_Slides.pdf        # Week 9 lab session slides
│   ├── PolyU_AAE4203_Lab_report_template.docx   # Lab report template
│   ├── WLSE_spp.py                              # Runnable weighted least-squares SPP (example and benchmark)
│   ├── pyubx2_csv_converter_gui.py              # UBX -> CSV converter (GUI)
│   ├── rinex2csv.m                              # RINEX -> CSV MATLAB function
│   ├── spp.conf                                 # RTKLIB SPP configuration
│   ├── example/                                 # Example dataset (Block X, 2025): UBX, RINEX, CSV inputs, receiver reference
│   ├── benchmark/                               # UrbanNav benchmark: OBS/NAV and GPS CSV inputs (no ground truth)
│   ├── tools/                                   # check_submission.py (CSV checker), verify_package.py (self-test)
│   ├── docs/                                    # Running instructions, verification and data sources
│   ├── legacy/                                  # Original 2025 and course scripts, for comparison
│   ├── matRTKLIB/                               # MATLAB interface for RTKLIB
│   └── figures/                                 # Example figures & workflow
├── Lab2/                                        # Lab 2: Inertial Navigation and PDR
│   ├── AAE4203 Lab2.pptx                        # Lab tutorial slides
│   ├── MobilePDR.py                             # PDR algorithm implementation
│   └── data/                                    # Sample datasets (RoundTrip, Square)
├── Week01/                                      # Week 1 slides & materials
├── Week02/                                      # Week 2 slides & materials (later weeks added progressively)
├── UIDemo.html                                  # Smartphone IMU data collection demo
└── AAE4203_Teaching_Plan_2026-27_Semester_1.pdf # Official teaching plan
```

---

## 🗂️ Weekly Materials

Weekly lecture slides and extra materials are released as the semester progresses:

- [Week 1](Week01/) — Mathematical Foundations for Guidance and Navigation - Part 1 (Sept 4)
- [Week 2](Week02/) — Mathematical Foundations for Guidance and Navigation - Part 2 (Sept 11)
- [Week 9](Lab1/) — Laboratory: GNSS Data Collection, Single Point Positioning and an Accuracy Challenge (Oct 30)

> Materials for upcoming weeks will be posted after each lecture.

---

## 📅 Teaching Schedule

| Week | Date    | Topics |
|------|---------|--------|
| 1    | Sept 4  | Mathematical Foundations for Guidance and Navigation - Part 1 |
| 2    | Sept 11 | Mathematical Foundations for Guidance and Navigation - Part 2 |
| 3    | Sept 18 | Reading Week - No Class |
| 4    | Sept 25 | GNSS Fundamentals - Part 1 |
| 5    | Oct 2   | GNSS Single Point Positioning (SPP) - Part 1 |
| 6    | Oct 9   | GNSS Single Point Positioning (SPP) - Part 2 |
| 7    | Oct 16  | GNSS Velocity Estimation |
| 8    | Oct 23  | **Mid-term Test** |
| 9    | Oct 30  | **Laboratory: GNSS Positioning and Velocity Estimation** |
| 10   | Nov 6   | GNSS RTK and High-Precision GNSS Positioning |
| 11   | Nov 13  | Group Project Presentations |
| 12   | Nov 20  | Group Project Presentations |
| 13   | Nov 27  | Tutorial and Revision |

> **Note**: Week 3 is designated as a reading week without an in-class meeting. Week 13 is reserved for tutorial and revision to consolidate the course material before the final examination.

---

## 📖 Required Textbooks

1. **Farrell, J.A. (2008)**. *Aided Navigation: GPS with High Rate Sensors*. McGraw-Hill.
2. **Groves, P.D. (2013)**. *Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems*. Artech House.

---

## 📜 Course Policies

- **Attendance**: Regular attendance is expected for all classes.
- **Late Submission**: A 10% penalty per day applies to late assignments.
- **Academic Integrity**: All submitted work must be original and comply with PolyU requirements.
- Students must attend the final examination in order to attain a passing grade in the overall result.

---

## 🔗 Resources

- **Lab 1 Guide**: [GNSS Single Point Positioning (SPP) Analysis Pipeline](Lab1/Readme.md) · [Tutorial 2026/27 (PDF)](Lab1/AAE4203_26_27_S1_Lab_Tutorial.pdf)
- **Lab 1 Leaderboard Platform**: [PolyU-TASLAB/AAE4203_ASN_2026](https://github.com/PolyU-TASLAB/AAE4203_ASN_2026)
- **Lab 2 Guide**: [Inertial Navigation and PDR](Lab2/readme.md)
- **Previous Year (2025-26 Semester 1)**: [AAE4203-S12526](https://github.com/weisongwen/AAE4203-S12526)
