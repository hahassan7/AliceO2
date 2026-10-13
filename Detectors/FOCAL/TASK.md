**Subject:** FOCAL-E pads digitizer – your first task (v1)

---

## 0. Before you start: use Claude as your assistant

O2 is a big framework, and nobody expects you to know it already. Our group has access to Claude (see the email from Nishida-kun a few weeks ago with the link for our group plan). I recommend you use it as an assistant for this task. It can read the O2 code and the papers, explain them to you, and write and fix code together with you.

**Setup:**

1. Join our group plan with the link in Nishida-kun's email, and install the **Claude desktop app** (claude.ai/download).
2. Make a working folder, e.g. `~/FOCAL-E_Pads_Digitizer/`, and put in it:
   - `AliceO2/`: a clone of my branch (link in section 2);
   - `TestbeamFOCAL_OldFramework/`: the old test-beam analysis repository ([https://gitlab.cern.ch/alice-focal-readout/focal_testbeamanalysis]), which contains the calibration code we used for the test beam;
   - `TestbeamPaper.pdf`: the FoCal prototype test-beam paper. Downloaded here: [https://arxiv.org/pdf/2311.07413];
   - `HGCROC3_datasheet_v1.pdf`: the HGCROC datasheet;
   - the two calibration files attached in the email;
   - `TASK.md`.
3. In the Claude app, create a **Project** called `FOCAL-E Pads Digitizer` and give it access to this folder. Then you can work on the code with Claude in that project.

**First request to Claude** (copy and paste):

> Please read TASK.md, then the FOCAL code in AliceO2/Detectors/FOCAL and AliceO2/DataFormats/Detectors/FOCAL (especially SDigitizer, SDigit, LabeledSDigit, PadIndex, PadMapper, FOCALDigitizerSpec, FOCALDigitWriterSpec), the EMCAL simulation code as an example, sections 5.5 and Appendix A of TestbeamPaper.pdf, sections 1.1–1.3 of the HGCROC datasheet, and FOCALEPadCalibration.h in the test-beam repository. Explain to me in simple words how the existing code works and what I have to implement. Then let's implement the task step by step, one class at a time, starting with PadSimParam and HgcrocResponse with a standalone test.

**How to use it well:**

- **You are responsible for the code, not the AI.** Ask it to explain every change until you understand it. Read the diffs before accepting them.
- **Work in small steps:** one class, compile, test, commit; then the next one.
- **Compile and run everything yourself** (O2 build, `o2-sim`, digitizer) and give Claude the error messages; it is very good at fixing them.
- **Check the physics against this email.** If Claude proposes a different equation or number than the ones below, ask why, and ask me before changing it.
- Don't let it push to GitHub or change anything outside our FOCAL directories without asking me.

---

## 1. What we want to achieve

When we simulate FOCAL with `o2-sim`, Geant4 tells us how much energy each particle deposits in each silicon pad. The real detector does not give us energies, though. Each sensor (9 × 8 = 72 pads) is read out by one **HGCROC** chip, which converts the charge of each pad into three integers:

- **ADC**: amplitude of the signal (10 bits, 0–1023). Used for small signals.
- **TOT** (time over threshold, 12 bits, 0–4095): used for large signals, where the ADC saturates.
- **TOA** (time of arrival, 10 bits): the time the signal crossed a threshold.

The **digitizer** is the code that converts the simulated energy into these numbers, so that the simulation looks like real data. **Your task (v1)** is to write this conversion for the FOCAL-E pads, run it on simulated electrons, and compare with the test beam.

**We ignore the following for now:**
- timing and pile-up (one collision = one readout);
- the readout electronics after the HGCROC (ECON-D, zero suppression, raw data);
- bad channels;
- pixels and HCAL.

These come later.

---

## 2. Getting started

**Code.** Start from this branch, which contains everything done so far (summed digitizer, pad indexing, workflow): [https://github.com/hahassan7/AliceO2/tree/Pads-Digitizer]

**Attached** are the two calibration files from the 2023 test beam, which you will need:

- `FOCALEPads_Calibration_20230208.txt`
- `calib_highinjection_pdfextracted.txt`

**Please read this existing code first.** Everything is in `AliceO2/Detectors/FOCAL` and `AliceO2/DataFormats/Detectors/FOCAL`.

1. `simulation/src/SDigitizer.cxx`: the **SDigitizer**. It sums all Geant4 hits of one collision per pad, giving one "summed digit" (SDigit) per pad.
2. `DataFormatsFOCAL/SDigit.h`: the SDigit object. Important getters:
   - `getEnergy()`: deposited energy in **eV**
   - `getIndex()`: unique pad index (`absID`)
   - `isPadDigit()`: we only digitize pads
3. `simulation/include/FOCALSimulation/LabeledSDigit.h`: an SDigit plus its **MC labels** (which simulated particles contributed, and with which fraction of the energy). You must carry these labels through to your digits.
4. `base/include/FOCALBase/PadIndex.h`: how the pads are numbered:
   - `absID = layer*7920 + row*90 + col`
   - `layer` 0–19, where 4 and 9 are pixel layers and never contain pads
   - `row` 0–87, bottom to top
   - `col` 0–89, left to right
   - `getLayerRowColFromIndex(absID)` / `getIndexFromLayerRowCol(...)` convert between `absID` and (layer, row, col)
   - `getPadLocation(absID)` gives the side, module, sensor and position inside the sensor
   - you get it configured for our geometry with `geometry->getPadIndex()`
5. `reconstruction/include/FOCALReconstruction/PadMapper.h`: which HGCROC channel (0–71) reads which pad of a sensor.
6. `workflow/src/FOCALDigitizerSpec.cxx`: the DPL task that runs everything. For each collision it calls the SDigitizer and currently writes the SDigits, the MC labels and one `TriggerRecord` per collision.
7. For comparison, look (just look) at EMCAL: `Detectors/EMCAL/simulation/SimParam.h/.cxx` and `Digitizer.h/.cxx`.

**Background reading:**

- **Test-beam paper:** section 5.5 (pad calibration), Appendix A (calibration of the pads, with Fig. 37), and section 5.8 for Fig. 27.

---

## 3. The code you have to write

You will write three small classes and modify the workflow:

```
SDigits (energy per pad) --PadsDigitizer--> charge per pad --HgcrocResponse--> Digit (ADC, TOA, TOT)
```

### 3.1 `PadSimParam.h/.cxx` (`simulation/`): all parameters in one place

This is an `O2ParamDef` configurable parameter class (see EMCAL `SimParam`), with key `FOCALPadSimParam`. It lets us change any number from the command line, e.g. `--configKeyValues "FOCALPadSimParam.mNoiseADC=2.0"`. **No numbers hard-coded anywhere else.**

| Parameter | Default | Meaning |
|---|---|---|
| `mEnergyPerPair` | 3.65 eV | energy to create one electron–hole pair in silicon |
| `mCCE` | 1.0 | charge collection efficiency |
| `mCalibFileLow` | path to `FOCALEPads_Calibration_20230208.txt` | ADC calibration and low-range TOT |
| `mCalibFileHigh` | path to `calib_highinjection_pdfextracted.txt` | high-range TOT |
| `mADCTOTGap` | 0.100 pC | gap between the end of the ADC dynamic range and the TOT start |
| `mTOTTransitionLow` / `mTOTTransitionHigh` | 0.5 / 2.0 pC | end of the TOT low range / start of the TOT high range |
| `mNoiseADC` | 1.5 ADC | incoherent noise per channel |
| `mCMNADC` | 4.5 ADC | common-mode noise per ASIC half |
| `mTOTSmear` | 2 TOT units | TOT resolution (to be tuned) |
| `mTOAThreshold` | 0.002 pC (2 fC) | TOA fires above this |
| `mReadoutAllSensors` | false | see section 3.3 |
| `mUseBadChannels` | false | keep false; not implemented in v1 |
| `mDebugStream` | false | write the debug file of section 3.5 |

### 3.2 `HgcrocResponse.h/.cxx` (`simulation/`): pad index + charge → ADC, TOT, TOA

This class does the whole electronics response for one pad.

**At initialisation** it reads the two calibration files once into arrays `[18 ASICs][72 channels]`.

Each line of `FOCALEPads_Calibration_20230208.txt` has the columns

`ASIC CHANNEL ADCLowFit ADCHighFit ADCConst ADCSlope TOTLowFit TOTHighFit TOTConst TOTSlope`

and each line of `calib_highinjection_pdfextracted.txt` has the columns

`ASIC CHANNEL TOTSlope_high TOTConst_high`

All charges are in **pC** and all slopes are per pC. Typical values (medians):

| Constant | Typical | Meaning |
|---|---|---|
| ADCConst | 68 | ADC at zero charge = **pedestal** |
| ADCSlope | 4050 ADC/pC | ADC gain |
| ADCHighFit | 0.20 pC | end of the linear ADC range |
| TOTConst / TOTSlope | −16 / 336 per pC | low-range TOT fit (valid ≈ 0.22–0.5 pC) |
| TOTConst_high / TOTSlope_high | 152 / 257 per pC | high-range TOT fit (valid above 2 pC) |

Some channels have nonsense fits (e.g. ADCSlope = 54, or a negative TOTSlope). If a constant is outside a sane window (ADCSlope 2500–4500, TOTSlope 200–550, TOTSlope_high 150–300), use the mean of the same ASIC half instead. For comparison, see `TestbeamFOCAL_OldFramework/focal_testbeamanalysis/include/FOCALEPadCalibration.h`.

**Which calibration line belongs to a pad?** The calibration is numbered by test-beam ASIC and HGCROC channel, while our pads are numbered by `absID`. Do the lookup inside `HgcrocResponse`, with the existing classes:

```cpp
auto loc = padIndex.getPadLocation(absID);               // layer, side, module, sensor, padX (0-8), padY (0-7)
int channel = padMapper.getChannelID(7 - loc.padY, loc.padX);  // HGCROC channel 0-71
// careful: PadMapper's first argument is the pad row counted from the TOP of the sensor, the second is the pad column
int half = (channel < 36) ? 0 : 1;                       // ASIC half (for the common-mode noise)
int tbAsic = loc.layer;                                  // test-beam ASIC: the test beam had only the 18 pad layers
for (int l = 0; l < loc.layer; ++l) {
  if (geometry->isPixelLayer(l)) {
    --tbAsic;                                            // layers 0,1,2,3,5,...,19 -> ASIC 0..17
  }
}
// -> use the constants of (tbAsic, channel)
```

- `PadMapper` lives in the reconstruction library, so add `O2::FOCALReconstruction` to the `PUBLIC_LINK_LIBRARIES` of `Detectors/FOCAL/simulation/CMakeLists.txt`.
- The sensor orientation and numbering are provisional until the engineers give us the final drawings. So never write your own formulas for this; always go through `PadIndex` and `PadMapper`.
- Every sensor of a layer uses the same 72 calibrated channels of "its" test-beam ASIC. This is a preliminary study; the real detector will have its own calibration.

**Precompute once per channel at initialisation** (only 18 × 72 channels, so this costs nothing):

```
G     = ADCSlope                       (ADC per pC)
Q_lin = ADCHighFit                     (end of the linear range)
D     = 1023 - ADCConst - G*Q_lin      (ADC counts left above the linear range; if D <= 0, treat the channel as "crazy" -> half mean)
Q_sat = Q_lin + 2*D/G                  (charge at which the ADC reaches 1023 = end of the ADC dynamic range, ~0.27 pC)
Q_TOT = Q_sat + mADCTOTGap             (TOT starts 100 fC after the end of the ADC dynamic range, ~0.37 pC)
```

**Per collision:** draw the common-mode noise, one value `CMN = Gaus(0, mCMNADC)` per (layer, side, module, sensor, half). Use, for example, a method `startCollision()`.

**Per pad**, given `absID` and charge `Q` in pC:

**Step A: noise.** One Gaussian value per pad, in charge, used for both the ADC and the TOT/TOA thresholds:

```
Qn = Q + Gaus(0, mNoiseADC / G)
```

**Step B: ADC (linear, non-linear, saturated).**

```
if   Qn <= Q_lin:   ADCsig = G*Qn                                          (linear)
elif Qn <  Q_sat:   x = Qn - Q_lin;   ADCsig = G*Q_lin + G*x - G*x*x/(2*(Q_sat - Q_lin))   (non-linear, bends smoothly)
else:               ADCsig = 1023 - ADCConst                               (saturated)

ADC = int(round(ADCConst + ADCsig + CMN)), clamped to [0, 1023]
```

**Step C: TOT (zero, low range, transition, high range, overflow).**

```
TOTlow(Q)  = TOTConst + TOTSlope*Q
TOThigh(Q) = TOTConst_high + TOTSlope_high*Q

if   Qn < Q_TOT:                 TOT = 0                                                      (no TOT)
elif Qn < mTOTTransitionLow:     TOT = TOTlow(Qn)                                             (low range)
elif Qn < mTOTTransitionHigh:    TOT = TOTlow(0.5) + (TOThigh(2.0) - TOTlow(0.5)) * (Qn - 0.5)/(2.0 - 0.5)   (transition, straight line)
else:                            TOT = TOThigh(Qn)                                            (high range)

if TOT > 0:  TOT = int(round(TOT + Gaus(0, mTOTSmear))), clamped to [1, 4095]   (4095 = overflow, Q above ~15 pC)
```

Between `Q_sat` and `Q_TOT` the ADC is saturated and there is still no TOT. This gap is real (test-beam paper, Appendix A), and we want to see its effect.

**Step D: TOA (placeholder for v1).** `TOA = 1` if `Qn > mTOAThreshold`, otherwise `TOA = 0`. A real time model comes later.

**Step E: flags.** Bit 0 = "TOT present" (TOT > 0), bit 1 = "ADC saturated" (ADC = 1023).

**Inverse (for the analysis).** Also add a method `getCharge(absID, ADC, TOT)` that converts a digit back into a charge in pC, using the same constants. You will need it in section 5.

```
if TOT == 0:                                       (use the ADC)
    y = ADC - ADCConst
    if y <= G*Q_lin:            Q = y / G
    elif ADC < 1023:            Q = Q_lin + (Q_sat - Q_lin) * (1 - sqrt(max(0, 1 - (y - G*Q_lin)/D)))   (max: noise can push y above the curve)
    else:                       Q = Q_sat          (saturated and no TOT: only a lower limit)
else:                                              (use the TOT)
    T1 = TOTlow(0.5);  T2 = TOThigh(2.0)
    if   TOT < T1:              Q = (TOT - TOTConst) / TOTSlope
    elif TOT < T2:              Q = 0.5 + (TOT - T1) * (2.0 - 0.5) / (T2 - T1)
    else:                       Q = (TOT - TOTConst_high) / TOTSlope_high     (TOT = 4095: overflow, lower limit)
```

### 3.3 `PadsDigitizer.h/.cxx` (`simulation/`): SDigits → digits

For each collision:

1. Take the `LabeledSDigit`s from the SDigitizer and keep only pads (`isPadDigit()`).
2. **Energy to charge:** `Q[pC] = mCCE * E[eV] / mEnergyPerPair / 6241.5 / 1000`, i.e. Q[fC] = E[eV] / 22780.
3. **Which pads to read out.** The real HGCROC sends **all 72 channels** when it is triggered, even channels with no signal (they show the pedestal). So for every sensor that has at least one SDigit in this collision, digitize all 72 pads of that sensor; pads without an SDigit get Q = 0. You can loop over the 72 pads of a sensor with `PadLocation` and `padIndex.getIndexFromPadLocation(loc)`. With `mReadoutAllSensors = true`, do this for all sensors of the detector (heavy: about 142 000 pads per collision; only for special studies).
4. Call `hgcrocResponse.startCollision()` once, then `HgcrocResponse` for every pad. There is **no zero suppression** in v1: every read-out pad gives a digit.
5. **MC labels:** pads with an SDigit get its labels. Pads without signal get a "noise" label, `MCLabel(true, 1.0)`, as EMCAL does.

### 3.4 The output object

On the branch, the energy-level digit is called **`SDigit`** (`SDigit.h/.cxx`) and the labeled one **`LabeledSDigit`**. The name **`Digit`** is kept free for the real, electronics-level digit, which you write in `DataFormats/Detectors/FOCAL/include/DataFormatsFOCAL/Digit.h` and `src/Digit.cxx`. It is small and simple:

```cpp
uint32_t mAbsID;    // pad index
uint16_t mADC;      // 0-1023
uint16_t mTOA;      // 0-1023
uint16_t mTOT;      // 0-4095
uint8_t  mFlags;    // bit 0: TOT present, bit 1: ADC saturated
// getters return int: getAbsID(), getADC(), getTOA(), getTOT(), getFlags()
```

The FOCAL-H readout chip (H2GCROC) gives exactly the same quantities, so the HCAL people will reuse this class later. Add it (and `std::vector<Digit>`) to `DataFormatsFOCALLinkDef.h` and to `CMakeLists.txt`.

### 3.5 Put it into the workflow

In `FOCALDigitizerSpec.cxx`, after the SDigitizer has run for a collision, call your `PadsDigitizer`, and write instead of the SDigits:

- the `Digit`s;
- one `TriggerRecord` per collision (first digit and number of digits of the collision);
- an `MCTruthContainer<MCLabel>` with the labels, filled with `labels.addElement(digitIndex, label)`.

`FOCALDigitWriterSpec.cxx` currently writes a `std::vector<SDigit>` in the branch `FOCALDigit`, plus the trigger records (`FOCALDigitTRGR`) and the labels (`FOCALDigitMCTruth`). Change only the type of the first branch from `SDigit` to your new `Digit` (and the `#include`). The branch name stays `FOCALDigit`.

**Debug output for validation.** We need the true energy of each pad next to its ADC/TOT, but we don't want SDigits in the final output. Do it like EMCAL, whose digitizer has a "debug stream" (look for `mDebugStream` / `TreeStreamRedirector` in `EMCAL/simulation/src/Digitizer.cxx`). If `mDebugStream` is true, `PadsDigitizer` writes a small tree to `focaldigitsDebug.root`, one entry per digitized pad, with:

`absID, layer, row, col, tbAsic, channel, E (eV), Q (pC), Qn (pC), ADC, TOT, TOA`

This needs no change to the data formats or the writer, and can be removed later.

---

## 4. How to run

**Simulation** (here 100 electrons of 100 GeV; note **PDG code 11** = electron, 211 would be a pion):

```bash
o2-sim -j8 -n100 -m FOC -e TGeant4 --detectorList ALICE2.1 -g external \
  --configKeyValues "GeneratorExternal.fileName=fixedEnergyGun.C;GeneratorExternal.funcName=fixedEnergyGun(100)" --noemptyevents
```

`fixedEnergyGun.C`:

```cpp
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <iostream>
#include "FairBoxGenerator.h"
#endif

FairGenerator* fixedEnergyGun(double energy)
{
  std::cout << "Single electron generator shooting at FOCAL with energy " << energy << " GeV/c" << std::endl;
  auto elecgen = new FairBoxGenerator(11, 1); // PDG 11 = electron, 1 particle per event
  elecgen->SetEtaRange(3.4, 5.5);
  elecgen->SetPhiRange(0., 360.);
  elecgen->SetPRange(energy, energy);
  return elecgen;
}
```

**Digitization:**

```bash
o2-sim-digitizer-workflow -b --onlyDet=FOC --configKeyValues "FOCALPadSimParam.mDebugStream=true"
```

This produces `focaldigits.root` (and `focaldigitsDebug.root` with the debug stream). The interaction rate does not matter for now, because we don't simulate timing.

**Reading the digits.** A skeleton macro, adapted from the one I use for EMCAL:

```cpp
TFile* f = TFile::Open("focaldigits.root");
TTreeReader reader("o2sim", f);
TTreeReaderValue<std::vector<o2::focal::Digit>> digits(reader, "FOCALDigit");
TTreeReaderValue<std::vector<o2::focal::TriggerRecord>> trgs(reader, "FOCALDigitTRGR");
while (reader.Next()) {
  for (auto& trg : *trgs) {                                   // one trigger record = one collision
    gsl::span<const o2::focal::Digit> evDigits(digits->data() + trg.getFirstEntry(), trg.getNumberOfObjects());
    double qSum = 0;
    for (auto& dig : evDigits) {
      qSum += hgcrocResponse.getCharge(dig.getAbsID(), dig.getADC(), dig.getTOT());   // inverse of section 3.2
    }
    hChargeSum->Fill(qSum);
  }
}
```

---

## 5. Validation and plots for v1

1. **Response curves** (standalone test of `HgcrocResponse`, no simulation needed). For one channel, scan Q from 0.001 to 15 pC, with noise switched off, and plot ADC and TOT vs Q. Compare with **Fig. 37** of the test-beam paper. You should see:
   - linear ADC up to ≈ 0.2 pC, then bending until 1023 at ≈ 0.27 pC;
   - the gap with no TOT, then TOT starting at ≈ 0.37 pC;
   - the TOT transition and the linear high range;
   - overflow at 4095.

   Also check that `getCharge(response(Q)) = Q`.
2. **Pedestal and noise.** In the simulated electron events, plot `ADC − ADCConst` for pads with Q = 0 (from the debug tree). It should be a Gaussian centred at 0 with σ ≈ √(1.5² + 4.5²) ≈ 4.7 ADC (incoherent + common-mode noise).
3. **Per-pad check.** Using the debug tree, plot the reconstructed charge `getCharge(ADC, TOT)` vs the true Q of each pad, and the ratio. It should be ≈ 1, with visible effects of the noise, of the ADC saturation/gap region, and of the TOT transition.
4. **Charge-sum distributions, as in Fig. 27.** For electrons of 20, 60, 100, 150, 200 and 300 GeV:
   - sum the reconstructed charge of all pads per event (the macro of section 4) and plot the distributions;
   - in the same plot, show the true charge sum from the debug tree (Q of the SDigits, no electronics).

   Expect a mean of roughly 0.32 pC per GeV.

Linearity and resolution (Figs. 28–29) come later, when we add the timing.

---

## 6. Practical advice

- Work on a git branch. Commit often, with small, understandable commits.
- First make it compile and run on 10 events, then worry about physics.
- Print a few digits (absID, Q, ADC, TOT) for one event and check them by hand against the formulas.
- Our simulation geometry is not the final one, and the calibration is from the test-beam prototype. **We are learning the method, not producing final numbers.**