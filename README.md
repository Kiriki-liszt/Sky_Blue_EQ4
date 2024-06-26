# Sky Blue EQ4

Sky Blue EQ4 is maag-like tone shaping eq.  
Runs in double precision 64-bit internal processing. Also double precision input / output if supported.  
Internal sample rated fixed at 192/176.2kHz in 48/44.1kHz, 96/88.1kHz, 192/176.2kHz sample rates.  
It does run under 44.1kHz, but may have some EQ curve cramping.  
At 48/44.1kHz - x4 oversampling is on with 24 sample latency, and at 96/88.1kHz - x2 oversampling is on with 12 sample latency.  

Windows and Mac, VST3 and AU.  

[![GitHub Release](https://img.shields.io/github/v/release/kiriki-liszt/Sky_Blue_EQ4?style=flat-square&label=Get%20latest%20Release)](https://github.com/Kiriki-liszt/Sky_Blue_EQ4/releases/latest)
[![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/kiriki-liszt/Sky_Blue_EQ4/total?style=flat-square&label=total%20downloads&color=blue)](https://github.com/Kiriki-liszt/Sky_Blue_EQ4/releases/latest)  

[![Static Badge](https://img.shields.io/badge/coffee%20maybe%3F%20%3D%5D%20-gray?style=for-the-badge&logo=buy-me-a-coffee)](https://buymeacoffee.com/kirikiaris)


<img src="https://github.com/Kiriki-liszt/Sky_Blue_EQ4/assets/107096260/97e4e687-3212-4a64-8bfa-03322c38140b"  width="600"/>  

## How to use  

1. Windows

Unzip Win.zip from latest release and copy to "C:\Program Files\Common Files\VST3".  

2. MacOS(Intel tested, Apple Silicon not tested).  

Unzip MacOS.zip from latest release and copy vst3 to "/Library/Audio/Plug-Ins/VST3" and component to "/Library/Audio/Plug-Ins/Components".  

> If it doesn't go well, configure security options in console as  
>  
> ``` console  
> sudo xattr -r -d com.apple.quarantine /Library/Audio/Plug-Ins/VST3/Sky_Blue_EQ4.vst3  
> sudo xattr -r -d com.apple.quarantine /Library/Audio/Plug-Ins/Components/Sky_Blue_EQ4.component
>
> sudo codesign --force --sign - /Library/Audio/Plug-Ins/VST3/Sky_Blue_EQ4.vst3  
> sudo codesign --force --sign - /Library/Audio/Plug-Ins/Components/Sky_Blue_EQ4.component
> ```  
>  
> tested by @jonasborneland [here](https://github.com/Kiriki-liszt/JS_Inflator_to_VST2_VST3/issues/12#issuecomment-1616671177)

## Licensing  

> Q: I would like to share the source code of my VST 3 plug-in/host on GitHub or other such platform.  
>
> * You can choose the GPLv3 license and feel free to share your plug-ins/host's source code including or referencing the VST 3 SDK's sources on GitHub.  
> * **You are allowed to provide a binary form of your plug-ins/host too, provided that you provide its source code as GPLv3 too.**  
> * Note that you have to follow the Steinberg VST usage guidelines.  
>  
> <https://steinbergmedia.github.io/vst3_dev_portal/pages/FAQ/Licensing.html>  

![120C7E3245DD5916ACD2E8E6AD51E8FD_snapshot](https://github.com/Kiriki-liszt/Sky_Blue_EQ4/assets/107096260/142e3c12-cd5f-415d-9b72-8b4f04419633)

VSTSDK 3.7.9 used  
VSTGUI 4.12 used  

## Version logs

v0.0.1   : intial try.  
v0.0.2   : GUI update.  
v1.0.0.b : AU added.  
v1.0.1   : AUv2 plist subtype fixed.  
v1.1.0   : 
  1. Oversampling method change - buffer block oversanpling to single sample oversampling.  
  2. FIR filter using Kaiser-Bessel window.  

## What I've learned

### 1. EQ cramp  

The Frequency Response of EQ is cramping near Nyquist, if not treated.  

One way to fix is adding high frequency contents, but it's not good at Phase Response.  

The other way, is by oversampling it.  
If we oversample internal IIR filter, it'll create Frequency Response at targeted Oversampled Frequency.  
Then EQ Shapes will be fixed, and Phase Response as well.  

<img src="https://raw.githubusercontent.com/Kiriki-liszt/Sky_Blue_EQ4/main/screenshot/box_tone.png"  width="600"/>  

Box tone, x4 oversampling at 48kHz to 192kHz.  

<img src="https://raw.githubusercontent.com/Kiriki-liszt/Sky_Blue_EQ4/main/screenshot/40k_5_FR.png"  width="600"/>  

Sky High Band +5.0 @ 40kHz, Frequency Response, x4 Oversampling at 48kHz to 192kHz.  

<img src="https://raw.githubusercontent.com/Kiriki-liszt/Sky_Blue_EQ4/main/screenshot/40k_5_PR.png"  width="600"/>  

Sky High Band +5.0 @ 40kHz, Phase Response, x4 oversampling at 48kHz to 192kHz.  

### 2. Oversampling without 'Nonlinearity'  

If we have any non-linear process, oversampling is used for antialiasing.  
It'll require massive amount of FIR filter taps for rejecting harmonic contents obove Nyquist.  

But, we're using oversampling only for Frequency Response without any harmonic content generated.  
So we don't need high tap FIRs.  
We just need some nice curve that folds back and become flat again.  
Also the Low Pass FIR filter is only applied once.  
No need for two because the EQ itself does not produce and harmonic content.  

The filters has 49 taps for x2 oversampling with 12 sample latency, and 193 taps for x4 oversampling with 24 sample latency.  


### 3. Parallel EQ and Band-pass filters

As I understand, the need of parallel EQ is for better resolution between two close EQ shape.  
It'll interfere less than Serial EQs.  

Also the original Maag EQs are using parallel EQ topology.  
Maag EQ is consisted of four Band-pass and one High-pass for re-constructing the frequency range, and one High-pass for adding 'Very High Frequency'.  

For example, we can view diffrences using two bands at once.  
I'm using Sonnox Oxford EQ as Serial EQ.  

<img src="https://raw.githubusercontent.com/Kiriki-liszt/Sky_Blue_EQ4/main/screenshot/Maag%20160Hz%20%2B2.5.png"  width="600"/>  

This is 160Hz +2.5.  

<img src="https://raw.githubusercontent.com/Kiriki-liszt/Sky_Blue_EQ4/main/screenshot/Maag%20650Hz%20%2B2.5.png"  width="600"/>  

This is 650Hz +2.5.  

<img src="https://raw.githubusercontent.com/Kiriki-liszt/Sky_Blue_EQ4/main/screenshot/Maag%20Parallel%20vs%20Serial.png"  width="600"/>  

This is 160Hz +2.5 and 650Hz +2.5 both.  
We can distinguish each peak better with Maag - Parallel topology.  

### 4. Different things about Maag EQ (and what I did)  

One thing before explaining what I did, there is something really not common about original Maag EQ.  
When we boost or cut any band, the whole curve moves along.  
This makes mixing choices hard as it boosts, sound become 'better' because of the that gain increase.  
I just want to hear what I do.  

Another thing that I noticed in plugin versions is that the low end is not going to -inf at 0 Hz.  
These does not happen with original gear, for example Maag EQ4M at Sound on Sound review of it.  

So, I designed my EQ in strictly Parallel topology with no input signal blended in, and without getting gain increase schema.  

### 5. EQ algo  

The EQ algorythm for Band-pass and High-pass filters is SVF sugested by Andy Cytomic.  
Used his 'trapezoidal SVF in tangent and state space form'.  

### 6. AU Wrapper

Be carefull of plist.info.  
subtype colision happens.  

## references

1. Sound on Sound - maag EQ4M  
<https://www.soundonsound.com/reviews/maag-audio-eq4m>  

2. IIR based EQ and distortions  
<https://vladgsound.wordpress.com/2013/06/01/iir-based-eq-and-distortions/>  

3. A classification of digital equalizers  
<https://vladgsound.wordpress.com/2015/01/12/a-classification-of-digital-equalizers-draft/>  

4. Nova-67P and parallel equalizers explained  
<https://vladgsound.wordpress.com/2014/08/09/nova-67p-and-parallel-equalizers-explained/>  
Dry mix is not used.

5. Simultaneous solving of all outputs of Linear SVF using trapezoidal integration in state space form  
<https://cytomic.com/files/dsp/SvfLinearTrapAllOutputs.pdf>

## Todo
