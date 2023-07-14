# Sky Blue EQ4

Sky Blue EQ4 is maag-like tone shaping eq.  
Internal sample rated fixed at 176.2/192kHz in 44.1/48kHz, 88.1/96kHz, 176.2/192kHz sample rates.  
It does run under 44.1kHz, but may have some EQ curve cramping.  
At 44.1/48kHz - x4 oversampling is on with 24 sample latency, and at 88.1/96kHz - x2 oversampling is on with 12 sample latency.  

<img src="https://github.com/Kiriki-liszt/Sky_Blue_EQ4/assets/107096260/97e4e687-3212-4a64-8bfa-03322c38140b"  width="600"/>  

Windows and MacOS(Intel tested, Apple Silicon not tested).  
I tried to Code Sign it, donn't know if it works.  
VST3.  

## Licensing  

> Q: I would like to share the source code of my VST 3 plug-in/host on GitHub or other such platform.  
>
> * You can choose the GPLv3 license and feel free to share your plug-ins/host's source code including or referencing the VST 3 SDK's sources on GitHub.  
> * **You are allowed to provide a binary form of your plug-ins/host too, provided that you provide its source code as GPLv3 too.**  
> * Note that you have to follow the Steinberg VST usage guidelines.  
>  
> <https://steinbergmedia.github.io/vst3_dev_portal/pages/FAQ/Licensing.html>  

![120C7E3245DD5916ACD2E8E6AD51E8FD_snapshot](https://github.com/Kiriki-liszt/Sky_Blue_EQ4/assets/107096260/142e3c12-cd5f-415d-9b72-8b4f04419633)


VSTSDK 3.7.7 used  
VSTGUI 4.12 used  

## Version logs

v0.0.1  : intial try.  
v0.0.2  : GUI update.  

## What I've learned

### 1. EQ cramp  

The Frequency Response of EQ is cramping near Nyquist, if not treated.  

One way to fix is adding high frequency contents, but it's not good at Phase Response.  

The other way, is by oversampling it.  
If we oversample internal IIR filter, it'll create Frequency Response at targeted Oversampled Frequency.  
Then EQ Shapes will be fixed, and Phase Response as well.  

<img src="https://github.com/Kiriki-liszt/Sky_Blue_EQ4/assets/107096260/86d167b5-8038-41fe-b4a3-d66b1166faa6"  width="600"/>  

Box tone, x4 oversampling at 48kHz to 196kHz.  

<img src="https://github.com/Kiriki-liszt/Sky_Blue_EQ4/assets/107096260/e2284855-84f1-40a3-aa0f-a3d18decb353"  width="600"/>  

Sky High Band +5.0 @ 40kHz, Frequency Response, x4 Oversampling at 48kHz to 196kHz.  

<img src="https://github.com/Kiriki-liszt/Sky_Blue_EQ4/assets/107096260/b1efedcb-13f6-4ba3-8e48-913f4c454e24"  width="600"/>  

Sky High Band +5.0 @ 40kHz, Phase Response, x4 oversampling at 48kHz to 196kHz.  

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

However, The Plugin-alliance version is not strictly Parallel topology.  
We can check this out with comparing with any Serial EQs like Pro-Q 3, as PA version Maag EQ curve can be recreated by simple Bell boosting.  

### 4. Different things about Maag EQ (and what I did)  

One thing before explaining what I did, there is something really not common about original Maag EQ.  
When we boost or cut any band, the whole curve moves along.  
This makes mixing choices hard as it boosts, sound become 'better' because of the that gain increase.  
I just want to hear what I do.  

And maybe because of this characteristic, the PA version recreated that 'gain increase as boosting' by blending in Original signal.  
This makes diffrence when boosting, Q for other band gets smaller, losing it's shape.  
Another thing that I noticed is the low end is not going to -inf at 0 Hz.  
These does not happen with original gear, for example Maag EQ4M at Sound on Sound review of it.  

So, I designed my EQ in strictly Parallel topology with no input signal blended in, and without getting gain increase schema.  

### EQ algo  

The EQ algorythm for Band-pass and High-pass filters is SVF sugested by Andy Cytomic.  
Used his 'trapezoidal SVF in tangent and state space form'.  

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
