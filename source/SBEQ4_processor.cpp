//------------------------------------------------------------------------
// Copyright(c) 2023 yg331.
//------------------------------------------------------------------------

#include "SBEQ4_processor.h"
// #include "SBEQ4_cids.h"

#include "base/source/fstreamer.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"

#include "public.sdk/source/vst/vstaudioprocessoralgo.h"
#include "public.sdk/source/vst/vsthelpers.h"
#include "pluginterfaces/base/futils.h"

#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/sse2.h"

using namespace Steinberg;

namespace yg331 {
	//------------------------------------------------------------------------
	// Sky_Blue_EQ4Processor
	//------------------------------------------------------------------------
	Sky_Blue_EQ4Processor::Sky_Blue_EQ4Processor() :
		fParamZoom(Init_Zoom),
		fParamInputAtten(Init_InputAtten),
		fParamSkyFreq(Init_SkyFreq),
		fParamSky(Init_Sky),
		fParam2k5(Init_2k5),
		fParam650(Init_650),
		fParam160(Init_160),
		fParam40(Init_40),
		fParamSub(Init_Sub),
		bParamBypass(Init_Bypass),
		fParamOS(overSample_1x)
	{
		//--- set the wanted controller for our processor
		setControllerClass(kSky_Blue_EQ4ControllerUID);
	}

	//------------------------------------------------------------------------
	Sky_Blue_EQ4Processor::~Sky_Blue_EQ4Processor()
	{}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Processor::initialize(FUnknown* context)
	{
		// Here the Plug-in will be instantiated

		//---always initialize the parent-------
		tresult result = AudioEffect::initialize(context);
		// if everything Ok, continue
		if (result != kResultOk)
		{
			return result;
		}

		//--- create Audio IO ------
        addAudioInput (STR16 ("Audio Input"), Vst::SpeakerArr::kStereo);
        addAudioOutput (STR16 ("Audio Output"), Vst::SpeakerArr::kStereo);

		/* If you don't need an event bus, you can remove the next line */
		// addEventInput(STR16("Event In"), 1);

		for (int channel = 0; channel < 2; channel++) {
			calcFilter(96000.0,  0.0, 24000.0, dnTap_21, 120.00, dnSample_21[channel].coef); //
			calcFilter(192000.0, 0.0, 24000.0, dnTap_42, 120.00, dnSample_42[channel].coef);

			for (int i = 0; i < dnTap_21; i++) dnSample_21[channel].coef[i] *= 2.0;
			for (int i = 0; i < dnTap_42; i++) dnSample_42[channel].coef[i] *= 4.0;
		}

		return kResultOk;
	}
 
 //------------------------------------------------------------------------
 tresult PLUGIN_API Sky_Blue_EQ4Processor::setBusArrangements(
     Vst::SpeakerArrangement* inputs, int32 numIns,
     Vst::SpeakerArrangement* outputs, int32 numOuts)
 {
     if (numIns == 1 && numOuts == 1)
     {
         // the host wants Mono => Mono (or 1 channel -> 1 channel)
         if (Vst::SpeakerArr::getChannelCount(inputs[0]) == 1 &&
             Vst::SpeakerArr::getChannelCount(outputs[0]) == 1)
         {
             auto* bus = FCast<Vst::AudioBus>(audioInputs.at(0));
             if (bus)
             {
                 // check if we are Mono => Mono, if not we need to recreate the busses
                 if (bus->getArrangement() != inputs[0])
                 {
                     getAudioInput(0)->setArrangement(inputs[0]);
                     getAudioInput(0)->setName(STR16("Mono In"));
                     getAudioOutput(0)->setArrangement(outputs[0]);
                     getAudioOutput(0)->setName(STR16("Mono Out"));
                 }
                 return kResultOk;
             }
         }
         // the host wants something else than Mono => Mono,
         // in this case we are always Stereo => Stereo
         else
         {
             auto* bus = FCast<Vst::AudioBus>(audioInputs.at(0));
             if (bus)
             {
                 tresult result = kResultFalse;

                 // the host wants 2->2 (could be LsRs -> LsRs)
                 if (Vst::SpeakerArr::getChannelCount(inputs[0]) == 2 &&
                     Vst::SpeakerArr::getChannelCount(outputs[0]) == 2)
                 {
                     getAudioInput(0)->setArrangement(inputs[0]);
                     getAudioInput(0)->setName(STR16("Stereo In"));
                     getAudioOutput(0)->setArrangement(outputs[0]);
                     getAudioOutput(0)->setName(STR16("Stereo Out"));
                     result = kResultTrue;
                 }
                 // the host want something different than 1->1 or 2->2 : in this case we want stereo
                 else if (bus->getArrangement() != Vst::SpeakerArr::kStereo)
                 {
                     getAudioInput(0)->setArrangement(Vst::SpeakerArr::kStereo);
                     getAudioInput(0)->setName(STR16("Stereo In"));
                     getAudioOutput(0)->setArrangement(Vst::SpeakerArr::kStereo);
                     getAudioOutput(0)->setName(STR16("Stereo Out"));
                     result = kResultFalse;
                 }

                 return result;
             }
         }
     }
     return kResultFalse;
 }


	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Processor::terminate()
	{
		// Here the Plug-in will be de-instantiated, last possibility to remove some memory!
		for (int channel = 0; channel < 2; channel++) {

		}
		//---do not forget to call parent ------
		return AudioEffect::terminate();
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Processor::setActive(TBool state)
	{
		//--- called when the Plug-in is enable/disable (On/Off) -----
		return AudioEffect::setActive(state);
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Processor::process(Vst::ProcessData& data)
	{
		Vst::IParameterChanges* paramChanges = data.inputParameterChanges;

		if (paramChanges)
		{
			int32 numParamsChanged = paramChanges->getParameterCount();

			for (int32 index = 0; index < numParamsChanged; index++)
			{
				Vst::IParamValueQueue* paramQueue = paramChanges->getParameterData(index);

				if (paramQueue)
				{
					Vst::ParamValue value;
					int32 sampleOffset;
					int32 numPoints = paramQueue->getPointCount();

					if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue) {
						switch (paramQueue->getParameterId()) {
						case kParamBypass:
							bParamBypass = (value > 0.5f);
							break;
						case kParamZoom:
							fParamZoom = value;
							break;
						case kParamInputAtten:
							fParamInputAtten = value;
							break;
						case kParamSkyFreq:
							fParamSkyFreq = value;
							break;
						case kParamSky:
							fParamSky = value;
							break;
						case kParam2k5:
							fParam2k5 = value;
							break;
						case kParam650:
							fParam650 = value;
							break;
						case kParam160:
							fParam160 = value;
							break;
						case kParam40:
							fParam40 = value;
							break;
						case kParamSub:
							fParamSub = value;
							break;
						}
					}
				}
			}
		}

		if (data.numInputs == 0 || data.numOutputs == 0) {
			return kResultOk;
		}

		// (simplification) we suppose in this example that we have the same input channel count than the output
		int32 numChannels = data.inputs[0].numChannels;

		//---get audio buffers----------------
		uint32 sampleFramesSize = getSampleFramesSizeInBytes(processSetup, data.numSamples);
		void** in = getChannelBuffersPointer(processSetup, data.inputs[0]);
		void** out = getChannelBuffersPointer(processSetup, data.outputs[0]);
		Vst::SampleRate getSampleRate = processSetup.sampleRate;

		//---check if silence---------------
		// check if all channel are silent then process silent
		if (data.inputs[0].silenceFlags == Vst::getChannelMask(data.inputs[0].numChannels))
		{
			// mark output silence too (it will help the host to propagate the silence)
			data.outputs[0].silenceFlags = data.inputs[0].silenceFlags;

			// the plug-in has to be sure that if it sets the flags silence that the output buffer are clear
			for (int32 channel = 0; channel < numChannels; channel++)
			{
				// do not need to be cleared if the buffers are the same (in this case input buffer are already cleared by the host)
				if (in[channel] != out[channel])
				{
					memset(out[channel], 0, sampleFramesSize);
				}
			}
			return kResultOk;
		}

		data.outputs[0].silenceFlags = data.inputs[0].silenceFlags;

		if (data.symbolicSampleSize == Vst::kSample32) {
			processSVF<Vst::Sample32>((Vst::Sample32**)in, (Vst::Sample32**)out, numChannels, getSampleRate, data.numSamples);
		}
		else {
			processSVF<Vst::Sample64>((Vst::Sample64**)in, (Vst::Sample64**)out, numChannels, getSampleRate, data.numSamples);
		}

		return kResultOk;
	}

	uint32 PLUGIN_API Sky_Blue_EQ4Processor::getLatencySamples()
	{
		if (fParamOS == overSample_1x) return 0;
		else if (fParamOS == overSample_2x) return latency_Fir_x2;
		else return latency_Fir_x4;
	}

	tresult PLUGIN_API Sky_Blue_EQ4Processor::setupProcessing(Vst::ProcessSetup& newSetup)
	{
		Vst::ParamValue OS_target = 0.0;
		if (newSetup.sampleRate <= 48000.0) {
			fParamOS = overSample_4x;
			OS_target = 4 * newSetup.sampleRate;
		}
		else if (newSetup.sampleRate <= 96000.0) {
			fParamOS = overSample_2x;
			OS_target = 2 * newSetup.sampleRate;
		}
		else {
			fParamOS = overSample_1x;
			OS_target = newSetup.sampleRate;
		}

		const Vst::ParamValue Fc_sub = 10.0;
		const Vst::ParamValue Fc_40  = 39.0;
		const Vst::ParamValue Fc_160 = 155.0;
		const Vst::ParamValue Fc_650 = 625.0;
		const Vst::ParamValue Fc_2k5 = 12.5;

		for (int channel = 0; channel < 2; channel++) {
			setSVF(&svf_Sub[channel], kBP, 0.465, 1.0, Fc_sub, OS_target);
			setSVF(&svf_40[channel], kBP, 0.465, 1.0, Fc_40, OS_target);
			setSVF(&svf_160[channel], kBP, 0.465, 1.0, Fc_160, OS_target);
			setSVF(&svf_650[channel], kBP, 0.45, 1.0, Fc_650, OS_target);
			setSVF(&svf_2k5[channel], kHP, 0.01, 1.0, Fc_2k5, OS_target);
			//setSVF(&svf_Sky[channel], kHP, 0.01, 1.0, Fc_air, OS_target);
			//setSVF(&svf_LC[channel], kHP, 0.72, 1.0, 0.3, OS_target);
			//setSVF(&svf_HC[channel], kLP, 0.1, 1.0, 72000, OS_target);
		}

		//--- called before any processing ----
		return AudioEffect::setupProcessing(newSetup);
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Processor::canProcessSampleSize(int32 symbolicSampleSize)
	{
		// by default kSample32 is supported
		if (symbolicSampleSize == Vst::kSample32)
			return kResultTrue;

		// disable the following comment if your processing support kSample64
		if (symbolicSampleSize == Vst::kSample64)
			return kResultTrue; 

		return kResultFalse;
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Processor::setState(IBStream* state)
	{
		// called when we load a preset, the model has to be reloaded
		IBStreamer streamer(state, kLittleEndian);

		int32           savedBypass = 0;
		Vst::ParamValue savedZoom = 0.0;

		Vst::ParamValue savedInputAtten = 0.0;
		Vst::ParamValue savedSkyFreq = 0.0;
		Vst::ParamValue savedSky = 0.0;
		Vst::ParamValue saved2k5 = 0.0;
		Vst::ParamValue saved650 = 0.0;
		Vst::ParamValue saved160 = 0.0;
		Vst::ParamValue saved40 = 0.0;
		Vst::ParamValue savedSub = 0.0;


		if (streamer.readInt32(savedBypass) == false) return kResultFalse;
		if (streamer.readDouble(savedZoom) == false) return kResultFalse;

		if (streamer.readDouble(savedInputAtten) == false) return kResultFalse;
		if (streamer.readDouble(savedSkyFreq) == false) return kResultFalse;
		if (streamer.readDouble(savedSky) == false) return kResultFalse;
		if (streamer.readDouble(saved2k5) == false) return kResultFalse;
		if (streamer.readDouble(saved650) == false) return kResultFalse;
		if (streamer.readDouble(saved160) == false) return kResultFalse;
		if (streamer.readDouble(saved40) == false) return kResultFalse;
		if (streamer.readDouble(savedSub) == false) return kResultFalse;


		bParamBypass = savedBypass > 0;
		fParamZoom = savedZoom;
		fParamInputAtten = savedInputAtten;
		fParamSkyFreq = savedSkyFreq;
		fParamSky = savedSky;
		fParam2k5 = saved2k5;
		fParam650 = saved650;
		fParam160 = saved160;
		fParam40 = saved40;
		fParamSub = savedSub;

		return kResultOk;
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Processor::getState(IBStream* state)
	{
		// here we need to save the model
		IBStreamer streamer(state, kLittleEndian);

		streamer.writeInt32(bParamBypass ? 1 : 0);
		streamer.writeDouble(fParamZoom);

		streamer.writeDouble(fParamInputAtten);
		streamer.writeDouble(fParamSkyFreq);
		streamer.writeDouble(fParamSky);
		streamer.writeDouble(fParam2k5);
		streamer.writeDouble(fParam650);
		streamer.writeDouble(fParam160);
		streamer.writeDouble(fParam40);
		streamer.writeDouble(fParamSub);

		return kResultOk;
	}




	// User made --------


	inline void Sky_Blue_EQ4Processor::setSVF
	(SVF_* svf_filter, filter_type kFilter, Vst::ParamValue q, Vst::ParamValue gain, Vst::Sample64 Fc, Vst::Sample64 Fs)
	{
		svf_filter->bellgaindB = gain;
		svf_filter->A = pow(10.0, svf_filter->bellgaindB / 40.0);

		svf_filter->w = Fc * M_PI / Fs;
		svf_filter->Q = q; // 0.5 == 6db/oct
		svf_filter->g0 = tan(svf_filter->w);
		svf_filter->k0 = 1 / svf_filter->Q;

		svf_filter->type = kFilter;
		svf_filter->g = svf_filter->g0;
		svf_filter->k = svf_filter->k0;

		switch (svf_filter->type)
		{
		case kLP:
			svf_filter->m0 = 0;
			svf_filter->m1 = 0;
			svf_filter->m2 = 1;
			break;
		case kHP:
			svf_filter->m0 = 1;
			svf_filter->m1 = 0;
			svf_filter->m2 = 0;
			break;
		case kBP:
			svf_filter->m0 = 0;
			svf_filter->m1 = 1;
			svf_filter->m2 = 0;
			break;
		case kNotch:
			svf_filter->m0 = 1;
			svf_filter->m1 = 0;
			svf_filter->m2 = 1;
			break;
		case kPeak:
			svf_filter->m0 = 1;
			svf_filter->m1 = 0;
			svf_filter->m2 = -1;
			break;
		case kBell:
			svf_filter->g = svf_filter->g0;
			svf_filter->k = svf_filter->k0 / svf_filter->A;
			svf_filter->m0 = 1;
			svf_filter->m1 = svf_filter->k0 * svf_filter->A;
			svf_filter->m2 = 1;
			break;
		case kLS:
			svf_filter->g = svf_filter->g0 / sqrt(svf_filter->A);
			svf_filter->k = svf_filter->k0;
			svf_filter->m0 = 1;
			svf_filter->m1 = svf_filter->k0 * svf_filter->A;
			svf_filter->m2 = svf_filter->A * svf_filter->A;
			break;
		case kHS:
			svf_filter->g = svf_filter->g0 * sqrt(svf_filter->A);
			svf_filter->k = svf_filter->k0;
			svf_filter->m0 = svf_filter->A * svf_filter->A;
			svf_filter->m1 = svf_filter->k0 * svf_filter->A;
			svf_filter->m2 = 1;
			break;
		default:
			break;
		}
		svf_filter->g_k = svf_filter->g + svf_filter->k;
		svf_filter->gt0 = 1 / (1 + svf_filter->g * svf_filter->g_k);
		svf_filter->gk0 = svf_filter->g_k * svf_filter->gt0;
		svf_filter->gt1 = svf_filter->g * svf_filter->gt0;
		svf_filter->gk1 = svf_filter->g * svf_filter->gk0;
		svf_filter->gt2 = svf_filter->g * svf_filter->gt1;
	}

	inline Vst::Sample64 Sky_Blue_EQ4Processor::computeSVF
	(SVF_* svf_filter, Vst::Sample64 input)
	{
		svf_filter->vin = input;
		svf_filter->t0 = svf_filter->vin - svf_filter->ic2eq;
		svf_filter->v0 = svf_filter->gt0 * svf_filter->t0 - svf_filter->gk0 * svf_filter->ic1eq;
		svf_filter->t1 = svf_filter->gt1 * svf_filter->t0 - svf_filter->gk1 * svf_filter->ic1eq;
		svf_filter->t2 = svf_filter->gt2 * svf_filter->t0 + svf_filter->gt1 * svf_filter->ic1eq;
		svf_filter->v1 = svf_filter->t1 + svf_filter->ic1eq;
		svf_filter->v2 = svf_filter->t2 + svf_filter->ic2eq;
		svf_filter->ic1eq += 2.0 * svf_filter->t1;
		svf_filter->ic2eq += 2.0 * svf_filter->t2;

		return svf_filter->m0 * svf_filter->v0 + svf_filter->m1 * svf_filter->v1 + svf_filter->m2 * svf_filter->v2;
	}

	inline void Sky_Blue_EQ4Processor::coefSVF
	(Vst::SampleRate Fs)
	{
		int32 iParamAirF = FromNormalized<Vst::ParamValue>(fParamSkyFreq, 6);
		double dAir_f[7] = { 0.0, 19.5, 39.0, 78.0, 110.0, 165.0, 350.0 };

		//double Fc_sub = 10.0;
		//double Fc_40 = 39.0;
		//double Fc_160 = 155.0;
		//double Fc_650 = 625.0;
		//double Fc_2k5 = 12.5;
		double Fc_air = dAir_f[iParamAirF];

		for (int channel = 0; channel < 2; channel++) {
			//setSVF(&svf_Sub[channel], kBP, 0.465, 1.0, Fc_sub, Fs);
			//setSVF(&svf_40[channel], kBP, 0.465, 1.0, Fc_40, Fs);
			//setSVF(&svf_160[channel], kBP, 0.465, 1.0, Fc_160, Fs);
			//setSVF(&svf_650[channel], kBP, 0.45, 1.0, Fc_650, Fs);
			//setSVF(&svf_2k5[channel], kHP, 0.01, 1.0, Fc_2k5, Fs);
			setSVF(&svf_Sky[channel], kHP, 0.01, 1.0, Fc_air, Fs);
			//setSVF(&svf_LC[channel], kHP, 0.72, 1.0, 0.3, Fs);
			//setSVF(&svf_HC[channel], kLP, 0.1, 1.0, 72000, Fs);
		}
		return;
	}


	template <typename SampleType>
	void Sky_Blue_EQ4Processor::processSVF(
		SampleType**     inputs, 
		SampleType**     outputs,
		Steinberg::int32 numChannels,
		Steinberg::Vst::SampleRate getSampleRate, 
		Steinberg::int32 sampleFrames)
	{

		Vst::Sample64 _db = (12.0 * fParamInputAtten - 12.0) - (3.0);

		Vst::Sample64 In_Atten = exp(log(10.0) * _db / 20.0);

		coefSVF(getSampleRate);

		int32 nParam_Sub = FromNormalized<Vst::ParamValue>(fParamSub, knob_stepCount);
		int32 nParam_40  = FromNormalized<Vst::ParamValue>(fParam40,  knob_stepCount);
		int32 nParam_160 = FromNormalized<Vst::ParamValue>(fParam160, knob_stepCount);
		int32 nParam_650 = FromNormalized<Vst::ParamValue>(fParam650, knob_stepCount);
		int32 nParam_2k5 = FromNormalized<Vst::ParamValue>(fParam2k5, knob_stepCount);
		int32 nParam_Sky = FromNormalized<Vst::ParamValue>(fParamSky, knob_stepCount);

		BP_gain_Sub = BP_arr[nParam_Sub];
		BP_gain_40  = BP_arr[nParam_40];
		BP_gain_160 = BP_arr[nParam_160];
		BP_gain_650 = BP_arr[nParam_650];
		HP_gain_2k5 = HP_arr[nParam_2k5];
		HP_gain_Sky = HP_Sky_arr[nParam_Sky];
		if (fParamSkyFreq == 0.0) HP_gain_Sky = 0.0;

		int32 oversampling = 1;
		if (fParamOS == overSample_2x) oversampling = 2;
		else if (fParamOS == overSample_4x) oversampling = 4;

		int32 latency = 0;
		if (fParamOS == overSample_2x) latency = latency_Fir_x2;
		else if (fParamOS == overSample_4x) latency = latency_Fir_x4;

		for (int32 channel = 0; channel < numChannels; channel++)
		{
			int32 samples = sampleFrames;

			SampleType* ptrIn  = (SampleType*)inputs[channel];
			SampleType* ptrOut = (SampleType*)outputs[channel];

			if (latency != latency_q[channel].size()) {
				int32 diff = latency - (int32)latency_q[channel].size();
				if   (diff > 0) for (int i = 0; i <  diff; i++) latency_q[channel].push(0.0);
				else            for (int i = 0; i < -diff; i++) latency_q[channel].pop();
			}

			while (--samples >= 0)
			{
				Vst::Sample64 inputSample = *ptrIn;
				ptrIn++;
				Vst::Sample64 drySample = inputSample;
				inputSample *= In_Atten;

				double up_x[4] = { 0.0, };
				double up_y[4] = { 0.0, };

				up_x[0] = inputSample;

				// Process
				for (int i = 0; i < oversampling; i++)
				{
					Vst::Sample64 overSampled = up_x[i];

					//inputSample = computeSVF(&svf_LC[channel], overSampled);

					svf_Sub[channel].vin = overSampled;
					svf_Sub[channel].t0 = svf_Sub[channel].vin - svf_Sub[channel].ic2eq;
					svf_Sub[channel].v0 = svf_Sub[channel].gt0 * svf_Sub[channel].t0 - svf_Sub[channel].gk0 * svf_Sub[channel].ic1eq;
					svf_Sub[channel].t1 = svf_Sub[channel].gt1 * svf_Sub[channel].t0 - svf_Sub[channel].gk1 * svf_Sub[channel].ic1eq;
					svf_Sub[channel].t2 = svf_Sub[channel].gt2 * svf_Sub[channel].t0 + svf_Sub[channel].gt1 * svf_Sub[channel].ic1eq;
					svf_Sub[channel].v1 = svf_Sub[channel].t1 + svf_Sub[channel].ic1eq;
					svf_Sub[channel].v2 = svf_Sub[channel].t2 + svf_Sub[channel].ic2eq;
					svf_Sub[channel].ic1eq += 2.0 * svf_Sub[channel].t1;
					svf_Sub[channel].ic2eq += 2.0 * svf_Sub[channel].t2;
					Vst::Sample64 tmp_sub = svf_Sub[channel].m0 * svf_Sub[channel].v0 
					                      + svf_Sub[channel].m1 * svf_Sub[channel].v1 
					                      + svf_Sub[channel].m2 * svf_Sub[channel].v2;

					svf_40[channel].vin = overSampled;
					svf_40[channel].t0 = svf_40[channel].vin - svf_40[channel].ic2eq;
					svf_40[channel].v0 = svf_40[channel].gt0 * svf_40[channel].t0 - svf_40[channel].gk0 * svf_40[channel].ic1eq;
					svf_40[channel].t1 = svf_40[channel].gt1 * svf_40[channel].t0 - svf_40[channel].gk1 * svf_40[channel].ic1eq;
					svf_40[channel].t2 = svf_40[channel].gt2 * svf_40[channel].t0 + svf_40[channel].gt1 * svf_40[channel].ic1eq;
					svf_40[channel].v1 = svf_40[channel].t1 + svf_40[channel].ic1eq;
					svf_40[channel].v2 = svf_40[channel].t2 + svf_40[channel].ic2eq;
					svf_40[channel].ic1eq += 2.0 * svf_40[channel].t1;
					svf_40[channel].ic2eq += 2.0 * svf_40[channel].t2;
					Vst::Sample64 tmp_40 = svf_40[channel].m0 * svf_40[channel].v0 
					                     + svf_40[channel].m1 * svf_40[channel].v1 
					                     + svf_40[channel].m2 * svf_40[channel].v2;

					svf_160[channel].vin = overSampled;
					svf_160[channel].t0 = svf_160[channel].vin - svf_160[channel].ic2eq;
					svf_160[channel].v0 = svf_160[channel].gt0 * svf_160[channel].t0 - svf_160[channel].gk0 * svf_160[channel].ic1eq;
					svf_160[channel].t1 = svf_160[channel].gt1 * svf_160[channel].t0 - svf_160[channel].gk1 * svf_160[channel].ic1eq;
					svf_160[channel].t2 = svf_160[channel].gt2 * svf_160[channel].t0 + svf_160[channel].gt1 * svf_160[channel].ic1eq;
					svf_160[channel].v1 = svf_160[channel].t1 + svf_160[channel].ic1eq;
					svf_160[channel].v2 = svf_160[channel].t2 + svf_160[channel].ic2eq;
					svf_160[channel].ic1eq += 2.0 * svf_160[channel].t1;
					svf_160[channel].ic2eq += 2.0 * svf_160[channel].t2;
					Vst::Sample64 tmp_160 = svf_160[channel].m0 * svf_160[channel].v0 
					                      + svf_160[channel].m1 * svf_160[channel].v1 
					                      + svf_160[channel].m2 * svf_160[channel].v2;

					svf_650[channel].vin = overSampled;
					svf_650[channel].t0 = svf_650[channel].vin - svf_650[channel].ic2eq;
					svf_650[channel].v0 = svf_650[channel].gt0 * svf_650[channel].t0 - svf_650[channel].gk0 * svf_650[channel].ic1eq;
					svf_650[channel].t1 = svf_650[channel].gt1 * svf_650[channel].t0 - svf_650[channel].gk1 * svf_650[channel].ic1eq;
					svf_650[channel].t2 = svf_650[channel].gt2 * svf_650[channel].t0 + svf_650[channel].gt1 * svf_650[channel].ic1eq;
					svf_650[channel].v1 = svf_650[channel].t1 + svf_650[channel].ic1eq;
					svf_650[channel].v2 = svf_650[channel].t2 + svf_650[channel].ic2eq;
					svf_650[channel].ic1eq += 2.0 * svf_650[channel].t1;
					svf_650[channel].ic2eq += 2.0 * svf_650[channel].t2;
					Vst::Sample64 tmp_650 = svf_650[channel].m0 * svf_650[channel].v0 
					                      + svf_650[channel].m1 * svf_650[channel].v1 
					                      + svf_650[channel].m2 * svf_650[channel].v2;

					svf_2k5[channel].vin = overSampled;
					svf_2k5[channel].t0 = svf_2k5[channel].vin - svf_2k5[channel].ic2eq;
					svf_2k5[channel].v0 = svf_2k5[channel].gt0 * svf_2k5[channel].t0 - svf_2k5[channel].gk0 * svf_2k5[channel].ic1eq;
					svf_2k5[channel].t1 = svf_2k5[channel].gt1 * svf_2k5[channel].t0 - svf_2k5[channel].gk1 * svf_2k5[channel].ic1eq;
					svf_2k5[channel].t2 = svf_2k5[channel].gt2 * svf_2k5[channel].t0 + svf_2k5[channel].gt1 * svf_2k5[channel].ic1eq;
					svf_2k5[channel].v1 = svf_2k5[channel].t1 + svf_2k5[channel].ic1eq;
					svf_2k5[channel].v2 = svf_2k5[channel].t2 + svf_2k5[channel].ic2eq;
					svf_2k5[channel].ic1eq += 2.0 * svf_2k5[channel].t1;
					svf_2k5[channel].ic2eq += 2.0 * svf_2k5[channel].t2;
					Vst::Sample64 tmp_2k5 = svf_2k5[channel].m0 * svf_2k5[channel].v0 
					                      + svf_2k5[channel].m1 * svf_2k5[channel].v1 
					                      + svf_2k5[channel].m2 * svf_2k5[channel].v2;

					svf_Sky[channel].vin = overSampled;
					svf_Sky[channel].t0 = svf_Sky[channel].vin - svf_Sky[channel].ic2eq;
					svf_Sky[channel].v0 = svf_Sky[channel].gt0 * svf_Sky[channel].t0 - svf_Sky[channel].gk0 * svf_Sky[channel].ic1eq;
					svf_Sky[channel].t1 = svf_Sky[channel].gt1 * svf_Sky[channel].t0 - svf_Sky[channel].gk1 * svf_Sky[channel].ic1eq;
					svf_Sky[channel].t2 = svf_Sky[channel].gt2 * svf_Sky[channel].t0 + svf_Sky[channel].gt1 * svf_Sky[channel].ic1eq;
					svf_Sky[channel].v1 = svf_Sky[channel].t1 + svf_Sky[channel].ic1eq;
					svf_Sky[channel].v2 = svf_Sky[channel].t2 + svf_Sky[channel].ic2eq;
					svf_Sky[channel].ic1eq += 2.0 * svf_Sky[channel].t1;
					svf_Sky[channel].ic2eq += 2.0 * svf_Sky[channel].t2;
					Vst::Sample64 tmp_sky = svf_Sky[channel].m0 * svf_Sky[channel].v0 
					                      + svf_Sky[channel].m1 * svf_Sky[channel].v1 
					                      + svf_Sky[channel].m2 * svf_Sky[channel].v2;

					Vst::Sample64 dataOut = 0.0;
					dataOut += tmp_sub * BP_gain_Sub;
					dataOut += tmp_40  * BP_gain_40;
					dataOut += tmp_160 * BP_gain_160;
					dataOut += tmp_650 * BP_gain_650;
					dataOut += tmp_2k5 * HP_gain_2k5;
					dataOut += tmp_sky * HP_gain_Sky;

					up_y[i] = dataOut; // * gain;
				}

				// Downsampling
				if (fParamOS == overSample_1x) inputSample = up_y[0];
				else if (fParamOS == overSample_2x) Fir_x2_dn(up_y, &inputSample, channel);
				else Fir_x4_dn(up_y, &inputSample, channel);

				// Latency compensate
				Vst::Sample64 delayed;
				if (fParamOS == overSample_1x) {
					delayed = drySample;
				}
				else {
					delayed = latency_q[channel].front();
					latency_q[channel].pop();
					latency_q[channel].push(drySample);
				}

				if (bParamBypass) {
					inputSample = delayed;
				}

				*ptrOut = (SampleType)inputSample;

				ptrOut++;
			}
		}
		return;
	}


	void Sky_Blue_EQ4Processor::Fir_x2_dn(
		Vst::Sample64* in, 
		Vst::Sample64* out, 
		Steinberg::int32 channel
	) 
	{
		double inter_21[2];

		const size_t dnTap_21_size = sizeof(double) * (dnTap_21);
		memmove(dnSample_21[channel].buff + 3, dnSample_21[channel].buff + 1, dnTap_21_size);
		dnSample_21[channel].buff[2] = in[0];
		dnSample_21[channel].buff[1] = in[1];
		__m128d _acc_out = _mm_setzero_pd();
		for (int i = 0; i < dnTap_21; i += 2) {
			__m128d coef = _mm_load_pd(&dnSample_21[channel].coef[i]);
			__m128d buff = _mm_load_pd(&dnSample_21[channel].buff[i + 2]);
			__m128d _mul = _mm_mul_pd(coef, buff);
			_acc_out = _mm_add_pd(_acc_out, _mul);
		}
		_mm_store_pd(inter_21, _acc_out);
		*out = inter_21[0] + inter_21[1];
		return;
	}

	void Sky_Blue_EQ4Processor::Fir_x4_dn(
		Vst::Sample64* in,
		Vst::Sample64* out,
		Steinberg::int32 channel
	)
	{
		double inter_42[2];

		int32 dnTap_42_size = sizeof(double) * (dnTap_42);
		memmove(dnSample_42[channel].buff + 5, dnSample_42[channel].buff + 1, dnTap_42_size);
		dnSample_42[channel].buff[4] = in[0]; // buff[3]
		dnSample_42[channel].buff[3] = in[1];
		dnSample_42[channel].buff[2] = in[2];
		dnSample_42[channel].buff[1] = in[3];
		__m128d _acc_out_a = _mm_setzero_pd();
		for (int i = 0; i < dnTap_42; i += 2) { // tt >= dnTap_42+3
			__m128d coef_a = _mm_load_pd(&dnSample_42[channel].coef[i]);
			__m128d buff_a = _mm_load_pd(&dnSample_42[channel].buff[i + 4]);
			__m128d _mul_a = _mm_mul_pd(coef_a, buff_a);
			_acc_out_a = _mm_add_pd(_acc_out_a, _mul_a);
		}
		_mm_store_pd(inter_42, _acc_out_a);
		*out = inter_42[0] + inter_42[1];
		return;
	}

	//------------------------------------------------------------------------
} // namespace yg331
