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
		fParamOS(overSample_4x)
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
		addAudioInput(STR16("Stereo In"), Steinberg::Vst::SpeakerArr::kStereo);
		addAudioOutput(STR16("Stereo Out"), Steinberg::Vst::SpeakerArr::kStereo);

		/* If you don't need an event bus, you can remove the next line */
		// addEventInput(STR16("Event In"), 1);

		out_0 = (Vst::Sample64**)malloc(sizeof(Vst::Sample64*) * 2);
		out_1 = (Vst::Sample64**)malloc(sizeof(Vst::Sample64*) * 2);
		out_2 = (Vst::Sample64**)malloc(sizeof(Vst::Sample64*) * 2);
		if ((out_0 == NULL) || (out_1 == NULL) || (out_2 == NULL) ) return kResultFalse;

		for (int c = 0; c < 2; c++) {
			out_0[c] = (Vst::Sample64*)malloc(sizeof(Vst::Sample64) * maxSample);
			out_1[c] = (Vst::Sample64*)malloc(sizeof(Vst::Sample64) * maxSample * 2);
			out_2[c] = (Vst::Sample64*)malloc(sizeof(Vst::Sample64) * maxSample * 4);
			if ((out_0[c] == NULL) || (out_1[c] == NULL) || (out_2[c] == NULL) ) return kResultFalse;
			memset(out_0[c], 0, maxSample);
			memset(out_1[c], 0, maxSample * 2);
			memset(out_2[c], 0, maxSample * 4);
		}

		return kResultOk;
	}

	tresult PLUGIN_API Sky_Blue_EQ4Processor::setBusArrangements(
		Vst::SpeakerArrangement* inputs, int32 numIns,
		Vst::SpeakerArrangement* outputs, int32 numOuts)
	{
		// we only support one in and output bus and these busses must have the same number of channels
		if (numIns == 1 &&
			numOuts == 1 &&
			Vst::SpeakerArr::getChannelCount(inputs[0]) == 2 &&
			Vst::SpeakerArr::getChannelCount(outputs[0]) == 2
			)	return AudioEffect::setBusArrangements(inputs, numIns, outputs, numOuts);

		return kResultFalse;
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Processor::terminate()
	{
		// Here the Plug-in will be de-instantiated, last possibility to remove some memory!
		for (int c = 0; c < 2; c++) {
			free(out_0[c]);
			free(out_1[c]);
			free(out_2[c]);
		}
		free(out_0);
		free(out_1);
		free(out_2);
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

		//---get audio buffers----------------
		uint32 sampleFramesSize = getSampleFramesSizeInBytes(processSetup, data.numSamples);
		void** in = getChannelBuffersPointer(processSetup, data.inputs[0]);
		void** out = getChannelBuffersPointer(processSetup, data.outputs[0]);
		Vst::SampleRate getSampleRate = processSetup.sampleRate;


		//---check if silence---------------
		if (data.inputs[0].silenceFlags != 0) // if flags is not zero => then it means that we have silent!
		{
			// As we know that the input is only filled with zero, the output will be then filled with zero too!

			data.outputs[0].silenceFlags = 0;

			if (data.inputs[0].silenceFlags & (uint64)1) { // Left
				if (in[0] != out[0]) memset(out[0], 0, sampleFramesSize);
				data.outputs[0].silenceFlags |= (uint64)1 << (uint64)0;
			}

			if (data.inputs[0].silenceFlags & (uint64)2) { // Right
				if (in[1] != out[1]) memset(out[1], 0, sampleFramesSize);
				data.outputs[0].silenceFlags |= (uint64)1 << (uint64)1;
			}

			if (data.inputs[0].silenceFlags & (uint64)3) {
				return kResultOk;
			}
		}

		data.outputs[0].silenceFlags = data.inputs[0].silenceFlags;

		//---in bypass mode outputs should be like inputs-----
		// handled with latency

		if (data.symbolicSampleSize == Vst::kSample32) {
			overSampling<Vst::Sample32>((Vst::Sample32**)in, (Vst::Sample32**)out, getSampleRate, data.numSamples);
		}
		else {
			overSampling<Vst::Sample64>((Vst::Sample64**)in, (Vst::Sample64**)out, getSampleRate, data.numSamples);
		}

		return kResultOk;
	}

	uint32 PLUGIN_API Sky_Blue_EQ4Processor::getLatencySamples()
	{
		if (fParamOS == overSample_1x) return 0;
		else if (fParamOS == overSample_2x) return 12;
		else return 24;
	}

	tresult PLUGIN_API Sky_Blue_EQ4Processor::setupProcessing(Vst::ProcessSetup& newSetup)
	{
		if (maxSample < newSetup.maxSamplesPerBlock) {
			for (int c = 0; c < 2; c++) {
				out_0[c] = (Vst::Sample64*)realloc(out_0[c], sizeof(Vst::Sample64) * newSetup.maxSamplesPerBlock);
				out_1[c] = (Vst::Sample64*)realloc(out_1[c], sizeof(Vst::Sample64) * newSetup.maxSamplesPerBlock * 2);
				out_2[c] = (Vst::Sample64*)realloc(out_2[c], sizeof(Vst::Sample64) * newSetup.maxSamplesPerBlock * 4);
				if ((out_0[c] == NULL) || (out_1[c] == NULL) || (out_2[c] == NULL) ) return kResultFalse;
				memset(out_0[c], 0, newSetup.maxSamplesPerBlock);
				memset(out_1[c], 0, newSetup.maxSamplesPerBlock * 2);
				memset(out_2[c], 0, newSetup.maxSamplesPerBlock * 4);
			}
		}

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

		double Fc_sub = 10.0;
		double Fc_40 = 39.0;
		double Fc_160 = 155.0;
		double Fc_650 = 625.0;
		double Fc_2k5 = 12.5;

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
#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* pi             */
#endif

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
	void Sky_Blue_EQ4Processor::processSVF(SampleType** inputs, Vst::Sample64** outputs, Steinberg::Vst::SampleRate getSampleRate, Steinberg::int32 sampleFrames)
	{

		SampleType* In_L = (SampleType*)inputs[0];
		SampleType* In_R = (SampleType*)inputs[1];
		Vst::Sample64* Out_L = (Vst::Sample64*)outputs[0];
		Vst::Sample64* Out_R = (Vst::Sample64*)outputs[1];

		Vst::Sample64 _db = (12.0 * fParamInputAtten - 12.0) - (3.0);

		Vst::Sample64 In_Atten = exp(log(10.0) * _db / 20.0);

		coefSVF(getSampleRate);

		int32 nParam_Sub = FromNormalized<Vst::ParamValue>(fParamSub, knob_stepCount);
		int32 nParam_40 = FromNormalized<Vst::ParamValue>(fParam40, knob_stepCount);
		int32 nParam_160 = FromNormalized<Vst::ParamValue>(fParam160, knob_stepCount);
		int32 nParam_650 = FromNormalized<Vst::ParamValue>(fParam650, knob_stepCount);
		int32 nParam_2k5 = FromNormalized<Vst::ParamValue>(fParam2k5, knob_stepCount);
		int32 nParam_Sky = FromNormalized<Vst::ParamValue>(fParamSky, knob_stepCount);

		BP_gain_Sub = BP_arr[nParam_Sub];
		BP_gain_40 = BP_arr[nParam_40];
		BP_gain_160 = BP_arr[nParam_160];
		BP_gain_650 = BP_arr[nParam_650];
		HP_gain_2k5 = HP_arr[nParam_2k5];
		HP_gain_Sky = HP_Sky_arr[nParam_Sky];
		if (fParamSkyFreq == 0.0) HP_gain_Sky = 0.0;

		int32 loop = 1;
		if (fParamOS == overSample_1x) loop = 1;
		else if (fParamOS == overSample_2x) loop = 2;
		else loop = 4;

		while (--sampleFrames >= 0)
		{
			Vst::Sample64 inputSampleL = *In_L;
			Vst::Sample64 inputSampleR = *In_R;
			In_L++;
			In_R++;

			inputSampleL *= In_Atten;
			inputSampleR *= In_Atten;

			for (int i = 0; i < loop; i++)
			{
				if (i > 0) {
					inputSampleL = 0.0;
					inputSampleR = 0.0;
				}

				//inputSampleL = computeSVF(&svf_LC[0], inputSampleL);
				//inputSampleR = computeSVF(&svf_LC[1], inputSampleR);

				Vst::Sample64 dataOutL = 0.0;
				Vst::Sample64 dataOutR = 0.0;

				svf_Sub[0].vin = inputSampleL;
				svf_Sub[0].t0 = svf_Sub[0].vin - svf_Sub[0].ic2eq;
				svf_Sub[0].v0 = svf_Sub[0].gt0 * svf_Sub[0].t0 - svf_Sub[0].gk0 * svf_Sub[0].ic1eq;
				svf_Sub[0].t1 = svf_Sub[0].gt1 * svf_Sub[0].t0 - svf_Sub[0].gk1 * svf_Sub[0].ic1eq;
				svf_Sub[0].t2 = svf_Sub[0].gt2 * svf_Sub[0].t0 + svf_Sub[0].gt1 * svf_Sub[0].ic1eq;
				svf_Sub[0].v1 = svf_Sub[0].t1 + svf_Sub[0].ic1eq;
				svf_Sub[0].v2 = svf_Sub[0].t2 + svf_Sub[0].ic2eq;
				svf_Sub[0].ic1eq += 2.0 * svf_Sub[0].t1;
				svf_Sub[0].ic2eq += 2.0 * svf_Sub[0].t2;
				Vst::Sample64 tmp_sub_L = svf_Sub[0].m0 * svf_Sub[0].v0 + svf_Sub[0].m1 * svf_Sub[0].v1 + svf_Sub[0].m2 * svf_Sub[0].v2;

				svf_40[0].vin = inputSampleL;
				svf_40[0].t0 = svf_40[0].vin - svf_40[0].ic2eq;
				svf_40[0].v0 = svf_40[0].gt0 * svf_40[0].t0 - svf_40[0].gk0 * svf_40[0].ic1eq;
				svf_40[0].t1 = svf_40[0].gt1 * svf_40[0].t0 - svf_40[0].gk1 * svf_40[0].ic1eq;
				svf_40[0].t2 = svf_40[0].gt2 * svf_40[0].t0 + svf_40[0].gt1 * svf_40[0].ic1eq;
				svf_40[0].v1 = svf_40[0].t1 + svf_40[0].ic1eq;
				svf_40[0].v2 = svf_40[0].t2 + svf_40[0].ic2eq;
				svf_40[0].ic1eq += 2.0 * svf_40[0].t1;
				svf_40[0].ic2eq += 2.0 * svf_40[0].t2;
				Vst::Sample64 tmp_40_L = svf_40[0].m0 * svf_40[0].v0 + svf_40[0].m1 * svf_40[0].v1 + svf_40[0].m2 * svf_40[0].v2;

				svf_160[0].vin = inputSampleL;
				svf_160[0].t0 = svf_160[0].vin - svf_160[0].ic2eq;
				svf_160[0].v0 = svf_160[0].gt0 * svf_160[0].t0 - svf_160[0].gk0 * svf_160[0].ic1eq;
				svf_160[0].t1 = svf_160[0].gt1 * svf_160[0].t0 - svf_160[0].gk1 * svf_160[0].ic1eq;
				svf_160[0].t2 = svf_160[0].gt2 * svf_160[0].t0 + svf_160[0].gt1 * svf_160[0].ic1eq;
				svf_160[0].v1 = svf_160[0].t1 + svf_160[0].ic1eq;
				svf_160[0].v2 = svf_160[0].t2 + svf_160[0].ic2eq;
				svf_160[0].ic1eq += 2.0 * svf_160[0].t1;
				svf_160[0].ic2eq += 2.0 * svf_160[0].t2;
				Vst::Sample64 tmp_160_L = svf_160[0].m0 * svf_160[0].v0 + svf_160[0].m1 * svf_160[0].v1 + svf_160[0].m2 * svf_160[0].v2;

				svf_650[0].vin = inputSampleL;
				svf_650[0].t0 = svf_650[0].vin - svf_650[0].ic2eq;
				svf_650[0].v0 = svf_650[0].gt0 * svf_650[0].t0 - svf_650[0].gk0 * svf_650[0].ic1eq;
				svf_650[0].t1 = svf_650[0].gt1 * svf_650[0].t0 - svf_650[0].gk1 * svf_650[0].ic1eq;
				svf_650[0].t2 = svf_650[0].gt2 * svf_650[0].t0 + svf_650[0].gt1 * svf_650[0].ic1eq;
				svf_650[0].v1 = svf_650[0].t1 + svf_650[0].ic1eq;
				svf_650[0].v2 = svf_650[0].t2 + svf_650[0].ic2eq;
				svf_650[0].ic1eq += 2.0 * svf_650[0].t1;
				svf_650[0].ic2eq += 2.0 * svf_650[0].t2;
				Vst::Sample64 tmp_650_L = svf_650[0].m0 * svf_650[0].v0 + svf_650[0].m1 * svf_650[0].v1 + svf_650[0].m2 * svf_650[0].v2;

				svf_2k5[0].vin = inputSampleL;
				svf_2k5[0].t0 = svf_2k5[0].vin - svf_2k5[0].ic2eq;
				svf_2k5[0].v0 = svf_2k5[0].gt0 * svf_2k5[0].t0 - svf_2k5[0].gk0 * svf_2k5[0].ic1eq;
				svf_2k5[0].t1 = svf_2k5[0].gt1 * svf_2k5[0].t0 - svf_2k5[0].gk1 * svf_2k5[0].ic1eq;
				svf_2k5[0].t2 = svf_2k5[0].gt2 * svf_2k5[0].t0 + svf_2k5[0].gt1 * svf_2k5[0].ic1eq;
				svf_2k5[0].v1 = svf_2k5[0].t1 + svf_2k5[0].ic1eq;
				svf_2k5[0].v2 = svf_2k5[0].t2 + svf_2k5[0].ic2eq;
				svf_2k5[0].ic1eq += 2.0 * svf_2k5[0].t1;
				svf_2k5[0].ic2eq += 2.0 * svf_2k5[0].t2;
				Vst::Sample64 tmp_2k5_L = svf_2k5[0].m0 * svf_2k5[0].v0 + svf_2k5[0].m1 * svf_2k5[0].v1 + svf_2k5[0].m2 * svf_2k5[0].v2;

				svf_Sky[0].vin = inputSampleL;
				svf_Sky[0].t0 = svf_Sky[0].vin - svf_Sky[0].ic2eq;
				svf_Sky[0].v0 = svf_Sky[0].gt0 * svf_Sky[0].t0 - svf_Sky[0].gk0 * svf_Sky[0].ic1eq;
				svf_Sky[0].t1 = svf_Sky[0].gt1 * svf_Sky[0].t0 - svf_Sky[0].gk1 * svf_Sky[0].ic1eq;
				svf_Sky[0].t2 = svf_Sky[0].gt2 * svf_Sky[0].t0 + svf_Sky[0].gt1 * svf_Sky[0].ic1eq;
				svf_Sky[0].v1 = svf_Sky[0].t1 + svf_Sky[0].ic1eq;
				svf_Sky[0].v2 = svf_Sky[0].t2 + svf_Sky[0].ic2eq;
				svf_Sky[0].ic1eq += 2.0 * svf_Sky[0].t1;
				svf_Sky[0].ic2eq += 2.0 * svf_Sky[0].t2;
				Vst::Sample64 tmp_sky_L = svf_Sky[0].m0 * svf_Sky[0].v0 + svf_Sky[0].m1 * svf_Sky[0].v1 + svf_Sky[0].m2 * svf_Sky[0].v2;



				svf_Sub[1].vin = inputSampleR;
				svf_Sub[1].t0 = svf_Sub[1].vin - svf_Sub[1].ic2eq;
				svf_Sub[1].v0 = svf_Sub[1].gt0 * svf_Sub[1].t0 - svf_Sub[1].gk0 * svf_Sub[1].ic1eq;
				svf_Sub[1].t1 = svf_Sub[1].gt1 * svf_Sub[1].t0 - svf_Sub[1].gk1 * svf_Sub[1].ic1eq;
				svf_Sub[1].t2 = svf_Sub[1].gt2 * svf_Sub[1].t0 + svf_Sub[1].gt1 * svf_Sub[1].ic1eq;
				svf_Sub[1].v1 = svf_Sub[1].t1 + svf_Sub[1].ic1eq;
				svf_Sub[1].v2 = svf_Sub[1].t2 + svf_Sub[1].ic2eq;
				svf_Sub[1].ic1eq += 2.0 * svf_Sub[1].t1;
				svf_Sub[1].ic2eq += 2.0 * svf_Sub[1].t2;
				Vst::Sample64 tmp_sub_R = svf_Sub[1].m0 * svf_Sub[1].v0 + svf_Sub[1].m1 * svf_Sub[1].v1 + svf_Sub[1].m2 * svf_Sub[1].v2;

				svf_40[1].vin = inputSampleR;
				svf_40[1].t0 = svf_40[1].vin - svf_40[1].ic2eq;
				svf_40[1].v0 = svf_40[1].gt0 * svf_40[1].t0 - svf_40[1].gk0 * svf_40[1].ic1eq;
				svf_40[1].t1 = svf_40[1].gt1 * svf_40[1].t0 - svf_40[1].gk1 * svf_40[1].ic1eq;
				svf_40[1].t2 = svf_40[1].gt2 * svf_40[1].t0 + svf_40[1].gt1 * svf_40[1].ic1eq;
				svf_40[1].v1 = svf_40[1].t1 + svf_40[1].ic1eq;
				svf_40[1].v2 = svf_40[1].t2 + svf_40[1].ic2eq;
				svf_40[1].ic1eq += 2.0 * svf_40[1].t1;
				svf_40[1].ic2eq += 2.0 * svf_40[1].t2;
				Vst::Sample64 tmp_40_R = svf_40[1].m0 * svf_40[1].v0 + svf_40[1].m1 * svf_40[1].v1 + svf_40[1].m2 * svf_40[1].v2;

				svf_160[1].vin = inputSampleR;
				svf_160[1].t0 = svf_160[1].vin - svf_160[1].ic2eq;
				svf_160[1].v0 = svf_160[1].gt0 * svf_160[1].t0 - svf_160[1].gk0 * svf_160[1].ic1eq;
				svf_160[1].t1 = svf_160[1].gt1 * svf_160[1].t0 - svf_160[1].gk1 * svf_160[1].ic1eq;
				svf_160[1].t2 = svf_160[1].gt2 * svf_160[1].t0 + svf_160[1].gt1 * svf_160[1].ic1eq;
				svf_160[1].v1 = svf_160[1].t1 + svf_160[1].ic1eq;
				svf_160[1].v2 = svf_160[1].t2 + svf_160[1].ic2eq;
				svf_160[1].ic1eq += 2.0 * svf_160[1].t1;
				svf_160[1].ic2eq += 2.0 * svf_160[1].t2;
				Vst::Sample64 tmp_160_R = svf_160[1].m0 * svf_160[1].v0 + svf_160[1].m1 * svf_160[1].v1 + svf_160[1].m2 * svf_160[1].v2;

				svf_650[1].vin = inputSampleR;
				svf_650[1].t0 = svf_650[1].vin - svf_650[1].ic2eq;
				svf_650[1].v0 = svf_650[1].gt0 * svf_650[1].t0 - svf_650[1].gk0 * svf_650[1].ic1eq;
				svf_650[1].t1 = svf_650[1].gt1 * svf_650[1].t0 - svf_650[1].gk1 * svf_650[1].ic1eq;
				svf_650[1].t2 = svf_650[1].gt2 * svf_650[1].t0 + svf_650[1].gt1 * svf_650[1].ic1eq;
				svf_650[1].v1 = svf_650[1].t1 + svf_650[1].ic1eq;
				svf_650[1].v2 = svf_650[1].t2 + svf_650[1].ic2eq;
				svf_650[1].ic1eq += 2.0 * svf_650[1].t1;
				svf_650[1].ic2eq += 2.0 * svf_650[1].t2;
				Vst::Sample64 tmp_650_R = svf_650[1].m0 * svf_650[1].v0 + svf_650[1].m1 * svf_650[1].v1 + svf_650[1].m2 * svf_650[1].v2;

				svf_2k5[1].vin = inputSampleR;
				svf_2k5[1].t0 = svf_2k5[1].vin - svf_2k5[1].ic2eq;
				svf_2k5[1].v0 = svf_2k5[1].gt0 * svf_2k5[1].t0 - svf_2k5[1].gk0 * svf_2k5[1].ic1eq;
				svf_2k5[1].t1 = svf_2k5[1].gt1 * svf_2k5[1].t0 - svf_2k5[1].gk1 * svf_2k5[1].ic1eq;
				svf_2k5[1].t2 = svf_2k5[1].gt2 * svf_2k5[1].t0 + svf_2k5[1].gt1 * svf_2k5[1].ic1eq;
				svf_2k5[1].v1 = svf_2k5[1].t1 + svf_2k5[1].ic1eq;
				svf_2k5[1].v2 = svf_2k5[1].t2 + svf_2k5[1].ic2eq;
				svf_2k5[1].ic1eq += 2.0 * svf_2k5[1].t1;
				svf_2k5[1].ic2eq += 2.0 * svf_2k5[1].t2;
				Vst::Sample64 tmp_2k5_R = svf_2k5[1].m0 * svf_2k5[1].v0 + svf_2k5[1].m1 * svf_2k5[1].v1 + svf_2k5[1].m2 * svf_2k5[1].v2;

				svf_Sky[1].vin = inputSampleR;
				svf_Sky[1].t0 = svf_Sky[1].vin - svf_Sky[1].ic2eq;
				svf_Sky[1].v0 = svf_Sky[1].gt0 * svf_Sky[1].t0 - svf_Sky[1].gk0 * svf_Sky[1].ic1eq;
				svf_Sky[1].t1 = svf_Sky[1].gt1 * svf_Sky[1].t0 - svf_Sky[1].gk1 * svf_Sky[1].ic1eq;
				svf_Sky[1].t2 = svf_Sky[1].gt2 * svf_Sky[1].t0 + svf_Sky[1].gt1 * svf_Sky[1].ic1eq;
				svf_Sky[1].v1 = svf_Sky[1].t1 + svf_Sky[1].ic1eq;
				svf_Sky[1].v2 = svf_Sky[1].t2 + svf_Sky[1].ic2eq;
				svf_Sky[1].ic1eq += 2.0 * svf_Sky[1].t1;
				svf_Sky[1].ic2eq += 2.0 * svf_Sky[1].t2;
				Vst::Sample64 tmp_sky_R = svf_Sky[1].m0 * svf_Sky[1].v0 + svf_Sky[1].m1 * svf_Sky[1].v1 + svf_Sky[1].m2 * svf_Sky[1].v2;


				dataOutL += tmp_sub_L * BP_gain_Sub;
				dataOutL += tmp_40_L * BP_gain_40;
				dataOutL += tmp_160_L * BP_gain_160;
				dataOutL += tmp_650_L * BP_gain_650;
				dataOutL += tmp_2k5_L * HP_gain_2k5;
				dataOutL += tmp_sky_L * HP_gain_Sky;

				dataOutR += tmp_sub_R * BP_gain_Sub;
				dataOutR += tmp_40_R * BP_gain_40;
				dataOutR += tmp_160_R * BP_gain_160;
				dataOutR += tmp_650_R * BP_gain_650;
				dataOutR += tmp_2k5_R * HP_gain_2k5;
				dataOutR += tmp_sky_R * HP_gain_Sky;

				inputSampleL = dataOutL; // * gain;
				inputSampleR = dataOutR; // * gain;

				*Out_L = inputSampleL;
				*Out_R = inputSampleR;

				Out_L++;
				Out_R++;
			}
		}
		return;
	}



	template <typename SampleType>
	void Sky_Blue_EQ4Processor::bypass_latency(SampleType** inputs, SampleType** outputs, Steinberg::int32 sampleFrames)
	{
		SampleType* In_L = (SampleType*)inputs[0];
		SampleType* In_R = (SampleType*)inputs[1];
		SampleType* Out_L = (SampleType*)outputs[0];
		SampleType* Out_R = (SampleType*)outputs[1];

		int32 latency = 0;
		if (fParamOS == overSample_1x) {
			memcpy(outputs[0], inputs[0], sizeof(SampleType) * sampleFrames);
			memcpy(outputs[1], inputs[1], sizeof(SampleType) * sampleFrames);
			return;
		}
		else if (fParamOS == overSample_2x) latency = 12;
		else latency = 24;

		while (--sampleFrames >= 0) {
			SampleType inputSampleL = *In_L;
			SampleType inputSampleR = *In_R;
			In_L++;
			In_R++;
			*Out_L = (SampleType)delay_buff[0][latency - 1];
			*Out_R = (SampleType)delay_buff[1][latency - 1];
			for (int i = latency - 1; i > 0; i--) {
				delay_buff[0][i] = delay_buff[0][i - 1];
				delay_buff[1][i] = delay_buff[1][i - 1];
			}
			delay_buff[0][0] = inputSampleL;
			delay_buff[1][0] = inputSampleR;
			Out_L++;
			Out_R++;
		}
		return;
	}

	template <typename SampleType>
	void Sky_Blue_EQ4Processor::proc_out(Vst::Sample64** inputs, SampleType** outputs, int32 sampleFrames)
	{
		Vst::Sample64* In_L = inputs[0];
		Vst::Sample64* In_R = inputs[1];
		SampleType* Out_L = (SampleType*)outputs[0];
		SampleType* Out_R = (SampleType*)outputs[1];

		while (--sampleFrames >= 0) {
			Vst::Sample64 inputSampleL = *In_L;
			Vst::Sample64 inputSampleR = *In_R;

			*Out_L = (SampleType)inputSampleL;
			*Out_R = (SampleType)inputSampleR;

			In_L++;
			In_R++;
			Out_L++;
			Out_R++;
		}
		return;
	}

	void Sky_Blue_EQ4Processor::FIR_dn_2_1(Vst::Sample64** inputs, Vst::Sample64** outputs, long sampleFrames) {
		Vst::Sample64* In_L = inputs[0];
		Vst::Sample64* In_R = inputs[1];
		Vst::Sample64* Out_L = outputs[0];
		Vst::Sample64* Out_R = outputs[1];

		int sf = sampleFrames;

		while (--sf >= 0) {

			Vst::Sample64 inputSampleL = *In_L;
			Vst::Sample64 inputSampleR = *In_R;


			double L = IR_coef_Dn_2[0] * inputSampleL;
			double R = IR_coef_Dn_2[0] * inputSampleR;
			for (int iter = 1, index = 0; iter < Num_coef_dn_2; iter += 1) {
				L += IR_coef_Dn_2[iter] * IR_buff_Dn_2[0][index];
				R += IR_coef_Dn_2[iter] * IR_buff_Dn_2[1][index];
				index++;
			}
			*Out_L = L;
			*Out_R = R;
			Out_L++;
			Out_R++;

			for (int k = Num_buff_dn_2; k > 2; k--) IR_buff_Dn_2[0][k - 1] = IR_buff_Dn_2[0][k - 3];
			for (int k = Num_buff_dn_2; k > 2; k--) IR_buff_Dn_2[1][k - 1] = IR_buff_Dn_2[1][k - 3];

			IR_buff_Dn_2[0][1] = inputSampleL;
			IR_buff_Dn_2[1][1] = inputSampleR;
			In_L++;
			In_R++;

			inputSampleL = *In_L;
			inputSampleR = *In_R;
			IR_buff_Dn_2[0][0] = inputSampleL;
			IR_buff_Dn_2[1][0] = inputSampleR;
			In_L++;
			In_R++;
		}
		return;
	}

	void Sky_Blue_EQ4Processor::FIR_dn_4_1(Vst::Sample64** inputs, Vst::Sample64** outputs, long sampleFrames) {
		Vst::Sample64* In_L = inputs[0];
		Vst::Sample64* In_R = inputs[1];
		Vst::Sample64* Out_L = outputs[0];
		Vst::Sample64* Out_R = outputs[1];

		int sf = sampleFrames;

		while (--sf >= 0) {

			Vst::Sample64 inputSampleL = *In_L;
			Vst::Sample64 inputSampleR = *In_R;


			double L = IR_coef_Dn_22[0] * inputSampleL;
			double R = IR_coef_Dn_22[0] * inputSampleR;
			for (int iter = 1, index = 0; iter < Num_coef_dn_22; iter += 1) {
				L += IR_coef_Dn_22[iter] * IR_buff_Dn_22[0][index];
				R += IR_coef_Dn_22[iter] * IR_buff_Dn_22[1][index];
				index++;
			}
			*Out_L = L;
			*Out_R = R;
			Out_L++;
			Out_R++;

			for (int k = Num_buff_dn_22; k > 4; k--) IR_buff_Dn_22[0][k - 1] = IR_buff_Dn_22[0][k - 5];
			for (int k = Num_buff_dn_22; k > 4; k--) IR_buff_Dn_22[1][k - 1] = IR_buff_Dn_22[1][k - 5];


			IR_buff_Dn_22[0][3] = inputSampleL;
			IR_buff_Dn_22[1][3] = inputSampleR;
			In_L++;
			In_R++;

			inputSampleL = *In_L;
			inputSampleR = *In_R;
			IR_buff_Dn_22[0][2] = inputSampleL;
			IR_buff_Dn_22[1][2] = inputSampleR;
			In_L++;
			In_R++;

			inputSampleL = *In_L;
			inputSampleR = *In_R;
			IR_buff_Dn_22[0][1] = inputSampleL;
			IR_buff_Dn_22[1][1] = inputSampleR;
			In_L++;
			In_R++;

			inputSampleL = *In_L;
			inputSampleR = *In_R;
			IR_buff_Dn_22[0][0] = inputSampleL;
			IR_buff_Dn_22[1][0] = inputSampleR;
			In_L++;
			In_R++;
		}
		return;
	}

	template <typename SampleType>
	void Sky_Blue_EQ4Processor::overSampling(
		SampleType** inputs,
		SampleType** outputs,
		Vst::SampleRate getSampleRate,
		int32 sampleFrames
	)
	{
		Vst::SampleRate SR_1 = getSampleRate;
		Vst::SampleRate SR_2 = getSampleRate * 2.0;
		Vst::SampleRate SR_4 = getSampleRate * 4.0;

		// memcpy(out_2[0], in_2[0], sizeof(SampleType) * len_4);
		// memcpy(out_2[1], in_2[1], sizeof(SampleType) * len_4);

		if (fParamOS == overSample_1x) {
			processSVF<SampleType>(inputs, out_0, SR_1, sampleFrames);
			if (bParamBypass) {
				bypass_latency<SampleType>(inputs, outputs, sampleFrames);
				return;
			}
		}
		else if (fParamOS == overSample_2x) {
			processSVF<SampleType>(inputs, out_1, SR_2, sampleFrames);
			if (bParamBypass) {
				bypass_latency<SampleType>(inputs, outputs, sampleFrames);
				return;
			}
			FIR_dn_2_1(out_1, out_0, sampleFrames);
		}
		else {
			processSVF<SampleType>(inputs, out_2, SR_4, sampleFrames);
			if (bParamBypass) {
				bypass_latency<SampleType>(inputs, outputs, sampleFrames);
				return;
			}
			FIR_dn_4_1(out_2, out_0, sampleFrames);
		}

		proc_out<SampleType>(out_0, (SampleType**)outputs, sampleFrames);

		return;
	}


	//------------------------------------------------------------------------
} // namespace yg331





/*

void Sky_Blue_EQ4Processor::FIR_up_x2(Vst::Sample64** inputs, Vst::Sample64** outputs, long sampleFrames) {
		Vst::Sample64* In_L = inputs[0];
		Vst::Sample64* In_R = inputs[1];
		Vst::Sample64* Out_L = outputs[0];
		Vst::Sample64* Out_R = outputs[1];

		int sf = sampleFrames;


		while (--sf >= 0) {
			Vst::Sample64 inputSampleL = *In_L;
			Vst::Sample64 inputSampleR = *In_R;
			In_L++;
			In_R++;

			for (int k = Num_buff_up; k > 1; k--) IR_buff_Up[0][k - 1] = IR_buff_Up[0][k - 2];
			for (int k = Num_buff_up; k > 1; k--) IR_buff_Up[1][k - 1] = IR_buff_Up[1][k - 2];
			IR_buff_Up[0][0] = inputSampleL;
			IR_buff_Up[1][0] = inputSampleR;

			for (int base = 0; base < 2; base++) {
				int index = 0;
				double L = 0.0, R = 0.0;
				for (int iter = base; iter < Num_coef_up; iter += 2) {
					L += IR_coef_Up[iter] * IR_buff_Up[0][index];
					R += IR_coef_Up[iter] * IR_buff_Up[1][index];
					index++;
				}
				*Out_L = L;
				*Out_R = R;
				Out_L++;
				Out_R++;
			}
		}
		return;
	}

	void Sky_Blue_EQ4Processor::FIR_up_x4(Vst::Sample64** inputs, Vst::Sample64** outputs, long sampleFrames) {
		Vst::Sample64* In_L = inputs[0];
		Vst::Sample64* In_R = inputs[1];
		Vst::Sample64* Out_L = outputs[0];
		Vst::Sample64* Out_R = outputs[1];

		int sf = sampleFrames;


		while (--sf >= 0) {
			Vst::Sample64 inputSampleL = *In_L;
			Vst::Sample64 inputSampleR = *In_R;
			In_L++;
			In_R++;

			for (int k = Num_buff_up; k > 1; k--) IR_buff_Up[0][k - 1] = IR_buff_Up[0][k - 2];
			for (int k = Num_buff_up; k > 1; k--) IR_buff_Up[1][k - 1] = IR_buff_Up[1][k - 2];
			IR_buff_Up[0][0] = inputSampleL;
			IR_buff_Up[1][0] = inputSampleR;



			inpt : [ s ]


			bffr : [ -0 -1 -2 -3 -4 ]
			bffr : [  0  1  2  3  4 ]
			vtrl : [ 4 0 0 0 3 0 0 0 2 0 0 0 1 0 0 0 s 0 0 0 ]

			coef : 19 tap
0.c18 + 0.c17 + 4.c16 + 0.c15 + 0.c14 + 0.c13 + 3.c12 + 0.c11 + 0.c10 + 0.c09 + 2.c08 + 0.c07 + 0.c06 + 0.c05 + 1.c04 + 0.c03 + 0.c02 + 0.c01 + s.c00
0.c18 + 4.c17 + 0.c16 + 0.c15 + 0.c14 + 3.c13 + 0.c12 + 0.c11 + 0.c10 + 2.c09 + 0.c08 + 0.c07 + 0.c06 + 1.c05 + 0.c04 + 0.c03 + 0.c02 + s.c01 + 0.c00
4.c18 + 0.c17 + 0.c16 + 0.c15 + 3.c14 + 0.c13 + 0.c12 + 0.c11 + 2.c10 + 0.c09 + 0.c08 + 0.c07 + 1.c06 + 0.c05 + 0.c04 + 0.c03 + s.c02 + 0.c01 + 0.c00
0.c18 + 0.c17 + 0.c16 + 3.c15 + 0.c14 + 0.c13 + 0.c12 + 2.c11 + 0.c10 + 0.c09 + 0.c08 + 1.c07 + 0.c06 + 0.c05 + 0.c04 + s.c03 + 0.c02 + 0.c01 + 0.c00


			a = 4.c16 + 3.c12 + 2.c08 + 1.c04 + s.c00
			b = 4.c17 + 3.c13 + 2.c09 + 1.c05 + s.c01
			c = 4.c18 + 3.c14 + 2.c10 + 1.c06 + s.c02
			d =         3.c15 + 2.c11 + 1.c07 + s.c03

			out  : [ a b c d ]


			for (int base = 0; base < 4; base++) {
				int index = 0;
				double L = 0.0, R = 0.0;
				for (int iter = base; iter < Num_coef_up; iter += 4) {
					L += IR_coef_Up[iter] * IR_buff_Up[0][index];
					R += IR_coef_Up[iter] * IR_buff_Up[1][index];
					index++;
				}
				*Out_L = L;
				*Out_R = R;
				Out_L++;
				Out_R++;
			}
		}
		return;
	}
	

void Sky_Blue_EQ4Processor::FIR_dn_x4(Vst::Sample64** inputs, Vst::Sample64** outputs, long sampleFrames) {
	Vst::Sample64* In_L = inputs[0];
	Vst::Sample64* In_R = inputs[1];
	Vst::Sample64* Out_L = outputs[0];
	Vst::Sample64* Out_R = outputs[1];

	int sf = sampleFrames;

	

	input : [ n 1 2 3 ]

	bffr : [ -4 -3 -2 -1 ]

	tap : 5
	c4 c3 c2 c1 c0

	out : -4.c4 + -3.c3 + -2.c2 + -1.c1 + n.c0


	

	while (--sf >= 0) {

		Vst::Sample64 inputSampleL = *In_L;
		Vst::Sample64 inputSampleR = *In_R;


		double L = IR_coef_Dn[0] * inputSampleL;
		double R = IR_coef_Dn[0] * inputSampleR;
		for (int iter = 1, index = 0; iter < Num_coef_dn; iter += 1) {
			L += IR_coef_Dn[iter] * IR_buff_Dn[0][index];
			R += IR_coef_Dn[iter] * IR_buff_Dn[1][index];
			index++;
		}
		*Out_L = L;
		*Out_R = R;
		Out_L++;
		Out_R++;

		for (int k = Num_buff_dn; k > 4; k--) IR_buff_Dn[0][k - 1] = IR_buff_Dn[0][k - 5];
		for (int k = Num_buff_dn; k > 4; k--) IR_buff_Dn[1][k - 1] = IR_buff_Dn[1][k - 5];


		IR_buff_Dn[0][3] = inputSampleL;
		IR_buff_Dn[1][3] = inputSampleR;
		In_L++;
		In_R++;

		inputSampleL = *In_L;
		inputSampleR = *In_R;
		IR_buff_Dn[0][2] = inputSampleL;
		IR_buff_Dn[1][2] = inputSampleR;
		In_L++;
		In_R++;

		inputSampleL = *In_L;
		inputSampleR = *In_R;
		IR_buff_Dn[0][1] = inputSampleL;
		IR_buff_Dn[1][1] = inputSampleR;
		In_L++;
		In_R++;

		inputSampleL = *In_L;
		inputSampleR = *In_R;
		IR_buff_Dn[0][0] = inputSampleL;
		IR_buff_Dn[1][0] = inputSampleR;
		In_L++;
		In_R++;
	}
	return;
}



	inline void Sky_Blue_EQ4Processor::setCoeffs(double Fs)
	{
		double Fc = 0.0, K = 0.0, norm = 0.0;
		double Q = 0.49763809;
		double PI_div_Fs = M_PI / Fs;

		int32 nParam_airf = FromNormalized<Vst::ParamValue>(fParamAirFreq, 6);
		double Air_f[7] = {0.0, 2000.0, 4000.0, 7800.0, 12200.0, 17000.0, 32500.0};
		Fc = Air_f[nParam_airf];

		// 20000Hz
		// Fc = 2200.0;
		K = tan(PI_div_Fs * Fc);
		norm = 1.0 / (K + 1.0);
		z_Air[0] = norm;
		z_Air[1] = -norm;
		z_Air[2] = 0.0;
		p_Air[0] = 1.0;
		p_Air[1] = (K - 1) * norm;
		p_Air[2] = 0.0;

		// 2500Hz
		Fc = 1200.0; // YES
		K = tan(PI_div_Fs * Fc);
		norm = 1.0 / (K + 1.0);
		z_2k5[0] = norm;
		z_2k5[1] = -norm;
		z_2k5[2] = 0.0;
		p_2k5[0] = 1.0;
		p_2k5[1] = (K - 1) * norm;
		p_2k5[2] = 0.0;

		Fc = 620.0;
		K = tan(PI_div_Fs * Fc);
		norm = 1.0 / (1.0 + K / Q + K * K);
		z_650[0] = K / Q * norm;
		z_650[1] = 0.0;
		z_650[2] = -z_650[0];
		p_650[0] = 1.0;
		p_650[1] = 2.0 * (K * K - 1.0) * norm;
		p_650[2] = (1.0 - K / Q + K * K) * norm;

		Fc = 155.0;
		K = tan(PI_div_Fs * Fc);
		norm = 1.0 / (1.0 + K / Q + K * K);
		z_160[0] = K / Q * norm;
		z_160[1] = 0.0;
		z_160[2] = -z_160[0];
		p_160[0] = 1.0;
		p_160[1] = 2.0 * (K * K - 1.0) * norm;
		p_160[2] = (1.0 - K / Q + K * K) * norm;

		Fc = 39.0;
		K = tan(PI_div_Fs * Fc);
		norm = 1.0 / (1.0 + K / Q + K * K);
		z_40[0] = K / Q * norm;
		z_40[1] = 0.0;
		z_40[2] = -z_40[0];
		p_40[0] = 1.0;
		p_40[1] = 2.0 * (K * K - 1.0) * norm;
		p_40[2] = (1.0 - K / Q + K * K) * norm;

		Fc = 9.5;
		K = tan(PI_div_Fs * Fc);
		norm = 1.0 / (1.0 + K / Q + K * K);
		z_10[0] = K / Q * norm;
		z_10[1] = 0.0;
		z_10[2] = -z_10[0];
		p_10[0] = 1.0;
		p_10[1] = 2.0 * (K * K - 1.0) * norm;
		p_10[2] = (1.0 - K / Q + K * K) * norm;

		for (int i = 0; i < 3; i++) {
			z_10[i] *= 5.85;
			z_40[i] *= 5.58;
			z_160[i] *= 5.40;
			z_650[i] *= 4.99;
			z_2k5[i] *= 8.174126;
			z_Air[i] *= 7.943282;
		}
	};

template <typename SampleType>
	void Sky_Blue_EQ4Processor::processAudio(Vst::Sample64** inputs, Vst::Sample64** outputs, Steinberg::Vst::SampleRate getSampleRate, Steinberg::int32 sampleFrames)
	{
		Vst::Sample64* In_L = (Vst::Sample64*)inputs[0];
		Vst::Sample64* In_R = (Vst::Sample64*)inputs[1];
		Vst::Sample64* Out_L = (Vst::Sample64*)outputs[0];
		Vst::Sample64* Out_R = (Vst::Sample64*)outputs[1];

		Vst::Sample64 In_Atten = exp(log(10.0) * (12.0 * fParamInputAtten - 12.0) / 20.0);

		Vst::Sample64 g_10, pg_10;
		Vst::Sample64 g_40, pg_40;
		Vst::Sample64 g_160, pg_160;
		Vst::Sample64 g_650, pg_650;
		Vst::Sample64 g_2k5, pg_2k5;
		Vst::Sample64 g_Air, pg_Air;

		double g_arr[21] = {
				0.02, // -5.0
				0.022, // -4.5
				0.024, // -4.0
				0.028, // -3.5
				0.036, // -3.0
				0.045, // -2.5
				0.05, // -2.0
				0.08,// -1.5
				0.13,// -1.0
				0.19, // -0.5
				0.25, // 0.0
				0.315, // +0.5
				0.38,  // +1.0
				0.48,  // +1.5
				0.7, // +2.0
				1.2, // +2.5
				1.42, // +3.0
				1.8, // +3.5
				2.35, // +4.0
				3.25, // +4.5
				3.3  // +5.0
		};
		double pg_arr[21] = {
				2.1, //  0/20  -5.0
				1.9, //  1/20  -4.5
				1.75, //  2/20  -4.0
				1.62, //  3/20  -3.5
				1.45, //  4/20  -3.0
				1.45, //  5/20  -2.5
				1.45, //  6/20  -2.0
				1.2, //  7/20  -1.5
				1.1, //  8/20  -1.0
				1.04, //  9/20  -0.5
				1.0, // 10/20   0.0
				0.98, // +0.5
				0.98, // +1.0
				0.94, // +1.5
				0.92, // +2.0
				0.92, // +2.5
				0.91, // +3.0
				0.9, // +3.5
				0.9, // +4.0
				0.9, // +4.5
				0.9  // +5.0
		};

		int32 nParam_10 = FromNormalized<Vst::ParamValue>(fParam10, knob_stepCount);
		g_10 = g_arr[nParam_10];
		pg_10 = pg_arr[nParam_10];
		int32 nParam_40 = FromNormalized<Vst::ParamValue>(fParam40, knob_stepCount);
		g_40 = g_arr[nParam_40];
		pg_40 = pg_arr[nParam_40];
		int32 nParam_160 = FromNormalized<Vst::ParamValue>(fParam160, knob_stepCount);
		g_160 = g_arr[nParam_160];
		pg_160 = pg_arr[nParam_160];
		int32 nParam_650 = FromNormalized<Vst::ParamValue>(fParam650, knob_stepCount);
		g_650 = g_arr[nParam_650];
		pg_650 = pg_arr[nParam_650];
		int32 nParam_2k5 = FromNormalized<Vst::ParamValue>(fParam2k5, knob_stepCount);
		g_2k5 = g_arr[nParam_2k5];
		pg_2k5 = pg_arr[nParam_2k5];

		double g_air[21] = {
				0.037, // 0/20
				0.0375, // 1/20
				0.0380, // 2/20
				0.042, // 3/20
				0.050, // 4/20
				0.052, // 5/20
				0.058, // 6/20
				0.072, // 7/20
				0.108, // 8/20
				0.148, // 9/20
				0.20, // 10/20
				0.25, // 11/20
				0.28,  // 12/20
				0.34,  // 13/20
				0.48, // 14/20
				1.04, // 15/20
				1.22, // 16/20
				1.52, // 17/20
				2.00, // 18/20
				2.78, // 19/20
				2.82  // 20/20
		};
		double pg_air[21] = {
				1.14, // 0/20
				1.14, // 1/20
				1.13, // 2/20
				1.10, // 3/20
				1.05, // 4/20
				1.158, // 5/20
				1.28, // 6/20
				1.368, // 7/20
				1.41, // 8/20
				1.43, // 9/20
				1.28, // 10/20
				1.18, // 11/20
				1.3,  // 12/20
				1.32,  // 13/20
				1.28, // 14/20
				1.275, // 15/20
				1.272, // 16/20
				1.28, // 17/20
				1.268, // 18/20
				1.275, // 19/20
				1.28  // 20/20
		};
		int32 nParam_air = FromNormalized<Vst::ParamValue>(fParamAir, knob_stepCount);
		g_Air = g_air[nParam_air];
		pg_Air = pg_air[nParam_air];

		if (fParamAirFreq == 0.0)
			g_Air = 0.0;


		Vst::Sample64 dcGain = g_10 + g_40 + g_160 + g_650 + g_2k5 + g_Air;
		Vst::Sample64 globalGain = 0.398 / dcGain;

		while (--sampleFrames >= 0)
		{
			Vst::Sample64 inputSampleL = *In_L;
			Vst::Sample64 inputSampleR = *In_R;

			inputSampleL *= In_Atten;
			inputSampleR *= In_Atten;

			Vst::Sample64 dataOutL = 0.0;
			Vst::Sample64 dataOutR = 0.0;

			x_10[0][0] = inputSampleL;
			x_40[0][0] = inputSampleL;
			x_160[0][0] = inputSampleL;
			x_650[0][0] = inputSampleL;
			x_2k5[0][0] = inputSampleL;
			x_Air[0][0] = inputSampleL;

			x_10[1][0] = inputSampleR;
			x_40[1][0] = inputSampleR;
			x_160[1][0] = inputSampleR;
			x_650[1][0] = inputSampleR;
			x_2k5[1][0] = inputSampleR;
			x_Air[1][0] = inputSampleR;

			// 10Hz
			y_10[0][0] = x_10[0][0] * z_10[0]
				+ x_10[0][1] * z_10[1]
				+ x_10[0][2] * z_10[2]
				- y_10[0][1] * p_10[1]
				- y_10[0][2] * p_10[2];
			x_10[0][2] = x_10[0][1];
			x_10[0][1] = x_10[0][0];
			y_10[0][2] = y_10[0][1];
			y_10[0][1] = y_10[0][0];

			y_10[1][0] = x_10[1][0] * z_10[0]
				+ x_10[1][1] * z_10[1]
				+ x_10[1][2] * z_10[2]
				- y_10[1][1] * p_10[1]
				- y_10[1][2] * p_10[2];
			x_10[1][2] = x_10[1][1];
			x_10[1][1] = x_10[1][0];
			y_10[1][2] = y_10[1][1];
			y_10[1][1] = y_10[1][0];

			// 40Hz
			y_40[0][0] = x_40[0][0] * z_40[0]
				+ x_40[0][1] * z_40[1]
				+ x_40[0][2] * z_40[2]
				- y_40[0][1] * p_40[1]
				- y_40[0][2] * p_40[2];
			x_40[0][2] = x_40[0][1];
			x_40[0][1] = x_40[0][0];
			y_40[0][2] = y_40[0][1];
			y_40[0][1] = y_40[0][0];

			y_40[1][0] = x_40[1][0] * z_40[0]
				+ x_40[1][1] * z_40[1]
				+ x_40[1][2] * z_40[2]
				- y_40[1][1] * p_40[1]
				- y_40[1][2] * p_40[2];
			x_40[1][2] = x_40[1][1];
			x_40[1][1] = x_40[1][0];
			y_40[1][2] = y_40[1][1];
			y_40[1][1] = y_40[1][0];

			// 160Hz
			y_160[0][0] = x_160[0][0] * z_160[0]
				+ x_160[0][1] * z_160[1]
				+ x_160[0][2] * z_160[2]
				- y_160[0][1] * p_160[1]
				- y_160[0][2] * p_160[2];
			x_160[0][2] = x_160[0][1];
			x_160[0][1] = x_160[0][0];
			y_160[0][2] = y_160[0][1];
			y_160[0][1] = y_160[0][0];

			y_160[1][0] = x_160[1][0] * z_160[0]
				+ x_160[1][1] * z_160[1]
				+ x_160[1][2] * z_160[2]
				- y_160[1][1] * p_160[1]
				- y_160[1][2] * p_160[2];
			x_160[1][2] = x_160[1][1];
			x_160[1][1] = x_160[1][0];
			y_160[1][2] = y_160[1][1];
			y_160[1][1] = y_160[1][0];


			// 640Hz
			y_650[0][0] = x_650[0][0] * z_650[0]
				+ x_650[0][1] * z_650[1]
				+ x_650[0][2] * z_650[2]
				- y_650[0][1] * p_650[1]
				- y_650[0][2] * p_650[2];
			x_650[0][2] = x_650[0][1];
			x_650[0][1] = x_650[0][0];
			y_650[0][2] = y_650[0][1];
			y_650[0][1] = y_650[0][0];

			y_650[1][0] = x_650[1][0] * z_650[0]
				+ x_650[1][1] * z_650[1]
				+ x_650[1][2] * z_650[2]
				- y_650[1][1] * p_650[1]
				- y_650[1][2] * p_650[2];
			x_650[1][2] = x_650[1][1];
			x_650[1][1] = x_650[1][0];
			y_650[1][2] = y_650[1][1];
			y_650[1][1] = y_650[1][0];

			// 2500Hz
			y_2k5[0][0] = x_2k5[0][0] * z_2k5[0]
				+ x_2k5[0][1] * z_2k5[1]
				+ x_2k5[0][2] * z_2k5[2]
				- y_2k5[0][1] * p_2k5[1]
				- y_2k5[0][2] * p_2k5[2];
			x_2k5[0][2] = x_2k5[0][1];
			x_2k5[0][1] = x_2k5[0][0];
			y_2k5[0][2] = y_2k5[0][1];
			y_2k5[0][1] = y_2k5[0][0];

			y_2k5[1][0] = x_2k5[1][0] * z_2k5[0]
				+ x_2k5[1][1] * z_2k5[1]
				+ x_2k5[1][2] * z_2k5[2]
				- y_2k5[1][1] * p_2k5[1]
				- y_2k5[1][2] * p_2k5[2];
			x_2k5[1][2] = x_2k5[1][1];
			x_2k5[1][1] = x_2k5[1][0];
			y_2k5[1][2] = y_2k5[1][1];
			y_2k5[1][1] = y_2k5[1][0];


			// 20kHz
			y_Air[0][0] = x_Air[0][0] * z_Air[0]
				+ x_Air[0][1] * z_Air[1]
				+ x_Air[0][2] * z_Air[2]
				- y_Air[0][1] * p_Air[1]
				- y_Air[0][2] * p_Air[2];
			x_Air[0][2] = x_Air[0][1];
			x_Air[0][1] = x_Air[0][0];
			y_Air[0][2] = y_Air[0][1];
			y_Air[0][1] = y_Air[0][0];

			y_Air[1][0] = x_Air[1][0] * z_Air[0]
				+ x_Air[1][1] * z_Air[1]
				+ x_Air[1][2] * z_Air[2]
				- y_Air[1][1] * p_Air[1]
				- y_Air[1][2] * p_Air[2];
			x_Air[1][2] = x_Air[1][1];
			x_Air[1][1] = x_Air[1][0];
			y_Air[1][2] = y_Air[1][1];
			y_Air[1][1] = y_Air[1][0];

			dataOutL += (y_10[0][0] * pg_10 + inputSampleL) * g_10;
			dataOutL += (y_40[0][0] * pg_40 + inputSampleL) * g_40;
			dataOutL += (y_160[0][0] * pg_160 + inputSampleL) * g_160;
			dataOutL += (y_650[0][0] * pg_650 + inputSampleL) * g_650;
			dataOutL += (y_2k5[0][0] * pg_2k5 + inputSampleL) * g_2k5;
			dataOutL += (y_Air[0][0] * pg_Air + inputSampleL) * g_Air;

			dataOutR += (y_10[1][0] * pg_10 + inputSampleR) * g_10;
			dataOutR += (y_40[1][0] * pg_40 + inputSampleR) * g_40;
			dataOutR += (y_160[1][0] * pg_160 + inputSampleR) * g_160;
			dataOutR += (y_650[1][0] * pg_650 + inputSampleR) * g_650;
			dataOutR += (y_2k5[1][0] * pg_2k5 + inputSampleR) * g_2k5;
			dataOutR += (y_Air[1][0] * pg_Air + inputSampleR) * g_Air;


			if (bParamAuto) {
				inputSampleL = dataOutL * globalGain;
				inputSampleR = dataOutR * globalGain;
			}
			else {
				inputSampleL = dataOutL * 0.316228; // -10dB
				inputSampleR = dataOutR * 0.316228; // -10dB
			}

			*Out_L = inputSampleL;
			*Out_R = inputSampleR;

			In_L++;
			In_R++;
			Out_L++;
			Out_R++;
		}
		return;
	}


*/