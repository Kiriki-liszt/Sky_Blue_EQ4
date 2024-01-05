//------------------------------------------------------------------------
// Copyright(c) 2023 yg331.
//------------------------------------------------------------------------

#pragma once

#include "SBEQ4_cids.h"
#include "public.sdk/source/vst/vstaudioeffect.h"

#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <queue>

namespace yg331 {

	//------------------------------------------------------------------------
	//  Sky_Blue_EQ4Processor
	//------------------------------------------------------------------------
	class Sky_Blue_EQ4Processor : public Steinberg::Vst::AudioEffect
	{
	public:
		Sky_Blue_EQ4Processor();
		~Sky_Blue_EQ4Processor() SMTG_OVERRIDE;

		// Create function
		static Steinberg::FUnknown* createInstance(void* /*context*/)
		{
			return (Steinberg::Vst::IAudioProcessor*)new Sky_Blue_EQ4Processor;
		}

		//--- ---------------------------------------------------------------------
		// AudioEffect overrides:
		//--- ---------------------------------------------------------------------
		/** Called at first after constructor */
		Steinberg::tresult PLUGIN_API initialize(Steinberg::FUnknown* context) SMTG_OVERRIDE;

		/** Called at the end before destructor */
		Steinberg::tresult PLUGIN_API terminate() SMTG_OVERRIDE;

		/** Switch the Plug-in on/off */
		Steinberg::tresult PLUGIN_API setActive(Steinberg::TBool state) SMTG_OVERRIDE;

		/** Will be called before any process call */
		Steinberg::tresult PLUGIN_API setupProcessing(Steinberg::Vst::ProcessSetup& newSetup) SMTG_OVERRIDE;

		/** Asks if a given sample size is supported see SymbolicSampleSizes. */
		Steinberg::tresult PLUGIN_API canProcessSampleSize(Steinberg::int32 symbolicSampleSize) SMTG_OVERRIDE;

		/** Here we go...the process call */
		Steinberg::tresult PLUGIN_API process(Steinberg::Vst::ProcessData& data) SMTG_OVERRIDE;

		/** For persistence */
		Steinberg::tresult PLUGIN_API setState(Steinberg::IBStream* state) SMTG_OVERRIDE;
		Steinberg::tresult PLUGIN_API getState(Steinberg::IBStream* state) SMTG_OVERRIDE;

		//------------------------------------------------------------------------

		Steinberg::tresult PLUGIN_API setBusArrangements(
			Steinberg::Vst::SpeakerArrangement* inputs, Steinberg::int32 numIns,
			Steinberg::Vst::SpeakerArrangement* outputs, Steinberg::int32 numOuts
		) SMTG_OVERRIDE;

		/** Gets the current Latency in samples. */
		Steinberg::uint32 PLUGIN_API getLatencySamples() SMTG_OVERRIDE;

	protected:
		inline void coefSVF ( Steinberg::Vst::SampleRate Fs );

		inline void setSVF
		(
			SVF_* svf_filter,
			filter_type kFilter,
			Steinberg::Vst::ParamValue q,
			Steinberg::Vst::ParamValue gain,
			Steinberg::Vst::Sample64 Fc,
			Steinberg::Vst::Sample64 Fs
		);

		inline Steinberg::Vst::Sample64 computeSVF
		(
			SVF_* svf_filter,
			Steinberg::Vst::Sample64 input
		);

		template <typename SampleType>
		void processSVF
		(
			SampleType** inputs,
			SampleType** outputs,
			Steinberg::int32 numChannels,
			Steinberg::Vst::SampleRate getSampleRate,
			Steinberg::int32 sampleFrames
		);

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* pi             */
#endif
#define maxTap 256

		inline double Ino(double x)
		{
			double d = 0, ds = 1, s = 1;
			do
			{
				d += 2;
				ds *= x * x / (d * d);
				s += ds;
			} while (ds > s * 1e-6);
			return s;
		}

		void calcFilter(double Fs, double Fa, double Fb, int M, double Att, double* dest)
		{
			// Kaiser windowed FIR filter "DIGITAL SIGNAL PROCESSING, II" IEEE Press pp 123-126.

			int Np = (M - 1) / 2;
			double A[maxTap] = { 0, };
			double Alpha;
			double Inoalpha;
			//double H[maxTap] = { 0, };

			A[0] = 2 * (Fb - Fa) / Fs;

			for (int j = 1; j <= Np; j++)
				A[j] = (sin(2.0 * j * M_PI * Fb / Fs) - sin(2.0 * j * M_PI * Fa / Fs)) / (j * M_PI);

			if (Att < 21.0)
				Alpha = 0;
			else if (Att > 50.0)
				Alpha = 0.1102 * (Att - 8.7);
			else
				Alpha = 0.5842 * pow((Att - 21), 0.4) + 0.07886 * (Att - 21);

			Inoalpha = Ino(Alpha);

			for (int j = 0; j <= Np; j++)
				dest[Np + j] = A[j] * Ino(Alpha * sqrt(1.0 - ((double)(j * j) / (double)(Np * Np)))) / Inoalpha;

			for (int j = 0; j < Np; j++)
				dest[j] = dest[M - 1 - j];
		}

		// Plugin controls ------------------------------------------------------------------
		Steinberg::Vst::ParamValue fParamZoom;
		Steinberg::Vst::ParamValue fParamInputAtten;
		Steinberg::Vst::ParamValue fParamSkyFreq;
		Steinberg::Vst::ParamValue fParamSky;
		Steinberg::Vst::ParamValue fParam2k5;
		Steinberg::Vst::ParamValue fParam650;
		Steinberg::Vst::ParamValue fParam160;
		Steinberg::Vst::ParamValue fParam40;
		Steinberg::Vst::ParamValue fParamSub;

		bool            bParamBypass;

		overSample      fParamOS;


		// Internal values --------------------------------------------------------------

		SVF_ svf_Sub[2], svf_40[2], svf_160[2], svf_650[2], svf_2k5[2], svf_Sky[2];
		SVF_ svf_LC[2], svf_HC[2];

		Steinberg::Vst::Sample64 BP_gain_Sub;
		Steinberg::Vst::Sample64 BP_gain_40;
		Steinberg::Vst::Sample64 BP_gain_160;
		Steinberg::Vst::Sample64 BP_gain_650;
		Steinberg::Vst::Sample64 HP_gain_2k5;
		Steinberg::Vst::Sample64 HP_gain_Sky;

		const Steinberg::Vst::Sample64 BP_arr[21] = {
			exp(log(10.0) * (6.0 - 7.5) / 20.0), //  0/20  -5.0
			exp(log(10.0) * (6.0 - 7.4) / 20.0), //  1/20  -4.5
			exp(log(10.0) * (6.0 - 7.3) / 20.0), //  2/20  -4.0
			exp(log(10.0) * (6.0 - 7.2) / 20.0), //  3/20  -3.5
			exp(log(10.0) * (6.0 - 6.8) / 20.0), //  4/20  -3.0
			exp(log(10.0) * (6.0 - 6.1) / 20.0), //  5/20  -2.5
			exp(log(10.0) * (6.0 - 5.7) / 20.0), //  6/20  -2.0
			exp(log(10.0) * (6.0 - 4.5) / 20.0), //  7/20  -1.5
			exp(log(10.0) * (6.0 - 3.0) / 20.0), //  8/20  -1.0
			exp(log(10.0) * (6.0 - 1.0) / 20.0), //  9/20  -0.5
			exp(log(10.0) * (6.0) / 20.0), // 10/20   0.0
			exp(log(10.0) * (6.0 + 1.2) / 20.0), // +0.5
			exp(log(10.0) * (6.0 + 2.3) / 20.0), // +1.0
			exp(log(10.0) * (6.0 + 3.5) / 20.0), // +1.5
			exp(log(10.0) * (6.0 + 5.7) / 20.0), // +2.0
			exp(log(10.0) * (6.0 + 8.75) / 20.0), // +2.5
			exp(log(10.0) * (6.0 + 9.65) / 20.0), // +3.0
			exp(log(10.0) * (6.0 + 10.8) / 20.0), // +3.5
			exp(log(10.0) * (6.0 + 12.05) / 20.0), // +4.0
			exp(log(10.0) * (6.0 + 13.4) / 20.0), // +4.5
			exp(log(10.0) * (6.0 + 13.5) / 20.0)  // +5.0
		};
		const Steinberg::Vst::Sample64 HP_arr[21] = {
			exp(log(10.0) * (3.15 - 7.5) / 20.0), //  0/20  -5.0
			exp(log(10.0) * (3.15 - 7.4) / 20.0), //  1/20  -4.5
			exp(log(10.0) * (3.15 - 7.3) / 20.0), //  2/20  -4.0
			exp(log(10.0) * (3.15 - 7.2) / 20.0), //  3/20  -3.5
			exp(log(10.0) * (3.15 - 6.8) / 20.0), //  4/20  -3.0
			exp(log(10.0) * (3.15 - 6.1) / 20.0), //  5/20  -2.5
			exp(log(10.0) * (3.15 - 5.7) / 20.0), //  6/20  -2.0
			exp(log(10.0) * (3.15 - 4.5) / 20.0), //  7/20  -1.5
			exp(log(10.0) * (3.15 - 3.0) / 20.0), //  8/20  -1.0
			exp(log(10.0) * (3.15 - 1.0) / 20.0), //  9/20  -0.5
			exp(log(10.0) * (3.15) / 20.0), // 10/20   0.0
			exp(log(10.0) * (3.15 + 1.2) / 20.0), // +0.5
			exp(log(10.0) * (3.15 + 2.3) / 20.0), // +1.0
			exp(log(10.0) * (3.15 + 3.5) / 20.0), // +1.5
			exp(log(10.0) * (3.15 + 5.7) / 20.0), // +2.0
			exp(log(10.0) * (3.15 + 8.75) / 20.0), // +2.5
			exp(log(10.0) * (3.15 + 9.65) / 20.0), // +3.0
			exp(log(10.0) * (3.15 + 10.8) / 20.0), // +3.5
			exp(log(10.0) * (3.15 + 12.05) / 20.0), // +4.0
			exp(log(10.0) * (3.15 + 13.4) / 20.0), // +4.5
			exp(log(10.0) * (3.15 + 13.5) / 20.0)  // +5.0
		};
		const Steinberg::Vst::Sample64 HP_Sky_arr[21] = {
			exp(log(10.0) * (0.0 - 16.0) / 20.0), //  0.0
			exp(log(10.0) * (0.0 - 15.8) / 20.0), //  0.5
			exp(log(10.0) * (0.0 - 15.6) / 20.0), //  1.0
			exp(log(10.0) * (0.0 - 15.0) / 20.0), //  1.5
			exp(log(10.0) * (0.0 - 14.5) / 20.0), //  2.0
			exp(log(10.0) * (0.0 - 13.0) / 20.0), //  2.5
			exp(log(10.0) * (0.0 - 11.8) / 20.0), //  3.0
			exp(log(10.0) * (0.0 - 9.4) / 20.0), //  3.5
			exp(log(10.0) * (0.0 - 5.9) / 20.0), //  4.0
			exp(log(10.0) * (0.0 - 3.1) / 20.0), //  4.5
			exp(log(10.0) * (0.0 - 1.5) / 20.0), // 5.0
			exp(log(10.0) * (0.0 - 0.3) / 20.0), // 5.5
			exp(log(10.0) * (0.0 + 1.2) / 20.0), // 6.0
			exp(log(10.0) * (0.0 + 2.8) / 20.0), // 6.5
			exp(log(10.0) * (0.0 + 5.2) / 20.0), // 7.0
			exp(log(10.0) * (0.0 + 10.7) / 20.0), // 7.5
			exp(log(10.0) * (0.0 + 11.7) / 20.0), // 8.0
			exp(log(10.0) * (0.0 + 13.1) / 20.0), // 8.5  <- EQ4M max
			exp(log(10.0) * (0.0 + 14.6) / 20.0), // 9.0
			exp(log(10.0) * (0.0 + 16.3) / 20.0), // 9.5
			exp(log(10.0) * (0.0 + 16.4) / 20.0)  // 10.
		};

		// Buffers ------------------------------------------------------------------
		typedef struct _Flt {
			double coef alignas(16)[maxTap] = { 0, };
			double buff alignas(16)[maxTap] = { 0, };
		} Flt;

		Flt dnSample_21[2];
		Flt dnSample_42[2];

		const Steinberg::int32 dnTap_21 = 49;
		const Steinberg::int32 dnTap_42 = 193;

		void Fir_x2_dn(Steinberg::Vst::Sample64* in, Steinberg::Vst::Sample64* out, Steinberg::int32 channel);
		void Fir_x4_dn(Steinberg::Vst::Sample64* in, Steinberg::Vst::Sample64* out, Steinberg::int32 channel);

		const Steinberg::int32 latency_Fir_x2 = 12;
		const Steinberg::int32 latency_Fir_x4 = 24;

		std::queue<double> latency_q[2];
	};
	//------------------------------------------------------------------------
} // namespace yg331
