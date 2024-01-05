//------------------------------------------------------------------------
// Copyright(c) 2023 yg331.
//------------------------------------------------------------------------

#pragma once

#include "SBEQ4_cids.h"
#include "public.sdk/source/vst/vstaudioeffect.h"

#include <cstdlib>
#include <stdlib.h>
#include <math.h>

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

		template <typename SampleType>
		void overSampling(SampleType** inputs, SampleType** outputs, Steinberg::int32 numChannels, Steinberg::Vst::SampleRate getSampleRate, Steinberg::int32 sampleFrames);

		template <typename SampleType>
		void proc_out(Steinberg::Vst::Sample64** inputs, SampleType** outputs, Steinberg::int32 numChannels, Steinberg::int32 sampleFrames);

		template <typename SampleType>
		void bypass_latency(SampleType** inputs, SampleType** outputs, Steinberg::int32 numChannels, Steinberg::int32 sampleFrames);

		void FIR_dn_2_1(Steinberg::Vst::Sample64** inputs, Steinberg::Vst::Sample64** outputs, Steinberg::int32 numChannels, long sampleFrames);
		void FIR_dn_4_1(Steinberg::Vst::Sample64** inputs, Steinberg::Vst::Sample64** outputs, Steinberg::int32 numChannels, long sampleFrames);

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
			Steinberg::Vst::Sample64** outputs,
			Steinberg::int32 numChannels,
			Steinberg::Vst::SampleRate getSampleRate,
			Steinberg::int32 sampleFrames
		);

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* pi             */
#endif
#define maxTap 128
#define halfTap 64

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
		Flt dnSample_41[2];
		Flt dnSample_42[2];

		const int32 dnTap_21 = 113;
		const int32 dnTap_41 = 113;
		const int32 dnTap_42 = 31;

		void Fir_x2_dn(Vst::Sample64* in, Vst::Sample64* out, int32 channel);
		void Fir_x4_dn(Vst::Sample64* in, Vst::Sample64* out, int32 channel);

		const int32 latency_Fir_x2 = 49;
		const int32 latency_Fir_x4 = 56;


#define maxSample 65536


#define Num_coef_dn_2 49 // 49 - 1 = 4 * 12
		Steinberg::Vst::Sample64 IR_coef_Dn_2[Num_coef_dn_2] =
		{
			 8.83644003748518e-008, -1.94193278915946e-005, -4.4985037174345e-007,   0.000125782066263887,
			 2.00732060683587e-006, -0.000498654865053702,  -6.13912560981722e-006,  0.00150429864464152,
			 1.49387715402172e-005, -0.00379464215391734,   -3.07373675173277e-005,  0.00839603954155323,
			 5.52843541114737e-005, -0.0168254657370688,    -8.87461409955332e-005,  0.0313574932644712,
			 0.000128900424304018,  -0.055906813113935,     -0.000170970557855707,   0.0994526426712894,
			 0.000208372850479974,  -0.1939599822113,       -0.000234286008046444,   0.628983177782361,
			 0.997874719778748,      0.628983177782361,     -0.000234286008046444,  -0.1939599822113,
			 0.000208372850479974,   0.0994526426712895,    -0.000170970557855707,  -0.055906813113935,
			 0.000128900424304018,   0.0313574932644712,    -8.87461409955332e-005, -0.0168254657370688,
			 5.52843541114738e-005,  0.00839603954155324,   -3.07373675173278e-005, -0.00379464215391735,
			 1.49387715402172e-005,  0.00150429864464152,   -6.13912560981721e-006, -0.000498654865053703,
			 2.00732060683587e-006,  0.000125782066263886,  -4.49850371743453e-007, -1.94193278915941e-005,
			 8.83644003748518e-008
		};


#define Num_buff_dn_2 49
		Steinberg::Vst::Sample64 IR_buff_Dn_2[2][Num_buff_dn_2 + 1] = { 0, };



#define Num_coef_dn_22 193 // 193 - 1 = 8 * 24
		Steinberg::Vst::Sample64 IR_coef_Dn_22[Num_coef_dn_22] =
		{
			 1.76310158052763e-007, -3.48368495625567e-006, -6.0099865839085e-006,  -5.50400561555831e-006,
			-3.41769436735522e-007,  8.58485797017538e-006,  1.62918687387007e-005,  1.52981747652374e-005,
			 8.97569494114445e-007, -2.28489834359722e-005, -4.16081552178807e-005, -3.74608091514729e-005,
			-2.02119091623382e-006,  5.20879426113635e-005,  9.16491242914849e-005,  7.99710103582398e-005,
			 4.0051311608801e-006,  -0.000105893439007756,  -0.000181840004621746,  -0.000155142609576021,
			-7.24972071716962e-006,  0.000198461085292744,   0.00033453139737088,    0.000280479212741291,
			 1.22491659761287e-005, -0.000349502190511278,  -0.000580293699380777,  -0.000479587230013075,
			-1.95681184070956e-005,  0.000585208004310622,   0.000959262920135358,   0.000783109163315969,
			 2.98067678861249e-005, -0.000939237697939664,  -0.00152250324427374,   -0.00122965912999226,
			-4.35538755036053e-005,  0.00145373317335439,    0.0023334096714335,     0.0018667937370317,
			 6.13291110686729e-005, -0.00218043248134293,   -0.00346928976792548,   -0.00275214765743677,
			-8.35182440276637e-005,  0.00318207702738813,    0.00502346550581221,    0.00395502767072196,
			 0.000110306788366018,  -0.00453451812070763,   -0.00710857794319149,   -0.00555903639382876,
			-0.000141619248313517,   0.0063302910754675,     0.0098623702270716,     0.00766678247956241,
			 0.000177071830727303,  -0.00868508179689356,   -0.0134583134092389,    -0.0104086408043634,
			-0.000215946167431355,   0.0117497875343864,     0.0181255975994597,     0.0139593594413516,
			 0.000257190158997307,  -0.0157335677135839,    -0.0241876646315684,    -0.0185703526566454,
			-0.000299449608033949,   0.0209495006629276,     0.0321394650660784,     0.0246353228973194,
			 0.000341131111058696,  -0.0279103612316827,    -0.042812571432583,     -0.0328334262542524,
			-0.000380493106911436,   0.037548028497909,      0.0577641053572444,     0.0444770383190673,
			 0.000415758496025339,  -0.0517864322124472,    -0.0803355888525183,    -0.0625034998056695,
			-0.000445239320218456,   0.0753663029401762,     0.119254735043198,      0.0951120232183091,
			 0.000467462042779561,  -0.123925877560685,     -0.206642025100748,     -0.176854177072887,
			-0.000481281292291734,   0.296719038135957,      0.632073495820056,      0.895856663391142,
			 0.995753897050922,      0.895856663391142,      0.632073495820056,      0.296719038135957,
			-0.000481281292291734,  -0.176854177072887,     -0.206642025100748,     -0.123925877560685,
			 0.000467462042779562,   0.0951120232183091,     0.119254735043198,      0.0753663029401762,
			-0.000445239320218456,  -0.0625034998056695,    -0.0803355888525183,    -0.0517864322124473,
			 0.000415758496025339,   0.0444770383190673,     0.0577641053572444,     0.037548028497909,
			-0.000380493106911436,  -0.0328334262542524,    -0.042812571432583,     -0.0279103612316827,
			 0.000341131111058696,   0.0246353228973194,     0.0321394650660784,     0.0209495006629276,
			-0.000299449608033949,  -0.0185703526566454,    -0.0241876646315684,    -0.0157335677135839,
			 0.000257190158997307,   0.0139593594413516,     0.0181255975994597,     0.0117497875343864,
			-0.000215946167431355,  -0.0104086408043634,    -0.0134583134092389,    -0.00868508179689357,
			 0.000177071830727303,   0.0076667824795624,     0.00986237022707161,    0.00633029107546751,
			-0.000141619248313517,  -0.00555903639382877,   -0.00710857794319148,   -0.00453451812070763,
			 0.000110306788366018,   0.00395502767072196,    0.00502346550581222,    0.00318207702738813,
			-8.35182440276638e-005, -0.00275214765743678,   -0.00346928976792548,   -0.00218043248134294,
			 6.13291110686729e-005,  0.0018667937370317,     0.0023334096714335,     0.00145373317335439,
			-4.35538755036053e-005, -0.00122965912999226,   -0.00152250324427374,   -0.000939237697939664,
			 2.98067678861249e-005,  0.000783109163315971,   0.00095926292013536,    0.000585208004310623,
			-1.95681184070957e-005, -0.000479587230013076,  -0.000580293699380777,  -0.000349502190511277,
			 1.22491659761287e-005,  0.000280479212741292,   0.00033453139737088,    0.000198461085292744,
			-7.24972071716964e-006, -0.000155142609576022,  -0.000181840004621746,  -0.000105893439007756,
			 4.0051311608801e-006,   7.99710103582402e-005,  9.16491242914852e-005,  5.20879426113636e-005,
			-2.02119091623381e-006, -3.74608091514732e-005, -4.1608155217881e-005,  -2.28489834359722e-005,
			 8.97569494114452e-007,  1.52981747652369e-005,  1.62918687387008e-005,  8.58485797017508e-006,
			-3.41769436735514e-007, -5.5040056155583e-006,  -6.0099865839085e-006,  -3.48368495625569e-006,
			 1.76310158052763e-007
		};



#define Num_buff_dn_22 193
		Steinberg::Vst::Sample64 IR_buff_Dn_22[2][Num_buff_dn_22 + 1] = { 0, };





	};
	//------------------------------------------------------------------------
} // namespace yg331
