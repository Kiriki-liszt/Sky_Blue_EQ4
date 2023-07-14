//------------------------------------------------------------------------
// Copyright(c) 2023 yg331.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/base/funknown.h"
#include "pluginterfaces/vst/vsttypes.h"

#include <vector>

namespace yg331 {
	//------------------------------------------------------------------------
	static const Steinberg::FUID kSky_Blue_EQ4ProcessorUID (0xDEFABAA2, 0x8A665F5E, 0x85CE3376, 0x9A96179B);
	static const Steinberg::FUID kSky_Blue_EQ4ControllerUID (0x3C9EE2F6, 0x5E085DBF, 0xA5664481, 0xFF226295);

#define Sky_Blue_EQ4VST3Category "Fx"

	struct ZoomFactor {
		const Steinberg::tchar* title;
		double factor;

		ZoomFactor(const Steinberg::tchar* title, double factor) : title(title), factor(factor) {}
	};
	typedef std::vector<ZoomFactor> ZoomFactorVector;
	typedef enum {
		overSample_1x,
		overSample_2x,
		overSample_4x,
		overSample_8x,
		overSample_num = 3
	} overSample;
	enum {
		kParamBypass = 0,
		kParamZoom,

		kParamInputAtten,
		kParamSkyFreq,
		kParamSky,
		kParam2k5,
		kParam650,
		kParam160,
		kParam40,
		kParamSub
	};
	enum filter_type
	{
		kLP,
		kHP,
		kBP,
		kNotch,
		kPeak,
		kBell,
		kHS,
		kLS
	};
	typedef struct SVF_st
	{
		double vin = 0.0;
		double v1 = 0.0;
		double v2 = 0.0;
		double v3 = 0.0;
		double t0 = 0.0;
		double v0 = 0.0;
		double t1 = 0.0;
		double t2 = 0.0;
		double ic1eq = 0.0;
		double ic2eq = 0.0;

		double bellgaindB = 0.0;
		double A = 0.0;

		double Q = 0.0;
		double w = 0.0;
		double g0 = 0.0;
		double k0 = 0.0;
		double g = 0.0;
		double k = 0.0;
		double g_k = 0.0;
		double gt0 = 0.0;
		double gk0 = 0.0;
		double gt1 = 0.0;
		double gk1 = 0.0;
		double gt2 = 0.0;

		double m0 = 0.0;
		double m1 = 0.0;
		double m2 = 0.0;

		filter_type type;
	} SVF_;
	constexpr auto knob_stepCount = 20;

	const bool
		Init_Bypass = false;

	const Steinberg::Vst::ParamValue
		Init_InputAtten = 1.0,
		Init_SkyFreq = 0.0,
		Init_Sky = 0.0,
		Init_2k5 = 0.5,
		Init_650 = 0.5,
		Init_160 = 0.5,
		Init_40 = 0.5,
		Init_Sub = 0.5,
		Init_Zoom = 0.0 / 6.0;

	


	//------------------------------------------------------------------------
} // namespace yg331
