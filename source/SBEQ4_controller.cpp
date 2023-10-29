//------------------------------------------------------------------------
// Copyright(c) 2023 yg331.
//------------------------------------------------------------------------

#include "SBEQ4_controller.h"
// #include "SBEQ4_cids.h"
#include "vstgui/plugin-bindings/vst3editor.h"

#include "pluginterfaces/base/ustring.h"
#include "base/source/fstreamer.h"
#include "public.sdk/source/vst/vsteditcontroller.h"

#include "vstgui/vstgui_uidescription.h"
#include "vstgui/uidescription/detail/uiviewcreatorattributes.h"

using namespace Steinberg;

namespace yg331 {


	//------------------------------------------------------------------------
	// Sky_Blue_EQ4Controller Implementation
	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Controller::initialize(FUnknown* context)
	{
		// Here the Plug-in will be instantiated

		//---do not forget to call parent ------
		tresult result = EditControllerEx1::initialize(context);
		if (result != kResultOk)
		{
			return result;
		}

		// Here you could register some parameters

		int32 stepCount;
		int32 flags;
		int32 tag;
		Vst::ParamValue defaultVal;
		Vst::ParamValue minPlain;
		Vst::ParamValue maxPlain;
		Vst::ParamValue defaultPlain;

		flags = Vst::ParameterInfo::kCanAutomate;

		tag = kParamInputAtten;
		auto* ParamInputAtten = new Vst::StringListParameter(STR16("InputAtten"), tag, STR16(""), flags);
		for (int i = 0; i < knob_stepCount + 1; i++) ParamInputAtten->appendString(Atten_knob[i]);
		//ParamInputAtten->setNormalized(ParamInputAtten->toNormalized(20));
		ParamInputAtten->getInfo().defaultNormalizedValue = 1.0;
		parameters.addParameter(ParamInputAtten);

		tag = kParamSkyFreq;
		auto* Air_Freq = new Vst::StringListParameter(STR("Sky Freq"), tag, STR16(""), flags);
		for (int i = 0; i < 7; i++) Air_Freq->appendString(Sky_freq[i]);
		//Air_Freq->setNormalized(Air_Freq->toNormalized(0));
		Air_Freq->getInfo().defaultNormalizedValue = 0.0;
		//Air_Freq->addDependent(this);
		parameters.addParameter(Air_Freq);

		tag = kParamSky;
		auto* ParamAir = new Vst::StringListParameter(STR("Sky"), tag, STR16(""), flags);
		for (int i = 0; i < knob_stepCount + 1; i++) ParamAir->appendString(Sky_band[i]);
		//ParamAir->setNormalized(ParamAir->toNormalized(0));
		ParamAir->getInfo().defaultNormalizedValue = 0.0;
		parameters.addParameter(ParamAir);


		stepCount = 20;
		minPlain = -5.0;
		maxPlain = 5.0;
		defaultPlain = 0.0;



		tag = kParam2k5;
		auto* Param2k5 = new Vst::StringListParameter(STR16("2k5"), tag, STR16(""), flags);
		for (int i = 0; i < knob_stepCount + 1; i++) Param2k5->appendString(Other_band[i]);
		//Param2k5->setNormalized(0.5);
		Param2k5->getInfo().defaultNormalizedValue = 0.5;
		//auto* Param2k5 = new Vst::RangeParameter(STR16("2k5"), tag, STR16(""), -5.0, +5.0, 0.0, 20, flags);
		parameters.addParameter(Param2k5);
		tag = kParam650;
		auto* Param650 = new Vst::StringListParameter(STR16("650"), tag, STR16(""), flags);
		for (int i = 0; i < knob_stepCount + 1; i++) Param650->appendString(Other_band[i]);
		//Param650->setNormalized(Param650->toNormalized(10));
		Param650->getInfo().defaultNormalizedValue = 0.5;
		parameters.addParameter(Param650);
		tag = kParam160;
		auto* Param160 = new Vst::StringListParameter(STR16("160"), tag, STR16(""), flags);
		for (int i = 0; i < knob_stepCount + 1; i++) Param160->appendString(Other_band[i]);
		//Param160->setNormalized(Param160->toNormalized(10));
		Param160->getInfo().defaultNormalizedValue = 0.5;
		parameters.addParameter(Param160);
		tag = kParam40;
		auto* Param40 = new Vst::StringListParameter(STR16("40"), tag, STR16(""), flags);
		for (int i = 0; i < knob_stepCount + 1; i++) Param40->appendString(Other_band[i]);
		//Param40->setNormalized(Param40->toNormalized(10));
		Param40->getInfo().defaultNormalizedValue = 0.5;
		parameters.addParameter(Param40);
		tag = kParamSub;
		auto* ParamSub = new Vst::StringListParameter(STR16("Sub"), tag, STR16(""), flags);
		for (int i = 0; i < knob_stepCount + 1; i++) ParamSub->appendString(Other_band[i]);
		//Param10->setNormalized(ParamSub->toNormalized(10));
		ParamSub->getInfo().defaultNormalizedValue = 0.5;
		parameters.addParameter(ParamSub);

		Vst::ParamValue zoom_coef = 1.0;

		if (zoomFactors.empty())
		{
			zoomFactors.push_back(ZoomFactor(STR("50%"), zoom_coef * 0.5));  // 0/6
			zoomFactors.push_back(ZoomFactor(STR("75%"), zoom_coef * 0.75)); // 1/6
			zoomFactors.push_back(ZoomFactor(STR("100%"), zoom_coef * 1.0));  // 2/6
			zoomFactors.push_back(ZoomFactor(STR("125%"), zoom_coef * 1.25)); // 3/6
			zoomFactors.push_back(ZoomFactor(STR("150%"), zoom_coef * 1.5));  // 4/6
			zoomFactors.push_back(ZoomFactor(STR("175%"), zoom_coef * 1.75)); // 5/6
			zoomFactors.push_back(ZoomFactor(STR("200%"), zoom_coef * 2.0));  // 6/6
		}

		Vst::StringListParameter* zoomParameter = new Vst::StringListParameter(STR("Zoom"), kParamZoom);
		for (ZoomFactorVector::const_iterator it = zoomFactors.begin(), end = zoomFactors.end(); it != end; ++it)
		{
			zoomParameter->appendString(it->title);
		}
		zoomParameter->setNormalized(zoomParameter->toNormalized(Init_Zoom)); // toNorm(2) == 100%
		zoomParameter->addDependent(this);
		parameters.addParameter(zoomParameter);


		tag = kParamBypass;
		stepCount = 1;
		defaultVal = Init_Bypass ? 1 : 0;
		flags = Vst::ParameterInfo::kCanAutomate | Vst::ParameterInfo::kIsBypass;
		parameters.addParameter(STR16("Bypass"), nullptr, stepCount, defaultVal, flags, tag);



		return result;
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Controller::terminate()
	{
		// Here the Plug-in will be de-instantiated, last possibility to remove some memory!

		//---do not forget to call parent ------
		return EditControllerEx1::terminate();
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Controller::setComponentState(IBStream* state)
	{
		// Here you get the state of the component (Processor part)
		if (!state)
			return kResultFalse;

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

		setParamNormalized(kParamBypass, savedBypass ? 1 : 0);
		setParamNormalized(kParamZoom, savedZoom);

		setParamNormalized(kParamInputAtten, savedInputAtten);
		setParamNormalized(kParamSkyFreq, savedSkyFreq);
		setParamNormalized(kParamSky, savedSky);
		setParamNormalized(kParam2k5, saved2k5);
		setParamNormalized(kParam650, saved650);
		setParamNormalized(kParam160, saved160);
		setParamNormalized(kParam40, saved40);
		setParamNormalized(kParamSub, savedSub);

		return kResultOk;
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Controller::setState(IBStream* state)
	{
		// Here you get the state of the controller

		return kResultTrue;
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Controller::getState(IBStream* state)
	{
		// Here you are asked to deliver the state of the controller (if needed)
		// Note: the real state of your plug-in is saved in the processor

		return kResultTrue;
	}

	//------------------------------------------------------------------------
	IPlugView* PLUGIN_API Sky_Blue_EQ4Controller::createView(FIDString name)
	{
		// Here the Host wants to open your editor (if you have one)
		if (FIDStringsEqual(name, Vst::ViewType::kEditor))
		{
			// create your editor here and return a IPlugView ptr of it
			auto* view = new VSTGUI::VST3Editor(this, "view", "SBEQ4_editor.uidesc");
			view->setZoomFactor(0.5);
			setKnobMode(Steinberg::Vst::KnobModes::kLinearMode);
			return view;
		}
		return nullptr;
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Controller::setParamNormalized(Vst::ParamID tag, Vst::ParamValue value)
	{
		// called by host to update your parameters
		tresult result = EditControllerEx1::setParamNormalized(tag, value);
		return result;
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Controller::getParamStringByValue(Vst::ParamID tag, Vst::ParamValue valueNormalized, Vst::String128 string)
	{
		// called by host to get a string for given normalized value of a specific parameter
		// (without having to set the value!)
		return EditControllerEx1::getParamStringByValue(tag, valueNormalized, string);
	}

	//------------------------------------------------------------------------
	tresult PLUGIN_API Sky_Blue_EQ4Controller::getParamValueByString(Vst::ParamID tag, Vst::TChar* string, Vst::ParamValue& valueNormalized)
	{
		// called by host to get a normalized value from a string representation of a specific parameter
		// (without having to set the value!)
		return EditControllerEx1::getParamValueByString(tag, string, valueNormalized);
	}

	//------------------------------------------------------------------------



	void PLUGIN_API Sky_Blue_EQ4Controller::update(FUnknown* changedUnknown, int32 message)
	{
		EditControllerEx1::update(changedUnknown, message);

		// GUI Resizing
		// check 'zoomtest' code at
		// https://github.com/steinbergmedia/vstgui/tree/vstgui4_10/vstgui/tests/uidescription%20vst3/source

		Vst::Parameter* param = FCast<Vst::Parameter>(changedUnknown);
		if (!param)
			return;

		if (param->getInfo().id == kParamZoom)
		{
			size_t index = static_cast<size_t> (param->toPlain(param->getNormalized()));

			if (index >= zoomFactors.size())
				return;

			for (EditorVector::const_iterator it = editors.begin(), end = editors.end(); it != end; ++it)
			{
				VSTGUI::VST3Editor* editor = dynamic_cast<VSTGUI::VST3Editor*>(*it);
				if (editor)
					editor->setZoomFactor(zoomFactors[index].factor);
			}
		}
	}
	//------------------------------------------------------------------------
	void Sky_Blue_EQ4Controller::editorAttached(Steinberg::Vst::EditorView* editor)
	{
		editors.push_back(editor);
	}

	//------------------------------------------------------------------------
	void Sky_Blue_EQ4Controller::editorRemoved(Steinberg::Vst::EditorView* editor)
	{
		editors.erase(std::find(editors.begin(), editors.end(), editor));
	}
} // namespace yg331
