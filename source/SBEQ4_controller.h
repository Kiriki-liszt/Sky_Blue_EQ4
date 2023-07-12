//------------------------------------------------------------------------
// Copyright(c) 2023 yg331.
//------------------------------------------------------------------------

#pragma once

#include "SBEQ4_cids.h"
#include "public.sdk/source/vst/vsteditcontroller.h"
#include "vstgui/vstgui.h"

namespace yg331 {
	//------------------------------------------------------------------------
	//  Sky_Blue_EQ4Controller
	//------------------------------------------------------------------------
	class Sky_Blue_EQ4Controller : public Steinberg::Vst::EditControllerEx1
	{
	public:
		//------------------------------------------------------------------------
		Sky_Blue_EQ4Controller() = default;
		~Sky_Blue_EQ4Controller() SMTG_OVERRIDE = default;

		// Create function
		static Steinberg::FUnknown* createInstance(void* /*context*/)
		{
			return (Steinberg::Vst::IEditController*)new Sky_Blue_EQ4Controller;
		}

		// IPluginBase
		Steinberg::tresult PLUGIN_API initialize(Steinberg::FUnknown* context) SMTG_OVERRIDE;
		Steinberg::tresult PLUGIN_API terminate() SMTG_OVERRIDE;

		// EditController
		Steinberg::tresult PLUGIN_API setComponentState(Steinberg::IBStream* state) SMTG_OVERRIDE;
		Steinberg::IPlugView* PLUGIN_API createView(Steinberg::FIDString name) SMTG_OVERRIDE;
		Steinberg::tresult PLUGIN_API setState(Steinberg::IBStream* state) SMTG_OVERRIDE;
		Steinberg::tresult PLUGIN_API getState(Steinberg::IBStream* state) SMTG_OVERRIDE;
		Steinberg::tresult PLUGIN_API setParamNormalized(Steinberg::Vst::ParamID tag,
			Steinberg::Vst::ParamValue value) SMTG_OVERRIDE;
		Steinberg::tresult PLUGIN_API getParamStringByValue(Steinberg::Vst::ParamID tag,
			Steinberg::Vst::ParamValue valueNormalized,
			Steinberg::Vst::String128 string) SMTG_OVERRIDE;
		Steinberg::tresult PLUGIN_API getParamValueByString(Steinberg::Vst::ParamID tag,
			Steinberg::Vst::TChar* string,
			Steinberg::Vst::ParamValue& valueNormalized) SMTG_OVERRIDE;

		//---Interface---------
		DEFINE_INTERFACES
			// Here you can add more supported VST3 interfaces
			// DEF_INTERFACE (Vst::IXXX)
			END_DEFINE_INTERFACES(EditController)
			DELEGATE_REFCOUNT(EditController)

		void PLUGIN_API update(Steinberg::FUnknown* changedUnknown, Steinberg::int32 message) SMTG_OVERRIDE;
		void editorAttached(Steinberg::Vst::EditorView* editor) SMTG_OVERRIDE;
		void editorRemoved(Steinberg::Vst::EditorView* editor) SMTG_OVERRIDE;


		//------------------------------------------------------------------------
	protected:
		typedef std::vector<Steinberg::Vst::EditorView*> EditorVector;
		EditorVector editors;

		ZoomFactorVector zoomFactors;

		Steinberg::tchar* Sky_freq[7] = {
			STR16("off"), STR16("2.5k"),
			STR16("5.0k"), STR16("10k"),
			STR16("15k"), STR16("20k"),
			STR16("40k")
		};

		Steinberg::tchar* Sky_band[knob_stepCount + 1] = {
			STR16("0.0"), STR16("0.5"),
			STR16("1.0"), STR16("1.5"),
			STR16("2.0"), STR16("2.5"),
			STR16("3.0"), STR16("3.5"),
			STR16("4.0"), STR16("4.5"),
			STR16("5.0"), STR16("5.5"),
			STR16("6.0"), STR16("6.5"),
			STR16("7.0"), STR16("7.5"),
			STR16("8.0"), STR16("8.5"),
			STR16("9.0"), STR16("9.5"),
			STR16("10.0")
		};

		Steinberg::tchar* Other_band[knob_stepCount + 1] = {
			STR16("-5.0"), STR16("-4.5"),
			STR16("-4.0"), STR16("-3.5"),
			STR16("-3.0"), STR16("-2.5"),
			STR16("-2.0"), STR16("-1.5"),
			STR16("-1.0"), STR16("-0.5"),
			STR16("0.0"), STR16("+0.5"),
			STR16("+1.0"), STR16("+1.5"),
			STR16("+2.0"), STR16("+2.5"),
			STR16("+3.0"), STR16("+3.5"),
			STR16("+4.0"), STR16("+4.5"),
			STR16("+5.0")
		};

		Steinberg::tchar* Atten_knob[knob_stepCount + 1] = {
			STR16("-10.0"), STR16("-9.5"),
			STR16("-9.0"), STR16("-8.5"),
			STR16("-8.0"), STR16("-7.5"),
			STR16("-7.0"), STR16("-6.5"),
			STR16("-6.0"), STR16("-5.5"),
			STR16("-5.0"), STR16("-4.5"),
			STR16("-4.0"), STR16("-3.5"),
			STR16("-3.0"), STR16("-2.5"),
			STR16("-2.0"), STR16("-1.5"),
			STR16("-1.0"), STR16("-0.5"),
			STR16("0.0")
		};
	};


	using namespace VSTGUI;
	class MyControl : public CParamDisplay
	{
	public:
		MyControl(const CRect& size);
		CLASS_METHODS(MyControl, CParamDisplay);

		CMouseEventResult onMouseDown(CPoint& where, const CButtonState& buttons) override
		{
			return kMouseEventHandled; // needed to get the onMouseUp call
		}
		CMouseEventResult onMouseUp(CPoint& where, const CButtonState& buttons) override
		{
			if (buttons.isLeftButton() && getViewSize().pointInside(where))
				doMouseClick();
			return kMouseEventHandled;
		}
		void doMouseClick() {}
	};


	//------------------------------------------------------------------------
} // namespace yg331

