cmake_minimum_required(VERSION 3.15.0)
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13 CACHE STRING "")

if(SMTG_MAC)
    set(vst3sdk_SOURCE_DIR "/Users/yoon/VST_SDK/vst3sdk")
elseif(SMTG_WIN)
    set(vst3sdk_SOURCE_DIR "C:/VST_SDK/vst3sdk")
endif(SMTG_MAC)

if(NOT vst3sdk_SOURCE_DIR)
    message(FATAL_ERROR "Path to VST3 SDK is empty!")
endif()

project(Sky_Blue_EQ4
    # This is your plug-in version number. Change it here only.
    # Version number symbols usable in C++ can be found in
    # source/version.h and ${PROJECT_BINARY_DIR}/projectversion.h.
    VERSION 1.0.0.0 
    DESCRIPTION "Sky_Blue_EQ4 VST 3 Plug-in"
)

set(SMTG_VSTGUI_ROOT "${vst3sdk_SOURCE_DIR}")

add_subdirectory(${vst3sdk_SOURCE_DIR} ${PROJECT_BINARY_DIR}/vst3sdk)
smtg_enable_vst3_sdk()

smtg_add_vst3plugin(Sky_Blue_EQ4
    source/version.h
    source/SBEQ4_cids.h
    source/SBEQ4_processor.h
    source/SBEQ4_processor.cpp
    source/SBEQ4_controller.h
    source/SBEQ4_controller.cpp
    source/SBEQ4_entry.cpp
)

#- VSTGUI Wanted ----
if(SMTG_ADD_VSTGUI)
    target_sources(Sky_Blue_EQ4
        PRIVATE
            resource/SBEQ4_editor.uidesc
    )
    target_link_libraries(Sky_Blue_EQ4
        PRIVATE
            vstgui_support
    )
    smtg_target_add_plugin_resources(Sky_Blue_EQ4
        RESOURCES
            "resource/SBEQ4_editor.uidesc"
            "resource/SB_background.png"
            "resource/SB_bypass.png"
            "resource/SB_knob_air.png"
            "resource/SB_knob_black.png"
            "resource/SB_knob_blue.png"
            "resource/SB_knob_green.png"
            "resource/SB_knob_magenta.png"
            "resource/SB_knob_orange.png"
            "resource/SB_knob_white.png"
            "resource/SB_knob_yellow.png"
    )
endif(SMTG_ADD_VSTGUI)
# -------------------

smtg_target_add_plugin_snapshots (Sky_Blue_EQ4
    RESOURCES
        resource/120C7E3245DD5916ACD2E8E6AD51E8FD_snapshot.png
        resource/120C7E3245DD5916ACD2E8E6AD51E8FD_snapshot_2.0x.png
)

target_link_libraries(Sky_Blue_EQ4
    PRIVATE
        sdk
)

smtg_target_configure_version_file(Sky_Blue_EQ4)

if(SMTG_MAC)
    smtg_target_set_bundle(Sky_Blue_EQ4
        BUNDLE_IDENTIFIER io.githuib.yg331.SkyBlueEQ4
        COMPANY_NAME "yg331"
    )
    smtg_target_set_debug_executable(Sky_Blue_EQ4
        "/Applications/VST3PluginTestHost.app"
        "--pluginfolder;$(BUILT_PRODUCTS_DIR)"
    )
elseif(SMTG_WIN)
    target_sources(Sky_Blue_EQ4 PRIVATE 
        resource/win32resource.rc
    )
    if(MSVC)
        set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Sky_Blue_EQ4)

        smtg_target_set_debug_executable(Sky_Blue_EQ4
            "$(ProgramW6432)/Steinberg/VST3PluginTestHost/VST3PluginTestHost.exe"
            "--pluginfolder \"$(OutDir)/\""
        )
    endif()
endif(SMTG_MAC)

# Add an AUv2 target
if (SMTG_MAC AND XCODE AND SMTG_COREAUDIO_SDK_PATH)
	smtg_target_add_auv2(Sky_Blue_EQ4-au
		BUNDLE_NAME Sky_Blue_EQ4
		BUNDLE_IDENTIFIER io.githuib.yg331.SkyBlueEQ4.audiounit
		INFO_PLIST_TEMPLATE resource/au-info.plist
		VST3_PLUGIN_TARGET Sky_Blue_EQ4)
endif(SMTG_MAC AND XCODE AND SMTG_COREAUDIO_SDK_PATH)
