/* Copyright (c) 2021 Evan Nemerson <evan@nemerson.com>
 *               2023 Michael R. Crusoe <crusoe@debian.org>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#define SIMDE_TEST_WASM_SIMD128_INSN make
#include "../../../simde/wasm/simd128.h"
#include "test-simd128.h"

static int
test_simde_wasm_i8x16_make(SIMDE_MUNIT_TEST_ARGS) {
  #if 1
    SIMDE_TEST_STRUCT_MODIFIERS struct {
      int8_t a[sizeof(simde_v128_t) / sizeof(int8_t)];
      int8_t r[sizeof(simde_v128_t) / sizeof(int8_t)];
    } test_vec[8] = {
      { { -INT8_C( 107),  INT8_C( 107), -INT8_C(  58), -INT8_C( 126), -INT8_C(  98), -INT8_C(  65),  INT8_C(  19), -INT8_C(  67),
          -INT8_C(  67), -INT8_C( 114), -INT8_C(  37), -INT8_C( 123), -INT8_C(  98),  INT8_C(  78), -INT8_C(  65), -INT8_C(  80) },
        { -INT8_C( 107),  INT8_C( 107), -INT8_C(  58), -INT8_C( 126), -INT8_C(  98), -INT8_C(  65),  INT8_C(  19), -INT8_C(  67),
          -INT8_C(  67), -INT8_C( 114), -INT8_C(  37), -INT8_C( 123), -INT8_C(  98),  INT8_C(  78), -INT8_C(  65), -INT8_C(  80) } },
      { {  INT8_C(  29),  INT8_C(  90),  INT8_C(  26),  INT8_C(  88), -INT8_C(  38), -INT8_C(  55),  INT8_C(   3), -INT8_C( 126),
           INT8_C( 113), -INT8_C(   3),  INT8_C(  55),  INT8_C(  65), -INT8_C(  12), -INT8_C( 124), -INT8_C(  62), -INT8_C(  82) },
        {  INT8_C(  29),  INT8_C(  90),  INT8_C(  26),  INT8_C(  88), -INT8_C(  38), -INT8_C(  55),  INT8_C(   3), -INT8_C( 126),
           INT8_C( 113), -INT8_C(   3),  INT8_C(  55),  INT8_C(  65), -INT8_C(  12), -INT8_C( 124), -INT8_C(  62), -INT8_C(  82) } },
      { {  INT8_C(  53), -INT8_C(  50), -INT8_C( 115), -INT8_C(  54), -INT8_C( 115), -INT8_C(  33), -INT8_C(   3),  INT8_C(  86),
          -INT8_C(  60),  INT8_C(  84), -INT8_C(  25), -INT8_C(  42), -INT8_C(  35),  INT8_C(  87), -INT8_C(   8),  INT8_C(  46) },
        {  INT8_C(  53), -INT8_C(  50), -INT8_C( 115), -INT8_C(  54), -INT8_C( 115), -INT8_C(  33), -INT8_C(   3),  INT8_C(  86),
          -INT8_C(  60),  INT8_C(  84), -INT8_C(  25), -INT8_C(  42), -INT8_C(  35),  INT8_C(  87), -INT8_C(   8),  INT8_C(  46) } },
      { { -INT8_C(  59),  INT8_C(  20),  INT8_C(  20), -INT8_C( 118),  INT8_C( 122), -INT8_C(  42),  INT8_C(  73),  INT8_C(  18),
          -INT8_C( 117), -INT8_C(   4),  INT8_C(  42), -INT8_C(  90),  INT8_C(  20),  INT8_C( 120),  INT8_C(  89),  INT8_C(  23) },
        { -INT8_C(  59),  INT8_C(  20),  INT8_C(  20), -INT8_C( 118),  INT8_C( 122), -INT8_C(  42),  INT8_C(  73),  INT8_C(  18),
          -INT8_C( 117), -INT8_C(   4),  INT8_C(  42), -INT8_C(  90),  INT8_C(  20),  INT8_C( 120),  INT8_C(  89),  INT8_C(  23) } },
      { {  INT8_C(  82),  INT8_C(  30), -INT8_C( 127), -INT8_C(  30),  INT8_C(  91), -INT8_C(  12),  INT8_C(  33),  INT8_C( 103),
          -INT8_C(  81),  INT8_C(  91),  INT8_C(  21), -INT8_C(  66),  INT8_C( 116), -INT8_C(  76),  INT8_C( 105),  INT8_C( 101) },
        {  INT8_C(  82),  INT8_C(  30), -INT8_C( 127), -INT8_C(  30),  INT8_C(  91), -INT8_C(  12),  INT8_C(  33),  INT8_C( 103),
          -INT8_C(  81),  INT8_C(  91),  INT8_C(  21), -INT8_C(  66),  INT8_C( 116), -INT8_C(  76),  INT8_C( 105),  INT8_C( 101) } },
      { { -INT8_C(  27),  INT8_C( 109),  INT8_C(  27),  INT8_C(  30), -INT8_C(   1),  INT8_C( 122),  INT8_C(  34),  INT8_C(  84),
           INT8_C(  42), -INT8_C(  83),  INT8_C(  44), -INT8_C(  56), -INT8_C(  41),  INT8_C(  94), -INT8_C(  44), -INT8_C( 111) },
        { -INT8_C(  27),  INT8_C( 109),  INT8_C(  27),  INT8_C(  30), -INT8_C(   1),  INT8_C( 122),  INT8_C(  34),  INT8_C(  84),
           INT8_C(  42), -INT8_C(  83),  INT8_C(  44), -INT8_C(  56), -INT8_C(  41),  INT8_C(  94), -INT8_C(  44), -INT8_C( 111) } },
      { {  INT8_C( 104), -INT8_C(  60), -INT8_C(  98), -INT8_C(  15),  INT8_C(  65), -INT8_C(  49),  INT8_C( 124),  INT8_C( 112),
           INT8_C(  39),  INT8_C( 103),  INT8_C(  36), -INT8_C(  79), -INT8_C( 124),  INT8_C(  23), -INT8_C( 106), -INT8_C( 110) },
        {  INT8_C( 104), -INT8_C(  60), -INT8_C(  98), -INT8_C(  15),  INT8_C(  65), -INT8_C(  49),  INT8_C( 124),  INT8_C( 112),
           INT8_C(  39),  INT8_C( 103),  INT8_C(  36), -INT8_C(  79), -INT8_C( 124),  INT8_C(  23), -INT8_C( 106), -INT8_C( 110) } },
      { { -INT8_C(   7),  INT8_C(  16),  INT8_C(  28),  INT8_C(  58),  INT8_C(  81), -INT8_C(  53),  INT8_C( 109),  INT8_C( 119),
          -INT8_C(  32), -INT8_C( 114),  INT8_C(  77),  INT8_C( 123), -INT8_C(  71),  INT8_C( 118),  INT8_C(  42), -INT8_C(  26) },
        { -INT8_C(   7),  INT8_C(  16),  INT8_C(  28),  INT8_C(  58),  INT8_C(  81), -INT8_C(  53),  INT8_C( 109),  INT8_C( 119),
          -INT8_C(  32), -INT8_C( 114),  INT8_C(  77),  INT8_C( 123), -INT8_C(  71),  INT8_C( 118),  INT8_C(  42), -INT8_C(  26) } }
    };

    for (size_t i = 0 ; i < (sizeof(test_vec) / sizeof(test_vec[0])) ; i++) {
      const int8_t* a = test_vec[i].a;
      simde_v128_t r =
        simde_wasm_i8x16_make(
            a[0], a[1], a[ 2], a[ 3], a[ 4], a[ 5], a[ 6], a[ 7],
            a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15]);
      simde_test_wasm_i8x16_assert_equal(r, simde_wasm_v128_load(test_vec[i].r));
    }
    return 0;
  #else
    fputc('\n', stdout);
    for (int i = 0 ; i < 8 ; i++) {
      int8_t a[sizeof(simde_v128_t) / sizeof(int8_t)];
      simde_v128_t r;

      simde_test_codegen_random_memory(sizeof(a), HEDLEY_REINTERPRET_CAST(uint8_t*, a));

      r = simde_wasm_i8x16_make(
          a[0], a[1], a[ 2], a[ 3], a[ 4], a[ 5], a[ 6], a[ 7],
          a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15]);

      simde_test_codegen_write_vi8(3, sizeof(a) / sizeof(a[0]), a, SIMDE_TEST_VEC_POS_FIRST);
      simde_test_wasm_i8x16_write(3, r, SIMDE_TEST_VEC_POS_LAST);
    }
    return 1;
  #endif
}

static int
test_simde_wasm_u8x16_make(SIMDE_MUNIT_TEST_ARGS) {
  #if 1
    SIMDE_TEST_STRUCT_MODIFIERS struct {
      uint8_t a[sizeof(simde_v128_t) / sizeof(uint8_t)];
      uint8_t r[sizeof(simde_v128_t) / sizeof(uint8_t)];
    } test_vec[] = {
      { { UINT8_C(  9), UINT8_C(209), UINT8_C(247), UINT8_C(100), UINT8_C(114), UINT8_C(139), UINT8_C(221), UINT8_C(224),
          UINT8_C(200), UINT8_C( 17), UINT8_C( 43), UINT8_C(165), UINT8_C( 92), UINT8_C(164), UINT8_C(111), UINT8_C( 35) },
        { UINT8_C(  9), UINT8_C(209), UINT8_C(247), UINT8_C(100), UINT8_C(114), UINT8_C(139), UINT8_C(221), UINT8_C(224),
          UINT8_C(200), UINT8_C( 17), UINT8_C( 43), UINT8_C(165), UINT8_C( 92), UINT8_C(164), UINT8_C(111), UINT8_C( 35) } },
      { { UINT8_C(147), UINT8_C(100), UINT8_C(235), UINT8_C( 87), UINT8_C(191), UINT8_C(237), UINT8_C(  6), UINT8_C( 93),
          UINT8_C( 85), UINT8_C(109), UINT8_C(140), UINT8_C(221), UINT8_C(126), UINT8_C( 96), UINT8_C( 36), UINT8_C(112) },
        { UINT8_C(147), UINT8_C(100), UINT8_C(235), UINT8_C( 87), UINT8_C(191), UINT8_C(237), UINT8_C(  6), UINT8_C( 93),
          UINT8_C( 85), UINT8_C(109), UINT8_C(140), UINT8_C(221), UINT8_C(126), UINT8_C( 96), UINT8_C( 36), UINT8_C(112) } },
      { { UINT8_C( 24), UINT8_C(115), UINT8_C(140), UINT8_C(171), UINT8_C(209),    UINT8_MAX, UINT8_C(145), UINT8_C(156),
          UINT8_C(219), UINT8_C(147), UINT8_C(191), UINT8_C( 89), UINT8_C(218), UINT8_C(190), UINT8_C(150), UINT8_C( 70) },
        { UINT8_C( 24), UINT8_C(115), UINT8_C(140), UINT8_C(171), UINT8_C(209),    UINT8_MAX, UINT8_C(145), UINT8_C(156),
          UINT8_C(219), UINT8_C(147), UINT8_C(191), UINT8_C( 89), UINT8_C(218), UINT8_C(190), UINT8_C(150), UINT8_C( 70) } },
      { { UINT8_C(223), UINT8_C(156), UINT8_C( 26), UINT8_C( 91), UINT8_C( 93), UINT8_C(111), UINT8_C(140), UINT8_C( 34),
          UINT8_C( 17), UINT8_C(212), UINT8_C(155), UINT8_C( 20), UINT8_C(125), UINT8_C(121), UINT8_C(225), UINT8_C( 92) },
        { UINT8_C(223), UINT8_C(156), UINT8_C( 26), UINT8_C( 91), UINT8_C( 93), UINT8_C(111), UINT8_C(140), UINT8_C( 34),
          UINT8_C( 17), UINT8_C(212), UINT8_C(155), UINT8_C( 20), UINT8_C(125), UINT8_C(121), UINT8_C(225), UINT8_C( 92) } },
      { { UINT8_C(203), UINT8_C(212), UINT8_C( 31), UINT8_C(184), UINT8_C(230), UINT8_C(233), UINT8_C( 33), UINT8_C( 22),
          UINT8_C( 86), UINT8_C(124), UINT8_C(111), UINT8_C( 86), UINT8_C(  3), UINT8_C( 87), UINT8_C(244), UINT8_C(165) },
        { UINT8_C(203), UINT8_C(212), UINT8_C( 31), UINT8_C(184), UINT8_C(230), UINT8_C(233), UINT8_C( 33), UINT8_C( 22),
          UINT8_C( 86), UINT8_C(124), UINT8_C(111), UINT8_C( 86), UINT8_C(  3), UINT8_C( 87), UINT8_C(244), UINT8_C(165) } },
      { { UINT8_C( 46), UINT8_C(233), UINT8_C( 13), UINT8_C( 87), UINT8_C( 97), UINT8_C(235), UINT8_C( 89), UINT8_C(169),
          UINT8_C( 50), UINT8_C( 49), UINT8_C(209), UINT8_C(140), UINT8_C(169), UINT8_C(207), UINT8_C(190), UINT8_C( 25) },
        { UINT8_C( 46), UINT8_C(233), UINT8_C( 13), UINT8_C( 87), UINT8_C( 97), UINT8_C(235), UINT8_C( 89), UINT8_C(169),
          UINT8_C( 50), UINT8_C( 49), UINT8_C(209), UINT8_C(140), UINT8_C(169), UINT8_C(207), UINT8_C(190), UINT8_C( 25) } },
      { { UINT8_C(110), UINT8_C(151), UINT8_C(159), UINT8_C(191), UINT8_C( 11), UINT8_C( 13), UINT8_C( 21), UINT8_C( 44),
          UINT8_C(129), UINT8_C(168), UINT8_C(249), UINT8_C(185), UINT8_C( 73), UINT8_C(225), UINT8_C( 63), UINT8_C(244) },
        { UINT8_C(110), UINT8_C(151), UINT8_C(159), UINT8_C(191), UINT8_C( 11), UINT8_C( 13), UINT8_C( 21), UINT8_C( 44),
          UINT8_C(129), UINT8_C(168), UINT8_C(249), UINT8_C(185), UINT8_C( 73), UINT8_C(225), UINT8_C( 63), UINT8_C(244) } },
      { { UINT8_C(  2), UINT8_C(194), UINT8_C(136), UINT8_C(110), UINT8_C(245), UINT8_C( 79), UINT8_C(159), UINT8_C( 82),
          UINT8_C(163), UINT8_C( 60), UINT8_C( 24), UINT8_C(240), UINT8_C(215), UINT8_C(179), UINT8_C(  4), UINT8_C(221) },
        { UINT8_C(  2), UINT8_C(194), UINT8_C(136), UINT8_C(110), UINT8_C(245), UINT8_C( 79), UINT8_C(159), UINT8_C( 82),
          UINT8_C(163), UINT8_C( 60), UINT8_C( 24), UINT8_C(240), UINT8_C(215), UINT8_C(179), UINT8_C(  4), UINT8_C(221) } },
    };

    for (size_t i = 0 ; i < (sizeof(test_vec) / sizeof(test_vec[0])) ; i++) {
      const uint8_t* a = test_vec[i].a;
      simde_v128_t r =
        simde_wasm_u8x16_make(
            a[0], a[1], a[ 2], a[ 3], a[ 4], a[ 5], a[ 6], a[ 7],
            a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15]);
      simde_test_wasm_i8x16_assert_equal(r, simde_wasm_v128_load(test_vec[i].r));
    }
    return 0;
  #else
    fputc('\n', stdout);
    for (int i = 0 ; i < 8 ; i++) {
      uint8_t a[sizeof(simde_v128_t) / sizeof(uint8_t)];
      simde_v128_t r;

      simde_test_codegen_random_memory(sizeof(a), HEDLEY_REINTERPRET_CAST(uint8_t*, a));

      r = simde_wasm_u8x16_make(
          a[0], a[1], a[ 2], a[ 3], a[ 4], a[ 5], a[ 6], a[ 7],
          a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15]);

      simde_test_codegen_write_vu8(3, sizeof(a) / sizeof(a[0]), a, SIMDE_TEST_VEC_POS_FIRST);
      simde_test_wasm_u8x16_write(3, r, SIMDE_TEST_VEC_POS_LAST);
    }
    return 1;
  #endif
}

static int
test_simde_wasm_i16x8_make(SIMDE_MUNIT_TEST_ARGS) {
  #if 1
    SIMDE_TEST_STRUCT_MODIFIERS struct {
      int16_t a[sizeof(simde_v128_t) / sizeof(int16_t)];
      int16_t r[sizeof(simde_v128_t) / sizeof(int16_t)];
    } test_vec[8] = {
      { {  INT16_C( 25276), -INT16_C( 21335),  INT16_C(  9881), -INT16_C( 22928),  INT16_C( 17610),  INT16_C( 25892), -INT16_C( 25524), -INT16_C( 25058) },
        {  INT16_C( 25276), -INT16_C( 21335),  INT16_C(  9881), -INT16_C( 22928),  INT16_C( 17610),  INT16_C( 25892), -INT16_C( 25524), -INT16_C( 25058) } },
      { {  INT16_C( 27009), -INT16_C( 14548),  INT16_C( 26207),  INT16_C( 31892), -INT16_C(  8887),  INT16_C( 25313),  INT16_C( 30268),  INT16_C( 32129) },
        {  INT16_C( 27009), -INT16_C( 14548),  INT16_C( 26207),  INT16_C( 31892), -INT16_C(  8887),  INT16_C( 25313),  INT16_C( 30268),  INT16_C( 32129) } },
      { {  INT16_C( 10968),  INT16_C( 23946),  INT16_C( 31163), -INT16_C( 11458), -INT16_C(  9405), -INT16_C( 11083), -INT16_C(   240), -INT16_C(  3905) },
        {  INT16_C( 10968),  INT16_C( 23946),  INT16_C( 31163), -INT16_C( 11458), -INT16_C(  9405), -INT16_C( 11083), -INT16_C(   240), -INT16_C(  3905) } },
      { {  INT16_C(  8561), -INT16_C( 13615), -INT16_C( 18924), -INT16_C( 30133), -INT16_C( 17435),  INT16_C(   448), -INT16_C( 11652),  INT16_C(  3750) },
        {  INT16_C(  8561), -INT16_C( 13615), -INT16_C( 18924), -INT16_C( 30133), -INT16_C( 17435),  INT16_C(   448), -INT16_C( 11652),  INT16_C(  3750) } },
      { {  INT16_C( 10059), -INT16_C( 19151), -INT16_C( 30928), -INT16_C(   586), -INT16_C(  3987),  INT16_C( 17896),  INT16_C( 30462), -INT16_C( 24335) },
        {  INT16_C( 10059), -INT16_C( 19151), -INT16_C( 30928), -INT16_C(   586), -INT16_C(  3987),  INT16_C( 17896),  INT16_C( 30462), -INT16_C( 24335) } },
      { { -INT16_C(  6214), -INT16_C( 27229), -INT16_C( 24158),  INT16_C( 22638), -INT16_C( 11044), -INT16_C(  6549), -INT16_C( 29906), -INT16_C( 20303) },
        { -INT16_C(  6214), -INT16_C( 27229), -INT16_C( 24158),  INT16_C( 22638), -INT16_C( 11044), -INT16_C(  6549), -INT16_C( 29906), -INT16_C( 20303) } },
      { {  INT16_C( 11305), -INT16_C( 19001), -INT16_C(  4764),  INT16_C( 16052),  INT16_C( 10294),  INT16_C(  9319),  INT16_C( 20218),  INT16_C( 13074) },
        {  INT16_C( 11305), -INT16_C( 19001), -INT16_C(  4764),  INT16_C( 16052),  INT16_C( 10294),  INT16_C(  9319),  INT16_C( 20218),  INT16_C( 13074) } },
      { {  INT16_C( 29737),  INT16_C( 12400), -INT16_C( 24056),  INT16_C( 23962), -INT16_C( 23877),  INT16_C(  1490),  INT16_C(  7748),  INT16_C( 22250) },
        {  INT16_C( 29737),  INT16_C( 12400), -INT16_C( 24056),  INT16_C( 23962), -INT16_C( 23877),  INT16_C(  1490),  INT16_C(  7748),  INT16_C( 22250) } }
    };

    for (size_t i = 0 ; i < (sizeof(test_vec) / sizeof(test_vec[0])) ; i++) {
      const int16_t* a = test_vec[i].a;
      simde_v128_t r = simde_wasm_i16x8_make(a[0], a[1], a[ 2], a[ 3], a[ 4], a[ 5], a[ 6], a[ 7]);
      simde_test_wasm_i16x8_assert_equal(r, simde_wasm_v128_load(test_vec[i].r));
    }
    return 0;
  #else
    fputc('\n', stdout);
    for (int i = 0 ; i < 8 ; i++) {
      int16_t a[sizeof(simde_v128_t) / sizeof(int16_t)];
      simde_v128_t r;

      simde_test_codegen_random_memory(sizeof(a), HEDLEY_REINTERPRET_CAST(uint8_t*, a));

      r = simde_wasm_i16x8_make(
          a[0], a[1], a[ 2], a[ 3], a[ 4], a[ 5], a[ 6], a[ 7]);

      simde_test_codegen_write_vi16(3, sizeof(a) / sizeof(a[0]), a, SIMDE_TEST_VEC_POS_FIRST);
      simde_test_wasm_i16x8_write(3, r, SIMDE_TEST_VEC_POS_LAST);
    }
    return 1;
  #endif
}

static int
test_simde_wasm_u16x8_make(SIMDE_MUNIT_TEST_ARGS) {
  #if 1
    SIMDE_TEST_STRUCT_MODIFIERS struct {
      uint16_t a[sizeof(simde_v128_t) / sizeof(uint16_t)];
      uint16_t r[sizeof(simde_v128_t) / sizeof(uint16_t)];
    } test_vec[] = {
      { { UINT16_C(20545), UINT16_C(44253), UINT16_C(58828), UINT16_C( 7787), UINT16_C(59469), UINT16_C(22081), UINT16_C(27088), UINT16_C(32395) },
        { UINT16_C(20545), UINT16_C(44253), UINT16_C(58828), UINT16_C( 7787), UINT16_C(59469), UINT16_C(22081), UINT16_C(27088), UINT16_C(32395) } },
      { { UINT16_C( 5775), UINT16_C(43217), UINT16_C(64334), UINT16_C( 1400), UINT16_C( 5738), UINT16_C(16205), UINT16_C(47899), UINT16_C( 2562) },
        { UINT16_C( 5775), UINT16_C(43217), UINT16_C(64334), UINT16_C( 1400), UINT16_C( 5738), UINT16_C(16205), UINT16_C(47899), UINT16_C( 2562) } },
      { { UINT16_C(29015), UINT16_C(25900), UINT16_C(65021), UINT16_C(47325), UINT16_C(54164), UINT16_C(53071), UINT16_C(52198), UINT16_C(13308) },
        { UINT16_C(29015), UINT16_C(25900), UINT16_C(65021), UINT16_C(47325), UINT16_C(54164), UINT16_C(53071), UINT16_C(52198), UINT16_C(13308) } },
      { { UINT16_C(62050), UINT16_C(62478), UINT16_C(55936), UINT16_C( 8179), UINT16_C(55718), UINT16_C(43779), UINT16_C(47349), UINT16_C(  892) },
        { UINT16_C(62050), UINT16_C(62478), UINT16_C(55936), UINT16_C( 8179), UINT16_C(55718), UINT16_C(43779), UINT16_C(47349), UINT16_C(  892) } },
      { { UINT16_C(44543), UINT16_C(33402), UINT16_C(50502), UINT16_C(  945), UINT16_C(  990), UINT16_C(15563), UINT16_C(27884), UINT16_C(14576) },
        { UINT16_C(44543), UINT16_C(33402), UINT16_C(50502), UINT16_C(  945), UINT16_C(  990), UINT16_C(15563), UINT16_C(27884), UINT16_C(14576) } },
      { { UINT16_C(42883), UINT16_C(43039), UINT16_C(63439), UINT16_C(58016), UINT16_C(38425), UINT16_C(63667), UINT16_C(12040), UINT16_C(30097) },
        { UINT16_C(42883), UINT16_C(43039), UINT16_C(63439), UINT16_C(58016), UINT16_C(38425), UINT16_C(63667), UINT16_C(12040), UINT16_C(30097) } },
      { { UINT16_C(57249), UINT16_C(38125), UINT16_C(59964), UINT16_C(29954), UINT16_C(59043), UINT16_C(17652), UINT16_C(28627), UINT16_C(63143) },
        { UINT16_C(57249), UINT16_C(38125), UINT16_C(59964), UINT16_C(29954), UINT16_C(59043), UINT16_C(17652), UINT16_C(28627), UINT16_C(63143) } },
      { { UINT16_C(29972), UINT16_C(31988), UINT16_C(38534), UINT16_C(25494), UINT16_C(57891), UINT16_C(22914), UINT16_C(20301), UINT16_C(10030) },
        { UINT16_C(29972), UINT16_C(31988), UINT16_C(38534), UINT16_C(25494), UINT16_C(57891), UINT16_C(22914), UINT16_C(20301), UINT16_C(10030) } }
    };

    for (size_t i = 0 ; i < (sizeof(test_vec) / sizeof(test_vec[0])) ; i++) {
      const uint16_t* a = test_vec[i].a;
      simde_v128_t r = simde_wasm_u16x8_make(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]);
      simde_test_wasm_u16x8_assert_equal(r, simde_wasm_v128_load(test_vec[i].r));
    }
    return 0;
  #else
    fputc('\n', stdout);
    for (int i = 0 ; i < 8 ; i++) {
      uint16_t a[sizeof(simde_v128_t) / sizeof(uint16_t)];
      simde_v128_t r;

      simde_test_codegen_random_memory(sizeof(a), HEDLEY_REINTERPRET_CAST(uint8_t*, a));

      r = simde_wasm_u16x8_make(
          a[0], a[1], a[ 2], a[ 3], a[ 4], a[ 5], a[ 6], a[ 7]);

      simde_test_codegen_write_vu16(3, sizeof(a) / sizeof(a[0]), a, SIMDE_TEST_VEC_POS_FIRST);
      simde_test_wasm_u16x8_write(3, r, SIMDE_TEST_VEC_POS_LAST);
    }
    return 1;
  #endif
}

static int
test_simde_wasm_i32x4_make(SIMDE_MUNIT_TEST_ARGS) {
  #if 1
    SIMDE_TEST_STRUCT_MODIFIERS struct {
      int32_t a[sizeof(simde_v128_t) / sizeof(int32_t)];
      int32_t r[sizeof(simde_v128_t) / sizeof(int32_t)];
    } test_vec[8] = {
      { {  INT32_C(   695075960),  INT32_C(  1003529442), -INT32_C(  1443479584), -INT32_C(  1440912186) },
        {  INT32_C(   695075960),  INT32_C(  1003529442), -INT32_C(  1443479584), -INT32_C(  1440912186) } },
      { {  INT32_C(   880505601),  INT32_C(   710643276),  INT32_C(  1076441165),  INT32_C(   843386822) },
        {  INT32_C(   880505601),  INT32_C(   710643276),  INT32_C(  1076441165),  INT32_C(   843386822) } },
      { { -INT32_C(  1059723883), -INT32_C(   586628182), -INT32_C(   963040743), -INT32_C(  1022298580) },
        { -INT32_C(  1059723883), -INT32_C(   586628182), -INT32_C(   963040743), -INT32_C(  1022298580) } },
      { { -INT32_C(   431686935), -INT32_C(  1520335465),  INT32_C(  1050402888), -INT32_C(   751240342) },
        { -INT32_C(   431686935), -INT32_C(  1520335465),  INT32_C(  1050402888), -INT32_C(   751240342) } },
      { { -INT32_C(   303617944),  INT32_C(        3995), -INT32_C(  1098262535),  INT32_C(   640209351) },
        { -INT32_C(   303617944),  INT32_C(        3995), -INT32_C(  1098262535),  INT32_C(   640209351) } },
      { {  INT32_C(  2027168138), -INT32_C(    56054172),  INT32_C(  1621230607),  INT32_C(  1551648147) },
        {  INT32_C(  2027168138), -INT32_C(    56054172),  INT32_C(  1621230607),  INT32_C(  1551648147) } },
      { {  INT32_C(   631271218), -INT32_C(   660470701),  INT32_C(  1855345560), -INT32_C(  1818340747) },
        {  INT32_C(   631271218), -INT32_C(   660470701),  INT32_C(  1855345560), -INT32_C(  1818340747) } },
      { {  INT32_C(   508548470),  INT32_C(  2044910564),  INT32_C(   842708574), -INT32_C(  1898864311) },
        {  INT32_C(   508548470),  INT32_C(  2044910564),  INT32_C(   842708574), -INT32_C(  1898864311) } }
    };

    for (size_t i = 0 ; i < (sizeof(test_vec) / sizeof(test_vec[0])) ; i++) {
      const int32_t* a = test_vec[i].a;
      simde_v128_t r = simde_wasm_i32x4_make(a[0], a[1], a[2], a[3]);
      simde_test_wasm_i32x4_assert_equal(r, simde_wasm_v128_load(test_vec[i].r));
    }
    return 0;
  #else
    fputc('\n', stdout);
    for (int i = 0 ; i < 8 ; i++) {
      int32_t a[sizeof(simde_v128_t) / sizeof(int32_t)];
      simde_v128_t r;

      simde_test_codegen_random_memory(sizeof(a), HEDLEY_REINTERPRET_CAST(uint8_t*, a));

      r = simde_wasm_i32x4_make(a[0], a[1], a[2], a[3]);

      simde_test_codegen_write_vi32(3, sizeof(a) / sizeof(a[0]), a, SIMDE_TEST_VEC_POS_FIRST);
      simde_test_wasm_i32x4_write(3, r, SIMDE_TEST_VEC_POS_LAST);
    }
    return 1;
  #endif
}

static int
test_simde_wasm_u32x4_make(SIMDE_MUNIT_TEST_ARGS) {
  #if 1
    SIMDE_TEST_STRUCT_MODIFIERS struct {
      uint32_t a[sizeof(simde_v128_t) / sizeof(uint32_t)];
      uint32_t r[sizeof(simde_v128_t) / sizeof(uint32_t)];
    } test_vec[] = {
      { { UINT32_C(3741086244), UINT32_C(2676567540), UINT32_C(3206543888), UINT32_C(2842651914) },
        { UINT32_C(3741086244), UINT32_C(2676567540), UINT32_C(3206543888), UINT32_C(2842651914) } },
      { { UINT32_C( 242331003), UINT32_C(4194749509), UINT32_C(3530874187), UINT32_C(1290889414) },
        { UINT32_C( 242331003), UINT32_C(4194749509), UINT32_C(3530874187), UINT32_C(1290889414) } },
      { { UINT32_C(2534034387), UINT32_C(1542371856), UINT32_C(3568344646), UINT32_C(1867754224) },
        { UINT32_C(2534034387), UINT32_C(1542371856), UINT32_C(3568344646), UINT32_C(1867754224) } },
      { { UINT32_C(3753926397), UINT32_C( 809341656), UINT32_C( 321631296), UINT32_C(1385685678) },
        { UINT32_C(3753926397), UINT32_C( 809341656), UINT32_C( 321631296), UINT32_C(1385685678) } },
      { { UINT32_C(2810325179), UINT32_C(2123103834), UINT32_C(2526303009), UINT32_C(3594571724) },
        { UINT32_C(2810325179), UINT32_C(2123103834), UINT32_C(2526303009), UINT32_C(3594571724) } },
      { { UINT32_C(3476255978), UINT32_C( 792158346), UINT32_C(3676351836), UINT32_C( 736736561) },
        { UINT32_C(3476255978), UINT32_C( 792158346), UINT32_C(3676351836), UINT32_C( 736736561) } },
      { { UINT32_C(3933311122), UINT32_C(  62652880), UINT32_C( 658539897), UINT32_C(3971227963) },
        { UINT32_C(3933311122), UINT32_C(  62652880), UINT32_C( 658539897), UINT32_C(3971227963) } },
      { { UINT32_C( 563898168), UINT32_C(2997174771), UINT32_C(3613371558), UINT32_C( 606428827) },
        { UINT32_C( 563898168), UINT32_C(2997174771), UINT32_C(3613371558), UINT32_C( 606428827) } },
    };

    for (size_t i = 0 ; i < (sizeof(test_vec) / sizeof(test_vec[0])) ; i++) {
      const uint32_t* a = test_vec[i].a;
      simde_v128_t r = simde_wasm_u32x4_make(a[0], a[1], a[2], a[3]);
      simde_test_wasm_u32x4_assert_equal(r, simde_wasm_v128_load(test_vec[i].r));
    }
    return 0;
  #else
    fputc('\n', stdout);
    for (int i = 0 ; i < 8 ; i++) {
      uint32_t a[sizeof(simde_v128_t) / sizeof(uint32_t)];
      simde_v128_t r;

      simde_test_codegen_random_memory(sizeof(a), HEDLEY_REINTERPRET_CAST(uint8_t*, a));

      r = simde_wasm_u32x4_make(a[0], a[1], a[2], a[3]);

      simde_test_codegen_write_vu32(3, sizeof(a) / sizeof(a[0]), a, SIMDE_TEST_VEC_POS_FIRST);
      simde_test_wasm_u32x4_write(3, r, SIMDE_TEST_VEC_POS_LAST);
    }
    return 1;
  #endif
}

static int
test_simde_wasm_i64x2_make(SIMDE_MUNIT_TEST_ARGS) {
  #if 1
    SIMDE_TEST_STRUCT_MODIFIERS struct {
      int64_t a[sizeof(simde_v128_t) / sizeof(int64_t)];
      int64_t r[sizeof(simde_v128_t) / sizeof(int64_t)];
    } test_vec[8] = {
      { {  INT64_C( 3706780716330511204), -INT64_C( 1056131380118931910) },
        {  INT64_C( 3706780716330511204), -INT64_C( 1056131380118931910) } },
      { {  INT64_C( 5808093934263827513),  INT64_C(  796694773045360021) },
        {  INT64_C( 5808093934263827513),  INT64_C(  796694773045360021) } },
      { { -INT64_C(  620641701573354864), -INT64_C( 6332908652456362593) },
        { -INT64_C(  620641701573354864), -INT64_C( 6332908652456362593) } },
      { { -INT64_C( 2825842138208504311),  INT64_C( 8727546400882524358) },
        { -INT64_C( 2825842138208504311),  INT64_C( 8727546400882524358) } },
      { {  INT64_C( 2466900102307772640), -INT64_C( 8981627329841242156) },
        {  INT64_C( 2466900102307772640), -INT64_C( 8981627329841242156) } },
      { {  INT64_C( 5976583528864353534),  INT64_C( 2340547445206133582) },
        {  INT64_C( 5976583528864353534),  INT64_C( 2340547445206133582) } },
      { {  INT64_C( 5001808799695288836),  INT64_C( 1033089374524368263) },
        {  INT64_C( 5001808799695288836),  INT64_C( 1033089374524368263) } },
      { {  INT64_C( 2505067071927086274),  INT64_C(  660034086012828132) },
        {  INT64_C( 2505067071927086274),  INT64_C(  660034086012828132) } }
    };

    for (size_t i = 0 ; i < (sizeof(test_vec) / sizeof(test_vec[0])) ; i++) {
      const int64_t* a = test_vec[i].a;
      simde_v128_t r = simde_wasm_i64x2_make(a[0], a[1]);
      simde_test_wasm_i64x2_assert_equal(r, simde_wasm_v128_load(test_vec[i].r));
    }
    return 0;
  #else
    fputc('\n', stdout);
    for (int i = 0 ; i < 8 ; i++) {
      int64_t a[sizeof(simde_v128_t) / sizeof(int64_t)];
      simde_v128_t r;

      simde_test_codegen_random_memory(sizeof(a), HEDLEY_REINTERPRET_CAST(uint8_t*, a));

      r = simde_wasm_i64x2_make(a[0], a[1]);

      simde_test_codegen_write_vi64(3, sizeof(a) / sizeof(a[0]), a, SIMDE_TEST_VEC_POS_FIRST);
      simde_test_wasm_i64x2_write(3, r, SIMDE_TEST_VEC_POS_LAST);
    }
    return 1;
  #endif
}

static int
test_simde_wasm_u64x2_make(SIMDE_MUNIT_TEST_ARGS) {
  #if 1
    SIMDE_TEST_STRUCT_MODIFIERS struct {
      uint64_t a[sizeof(simde_v128_t) / sizeof(uint64_t)];
      uint64_t r[sizeof(simde_v128_t) / sizeof(uint64_t)];
    } test_vec[] = {
      { { UINT64_C(13667080432244765976), UINT64_C(14908660942780224231) },
        { UINT64_C(13667080432244765976), UINT64_C(14908660942780224231) } },
      { { UINT64_C( 8032851988975082921), UINT64_C( 3524443896702843479) },
        { UINT64_C( 8032851988975082921), UINT64_C( 3524443896702843479) } },
      { { UINT64_C( 9626695193294003197), UINT64_C(13219513155173601028) },
        { UINT64_C( 9626695193294003197), UINT64_C(13219513155173601028) } },
      { { UINT64_C(  933752945553757041), UINT64_C(13014576830931341019) },
        { UINT64_C(  933752945553757041), UINT64_C(13014576830931341019) } },
      { { UINT64_C(16763484142101381182), UINT64_C(10273084224497800495) },
        { UINT64_C(16763484142101381182), UINT64_C(10273084224497800495) } },
      { { UINT64_C(10100836791351952211), UINT64_C(15952138415039713109) },
        { UINT64_C(10100836791351952211), UINT64_C(15952138415039713109) } },
      { { UINT64_C( 6493330894824637967), UINT64_C(17916322195326844626) },
        { UINT64_C( 6493330894824637967), UINT64_C(17916322195326844626) } },
      { { UINT64_C( 4331549136355019593), UINT64_C( 7743807550277083556) },
        { UINT64_C( 4331549136355019593), UINT64_C( 7743807550277083556) } },
    };

    for (size_t i = 0 ; i < (sizeof(test_vec) / sizeof(test_vec[0])) ; i++) {
      const uint64_t* a = test_vec[i].a;
      simde_v128_t r = simde_wasm_u64x2_make(a[0], a[1]);
      simde_test_wasm_u64x2_assert_equal(r, simde_wasm_v128_load(test_vec[i].r));
    }
    return 0;
  #else
    fputc('\n', stdout);
    for (int i = 0 ; i < 8 ; i++) {
      uint64_t a[sizeof(simde_v128_t) / sizeof(uint64_t)];
      simde_v128_t r;

      simde_test_codegen_random_memory(sizeof(a), HEDLEY_REINTERPRET_CAST(uint8_t*, a));

      r = simde_wasm_u64x2_make(a[0], a[1]);

      simde_test_codegen_write_vu64(3, sizeof(a) / sizeof(a[0]), a, SIMDE_TEST_VEC_POS_FIRST);
      simde_test_wasm_u64x2_write(3, r, SIMDE_TEST_VEC_POS_LAST);
    }
    return 1;
  #endif
}

static int
test_simde_wasm_f32x4_make(SIMDE_MUNIT_TEST_ARGS) {
  #if 1
    SIMDE_TEST_STRUCT_MODIFIERS struct {
      simde_float32 a[sizeof(simde_v128_t) / sizeof(simde_float32)];
      simde_float32 r[sizeof(simde_v128_t) / sizeof(simde_float32)];
    } test_vec[8] = {
      { { SIMDE_FLOAT32_C(  -591.06), SIMDE_FLOAT32_C(  -843.73), SIMDE_FLOAT32_C(   664.08), SIMDE_FLOAT32_C(  -354.19) },
        { SIMDE_FLOAT32_C(  -591.06), SIMDE_FLOAT32_C(  -843.73), SIMDE_FLOAT32_C(   664.08), SIMDE_FLOAT32_C(  -354.19) } },
      { { SIMDE_FLOAT32_C(  -257.96), SIMDE_FLOAT32_C(  -389.78), SIMDE_FLOAT32_C(  -814.91), SIMDE_FLOAT32_C(   348.16) },
        { SIMDE_FLOAT32_C(  -257.96), SIMDE_FLOAT32_C(  -389.78), SIMDE_FLOAT32_C(  -814.91), SIMDE_FLOAT32_C(   348.16) } },
      { { SIMDE_FLOAT32_C(  -599.02), SIMDE_FLOAT32_C(   786.02), SIMDE_FLOAT32_C(  -783.49), SIMDE_FLOAT32_C(  -830.91) },
        { SIMDE_FLOAT32_C(  -599.02), SIMDE_FLOAT32_C(   786.02), SIMDE_FLOAT32_C(  -783.49), SIMDE_FLOAT32_C(  -830.91) } },
      { { SIMDE_FLOAT32_C(  -437.22), SIMDE_FLOAT32_C(  -522.93), SIMDE_FLOAT32_C(   494.10), SIMDE_FLOAT32_C(  -978.90) },
        { SIMDE_FLOAT32_C(  -437.22), SIMDE_FLOAT32_C(  -522.93), SIMDE_FLOAT32_C(   494.10), SIMDE_FLOAT32_C(  -978.90) } },
      { { SIMDE_FLOAT32_C(   161.42), SIMDE_FLOAT32_C(   238.20), SIMDE_FLOAT32_C(   569.29), SIMDE_FLOAT32_C(    42.14) },
        { SIMDE_FLOAT32_C(   161.42), SIMDE_FLOAT32_C(   238.20), SIMDE_FLOAT32_C(   569.29), SIMDE_FLOAT32_C(    42.14) } },
      { { SIMDE_FLOAT32_C(   469.46), SIMDE_FLOAT32_C(   258.39), SIMDE_FLOAT32_C(   977.69), SIMDE_FLOAT32_C(  -720.50) },
        { SIMDE_FLOAT32_C(   469.46), SIMDE_FLOAT32_C(   258.39), SIMDE_FLOAT32_C(   977.69), SIMDE_FLOAT32_C(  -720.50) } },
      { { SIMDE_FLOAT32_C(  -450.95), SIMDE_FLOAT32_C(  -358.31), SIMDE_FLOAT32_C(   417.53), SIMDE_FLOAT32_C(  -526.21) },
        { SIMDE_FLOAT32_C(  -450.95), SIMDE_FLOAT32_C(  -358.31), SIMDE_FLOAT32_C(   417.53), SIMDE_FLOAT32_C(  -526.21) } },
      { { SIMDE_FLOAT32_C(   659.89), SIMDE_FLOAT32_C(  -270.01), SIMDE_FLOAT32_C(   814.76), SIMDE_FLOAT32_C(  -764.06) },
        { SIMDE_FLOAT32_C(   659.89), SIMDE_FLOAT32_C(  -270.01), SIMDE_FLOAT32_C(   814.76), SIMDE_FLOAT32_C(  -764.06) } }
    };

    for (size_t i = 0 ; i < (sizeof(test_vec) / sizeof(test_vec[0])) ; i++) {
      const simde_float32* a = test_vec[i].a;
      simde_v128_t r = simde_wasm_f32x4_make(a[0], a[1], a[2], a[3]);
      simde_test_wasm_i32x4_assert_equal(r, simde_wasm_v128_load(test_vec[i].r));
    }
    return 0;
  #else
    fputc('\n', stdout);
    for (int i = 0 ; i < 8 ; i++) {
      simde_float32 a[sizeof(simde_v128_t) / sizeof(simde_float32)];
      simde_v128_t r;

      simde_test_codegen_random_vf32(sizeof(a) / sizeof(a[0]), a, -SIMDE_FLOAT32_C(1000.0), SIMDE_FLOAT32_C(1000.0));

      r = simde_wasm_f32x4_make(a[0], a[1], a[2], a[3]);

      simde_test_codegen_write_vf32(3, sizeof(a) / sizeof(a[0]), a, SIMDE_TEST_VEC_POS_FIRST);
      simde_test_wasm_f32x4_write(3, r, SIMDE_TEST_VEC_POS_LAST);
    }
    return 1;
  #endif
}

static int
test_simde_wasm_f64x2_make(SIMDE_MUNIT_TEST_ARGS) {
  #if 1
    SIMDE_TEST_STRUCT_MODIFIERS struct {
      simde_float64 a[sizeof(simde_v128_t) / sizeof(simde_float64)];
      simde_float64 r[sizeof(simde_v128_t) / sizeof(simde_float64)];
    } test_vec[8] = {
      { { SIMDE_FLOAT64_C(   639.11), SIMDE_FLOAT64_C(  -669.03) },
        { SIMDE_FLOAT64_C(   639.11), SIMDE_FLOAT64_C(  -669.03) } },
      { { SIMDE_FLOAT64_C(  -961.87), SIMDE_FLOAT64_C(  -484.03) },
        { SIMDE_FLOAT64_C(  -961.87), SIMDE_FLOAT64_C(  -484.03) } },
      { { SIMDE_FLOAT64_C(  -709.43), SIMDE_FLOAT64_C(  -546.13) },
        { SIMDE_FLOAT64_C(  -709.43), SIMDE_FLOAT64_C(  -546.13) } },
      { { SIMDE_FLOAT64_C(  -223.36), SIMDE_FLOAT64_C(  -698.92) },
        { SIMDE_FLOAT64_C(  -223.36), SIMDE_FLOAT64_C(  -698.92) } },
      { { SIMDE_FLOAT64_C(   616.46), SIMDE_FLOAT64_C(   -91.93) },
        { SIMDE_FLOAT64_C(   616.46), SIMDE_FLOAT64_C(   -91.93) } },
      { { SIMDE_FLOAT64_C(   779.18), SIMDE_FLOAT64_C(   137.61) },
        { SIMDE_FLOAT64_C(   779.18), SIMDE_FLOAT64_C(   137.61) } },
      { { SIMDE_FLOAT64_C(  -932.66), SIMDE_FLOAT64_C(   871.28) },
        { SIMDE_FLOAT64_C(  -932.66), SIMDE_FLOAT64_C(   871.28) } },
      { { SIMDE_FLOAT64_C(   895.78), SIMDE_FLOAT64_C(  -537.98) },
        { SIMDE_FLOAT64_C(   895.78), SIMDE_FLOAT64_C(  -537.98) } }
    };

    for (size_t i = 0 ; i < (sizeof(test_vec) / sizeof(test_vec[0])) ; i++) {
      const simde_float64* a = test_vec[i].a;
      simde_v128_t r = simde_wasm_f64x2_make(a[0], a[1]);
      simde_test_wasm_i32x4_assert_equal(r, simde_wasm_v128_load(test_vec[i].r));
    }
    return 0;
  #else
    fputc('\n', stdout);
    for (int i = 0 ; i < 8 ; i++) {
      simde_float64 a[sizeof(simde_v128_t) / sizeof(simde_float64)];
      simde_v128_t r;

      simde_test_codegen_random_vf64(sizeof(a) / sizeof(a[0]), a, -SIMDE_FLOAT64_C(1000.0), SIMDE_FLOAT64_C(1000.0));

      r = simde_wasm_f64x2_make(a[0], a[1]);

      simde_test_codegen_write_vf64(3, sizeof(a) / sizeof(a[0]), a, SIMDE_TEST_VEC_POS_FIRST);
      simde_test_wasm_f64x2_write(3, r, SIMDE_TEST_VEC_POS_LAST);
    }
    return 1;
  #endif
}

SIMDE_TEST_FUNC_LIST_BEGIN
  SIMDE_TEST_FUNC_LIST_ENTRY(wasm_i8x16_make)
  SIMDE_TEST_FUNC_LIST_ENTRY(wasm_u8x16_make)
  SIMDE_TEST_FUNC_LIST_ENTRY(wasm_i16x8_make)
  SIMDE_TEST_FUNC_LIST_ENTRY(wasm_u16x8_make)
  SIMDE_TEST_FUNC_LIST_ENTRY(wasm_i32x4_make)
  SIMDE_TEST_FUNC_LIST_ENTRY(wasm_u32x4_make)
  SIMDE_TEST_FUNC_LIST_ENTRY(wasm_i64x2_make)
  SIMDE_TEST_FUNC_LIST_ENTRY(wasm_u64x2_make)
  SIMDE_TEST_FUNC_LIST_ENTRY(wasm_f32x4_make)
  SIMDE_TEST_FUNC_LIST_ENTRY(wasm_f64x2_make)
SIMDE_TEST_FUNC_LIST_END

#include "test-simd128-footer.h"
