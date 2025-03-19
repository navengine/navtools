/**
 * *binary-ops.hpp*
 *
 * =======  ========================================================================================
 * @file    include/navtools/binary-ops.hpp
 * @brief   Useful binary operations.
 * @author  Blake Baker, Daniel Sturdivant
 * @date    March 2025
 * =======  ========================================================================================
 */

#ifndef NAVTOOLS_BINARY_OPS_HPP
#define NAVTOOLS_BINARY_OPS_HPP

#include <array>
#include <cassert>
#include <cstdint>

namespace navtools {

/**
 * *=== SetBit ===*
 * @brief Set a data bit to 1
 * @param x   Number to modify
 * @param n   Position of bit to set (Position 0 is MSB and 31 is LSB by default)
 */
template <bool LsbIsZero = false>
void SetBit(uint32_t &x, const uint8_t n) {
  assert(!(n > 31));
  if constexpr (LsbIsZero) {
    x |= (0x00000001 << n);
  } else {
    x |= (0x80000000 >> n);
  }
}

/**
 * *=== UnsetBit ===*
 * @brief Set a data bit to 0
 * @param x   Number to modify
 * @param n   Position of bit to set (Position 0 is MSB and 31 is LSB by default)
 */
template <bool LsbIsZero = false>
void UnsetBit(uint32_t &x, const uint8_t n) {
  assert(!(n > 31));
  if constexpr (LsbIsZero) {
    x &= ~(0x00000001 << n);
  } else {
    x &= ~(0x80000000 >> n);
  }
}

/**
 * *=== SetBitTo ===
 * @brief Set a data bit to "b"
 * @param x   Number to modify
 * @param n   Position of bit to set (Position 0 is MSB and 31 is LSB by default)
 * @param b   New value of bit you want to modify
 */
template <bool LsbIsZero = false>
void SetBitTo(uint32_t &x, const uint8_t n, const bool b) {
  assert(!(n > 31));
  if constexpr (LsbIsZero) {
    x = (x & ~(0x00000001 << n)) | (b << n);
    // x ^= (-(uint32_t)b ^ x) & (0x00000001 << n);
  } else {
    x = (x & ~(0x80000000 >> n)) | (b << (31 - n));
    // x ^= (-(uint32_t)b ^ x) & (0x80000000 >> n);
  }
}

/**
 * *=== GetBit ===*
 * @brief Return the value of a bit
 * @param x   Number to modify
 * @param n   Position of bit to check (Position 0 is MSB and 31 is LSB by default)
 * @return Inspected bit
 */
template <bool LsbIsZero = false>
bool GetBit(const uint32_t &x, const uint8_t n) {
  assert(!(n > 31));
  if constexpr (LsbIsZero) {
    return static_cast<bool>((x >> n) & 0x00000001);
  } else {
    return static_cast<bool>(x & (0x80000000 >> n));
  }
}

/**
 * *=== GetBits ===*
 * @brief Return the value of multiple bits in series
 * @param x   Number to modify
 * @param b   Position of first bit to check (Position 0 is MSB and 31 is LSB by default)
 * @param e   Position of last bit to check (Position 0 is MSB and 31 is LSB by default)
 * @return Inspected bits shifted to the LSB end
 */
template <bool LsbIsZero = false>
uint32_t GetBits(const uint32_t &x, const uint8_t b, const uint8_t e) {
  assert(!(b > 31));
  assert(!(e > 31));
  if constexpr (LsbIsZero) {
    // TODO: figure this out later
    return x;
  } else {
    return (x >> (31 - e)) & ((1 << (e - b + 1)) - 1);
  }
}

/**
 * *=== MultiXor ===*
 * @param x   Number to modify
 * @param n   Positions of bits to check (Position 0 is MSB and 31 is LSB by default)
 * @return XOR'd number
 */
template <int Size, bool LsbIsZero = false>
bool MultiXor(const uint32_t &x, const uint8_t n[]) {
  bool r = GetBit<LsbIsZero>(x, n[0]);
  for (uint8_t i = 1; i < Size; i++) {
    r ^= GetBit<LsbIsZero>(x, n[i]);
  }
  return r;
}
template <int Size, bool LsbIsZero = false>
bool MultiXor(const uint32_t &x, const std::array<uint8_t, Size> n[Size]) {
  bool r = GetBit<LsbIsZero>(x, n[0]);
  for (uint8_t i = 1; i < Size; i++) {
    r ^= GetBit<LsbIsZero>(x, n[i]);
  }
  return r;
}
template <bool LsbIsZero = false>
bool MultiXor(const uint32_t &x, const uint8_t* n, const int Size) {
  bool r = GetBit<LsbIsZero>(x, n[0]);
  for (uint8_t i = 1; i < Size; i++) {
    r ^= GetBit<LsbIsZero>(x, n[i]);
  }
  return r;
}

/**
 * *=== TwosComp ===*
 * @brief Two's compliment to create a signed integer
 * @param x   Number to use
 * @param n   Number of bits in the integer
 * @return signed integer
 */
inline double TwosComp(uint32_t &x, const uint8_t n) {
  assert(!(n > 32));
  if (GetBit<true>(x, n - 1)) {
    uint32_t sgn_b = 0x00000001 << (n - 1);
    x &= ~sgn_b;
    return -static_cast<double>(sgn_b - x);
  } else {
    return static_cast<double>(x);
  }
  // if ((x & (1 << (n - 1))) != 0) {
  //   return static_cast<double>(x - (1 << n));
  // }
  // return static_cast<double>(x);
}

// // Obtains the value of a bit in num. The bit position is chosen with pos.
// // LsbIsZero determines whether (pos = 0) signifies the LSB or the MSB
// template <bool LsbIsZero = true>
// bool BitVal(const uint32_t& num, const uint8_t pos) {
//     assert(!(pos > 31));

//     if constexpr (LsbIsZero) {
//         return (num >> pos) & 1;
//     } else {
//         return num & (0x80000000 >> pos);
//     }
// }

// // Sets the value of a bit in num to true. The bit position is chosen with pos.
// // LsbIsZero determines whether (pos = 0) signifies the LSB or the MSB
// template <bool LsbIsZero = true>
// void BitSet(uint32_t& num, const uint8_t pos) {
//     assert(!(pos > 31));
//     if constexpr (LsbIsZero) {
//         num |= 1 << pos;
//     } else {
//         num |= (0x80000000 >> pos);
//     }
// }

// // Sets the value of a bit in num to false. The bit position is chosen with pos.
// // LsbIsZero determines whether (pos = 0) signifies the LSB or the MSB
// inline void BitUnset(uint32_t& num, const uint8_t pos) {
//     assert(!(pos > 31));
//     num &= ~(1 << pos);
// }

// // Toggles the value of a bit in num. The bit position is chosen with pos.
// // LsbIsZero determines whether (pos = 0) signifies the LSB or the MSB
// inline void BitToggle(uint32_t& num, const uint8_t pos) {
//     assert(!(pos > 31));
//     num ^= 1 << pos;
// }

// // Sets the value of a bit in num to the value val. The bit position is chosen with pos.
// // LsbIsZero determines whether (pos = 0) signifies the LSB or the MSB
// template <bool LsbIsZero = true>
// void BitEqu(uint16_t& num, const uint8_t pos, bool val) {
//     assert(!(pos > 15));
//     if constexpr (LsbIsZero) {
//         num ^= (-(uint16_t)val ^ num) & (1 << pos);
//     } else {
//         num ^= (-(uint16_t)val ^ num) & (0x8000 >> pos);
//     }
// }

// // Sets the value of a bit in num to the value val. The bit position is chosen with pos.
// // LsbIsZero determines whether (pos = 0) signifies the LSB or the MSB
// template <bool LsbIsZero = true>
// void BitEqu(uint32_t& num, const uint8_t pos, bool val) {
//     assert(!(pos > 31));
//     if constexpr (LsbIsZero) {
//         num ^= (-(uint32_t)val ^ num) & (1 << pos);
//     } else {
//         num ^= (-(uint32_t)val ^ num) & (0x80000000 >> pos);
//     }
// }

// // prints MSB to LSB
// template <bool LsbFirst = true>
// void PrintBinary(const uint8_t num) {
//     for (uint8_t i = 0; i < 8; i++) {
//         std::cout << BitVal<LsbFirst>(num, i);
//     }
//     std::cout << '\n';
// }

// template <bool LsbFirst = true>
// void PrintBinary(const uint16_t num) {
//     for (uint8_t i = 0; i < 16; i++) {
//         std::cout << BitVal<LsbFirst>(num, i);
//     }
//     std::cout << '\n';
// }

// template <bool LsbFirst = true>
// void PrintBinary(const uint32_t num) {
//     for (uint8_t i = 0; i < 32; i++) {
//         std::cout << BitVal<LsbFirst>(num, i);
//     }
//     std::cout << '\n';
// }

// // XOR a set of binary values, determined by the bit values in num. The set of positions XOR'd
// // together is given with the positions variable. LsbIsZero determines whether (pos = 0)
// signifies
// // the LSB or the MSB
// template <int Size, bool IsLsbZero = true>
// bool MultiXor(const uint32_t& num, const std::array<uint8_t, Size> positions) {
//     bool result = BitVal<IsLsbZero>(num, positions[0]);
//     for (uint8_t i = 1; i < Size; i++) {
//         result ^= BitVal<IsLsbZero>(num, positions[i]);
//     }
//     return result;
// }

// // XOR a set of binary values, determined by the bit values in num. The set of positions XOR'd
// // together is given with the positions variable. LsbIsZero determines whether (pos = 0)
// signifies
// // the LSB or the MSB
// template <int Size, bool IsLsbZero = true>
// bool MultiXor(const uint32_t& num, const uint8_t positions[]) {
//     bool result = BitVal<IsLsbZero>(num, positions[0]);
//     for (uint8_t i = 1; i < Size; i++) {
//         result ^= BitVal<IsLsbZero>(num, positions[i]);
//     }
//     return result;
// }

}  // namespace navtools
#endif
