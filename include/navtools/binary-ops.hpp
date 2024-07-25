/**
|======================================== binary-ops.hpp ==========================================|
|                                                                                                  |
|   @file     include/navtools/binary-ops.hpp                                                      |
|   @brief    Useful binary operations.                                                            |
|   @date     July 2024                                                                            |
|                                                                                                  |
|==================================================================================================|
*/

#ifndef NAVTOOLS_BINARY_OPS_HPP
#define NAVTOOLS_BINARY_OPS_HPP

#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>

namespace navtools {

// Obtains the value of a bit in num. The bit position is chosen with pos.
// LsbIsZero determines whether (pos = 0) signifies the LSB or the MSB
template <bool LsbIsZero = true>
bool BitVal(const uint32_t& num, const uint8_t pos) {
    assert(!(pos > 31));

    if constexpr (LsbIsZero) {
        return (num >> pos) & 1;
    } else {
        return num & (0x80000000 >> pos);
    }
}

// Sets the value of a bit in num to true. The bit position is chosen with pos.
// LsbIsZero determines whether (pos = 0) signifies the LSB or the MSB
template <bool LsbIsZero = true>
void BitSet(uint32_t& num, const uint8_t pos) {
    assert(!(pos > 31));
    if constexpr (LsbIsZero) {
        num |= 1 << pos;
    } else {
        num |= (0x80000000 >> pos);
    }
}

// Sets the value of a bit in num to false. The bit position is chosen with pos.
// LsbIsZero determines whether (pos = 0) signifies the LSB or the MSB
inline void BitUnset(uint32_t& num, const uint8_t pos) {
    assert(!(pos > 31));
    num &= ~(1 << pos);
}

// Toggles the value of a bit in num. The bit position is chosen with pos.
// LsbIsZero determines whether (pos = 0) signifies the LSB or the MSB
inline void BitToggle(uint32_t& num, const uint8_t pos) {
    assert(!(pos > 31));
    num ^= 1 << pos;
}

// Sets the value of a bit in num to the value val. The bit position is chosen with pos.
// LsbIsZero determines whether (pos = 0) signifies the LSB or the MSB
template <bool LsbIsZero = true>
void BitEqu(uint16_t& num, const uint8_t pos, bool val) {
    assert(!(pos > 15));
    if constexpr (LsbIsZero) {
        num ^= (-(uint16_t)val ^ num) & (1 << pos);
    } else {
        num ^= (-(uint16_t)val ^ num) & (0x8000 >> pos);
    }
}

// Sets the value of a bit in num to the value val. The bit position is chosen with pos.
// LsbIsZero determines whether (pos = 0) signifies the LSB or the MSB
template <bool LsbIsZero = true>
void BitEqu(uint32_t& num, const uint8_t pos, bool val) {
    assert(!(pos > 31));
    if constexpr (LsbIsZero) {
        num ^= (-(uint32_t)val ^ num) & (1 << pos);
    } else {
        num ^= (-(uint32_t)val ^ num) & (0x80000000 >> pos);
    }
}

// prints MSB to LSB
template <bool LsbFirst = true>
void PrintBinary(const uint8_t num) {
    for (uint8_t i = 0; i < 8; i++) {
        std::cout << BitVal<LsbFirst>(num, i);
    }
    std::cout << '\n';
}

template <bool LsbFirst = true>
void PrintBinary(const uint16_t num) {
    for (uint8_t i = 0; i < 16; i++) {
        std::cout << BitVal<LsbFirst>(num, i);
    }
    std::cout << '\n';
}

template <bool LsbFirst = true>
void PrintBinary(const uint32_t num) {
    for (uint8_t i = 0; i < 32; i++) {
        std::cout << BitVal<LsbFirst>(num, i);
    }
    std::cout << '\n';
}

// XOR a set of binary values, determined by the bit values in num. The set of positions XOR'd
// together is given with the positions variable. LsbIsZero determines whether (pos = 0) signifies
// the LSB or the MSB
template <int Size, bool IsLsbZero = true>
bool MultiXor(const uint32_t& num, const std::array<uint8_t, Size> positions) {
    bool result = BitVal<IsLsbZero>(num, positions[0]);
    for (uint8_t i = 1; i < Size; i++) {
        result ^= BitVal<IsLsbZero>(num, positions[i]);
    }
    return result;
}

// XOR a set of binary values, determined by the bit values in num. The set of positions XOR'd
// together is given with the positions variable. LsbIsZero determines whether (pos = 0) signifies
// the LSB or the MSB
template <int Size, bool IsLsbZero = true>
bool MultiXor(const uint32_t& num, const uint8_t positions[]) {
    bool result = BitVal<IsLsbZero>(num, positions[0]);
    for (uint8_t i = 1; i < Size; i++) {
        result ^= BitVal<IsLsbZero>(num, positions[i]);
    }
    return result;
}

}  // namespace navtools
#endif
