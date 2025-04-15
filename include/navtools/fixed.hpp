#ifndef NAVTOOLS_FIXED_HPP
#define NAVTOOLS_FIXED_HPP

#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <type_traits>
#include <cassert>
#include <array>
#include <vector>
#include <algorithm>
#include <cctype>
#include <climits>
#include <ios>


namespace navtools
{
// This is a modified version of the fixed class from https://github.com/MikeLankamp/navtools

//! Fixed-point number type
//! \tparam BaseType         the base integer type used to store the fixed-point number. This can be a signed or unsigned type.
//! \tparam IntermediateType the integer type used to store intermediate results during calculations.
//! \tparam FractionBits     the number of bits of the BaseType used to store the fraction
//! \tparam EnableRounding   enable rounding of LSB for multiplication, division, and type conversion
template <typename BaseType, typename IntermediateType, unsigned int FractionBits, bool EnableRounding = true>
class fixed
{
    static_assert(std::is_integral<BaseType>::value, "BaseType must be an integral type");
    static_assert(FractionBits > 0, "FractionBits must be greater than zero");
    static_assert(FractionBits <= sizeof(BaseType) * 8 - 1, "BaseType must at least be able to contain entire fraction, with space for at least one integral bit");
    static_assert((sizeof(IntermediateType) > sizeof(BaseType)) || (sizeof(BaseType)*8 - FractionBits > 1), "IntermediateType must be larger than BaseType");
    static_assert(std::is_signed<IntermediateType>::value == std::is_signed<BaseType>::value, "IntermediateType must have same signedness as BaseType");

    // Although this value fits in the BaseType in terms of bits, if there's only one integral bit, this value
    // is incorrect (flips from positive to negative), so we must extend the size to IntermediateType.
    static constexpr IntermediateType FRACTION_MULT = IntermediateType(1) << FractionBits;
    static constexpr bool DISABLE_MULT = (sizeof(BaseType) == sizeof(IntermediateType));

    struct raw_construct_tag {};
    constexpr inline fixed(BaseType val, raw_construct_tag) noexcept : m_value(val) {}

public:
    inline fixed() noexcept = default;

    // Converts an integral number to the fixed-point type.
    // Like static_cast, this truncates bits that don't fit.
    template <typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
    constexpr inline explicit fixed(T val) noexcept
        : m_value(static_cast<BaseType>(val * FRACTION_MULT))
    {}

    // Converts an floating-point number to the fixed-point type.
    // Like static_cast, this truncates bits that don't fit.
    template <typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
    constexpr inline explicit fixed(T val) noexcept
        : m_value(static_cast<BaseType>((EnableRounding) ?
		       (val >= 0.0) ? (val * FRACTION_MULT + T{0.5}) : (val * FRACTION_MULT - T{0.5})
		      : (val * FRACTION_MULT)))
    {}

    // Constructs from another fixed-point type with possibly different underlying representation.
    // Like static_cast, this truncates bits that don't fit.
    template <typename B, typename I, unsigned int F, bool R>
    constexpr inline explicit fixed(fixed<B,I,F,R> val) noexcept
        : m_value(from_fixed_point<F>(val.raw_value()).raw_value())
    {}

    // Explicit conversion to a floating-point type
    template <typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
    constexpr inline explicit operator T() const noexcept
    {
      return static_cast<T>(m_value) / FRACTION_MULT;
    }

    // Explicit conversion to an integral type
    template <typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
    constexpr inline explicit operator T() const noexcept
    {
      return static_cast<T>(m_value / FRACTION_MULT);
    }

    // Returns the raw underlying value of this type.
    // Do not use this unless you know what you're doing.
    constexpr inline BaseType raw_value() const noexcept
    {
      return m_value;
    }

    //! Constructs a fixed-point number from another fixed-point number.
    //! \tparam NumFractionBits the number of bits used by the fraction in \a value.
    //! \param value the integer fixed-point number
    template <unsigned int NumFractionBits, typename T, typename std::enable_if<(NumFractionBits > FractionBits)>::type* = nullptr>
    static constexpr inline fixed from_fixed_point(T value) noexcept
    {
      // To correctly round the last bit in the result, we need one more bit of information.
      // We do this by multiplying by two before dividing and adding the LSB to the real result.
      return (EnableRounding) ? fixed(static_cast<BaseType>(
             value / (T(1) << (NumFractionBits - FractionBits)) +
            (value / (T(1) << (NumFractionBits - FractionBits - 1)) % 2)),
            raw_construct_tag{}) :
            fixed(static_cast<BaseType>(value / (T(1) << (NumFractionBits - FractionBits))),
            raw_construct_tag{});
    }

    template <unsigned int NumFractionBits, typename T, typename std::enable_if<(NumFractionBits <= FractionBits)>::type* = nullptr>
    static constexpr inline fixed from_fixed_point(T value) noexcept
    {
      return fixed(static_cast<BaseType>(
          value * (T(1) << (FractionBits - NumFractionBits))),
          raw_construct_tag{});
    }

    // Constructs a fixed-point number from its raw underlying value.
    // Do not use this unless you know what you're doing.
    static constexpr inline fixed from_raw_value(BaseType value) noexcept
    {
      return fixed(value, raw_construct_tag{});
    }

    //
    // Constants
    //
    static constexpr fixed e() { return from_fixed_point<61>(6267931151224907085ll); }
    static constexpr fixed pi() { return from_fixed_point<61>(7244019458077122842ll); }
    static constexpr fixed half_pi() { return from_fixed_point<62>(7244019458077122842ll); }
    static constexpr fixed two_pi() { return from_fixed_point<60>(7244019458077122842ll); }

    //
    // Arithmetic member operators
    //

    constexpr inline fixed operator-() const noexcept
    {
      return fixed::from_raw_value(-m_value);
    }

    inline fixed& operator+=(const fixed& y) noexcept
    {
      m_value += y.m_value;
      return *this;
    }

    template <typename I, typename std::enable_if<std::is_integral<I>::value>::type* = nullptr>
    inline fixed& operator+=(I y) noexcept
    {
      m_value += y * FRACTION_MULT;
      return *this;
    }

    inline fixed& operator-=(const fixed& y) noexcept
    {
      m_value -= y.m_value;
      return *this;
    }

    template <typename I, typename std::enable_if<std::is_integral<I>::value>::type* = nullptr>
    inline fixed& operator-=(I y) noexcept
    {
      m_value -= y * FRACTION_MULT;
      return *this;
    }


    inline fixed& operator*=(const fixed& y) noexcept
    {
      static_assert(!DISABLE_MULT, "This type does not support multiplication since the base and intermediate types are the same size.");
      if constexpr (EnableRounding){
        // Normal fixed-point multiplication is: x * y / 2**FractionBits.
        // To correctly round the last bit in the result, we need one more bit of information.
        // We do this by multiplying by two before dividing and adding the LSB to the real result.
        auto value = (static_cast<IntermediateType>(m_value) * y.m_value) / (FRACTION_MULT / 2);
        m_value = static_cast<BaseType>((value / 2) + (value % 2));
      } else {
        auto value = (static_cast<IntermediateType>(m_value) * y.m_value) / FRACTION_MULT;
        m_value = static_cast<BaseType>(value);
      }
      return *this;
    }

    template <typename I, typename std::enable_if<std::is_integral<I>::value>::type* = nullptr>
    inline fixed& operator*=(I y) noexcept
    {
      static_assert(!DISABLE_MULT, "This type does not support multiplication since the base and intermediate types are the same size.");
      m_value *= y;
      return *this;
    }

    inline fixed& operator/=(const fixed& y) noexcept
    {
      static_assert(!DISABLE_MULT, "This type does not support division since the base and intermediate types are the same size.");
      assert(y.m_value != 0);
      if constexpr (EnableRounding){
        // Normal fixed-point division is: x * 2**FractionBits / y.
        // To correctly round the last bit in the result, we need one more bit of information.
        // We do this by multiplying by two before dividing and adding the LSB to the real result.
        auto value = (static_cast<IntermediateType>(m_value) * FRACTION_MULT * 2) / y.m_value;
        m_value = static_cast<BaseType>((value / 2) + (value % 2));
      } else {
        auto value = (static_cast<IntermediateType>(m_value) * FRACTION_MULT) / y.m_value;
        m_value = static_cast<BaseType>(value);
      }
      return *this;
    }

    template <typename I, typename std::enable_if<std::is_integral<I>::value>::type* = nullptr>
    inline fixed& operator/=(I y) noexcept
    {
      static_assert(!DISABLE_MULT, "This type does not support division since the base and intermediate types are the same size.");
      m_value /= y;
      return *this;
    }

private:
    BaseType m_value;
};

//
// Convenience typedefs
//
using fixed_16_16 = fixed<std::int32_t, std::int64_t, 16>;
using fixed_24_8 = fixed<std::int32_t, std::int64_t, 8>;
using fixed_8_24 = fixed<std::int32_t, std::int64_t, 24>;

//
// Addition
//

template <typename B, typename I, unsigned int F, bool R>
constexpr inline fixed<B, I, F, R> operator+(const fixed<B, I, F, R>& x, const fixed<B, I, F, R>& y) noexcept
{
  return fixed<B, I, F, R>(x) += y;
}

template <typename B, typename I, unsigned int F, bool R, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
constexpr inline fixed<B, I, F, R> operator+(const fixed<B, I, F, R>& x, T y) noexcept
{
  return fixed<B, I, F, R>(x) += y;
}

template <typename B, typename I, unsigned int F, bool R, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
constexpr inline fixed<B, I, F, R> operator+(T x, const fixed<B, I, F, R>& y) noexcept
{
  return fixed<B, I, F, R>(y) += x;
}

//
// Subtraction
//

template <typename B, typename I, unsigned int F, bool R>
constexpr inline fixed<B, I, F, R> operator-(const fixed<B, I, F, R>& x, const fixed<B, I, F, R>& y) noexcept
{
  return fixed<B, I, F, R>(x) -= y;
}

template <typename B, typename I, unsigned int F, bool R, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
constexpr inline fixed<B, I, F, R> operator-(const fixed<B, I, F, R>& x, T y) noexcept
{
  return fixed<B, I, F, R>(x) -= y;
}

template <typename B, typename I, unsigned int F, bool R, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
constexpr inline fixed<B, I, F, R> operator-(T x, const fixed<B, I, F, R>& y) noexcept
{
  return fixed<B, I, F, R>(x) -= y;
}


//
// Multiplication
//

template <typename B, typename I, unsigned int F, bool R>
constexpr inline fixed<B, I, F, R> operator*(const fixed<B, I, F, R>& x, const fixed<B, I, F, R>& y) noexcept
{
  static_assert(sizeof(B)==sizeof(I), "This fixed type does not support multiplication since the base and intermediate types are the same size.");
  return fixed<B, I, F, R>(x) *= y;
}

template <typename B, typename I, unsigned int F, bool R, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
constexpr inline fixed<B, I, F, R> operator*(const fixed<B, I, F, R>& x, T y) noexcept
{
  static_assert(sizeof(B)==sizeof(I), "This fixed type does not support multiplication since the base and intermediate types are the same size.");
  return fixed<B, I, F, R>(x) *= y;
}

template <typename B, typename I, unsigned int F, bool R, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
constexpr inline fixed<B, I, F, R> operator*(T x, const fixed<B, I, F, R>& y) noexcept
{
  static_assert(sizeof(B)==sizeof(I), "This fixed type does not support multiplication since the base and intermediate types are the same size.");
  return fixed<B, I, F, R>(y) *= x;
}

//
// Division
//

template <typename B, typename I, unsigned int F, bool R>
constexpr inline fixed<B, I, F, R> operator/(const fixed<B, I, F, R>& x, const fixed<B, I, F, R>& y) noexcept
{
  static_assert(sizeof(B)==sizeof(I), "This fixed type does not support division since the base and intermediate types are the same size.");
  return fixed<B, I, F, R>(x) /= y;
}

template <typename B, typename I, unsigned int F, typename T, bool R, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
constexpr inline fixed<B, I, F, R> operator/(const fixed<B, I, F, R>& x, T y) noexcept
{
  static_assert(sizeof(B)==sizeof(I), "This fixed type does not support division since the base and intermediate types are the same size.");
  return fixed<B, I, F, R>(x) /= y;
}

template <typename B, typename I, unsigned int F, typename T, bool R, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
constexpr inline fixed<B, I, F, R> operator/(T x, const fixed<B, I, F, R>& y) noexcept
{
  static_assert(sizeof(B)==sizeof(I), "This fixed type does not support division since the base and intermediate types are the same size.");
  return fixed<B, I, F, R>(x) /= y;
}


//
// Comparison operators
//

template <typename B, typename I, unsigned int F, bool R>
constexpr inline bool operator==(const fixed<B, I, F, R>& x, const fixed<B, I, F, R>& y) noexcept
{
  return x.raw_value() == y.raw_value();
}

template <typename B, typename I, unsigned int F, bool R>
constexpr inline bool operator!=(const fixed<B, I, F, R>& x, const fixed<B, I, F, R>& y) noexcept
{
  return x.raw_value() != y.raw_value();
}

template <typename B, typename I, unsigned int F, bool R>
constexpr inline bool operator<(const fixed<B, I, F, R>& x, const fixed<B, I, F, R>& y) noexcept
{
  return x.raw_value() < y.raw_value();
}

template <typename B, typename I, unsigned int F, bool R>
constexpr inline bool operator>(const fixed<B, I, F, R>& x, const fixed<B, I, F, R>& y) noexcept
{
  return x.raw_value() > y.raw_value();
}

template <typename B, typename I, unsigned int F, bool R>
constexpr inline bool operator<=(const fixed<B, I, F, R>& x, const fixed<B, I, F, R>& y) noexcept
{
  return x.raw_value() <= y.raw_value();
}

template <typename B, typename I, unsigned int F, bool R>
constexpr inline bool operator>=(const fixed<B, I, F, R>& x, const fixed<B, I, F, R>& y) noexcept
{
  return x.raw_value() >= y.raw_value();
}


// ------------------------------ MATH FUNCTIONS ------------------------------
template <typename B, typename I, unsigned int F, bool R>
inline fixed<B, I, F, R> ceil(fixed<B, I, F, R> x) noexcept
{
    constexpr auto FRAC = B(1) << F;
    auto value = x.raw_value();
    if (value > 0) value += FRAC - 1;
    return fixed<B, I, F, R>::from_raw_value(value / FRAC * FRAC);
}

template <typename B, typename I, unsigned int F, bool R>
inline fixed<B, I, F, R> floor(fixed<B, I, F, R> x) noexcept
{
    constexpr auto FRAC = B(1) << F;
    auto value = x.raw_value();
    if (value < 0) value -= FRAC - 1;
    return fixed<B, I, F, R>::from_raw_value(value / FRAC * FRAC);
}

template <typename B, typename I, unsigned int F, bool R>
inline fixed<B, I, F, R> round(fixed<B, I, F, R> x) noexcept
{
    constexpr auto FRAC = B(1) << F;
    auto value = x.raw_value() / (FRAC / 2);
    return fixed<B, I, F, R>::from_raw_value(((value / 2) + (value % 2)) * FRAC);
}

template <typename B, typename I, unsigned int F, bool R>
constexpr inline fixed<B, I, F, R> CircMod(fixed<B, I, F, R>& x, fixed<B, I, F, R> y) 
{
  assert(y.raw_value() > 0);
  x =
    fixed<B, I, F, R>::from_raw_value(
      (x.raw_value() >= 0) ? (x.raw_value() % y.raw_value())
      : ((x.raw_value() + 1) % y.raw_value() + (y.raw_value() - 1))
    );
}

template <typename B, typename I, unsigned int F, bool R, typename C, typename J, unsigned int G, bool S>
constexpr inline fixed<B, I, F, R> copysign(fixed<B, I, F, R> x, fixed<C, J, G, S> y) noexcept
{
    return
        x = abs(x),
        (y >= fixed<C, J, G, S>{0}) ? x : -x;
}



namespace detail
{
// Number of base-10 digits required to fully represent a number of bits
static constexpr int max_digits10(int bits)
{
  // 8.24 fixed-point equivalent of (int)ceil(bits * std::log10(2));
  using T = long long;
  return static_cast<int>((T{bits} * 5050445 + (T{1} << 24) - 1) >> 24);
}

// Number of base-10 digits that can be fully represented by a number of bits
static constexpr int digits10(int bits)
{
  // 8.24 fixed-point equivalent of (int)(bits * std::log10(2));
  using T = long long;
  return static_cast<int>((T{bits} * 5050445) >> 24);
}

// Returns the index of the most-signifcant set bit
inline long find_highest_bit(unsigned long long value) noexcept
{
    assert(value != 0);
#if defined(_MSC_VER)
    unsigned long index;
#if defined(_WIN64)
    _BitScanReverse64(&index, value);
#else
    if (_BitScanReverse(&index, static_cast<unsigned long>(value >> 32)) != 0) {
        index += 32;
    } else {
        _BitScanReverse(&index, static_cast<unsigned long>(value & 0xfffffffflu));
    }
#endif
    return index;
#elif defined(__GNUC__) || defined(__clang__)
    return sizeof(value) * 8 - 1 - __builtin_clzll(value);
#else
#   error "your platform does not support find_highest_bit()"
#endif
}

} // namespace detail
} // namespace navtools

// Specializations for customization points
namespace std
{

template <typename B, typename I, unsigned int F, bool R>
struct hash<navtools::fixed<B,I,F,R>>
{
  using argument_type = navtools::fixed<B, I, F, R>;
  using result_type = std::size_t;

  result_type operator()(argument_type arg) const noexcept(noexcept(std::declval<std::hash<B>>()(arg.raw_value()))) {
    return m_hash(arg.raw_value());
  }

private:
  std::hash<B> m_hash;
};

template <typename B, typename I, unsigned int F, bool R>
struct numeric_limits<navtools::fixed<B,I,F,R>>
{
  static constexpr bool is_specialized = true;
  static constexpr bool is_signed = std::numeric_limits<B>::is_signed;
  static constexpr bool is_integer = false;
  static constexpr bool is_exact = true;
  static constexpr bool has_infinity = false;
  static constexpr bool has_quiet_NaN = false;
  static constexpr bool has_signaling_NaN = false;
  static constexpr std::float_denorm_style has_denorm = std::denorm_absent;
  static constexpr bool has_denorm_loss = false;
  static constexpr std::float_round_style round_style = std::round_to_nearest;
  static constexpr bool is_iec559 = false;
  static constexpr bool is_bounded = true;
  static constexpr bool is_modulo = std::numeric_limits<B>::is_modulo;
  static constexpr int digits = std::numeric_limits<B>::digits;

  // Any number with `digits10` significant base-10 digits (that fits in
  // the range of the type) is guaranteed to be convertible from text and
  // back without change. Worst case, this is 0.000...001, so we can only
  // guarantee this case. Nothing more.
  static constexpr int digits10 = 1;

  // This is equal to max_digits10 for the integer and fractional part together.
  static constexpr int max_digits10 =
      navtools::detail::max_digits10(std::numeric_limits<B>::digits - F) + navtools::detail::max_digits10(F);

  static constexpr int radix = 2;
  static constexpr int min_exponent = 1 - F;
  static constexpr int min_exponent10 = -navtools::detail::digits10(F);
  static constexpr int max_exponent = std::numeric_limits<B>::digits - F;
  static constexpr int max_exponent10 = navtools::detail::digits10(std::numeric_limits<B>::digits - F);
  static constexpr bool traps = true;
  static constexpr bool tinyness_before = false;

  static constexpr navtools::fixed<B,I,F,R> lowest() noexcept {
    return navtools::fixed<B,I,F,R>::from_raw_value(std::numeric_limits<B>::lowest());
  };

  static constexpr navtools::fixed<B,I,F,R> min() noexcept {
    return lowest();
  }

  static constexpr navtools::fixed<B,I,F,R> max() noexcept {
    return navtools::fixed<B,I,F,R>::from_raw_value(std::numeric_limits<B>::max());
  };

  static constexpr navtools::fixed<B,I,F,R> epsilon() noexcept {
    return navtools::fixed<B,I,F,R>::from_raw_value(1);
  };

  static constexpr navtools::fixed<B,I,F,R> round_error() noexcept {
    return navtools::fixed<B,I,F,R>(1) / 2;
  };

  static constexpr navtools::fixed<B,I,F,R> denorm_min() noexcept {
    return min();
  }
};

template <typename B, typename I, unsigned int F, bool R>
constexpr bool numeric_limits<navtools::fixed<B,I,F,R>>::is_specialized;
template <typename B, typename I, unsigned int F, bool R>
constexpr bool numeric_limits<navtools::fixed<B,I,F,R>>::is_signed;
template <typename B, typename I, unsigned int F, bool R>
constexpr bool numeric_limits<navtools::fixed<B,I,F,R>>::is_integer;
template <typename B, typename I, unsigned int F, bool R>
constexpr bool numeric_limits<navtools::fixed<B,I,F,R>>::is_exact;
template <typename B, typename I, unsigned int F, bool R>
constexpr bool numeric_limits<navtools::fixed<B,I,F,R>>::has_infinity;
template <typename B, typename I, unsigned int F, bool R>
constexpr bool numeric_limits<navtools::fixed<B,I,F,R>>::has_quiet_NaN;
template <typename B, typename I, unsigned int F, bool R>
constexpr bool numeric_limits<navtools::fixed<B,I,F,R>>::has_signaling_NaN;
template <typename B, typename I, unsigned int F, bool R>
constexpr std::float_denorm_style numeric_limits<navtools::fixed<B,I,F,R>>::has_denorm;
template <typename B, typename I, unsigned int F, bool R>
constexpr bool numeric_limits<navtools::fixed<B,I,F,R>>::has_denorm_loss;
template <typename B, typename I, unsigned int F, bool R>
constexpr std::float_round_style numeric_limits<navtools::fixed<B,I,F,R>>::round_style;
template <typename B, typename I, unsigned int F, bool R>
constexpr bool numeric_limits<navtools::fixed<B,I,F,R>>::is_iec559;
template <typename B, typename I, unsigned int F, bool R>
constexpr bool numeric_limits<navtools::fixed<B,I,F,R>>::is_bounded;
template <typename B, typename I, unsigned int F, bool R>
constexpr bool numeric_limits<navtools::fixed<B,I,F,R>>::is_modulo;
template <typename B, typename I, unsigned int F, bool R>
constexpr int numeric_limits<navtools::fixed<B,I,F,R>>::digits;
template <typename B, typename I, unsigned int F, bool R>
constexpr int numeric_limits<navtools::fixed<B,I,F,R>>::digits10;
template <typename B, typename I, unsigned int F, bool R>
constexpr int numeric_limits<navtools::fixed<B,I,F,R>>::max_digits10;
template <typename B, typename I, unsigned int F, bool R>
constexpr int numeric_limits<navtools::fixed<B,I,F,R>>::radix;
template <typename B, typename I, unsigned int F, bool R>
constexpr int numeric_limits<navtools::fixed<B,I,F,R>>::min_exponent;
template <typename B, typename I, unsigned int F, bool R>
constexpr int numeric_limits<navtools::fixed<B,I,F,R>>::min_exponent10;
template <typename B, typename I, unsigned int F, bool R>
constexpr int numeric_limits<navtools::fixed<B,I,F,R>>::max_exponent;
template <typename B, typename I, unsigned int F, bool R>
constexpr int numeric_limits<navtools::fixed<B,I,F,R>>::max_exponent10;
template <typename B, typename I, unsigned int F, bool R>
constexpr bool numeric_limits<navtools::fixed<B,I,F,R>>::traps;
template <typename B, typename I, unsigned int F, bool R>
constexpr bool numeric_limits<navtools::fixed<B,I,F,R>>::tinyness_before;

}

namespace navtools {

template<typename T>
struct is_fixed : std::false_type {};

template<typename BaseType, typename IntermediateType, unsigned int FractionBits, bool EnableRounding>
struct is_fixed<fixed<BaseType, IntermediateType, FractionBits, EnableRounding>> : std::true_type {};

#if  __cplusplus >= 201703L
template<typename T>
inline constexpr bool is_fixed_v = is_fixed<T>::value;
#endif

typedef navtools::fixed<int64_t,int64_t,42> Fixed42;


// ---------------------------- IOS IMPLEMENTATION ----------------------------
template <typename CharT, typename B, typename I, unsigned int F, bool R>
std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os, fixed<B, I, F, R> x) noexcept
{
    const auto uppercase = ((os.flags() & std::ios_base::uppercase) != 0);
    const auto showpoint = ((os.flags() & std::ios_base::showpoint) != 0);
    const auto adjustfield = (os.flags() & std::ios_base::adjustfield);
    const auto width = os.width();
    const auto& ctype = std::use_facet<std::ctype<CharT>>(os.getloc());
    const auto& numpunct = std::use_facet<std::numpunct<CharT>>(os.getloc());

    auto floatfield = (os.flags() & std::ios_base::floatfield);
    auto precision = os.precision();
    auto show_trailing_zeros = true;
    auto use_significant_digits = false;

    // Invalid precision? Reset to the default
    if (precision < 0)
    {
        precision = 6;
    }

    // Output buffer. Needs to be big enough for the formatted number without padding.
    // Optional prefixes (i.e. "+"/"-", decimal separator, exponent "e+/-" and/or "0x").
    constexpr auto worst_case_constant_size = 6;
    // Maximum number of digits from the base type (covers integral + fractional digits)
    constexpr auto worst_case_digit_count = std::numeric_limits<B>::digits10 + 2;
    // Exponent suffixes (i.e. maximum digits based on log of the base type size).
    // Needs a log10, but that isn't constexpr, so we're over-allocating on the stack. Can't hurt.
    constexpr auto worst_case_suffix_size = std::numeric_limits<B>::digits;
    // Double the digit count: in the worst case the thousands grouping add a character per digit.
    using buffer_t = std::array<CharT, worst_case_constant_size + worst_case_digit_count * 2 + worst_case_suffix_size>;
    buffer_t buffer;

    // Output cursor
    auto end = buffer.begin();

    // Keep track of the start of "internal" padding
    typename buffer_t::iterator internal_pad = buffer.end();

    // Representation of a number.
    // The value of the number is: raw / divisor * (10|2) ^ exponent
    // The base of the exponent is 2 in hexfloat mode, or 10 otherwise.
    struct number_t {
        I raw;          // raw fixed-point value
        I divisor;      // the divisor indicating the place of the decimal point
        int exponent;   // the exponent applied
    };

    // Convert a value without exponent to scientific representation
    // where the part before the decimal point is less than 10.
    const auto as_scientific = [](number_t value) {
        assert(value.exponent == 0);
        if (value.raw > 0)
        {
            while (value.raw / 10 >= value.divisor) {
                value.divisor *= 10;
                ++value.exponent;
            }
            while (value.raw < value.divisor) {
                 value.raw *= 10;
                --value.exponent;
            }
        }
        return value;
    };

    number_t value = { x.raw_value(), I{1} << F, 0};

    auto base = B{10};

    // First write the sign
    if (value.raw < 0)
    {
        *end++ = ctype.widen('-');
        value.raw = -value.raw;
        internal_pad = end;
    }
    else if (os.flags() & std::ios_base::showpos)
    {
        *end++ = ctype.widen('+');
        internal_pad = end;
    }
    assert(value.raw >= 0);

    switch (floatfield)
    {
    case std::ios_base::fixed | std::ios_base::scientific:
        // Hexadecimal mode: figure out the hexadecimal exponent and write "0x"
        if (value.raw > 0)
        {
            auto bit  = detail::find_highest_bit(value.raw);
            value.exponent = bit - F;    // exponent is applied to base 2
            value.divisor = I{1} << bit; // divisor is at the highest bit, ensuring it starts with "1."
            precision = (bit + 3) / 4;   // precision is number of nibbles, so we show all of them
        }
        base = 16;
        show_trailing_zeros = false; // Always strip trailing zeros in hexfloat mode

        *end++ = ctype.widen('0');
        *end++ = ctype.widen(uppercase ? 'X' : 'x');
        break;

    case std::ios_base::scientific:
        // Scientific mode, normalize value to scientific notation
        value = as_scientific(value);
        break;

    case std::ios_base::fixed:
        // Fixed mode. Nothing to do.
        break;

    default:
    {
        // "auto" mode: figure out the exponent
        const number_t sci_value = as_scientific(value);

        // Now `precision` indicates the number of *significant digits* (not fractional digits).
        use_significant_digits = true;
        precision = std::max<std::streamsize>(precision, 1);

        if (sci_value.exponent >= precision || sci_value.exponent < -4) {
            // Display as scientific format
            floatfield = std::ios_base::scientific;
            value = sci_value;
        } else {
            // Display as fixed format.
            // "showpoint" indicates whether or not we show trailing zeros
            floatfield = std::ios_base::fixed;
            show_trailing_zeros = showpoint;
        }
        break;
    }
    };

    // If we didn't write a sign, any internal padding starts here
    // (after a potential "0x" for hexfloats).
    if (internal_pad == buffer.end()) {
        internal_pad = end;
    }

    // Separate out the integral part of the number
    I integral = value.raw / value.divisor;
    value.raw %= value.divisor;

    // Here we start printing the number itself
    const char* const digits = uppercase ? "0123456789ABCDEF" : "0123456789abcdef";
    const auto digits_start = end;

    // Are we already printing significant digits? (yes if we're not counting significant digits)
    bool significant_digits = !use_significant_digits;

    // Print the integral part
    int last_digit = 0;
    if (integral == 0) {
        *end++ = ctype.widen('0');
        if (value.raw == 0) {
            // If the fraction is zero too, all zeros including the integral count
            // as significant digits.
            significant_digits = true;
        }
    } else {
        while (integral > 0) {
            last_digit = integral % base;
            *end++ = ctype.widen(digits[last_digit]);
            integral /= base;
        }
        std::reverse(digits_start, end);
        significant_digits = true;
    }

    if (use_significant_digits && significant_digits)
    {
        // Apparently the integral part was significant; subtract its
        // length from the remaining significant digits.
        precision -= (end - digits_start);
    }

    // At this point, `value` contains only the fraction and
    // `precision` holds the number of digits to print.
    assert(value.raw < value.divisor);
    assert(precision >= 0);

    // Location of decimal point
    typename buffer_t::iterator point = buffer.end();

    // Start (and length) of the trailing zeros to insert while printing
    // By tracking this to print them later instead of actually printing them now,
    // we can support large precisions with a small printing buffer.
    typename buffer_t::iterator trailing_zeros_start = buffer.end();
    std::streamsize trailing_zeros_count = 0;

    if (precision > 0)
    {
        // Print the fractional part
        *(point = end++) = numpunct.decimal_point();

        for (int i = 0; i < precision; ++i)
        {
            if (value.raw == 0)
            {
                // The rest of the digits are all zeros, mark them
                // to be printed in this spot.
                trailing_zeros_start = end;
                trailing_zeros_count = precision - i;
                break;
            }

            // Shift the divisor if we can to avoid overflow on the value
            if (value.divisor % base == 0) {
                value.divisor /= base;
            } else {
                value.raw *= base;
            }
            assert(value.divisor > 0);
            assert(value.raw >= 0);
            last_digit = (value.raw / value.divisor) % base;
            value.raw %= value.divisor;
            *end++ = ctype.widen(digits[last_digit]);

            if (!significant_digits) {
                // We're still finding the first significant digit
                if (last_digit != 0) {
                    // Found it
                    significant_digits = true;
                } else {
                    // Not yet; increment number of digits to print
                    ++precision;
                }
            }
        }
    }
    else if (showpoint)
    {
        // No fractional part to print, but we still want the point
        *(point = end++) = numpunct.decimal_point();
    }

    // Insert `ch` into the output at `position`, updating all references accordingly
    const auto insert_character = [&](typename buffer_t::iterator position, CharT ch) {
        assert(position >= buffer.begin() && position < end);
        std::move_backward(position, end, end + 1);
        if (point != buffer.end() && position < point) {
            ++point;
        }
        if (trailing_zeros_start != buffer.end() && position < trailing_zeros_start) {
            ++trailing_zeros_start;
        }
        ++end;
        *position = ch;
    };

    // Round the number: round to nearest
    bool increment = false;
    if (value.raw > value.divisor / 2) {
        // Round up
        increment = true;
    } else if (value.raw == value.divisor / 2) {
        // It's a tie (i.e. "xyzw.5"): round to even
        increment = ((last_digit % 2) == 1);
    }

    if (increment)
    {
        auto p = end - 1;
        // Increment all digits backwards while we see "9"
        while (p >= digits_start) {
            if (p == point) {
                // Skip over the decimal point
                --p;
            }
            if ((*p)++ != ctype.widen('9')) {
                break;
            }
            *p-- = ctype.widen('0');
        }

        if (p < digits_start) {
            // We've incremented all the way to the start (all 9's), we need to insert the
            // carried-over 1 from incrementing the last 9.
            assert(p == digits_start - 1);
            insert_character(++p, ctype.widen('1'));

            if (floatfield == std::ios::scientific)
            {
                // We just made the integral part equal to 10, so we shift the decimal point
                // back one place (if any) and tweak the exponent, so that we keep the integer part
                // less than 10.
                if (point != buffer.end()) {
                    assert(p + 2 == point);
                    std::swap(*(point - 1), *point);
                    --point;
                }
                ++value.exponent;

                // We've introduced an extra digit so we need to strip the last digit
                // to maintain the same precision
                --end;
            }
        }

        if (use_significant_digits && *p == ctype.widen('1') && point != buffer.end()) {
            // We've converted a leading zero to a 1 so we need to strip the last digit
            // (behind the decimal point) to maintain the same significant digit count.
            --end;
        }
    }

    if (point != buffer.end())
    {
        if (!show_trailing_zeros)
        {
            // Remove trailing zeros
            while (*(end - 1) == ctype.widen('0')) {
                --end;
            }

            // Also clear the "trailing zeros to append during printing" range
            trailing_zeros_start = buffer.end();
            trailing_zeros_count = 0;
        }

        if (end - 1 == point && trailing_zeros_count == 0 && !showpoint) {
            // Remove the decimal point, too
            --end;
        }
    }

    // Apply thousands grouping
    const auto& grouping = numpunct.grouping();
    if (!grouping.empty())
    {
        // Step backwards from the end or decimal point, inserting the
        // thousands separator at every group interval.
        const CharT thousands_sep = ctype.widen(numpunct.thousands_sep());
        std::size_t group = 0;
        auto p = point != buffer.end() ? point : end;
        auto size = static_cast<int>(grouping[group]);
        while (size > 0 && size < CHAR_MAX && p - digits_start > size) {
            p -= size;
            insert_character(p, thousands_sep);
            if (group < grouping.size() - 1) {
                size = static_cast<int>(grouping[++group]);
            }
        }
    }

    // Print the exponent if required
    assert(floatfield != 0);
    if (floatfield & std::ios_base::scientific)
    {
        // Hexadecimal (%a/%A) or decimal (%e/%E) scientific notation
        if (floatfield & std::ios_base::fixed) {
            *end++ = ctype.widen(uppercase ? 'P' : 'p');
        } else {
            *end++ = ctype.widen(uppercase ? 'E' : 'e');
        }

        if (value.exponent < 0) {
            *end++ = ctype.widen('-');
            value.exponent = -value.exponent;
        } else {
            *end++ = ctype.widen('+');
        }

        if (floatfield == std::ios_base::scientific) {
            // In decimal scientific notation (%e/%E), the exponent is at least two digits
            if (value.exponent < 10) {
                *end++ = ctype.widen('0');
            }
        }

        const auto exponent_start = end;
        if (value.exponent == 0) {
            *end++ = ctype.widen('0');
        } else while (value.exponent > 0) {
            *end++ = ctype.widen(digits[value.exponent % 10]);
            value.exponent /= 10;
        }
        std::reverse(exponent_start, end);
    }

    // Write character `ch` `count` times to the stream
    const auto sputcn = [&](CharT ch, std::streamsize count){
        // Fill a buffer to output larger chunks
        constexpr std::streamsize chunk_size = 64;
        std::array<CharT, chunk_size> fill_buffer;
        std::fill_n(fill_buffer.begin(), std::min(count, chunk_size), ch);

        for (std::streamsize size, left = count; left > 0; left -= size) {
            size = std::min(chunk_size, left);
            os.rdbuf()->sputn(&fill_buffer[0], size);
        }
    };

    // Outputs a range of characters, making sure to output the trailing zeros range
    // if it lies in the specified range
    const auto put_range = [&](typename buffer_t::const_iterator begin, typename buffer_t::const_iterator end) {
        assert(end >= begin);
        if (trailing_zeros_start >= begin && trailing_zeros_start <= end) {
            // Print range with trailing zeros range in the middle
            assert(trailing_zeros_count > 0);
            os.rdbuf()->sputn(&*begin, trailing_zeros_start - begin);
            sputcn(ctype.widen('0'), trailing_zeros_count);
            os.rdbuf()->sputn(&*trailing_zeros_start, end - trailing_zeros_start);
        } else {
            // Print range as-is
            os.rdbuf()->sputn(&*begin, end - begin);
        }
    };

    // Pad the buffer if necessary.
    // Note that the length of trailing zeros is counted towards the length of the content.
    const auto content_size = end - buffer.begin() + trailing_zeros_count;
    if (content_size >= width)
    {
        // Buffer needs no padding, output as-is
        put_range(buffer.begin(), end);
    }
    else
    {
        const auto pad_size = width - content_size;
        switch (adjustfield)
        {
        case std::ios_base::left:
            // Content is left-aligned, so output the buffer, followed by the padding
            put_range(buffer.begin(), end);
            sputcn(os.fill(), pad_size);
            break;
        case std::ios_base::internal:
            // Content is internally aligned, so output the buffer up to the "internal pad"
            // point, followed by the padding, followed by the remainder of the buffer.
            put_range(buffer.begin(), internal_pad);
            sputcn(os.fill(), pad_size);
            put_range(internal_pad, end);
            break;
        default:
            // Content is right-aligned, so output the padding, followed by the buffer
            sputcn(os.fill(), pad_size);
            put_range(buffer.begin(), end);
            break;
        }
    }

    // Width is reset after every write
    os.width(0);

    return os;
}


template <typename CharT, class Traits, typename B, typename I, unsigned int F, bool R>
std::basic_istream<CharT, Traits>& operator>>(std::basic_istream<CharT, Traits>& is, fixed<B, I, F, R>& x)
{
    typename std::basic_istream<CharT, Traits>::sentry sentry(is);
    if (!sentry)
    {
        return is;
    }

    const auto& ctype = std::use_facet<std::ctype<CharT>>(is.getloc());
    const auto& numpunct = std::use_facet<std::numpunct<CharT>>(is.getloc());

    bool thousands_separator_allowed = false;
    const bool supports_thousands_separators = !numpunct.grouping().empty();

    const auto& is_valid_character = [](char ch) {
        // Note: allowing ['p', 'i', 'n', 't', 'y'] is technically in violation of the spec (we are emulating std::num_get),
        // but otherwise we cannot parse hexfloats and "infinity". This is a known issue with the spec (LWG #2381).
        return std::isxdigit(ch) ||
            ch == 'x' || ch == 'X' || ch == 'p' || ch == 'P' ||
            ch == 'i' || ch == 'I' || ch == 'n' || ch == 'N' ||
            ch == 't' || ch == 'T' || ch == 'y' || ch == 'Y' ||
            ch == '-' || ch == '+';
    };

    const auto& peek = [&]() {
        for(;;) {
            auto ch = is.rdbuf()->sgetc();
            if (ch == Traits::eof()) {
                is.setstate(std::ios::eofbit);
                return '\0';
            }
            if (ch == numpunct.decimal_point()) {
                return '.';
            }
            if (ch == numpunct.thousands_sep())
            {
                if (!supports_thousands_separators || !thousands_separator_allowed) {
                    return '\0';
                }
                // Ignore valid thousands separators
                is.rdbuf()->sbumpc();
                continue;
            }
            auto res = ctype.narrow(ch, 0);
            if (!is_valid_character(res)) {
                // Invalid character: end input
                return '\0';
            }
            return res;
        }
    };

    const auto& bump = [&]() {
        is.rdbuf()->sbumpc();
    };

    const auto& next = [&]() {
        bump();
        return peek();
    };

    bool negate = false;
    auto ch = peek();
    if (ch == '-') {
        negate = true;
        ch = next();
    } else if (ch == '+') {
        ch = next();
    }

    const char infinity[] = "infinity";
    // Must be "inf" or "infinity"
    int i = 0;
    while (i < 8 && ch == infinity[i]) {
        ++i;
        ch = next();
    }

    if (i > 0) {
        if (i == 3 || i == 8) {
            x = negate ? std::numeric_limits<fixed<B, I, F, R>>::min() : std::numeric_limits<fixed<B, I, F, R>>::max();
        } else {
            is.setstate(std::ios::failbit);
        }
        return is;
    }

    char exponent_char = 'e';
    int base = 10;

    constexpr auto NoFraction = std::numeric_limits<std::size_t>::max();
    std::size_t fraction_start = NoFraction;
    std::vector<unsigned char> significand;

    if (ch == '0') {
        ch = next();
        if (ch == 'x' || ch == 'X') {
            // Hexfloat
            exponent_char = 'p';
            base = 16;
            ch = next();
        } else {
            significand.push_back(0);
        }
    }

    // Parse the significand
    thousands_separator_allowed = true;
    for (;; ch = next()) {
        if (ch == '.') {
            if (fraction_start != NoFraction) {
                // Double decimal point. Stop parsing.
                break;
            }
            fraction_start = significand.size();
            thousands_separator_allowed = false;
        } else {
            unsigned char val = base;
            if (ch >= '0' && ch <= '9') {
                val = ch - '0';
            } else if (ch >= 'a' && ch <= 'f') {
                val = ch - 'a' + 10;
            } else if (ch >= 'A' && ch <= 'F') {
                val = ch - 'A' + 10;
            }
            if (val < 0 || val >= base) {
                break;
            }
            significand.push_back(val);
        }
    }
    if (significand.empty()) {
        // We need a significand
        is.setstate(std::ios::failbit);
        return is;
    }
    thousands_separator_allowed = false;

    if (fraction_start == NoFraction) {
        // If we haven't seen a fraction yet, place it at the end of the significand
        fraction_start = significand.size();
    }

    // Parse the exponent
    bool exponent_overflow = false;
    std::size_t exponent = 0;
    bool exponent_negate = false;
    if (std::tolower(ch) == exponent_char)
    {
        ch = next();
        if (ch == '-') {
            exponent_negate = true;
            ch = next();
        } else if (ch == '+') {
            ch = next();
        }

        bool parsed = false;
        while (std::isdigit(ch)) {
            if (exponent <= std::numeric_limits<int>::max() / 10) {
                exponent = exponent * 10 + (ch - '0');
            } else {
                exponent_overflow = true;
            }
            parsed = true;
            ch = next();
        }
        if (!parsed) {
            // If the exponent character is given, the exponent value may not be empty
            is.setstate(std::ios::failbit);
            return is;
        }
    }

    // We've parsed all we need. Construct the value.
    if (exponent_overflow) {
        // Absolute exponent is too large
        if (std::all_of(significand.begin(), significand.end(), [](unsigned char x){ return x == 0; })) {
            // Significand is zero. Exponent doesn't matter.
            x = fixed<B, I, F, R>(0);
        } else if (exponent_negate) {
            // A huge negative exponent approaches 0.
            x = fixed<B, I, F, R>::from_raw_value(0);
        } else {
            // A huge positive exponent approaches infinity.
            x = std::numeric_limits<fixed<B, I, F, R>>::max();
        }
        return is;
    }

    // Shift the fraction offset according to exponent
    {
        const auto exponent_mult = (base == 10) ? 1: 4;
        if (exponent_negate) {
            const auto adjust = std::min(exponent / exponent_mult, fraction_start);
            fraction_start -= adjust;
            exponent -= adjust * exponent_mult;
        } else {
            const auto adjust = std::min(exponent / exponent_mult, significand.size() - fraction_start);
            fraction_start += adjust;
            exponent -= adjust * exponent_mult;
        }
    }

    constexpr auto IsSigned = std::is_signed<B>::value;
    constexpr auto IntBits = sizeof(B) * 8 - F - (IsSigned ? 1 : 0);
    constexpr auto MaxInt = (I{1} << IntBits) - 1;
    constexpr auto MaxFraction = (I{1} << F) - 1;
    constexpr auto MaxValue = (I{1} << sizeof(B) * 8) - 1;

    // Parse the integer part
    I integer = 0;
    for (std::size_t i = 0; i < fraction_start; ++i) {
        if (integer > MaxInt / base) {
            // Overflow
            x = negate ? std::numeric_limits<fixed<B, I, F, R>>::min() : std::numeric_limits<fixed<B, I, F, R>>::max();
            return is;
        }
        assert(significand[i] < base);
        integer = integer * base + significand[i];
    }

    // Parse the fractional part
    I fraction = 0;
    I divisor = 1;
    for (std::size_t i = fraction_start; i < significand.size(); ++i) {
        assert(significand[i] < base);
        if (divisor > MaxFraction / base) {
            // We're done
            break;
        }
        fraction = fraction * base + significand[i];
        divisor *= base;
    }

    // Construct the value from the parsed parts
    I raw_value = (integer << F) + (fraction << F) / divisor;

    // Apply remaining exponent
    if (exponent_char == 'p') {
        // Base-2 exponent
        if (exponent_negate) {
            raw_value >>= exponent;
        } else {
            raw_value <<= exponent;
        }
    } else {
        // Base-10 exponent
        if (exponent_negate) {
            I remainder = 0;
            for (std::size_t e = 0; e < exponent; ++e) {
                remainder = raw_value % 10;
                raw_value /= 10;
            }
            raw_value += remainder / 5;
        } else {
            for (std::size_t e = 0; e < exponent; ++e) {
                if (raw_value > MaxValue / 10) {
                    // Overflow
                    x = negate ? std::numeric_limits<fixed<B, I, F, R>>::min() : std::numeric_limits<fixed<B, I, F, R>>::max();
                    return is;
                }
                raw_value *= 10;
            }
        }
    }
    x = fixed<B, I, F, R>::from_raw_value(static_cast<B>(negate ? -raw_value : raw_value));
    return is;
}


} // namespace navtools
#endif

