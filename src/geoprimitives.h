/* Copyright (c) 2013 Kevin L. Stern
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#pragma once

#include <bitset>
#include <cmath>
#include <cstdint>
#include <istream>
#include <limits>
#include <vector>
#include <stdexcept>

// -------------------------------------------------------------------------------------------------

// A CoordinateSystem is a Cartesian coordinate system with custom x and y coordinate extents.
// Following are pre-defined CoordinateSystems. A custom CoordinateSystem can be defined by creating
// a class with a static function GetRanges which initializes an x and y coordinate extent array.

// The GeoHash and related primitives are templated so that their CoordinateSystem is part of their
// type. This CoordinateSystem gives latitude/longitude in degrees.
class GeoDegreeCoordinateSystem {
public:
  static void GetRanges(double xrange[2], double yrange[2]) {
    xrange[0] = -180;
    xrange[1] = 180;
    yrange[0] = -90;
    yrange[1] = 90;
  }
};

// The GeoHash and related primitives are templated so that their CoordinateSystem is part of their
// type. This CoordinateSystem gives latitude/longitude in radians.
class GeoRadianCoordinateSystem {
public:
  static void GetRanges(double xrange[2], double yrange[2]) {
    xrange[0] = -M_PI;
    xrange[1] = M_PI;
    yrange[0] = -M_PI_2;
    yrange[1] = M_PI_2;
  }
};

// -------------------------------------------------------------------------------------------------

// A GeoPoint is a two dimensional vector of type double which represents a point in
// CoordinateSystem.
template <typename CoordinateSystem>
class GeoPoint {
public:
  GeoPoint() : x_(0), y_(0) {}

  GeoPoint(double x, double y) : x_(x), y_(y) {}

  GeoPoint(const GeoPoint<CoordinateSystem>& point) : x_(point.x_), y_(point.y_) {}

  GeoPoint<CoordinateSystem>& operator=(const GeoPoint<CoordinateSystem>& point) {
    x_ = point.x_;
    y_ = point.y_;
    return *this;
  }

  bool operator==(const GeoPoint<CoordinateSystem>& point) const {
    return x_ == point.x_ && y_ == point.y_;
  }

  // Get the x coordinate.
  double x() const {
    return x_;
  }

  // Get the y coordinate.
  double y() const {
    return y_;
  }

private:
  double x_, y_;

  // For pretty printing.
  friend std::ostream& operator<<(std::ostream& out, const GeoPoint<CoordinateSystem>& point) {
    out << "(" << point.x_ << "," << point.y_ << ")";
    return out;
  }
};

namespace std {
template <typename CoordinateSystem>
class hash<GeoPoint<CoordinateSystem>> {
public:
  size_t operator()(const GeoPoint<CoordinateSystem>& point) const {
    size_t result = 31;
    result = 31 * result + std::hash<double>()(point.x());
    result = 31 * result + std::hash<double>()(point.y());
    return result;
  }
};
}

// -------------------------------------------------------------------------------------------------

// A GeoRectangle represents the rectangle defined by a lower-left and an upper-right GeoPoint
// within CoordinateSystem.
template <typename CoordinateSystem>
class GeoRectangle {
public:
  GeoRectangle() {}

  GeoRectangle(double ll_x, double ll_y, double ur_x, double ur_y)
      : ll_(ll_x, ll_y), ur_(ur_x, ur_y) {}

  GeoRectangle(GeoPoint<CoordinateSystem> ll, GeoPoint<CoordinateSystem> ur) : ll_(ll), ur_(ur) {}

  GeoRectangle(const GeoRectangle<CoordinateSystem>& rect) : ll_(rect.ll_), ur_(rect.ur_) {}

  GeoRectangle<CoordinateSystem>& operator=(const GeoRectangle<CoordinateSystem>& rect) {
    ll_ = rect.ll_;
    ur_ = rect.ur_;
    return *this;
  }

  bool operator==(const GeoRectangle<CoordinateSystem>& rect) const {
    return ll_ == rect.ll_ && ur_ == rect.ur_;
  }

  // Get the lower-left GeoPoint.
  const GeoPoint<CoordinateSystem>& ll() const {
    return ll_;
  }

  // Get the upper-right GeoPoint.
  const GeoPoint<CoordinateSystem>& ur() const {
    return ur_;
  }

  // Test whether or not this rectangle contains point.
  bool contains(const GeoPoint<CoordinateSystem>& point) const {
    return ll_.x() <= point.x() && point.x() <= ur_.x()
        && ll_.y() <= point.y() && point.y() <= ur_.y();
  }

  // Test whether or not this rectangle intersects rect.
  bool intersects(const GeoRectangle<CoordinateSystem> rect) const {
    return !(rect.ll_.x() > ur_.x() || rect.ur_.x() < ll_.x()
        || rect.ll_.y() > ur_.y() || rect.ur_.y() < ll_.y());
  }

private:
  GeoPoint<CoordinateSystem> ll_, ur_;

  // For pretty printing.
  friend std::ostream& operator<<(std::ostream& out, const GeoRectangle<CoordinateSystem>& rect) {
    out << "Rectangle[" << rect.ll_ << "," << rect.ur_ << "]";
    return out;
  }
};

namespace std {
template <typename CoordinateSystem>
class hash<GeoRectangle<CoordinateSystem>> {
public:
  size_t operator()(const GeoRectangle<CoordinateSystem>& rect) const {
    size_t result = 31;
    result = 31 * result + std::hash<GeoPoint<CoordinateSystem>>()(rect.ll());
    result = 31 * result + std::hash<GeoPoint<CoordinateSystem>>()(rect.ur());
    return result;
  }
};
}

// -------------------------------------------------------------------------------------------------

// GeoHash wraps a HashType geohash value within CoordinateSystem.
template <typename CoordinateSystem, typename HashType>
class GeoHash {
public:
  GeoHash(HashType value) : value_(value) {}

  GeoHash(const GeoHash& hash) : value_(hash.value_) {}

  GeoHash<CoordinateSystem, HashType>& operator=(const GeoHash<CoordinateSystem, HashType>& hash) {
    value_ = hash.value_;
    return *this;
  }

  bool operator==(const GeoHash<CoordinateSystem, HashType>& hash) const {
    return value_ == hash.value_;
  }

  // Get the hash value.
  HashType value() const {
    return value_;
  }

private:
  HashType value_;
};

namespace std {
template <typename CoordinateSystem, typename HashType>
class hash<GeoHash<CoordinateSystem, HashType>> {
public:
  size_t operator()(const GeoHash<CoordinateSystem, HashType>& hash) const {
    size_t result = 31;
    result = 31 * result + std::hash<HashType>()(hash.value());
    return result;
  }
};
}

// Used for mapping from hash type to coordinate type and for compile-time validation of the hash
// type:
// A 64 bit hash type carries 32 bits from each coordinate.
// A 32 bit hash type carries 16 bits from each coordinate.
// A 16 bit hash type carries 8 bits from each coordinate.
template <typename HashType>
struct HashTypeInterpreter;

template <>
struct HashTypeInterpreter<uint64_t> {
  typedef uint32_t CoordinateType;
};

template <>
struct HashTypeInterpreter<uint32_t> {
  typedef uint16_t CoordinateType;
};

template <>
struct HashTypeInterpreter<uint16_t> {
  typedef uint8_t CoordinateType;
};

// Maps and unmaps points to distances down the Hilbert curve of #bits(sizeof(CoordinateType))
// iterations.
//
// To map a point to its Hilbert value, this implementation operates iteratively by dividing the
// current area into quadrants, calculating the quadrant into which the point falls, and rotating
// the quadrant ordering to match the current Hilbert curve iteration. This process repeats with the
// calculated quadrant for #bits(sizeof(CoordinateType)) iterations.
//
// To map a Hilbert value to a point, this implementation operates iteratively by extracting the
// rotated quadrant from the appropriate two-bits of the Hilbert value and un-rotating based upon
// the Hilbert curve iteration which produced the bits. This process repeats for
// #bits(sizeof(CoordinateType)) iterations.
//
// See http://en.wikipedia.org/wiki/Hilbert_curve for the construction of a Hilbert curve.
template <typename HashType>
class Hilbert {
private:
  typedef typename HashTypeInterpreter<HashType>::CoordinateType CoordinateType;

public:
  enum Rotation {Up, Right, Down, Left};

  HashType hilbert(CoordinateType x, CoordinateType y) const {
    return hilbert(x, y, nullptr);
  }

  // Map from the point (x, y) to the length down the Hilbert curve of #bits(sizeof(CoordinateType))
  // iterations.
  // rotations: #bits(sizeof(CoordinateType)) Rotation wide array which will be populated with the
  //            quadrant rotations at each iteration of hilbert curve generation; ignored if
  //            nullptr.
  HashType hilbert(CoordinateType x, CoordinateType y, Rotation* rotations) const {
    HashType result = 0;
    Rotation rotation = Up;
    for (int i = (sizeof(CoordinateType) << 3) - 1; i >= 0; --i) {
      uint8_t quadrant = 0;
      // | 0 | 3 |
      //  -------
      // | 1 | 2 |
      if (x & (1 << i)) {
        if (y & (1 << i)) {
          quadrant = 3;
        } else {
          quadrant = 2;
        }
      } else {
        if (y & (1 << i)) {
          quadrant = 0;
        } else {
          quadrant = 1;
        }
      }

      if (rotations != nullptr) {
        rotations[i] = rotation;
      }
      HashType rotated = rotate(rotation, quadrant);
      rotation = new_rotation(rotation, rotated);
      result |= (rotated << (i << 1));
    }
    return result;
  }

  void unhilbert(HashType d, CoordinateType& x, CoordinateType& y) const {
    unhilbert(d, x, y, nullptr);
  }

  // Map from d, the output of the hilbert function, to the point (x, y).
  // rotations: #bits(sizeof(CoordinateType)) Rotation wide array which will be populated with the
  //            quadrant rotations at each iteration of hilbert curve generation; ignored if
  //            nullptr.
  void unhilbert(HashType d, CoordinateType& x, CoordinateType& y, Rotation* rotations) const {
    x = 0;
    y = 0;
    Rotation rotation = Up;
    for (int i = (sizeof(CoordinateType) << 3) - 1; i >= 0; --i) {
      if (rotations != nullptr) {
        rotations[i] = rotation;
      }
      uint8_t quadrant = ((d & (3ll << (i << 1))) >> (i << 1));
      uint8_t unrotated = rotate(rotation, quadrant);
      rotation = new_rotation(rotation, quadrant);
      // | 0 | 3 |
      //  -------
      // | 1 | 2 |
      switch (unrotated) {
        case 0:
          y |= (1 << i);
          break;
        case 1:
          // no op
          break;
        case 2:
          x |= (1 << i);
          break;
        case 3:
          x |= (1 << i);
          y |= (1 << i);
          break;
        default:
          throw std::runtime_error("internal error");
      }
    }
  }

private:
  // Rotate quadrant by rotation. The provided quadrant is with respect to the Up rotation.
  //
  // Note that this function is self-inverting in the sense that rotate(r, rotate(r, x)) == x.
  uint8_t rotate(Rotation rotation, uint8_t quadrant) const {
    // | 2 | 3 |
    //  -------
    // | 1 | 0 |
    static uint8_t right[4] = {2, 1, 0, 3};
    // | 2 | 1 |
    //  -------
    // | 3 | 0 |
    static uint8_t down[4] = {2, 3, 0, 1};
    // | 0 | 1 |
    //  -------
    // | 3 | 2 |
    static uint8_t left[4] = {0, 3, 2, 1};

    switch (rotation) {
      case Up:
        return quadrant;
      case Right:
        return right[quadrant];
      case Down:
        return down[quadrant];
      case Left:
        return left[quadrant];
      default:
        throw std::runtime_error("internal error");
    }
  }

  // Compute the new rotation for quadrant and rotation. The provided quadrant is with respect to
  // rotation.
  Rotation new_rotation(Rotation rotation, uint8_t quadrant) const {
    Rotation result = rotation;
    // | Left | Right |
    //  --------------
    // |  Up  |   Up  |
    if (quadrant == 0) {
      switch (rotation) {
        case Up:
          result = Left;
          break;
        case Right:
          result = Down;
          break;
        case Down:
          result = Right;
          break;
        case Left:
          result = Up;
          break;
      }
    } else if (quadrant == 3) {
      switch (rotation) {
        case Up:
          result = Right;
          break;
        case Right:
          result = Up;
          break;
        case Down:
          result = Left;
          break;
        case Left:
          result = Down;
          break;
      }
    }
    return result;
  }
};

// GeoHasher hashes and unhashes GeoPoints. The GeoHash produced by geohashing a point is computed
// as the distance down the Hilbert curve of #bits(sizeof(CoordinateType)) iterations overlaid onto
// the CoordinateSystem's range. Each coordinate is first mapped to a #bits(sizeof(CoordinateType))
// bit integer which represents the coordinate's location within the CoordinateSystem's range as
// follows: The high-order bit gives the coordinate's location relative to the midpoint of the
// range, the next bit gives the coordinate's location relative to the the midpoint of the relevant
// half of the range, and so on. The results of this process are the coordinates within an integer
// coordinate system that spans the CoordinateSystem's range onto which a Hilbert curve can be
// overlaid.
//
// The construction of the Hilbert curve gives the nice property that points which are near one
// another spatially tend to fall near one another on the curve; hence, points that are near one
// another spatially tend to produce GeoHashes that are close in magnitude.
template <typename HashType>
class GeoHasher {
private:
  typedef typename HashTypeInterpreter<HashType>::CoordinateType CoordinateType;

public: 
  // Compute the GeoHash for point.
  template <typename CoordinateSystem>
  GeoHash<CoordinateSystem, HashType> hash(const GeoPoint<CoordinateSystem>& point) const {
    double x_range[2];
    double y_range[2];
    CoordinateSystem::GetRanges(x_range, y_range);

    CoordinateType x = 0;
    CoordinateType y = 0;
    for (int i = (sizeof(CoordinateType) << 3) - 1; i >= 0; --i) {
      double y_mid = (y_range[0] + y_range[1]) / 2.;
      if (point.y() >= y_mid) {
        y |= 1 << i;
        y_range[0] = y_mid;
      } else {
        y_range[1] = y_mid;
      }
      double x_mid = (x_range[0] + x_range[1]) / 2.;
      if (point.x() >= x_mid) {
        x |= 1 << i;
        x_range[0] = x_mid;
      } else {
        x_range[1] = x_mid;
      }
    }
    return hilbert_.hilbert(x, y);
  }

  // Compute the GeoPoint for hash.
  template <typename CoordinateSystem>
  GeoPoint<CoordinateSystem> unhash(const GeoHash<CoordinateSystem, HashType>& hash) const {
    double x_range[2];
    double y_range[2];
    CoordinateSystem::GetRanges(x_range, y_range);

    CoordinateType x = 0;
    CoordinateType y = 0;
    hilbert_.unhilbert(hash.value(), x, y);

    for (int i = (sizeof(CoordinateType) << 3) - 1; i >= 0; --i) {
      double y_mid = (y_range[0] + y_range[1]) / 2.;
      if (y & (1 << i)) {
        y_range[0] = y_mid;
      } else {
        y_range[1] = y_mid;
      }
      double x_mid = (x_range[0] + x_range[1]) / 2.;
      if (x & (1 << i)) {
        x_range[0] = x_mid;
      } else {
        x_range[1] = x_mid;
      }
    }
    return GeoPoint<CoordinateSystem>(x_range[0], y_range[0]);
  }

private:
  Hilbert<HashType> hilbert_;
};

// GeoPrefix represents a GeoHash prefix with a hash value and an offset which gives the number of
// low-order bits that are not part of the prefix.
// The GeoHash itself represents a distance down the Hilbert curve of #bits(sizeof(CoordinateType))
// iterations within CoordinateSystem; distance values that share a prefix are no farther from one
// another than the width of the prefix range.
//
// GeoPrefix provides a convenience method for finding the joint prefix for two GeoHashes and other
// functionality such as giving the bounding box associated with the prefix.
template <typename CoordinateSystem, typename HashType>
class GeoPrefix {
private:
  typedef typename HashTypeInterpreter<HashType>::CoordinateType CoordinateType;

public:
  // Calculate the joint GeoPrefix for h1 and h2.
  static GeoPrefix<CoordinateSystem, HashType> joint_prefix(
      GeoHash<CoordinateSystem, HashType> h1, GeoHash<CoordinateSystem, HashType> h2) {
    HashType value = 0;
    uint8_t offset = (sizeof(HashType) << 3);
    for (int i = 0; i < (sizeof(HashType) << 3); ++i) {
      HashType mask = 1ll << ((sizeof(HashType) << 3) - 1 - i);
      if ((h1.value() & mask) != (h2.value() & mask)) {
        break;
      }
      value |= (h1.value() & mask);
      --offset;
    }
    return GeoPrefix<CoordinateSystem, HashType>(value, offset);
  }

  // Compute a collection of prefixes such that any point within bounds will exhibit a prefix from
  // the collection.
  // prefix_refinement_stop_size: The prefix count at which to stop refining query prefixes.
  static void search_prefixes(const GeoRectangle<CoordinateSystem>& bounds,
                              std::vector<GeoPrefix<CoordinateSystem, HashType>>& prefixes,
                              size_t prefix_refinement_stop_size) {
    GeoHasher<HashType> hasher;
    // The joint prefix between bounds.ll() and bounds.ur() gives a quadrant that contains both
    // points and, therefore, all points within bounds.
    prefixes.push_back(joint_prefix(hasher.hash(bounds.ll()), hasher.hash(bounds.ur())));
    std::vector<GeoPrefix<CoordinateSystem, HashType>> workspace;
    // Split the joint prefix and filter out unnecessary search areas.
    bool stop_splitting = false;
    for (int i = 0; !stop_splitting; ++i) {
      size_t prefix_count = prefixes.size();
      for (auto& prefix : prefixes) {
        if (stop_splitting || prefix.offset() == 0) {
          workspace.push_back(prefix);
          // If we're here because prefix.offset() == 0, then no additional splitting will occur.
          stop_splitting = true;
        } else {
          uint8_t new_offset = prefix.offset() - 1;
          --prefix_count;
          // Extend prefix by one bit and select the forks that intersect bounds. Heuristically,
          // this will eliminate an unnecessary half of a prefix that is too coarse.
          for (HashType i : {0, 1}) {
            GeoPrefix<CoordinateSystem, HashType> split(prefix.value() | (i << new_offset),
                                                        new_offset);
            if (bounds.intersects(split.bounds())) {
              workspace.push_back(split);
              ++prefix_count;
            }
          }
          if (prefix_count >= prefix_refinement_stop_size) {
            stop_splitting = true;
          }
        }
      }
      prefixes.swap(workspace);
      workspace.clear();
    }
    // Recombine prefixes that were split but not filtered due to both halves intersecting bounds.
    while (true) {
      int i;
      for (i = 0; i < prefixes.size() - 1; ++i) {
        // {prefix}0 and {prefix}1 will be adjacent, if they are both present.
        GeoPrefix<CoordinateSystem, HashType> prefix(prefixes[i].value(), prefixes[i].offset() + 1);
        if (prefixes[i].offset() == prefixes[i + 1].offset() && prefix.applies_to(prefixes[i + 1])) {
          workspace.push_back(prefix);
          ++i;
        } else {
          workspace.push_back(prefixes[i]);
        }
      }
      if (i < prefixes.size()) {
        workspace.push_back(prefixes[i]);
      }
      if (prefixes.size() == workspace.size()) {
        break;
      }
      prefixes.swap(workspace);
      workspace.clear();
    }
  }

  // value is a GeoHash value and offset gives the number of low-order bits to ignore.
  GeoPrefix(HashType value, uint8_t offset) : value_(value), offset_(offset) {
    if (offset == (sizeof(HashType) << 3)) {
      value_ = 0;
    } else {
      value_ &= static_cast<HashType>(0xffffffffffffffffll << offset_);
    }
  }

  GeoPrefix(const GeoPrefix<CoordinateSystem, HashType>& geo_prefix)
      : value_(geo_prefix.value_), offset_(geo_prefix.offset_) {}

  bool operator==(const GeoPrefix<CoordinateSystem, HashType>& prefix) const {
    return value_ == prefix.value_ && offset_ == prefix.offset_;
  }

  // Get the value.
  HashType value() const {
    return value_;
  }

  // Get the offset.
  uint8_t offset() const {
    return offset_;
  }

  // Test whether or not this prefix is a prefix of hash.
  bool applies_to(GeoHash<CoordinateSystem, HashType> hash) const {
    if (offset_ == (sizeof(HashType) << 3)) {
      // The C++ standard leaves bit-shifting by a value greater than the type's width undefined;
      // therefore, explicitly defining this behavior is required.
      return true;
    }
    return (hash.value() & static_cast<HashType>(0xffffffffffffffffll << offset_)) == value_;
  }

  // Test whether or not this prefix is a prefix of prefix.
  bool applies_to(const GeoPrefix<CoordinateSystem, HashType>& prefix) const {
    if (offset_ == (sizeof(HashType) << 3)) {
      // The C++ standard leaves bit-shifting by a value greater than the type's width undefined;
      // therefore, explicitly defining this behavior is required.
      return true;
    }
    return offset_ >= prefix.offset_ &&
           (prefix.value() & static_cast<HashType>(0xffffffffffffffffll << offset_)) == value_;
  }

  // Calculate the bounds for this prefix. The bounds for a prefix is the rectangle that delimits
  // the area in which all points have this prefix and out of which no points have this prefix.
  GeoRectangle<CoordinateSystem> bounds() const {
    double x_range[2];
    double y_range[2];
    CoordinateSystem::GetRanges(x_range, y_range);
    if (offset_ < (sizeof(HashType) << 3)) {
      CoordinateType x = 0;
      CoordinateType y = 0;
      typename TypedHilbert::Rotation rotations[32];
      hilbert_.unhilbert(value_, x, y, rotations);
      int i = (sizeof(HashType) << 3) - 1;
      while (i > offset_) {
        if (y & (1 << (i >> 1))) {
          y_range[0] = (y_range[0] + y_range[1]) / 2.;
        } else {
          y_range[1] = (y_range[0] + y_range[1]) / 2.;
        }
        --i;
        if (x & (1 << (i >> 1))) {
          x_range[0] = (x_range[0] + x_range[1]) / 2.;
        } else {
          x_range[1] = (x_range[0] + x_range[1]) / 2.;
        }
        --i;
      }
      if (i == offset_) {
        switch (rotations[(i >> 1)]) {
          case TypedHilbert::Up:
          case TypedHilbert::Down:
            if (x & (1 << (i >> 1))) {
              x_range[0] = (x_range[0] + x_range[1]) / 2.;
            } else {
              x_range[1] = (x_range[0] + x_range[1]) / 2.;
            }
            break;
          case TypedHilbert::Left:
          case TypedHilbert::Right:
            if (y & (1 << (i >> 1))) {
              y_range[0] = (y_range[0] + y_range[1]) / 2.;
            } else {
              y_range[1] = (y_range[0] + y_range[1]) / 2.;
            }
            break;
          default:
            throw std::runtime_error("internal error");
        }
      }
    }
    return GeoRectangle<CoordinateSystem>(x_range[0], y_range[0], x_range[1], y_range[1]);
  }

private:
  typedef Hilbert<HashType> TypedHilbert;

  TypedHilbert hilbert_;
  HashType value_;
  uint8_t offset_;

  // For pretty printing.
  friend std::ostream& operator<<(std::ostream& out,
                                  const GeoPrefix<CoordinateSystem, HashType>& prefix) {
    out << "value=" << std::bitset<(sizeof(HashType) << 3)>(prefix.value_) << " offset="
        << (int) prefix.offset_;
    return out;
  }
};

namespace std {
template <typename CoordinateSystem, typename HashType>
class hash<GeoPrefix<CoordinateSystem, HashType>> {
public:
  size_t operator()(const GeoPrefix<CoordinateSystem, HashType>& prefix) const {
    size_t result = 31;
    result = 31 * result + std::hash<HashType>()(prefix.value());
    result = 31 * result + std::hash<uint8_t>()(prefix.offset());
    return result;
  }
};
}
