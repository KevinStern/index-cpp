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

#include "test.h"

#include "geoprimitives.h"

typedef GeoPoint<GeoDegreeCoordinateSystem> AnchoredGeoPoint;
typedef GeoRectangle<GeoDegreeCoordinateSystem> AnchoredGeoRectangle;
typedef GeoHash<GeoDegreeCoordinateSystem> AnchoredGeoHash;
typedef GeoPrefix<GeoDegreeCoordinateSystem> AnchoredGeoPrefix;

TEST(GeoHashUnhash) {
  GeoHasher hasher;
  AnchoredGeoPoint p, q;
  p = AnchoredGeoPoint(0, 0);
  q = hasher.unhash(hasher.hash(p));
  ASSERT_EQ(p.x(), q.x(), 0.000000001);
  ASSERT_EQ(p.y(), q.y(), 0.000000001);

  p = AnchoredGeoPoint(90, 90);
  q = hasher.unhash(hasher.hash(p));
  ASSERT_EQ(p.x(), q.x(), 0.000000001);
  ASSERT_EQ(p.y(), q.y(), 0.000000001);

  p = AnchoredGeoPoint(-90, 90);
  q = hasher.unhash(hasher.hash(p));
  ASSERT_EQ(p.x(), q.x(), 0.000000001);
  ASSERT_EQ(p.y(), q.y(), 0.000000001);

  p = AnchoredGeoPoint(90, -90);
  q = hasher.unhash(hasher.hash(p));
  ASSERT_EQ(p.x(), q.x(), 0.000000001);
  ASSERT_EQ(p.y(), q.y(), 0.000000001);

  p = AnchoredGeoPoint(-90, -90);
  q = hasher.unhash(hasher.hash(p));
  ASSERT_EQ(p.x(), q.x(), 0.000000001);
  ASSERT_EQ(p.y(), q.y(), 0.000000001);
}

TEST(GeoPrefixJointPrefix) {
  ASSERT_EQ(AnchoredGeoPrefix(8, 3), AnchoredGeoPrefix::joint_prefix(8, 12));
  ASSERT_EQ(AnchoredGeoPrefix(0, 64), AnchoredGeoPrefix::joint_prefix(0, 1ll << 63));
}

TEST(GeoPrefixAppliesTo) {
  for (int i = 8; i < 16; ++i) {
    ASSERT_TRUE(AnchoredGeoPrefix(8, 3).applies_to(i));
  }
  ASSERT_FALSE(AnchoredGeoPrefix(8, 3).applies_to(16));
  ASSERT_TRUE(AnchoredGeoPrefix(0, 64).applies_to(16));
}

TEST(Hilbert) {
  Hilbert hilbert;
  // |  0 |  1 | 14 | 15 |
  //  -------------------
  // |  3 |  2 | 13 | 12 |
  //  -------------------
  // |  4 |  7 |  8 | 11 |
  //  -------------------
  // |  5 |  6 |  9 | 10 |
  ASSERT_EQ(0, hilbert.hilbert(0, 0xffffffff));
  ASSERT_EQ(1, hilbert.hilbert(0x1, 0xffffffff));
  ASSERT_EQ(2, hilbert.hilbert(0x1, 0xfffffffe));
  ASSERT_EQ(3, hilbert.hilbert(0, 0xfffffffe));
  ASSERT_EQ(4, hilbert.hilbert(0, 0xfffffffd));
  ASSERT_EQ(5, hilbert.hilbert(0, 0xfffffffc));
  ASSERT_EQ(6, hilbert.hilbert(0x1, 0xfffffffc));
  ASSERT_EQ(7, hilbert.hilbert(0x1, 0xfffffffd));
  ASSERT_EQ(8, hilbert.hilbert(0x2, 0xfffffffd));
  ASSERT_EQ(9, hilbert.hilbert(0x2, 0xfffffffc));
  ASSERT_EQ(10, hilbert.hilbert(0x3, 0xfffffffc));
  ASSERT_EQ(11, hilbert.hilbert(0x3, 0xfffffffd));
  ASSERT_EQ(12, hilbert.hilbert(0x3, 0xfffffffe));
  ASSERT_EQ(13, hilbert.hilbert(0x2, 0xfffffffe));
  ASSERT_EQ(14, hilbert.hilbert(0x2, 0xffffffff));
  ASSERT_EQ(15, hilbert.hilbert(0x3, 0xffffffff));
}

TEST(UnHilbert) {
  Hilbert hilbert;
  // Try out a bunch of points.
  for (int i = 0; i < 0x000000ff; ++i) {
    for (int j = 0; j < 0x000000ff; ++j) {
      uint32_t x, y;
      hilbert.unhilbert(hilbert.hilbert(i, j), x, y);
      ASSERT_EQ(i, x);
      ASSERT_EQ(j, y);
    }
  }
  // Try out a bunch more points.
  for (int i = 0; i < 0x000000ff; ++i) {
    for (int j = 0xffffff00; j < 0xffffffff; ++j) {
      uint32_t x, y;
      hilbert.unhilbert(hilbert.hilbert(i, j), x, y);
      ASSERT_EQ(i, x);
      ASSERT_EQ(j, y);
    }
  }
  // Even more ...
  for (int i = 0xffffff00; i < 0xffffffff; ++i) {
    for (int j = 0; j < 0x000000ff; ++j) {
      uint32_t x, y;
      hilbert.unhilbert(hilbert.hilbert(i, j), x, y);
      ASSERT_EQ(i, x);
      ASSERT_EQ(j, y);
    }
  }
  // This is it, really ...
  for (int i = 0xffffff00; i < 0xffffffff; ++i) {
    for (int j = 0xffffff00; j < 0xffffffff; ++j) {
      uint32_t x, y;
      hilbert.unhilbert(hilbert.hilbert(i, j), x, y);
      ASSERT_EQ(i, x);
      ASSERT_EQ(j, y);
    }
  }
}
