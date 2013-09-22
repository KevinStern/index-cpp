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

#include <algorithm>
#include <bitset>
#include <random>
#include <unordered_set>

#include "geoindex.h"

static double XMIN = -180;
static double XMAX = 180;
static double YMIN = -90;
static double YMAX = 90;

typedef GeoPoint<GeoDegreeCoordinateSystem> AnchoredGeoPoint;
typedef GeoRectangle<GeoDegreeCoordinateSystem> AnchoredGeoRectangle;
typedef GeoHash<GeoDegreeCoordinateSystem> AnchoredGeoHash;
typedef GeoPrefix<GeoDegreeCoordinateSystem> AnchoredGeoPrefix;

static const char TEST_FILENAME[] = "GeoIndexTest.index";
static std::default_random_engine E;
static std::uniform_real_distribution<double> X{XMIN, XMAX};
static std::uniform_real_distribution<double> Y{YMIN, YMAX};

struct IndexRecord {
  AnchoredGeoPoint point;

  IndexRecord() {}
  IndexRecord(const AnchoredGeoPoint& point) : point(point) {}
};

typedef GeoIndex<IndexRecord, GeoDegreeCoordinateSystem> AnchoredGeoIndex;

class IndexRecordSerializer {
public:
  static constexpr uint32_t ENCODED_SIZE = 16;

  static void parse(Buffer& buffer, IndexRecord& record) {
    double x = buffer.get_double();
    double y = buffer.get_double();
    record.point = AnchoredGeoPoint(x, y);
  }

  static void write(const IndexRecord& record, Buffer& buffer) {
    buffer.put_double(record.point.x());
    buffer.put_double(record.point.y());
  }

  static AnchoredGeoPoint geopoint(const IndexRecord& record) {
    return record.point;
  }
};

namespace {
AnchoredGeoRectangle MakeRandomGeoRectangle(double min_dimension, double max_dimension) {
  std::uniform_real_distribution<double> d{min_dimension, max_dimension};
  double x1 = X(E);
  double y1 = Y(E);
  double x2 = x1 + d(E);
  double y2 = y1 + d(E);
  return AnchoredGeoRectangle(x1, y1, x2, y2);
}
}

TEST(GeoIndexBasic) {
  E.seed(time(NULL));
  Cleaner cleaner([&](){remove(TEST_FILENAME);});
  std::vector<IndexRecord> points;
  for (int i = 0; i < 1000000; ++i) {
    points.push_back(IndexRecord(AnchoredGeoPoint(X(E), Y(E))));
  }
  AnchoredGeoIndex::order(points.begin(), points.end(), &IndexRecordSerializer::geopoint, points);
  Buffer buffer(Buffer::LITTLE, points.size() * IndexRecordSerializer::ENCODED_SIZE);
  for (const auto& next : points) {
    IndexRecordSerializer::write(next, buffer);
  }   
  buffer.switch_mode();
  FILE* file = fopen(TEST_FILENAME, "wb");
  fwrite(buffer.data(), 1, buffer.remaining(), file);
  fclose(file);
  AnchoredGeoIndex::FileDataSource file_source(
      TEST_FILENAME, &IndexRecordSerializer::parse, IndexRecordSerializer::ENCODED_SIZE,
      Buffer::LITTLE, 100);
  AnchoredGeoIndex::VectorDataSource vector_source(points);
  for (int i = 0; i < 100; ++i) {
    AnchoredGeoIndex file_index(&file_source, &IndexRecordSerializer::geopoint);
    AnchoredGeoIndex vector_index(&vector_source, &IndexRecordSerializer::geopoint);
    AnchoredGeoRectangle rect = MakeRandomGeoRectangle(0.01, 2);

    std::unordered_set<AnchoredGeoPoint> expected;
    for (auto& next : points) {
      if (rect.contains(next.point)) {
        expected.insert(next.point);
      }
    }

    std::unordered_set<AnchoredGeoPoint> file_actual;
    std::unordered_set<AnchoredGeoPoint> vector_actual;
    file_index.query(rect, [&](const IndexRecord& record){file_actual.insert(record.point);});
    vector_index.query(rect, [&](const IndexRecord& record){vector_actual.insert(record.point);});
    ASSERT_TRUE(expected == file_actual);
    ASSERT_TRUE(expected == vector_actual);
  }
}

TEST(GeoIndexWorstCase) {
  E.seed(time(NULL));
  Cleaner cleaner([&](){remove(TEST_FILENAME);});
  std::vector<IndexRecord> points;
  for (int i = 0; i < 1000000; ++i) {
    points.push_back(IndexRecord(AnchoredGeoPoint(X(E), Y(E))));
  }
  AnchoredGeoIndex::order(points.begin(), points.end(), &IndexRecordSerializer::geopoint, points);
  Buffer buffer(Buffer::LITTLE, points.size() * IndexRecordSerializer::ENCODED_SIZE);
  for (const auto& next : points) {
    IndexRecordSerializer::write(next, buffer);
  }   
  buffer.switch_mode();
  FILE* file = fopen(TEST_FILENAME, "wb");
  fwrite(buffer.data(), 1, buffer.remaining(), file);
  fclose(file);
  AnchoredGeoIndex::FileDataSource source(
      TEST_FILENAME, &IndexRecordSerializer::parse, IndexRecordSerializer::ENCODED_SIZE,
      Buffer::LITTLE, 100);
  AnchoredGeoIndex index(&source, &IndexRecordSerializer::geopoint);
  AnchoredGeoRectangle rect(-1, -1, 1, 1);

  std::unordered_set<AnchoredGeoPoint> expected;
  for (auto& next : points) {
    if (rect.contains(next.point)) {
      expected.insert(next.point);
    }
  }

  std::unordered_set<AnchoredGeoPoint> actual;
  index.query(rect, [&](const IndexRecord& record){actual.insert(record.point);});
  ASSERT_TRUE(expected == actual);
}
